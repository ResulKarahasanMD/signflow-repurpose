"""Remote resource verification for Phase 1."""

from __future__ import annotations

import time
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from typing import Any

import requests

DEFAULT_TIMEOUT = 20
PACKAGE_VERSION_NAMES = [
    "decoupler",
    "gseapy",
    "matplotlib",
    "networkx",
    "numpy",
    "pandas",
    "pyarrow",
    "pydeseq2",
    "pyyaml",
    "requests",
    "scikit-learn",
    "scipy",
    "statsmodels",
    "torch",
    "torch-geometric",
]


@dataclass
class ResourceRecord:
    name: str
    ok: bool
    status: str
    endpoint: str
    checked_at_utc: str
    details: dict[str, Any] = field(default_factory=dict)
    warning: str | None = None
    error: str | None = None
    latency_s: float | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def utc_now_iso() -> str:
    return datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def _request(
    session: requests.Session,
    method: str,
    url: str,
    *,
    retries: int = 3,
    backoff_seconds: float = 1.0,
    timeout: int = DEFAULT_TIMEOUT,
    **kwargs: Any,
) -> requests.Response:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            response = session.request(method, url, timeout=timeout, **kwargs)
            response.raise_for_status()
            return response
        except requests.RequestException as exc:
            last_error = exc
            if attempt == retries:
                break
            time.sleep(backoff_seconds * attempt)
    assert last_error is not None
    raise last_error


def _record(name: str, endpoint: str, ok: bool, status: str, **kwargs: Any) -> ResourceRecord:
    return ResourceRecord(
        name=name,
        endpoint=endpoint,
        ok=ok,
        status=status,
        checked_at_utc=utc_now_iso(),
        **kwargs,
    )


def _graphql(
    session: requests.Session,
    endpoint: str,
    query: str,
    variables: dict[str, Any],
) -> dict[str, Any]:
    response = _request(session, "POST", endpoint, json={"query": query, "variables": variables})
    payload = response.json()
    if payload.get("errors"):
        raise RuntimeError(str(payload["errors"]))
    return payload["data"]


def _ot_label_score(row: dict[str, Any], allowlist: set[str]) -> float | None:
    scores = row.get("datatypeScores") or []
    allowed = [item["score"] for item in scores if item.get("id") in allowlist]
    if not allowed:
        return None
    return max(allowed)


def count_filtered_targets(rows: list[dict[str, Any]], allowlist: set[str], threshold: float) -> int:
    count = 0
    for row in rows:
        label_score = _ot_label_score(row, allowlist)
        if label_score is not None and label_score >= threshold:
            count += 1
    return count


def fetch_ot_associated_targets(
    session: requests.Session,
    endpoint: str,
    disease_id: str,
    page_size: int = 500,
) -> tuple[int, list[dict[str, Any]]]:
    query = """
    query AssociatedTargets($diseaseId: String!, $index: Int!, $size: Int!) {
      disease(efoId: $diseaseId) {
        associatedTargets(page: {index: $index, size: $size}) {
          count
          rows {
            score
            datatypeScores { id score }
            datasourceScores { id score }
            target { id approvedSymbol }
          }
        }
      }
    }
    """
    rows: list[dict[str, Any]] = []
    index = 0
    total_count = 0
    while True:
        data = _graphql(
            session,
            endpoint,
            query,
            {"diseaseId": disease_id, "index": index, "size": page_size},
        )
        page = data["disease"]["associatedTargets"]
        total_count = int(page["count"])
        page_rows = page["rows"]
        rows.extend(page_rows)
        if len(rows) >= total_count or not page_rows:
            break
        index += 1
    return total_count, rows


def verify_open_targets(config: dict[str, Any], session: requests.Session) -> ResourceRecord:
    endpoint = config["resources"]["open_targets"]["endpoint"]
    started = time.perf_counter()
    try:
        uc_count, uc_rows = fetch_ot_associated_targets(session, endpoint, config["plan"]["benchmark_disease_id"])
        ic_count, ic_rows = fetch_ot_associated_targets(session, endpoint, config["plan"]["target_disease_id"])

        threshold = config["plan"]["label_score_threshold"]
        low_threshold, high_threshold = config["plan"]["label_score_sensitivities"]

        uc_primary = count_filtered_targets(uc_rows, {"genetic_association"}, threshold)
        uc_low = count_filtered_targets(uc_rows, {"genetic_association"}, low_threshold)
        uc_high = count_filtered_targets(uc_rows, {"genetic_association"}, high_threshold)
        uc_union_primary = count_filtered_targets(uc_rows, {"genetic_association", "somatic_mutation"}, threshold)

        ic_primary = count_filtered_targets(ic_rows, {"animal_model"}, threshold)
        ic_low = count_filtered_targets(ic_rows, {"animal_model"}, low_threshold)
        ic_high = count_filtered_targets(ic_rows, {"animal_model"}, high_threshold)

        drift: dict[str, Any] = {}
        ok = True
        status = "ok"
        if uc_primary != config["plan"]["frozen_counts"]["uc_genetics_primary"]:
            ok = False
            status = "drift"
            drift["uc_genetics_primary"] = {
                "expected": config["plan"]["frozen_counts"]["uc_genetics_primary"],
                "observed": uc_primary,
            }
        if ic_primary != config["plan"]["frozen_counts"]["ic_animal_model_primary"]:
            ok = False
            status = "drift"
            drift["ic_animal_model_primary"] = {
                "expected": config["plan"]["frozen_counts"]["ic_animal_model_primary"],
                "observed": ic_primary,
                "note": "Historical count is still verified even though TARGET animal_model is dropped.",
            }

        return _record(
            "open_targets",
            endpoint,
            ok,
            status,
            latency_s=round(time.perf_counter() - started, 3),
            warning="Drift detected in frozen Open Targets counts." if drift else None,
            details={
                "download_policy": "GraphQL metadata only; no patient matrices downloaded.",
                "uc_total_associations": uc_count,
                "ic_total_associations": ic_count,
                "uc_counts": {
                    "genetics_only": {"thr_0_20": uc_low, "thr_0_30": uc_primary, "thr_0_40": uc_high},
                    "genetics_plus_somatic_thr_0_30": uc_union_primary,
                },
                "ic_counts": {
                    "animal_model_only": {"thr_0_20": ic_low, "thr_0_30": ic_primary, "thr_0_40": ic_high},
                    "channel_status": "dropped_for_target_labels",
                },
                "drift": drift,
            },
        )
    except Exception as exc:  # pragma: no cover - exercised by live runs
        return _record(
            "open_targets",
            endpoint,
            False,
            "error",
            latency_s=round(time.perf_counter() - started, 3),
            error=str(exc),
        )


def verify_geo(config: dict[str, Any], session: requests.Session) -> ResourceRecord:
    endpoint = config["resources"]["geo"]["esearch_endpoint"]
    started = time.perf_counter()
    term = config["resources"]["geo"]["target_search_term"]
    try:
        response = _request(
            session,
            "GET",
            endpoint,
            params={"db": "gds", "term": term, "retmode": "json"},
        ).json()
        result = response["esearchresult"]
        accessions = config["plan"]["target_primary_cohorts"] + [
            config["plan"]["target_inflammation_control"],
            config["plan"]["target_scrna"],
        ]
        return _record(
            "geo_eutils",
            endpoint,
            True,
            "ok",
            latency_s=round(time.perf_counter() - started, 3),
            details={
                "search_term": term,
                "count": int(result["count"]),
                "required_accessions": accessions,
                "download_policy": "E-utilities metadata only; no expression matrices fetched.",
            },
        )
    except Exception as exc:  # pragma: no cover - exercised by live runs
        return _record(
            "geo_eutils",
            endpoint,
            False,
            "error",
            latency_s=round(time.perf_counter() - started, 3),
            error=str(exc),
        )


def verify_omnipath(config: dict[str, Any], session: requests.Session) -> ResourceRecord:
    endpoint = config["resources"]["omnipath"]["endpoint"]
    started = time.perf_counter()
    try:
        payload = _request(session, "GET", endpoint, params=config["resources"]["omnipath"]["params"]).json()
        sample_fields: list[str] = []
        first_row = payload[0] if isinstance(payload, list) and payload else payload
        if isinstance(first_row, dict):
            sample_fields = sorted(first_row.keys())
        elif isinstance(first_row, list) and first_row and isinstance(first_row[0], dict):
            sample_fields = sorted(first_row[0].keys())
        return _record(
            "omnipath",
            endpoint,
            bool(payload),
            "ok",
            latency_s=round(time.perf_counter() - started, 3),
            details={
                "sample_fields": sample_fields,
                "row_count": len(payload) if isinstance(payload, list) else 1,
                "payload_type": type(payload).__name__,
            },
        )
    except Exception as exc:  # pragma: no cover - exercised by live runs
        return _record("omnipath", endpoint, False, "error", error=str(exc), latency_s=round(time.perf_counter() - started, 3))


def verify_signor(config: dict[str, Any], session: requests.Session) -> ResourceRecord:
    endpoint = config["resources"]["signor"]["endpoint"]
    started = time.perf_counter()
    try:
        response = _request(session, "GET", endpoint)
        return _record(
            "signor",
            endpoint,
            True,
            "not_verified",
            latency_s=round(time.perf_counter() - started, 3),
            warning="Portal reachability confirmed only; no pinned SIGNOR data-file endpoint is verified yet.",
            details={"http_status": response.status_code, "verification_scope": "download_portal_only"},
        )
    except Exception as exc:  # pragma: no cover - exercised by live runs
        return _record("signor", endpoint, False, "error", error=str(exc), latency_s=round(time.perf_counter() - started, 3))


def verify_string(config: dict[str, Any], session: requests.Session) -> ResourceRecord:
    endpoint = config["resources"]["string"]["endpoint"]
    started = time.perf_counter()
    try:
        payload = _request(session, "GET", endpoint).json()
        version = payload[0].get("string_version") if payload else None
        return _record(
            "string",
            endpoint,
            bool(version),
            "ok" if version else "warning",
            latency_s=round(time.perf_counter() - started, 3),
            details={"string_version": version},
        )
    except Exception as exc:  # pragma: no cover - exercised by live runs
        return _record("string", endpoint, False, "error", error=str(exc), latency_s=round(time.perf_counter() - started, 3))


def verify_dgidb(config: dict[str, Any], session: requests.Session) -> ResourceRecord:
    endpoint = config["resources"]["dgidb"]["endpoint"]
    started = time.perf_counter()
    try:
        payload = _request(session, "POST", endpoint, json={"query": "query { __typename }"}).json()
        return _record(
            "dgidb_graphql",
            endpoint,
            "data" in payload,
            "ok" if "data" in payload else "warning",
            latency_s=round(time.perf_counter() - started, 3),
            details={"keys": sorted(payload.keys())},
        )
    except Exception as exc:  # pragma: no cover - exercised by live runs
        return _record("dgidb_graphql", endpoint, False, "error", error=str(exc), latency_s=round(time.perf_counter() - started, 3))


def verify_chembl(config: dict[str, Any], session: requests.Session) -> ResourceRecord:
    endpoint = config["resources"]["chembl"]["endpoint"]
    started = time.perf_counter()
    try:
        payload = _request(session, "GET", endpoint).json()
        return _record(
            "chembl",
            endpoint,
            True,
            "ok",
            latency_s=round(time.perf_counter() - started, 3),
            details={"keys": sorted(payload.keys())},
        )
    except Exception as exc:  # pragma: no cover - exercised by live runs
        return _record("chembl", endpoint, False, "error", error=str(exc), latency_s=round(time.perf_counter() - started, 3))


def verify_gtex_hpa(config: dict[str, Any], session: requests.Session) -> list[ResourceRecord]:
    records: list[ResourceRecord] = []
    for name, endpoint in config["resources"]["best_effort"].items():
        started = time.perf_counter()
        try:
            response = session.request("HEAD", endpoint, timeout=DEFAULT_TIMEOUT, allow_redirects=True)
            if response.status_code >= 400:
                response = _request(session, "GET", endpoint)
            records.append(
                _record(
                    name,
                    endpoint,
                    True,
                    "not_verified",
                    latency_s=round(time.perf_counter() - started, 3),
                    warning="Site reachability confirmed only; no pinned GTEx/HPA data endpoint is verified yet.",
                    details={"http_status": response.status_code, "verification_scope": "site_only"},
                )
            )
        except Exception as exc:  # pragma: no cover - exercised by live runs
            records.append(
                _record(
                    name,
                    endpoint,
                    False,
                    "error",
                    latency_s=round(time.perf_counter() - started, 3),
                    error=str(exc),
                    details={"verification_scope": "site_only"},
                )
            )
    return records


def verify_clinicaltrials(config: dict[str, Any], session: requests.Session) -> ResourceRecord:
    endpoint = config["resources"]["clinicaltrials"]["endpoint"]
    started = time.perf_counter()
    try:
        payload = _request(
            session,
            "GET",
            endpoint,
            params={"query.cond": "Interstitial Cystitis", "pageSize": 1, "countTotal": "true"},
        ).json()
        return _record(
            "clinicaltrials_v2",
            endpoint,
            "totalCount" in payload,
            "ok" if "totalCount" in payload else "warning",
            latency_s=round(time.perf_counter() - started, 3),
            details={"totalCount": payload.get("totalCount")},
        )
    except Exception as exc:  # pragma: no cover - exercised by live runs
        return _record("clinicaltrials_v2", endpoint, False, "error", error=str(exc), latency_s=round(time.perf_counter() - started, 3))


def verify_decoupler_installation(config: dict[str, Any]) -> ResourceRecord:
    del config
    endpoint = "python:decoupler"
    started = time.perf_counter()
    try:
        import decoupler  # type: ignore

        has_show_resources = hasattr(decoupler, "show_resources")
        return _record(
            "decoupler",
            endpoint,
            True,
            "ok",
            latency_s=round(time.perf_counter() - started, 3),
            details={
                "module": getattr(decoupler, "__file__", "unknown"),
                "has_show_resources": has_show_resources,
                "download_policy": "Import/liveness only; no resource downloads performed.",
            },
        )
    except Exception as exc:
        return _record("decoupler", endpoint, False, "error", error=str(exc), latency_s=round(time.perf_counter() - started, 3))


def verify_resources(config: dict[str, Any]) -> list[dict[str, Any]]:
    session = requests.Session()
    session.headers.update({"User-Agent": "signflow-phase1/0.1.0"})

    records = [
        verify_open_targets(config, session),
        verify_geo(config, session),
        verify_omnipath(config, session),
        verify_signor(config, session),
        verify_string(config, session),
        verify_dgidb(config, session),
        verify_chembl(config, session),
        verify_clinicaltrials(config, session),
        verify_decoupler_installation(config),
    ]
    records.extend(verify_gtex_hpa(config, session))
    return [record.to_dict() for record in records]
