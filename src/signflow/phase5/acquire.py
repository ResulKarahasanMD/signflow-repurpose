"""Acquire drug-target resources with SHA-256 cached manifests."""

from __future__ import annotations

import json
import os
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd
import requests
import yaml

from signflow.config import load_config
from signflow.phase2 import PROJECT_ROOT

from . import BENCHMARK_DRUGS, FROZEN_PLAN, PHASE5_RAW, StopAndReportError, ensure_directory, sha256_path

REQUEST_TIMEOUT = 120
CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
DGIDB_ENDPOINT = "https://dgidb.org/api/graphql"
OT_ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"


def _utc_now() -> str:
    return datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def _request_json(url: str, *, method: str = "GET", **kwargs: Any) -> Any:
    last: Exception | None = None
    for attempt in range(1, 4):
        try:
            response = requests.request(method, url, timeout=REQUEST_TIMEOUT, **kwargs)
            response.raise_for_status()
            return response.json()
        except requests.RequestException as exc:
            last = exc
            time.sleep(attempt)
    assert last is not None
    raise StopAndReportError(f"HTTP request failed: {url}", reason="network-resource-unavailable", details={"error": str(last)})


def _graphql(endpoint: str, query: str, variables: dict[str, Any] | None = None) -> dict[str, Any]:
    payload = {"query": query, "variables": variables or {}}
    data = _request_json(endpoint, method="POST", json=payload)
    if data.get("errors"):
        raise StopAndReportError(str(data["errors"]), reason="graphql-error", details={"endpoint": endpoint})
    return data["data"]


def verify_benchmark_drugs_sha() -> dict[str, object]:
    plan = yaml.safe_load(FROZEN_PLAN.read_text())["plan"]
    expected = str(plan["benchmark_gold_drugs_sha256"])
    observed = sha256_path(BENCHMARK_DRUGS)
    if observed != expected:
        raise StopAndReportError(
            f"benchmark_drugs.yaml SHA mismatch: {observed} != {expected}",
            reason="gold-drugs-sha-mismatch",
            details={"expected": expected, "observed": observed, "path": str(BENCHMARK_DRUGS)},
        )
    return {"path": str(BENCHMARK_DRUGS), "sha256": observed, "verified": True}


def load_gold_drugs() -> pd.DataFrame:
    verify_benchmark_drugs_sha()
    doc = yaml.safe_load(BENCHMARK_DRUGS.read_text())
    rows = []
    for item in doc.get("drugs", []):
        if not item.get("included", True):
            continue
        rows.append(
            {
                "chembl_id": item["chembl_id"],
                "pref_name": item.get("pref_name", ""),
                "canonical_group": item.get("canonical_group", ""),
                "is_gold": True,
            }
        )
    return pd.DataFrame(rows)


def _cache_manifest(raw_root: Path, name: str, record: dict[str, Any]) -> None:
    manifest_dir = ensure_directory(raw_root / "manifests")
    (manifest_dir / f"{name}.json").write_text(json.dumps(record, indent=2) + "\n")


def fetch_opentargets_candidates(disease_id: str, raw_root: Path) -> dict[str, object]:
    cache = raw_root / "opentargets" / f"{disease_id}_drug_candidates.json"
    ensure_directory(cache.parent)
    query = """
    query Cand($efoId: String!) {
      disease(efoId: $efoId) {
        id
        name
        drugAndClinicalCandidates {
          count
          rows {
            maxClinicalStage
            drug { id name }
          }
        }
      }
    }
    """
    if cache.exists() and cache.stat().st_size > 0:
        record = {
            "source": "opentargets",
            "disease_id": disease_id,
            "endpoint": OT_ENDPOINT,
            "path": str(cache),
            "sha256": sha256_path(cache),
            "status": "reused",
            "accessed_at_utc": _utc_now(),
            "note": "Open Targets v4 exposes drugAndClinicalCandidates (knownDrugs unavailable).",
        }
        _cache_manifest(raw_root, f"opentargets_{disease_id}", record)
        return record

    data = _graphql(OT_ENDPOINT, query, {"efoId": disease_id})
    payload = {
        "accessed_at_utc": _utc_now(),
        "endpoint": OT_ENDPOINT,
        "disease_id": disease_id,
        "response": data,
        "note": "knownDrugs field absent in OT v4; using drugAndClinicalCandidates.",
    }
    cache.write_text(json.dumps(payload, indent=2) + "\n")
    record = {
        "source": "opentargets",
        "disease_id": disease_id,
        "endpoint": OT_ENDPOINT,
        "path": str(cache),
        "sha256": sha256_path(cache),
        "status": "downloaded",
        "accessed_at_utc": payload["accessed_at_utc"],
    }
    _cache_manifest(raw_root, f"opentargets_{disease_id}", record)
    return record



def ot_all_candidates(disease_id: str, raw_root: Path) -> pd.DataFrame:
    """All Open Targets drugAndClinicalCandidates for the disease (ranking background).

    Gold rediscovery uses the frozen APPROVAL set; the ranking universe must be broader
    than gold alone or enrichment@k is structurally uninformative.
    """

    cache = raw_root / "opentargets" / f"{disease_id}_drug_candidates.json"
    if not cache.exists():
        fetch_opentargets_candidates(disease_id, raw_root)
    payload = json.loads(cache.read_text())
    rows = payload["response"]["disease"]["drugAndClinicalCandidates"]["rows"]
    out = []
    for row in rows:
        drug = row["drug"]
        out.append(
            {
                "chembl_id": drug["id"],
                "pref_name": drug.get("name", ""),
                "max_clinical_stage": row.get("maxClinicalStage"),
                "is_gold": False,
            }
        )
    return pd.DataFrame(out).drop_duplicates(subset=["chembl_id"])

def ot_approval_drugs(disease_id: str, raw_root: Path) -> pd.DataFrame:
    cache = raw_root / "opentargets" / f"{disease_id}_drug_candidates.json"
    if not cache.exists():
        fetch_opentargets_candidates(disease_id, raw_root)
    payload = json.loads(cache.read_text())
    rows = payload["response"]["disease"]["drugAndClinicalCandidates"]["rows"]
    out = []
    for row in rows:
        if row.get("maxClinicalStage") != "APPROVAL":
            continue
        drug = row["drug"]
        out.append({"chembl_id": drug["id"], "pref_name": drug.get("name", ""), "is_gold": False})
    return pd.DataFrame(out).drop_duplicates(subset=["chembl_id"])


def _chembl_get_json(session: requests.Session, url: str, *, retries: int = 5) -> dict:
    last_err: Exception | None = None
    for attempt in range(retries):
        try:
            resp = session.get(url, timeout=REQUEST_TIMEOUT)
            if resp.status_code in {429, 500, 502, 503, 504} or not resp.text.strip():
                time.sleep(min(2 ** attempt, 30))
                continue
            resp.raise_for_status()
            return resp.json()
        except Exception as exc:  # noqa: BLE001
            last_err = exc
            time.sleep(min(2 ** attempt, 30))
    raise RuntimeError(f"ChEMBL request failed after {retries} retries: {url}") from last_err


def _chembl_target_genes(target_chembl_id: str, session: requests.Session, raw_root: Path) -> list[str]:
    if not target_chembl_id or target_chembl_id == "None":
        return []
    cache = raw_root / "chembl" / "targets" / f"{target_chembl_id}.json"
    ensure_directory(cache.parent)
    if cache.exists() and cache.stat().st_size > 0:
        payload = json.loads(cache.read_text())
    else:
        url = f"{CHEMBL_BASE}/target/{target_chembl_id}.json"
        payload = _chembl_get_json(session, url)
        cache.write_text(json.dumps(payload, indent=2) + "\n")
    genes: list[str] = []
    for comp in payload.get("target_components", []) or []:
        for syn in comp.get("target_component_synonyms", []) or []:
            if syn.get("syn_type") == "GENE_SYMBOL":
                genes.append(str(syn["component_synonym"]))
    return sorted(set(genes))


def fetch_chembl_mechanisms(chembl_ids: list[str], raw_root: Path) -> dict[str, object]:
    ensure_directory(raw_root / "chembl")
    session = requests.Session()
    session.headers.update({"User-Agent": "signflow-phase5/0.1.0"})
    all_rows: list[dict[str, object]] = []
    n_ids = len(set(chembl_ids))
    for i, chembl_id in enumerate(sorted(set(chembl_ids)), start=1):
        if i == 1 or i % 25 == 0 or i == n_ids:
            print(f"[phase5] chembl mechanisms {i}/{n_ids}", flush=True)
        cache = raw_root / "chembl" / f"{chembl_id}_mechanisms.json"
        if cache.exists() and cache.stat().st_size > 0:
            mechs = json.loads(cache.read_text()).get("mechanisms", [])
        else:
            url = f"{CHEMBL_BASE}/mechanism.json?molecule_chembl_id={chembl_id}&limit=1000"
            payload = _chembl_get_json(session, url)
            mechs = payload.get("mechanisms", [])
            cache.write_text(json.dumps({"mechanisms": mechs}, indent=2) + "\n")
            time.sleep(0.05)
        for mech in mechs:
            tid = mech.get("target_chembl_id")
            if not tid:
                continue
            genes = _chembl_target_genes(str(tid), session, raw_root)
            for gene in genes:
                all_rows.append(
                    {
                        "chembl_id": chembl_id,
                        "gene_symbol": gene,
                        "raw_action": mech.get("action_type") or mech.get("mechanism_of_action"),
                        "source": "chembl_mechanism",
                        "target_chembl_id": tid,
                    }
                )
    frame = pd.DataFrame(all_rows)
    out_path = raw_root / "chembl" / "mechanisms_merged.tsv"
    frame.to_csv(out_path, sep="\t", index=False)
    record = {
        "source": "chembl",
        "endpoint": CHEMBL_BASE,
        "path": str(out_path),
        "sha256": sha256_path(out_path),
        "n_rows": int(len(frame)),
        "n_drugs": int(frame["chembl_id"].nunique()) if len(frame) else 0,
        "accessed_at_utc": _utc_now(),
    }
    _cache_manifest(raw_root, "chembl_mechanisms", record)
    return record


def _dgidb_chunk(names: list[str]) -> list[dict[str, object]]:
    query = """
    query DrugIx($names: [String!]!) {
      drugs(names: $names) {
        nodes {
          name
          conceptId
          interactions {
            gene { name }
            interactionTypes { type directionality }
          }
        }
      }
    }
    """
    data = _graphql(DGIDB_ENDPOINT, query, {"names": names})
    rows: list[dict[str, object]] = []
    for node in data.get("drugs", {}).get("nodes", []) or []:
        drug_name = node.get("name")
        for ix in node.get("interactions", []) or []:
            gene = (ix.get("gene") or {}).get("name")
            if not gene:
                continue
            raw_action = None
            types = ix.get("interactionTypes") or []
            if types:
                raw_action = types[0].get("type") or types[0].get("directionality")
            rows.append(
                {
                    "drug_name": drug_name,
                    "gene_symbol": gene,
                    "raw_action": raw_action,
                    "source": "dgidb",
                }
            )
    return rows


def fetch_dgidb_by_names(drug_names: list[str], raw_root: Path) -> dict[str, object]:
    ensure_directory(raw_root / "dgidb")
    names = sorted({n.strip().upper() for n in drug_names if n and str(n).strip()})
    rows: list[dict[str, object]] = []
    chunk = 25
    for i in range(0, len(names), chunk):
        rows.extend(_dgidb_chunk(names[i : i + chunk]))
    frame = pd.DataFrame(rows)
    out_path = raw_root / "dgidb" / "interactions_by_name.tsv"
    frame.to_csv(out_path, sep="\t", index=False)
    record = {
        "source": "dgidb",
        "endpoint": DGIDB_ENDPOINT,
        "path": str(out_path),
        "sha256": sha256_path(out_path),
        "n_rows": int(len(frame)),
        "n_drugs_queried": len(names),
        "accessed_at_utc": _utc_now(),
    }
    _cache_manifest(raw_root, "dgidb_interactions", record)
    return record


def maybe_fetch_drugbank(raw_root: Path) -> dict[str, object] | None:
    xml_path = os.environ.get("DRUGBANK_XML")
    if not xml_path:
        return None
    path = Path(xml_path)
    if not path.exists():
        raise StopAndReportError(f"DRUGBANK_XML set but missing: {path}", reason="drugbank-missing")
    record = {
        "source": "drugbank",
        "path": str(path),
        "sha256": sha256_path(path),
        "status": "present_not_parsed",
        "accessed_at_utc": _utc_now(),
        "note": "Parsing deferred; ChEMBL+DGIdb primary.",
    }
    _cache_manifest(raw_root, "drugbank", record)
    return record


def fetch_chembl_approved_molecules(raw_root: Path, *, max_phase: int = 4) -> dict[str, object]:
    """Disease-agnostic approved-drug universe from ChEMBL (max_phase filter; no disease ID)."""

    pool_dir = ensure_directory(raw_root / "pool")
    cache = pool_dir / "chembl_max_phase4_molecules.json"
    out_tsv = pool_dir / "approved_drug_universe.tsv"
    session = requests.Session()
    session.headers.update({"User-Agent": "signflow-phase5/0.1.0"})

    if cache.exists() and cache.stat().st_size > 0 and out_tsv.exists():
        molecules = json.loads(cache.read_text()).get("molecules", [])
    else:
        molecules = []
        limit = 1000
        offset = 0
        total = None
        while True:
            url = f"{CHEMBL_BASE}/molecule.json?max_phase={max_phase}&limit={limit}&offset={offset}"
            print(f"[phase5] chembl approved molecules offset={offset}", flush=True)
            payload = _chembl_get_json(session, url)
            batch = payload.get("molecules", []) or []
            molecules.extend(batch)
            meta = payload.get("page_meta") or {}
            total = meta.get("total_count", total)
            offset += limit
            if not batch or (total is not None and offset >= int(total)):
                break
            time.sleep(0.05)
        cache.write_text(
            json.dumps(
                {
                    "accessed_at_utc": _utc_now(),
                    "endpoint": CHEMBL_BASE,
                    "query": f"molecule.json?max_phase={max_phase}",
                    "disease_conditioning": None,
                    "n_molecules": len(molecules),
                    "molecules": molecules,
                },
                indent=2,
            )
            + "\n"
        )

    rows = []
    for mol in molecules:
        chembl_id = mol.get("molecule_chembl_id")
        if not chembl_id:
            continue
        pref = mol.get("pref_name") or ""
        rows.append(
            {
                "chembl_id": str(chembl_id),
                "pref_name": str(pref),
                "max_phase": mol.get("max_phase", max_phase),
                "source": "chembl_max_phase_4",
                "is_gold": False,
            }
        )
    frame = pd.DataFrame(rows).drop_duplicates(subset=["chembl_id"])
    # Union frozen gold (labels only; must remain a small minority)
    gold = load_gold_drugs()
    gold = gold.copy()
    gold["max_phase"] = 4
    gold["source"] = "frozen_uc_gold"
    gold_ids = set(gold["chembl_id"].astype(str))
    frame = frame.loc[~frame["chembl_id"].astype(str).isin(gold_ids)]
    universe = pd.concat([gold, frame], ignore_index=True).drop_duplicates(subset=["chembl_id"], keep="first")
    if "canonical_group" not in universe.columns:
        universe["canonical_group"] = ""
    missing_group = universe["canonical_group"].isna() | (universe["canonical_group"].astype(str).str.strip() == "")
    universe.loc[missing_group, "canonical_group"] = (
        universe.loc[missing_group, "pref_name"].astype(str).str.lower().str.replace(" ", "_")
    )
    universe.to_csv(out_tsv, sep="\t", index=False)
    digest = sha256_path(out_tsv)
    record = {
        "source": "chembl_max_phase_4",
        "endpoint": CHEMBL_BASE,
        "query": f"molecule.json?max_phase={max_phase}",
        "disease_conditioning": None,
        "path": str(out_tsv),
        "cache_path": str(cache),
        "sha256": digest,
        "n_molecules_raw": int(len(molecules)),
        "n_universe": int(len(universe)),
        "n_gold": int(universe["is_gold"].astype(bool).sum()) if "is_gold" in universe.columns else 0,
        "accessed_at_utc": _utc_now(),
    }
    _cache_manifest(raw_root, "approved_drug_universe", record)
    return record


def load_approved_drug_universe(raw_root: Path | None = None) -> pd.DataFrame:
    raw_root = Path(raw_root or PHASE5_RAW)
    path = raw_root / "pool" / "approved_drug_universe.tsv"
    if not path.exists():
        fetch_chembl_approved_molecules(raw_root)
    return pd.read_csv(path, sep="\t")


def assert_pool_guards(
    universe: pd.DataFrame,
    *,
    n_in_substrate: int | None = None,
    plan: dict[str, object] | None = None,
) -> dict[str, object]:
    """Stop-and-report if pool floor/ceiling violated (failure mode that voided the 30-drug run)."""

    plan = plan or yaml.safe_load(FROZEN_PLAN.read_text())["plan"]
    n_total = int(len(universe))
    n_gold = int(universe["is_gold"].astype(bool).sum()) if "is_gold" in universe.columns else 0
    prevalence = float(n_gold / max(n_total, 1))
    floor = int(plan["phase5_candidate_pool_min_in_substrate"])
    ceiling = float(plan["phase5_candidate_pool_gold_prevalence_ceiling"])
    report = {
        "n_universe": n_total,
        "n_gold": n_gold,
        "gold_prevalence": prevalence,
        "n_in_substrate": n_in_substrate,
        "floor_in_substrate": floor,
        "gold_prevalence_ceiling": ceiling,
    }
    if prevalence > ceiling:
        raise StopAndReportError(
            f"Gold prevalence {prevalence:.3f} exceeds ceiling {ceiling:.3f} "
            f"(n_gold={n_gold}, n_universe={n_total}) — candidate-pool failure mode.",
            reason="gold-prevalence-ceiling-exceeded",
            details=report,
        )
    if n_in_substrate is not None and n_in_substrate < floor:
        raise StopAndReportError(
            f"In-substrate pool size {n_in_substrate} below floor {floor}.",
            reason="candidate-pool-below-floor",
            details=report,
        )
    return report


def fetch_chembl_mechanisms_for_universe(chembl_ids: list[str], raw_root: Path) -> dict[str, object]:
    """Prefer bulk mechanism pagination filtered to the universe; fall back to per-drug cache."""

    ensure_directory(raw_root / "chembl")
    session = requests.Session()
    session.headers.update({"User-Agent": "signflow-phase5/0.1.0"})
    wanted = set(map(str, chembl_ids))
    bulk_cache = raw_root / "chembl" / "mechanisms_all_pages.json"
    all_mechs: list[dict[str, object]] = []
    if bulk_cache.exists() and bulk_cache.stat().st_size > 0:
        all_mechs = json.loads(bulk_cache.read_text()).get("mechanisms", [])
    else:
        limit = 1000
        offset = 0
        total = None
        while True:
            url = f"{CHEMBL_BASE}/mechanism.json?limit={limit}&offset={offset}"
            print(f"[phase5] chembl mechanisms bulk offset={offset}", flush=True)
            payload = _chembl_get_json(session, url)
            batch = payload.get("mechanisms", []) or []
            all_mechs.extend(batch)
            meta = payload.get("page_meta") or {}
            total = meta.get("total_count", total)
            offset += limit
            if not batch or (total is not None and offset >= int(total)):
                break
            time.sleep(0.05)
        bulk_cache.write_text(json.dumps({"mechanisms": all_mechs, "accessed_at_utc": _utc_now()}, indent=2) + "\n")

    filtered = [m for m in all_mechs if str(m.get("molecule_chembl_id")) in wanted]
    # Also fill missing drugs via per-id cache (covers IDs absent from bulk dump)
    have = {str(m.get("molecule_chembl_id")) for m in filtered}
    missing = sorted(wanted - have)
    if missing:
        print(f"[phase5] chembl per-drug fill for {len(missing)} ids absent from bulk", flush=True)
        for i, chembl_id in enumerate(missing, start=1):
            if i == 1 or i % 50 == 0 or i == len(missing):
                print(f"[phase5] chembl fill {i}/{len(missing)}", flush=True)
            cache = raw_root / "chembl" / f"{chembl_id}_mechanisms.json"
            if cache.exists() and cache.stat().st_size > 0:
                mechs = json.loads(cache.read_text()).get("mechanisms", [])
            else:
                url = f"{CHEMBL_BASE}/mechanism.json?molecule_chembl_id={chembl_id}&limit=1000"
                payload = _chembl_get_json(session, url)
                mechs = payload.get("mechanisms", [])
                cache.write_text(json.dumps({"mechanisms": mechs}, indent=2) + "\n")
                time.sleep(0.05)
            filtered.extend(mechs)

    all_rows: list[dict[str, object]] = []
    target_ids = sorted({str(m.get("target_chembl_id")) for m in filtered if m.get("target_chembl_id")})
    gene_cache: dict[str, list[str]] = {}
    for i, tid in enumerate(target_ids, start=1):
        if i == 1 or i % 100 == 0 or i == len(target_ids):
            print(f"[phase5] chembl targets {i}/{len(target_ids)}", flush=True)
        gene_cache[tid] = _chembl_target_genes(tid, session, raw_root)
    for mech in filtered:
        chembl_id = str(mech.get("molecule_chembl_id") or "")
        tid = mech.get("target_chembl_id")
        if not chembl_id or not tid:
            continue
        for gene in gene_cache.get(str(tid), []):
            all_rows.append(
                {
                    "chembl_id": chembl_id,
                    "gene_symbol": gene,
                    "raw_action": mech.get("action_type") or mech.get("mechanism_of_action"),
                    "source": "chembl_mechanism",
                    "target_chembl_id": tid,
                }
            )
    frame = pd.DataFrame(all_rows)
    out_path = raw_root / "chembl" / "mechanisms_merged.tsv"
    frame.to_csv(out_path, sep="\t", index=False)
    record = {
        "source": "chembl",
        "endpoint": CHEMBL_BASE,
        "path": str(out_path),
        "sha256": sha256_path(out_path),
        "n_rows": int(len(frame)),
        "n_drugs": int(frame["chembl_id"].nunique()) if len(frame) else 0,
        "n_mechanisms_filtered": int(len(filtered)),
        "accessed_at_utc": _utc_now(),
    }
    _cache_manifest(raw_root, "chembl_mechanisms", record)
    return record


def acquire_drug_resources(drug_universe: pd.DataFrame, *, raw_root: str | Path | None = None) -> dict[str, object]:
    verify_benchmark_drugs_sha()
    raw_root = ensure_directory(raw_root or PHASE5_RAW)
    plan = yaml.safe_load(FROZEN_PLAN.read_text())["plan"]
    config = load_config(PROJECT_ROOT / "configs" / "default.yaml")

    if str(plan.get("phase5_candidate_pool_disease_conditioning")) != "forbidden":
        raise StopAndReportError(
            "phase5_candidate_pool_disease_conditioning must be forbidden",
            reason="policy-drift",
        )

    pool_record = fetch_chembl_approved_molecules(raw_root)
    # Re-load canonical universe (gold ∪ approved) after freeze write
    drug_universe = load_approved_drug_universe(raw_root)
    assert_pool_guards(drug_universe, plan=plan)

    chembl_ids = drug_universe["chembl_id"].astype(str).tolist()
    names = drug_universe.get("pref_name", pd.Series(dtype=str)).astype(str).tolist()

    manifests = {
        "benchmark_drugs_sha": verify_benchmark_drugs_sha(),
        "approved_drug_universe": pool_record,
        # Disease-conditioned OT retained only as annotation/audit, never as ranking universe
        "opentargets_benchmark_annotation_only": fetch_opentargets_candidates(
            str(plan["benchmark_disease_id"]), raw_root
        ),
        "chembl": fetch_chembl_mechanisms_for_universe(chembl_ids, raw_root),
        "dgidb": fetch_dgidb_by_names(names, raw_root),
        "drugbank": maybe_fetch_drugbank(raw_root),
        "config_dgidb_endpoint": config["resources"]["dgidb"]["endpoint"],
        "config_chembl_endpoint": config["resources"]["chembl"]["endpoint"],
        "drug_data_touched": True,
    }
    summary_path = raw_root / "acquisition_summary.json"
    summary_path.write_text(json.dumps(manifests, indent=2, default=str) + "\n")
    return {
        "raw_root": str(raw_root),
        "manifests": manifests,
        "summary_path": str(summary_path),
        "drug_data_touched": True,
        "drug_universe_path": pool_record["path"],
        "drug_universe_sha256": pool_record["sha256"],
    }
