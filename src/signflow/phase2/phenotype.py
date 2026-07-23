"""Phenotype parsing and frozen-scheme sample selection for Phase 2."""

from __future__ import annotations

import csv
import gzip
import re
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from . import PHENOTYPE_MAP_PATH, StopAndReportError


def _read_text(path: str | Path) -> str:
    payload = Path(path).read_bytes()
    if str(path).endswith(".gz"):
        return gzip.decompress(payload).decode("utf-8", errors="replace")
    return payload.decode("utf-8", errors="replace")


def _parse_tsv_line(line: str) -> list[str]:
    return next(csv.reader([line], delimiter="\t", quotechar='"'))


def load_phenotype_map(path: str | Path = PHENOTYPE_MAP_PATH) -> dict[str, Any]:
    return yaml.safe_load(Path(path).read_text())


def parse_series_matrix_samples(series_matrix_path: str | Path) -> pd.DataFrame:
    """Parse GEO sample metadata from a series matrix file."""

    sample_fields: dict[str, list[list[str]]] = {}
    sample_ids: list[str] | None = None
    sample_titles: list[str] | None = None

    for line in _read_text(series_matrix_path).splitlines():
        if line.startswith("!series_matrix_table_begin"):
            break
        if not line.startswith("!Sample_"):
            continue

        values = _parse_tsv_line(line)
        key = values[0][1:]
        field_values = values[1:]
        if key == "Sample_geo_accession":
            sample_ids = field_values
        if key == "Sample_title":
            sample_titles = field_values
        sample_fields.setdefault(key, []).append(field_values)

    if not sample_ids:
        raise StopAndReportError("Series matrix is missing !Sample_geo_accession metadata", reason="unavailable-files")

    rows = []
    for index, sample_id in enumerate(sample_ids):
        row: dict[str, Any] = {"sample_id": sample_id}
        if sample_titles:
            row["sample_title"] = sample_titles[index]
        for key, value_lists in sample_fields.items():
            row[key] = [values[index] for values in value_lists if index < len(values) and values[index]]
        rows.append(row)
    return pd.DataFrame(rows)


def _normalized_strings(record: pd.Series) -> list[str]:
    texts: list[str] = []
    for key, value in record.items():
        if key == "sample_id" or value is None:
            continue
        if isinstance(value, list):
            texts.extend(str(item).strip() for item in value if str(item).strip())
        elif str(value).strip():
            texts.append(str(value).strip())
    return texts


def _extract_title_metadata(accession: str, sample_title: str, study_map: dict[str, Any]) -> dict[str, str]:
    regex = study_map.get("title_regex")
    if not regex:
        return {}
    match = re.search(regex, sample_title)
    if not match:
        raise StopAndReportError(
            f"{accession} title did not match the frozen parsing rule: {sample_title}",
            accession=accession,
            reason="ambiguous-phenotype",
            details={"sample_title": sample_title},
        )
    extracted = {key: value for key, value in match.groupdict().items() if value}
    site_aliases = study_map.get("site_aliases", {})
    if "site_code" in extracted:
        extracted["site"] = site_aliases.get(extracted["site_code"], extracted["site_code"])
    return extracted


def _match_rule(normalized: list[str], rule: dict) -> bool:
    all_of = [s.lower() for s in rule.get("all_of", [])]
    any_of = [s.lower() for s in rule.get("any_of", [])]
    if not all_of and not any_of:
        return False
    if all_of and not all(any(req in text for text in normalized) for req in all_of):
        return False
    if any_of and not any(any(req in text for text in normalized) for req in any_of):
        return False
    return True


def classify_phenotypes(
    accession: str,
    series_matrix_path: str | Path,
    mapping_path: str | Path = PHENOTYPE_MAP_PATH,
) -> pd.DataFrame:
    """Classify each sample as case, control, or excluded."""

    mappings = load_phenotype_map(mapping_path)
    if accession not in mappings:
        raise StopAndReportError(
            f"No frozen phenotype map exists for {accession}",
            accession=accession,
            reason="ambiguous-phenotype",
        )

    study_map = mappings[accession]
    samples = parse_series_matrix_samples(series_matrix_path)
    rows = []
    for _, record in samples.iterrows():
        raw_texts = _normalized_strings(record)
        normalized = [text.lower() for text in raw_texts]
        matches = []

        if "rules" in study_map:
            for rule in study_map["rules"]:
                if _match_rule(normalized, rule):
                    label = "+".join(rule.get("all_of") or rule.get("any_of") or ["rule"])
                    matches.append((label, rule))
        else:
            raw_map = study_map["raw_string_map"]
            for raw_key, payload in raw_map.items():
                raw_key_normalized = raw_key.lower()
                if any(raw_key_normalized in text for text in normalized):
                    matches.append((raw_key, payload))

        groups = {payload["group"] for _, payload in matches}
        if len(groups) != 1 or not matches:
            raise StopAndReportError(
                f"{accession} sample {record['sample_id']} could not be classified unambiguously",
                accession=accession,
                reason="ambiguous-phenotype",
                details={"sample_id": record["sample_id"], "matches": [match[0] for match in matches]},
            )

        metadata = _extract_title_metadata(accession, str(record.get("sample_title", "")), study_map)
        matched_key, matched_payload = matches[0]
        metadata.setdefault("site", matched_payload.get("site"))
        rows.append(
            {
                "sample_id": record["sample_id"],
                "sample_title": record.get("sample_title", ""),
                "group": matched_payload["group"],
                "patient_id": metadata.get("patient_id"),
                "site": metadata.get("site"),
                "raw_match": matched_key,
                "raw_blob": " | ".join(raw_texts),
            }
        )

    return pd.DataFrame(rows).sort_values(["group", "sample_id"]).reset_index(drop=True)


def phenotype_mapping_summary(phenotypes: pd.DataFrame) -> str:
    """Return a text summary of mapping decisions and counts."""

    counts = phenotypes["group"].value_counts().to_dict()
    lines = ["Phenotype mapping summary:"]
    lines.extend(f"- {group}: {count}" for group, count in sorted(counts.items()))
    for raw_value, count in Counter(phenotypes["raw_match"]).items():
        lines.append(f"- {raw_value} -> {count}")
    return "\n".join(lines)


def select_uc_only(phenotypes: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    """Keep UC cases and controls while counting excluded samples."""

    counts = phenotypes["group"].value_counts().to_dict()
    filtered = phenotypes[phenotypes["group"] != "excluded"].reset_index(drop=True)
    return filtered, {"excluded": int(counts.get("excluded", 0))}


def select_gse11783_primary_samples(
    phenotypes: pd.DataFrame,
    primary_site: str = "non-ulcer",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Enforce one sample per patient with the frozen primary site rule."""

    controls = phenotypes[phenotypes["group"] == "control"]
    cases = phenotypes[phenotypes["group"] == "case"]
    selected_cases = []
    mapping_rows = []
    for patient_id, patient_rows in cases.groupby("patient_id", dropna=False):
        if not patient_id:
            raise StopAndReportError(
                "GSE11783 case sample is missing a parsed patient identifier",
                accession="GSE11783",
                reason="ambiguous-phenotype",
            )
        matches = patient_rows[patient_rows["site"] == primary_site]
        if len(matches) != 1:
            raise StopAndReportError(
                f"GSE11783 patient {patient_id} did not resolve to exactly one {primary_site} sample",
                accession="GSE11783",
                reason="ambiguous-phenotype",
                details={"patient_id": patient_id, "available_sites": sorted(patient_rows['site'].dropna().unique())},
            )
        chosen = matches.iloc[0]
        selected_cases.append(chosen)
        mapping_rows.append({"patient_id": patient_id, "selected_sample_id": chosen["sample_id"], "selected_site": primary_site})

    selected = pd.concat([controls, pd.DataFrame(selected_cases)], ignore_index=True)
    selected = selected.sort_values(["group", "sample_id"]).reset_index(drop=True)
    mapping = pd.DataFrame(mapping_rows).sort_values("patient_id").reset_index(drop=True)
    return selected, mapping
