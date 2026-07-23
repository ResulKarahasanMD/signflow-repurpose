"""Build signed directed causal edges and separate STRING undirected layer."""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from . import FROZEN_PLAN, StopAndReportError


def load_conflict_rule() -> str:
    return str(yaml.safe_load(FROZEN_PLAN.read_text())["plan"]["edge_sign_conflict_rule"])


def load_string_score_min() -> int:
    return int(yaml.safe_load(FROZEN_PLAN.read_text())["plan"]["string_score_min"])


def _ref_count(value: object) -> int:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return 0
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none"}:
        return 0
    # OmniPath references like "TRIP:22226147;Signor:123"
    parts = [p for p in text.replace("|", ";").split(";") if p.strip()]
    return len(parts)



def parse_omnipath(path: str | Path, *, conflict_rule: str) -> tuple[pd.DataFrame, dict[str, int]]:
    frame = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    for col in ("is_directed", "is_stimulation", "is_inhibition"):
        if col not in frame.columns:
            raise StopAndReportError(f"OmniPath missing {col}", reason="omnipath-fields-absent")

    def _bool(series: pd.Series) -> pd.Series:
        return series.astype(str).str.lower().isin({"true", "1", "t", "yes"})

    directed = _bool(frame["is_directed"])
    stim = _bool(frame["is_stimulation"])
    inhib = _bool(frame["is_inhibition"])
    frame = frame.loc[directed].copy()
    stim = stim.loc[frame.index]
    inhib = inhib.loc[frame.index]
    conflict = stim & inhib
    n_conflict = int(conflict.sum())

    relation = np.where(
        stim & inhib,
        "ambiguous",
        np.where(stim, "stimulation", np.where(inhib, "inhibition", "unsigned_directed")),
    )
    sign = np.where(relation == "stimulation", 1, np.where(relation == "inhibition", -1, 0))
    if conflict_rule == "drop":
        keep = ~(stim & inhib)
        frame = frame.loc[keep]
        relation = relation[keep.to_numpy()]
        sign = sign[keep.to_numpy()]
        stim = stim.loc[frame.index]
        inhib = inhib.loc[frame.index]

    src = frame["source_genesymbol"].astype(str).str.strip()
    tgt = frame["target_genesymbol"].astype(str).str.strip()
    ok = (src != "") & (tgt != "") & (src != "nan") & (tgt != "nan")
    refs = frame.get("references", pd.Series([""] * len(frame)))
    n_ref = refs.fillna("").astype(str).map(_ref_count)
    edges = pd.DataFrame(
        {
            "source": src[ok].to_numpy(),
            "target": tgt[ok].to_numpy(),
            "sign": sign[ok.to_numpy()],
            "relation_type": relation[ok.to_numpy()],
            "confidence": n_ref[ok].astype(float).to_numpy(),
            "provenance": "omnipath",
            "source_datasets": frame.get("sources", pd.Series([""] * len(frame)))[ok].astype(str).to_numpy(),
            "n_references": n_ref[ok].astype(int).to_numpy(),
            "source_uniprot": frame.get("source", pd.Series([""] * len(frame)))[ok].astype(str).to_numpy(),
            "target_uniprot": frame.get("target", pd.Series([""] * len(frame)))[ok].astype(str).to_numpy(),
        }
    )
    stats = {
        "n_input_directed": int(directed.sum()),
        "n_conflict_stim_inhib": n_conflict,
        "n_edges_out": int(len(edges)),
        "n_conflict_retained_or_dropped": n_conflict,
    }
    return edges, stats


def _signor_sign(effect: str) -> tuple[int, str]:
    e = effect.strip().strip('"').lower()
    if e.startswith("up-regulates"):
        return 1, "stimulation"
    if e.startswith("down-regulates"):
        return -1, "inhibition"
    if e in {"unknown", "form complex"} or "complex" in e:
        return 0, "ambiguous"
    return 0, "ambiguous"


def parse_signor(path: str | Path) -> tuple[pd.DataFrame, dict[str, int]]:
    frame = pd.read_csv(path, sep="\t", dtype=str, low_memory=False, quoting=0)
    # Human only
    tax = frame["TAX_ID"].astype(str).str.strip()
    frame = frame.loc[tax == "9606"].copy()
    # Protein–protein only
    frame = frame.loc[
        (frame["TYPEA"].astype(str).str.lower() == "protein")
        & (frame["TYPEB"].astype(str).str.lower() == "protein")
    ].copy()
    rows = []
    for _, row in frame.iterrows():
        src = str(row.get("ENTITYA") or "").strip()
        tgt = str(row.get("ENTITYB") or "").strip()
        if not src or not tgt:
            continue
        sign, relation = _signor_sign(str(row.get("EFFECT") or ""))
        pmid = str(row.get("PMID") or "").strip()
        n_ref = len([p for p in pmid.replace("|", ";").split(";") if p.strip()]) if pmid else 0
        rows.append(
            {
                "source": src,
                "target": tgt,
                "sign": sign,
                "relation_type": relation,
                "confidence": float(n_ref),
                "provenance": "signor",
                "source_datasets": "SIGNOR",
                "n_references": n_ref,
                "source_uniprot": str(row.get("IDA") or ""),
                "target_uniprot": str(row.get("IDB") or ""),
            }
        )
    edges = pd.DataFrame(rows)
    return edges, {"n_human_protein_edges": int(len(edges))}


def build_causal_edge_table(
    omnipath_path: str | Path,
    signor_path: str | Path,
    *,
    conflict_rule: str | None = None,
) -> tuple[pd.DataFrame, dict[str, object]]:
    conflict_rule = conflict_rule or load_conflict_rule()
    if conflict_rule not in {"retain_as_ambiguous", "drop"}:
        raise StopAndReportError(f"Unknown conflict rule: {conflict_rule}", reason="bad-conflict-rule")
    omni, omni_stats = parse_omnipath(omnipath_path, conflict_rule=conflict_rule)
    signor, signor_stats = parse_signor(signor_path)
    combined = pd.concat([omni, signor], ignore_index=True)
    # Deduplicate identical source/target/sign/relation keeping max confidence
    combined = (
        combined.sort_values("confidence", ascending=False)
        .drop_duplicates(subset=["source", "target", "sign", "relation_type", "provenance"], keep="first")
        .reset_index(drop=True)
    )
    summary = {
        "conflict_rule": conflict_rule,
        "omnipath": omni_stats,
        "signor": signor_stats,
        "n_causal_edges": int(len(combined)),
        "sign_balance": combined["sign"].value_counts().to_dict(),
        "relation_type_distribution": combined["relation_type"].value_counts().to_dict(),
        "provenance_counts": combined["provenance"].value_counts().to_dict(),
    }
    return combined, summary


def build_string_edge_table(
    links_gz: str | Path,
    info_gz: str | Path,
    *,
    score_min: int | None = None,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Undirected auxiliary STRING layer — never merged into causal signed edges."""

    score_min = score_min if score_min is not None else load_string_score_min()
    info = pd.read_csv(info_gz, sep="\t", compression="gzip", dtype=str)
    # columns: #string_protein_id preferred_name ...
    id_col = [c for c in info.columns if "string_protein_id" in c.lower() or c.startswith("#")][0]
    name_col = "preferred_name" if "preferred_name" in info.columns else info.columns[1]
    id_to_symbol = dict(zip(info[id_col].astype(str), info[name_col].astype(str), strict=False))

    rows = []
    with gzip.open(links_gz, "rt") as handle:
        header = handle.readline().strip().split()
        # protein1 protein2 combined_score
        for line in handle:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            p1, p2, score_s = parts[0], parts[1], parts[2]
            score = int(score_s)
            if score < score_min:
                continue
            s1 = id_to_symbol.get(p1)
            s2 = id_to_symbol.get(p2)
            if not s1 or not s2:
                continue
            rows.append(
                {
                    "source": s1,
                    "target": s2,
                    "sign": 0,
                    "relation_type": "string_physical",
                    "confidence": score / 1000.0,
                    "provenance": "string_v12_physical",
                    "undirected": True,
                }
            )
    edges = pd.DataFrame(rows)
    summary = {
        "score_min": score_min,
        "n_edges": int(len(edges)),
        "layer": "auxiliary_undirected_separate",
    }
    return edges, summary
