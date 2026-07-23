"""BENCHMARK UC labels via Phase 1 Open Targets verifier path; TARGET emits no labels."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import requests
import yaml

from signflow.config import load_config
from signflow.phase2 import PROJECT_ROOT
from signflow.resources import count_filtered_targets, fetch_ot_associated_targets

from . import FROZEN_PLAN, StopAndReportError, ensure_directory


def load_uc_positives(session: requests.Session | None = None) -> pd.DataFrame:
    """Load UC genetics positives; hard-fail if count drifts from frozen 396."""

    config = load_config(PROJECT_ROOT / "configs" / "default.yaml")
    plan = yaml.safe_load(FROZEN_PLAN.read_text())["plan"]
    allowlist = set(plan["uc_genetics_allowlist"])
    threshold = float(plan["label_score_threshold"])
    expected = int(plan["frozen_counts"]["uc_genetics_primary"])
    session = session or requests.Session()
    endpoint = config["resources"]["open_targets"]["endpoint"]
    _, rows = fetch_ot_associated_targets(session, endpoint, plan["benchmark_disease_id"])
    observed = count_filtered_targets(rows, allowlist, threshold)
    if observed != expected:
        raise StopAndReportError(
            f"UC label count drifted: expected {expected}, observed {observed}",
            reason="uc-label-drift",
            details={"expected": expected, "observed": observed},
        )
    out_rows = []
    for row in rows:
        scores = row.get("datatypeScores") or []
        allowed = [item["score"] for item in scores if item.get("id") in allowlist]
        if not allowed:
            continue
        score = max(allowed)
        if score < threshold:
            continue
        target = row.get("target") or {}
        out_rows.append(
            {
                "gene_symbol": target.get("approvedSymbol"),
                "ensembl_id": target.get("id"),
                "genetics_score": score,
                "label": 1,
                "pu_status": "positive",
            }
        )
    frame = pd.DataFrame(out_rows).dropna(subset=["gene_symbol"])
    if len(frame) != expected:
        raise StopAndReportError(
            f"UC positive table length {len(frame)} != frozen {expected}",
            reason="uc-label-drift",
            details={"table_n": len(frame), "expected": expected},
        )
    return frame


def degree_matched_negatives(
    positive_genes: list[str],
    candidate_genes: list[str],
    degrees: pd.Series,
    *,
    seed: int,
    tolerance: float,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Sample unlabeled genes with degree matched to positives (PU: unlabeled ≠ negative for training framing)."""

    rng = np.random.default_rng(seed)
    pos = [g for g in positive_genes if g in degrees.index]
    cand = [g for g in candidate_genes if g not in set(pos) and g in degrees.index]
    if not pos:
        raise StopAndReportError("No positives in degree index", reason="label-degree-empty")
    pos_deg = degrees.loc[pos].astype(float)
    selected = []
    for g in pos:
        target = float(degrees.loc[g])
        lo, hi = target * (1 - tolerance), target * (1 + tolerance)
        pool = [c for c in cand if c not in selected and lo <= float(degrees.loc[c]) <= hi]
        if not pool:
            # widen to nearest by absolute degree distance
            dists = sorted(cand, key=lambda c: abs(float(degrees.loc[c]) - target))
            pool = [c for c in dists if c not in selected][:20]
        if not pool:
            continue
        pick = str(rng.choice(pool))
        selected.append(pick)
    neg = pd.DataFrame(
        {
            "gene_symbol": selected,
            "label": 0,
            "pu_status": "degree_matched_unlabeled_control",
            "degree": [float(degrees.loc[g]) for g in selected],
        }
    )
    qc = {
        "n_positives": len(pos),
        "n_negatives": len(selected),
        "pos_degree_mean": float(pos_deg.mean()),
        "neg_degree_mean": float(neg["degree"].mean()) if len(neg) else float("nan"),
        "pos_degree_median": float(pos_deg.median()),
        "neg_degree_median": float(neg["degree"].median()) if len(neg) else float("nan"),
        "tolerance": tolerance,
        "note": "PU framing: these controls are for evaluation/baselines; unlabeled ≠ biological negative.",
    }
    return neg, qc


def write_benchmark_labels(
    output_dir: str | Path,
    node_genes: list[str],
    degrees: pd.Series,
    *,
    seed: int,
    tolerance: float,
) -> dict[str, object]:
    output_dir = ensure_directory(output_dir)
    positives = load_uc_positives()
    positives.to_csv(Path(output_dir) / "uc_positives.tsv", sep="\t", index=False)
    in_graph = positives[positives["gene_symbol"].isin(node_genes)].copy()
    in_graph.to_csv(Path(output_dir) / "uc_positives_in_graph.tsv", sep="\t", index=False)
    neg, qc = degree_matched_negatives(
        in_graph["gene_symbol"].tolist(),
        node_genes,
        degrees,
        seed=seed,
        tolerance=tolerance,
    )
    neg.to_csv(Path(output_dir) / "degree_matched_controls.tsv", sep="\t", index=False)
    pd.DataFrame([qc]).to_csv(Path(output_dir) / "degree_match_qc.tsv", sep="\t", index=False)
    meta = {
        "genetics_channel_excluded_from_convergence_rank": True,
        "a6_2": "genetics channel used for labels excluded from later triangulation/convergence rank",
        "n_positives_total": int(len(positives)),
        "n_positives_in_graph": int(len(in_graph)),
        "degree_match_qc": qc,
    }
    return meta
