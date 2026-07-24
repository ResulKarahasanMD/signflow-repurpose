#!/usr/bin/env python3
"""Phase 5 §A1 — BENCHMARK control-graph AUPRC/AUROC (decisive diagnostic).

Control-graph scoring ran for TARGET only; the benchmark control-graph AUPRC was
never reported. This is the number that determines what the Phase 5 result says.

Scores the frozen disease-agnostic drug pool on the benchmark real signed graph and
its three degree/structure-preserving controls (rewired, edge_shuffled, sign_shuffled),
holding the module and disease direction fixed so only the network changes. Reports
AUPRC and AUROC with n_boot=1000 bootstrap CIs, the chance baseline beside each, and
the Spearman correlation of drug scores between the real graph and each control.

Pre-committed interpretation (frozen before running, per prompt §A1):
  - If rewired/edge_shuffled AUPRC ~= real (CIs overlap / difference below the
    detectable threshold), the score carries NO network information: the modest
    discrimination comes from drug target sets and the disease signature, not from
    propagation.
  - If real materially exceeds rewired, the network contributes and the current
    framing stands.
Either outcome is reportable; the script commits to reporting whichever occurs.
"""

from __future__ import annotations

import argparse
import json
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd

from signflow.phase4.metrics import bootstrap_ci, ci_overlap, minimum_detectable_difference
from signflow.phase5 import PHASE3_LATEST, PHASE5_OUT, PHASE5_RAW, ensure_directory
from signflow.phase5.evaluate import (
    build_drug_universe,
    compute_ranking_metrics,
    prepare_drug_targets,
)
from signflow.phase5.drug_targets import drug_level_target_summary
from signflow.phase5.module import build_benchmark_modules, write_module_artifacts
from signflow.phase5.score import CachedFlowScorer, node_degrees
from signflow.phase5.substrate import load_plan, signed_causal_edges


def _score_pool_on_edges(
    *,
    variant: str,
    raw_edges: pd.DataFrame,
    targets: pd.DataFrame,
    drug_summary: pd.DataFrame,
    module_genes: set[str],
    disease_direction: pd.Series,
    plan: dict[str, object],
    out_path: Path,
) -> pd.DataFrame:
    """Score every eligible drug on `raw_edges` after the frozen signed filter.

    Node universe and eligibility (>= min targets in substrate) are recomputed on
    the control graph's own signed node set, exactly as the real substrate is built,
    so the comparison isolates network structure with module + direction held fixed.
    """

    signed_edges = signed_causal_edges(raw_edges, plan)
    nodes = sorted(set(signed_edges["source"].astype(str)) | set(signed_edges["target"].astype(str)))
    node_set = set(nodes)
    degrees = pd.Series(node_degrees(signed_edges, nodes), index=nodes)
    scorer = CachedFlowScorer(
        nodes,
        signed_edges,
        module_genes,
        disease_direction,
        alpha=float(plan["phase5_flow_damping"]),
        max_iter=int(plan["phase5_flow_max_iter"]),
        tol=float(plan["phase5_flow_tol"]),
        ignore_sign=False,
        variant=variant,
    )
    meta = drug_summary.set_index("chembl_id")
    min_targets = int(plan["phase5_min_drug_targets_in_substrate"])
    rows = []
    groups = list(targets.groupby("chembl_id"))
    for i, (chembl_id, group) in enumerate(groups, start=1):
        if i == 1 or i % 400 == 0 or i == len(groups):
            print(f"[a1] score {variant} {i}/{len(groups)}", flush=True)
        if chembl_id not in meta.index:
            continue
        m = meta.loc[chembl_id]
        in_sub = sorted(set(group["gene_symbol"].astype(str)) & node_set)
        if len(in_sub) < min_targets:
            continue
        score_val = scorer.score_drug_targets(group, unknown_policy=str(plan["phase5_unknown_action_policy"]))
        rows.append(
            {
                "chembl_id": chembl_id,
                "pref_name": m.get("pref_name", ""),
                "canonical_group": m.get("canonical_group", ""),
                "is_gold": bool(m.get("is_gold", False)),
                "gold_eligible": bool(m.get("is_gold", False)),
                "score": score_val,
                "score_variant": variant,
                "n_targets_in_substrate": len(in_sub),
                "mean_target_degree": float(degrees.loc[in_sub].mean()) if in_sub else float("nan"),
            }
        )
    ranking = pd.DataFrame(rows).sort_values(["score", "chembl_id"], ascending=[False, True]).reset_index(drop=True)
    ranking.to_csv(out_path, sep="\t", index=False)
    return ranking


def _y_s(ranking: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    ranked = ranking.sort_values(["score", "chembl_id"], ascending=[False, True]).reset_index(drop=True)
    is_gold = ranked["is_gold"].astype(bool)
    gold_eligible = ranked.get("gold_eligible", is_gold).astype(bool)
    y = (is_gold & gold_eligible).astype(int).to_numpy()
    s = ranked["score"].astype(float).to_numpy()
    return y, s


def _spearman_vs_real(real: pd.DataFrame, ctrl: pd.DataFrame) -> dict[str, object]:
    merged = real[["chembl_id", "score"]].merge(
        ctrl[["chembl_id", "score"]], on="chembl_id", suffixes=("_real", "_ctrl")
    )
    if len(merged) < 3:
        return {"rho": float("nan"), "n": int(len(merged))}
    rho = float(merged["score_real"].corr(merged["score_ctrl"], method="spearman"))
    return {"rho": rho, "n": int(len(merged))}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase3-root", default=None)
    parser.add_argument("--prior-stamp", default="20260724T075949Z")
    parser.add_argument("--n-boot", type=int, default=1000)
    args = parser.parse_args()

    plan = load_plan()
    p3 = Path(args.phase3_root or PHASE3_LATEST).resolve()
    bench_dir = p3 / "benchmark"
    ctrl_dir = bench_dir / "controls"

    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    out = ensure_directory(PHASE5_OUT / stamp / "scores" / "benchmark" / "control_graphs")
    print(f"[a1] out={out}", flush=True)
    print(f"[a1] phase3_root={p3}", flush=True)

    drug_universe = build_drug_universe(PHASE5_RAW, plan=plan)
    targets = prepare_drug_targets(PHASE5_RAW, drug_universe)

    # Frozen module + disease direction (held fixed across all graphs).
    module_bundle = build_benchmark_modules(p3, plan=plan)
    write_module_artifacts(module_bundle, out / "modules")
    module_genes = set(module_bundle["modules"]["conformal_a0.1"])
    disease_direction = module_bundle["disease_direction"]

    # Drug summary keyed to the REAL signed node set (eligibility base for reporting);
    # each graph recomputes its own eligibility internally.
    real_raw = pd.read_csv(bench_dir / "graph_edges.tsv", sep="\t")
    real_signed = signed_causal_edges(real_raw, plan)
    real_nodes = set(real_signed["source"].astype(str)) | set(real_signed["target"].astype(str))
    summary = drug_level_target_summary(targets, real_nodes)
    meta = drug_universe[["chembl_id", "is_gold", "canonical_group", "pref_name"]].drop_duplicates("chembl_id")
    summary = summary.drop(columns=[c for c in ("is_gold", "canonical_group", "pref_name") if c in summary.columns])
    summary = summary.merge(meta, on="chembl_id", how="left")
    summary["is_gold"] = summary["is_gold"].fillna(False).astype(bool)

    graphs = {
        "real": bench_dir / "graph_edges.tsv",
        "rewired": ctrl_dir / "edges_rewired.tsv",
        "edge_shuffled": ctrl_dir / "edges_edge_shuffled.tsv",
        "sign_shuffled": ctrl_dir / "edges_sign_shuffled.tsv",
    }

    rankings: dict[str, pd.DataFrame] = {}
    arms: dict[str, dict[str, object]] = {}
    for name, path in graphs.items():
        if not path.exists():
            arms[name] = {"status": "missing", "path": str(path)}
            continue
        raw = pd.read_csv(path, sep="\t")
        ranking = _score_pool_on_edges(
            variant=f"benchmark_{name}",
            raw_edges=raw,
            targets=targets,
            drug_summary=summary,
            module_genes=module_genes,
            disease_direction=disease_direction,
            plan=plan,
            out_path=out / f"ranking_{name}.tsv",
        )
        rankings[name] = ranking
        y, s = _y_s(ranking)
        metrics = compute_ranking_metrics(ranking, n_boot=args.n_boot)
        arms[name] = {
            "n_scored": int(len(ranking)),
            "n_gold_scored": int(y.sum()),
            "auprc": metrics.get("auprc"),
            "auprc_ci": metrics.get("auprc_ci"),
            "auroc": metrics.get("auroc"),
            "auroc_ci": metrics.get("auroc_ci"),
            "chance_auprc": metrics.get("chance_auprc"),
            "enrichment_at_20": metrics.get("enrichment_at_20"),
            "enrichment_at_100": metrics.get("enrichment_at_100"),
            "enrichment_at_250": metrics.get("enrichment_at_250"),
        }

    # Spearman of scores + detectable-difference interpretation vs the real graph.
    real = rankings.get("real")
    interpretation: dict[str, object] = {}
    if real is not None:
        real_y, real_s = _y_s(real)
        real_auprc_ci = arms["real"]["auprc_ci"]
        for name in ("rewired", "edge_shuffled", "sign_shuffled"):
            if name not in rankings:
                continue
            arms[name]["spearman_vs_real"] = _spearman_vs_real(real, rankings[name])
            ctrl_ci = arms[name]["auprc_ci"]
            overlaps = ci_overlap(real_auprc_ci, ctrl_ci)
            mdd = minimum_detectable_difference(
                float(real_auprc_ci.get("width", float("nan"))),
                float(ctrl_ci.get("width", float("nan"))),
            )
            observed_gap = float(arms["real"]["auprc"]) - float(arms[name]["auprc"])
            interpretation[name] = {
                "auprc_ci_overlaps_real": bool(overlaps),
                "observed_auprc_gap_real_minus_control": observed_gap,
                "minimum_detectable_auprc_difference": mdd,
                "real_materially_exceeds_control": bool((not overlaps) and observed_gap > 0),
            }

    # Pre-committed verdict: network information present only if the real graph's AUPRC
    # materially exceeds the structure-preserving controls (rewired / edge_shuffled),
    # which destroy real connectivity while preserving degree / edge multiset.
    structural_controls = [c for c in ("rewired", "edge_shuffled") if c in interpretation]
    network_info = any(interpretation[c]["real_materially_exceeds_control"] for c in structural_controls)
    verdict = {
        "network_information_present": bool(network_info),
        "committed_reading": (
            "Real signed graph materially exceeds structure-preserving controls: the "
            "network contributes and the current framing stands."
            if network_info
            else "Real signed graph does NOT materially exceed rewired / edge_shuffled "
            "controls: the control-graph score carries no network information. The modest "
            "discrimination comes from the drug target sets and the disease signature, not "
            "from propagation over the molecular network."
        ),
        "structural_controls_considered": structural_controls,
        "note": (
            "Pre-committed before running (prompt §A1). Detectable threshold = AUPRC "
            "bootstrap-CI overlap (n_boot=%d) plus positive observed gap." % args.n_boot
        ),
    }

    result = {
        "diagnostic": "A1_benchmark_control_graph_auprc",
        "started_at_utc": stamp,
        "phase3_root": str(p3),
        "output_dir": str(out),
        "n_boot": args.n_boot,
        "chance_auprc": arms.get("real", {}).get("chance_auprc"),
        "arms": arms,
        "interpretation_vs_real": interpretation,
        "verdict": verdict,
    }
    (out / "a1_control_graph_metrics.json").write_text(json.dumps(result, indent=2, default=str) + "\n")
    print(json.dumps({"verdict": verdict, "arms": {k: {"auprc": v.get("auprc"), "auroc": v.get("auroc")} for k, v in arms.items()}}, indent=2, default=str))
    print(f"[a1] wrote {out / 'a1_control_graph_metrics.json'}", flush=True)


if __name__ == "__main__":
    main()
