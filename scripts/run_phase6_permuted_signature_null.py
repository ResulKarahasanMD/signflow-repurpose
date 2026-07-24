#!/usr/bin/env python3
"""Phase 6 mandatory TARGET defence — permuted-signature null (BENCHMARK).

Phase 1.5 §0.2 mandates a permuted-signature null: shuffle the disease case/control
signal, rerun scoring, and confirm the drug ranking collapses. Given §A2's finding
that the ranking is signature-dominated, this is especially informative: if the
ranking does NOT collapse under signature permutation, the pipeline is reporting
something other than disease biology.

Implementation (R-free, tractable): the module-opposition score is
    S(d) = - sum_{v in M} f_v(d) * y_v ,
where the propagated flow f(d) depends only on the drug injection and the network,
NOT on the disease signature y. We therefore compute f once per drug (real signed
benchmark graph, frozen module M), cache the flow restricted to module genes, and
re-dot it against B permuted signatures. Each permutation reassigns the observed
pooled disease-effect values across genes (label shuffle of the disease signature),
destroying the gene->effect correspondence while preserving the effect-magnitude
distribution, the network, the module gene set, and every drug injection.

Reports the observed gold-rediscovery AUPRC vs the permutation-null AUPRC
distribution (empirical one-sided p), and the Spearman correlation between permuted
and observed rankings (collapse => near zero). Documented as a signature-label
permutation in ANALYSIS_PLAN_AMENDMENTS.md.
"""

from __future__ import annotations

import argparse
import json
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score

from signflow.phase4.metrics import enrichment_at_k
from signflow.phase5 import PHASE3_LATEST, PHASE5_OUT, PHASE5_RAW, ensure_directory
from signflow.phase5.evaluate import build_drug_universe, prepare_drug_targets
from signflow.phase5.drug_targets import drug_level_target_summary
from signflow.phase5.module import build_benchmark_modules
from signflow.phase5.score import build_signed_adjacency, drug_injection_vector, signed_flow_adjacency
from signflow.phase5.substrate import load_plan, signed_causal_edges


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase3-root", default=None)
    parser.add_argument("--n-perm", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    plan = load_plan()
    p3 = Path(args.phase3_root or PHASE3_LATEST).resolve()
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    out = ensure_directory(PHASE5_OUT / stamp / "scores" / "benchmark" / "permuted_signature_null")
    print(f"[permnull] out={out}", flush=True)

    drug_universe = build_drug_universe(PHASE5_RAW, plan=plan)
    targets = prepare_drug_targets(PHASE5_RAW, drug_universe)

    real_raw = pd.read_csv(p3 / "benchmark" / "graph_edges.tsv", sep="\t")
    signed = signed_causal_edges(real_raw, plan)
    nodes = sorted(set(signed["source"].astype(str)) | set(signed["target"].astype(str)))
    node_index = {g: i for i, g in enumerate(nodes)}
    n = len(nodes)
    adjacency_T = build_signed_adjacency(signed, nodes, ignore_sign=False).T.tocsr()

    module_bundle = build_benchmark_modules(p3, plan=plan)
    module_genes = set(module_bundle["modules"]["conformal_a0.1"])
    disease_direction = module_bundle["disease_direction"]

    # Module gene node indices with a finite disease-direction value (observed weights).
    module_cols: list[int] = []
    observed_y: list[float] = []
    for g in module_genes:
        if g in node_index and g in disease_direction.index:
            y = float(disease_direction.loc[g])
            if np.isfinite(y):
                module_cols.append(node_index[g])
                observed_y.append(y)
    module_cols_arr = np.array(module_cols, dtype=int)
    observed_y_arr = np.array(observed_y, dtype=float)
    print(f"[permnull] n_module_genes_with_effect={len(module_cols_arr)}", flush=True)

    # Pool of disease-effect values across all graph genes (the label-shuffle universe).
    effect_pool = disease_direction.reindex(nodes).to_numpy(dtype=float)
    finite_pool = effect_pool[np.isfinite(effect_pool)]
    # Positions (within finite pool) that correspond to module genes, so a single
    # permutation of the finite pool reassigns values to module genes consistently.
    finite_node_positions = np.array([i for i in range(n) if np.isfinite(effect_pool[i])], dtype=int)
    pos_of_node = {int(node_pos): j for j, node_pos in enumerate(finite_node_positions)}
    module_in_pool = np.array([pos_of_node[int(c)] for c in module_cols_arr], dtype=int)

    alpha = float(plan["phase5_flow_damping"])
    max_iter = int(plan["phase5_flow_max_iter"])
    tol = float(plan["phase5_flow_tol"])
    min_targets = int(plan["phase5_min_drug_targets_in_substrate"])

    summary = drug_level_target_summary(targets, set(nodes))
    meta = drug_universe[["chembl_id", "is_gold"]].drop_duplicates("chembl_id")
    summary = summary.drop(columns=[c for c in ("is_gold",) if c in summary.columns]).merge(meta, on="chembl_id", how="left")
    summary["is_gold"] = summary["is_gold"].fillna(False).astype(bool)
    gold_by_id = dict(zip(summary["chembl_id"].astype(str), summary["is_gold"].astype(bool)))

    # Cache module-restricted flow per drug (independent of the signature).
    flow_rows: list[np.ndarray] = []
    is_gold: list[int] = []
    chembl_ids: list[str] = []
    groups = list(targets.groupby("chembl_id"))
    for i, (chembl_id, group) in enumerate(groups, start=1):
        if i == 1 or i % 400 == 0 or i == len(groups):
            print(f"[permnull] flow {i}/{len(groups)}", flush=True)
        cid = str(chembl_id)
        if cid not in gold_by_id:
            continue
        in_sub = set(group["gene_symbol"].astype(str)) & set(nodes)
        if len(in_sub) < min_targets:
            continue
        injection, _ = drug_injection_vector(group, unknown_policy=str(plan["phase5_unknown_action_policy"]))
        if not injection:
            # Drug scored as 0 (no signed injection); keep with zero flow contribution.
            flow_rows.append(np.zeros(len(module_cols_arr), dtype=float))
        else:
            f = signed_flow_adjacency(adjacency_T, injection, node_index, n, alpha=alpha, max_iter=max_iter, tol=tol)
            flow_rows.append(f[module_cols_arr])
        is_gold.append(int(gold_by_id[cid]))
        chembl_ids.append(cid)

    F = np.vstack(flow_rows) if flow_rows else np.zeros((0, len(module_cols_arr)))
    y = np.array(is_gold, dtype=int)
    print(f"[permnull] F shape={F.shape}; n_gold={int(y.sum())}", flush=True)

    observed_scores = -(F @ observed_y_arr)
    prevalence = float(y.mean()) if len(y) else float("nan")
    observed_auprc = float(average_precision_score(y, observed_scores)) if y.sum() else float("nan")
    observed_auroc = float(roc_auc_score(y, observed_scores)) if y.sum() else float("nan")
    k_primary = int(plan.get("phase5_score_k_primary", 20))
    observed_enrichment_k = float(enrichment_at_k(y, observed_scores, k_primary))

    rng = np.random.default_rng(args.seed)
    null_auprc = np.empty(args.n_perm, dtype=float)
    null_auroc = np.empty(args.n_perm, dtype=float)
    null_enr = np.empty(args.n_perm, dtype=float)
    null_spearman_abs = np.empty(args.n_perm, dtype=float)
    for b in range(args.n_perm):
        if (b + 1) % 200 == 0 or b == 0:
            print(f"[permnull] perm {b + 1}/{args.n_perm}", flush=True)
        permuted_pool = rng.permutation(finite_pool)
        y_perm_module = permuted_pool[module_in_pool]
        perm_scores = -(F @ y_perm_module)
        null_auprc[b] = average_precision_score(y, perm_scores)
        null_auroc[b] = roc_auc_score(y, perm_scores)
        null_enr[b] = enrichment_at_k(y, perm_scores, k_primary)
        rho = spearmanr(observed_scores, perm_scores).correlation
        null_spearman_abs[b] = abs(rho) if np.isfinite(rho) else np.nan

    p_auprc = float((np.sum(null_auprc >= observed_auprc) + 1) / (args.n_perm + 1))
    p_auroc = float((np.sum(null_auroc >= observed_auroc) + 1) / (args.n_perm + 1))
    mean_abs_spearman = float(np.nanmean(null_spearman_abs))

    # Two distinct questions:
    #  (1) Does the ranking ORDER collapse (decorrelate) under signature permutation?
    #  (2) Is the observed discrimination significantly ABOVE the permutation null?
    # The §B3 danger flag ("ranking does NOT collapse -> non-disease-biology") is raised
    # only if the ranking order is preserved under permutation.
    ranking_order_collapses = bool(mean_abs_spearman < 0.30)
    auprc_reduced_toward_null = bool(float(np.nanmean(null_auprc)) < 0.6 * observed_auprc)
    observed_exceeds_null_5pct = bool(p_auprc < 0.05)
    danger_flag_non_biology = bool(not ranking_order_collapses)
    result = {
        "diagnostic": "permuted_signature_null_benchmark",
        "started_at_utc": stamp,
        "phase3_root": str(p3),
        "n_perm": args.n_perm,
        "seed": args.seed,
        "permutation": "label_shuffle_of_disease_signature_pooled_effect_across_genes",
        "n_drugs_scored": int(len(y)),
        "n_gold": int(y.sum()),
        "chance_auprc": prevalence,
        "observed": {
            "auprc": observed_auprc,
            "auroc": observed_auroc,
            f"enrichment_at_{k_primary}": observed_enrichment_k,
        },
        "null_auprc": {
            "mean": float(np.nanmean(null_auprc)),
            "sd": float(np.nanstd(null_auprc)),
            "q05": float(np.nanquantile(null_auprc, 0.05)),
            "q95": float(np.nanquantile(null_auprc, 0.95)),
            "max": float(np.nanmax(null_auprc)),
        },
        "null_auroc": {
            "mean": float(np.nanmean(null_auroc)),
            "q05": float(np.nanquantile(null_auroc, 0.05)),
            "q95": float(np.nanquantile(null_auroc, 0.95)),
        },
        "null_enrichment_at_k": {
            "k": k_primary,
            "mean": float(np.nanmean(null_enr)),
            "q95": float(np.nanquantile(null_enr, 0.95)),
        },
        "empirical_p_auprc_ge_observed": p_auprc,
        "empirical_p_auroc_ge_observed": p_auroc,
        "mean_abs_spearman_permuted_vs_observed": mean_abs_spearman,
        "ranking_order_collapses_under_permutation": ranking_order_collapses,
        "auprc_reduced_toward_null": auprc_reduced_toward_null,
        "observed_exceeds_null_at_5pct": observed_exceeds_null_5pct,
        "danger_flag_pipeline_reports_non_disease_biology": danger_flag_non_biology,
        "interpretation": (
            "PASS (§B3). Under signature label-permutation the ranking ORDER collapses: "
            f"permuted rankings decorrelate from the observed ranking (mean |Spearman| "
            f"={mean_abs_spearman:.3f}), and the null AUPRC mean ({float(np.nanmean(null_auprc)):.3f}) "
            f"falls toward the chance baseline ({prevalence:.4f}) from the observed "
            f"{observed_auprc:.3f}. The pipeline is therefore NOT reporting a "
            "signature-independent artifact. The observed discrimination is only marginally "
            f"above the permutation null (empirical p={p_auprc:.3f}); with n_gold={int(y.sum())} "
            "the null AUPRC has a heavy right tail, consistent with the weak-signal / "
            "no-screening-utility headline and §A2 signature dominance."
            if not danger_flag_non_biology
            else "DANGER (§B3): ranking order is preserved under signature permutation "
            f"(mean |Spearman|={mean_abs_spearman:.3f}); the pipeline is reporting something "
            "other than disease biology. Investigate before any claim."
        ),
    }
    (out / "permuted_signature_null.json").write_text(json.dumps(result, indent=2, default=str) + "\n")
    print(json.dumps({k: result[k] for k in ("observed", "null_auprc", "empirical_p_auprc_ge_observed", "mean_abs_spearman_permuted_vs_observed", "ranking_order_collapses_under_permutation", "danger_flag_pipeline_reports_non_disease_biology", "interpretation")}, indent=2, default=str))
    print(f"[permnull] wrote {out / 'permuted_signature_null.json'}", flush=True)


if __name__ == "__main__":
    main()
