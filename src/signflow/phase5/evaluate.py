"""Score drugs, nulls, rankings, and gold rediscovery metrics."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import average_precision_score, roc_auc_score

from signflow.phase4.metrics import bootstrap_ci, enrichment_at_k, paired_difference_bootstrap, precision_at_k

from . import ensure_directory
from .acquire import assert_pool_guards, load_approved_drug_universe, load_gold_drugs
from .controls import select_unrelated_negative_drugs
from .drug_targets import build_drug_target_table, drug_level_target_summary
from .module import build_benchmark_modules, load_target_modules, write_module_artifacts
from .nulls import attach_nulls_and_fdr, spearman_score_degree
from .score import CachedFlowScorer, drug_injection_vector, module_opposition_score, node_degrees, signed_flow
from .substrate import load_phase3_substrates, load_plan


def build_drug_universe(raw_root: str | Path, plan: dict[str, object] | None = None) -> pd.DataFrame:
    """Disease-agnostic approved-drug universe ∪ frozen gold labels.

    Disease-conditioned OT candidates are never used as the ranking pool.
    """

    plan = plan or load_plan()
    if str(plan.get("phase5_candidate_pool_disease_conditioning")) != "forbidden":
        from . import StopAndReportError

        raise StopAndReportError(
            "phase5_candidate_pool_disease_conditioning must be forbidden",
            reason="policy-drift",
        )
    uni = load_approved_drug_universe(Path(raw_root))
    if "canonical_group" not in uni.columns:
        uni["canonical_group"] = uni["pref_name"].astype(str).str.lower().str.replace(" ", "_")
    if "is_gold" not in uni.columns:
        gold_ids = set(load_gold_drugs()["chembl_id"].astype(str))
        uni["is_gold"] = uni["chembl_id"].astype(str).isin(gold_ids)
    assert_pool_guards(uni, plan=plan)
    return uni


def _score_one_drug(
    drug_targets: pd.DataFrame,
    *,
    nodes: list[str],
    signed_edges: pd.DataFrame,
    module_genes: set[str],
    disease_direction: pd.Series,
    plan: dict[str, object],
    ignore_sign: bool = False,
) -> dict[str, object]:
    injection, diag = drug_injection_vector(
        drug_targets, unknown_policy=str(plan["phase5_unknown_action_policy"])
    )
    flow = signed_flow(
        nodes,
        signed_edges,
        injection,
        alpha=float(plan["phase5_flow_damping"]),
        max_iter=int(plan["phase5_flow_max_iter"]),
        tol=float(plan["phase5_flow_tol"]),
        ignore_sign=ignore_sign,
    )
    score = module_opposition_score(flow, nodes, module_genes, disease_direction)
    in_sub = sorted(set(drug_targets["gene_symbol"].astype(str)) & set(nodes))
    signed_genes = drug_targets.loc[drug_targets["is_signed"], "gene_symbol"].astype(str).tolist()
    signed_in_sub = sorted(set(signed_genes) & set(nodes))
    return {
        "score": float(score),
        "injection_diag": diag,
        "targets_in_substrate": in_sub,
        "n_targets_in_substrate": len(in_sub),
        "n_targets_signed_in_substrate": len(signed_in_sub),
    }


def compute_ranking_metrics(
    ranking: pd.DataFrame,
    *,
    k_values: list[int] | None = None,
    n_boot: int = 1000,
    score_col: str = "score",
) -> dict[str, object]:
    """Enrichment/recall@k plus pool-wide AUROC/AUPRC with chance baselines and CIs."""

    plan = load_plan()
    k_values = k_values or list(plan.get("phase5_score_k_values", [20, 50, 100, 250]))
    ranked = ranking.sort_values([score_col, "chembl_id"], ascending=[False, True]).reset_index(drop=True)
    is_gold = ranked.get("is_gold", False)
    if not isinstance(is_gold, pd.Series):
        is_gold = pd.Series([False] * len(ranked))
    is_gold = is_gold.astype(bool)
    gold_eligible = ranked.get("gold_eligible", is_gold)
    if not isinstance(gold_eligible, pd.Series):
        gold_eligible = is_gold
    gold_eligible = gold_eligible.astype(bool)

    y = (is_gold & gold_eligible).astype(int).to_numpy()
    s = ranked[score_col].astype(float).to_numpy()
    prevalence = float(y.mean()) if len(y) else float("nan")
    metrics: dict[str, object] = {
        "n_scored_universe": int(len(ranked)),
        "n_gold_total": int(is_gold.sum()),
        "n_gold_eligible_denominator": int((is_gold & gold_eligible).sum()),
        "n_gold_coverage_misses": int((is_gold & ~gold_eligible).sum()),
        "n_non_gold_scored": int((~is_gold).sum()),
        "gold_prevalence": prevalence,
        "chance_precision": prevalence,
        "chance_auprc": prevalence,
        "chance_auroc": 0.5,
        "note": (
            "Enrichment@k / precision@k / AUROC/AUPRC use the full scored universe as background. "
            "Chance baselines reported beside every figure."
        ),
    }
    if y.sum() == 0 or y.sum() == len(y):
        metrics["warning"] = "Degenerate gold labels for rediscovery metrics"
        return metrics

    for k in k_values:
        metrics[f"precision_at_{k}"] = float(precision_at_k(y, s, k))
        metrics[f"enrichment_at_{k}"] = float(enrichment_at_k(y, s, k))
        top_k = np.argsort(s)[::-1][: min(k, len(s))]
        metrics[f"recall_at_{k}"] = float(y[top_k].sum() / max(y.sum(), 1))
        metrics[f"expected_gold_in_top_{k}_at_chance"] = float(prevalence * min(k, len(s)))

    metrics["auroc"] = float(roc_auc_score(y, s))
    metrics["auprc"] = float(average_precision_score(y, s))
    metrics["auroc_ci"] = bootstrap_ci(y, s, metric="auroc", seed=0, n_boot=n_boot)
    metrics["auprc_ci"] = bootstrap_ci(y, s, metric="auprc", seed=0, n_boot=n_boot)

    gold_groups = set(ranked.loc[is_gold & gold_eligible, "canonical_group"].astype(str))
    for k in k_values:
        top = ranked.head(k)
        hit_groups = set(top.loc[top["is_gold"].astype(bool), "canonical_group"].astype(str))
        metrics[f"canonical_group_hits_at_{k}"] = sorted(hit_groups)
        metrics[f"n_canonical_group_hits_at_{k}"] = len(hit_groups)
        metrics[f"canonical_group_recall_at_{k}"] = float(len(hit_groups & gold_groups) / max(len(gold_groups), 1))

    gold_ranks = []
    for i, row in ranked.iterrows():
        if bool(row.get("is_gold")):
            gold_ranks.append(
                {
                    "rank": int(i) + 1,
                    "chembl_id": str(row["chembl_id"]),
                    "pref_name": str(row.get("pref_name", "")),
                    "canonical_group": str(row.get("canonical_group", "")),
                    "score": float(row[score_col]),
                    "bh_fdr": float(row["bh_fdr"])
                    if "bh_fdr" in ranked.columns and pd.notna(row.get("bh_fdr"))
                    else None,
                    "gold_eligible": bool(row.get("gold_eligible", True)),
                }
            )
    metrics["gold_ranks"] = gold_ranks
    return metrics


# Back-compat alias
def compute_gold_rediscovery(ranking: pd.DataFrame, *, k_values: list[int] | None = None) -> dict[str, object]:
    return compute_ranking_metrics(ranking, k_values=k_values)


def _rank_variant(
    *,
    variant: str,
    ignore_sign: bool,
    targets: pd.DataFrame,
    drug_summary: pd.DataFrame,
    nodes: list[str],
    signed_edges: pd.DataFrame,
    module_genes: set[str],
    disease_direction: pd.Series,
    degrees: pd.Series,
    plan: dict[str, object],
    out_dir: Path,
) -> tuple[pd.DataFrame, dict[str, object]]:
    scorer = CachedFlowScorer(
        nodes,
        signed_edges,
        module_genes,
        disease_direction,
        alpha=float(plan["phase5_flow_damping"]),
        max_iter=int(plan["phase5_flow_max_iter"]),
        tol=float(plan["phase5_flow_tol"]),
        ignore_sign=ignore_sign,
        variant=variant,
    )
    rows = []
    groups = list(targets.groupby("chembl_id"))
    for i, (chembl_id, group) in enumerate(groups, start=1):
        if i == 1 or i % 200 == 0 or i == len(groups):
            print(f"[phase5] score {variant} observed {i}/{len(groups)}", flush=True)
        meta_rows = drug_summary.loc[drug_summary["chembl_id"] == chembl_id]
        if meta_rows.empty:
            continue
        meta = meta_rows.iloc[0]
        if int(meta["n_targets_in_substrate"]) < int(plan["phase5_min_drug_targets_in_substrate"]):
            continue
        in_sub = sorted(set(group["gene_symbol"].astype(str)) & set(nodes))
        score_val = scorer.score_drug_targets(group, unknown_policy=str(plan["phase5_unknown_action_policy"]))
        signed_in_sub = sorted(set(group.loc[group["is_signed"], "gene_symbol"].astype(str)) & set(nodes))
        rows.append(
            {
                "chembl_id": chembl_id,
                "pref_name": meta.get("pref_name", ""),
                "canonical_group": meta.get("canonical_group", ""),
                "is_gold": bool(meta.get("is_gold", False)),
                "gold_eligible": bool(meta.get("is_gold", False)),
                "score": score_val,
                "score_variant": variant,
                "n_targets_in_substrate": len(in_sub),
                "n_targets_signed_in_substrate": len(signed_in_sub),
                "targets_in_substrate": in_sub,
                "mean_target_degree": float(np.mean([degrees[g] for g in in_sub])) if in_sub else float("nan"),
            }
        )
    ranking = pd.DataFrame(rows)

    def score_fn(gene_list: list[str]) -> float:
        return scorer.score_genes(gene_list)

    if len(ranking) and int(plan.get("null_permutations", 0)) > 0:
        print(f"[phase5] nulls {variant} n_drugs={len(ranking)} n_perm={plan['null_permutations']}", flush=True)
        ranking = attach_nulls_and_fdr(
            ranking,
            degrees=degrees,
            score_fn_for_genes=score_fn,
            n_perm=int(plan["null_permutations"]),
            tolerance=float(plan["degree_match_tolerance"]),
            seed=0,
        )
        sp = spearman_score_degree(ranking["score"].to_numpy(), ranking["mean_target_degree"].to_numpy())
    elif len(ranking):
        ranking = ranking.copy()
        ranking["bh_fdr"] = np.nan
        ranking["z_score"] = np.nan
        ranking["null_mean"] = np.nan
        ranking["null_sd"] = np.nan
        sp = spearman_score_degree(ranking["score"].to_numpy(), ranking["mean_target_degree"].to_numpy())
    else:
        sp = {"rho": float("nan"), "p": float("nan")}

    ranking = ranking.sort_values("score", ascending=False)
    rank_path = out_dir / f"ranking_{variant}.tsv"
    ranking.to_csv(rank_path, sep="\t", index=False)
    rediscovery = compute_ranking_metrics(
        ranking,
        k_values=list(plan.get("phase5_score_k_values", [20, 50, 100, 250])),
        n_boot=int(plan.get("phase5_bootstrap_n", 1000)),
    )
    n_bh = (
        int((ranking["bh_fdr"] < 0.05).sum())
        if len(ranking) and "bh_fdr" in ranking.columns and ranking["bh_fdr"].notna().any()
        else 0
    )
    summary = {
        "variant": variant,
        "ignore_sign": ignore_sign,
        "ranking_path": str(rank_path),
        "n_scored_drugs": int(len(ranking)),
        "spearman_score_degree": sp,
        "gold_rediscovery": rediscovery,
        "n_bh_fdr_survivors": n_bh,
    }
    return ranking, summary


def compare_signed_unsigned(
    signed_ranking: pd.DataFrame,
    unsigned_ranking: pd.DataFrame,
    *,
    n_boot: int = 1000,
    unsigned_label: str = "unsigned",
) -> dict[str, object]:
    """Paired bootstrap over gold set for enrichment@k and AUROC/AUPRC; decision rule applied."""

    plan = load_plan()
    k_primary = int(plan.get("phase5_score_k_primary", 20))
    merged = signed_ranking[["chembl_id", "is_gold", "score"]].merge(
        unsigned_ranking[["chembl_id", "score"]],
        on="chembl_id",
        suffixes=("_signed", f"_{unsigned_label}"),
    )
    if merged.empty:
        return {"status": "empty", "decision": "no_measurable_sign_benefit"}
    y = merged["is_gold"].astype(bool).astype(int).to_numpy()
    s_signed = merged["score_signed"].astype(float).to_numpy()
    s_unsigned = merged[f"score_{unsigned_label}"].astype(float).to_numpy()
    enr_signed = float(enrichment_at_k(y, s_signed, k_primary))
    enr_unsigned = float(enrichment_at_k(y, s_unsigned, k_primary))
    paired_auprc = paired_difference_bootstrap(y, s_signed, s_unsigned, metric="auprc", seed=0, n_boot=n_boot)
    # Detectable threshold: lower CI of signed−unsigned enrichment proxy via paired AUPRC
    exceeds = bool(paired_auprc.get("lo", 0) > 0) and enr_signed > enr_unsigned
    decision = (
        "signed_exceeds_unsigned"
        if exceeds
        else "sign_contributed_no_measurable_benefit_to_drug_scoring"
    )
    return {
        "k_primary": k_primary,
        "unsigned_arm": unsigned_label,
        "enrichment_signed": enr_signed,
        "enrichment_unsigned": enr_unsigned,
        "auroc_signed": float(roc_auc_score(y, s_signed)) if len(np.unique(y)) > 1 else float("nan"),
        "auroc_unsigned": float(roc_auc_score(y, s_unsigned)) if len(np.unique(y)) > 1 else float("nan"),
        "auprc_signed": float(average_precision_score(y, s_signed)) if len(np.unique(y)) > 1 else float("nan"),
        "auprc_unsigned": float(average_precision_score(y, s_unsigned)) if len(np.unique(y)) > 1 else float("nan"),
        "paired_auprc_signed_minus_unsigned": paired_auprc,
        "decision": decision,
        "primary_variant_retained": "signed",
        "note": (
            "Signed remains the pre-registered primary even if it underperforms. "
            "Decision reports whether sign contributed measurable benefit. "
            "Binding sign contrast uses unsigned_matched (same reduced edge set)."
        ),
    }


def compare_connectivity_cost(
    unsigned_matched: pd.DataFrame,
    unsigned_full: pd.DataFrame,
    *,
    n_boot: int = 1000,
) -> dict[str, object]:
    """Quantify the connectivity cost of exclude_from_signed_propagation."""

    plan = load_plan()
    k_primary = int(plan.get("phase5_score_k_primary", 20))
    merged = unsigned_matched[["chembl_id", "is_gold", "score"]].merge(
        unsigned_full[["chembl_id", "score"]],
        on="chembl_id",
        suffixes=("_matched", "_full"),
    )
    if merged.empty:
        return {"status": "empty"}
    y = merged["is_gold"].astype(bool).astype(int).to_numpy()
    s_m = merged["score_matched"].astype(float).to_numpy()
    s_f = merged["score_full"].astype(float).to_numpy()
    paired = paired_difference_bootstrap(y, s_m, s_f, metric="auprc", seed=1, n_boot=n_boot)
    return {
        "k_primary": k_primary,
        "enrichment_unsigned_matched": float(enrichment_at_k(y, s_m, k_primary)),
        "enrichment_unsigned_full": float(enrichment_at_k(y, s_f, k_primary)),
        "auroc_unsigned_matched": float(roc_auc_score(y, s_m)) if len(np.unique(y)) > 1 else float("nan"),
        "auroc_unsigned_full": float(roc_auc_score(y, s_f)) if len(np.unique(y)) > 1 else float("nan"),
        "auprc_unsigned_matched": float(average_precision_score(y, s_m)) if len(np.unique(y)) > 1 else float("nan"),
        "auprc_unsigned_full": float(average_precision_score(y, s_f)) if len(np.unique(y)) > 1 else float("nan"),
        "paired_auprc_matched_minus_full": paired,
        "interpretation": (
            "unsigned_matched vs unsigned_full isolates the connectivity cost of "
            "excluding ambiguous + unsigned_directed edges from the signed substrate."
        ),
    }


def acquired_pool_sensitivity(
    scored_ranking: pd.DataFrame,
    drug_universe: pd.DataFrame,
    *,
    k_values: list[int] | None = None,
) -> dict[str, object]:
    """Enrichment against the acquired pool; uncovered drugs counted as misses (score=-inf)."""

    plan = load_plan()
    k_values = k_values or list(plan.get("phase5_score_k_values", [20, 50, 100, 250]))
    uni = drug_universe[["chembl_id", "is_gold", "canonical_group", "pref_name"]].drop_duplicates("chembl_id")
    scored = scored_ranking[["chembl_id", "score"]].copy()
    merged = uni.merge(scored, on="chembl_id", how="left")
    merged["uncovered"] = merged["score"].isna()
    # Finite sentinel so AUROC/AUPRC remain defined; uncovered drugs sort last.
    miss_score = float(np.nanmin(merged["score"].to_numpy(dtype=float))) - 1.0 if merged["score"].notna().any() else -1.0
    if not np.isfinite(miss_score):
        miss_score = -1.0
    merged["score"] = merged["score"].fillna(miss_score)
    merged["gold_eligible"] = merged["is_gold"].astype(bool)
    metrics = compute_ranking_metrics(merged, k_values=k_values)
    metrics["n_acquired_universe"] = int(len(merged))
    metrics["n_uncovered"] = int(merged["uncovered"].sum())
    metrics["n_uncovered_gold"] = int((merged["uncovered"] & merged["is_gold"].astype(bool)).sum())
    metrics["note"] = (
        "Sensitivity vs acquired pool (n includes uncovered drugs as score=-inf misses). "
        "Primary rediscovery metrics remain scored-pool only."
    )
    return metrics


def differential_coverage_rates(
    coverage_report: dict[str, object],
    drug_summary: pd.DataFrame,
    *,
    n_acquired_universe: int | None = None,
) -> dict[str, object]:
    """Gold vs non-gold substrate coverage rates (A3).

    Relative inflation uses the acquired ranking universe (ChEMBL max_phase=4 ∪ gold),
    not the mapped drug-summary size.
    """

    gold = drug_summary["is_gold"].astype(bool)
    eligible = drug_summary["n_targets_in_substrate"] >= 1
    n_gold = int(gold.sum())
    n_non = int((~gold).sum())
    gold_cov = float(eligible[gold].mean()) if n_gold else float("nan")
    non_cov = float(eligible[~gold].mean()) if n_non else float("nan")
    scored_n = int(eligible.sum())
    gold_in_scored = int((gold & eligible).sum())
    acquired_n = int(n_acquired_universe) if n_acquired_universe is not None else int(
        coverage_report.get("n_acquired_universe")
        or coverage_report.get("n_universe")
        or coverage_report.get("n_drugs_total")
        or len(drug_summary)
    )
    gold_prev_acquired = float(n_gold / max(acquired_n, 1))
    gold_prev_scored = float(gold_in_scored / max(scored_n, 1))
    return {
        "n_gold": n_gold,
        "n_non_gold": n_non,
        "n_acquired_universe": acquired_n,
        "gold_coverage_rate": gold_cov,
        "non_gold_coverage_rate": non_cov,
        "gold_prevalence_acquired": gold_prev_acquired,
        "gold_prevalence_scored": gold_prev_scored,
        "relative_inflation": float(gold_prev_scored / max(gold_prev_acquired, 1e-12) - 1.0),
        "note": (
            "Gold prevalence rises from acquired universe to scored pool when substrate "
            "coverage correlates with approval status; report rates separately. "
            "Acquired denominator = disease-agnostic approved universe ∪ gold."
        ),
    }


def run_benchmark_evaluation(
    *,
    phase3_root: str | Path | None,
    output_dir: str | Path,
    drug_universe: pd.DataFrame,
    targets: pd.DataFrame,
    substrates: dict[str, object],
    acquire_info: dict[str, object] | None = None,
    skip_nulls_for_unsigned_full: bool = True,
) -> dict[str, object]:
    plan = load_plan()
    output_dir = ensure_directory(output_dir)
    bench_out = ensure_directory(output_dir / "benchmark")

    primary = substrates["primary"]
    nodes = sorted(primary["nodes"])
    signed_edges = primary["signed_edges"]
    full_unsigned_edges = primary.get("full_unsigned_edges")
    if full_unsigned_edges is None:
        from .substrate import unsigned_full_causal_edges

        full_unsigned_edges = unsigned_full_causal_edges(primary.get("causal_edges", signed_edges))
    full_nodes = sorted(set(primary.get("full_unsigned_nodes") or nodes) | set(nodes))
    degrees = pd.Series(node_degrees(signed_edges, nodes), index=nodes)

    # Pool guards after substrate mapping
    drug_summary = drug_level_target_summary(targets, set(nodes))
    meta = drug_universe[["chembl_id", "is_gold", "canonical_group", "pref_name"]].drop_duplicates("chembl_id")
    if "is_gold" in drug_summary.columns:
        drug_summary = drug_summary.drop(
            columns=[c for c in ("is_gold", "canonical_group", "pref_name") if c in drug_summary.columns]
        )
    drug_summary = drug_summary.merge(meta, on="chembl_id", how="left")
    n_in_sub = int((drug_summary["n_targets_in_substrate"] >= int(plan["phase5_min_drug_targets_in_substrate"])).sum())
    pool_guard = assert_pool_guards(drug_universe, n_in_substrate=n_in_sub, plan=plan)

    module_bundle = build_benchmark_modules(phase3_root, plan=plan)
    mod_paths = write_module_artifacts(module_bundle, bench_out / "modules")
    module_genes = set(module_bundle["modules"]["conformal_a0.1"])
    disease_direction = module_bundle["disease_direction"]

    signed_ranking, signed_summary = _rank_variant(
        variant="signed",
        ignore_sign=False,
        targets=targets,
        drug_summary=drug_summary,
        nodes=nodes,
        signed_edges=signed_edges,
        module_genes=module_genes,
        disease_direction=disease_direction,
        degrees=degrees,
        plan=plan,
        out_dir=bench_out,
    )
    # unsigned_matched = co-primary on the SAME reduced signed-only edge set
    unsigned_matched_ranking, unsigned_matched_summary = _rank_variant(
        variant="unsigned_matched",
        ignore_sign=True,
        targets=targets,
        drug_summary=drug_summary,
        nodes=nodes,
        signed_edges=signed_edges,
        module_genes=module_genes,
        disease_direction=disease_direction,
        degrees=degrees,
        plan=plan,
        out_dir=bench_out,
    )
    # Back-compat alias path expected by prior stamps / docs
    unsigned_alias = bench_out / "ranking_unsigned.tsv"
    unsigned_matched_ranking.to_csv(unsigned_alias, sep="\t", index=False)

    # unsigned_full = connectivity diagnostic on stim+inhib+unsigned_directed
    degrees_full = pd.Series(node_degrees(full_unsigned_edges, full_nodes), index=full_nodes)
    unsigned_full_plan = dict(plan)
    if skip_nulls_for_unsigned_full:
        # Comparison metrics (AUROC/AUPRC/enrichment) do not require degree-matched nulls.
        unsigned_full_plan = {**plan, "null_permutations": 0}
    unsigned_full_ranking, unsigned_full_summary = _rank_variant(
        variant="unsigned_full",
        ignore_sign=True,
        targets=targets,
        drug_summary=drug_summary,
        nodes=full_nodes,
        signed_edges=full_unsigned_edges,
        module_genes=module_genes,
        disease_direction=disease_direction,
        degrees=degrees_full,
        plan=unsigned_full_plan,
        out_dir=bench_out,
    )

    sign_comparison = compare_signed_unsigned(
        signed_ranking,
        unsigned_matched_ranking,
        n_boot=int(plan.get("phase5_bootstrap_n", 1000)),
        unsigned_label="unsigned_matched",
    )
    connectivity_comparison = compare_connectivity_cost(
        unsigned_matched_ranking,
        unsigned_full_ranking,
        n_boot=int(plan.get("phase5_bootstrap_n", 1000)),
    )
    three_arm = {
        "arms": [
            {
                "arm": "signed",
                "edge_set": "signed_only_reduced",
                "n_edges": int(len(signed_edges)),
                "sign_used": True,
                "auroc": signed_summary["gold_rediscovery"].get("auroc"),
                "auprc": signed_summary["gold_rediscovery"].get("auprc"),
                "chance_auprc": signed_summary["gold_rediscovery"].get("chance_auprc"),
                "enrichment_at_20": signed_summary["gold_rediscovery"].get("enrichment_at_20"),
            },
            {
                "arm": "unsigned_matched",
                "edge_set": "signed_only_reduced",
                "n_edges": int(len(signed_edges)),
                "sign_used": False,
                "auroc": unsigned_matched_summary["gold_rediscovery"].get("auroc"),
                "auprc": unsigned_matched_summary["gold_rediscovery"].get("auprc"),
                "chance_auprc": unsigned_matched_summary["gold_rediscovery"].get("chance_auprc"),
                "enrichment_at_20": unsigned_matched_summary["gold_rediscovery"].get("enrichment_at_20"),
            },
            {
                "arm": "unsigned_full",
                "edge_set": "full_causal_directed_excluding_ambiguous",
                "n_edges": int(len(full_unsigned_edges)),
                "sign_used": False,
                "auroc": unsigned_full_summary["gold_rediscovery"].get("auroc"),
                "auprc": unsigned_full_summary["gold_rediscovery"].get("auprc"),
                "chance_auprc": unsigned_full_summary["gold_rediscovery"].get("chance_auprc"),
                "enrichment_at_20": unsigned_full_summary["gold_rediscovery"].get("enrichment_at_20"),
            },
        ],
        "signed_vs_unsigned_matched": sign_comparison,
        "unsigned_matched_vs_unsigned_full": connectivity_comparison,
        "binding_sign_decision": sign_comparison.get("decision"),
        "note": plan.get("phase5_score_three_arm_note"),
    }

    ensure_directory(bench_out / "controls")
    negatives = select_unrelated_negative_drugs(drug_universe, plan=plan)
    negatives.to_csv(bench_out / "controls" / "unrelated_negatives.tsv", sep="\t", index=False)

    # Primary ranking alias (signed retained even if unsigned wins)
    primary_path = bench_out / "ranking_primary.tsv"
    signed_ranking.to_csv(primary_path, sep="\t", index=False)

    acquired_sens = acquired_pool_sensitivity(signed_ranking, drug_universe)
    diff_cov = differential_coverage_rates(
        {"n_drugs_total": int(len(drug_summary))},
        drug_summary,
        n_acquired_universe=int(len(drug_universe.drop_duplicates("chembl_id"))),
    )

    result = {
        "module_paths": mod_paths,
        "ranking_primary": str(primary_path),
        "ranking_signed": signed_summary["ranking_path"],
        "ranking_unsigned": str(unsigned_alias),
        "ranking_unsigned_matched": unsigned_matched_summary["ranking_path"],
        "ranking_unsigned_full": unsigned_full_summary["ranking_path"],
        "n_scored_drugs": signed_summary["n_scored_drugs"],
        "spearman_score_degree": signed_summary["spearman_score_degree"],
        "gold_rediscovery": signed_summary["gold_rediscovery"],
        "gold_rediscovery_unsigned": unsigned_matched_summary["gold_rediscovery"],
        "gold_rediscovery_unsigned_matched": unsigned_matched_summary["gold_rediscovery"],
        "gold_rediscovery_unsigned_full": unsigned_full_summary["gold_rediscovery"],
        "gold_rediscovery_acquired_pool_sensitivity": acquired_sens,
        "differential_coverage": diff_cov,
        "signed_summary": signed_summary,
        "unsigned_summary": unsigned_matched_summary,
        "unsigned_matched_summary": unsigned_matched_summary,
        "unsigned_full_summary": unsigned_full_summary,
        "signed_vs_unsigned": sign_comparison,
        "three_arm_comparison": three_arm,
        "n_bh_fdr_survivors": signed_summary["n_bh_fdr_survivors"],
        "primary_module": "conformal_a0.1",
        "substrate": primary["name"],
        "substrate_connectivity": primary.get("connectivity"),
        "substrate_connectivity_unsigned_full": primary.get("connectivity_unsigned_full"),
        "sensitivity_connectivity": substrates.get("sensitivity", {}).get("connectivity"),
        "edge_set_by_arm": {
            "signed": "signed_only_reduced",
            "unsigned_matched": "signed_only_reduced",
            "unsigned_full": "full_causal_directed_excluding_ambiguous",
        },
        "pool_guard": pool_guard,
        "acquire_info": acquire_info,
    }
    (bench_out / "benchmark_eval_summary.json").write_text(json.dumps(result, indent=2, default=str) + "\n")
    return result


def _load_soft_positives(plan: dict[str, object]) -> list[dict[str, object]]:
    path = Path(plan.get("phase5_soft_positives_path", "configs/target_soft_positives.yaml"))
    if not path.is_absolute():
        from . import PROJECT_ROOT

        path = PROJECT_ROOT / path
    if not path.exists():
        return []
    payload = yaml.safe_load(path.read_text()) or {}
    return list(payload.get("drugs", []) or [])


def report_soft_positive_ranks(ranking: pd.DataFrame, soft_positives: list[dict[str, object]]) -> list[dict[str, object]]:
    ranked = ranking.reset_index(drop=True)
    ranked["rank"] = np.arange(1, len(ranked) + 1)
    by_id = ranked.set_index("chembl_id")
    rows = []
    for sp in soft_positives:
        ids = [str(x) for x in sp.get("chembl_ids", [])]
        hit = None
        for cid in ids:
            if cid in by_id.index:
                hit = cid
                break
        if hit is None:
            rows.append(
                {
                    "name": sp.get("name"),
                    "chembl_ids": ids,
                    "status": "not_in_scored_pool",
                    "rank": None,
                    "score": None,
                    "note": sp.get("note"),
                }
            )
            continue
        row = by_id.loc[hit]
        rows.append(
            {
                "name": sp.get("name"),
                "chembl_id": hit,
                "chembl_ids": ids,
                "status": "scored",
                "rank": int(row["rank"]),
                "score": float(row["score"]),
                "pref_name": str(row.get("pref_name", "")),
                "note": sp.get("note"),
            }
        )
    return rows


def _score_ranking_on_edges(
    *,
    variant: str,
    targets: pd.DataFrame,
    drug_universe: pd.DataFrame,
    nodes: list[str],
    edges: pd.DataFrame,
    module_genes: set[str],
    disease_direction: pd.Series,
    plan: dict[str, object],
    out_path: Path,
    ignore_sign: bool = False,
) -> pd.DataFrame:
    scorer = CachedFlowScorer(
        nodes,
        edges,
        module_genes,
        disease_direction,
        alpha=float(plan["phase5_flow_damping"]),
        max_iter=int(plan["phase5_flow_max_iter"]),
        tol=float(plan["phase5_flow_tol"]),
        ignore_sign=ignore_sign,
        variant=variant,
    )
    meta = drug_universe.drop_duplicates("chembl_id").set_index("chembl_id")
    rows = []
    groups = list(targets.groupby("chembl_id"))
    for i, (chembl_id, group) in enumerate(groups, start=1):
        if i == 1 or i % 400 == 0 or i == len(groups):
            print(f"[phase5] target score {variant} {i}/{len(groups)}", flush=True)
        genes_in = set(group["gene_symbol"].astype(str)) & set(nodes)
        if len(genes_in) < int(plan["phase5_min_drug_targets_in_substrate"]):
            continue
        if chembl_id not in meta.index:
            continue
        score_val = scorer.score_drug_targets(group, unknown_policy=str(plan["phase5_unknown_action_policy"]))
        m = meta.loc[chembl_id]
        rows.append(
            {
                "chembl_id": chembl_id,
                "pref_name": m.get("pref_name", ""),
                "score": score_val,
                "n_targets_in_substrate": len(genes_in),
                "score_variant": variant,
            }
        )
    ranking = pd.DataFrame(rows).sort_values(["score", "chembl_id"], ascending=[False, True]).reset_index(drop=True)
    ranking.to_csv(out_path, sep="\t", index=False)
    return ranking


def _top_k_jaccard(a: pd.DataFrame, b: pd.DataFrame, k: int = 20) -> float:
    sa = set(a.head(k)["chembl_id"].astype(str))
    sb = set(b.head(k)["chembl_id"].astype(str))
    union = len(sa | sb)
    return float(len(sa & sb) / union) if union else 1.0


def run_target_evaluation(
    *,
    phase3_root: str | Path | None,
    output_dir: str | Path,
    drug_universe: pd.DataFrame,
    targets: pd.DataFrame,
    substrates: dict[str, object],
    run_h7: bool = True,
    run_control_graphs: bool = True,
) -> dict[str, object]:
    from .h7 import run_h7_futility
    from .substrate import signed_causal_edges

    plan = load_plan()
    output_dir = ensure_directory(output_dir)
    target_out = ensure_directory(output_dir / "target")
    root = Path(phase3_root).resolve() if phase3_root else None
    tmods = load_target_modules(root)
    write_module_artifacts(tmods, target_out / "modules")
    soft_positives = _load_soft_positives(plan)

    primary = substrates["primary"]
    sensitivity = substrates.get("sensitivity", primary)
    nodes_primary = sorted(primary["nodes"])
    signed_edges = primary["signed_edges"]
    disease_direction = tmods["disease_direction"]

    seed_sets = {
        "primary": set(tmods["modules"]["primary"]),
        "sensitivity_concordant": set(tmods["modules"]["sensitivity_concordant"]),
    }
    # Carry per-seed replication annotation through module meta
    seed_annotation = {}
    for name, frame in tmods.get("seed_frames", {}).items():
        seed_annotation[name] = {
            "n_seeds": int(len(frame)),
            "n_direction_agree": int(frame["direction_agree"].sum()) if "direction_agree" in frame.columns else None,
            "median_i2": float(frame["i2"].median()) if "i2" in frame.columns else None,
            "median_weight_share_GSE11783": float(frame["weight_share_GSE11783"].median())
            if "weight_share_GSE11783" in frame.columns
            else None,
            "genes": frame["gene_symbol"].astype(str).tolist(),
        }

    rankings: dict[str, pd.DataFrame] = {}
    soft_reports: dict[str, list] = {}
    for seed_name, module_genes in seed_sets.items():
        for sub_name, sub in (("primary_substrate", primary), ("sensitivity_substrate", sensitivity)):
            nodes = sorted(sub["nodes"])
            edges = sub["signed_edges"]
            variant = f"{seed_name}__{sub_name}"
            path = target_out / f"ranking_{variant}.tsv"
            ranking = _score_ranking_on_edges(
                variant=variant,
                targets=targets,
                drug_universe=drug_universe,
                nodes=nodes,
                edges=edges,
                module_genes=module_genes,
                disease_direction=disease_direction,
                plan=plan,
                out_path=path,
                ignore_sign=False,
            )
            rankings[variant] = ranking
            soft_reports[variant] = report_soft_positive_ranks(ranking, soft_positives)

    # Primary alias
    primary_ranking = rankings["primary__primary_substrate"]
    primary_path = target_out / "ranking_primary.tsv"
    primary_ranking.to_csv(primary_path, sep="\t", index=False)

    # Axis disagreement between seed sets / substrates (descriptive)
    axis_disagreement = {
        "top20_jaccard_primary_vs_sensitivity_seed_on_primary_substrate": _top_k_jaccard(
            rankings["primary__primary_substrate"],
            rankings["sensitivity_concordant__primary_substrate"],
            20,
        ),
        "top20_jaccard_primary_vs_sensitivity_substrate_on_primary_seeds": _top_k_jaccard(
            rankings["primary__primary_substrate"],
            rankings["primary__sensitivity_substrate"],
            20,
        ),
    }

    control_summaries = {}
    if run_control_graphs and root is not None:
        ctrl_dir = root / "target" / "controls"
        ctrl_out = ensure_directory(target_out / "controls")
        for variant in plan.get("phase5_control_graph_variants", ["rewired", "edge_shuffled", "sign_shuffled"]):
            edge_path = ctrl_dir / f"edges_{variant}.tsv"
            if not edge_path.exists():
                control_summaries[variant] = {"status": "missing", "path": str(edge_path)}
                continue
            raw_edges = pd.read_csv(edge_path, sep="\t")
            ctrl_edges = signed_causal_edges(raw_edges, plan)
            ctrl_nodes = sorted(set(ctrl_edges["source"].astype(str)) | set(ctrl_edges["target"].astype(str)))
            out_path = ctrl_out / f"ranking_primary_seeds_{variant}.tsv"
            ctrl_rank = _score_ranking_on_edges(
                variant=f"control_{variant}",
                targets=targets,
                drug_universe=drug_universe,
                nodes=ctrl_nodes,
                edges=ctrl_edges,
                module_genes=seed_sets["primary"],
                disease_direction=disease_direction,
                plan=plan,
                out_path=out_path,
                ignore_sign=False,
            )
            control_summaries[variant] = {
                "ranking_path": str(out_path),
                "n_scored": int(len(ctrl_rank)),
                "n_edges": int(len(ctrl_edges)),
                "top20_jaccard_vs_real": _top_k_jaccard(primary_ranking, ctrl_rank, 20),
                "spearman_vs_real": float(
                    primary_ranking.merge(ctrl_rank, on="chembl_id", suffixes=("_real", "_ctrl"))
                    .pipe(lambda d: d["score_real"].corr(d["score_ctrl"], method="spearman"))
                )
                if len(ctrl_rank)
                else float("nan"),
            }

    h7 = None
    if run_h7:
        h7 = run_h7_futility(
            nodes=nodes_primary,
            signed_edges=signed_edges,
            targets=targets,
            scored_drug_ids=primary_ranking["chembl_id"].astype(str).tolist(),
            output_dir=target_out,
            plan=plan,
        )

    summary = {
        "ranking_primary": str(primary_path),
        "n_scored": int(len(primary_ranking)),
        "module": "primary_seeds",
        "seed_annotation": seed_annotation,
        "rankings": {k: str(target_out / f"ranking_{k}.tsv") for k in rankings},
        "soft_positives": soft_reports,
        "axis_disagreement": axis_disagreement,
        "control_graphs": control_summaries,
        "h7_futility": h7,
        "substrate_connectivity": primary.get("connectivity"),
        "present_ranked_list_as_finding": bool(h7 and not h7.get("futile", True)),
        "outstanding_mandatory_defences": [
            "GSE238208_inflammation_specificity",
            "permuted_signature_null",
        ],
        "note": (
            "TARGET remains unsupervised / hypothesis-generating. If H7 fails, report as "
            "non-informative — stability stats and disagreements only; no ranked candidate list "
            "as a finding. Soft positives reported at their ranks, never forced."
        ),
    }
    (target_out / "target_eval_summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    return summary


def prepare_drug_targets(raw_root: str | Path, drug_universe: pd.DataFrame) -> pd.DataFrame:
    raw_root = Path(raw_root)
    chembl = pd.read_csv(raw_root / "chembl" / "mechanisms_merged.tsv", sep="\t")
    dgidb_path = raw_root / "dgidb" / "interactions_by_name.tsv"
    dgidb = pd.read_csv(dgidb_path, sep="\t") if dgidb_path.exists() else pd.DataFrame()
    return build_drug_target_table(chembl, dgidb, drug_universe)
