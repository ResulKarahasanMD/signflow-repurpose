"""Score drugs, nulls, rankings, and gold rediscovery metrics."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
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
    for chembl_id, group in targets.groupby("chembl_id"):
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

    if len(ranking):
        ranking = attach_nulls_and_fdr(
            ranking,
            degrees=degrees,
            score_fn_for_genes=score_fn,
            n_perm=int(plan["null_permutations"]),
            tolerance=float(plan["degree_match_tolerance"]),
            seed=0,
        )
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
    n_bh = int((ranking["bh_fdr"] < 0.05).sum()) if len(ranking) and "bh_fdr" in ranking.columns else 0
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
) -> dict[str, object]:
    """Paired bootstrap over gold set for enrichment@k and AUROC/AUPRC; decision rule applied."""

    plan = load_plan()
    k_primary = int(plan.get("phase5_score_k_primary", 20))
    merged = signed_ranking[["chembl_id", "is_gold", "score"]].merge(
        unsigned_ranking[["chembl_id", "score"]],
        on="chembl_id",
        suffixes=("_signed", "_unsigned"),
    )
    if merged.empty:
        return {"status": "empty", "decision": "no_measurable_sign_benefit"}
    y = merged["is_gold"].astype(bool).astype(int).to_numpy()
    s_signed = merged["score_signed"].astype(float).to_numpy()
    s_unsigned = merged["score_unsigned"].astype(float).to_numpy()
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
            "Decision reports whether sign contributed measurable benefit."
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
) -> dict[str, object]:
    plan = load_plan()
    output_dir = ensure_directory(output_dir)
    bench_out = ensure_directory(output_dir / "benchmark")

    primary = substrates["primary"]
    nodes = sorted(primary["nodes"])
    signed_edges = primary["signed_edges"]
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
    unsigned_ranking, unsigned_summary = _rank_variant(
        variant="unsigned",
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
    comparison = compare_signed_unsigned(
        signed_ranking,
        unsigned_ranking,
        n_boot=int(plan.get("phase5_bootstrap_n", 1000)),
    )

    ensure_directory(bench_out / "controls")
    negatives = select_unrelated_negative_drugs(drug_universe, plan=plan)
    negatives.to_csv(bench_out / "controls" / "unrelated_negatives.tsv", sep="\t", index=False)

    # Primary ranking alias (signed retained even if unsigned wins)
    primary_path = bench_out / "ranking_primary.tsv"
    signed_ranking.to_csv(primary_path, sep="\t", index=False)

    result = {
        "module_paths": mod_paths,
        "ranking_primary": str(primary_path),
        "ranking_signed": signed_summary["ranking_path"],
        "ranking_unsigned": unsigned_summary["ranking_path"],
        "n_scored_drugs": signed_summary["n_scored_drugs"],
        "spearman_score_degree": signed_summary["spearman_score_degree"],
        "gold_rediscovery": signed_summary["gold_rediscovery"],
        "gold_rediscovery_unsigned": unsigned_summary["gold_rediscovery"],
        "signed_summary": signed_summary,
        "unsigned_summary": unsigned_summary,
        "signed_vs_unsigned": comparison,
        "n_bh_fdr_survivors": signed_summary["n_bh_fdr_survivors"],
        "primary_module": "conformal_a0.1",
        "substrate": primary["name"],
        "substrate_connectivity": primary.get("connectivity"),
        "sensitivity_connectivity": substrates.get("sensitivity", {}).get("connectivity"),
        "pool_guard": pool_guard,
        "acquire_info": acquire_info,
    }
    (bench_out / "benchmark_eval_summary.json").write_text(json.dumps(result, indent=2, default=str) + "\n")
    return result


def run_target_evaluation(
    *,
    phase3_root: str | Path | None,
    output_dir: str | Path,
    drug_universe: pd.DataFrame,
    targets: pd.DataFrame,
    substrates: dict[str, object],
) -> dict[str, object]:
    plan = load_plan()
    output_dir = ensure_directory(output_dir)
    target_out = ensure_directory(output_dir / "target")
    root = Path(phase3_root) if phase3_root else None
    tmods = load_target_modules(root)
    write_module_artifacts(tmods, target_out / "modules")

    primary = substrates["primary"]
    nodes = sorted(primary["nodes"])
    signed_edges = primary["signed_edges"]
    module_genes = set(tmods["modules"]["primary"])
    disease_direction = tmods["disease_direction"]

    rows = []
    for chembl_id, group in targets.groupby("chembl_id"):
        genes_in = set(group["gene_symbol"].astype(str)) & set(nodes)
        if len(genes_in) < int(plan["phase5_min_drug_targets_in_substrate"]):
            continue
        scored = _score_one_drug(
            group,
            nodes=nodes,
            signed_edges=signed_edges,
            module_genes=module_genes,
            disease_direction=disease_direction,
            plan=plan,
            ignore_sign=False,
        )
        meta = drug_universe.loc[drug_universe["chembl_id"] == chembl_id].iloc[0]
        rows.append(
            {
                "chembl_id": chembl_id,
                "pref_name": meta.get("pref_name", ""),
                "score": scored["score"],
                "n_targets_in_substrate": scored["n_targets_in_substrate"],
            }
        )
    ranking = pd.DataFrame(rows).sort_values("score", ascending=False)
    rank_path = target_out / "ranking_primary.tsv"
    ranking.to_csv(rank_path, sep="\t", index=False)
    summary = {
        "ranking_primary": str(rank_path),
        "n_scored": int(len(ranking)),
        "module": "primary_seeds",
        "substrate_connectivity": primary.get("connectivity"),
    }
    (target_out / "target_eval_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def prepare_drug_targets(raw_root: str | Path, drug_universe: pd.DataFrame) -> pd.DataFrame:
    raw_root = Path(raw_root)
    chembl = pd.read_csv(raw_root / "chembl" / "mechanisms_merged.tsv", sep="\t")
    dgidb_path = raw_root / "dgidb" / "interactions_by_name.tsv"
    dgidb = pd.read_csv(dgidb_path, sep="\t") if dgidb_path.exists() else pd.DataFrame()
    return build_drug_target_table(chembl, dgidb, drug_universe)
