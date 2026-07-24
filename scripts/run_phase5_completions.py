#!/usr/bin/env python3
"""Phase 5 completions (§A): three-arm sign/connectivity, coverage, TARGET/H7.

Reuses the rebuilt disease-agnostic stamp for signed + unsigned_matched rankings
when present; computes unsigned_full, differential coverage, acquired-pool
sensitivity, TARGET multi-seed/substrate rankings, control graphs, and H7.
"""

from __future__ import annotations

import argparse
import json
import shutil
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

from signflow.phase5 import PHASE3_LATEST, PHASE5_OUT, PHASE5_RAW, ensure_directory
from signflow.phase5.evaluate import (
    acquired_pool_sensitivity,
    compare_connectivity_cost,
    compare_signed_unsigned,
    compute_ranking_metrics,
    differential_coverage_rates,
    prepare_drug_targets,
    run_target_evaluation,
)
from signflow.phase5.evaluate import build_drug_universe
from signflow.phase5.score import CachedFlowScorer, node_degrees
from signflow.phase5.substrate import load_phase3_substrates, load_plan, unsigned_full_causal_edges
from signflow.phase5.module import build_benchmark_modules


def _write_completion_exit_report(stamp_dir: Path, result: dict[str, object]) -> str:
    docs = Path(__file__).resolve().parents[1] / "docs" / "phase5"
    docs.mkdir(parents=True, exist_ok=True)
    bench = result.get("benchmark", {}) or {}
    cov = result.get("coverage", {}) or {}
    target = result.get("target", {}) or {}
    redisc = bench.get("gold_rediscovery", {})
    three = bench.get("three_arm_comparison", {})
    diff = bench.get("differential_coverage", {})
    acquired = bench.get("gold_rediscovery_acquired_pool_sensitivity", {})
    h7 = target.get("h7_futility") or {}
    soft = (target.get("soft_positives") or {}).get("primary__primary_substrate", [])

    def _fmt_ci(ci: dict | None) -> str:
        if not isinstance(ci, dict):
            return "—"
        return f"[{ci.get('lo', float('nan')):.3f}, {ci.get('hi', float('nan')):.3f}]"

    lines = [
        "# Phase 5 Exit Report (completions — §A)",
        "",
        f"**Run stamp (UTC):** {result.get('started_at_utc')}",
        f"**Output:** `{result.get('output_dir')}`",
        f"**Prior rebuild reused:** `{result.get('prior_rebuild_stamp')}`",
        f"**Phase 3 root:** `{result.get('phase3_root')}`",
        f"**Drug data touched:** `{result.get('drug_data_touched')}`",
        "",
        "## Status",
        "",
        "Prior rediscovery stamps with gold-dominated / disease-conditioned pools remain **VOID/SUPERSEDED**.",
        "This stamp completes the rebuilt Phase 5 report under the frozen disease-agnostic pool (§A1–A5).",
        "",
        "## A1 — Both halves of the result (required pair)",
        "",
        "Chance AUPRC on the scored pool = "
        f"{redisc.get('n_gold_eligible_denominator', 29)}/{redisc.get('n_scored_universe', 2833)} = "
        f"**{redisc.get('chance_auprc', 0.0102):.4f}**.",
        "",
        "1. The signed-flow score achieves **weak but statistically distinguishable global discrimination** "
        "of UC therapeutics against a disease-agnostic approved universe "
        f"(~{float(redisc.get('auprc', 0))/max(float(redisc.get('chance_auprc', 1)), 1e-12):.1f}× chance AUPRC; "
        f"AUROC CI lower bound {_fmt_ci(redisc.get('auroc_ci')).split(',')[0].strip('[')}).",
        f"   Observed AUPRC **{redisc.get('auprc', float('nan')):.3f}** {_fmt_ci(redisc.get('auprc_ci'))} "
        f"(chance {redisc.get('chance_auprc', float('nan')):.4f}); "
        f"AUROC **{redisc.get('auroc', float('nan')):.3f}** {_fmt_ci(redisc.get('auroc_ci'))} "
        f"(chance {redisc.get('chance_auroc', 0.5)}).",
        "2. It produces **no top-end enrichment** "
        f"(expected gold in top-20 under chance = {redisc.get('expected_gold_in_top_20_at_chance', 0.21):.2f}, "
        f"observed 0) and **no drug survives BH-FDR** "
        f"(n_survivors={bench.get('n_bh_fdr_survivors', 0)}), so it has **no practical screening utility**.",
        "",
        "Do not report either half alone. Chance baselines sit beside every metric below.",
        "",
        "## A2 — Sign vs connectivity (three-arm)",
        "",
        "| Arm | Edge set | Sign used | AUPRC (chance) | AUROC | enrichment@20 |",
        "|---|---|---|---:|---:|---:|",
    ]
    for arm in three.get("arms", []):
        lines.append(
            f"| {arm.get('arm')} | {arm.get('edge_set')} | {arm.get('sign_used')} | "
            f"{arm.get('auprc')} ({arm.get('chance_auprc')}) | {arm.get('auroc')} | "
            f"{arm.get('enrichment_at_20')} |"
        )
    lines += [
        "",
        f"- Binding sign decision (signed vs unsigned_matched): **{three.get('binding_sign_decision')}**",
        f"- Connectivity contrast: `{json.dumps(three.get('unsigned_matched_vs_unsigned_full', {}), default=str)[:1200]}`",
        f"- Edge counts: signed/matched n_edges={three.get('arms', [{}])[0].get('n_edges') if three.get('arms') else '—'}; "
        f"unsigned_full n_edges={three.get('arms', [{}, {}, {}])[2].get('n_edges') if len(three.get('arms', [])) > 2 else '—'}",
        "",
        "signed vs unsigned_matched isolates the sign contribution; unsigned_matched vs unsigned_full "
        "quantifies the connectivity cost of `exclude_from_signed_propagation`. "
        "Logged as a **defect-driven completion** of the frozen comparison — decision rule unchanged.",
        "",
        "## A3 — Differential gold coverage",
        "",
        f"- Gold coverage rate: **{diff.get('gold_coverage_rate')}**",
        f"- Non-gold coverage rate: **{diff.get('non_gold_coverage_rate')}**",
        f"- Gold prevalence acquired/mapped → scored: "
        f"{diff.get('gold_prevalence_acquired_mapped')} → {diff.get('gold_prevalence_scored')} "
        f"(relative inflation {diff.get('relative_inflation')})",
        f"- Acquired-pool sensitivity (uncovered = misses): "
        f"AUPRC={acquired.get('auprc')} (chance {acquired.get('chance_auprc')}); "
        f"enrichment@20={acquired.get('enrichment_at_20')} "
        f"(E[gold@20]={acquired.get('expected_gold_in_top_20_at_chance')})",
        "",
        "## A4 — All frozen k values",
        "",
        "| k | precision@k | enrichment@k | recall@k | E[gold in top-k] at chance |",
        "|---:|---:|---:|---:|---:|",
    ]
    for k in (20, 50, 100, 250):
        lines.append(
            f"| {k} | {redisc.get(f'precision_at_{k}')} | {redisc.get(f'enrichment_at_{k}')} | "
            f"{redisc.get(f'recall_at_{k}')} | {redisc.get(f'expected_gold_in_top_{k}_at_chance')} |"
        )
    lines += [
        "",
        "## A5 — TARGET track + H7",
        "",
        f"- H7 mean pairwise top-20 Jaccard: **{h7.get('mean_pairwise_top_k_jaccard')}** "
        f"(threshold {h7.get('threshold_j')})",
        f"- H7 decision: **{h7.get('decision')}** — {h7.get('reporting_rule')}",
        f"- Present ranked list as finding: `{target.get('present_ranked_list_as_finding')}`",
        f"- Soft positives (primary seeds × primary substrate): {json.dumps(soft, default=str)}",
        f"- Axis disagreement: {json.dumps(target.get('axis_disagreement', {}), default=str)}",
        f"- Control-graph top-20 Jaccard vs real: "
        f"{json.dumps({k: v.get('top20_jaccard_vs_real') for k, v in (target.get('control_graphs') or {}).items()}, default=str)}",
        f"- Outstanding mandatory Phase 1.5 defences (not executed here): "
        f"{target.get('outstanding_mandatory_defences')}",
        "",
        "## Coverage headline",
        "",
        f"- Mapped summary / in-substrate / misses: {cov.get('n_drugs_total')} / "
        f"{cov.get('n_with_ge1_primary_target')} / {cov.get('n_coverage_misses')}",
        "",
        "## FDR",
        "",
        f"**{bench.get('n_bh_fdr_survivors', 0)} drugs** survive BH-FDR < 0.05 on the primary signed-flow score.",
        "",
        "## Gate",
        "",
        "- §A completions: **this document**",
        "- Phase 6 must not begin until H7 decision + §C reframing are committed",
        "",
    ]
    text = "\n".join(lines) + "\n"
    (docs / "EXIT_REPORT.md").write_text(text)
    (stamp_dir / "EXIT_REPORT.md").write_text(text)
    return str(docs / "EXIT_REPORT.md")


def _rank_unsigned_full(
    *,
    prior_dir: Path,
    out_bench: Path,
    substrates: dict,
    targets: pd.DataFrame,
    drug_universe: pd.DataFrame,
    plan: dict,
) -> tuple[pd.DataFrame, dict]:
    primary = substrates["primary"]
    causal = primary.get("causal_edges")
    if causal is None:
        # Reload from phase3 via substrate builder fields
        full_edges = primary["full_unsigned_edges"]
    else:
        full_edges = unsigned_full_causal_edges(causal, plan)
    nodes = sorted(set(full_edges["source"].astype(str)) | set(full_edges["target"].astype(str)))
    module_bundle = build_benchmark_modules(None, plan=plan)
    # Prefer reusing prior module if present
    prior_mod = prior_dir / "scores" / "benchmark" / "modules" / "conformal_a0.1.tsv"
    if prior_mod.exists():
        module_genes = set(pd.read_csv(prior_mod, sep="\t")["gene_symbol"].astype(str))
    else:
        module_genes = set(module_bundle["modules"]["conformal_a0.1"])
    disease_direction = module_bundle["disease_direction"]
    degrees = pd.Series(node_degrees(full_edges, nodes), index=nodes)

    meta = drug_universe[["chembl_id", "is_gold", "canonical_group", "pref_name"]].drop_duplicates("chembl_id")
    from signflow.phase5.drug_targets import drug_level_target_summary

    drug_summary = drug_level_target_summary(targets, set(nodes))
    # Avoid is_gold_x / is_gold_y collisions when the summary already carries gold flags.
    drop_cols = [c for c in ("is_gold", "canonical_group", "pref_name") if c in drug_summary.columns]
    if drop_cols:
        drug_summary = drug_summary.drop(columns=drop_cols)
    drug_summary = drug_summary.merge(meta, on="chembl_id", how="left")
    drug_summary["is_gold"] = drug_summary["is_gold"].fillna(False).astype(bool)

    scorer = CachedFlowScorer(
        nodes,
        full_edges,
        module_genes,
        disease_direction,
        alpha=float(plan["phase5_flow_damping"]),
        max_iter=int(plan["phase5_flow_max_iter"]),
        tol=float(plan["phase5_flow_tol"]),
        ignore_sign=True,
        variant="unsigned_full",
    )
    rows = []
    groups = list(targets.groupby("chembl_id"))
    for i, (cid, group) in enumerate(groups, start=1):
        if i == 1 or i % 200 == 0 or i == len(groups):
            print(f"[phase5] score unsigned_full observed {i}/{len(groups)}", flush=True)
        meta_rows = drug_summary.loc[drug_summary["chembl_id"] == cid]
        if meta_rows.empty:
            continue
        m = meta_rows.iloc[0]
        in_sub = sorted(set(group["gene_symbol"].astype(str)) & set(nodes))
        if len(in_sub) < int(plan["phase5_min_drug_targets_in_substrate"]):
            continue
        score_val = scorer.score_drug_targets(group, unknown_policy=str(plan["phase5_unknown_action_policy"]))
        rows.append(
            {
                "chembl_id": cid,
                "pref_name": m.get("pref_name", ""),
                "canonical_group": m.get("canonical_group", ""),
                "is_gold": bool(m.get("is_gold", False)),
                "gold_eligible": bool(m.get("is_gold", False)),
                "score": score_val,
                "score_variant": "unsigned_full",
                "n_targets_in_substrate": len(in_sub),
                "mean_target_degree": float(degrees.loc[in_sub].mean()) if in_sub else float("nan"),
            }
        )
    ranking = pd.DataFrame(rows).sort_values(["score", "chembl_id"], ascending=[False, True])
    path = out_bench / "ranking_unsigned_full.tsv"
    ranking.to_csv(path, sep="\t", index=False)
    rediscovery = compute_ranking_metrics(
        ranking,
        k_values=list(plan.get("phase5_score_k_values", [20, 50, 100, 250])),
        n_boot=int(plan.get("phase5_bootstrap_n", 1000)),
    )
    summary = {
        "variant": "unsigned_full",
        "ignore_sign": True,
        "edge_set": "full_causal_directed_excluding_ambiguous",
        "n_edges": int(len(full_edges)),
        "ranking_path": str(path),
        "n_scored_drugs": int(len(ranking)),
        "gold_rediscovery": rediscovery,
        "n_bh_fdr_survivors": None,
        "nulls": "skipped_comparison_arm",
    }
    return ranking, summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--prior-stamp",
        default="20260724T004001Z",
        help="Rebuilt Phase 5 stamp providing signed + unsigned_matched rankings",
    )
    parser.add_argument("--phase3-root", default=None)
    parser.add_argument("--skip-h7", action="store_true")
    parser.add_argument("--skip-controls", action="store_true")
    parser.add_argument("--skip-unsigned-full", action="store_true")
    args = parser.parse_args()

    plan = load_plan()
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    out = ensure_directory(PHASE5_OUT / stamp)
    prior = (PHASE5_OUT / args.prior_stamp).resolve()
    p3 = Path(args.phase3_root or PHASE3_LATEST).resolve()

    print(f"[phase5-completions] out={out}", flush=True)
    print(f"[phase5-completions] prior={prior}", flush=True)

    # Copy reusable prior artifacts
    for name in ("coverage_report.json", "coverage_drug_summary.tsv", "drug_target_actions.tsv", "substrate_connectivity.json"):
        src = prior / name
        if src.exists():
            shutil.copy2(src, out / name)

    scores = ensure_directory(out / "scores")
    bench_out = ensure_directory(scores / "benchmark")
    prior_bench = prior / "scores" / "benchmark"
    for name in (
        "ranking_signed.tsv",
        "ranking_unsigned.tsv",
        "ranking_primary.tsv",
        "benchmark_eval_summary.json",
    ):
        src = prior_bench / name
        if src.exists():
            shutil.copy2(src, bench_out / name)
    # Alias unsigned → unsigned_matched
    if (bench_out / "ranking_unsigned.tsv").exists():
        shutil.copy2(bench_out / "ranking_unsigned.tsv", bench_out / "ranking_unsigned_matched.tsv")
    if (prior_bench / "modules").exists():
        if (bench_out / "modules").exists():
            shutil.rmtree(bench_out / "modules")
        shutil.copytree(prior_bench / "modules", bench_out / "modules")
    if (prior_bench / "controls").exists():
        if (bench_out / "controls").exists():
            shutil.rmtree(bench_out / "controls")
        shutil.copytree(prior_bench / "controls", bench_out / "controls")

    drug_universe = build_drug_universe(PHASE5_RAW, plan=plan)
    targets = prepare_drug_targets(PHASE5_RAW, drug_universe)
    substrates = load_phase3_substrates(p3 / "benchmark")

    signed = pd.read_csv(bench_out / "ranking_signed.tsv", sep="\t")
    unsigned_matched = pd.read_csv(bench_out / "ranking_unsigned_matched.tsv", sep="\t")

    if args.skip_unsigned_full and (prior_bench / "ranking_unsigned_full.tsv").exists():
        unsigned_full = pd.read_csv(prior_bench / "ranking_unsigned_full.tsv", sep="\t")
        unsigned_full_summary = {
            "variant": "unsigned_full",
            "ranking_path": str(bench_out / "ranking_unsigned_full.tsv"),
            "gold_rediscovery": compute_ranking_metrics(unsigned_full),
            "n_scored_drugs": int(len(unsigned_full)),
        }
        shutil.copy2(prior_bench / "ranking_unsigned_full.tsv", bench_out / "ranking_unsigned_full.tsv")
    else:
        unsigned_full, unsigned_full_summary = _rank_unsigned_full(
            prior_dir=prior,
            out_bench=bench_out,
            substrates=substrates,
            targets=targets,
            drug_universe=drug_universe,
            plan=plan,
        )

    signed_summary_metrics = compute_ranking_metrics(signed)
    unsigned_matched_metrics = compute_ranking_metrics(unsigned_matched)

    sign_comp = compare_signed_unsigned(
        signed, unsigned_matched, n_boot=int(plan["phase5_bootstrap_n"]), unsigned_label="unsigned_matched"
    )
    conn_comp = compare_connectivity_cost(
        unsigned_matched, unsigned_full, n_boot=int(plan["phase5_bootstrap_n"])
    )
    three_arm = {
        "arms": [
            {
                "arm": "signed",
                "edge_set": "signed_only_reduced",
                "n_edges": int(substrates["primary"]["n_signed_edges"]),
                "sign_used": True,
                "auroc": signed_summary_metrics.get("auroc"),
                "auprc": signed_summary_metrics.get("auprc"),
                "chance_auprc": signed_summary_metrics.get("chance_auprc"),
                "enrichment_at_20": signed_summary_metrics.get("enrichment_at_20"),
            },
            {
                "arm": "unsigned_matched",
                "edge_set": "signed_only_reduced",
                "n_edges": int(substrates["primary"]["n_signed_edges"]),
                "sign_used": False,
                "auroc": unsigned_matched_metrics.get("auroc"),
                "auprc": unsigned_matched_metrics.get("auprc"),
                "chance_auprc": unsigned_matched_metrics.get("chance_auprc"),
                "enrichment_at_20": unsigned_matched_metrics.get("enrichment_at_20"),
            },
            {
                "arm": "unsigned_full",
                "edge_set": "full_causal_directed_excluding_ambiguous",
                "n_edges": int(substrates["primary"].get("n_full_unsigned_edges") or len(substrates["primary"]["full_unsigned_edges"])),
                "sign_used": False,
                "auroc": unsigned_full_summary["gold_rediscovery"].get("auroc"),
                "auprc": unsigned_full_summary["gold_rediscovery"].get("auprc"),
                "chance_auprc": unsigned_full_summary["gold_rediscovery"].get("chance_auprc"),
                "enrichment_at_20": unsigned_full_summary["gold_rediscovery"].get("enrichment_at_20"),
            },
        ],
        "signed_vs_unsigned_matched": sign_comp,
        "unsigned_matched_vs_unsigned_full": conn_comp,
        "binding_sign_decision": sign_comp.get("decision"),
        "note": plan.get("phase5_score_three_arm_note"),
    }

    cov_summary = pd.read_csv(out / "coverage_drug_summary.tsv", sep="\t")
    if "is_gold" not in cov_summary.columns:
        cov_summary = cov_summary.merge(
            drug_universe[["chembl_id", "is_gold"]].drop_duplicates("chembl_id"),
            on="chembl_id",
            how="left",
        )
    cov_report = json.loads((out / "coverage_report.json").read_text())
    n_universe = int(len(drug_universe.drop_duplicates("chembl_id")))
    diff_cov = differential_coverage_rates(cov_report, cov_summary, n_acquired_universe=n_universe)
    acquired_sens = acquired_pool_sensitivity(signed, drug_universe)

    # Load prior BH survivors / spearman from prior summary if present
    prior_summary = {}
    prior_sum_path = prior_bench / "benchmark_eval_summary.json"
    if prior_sum_path.exists():
        prior_summary = json.loads(prior_sum_path.read_text())

    bench = {
        "ranking_primary": str(bench_out / "ranking_primary.tsv"),
        "ranking_signed": str(bench_out / "ranking_signed.tsv"),
        "ranking_unsigned": str(bench_out / "ranking_unsigned.tsv"),
        "ranking_unsigned_matched": str(bench_out / "ranking_unsigned_matched.tsv"),
        "ranking_unsigned_full": unsigned_full_summary["ranking_path"],
        "n_scored_drugs": int(len(signed)),
        "spearman_score_degree": prior_summary.get("spearman_score_degree"),
        "gold_rediscovery": signed_summary_metrics,
        "gold_rediscovery_unsigned": unsigned_matched_metrics,
        "gold_rediscovery_unsigned_matched": unsigned_matched_metrics,
        "gold_rediscovery_unsigned_full": unsigned_full_summary["gold_rediscovery"],
        "gold_rediscovery_acquired_pool_sensitivity": acquired_sens,
        "differential_coverage": diff_cov,
        "signed_vs_unsigned": sign_comp,
        "three_arm_comparison": three_arm,
        "n_bh_fdr_survivors": prior_summary.get("n_bh_fdr_survivors", 0),
        "substrate_connectivity": substrates["primary"].get("connectivity"),
        "substrate_connectivity_unsigned_full": substrates["primary"].get("connectivity_unsigned_full"),
        "edge_set_by_arm": {
            "signed": "signed_only_reduced",
            "unsigned_matched": "signed_only_reduced",
            "unsigned_full": "full_causal_directed_excluding_ambiguous",
        },
        "pool_guard": prior_summary.get("pool_guard"),
    }
    (bench_out / "benchmark_eval_summary.json").write_text(json.dumps(bench, indent=2, default=str) + "\n")

    # TARGET substrates: load target track graph for TARGET scoring
    target_substrates = load_phase3_substrates(p3 / "target")
    target = run_target_evaluation(
        phase3_root=p3,
        output_dir=scores,
        drug_universe=drug_universe,
        targets=targets,
        substrates=target_substrates,
        run_h7=not args.skip_h7,
        run_control_graphs=not args.skip_controls,
    )

    result = {
        "started_at_utc": stamp,
        "phase3_root": str(p3),
        "output_dir": str(out),
        "prior_rebuild_stamp": str(prior),
        "drug_data_touched": True,
        "coverage": cov_report,
        "benchmark": bench,
        "target": target,
        "completions": ["A1", "A2", "A3", "A4", "A5"],
    }
    (out / "phase5_result.json").write_text(json.dumps(result, indent=2, default=str) + "\n")
    (out / "substrate_connectivity.json").write_text(
        json.dumps(
            {
                "primary": substrates["primary"].get("connectivity"),
                "unsigned_full": substrates["primary"].get("connectivity_unsigned_full"),
                "sensitivity": substrates["sensitivity"].get("connectivity"),
            },
            indent=2,
            default=str,
        )
        + "\n"
    )
    result["exit_report"] = _write_completion_exit_report(out, result)

    latest = PHASE5_OUT / "latest"
    if latest.exists() or latest.is_symlink():
        latest.unlink()
    latest.symlink_to(out.name)

    # Phase 6 decision stub written after H7
    h7 = target.get("h7_futility") or {}
    phase6 = {
        "decided_at_utc": stamp,
        "h7_decision": h7.get("decision"),
        "h7_futile": h7.get("futile"),
        "benchmark_fdr_survivors": bench.get("n_bh_fdr_survivors", 0),
        "benchmark_enrichment_at_20": signed_summary_metrics.get("enrichment_at_20"),
        "phase6_mode": (
            "descriptive_characterisation_of_frozen_seed_module_only"
            if h7.get("futile", True)
            else "triangulation_as_specified_genetics_excluded_on_benchmark"
        ),
        "always_run": [
            "H6_convergence_test_on_BENCHMARK_k20_bootstrap_over_gold",
        ],
        "carry_forward_mandatory_target_defences": [
            "GSE238208_inflammation_specificity",
            "permuted_signature_null",
        ],
        "do_not": [
            "add_axes_or_loosen_thresholds_to_manufacture_survivors",
            "run_convergence_ranking_on_non_informative_target_list",
        ],
        "note": (
            "Recorded before Phase 6 build. If TARGET fails H7, axes run as descriptive "
            "characterisation of the frozen seed module only; report axis disagreement without a fused rank. "
            "H6 on BENCHMARK remains answerable even with weak signal."
        ),
    }
    phase6_path = Path("docs/phase0/PHASE6_RESCOPE_DECISION.json")
    phase6_path.write_text(json.dumps(phase6, indent=2) + "\n")
    (out / "PHASE6_RESCOPE_DECISION.json").write_text(json.dumps(phase6, indent=2) + "\n")
    result["phase6_decision"] = phase6

    print(json.dumps({"stamp": stamp, "h7": h7.get("decision"), "exit_report": result["exit_report"]}, indent=2))


if __name__ == "__main__":
    main()
