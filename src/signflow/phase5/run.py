"""Phase 5 orchestrator."""

from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

from . import FROZEN_PLAN, PHASE3_LATEST, PHASE5_OUT, PHASE5_RAW, StopAndReportError, ensure_directory
from .acquire import acquire_drug_resources, load_approved_drug_universe, verify_benchmark_drugs_sha
from .coverage import write_coverage_report
from .evaluate import (
    build_drug_universe,
    prepare_drug_targets,
    run_benchmark_evaluation,
    run_target_evaluation,
)
from .substrate import load_phase3_substrates, load_plan


def _mark_drug_data_touched(result: dict[str, object], *, resources_read: bool) -> None:
    """Provenance flag must be true whenever any drug resource is read."""

    if resources_read:
        result["drug_data_touched"] = True


def assert_drug_data_touched(result: dict[str, object], *, resources_read: bool) -> None:
    """Stop-and-report if the provenance flag inverts reality."""

    if resources_read and result.get("drug_data_touched") is not True:
        raise StopAndReportError(
            "drug_data_touched is false while drug resources were read",
            reason="drug-data-touched-false",
            details={"drug_data_touched": result.get("drug_data_touched")},
        )


def _write_exit_report(stamp_dir: Path, result: dict[str, object]) -> str:
    docs = Path(__file__).resolve().parents[3] / "docs" / "phase5"
    docs.mkdir(parents=True, exist_ok=True)
    cov = result.get("coverage", {})
    bench = result.get("benchmark", {}) or {}
    redisc = bench.get("gold_rediscovery", {})
    unsigned = bench.get("gold_rediscovery_unsigned", {})
    comparison = bench.get("signed_vs_unsigned", {})
    conn = bench.get("substrate_connectivity", {})
    pool = bench.get("pool_guard", {})
    lines = [
        "# Phase 5 Exit Report",
        "",
        f"**Run stamp (UTC):** {result.get('started_at_utc')}",
        f"**Output:** `{result.get('output_dir')}`",
        f"**Phase 3 root:** `{result.get('phase3_root')}`",
        f"**Drug data touched:** `{result.get('drug_data_touched')}`",
        f"**Pool SHA-256:** `{result.get('pool_sha256')}`",
        "",
        "## Status",
        "",
        "- Prior rediscovery stamps with gold-dominated / disease-conditioned pools are **VOID/SUPERSEDED** "
        "(see `VOID.md` under those run directories).",
        "- Ranking universe = disease-agnostic approved drugs (ChEMBL `max_phase=4`) ∪ frozen UC gold labels.",
        "- Signed is the pre-registered primary; unsigned is co-primary. Signed is never dropped if it underperforms.",
        "- No clinical recommendations.",
        "",
        "## Pool construction",
        "",
        f"- Rule: disease-agnostic approved drugs; disease conditioning forbidden.",
        f"- Universe size / gold / prevalence: {pool}",
        f"- In-substrate floor / gold-prevalence ceiling enforced (stop-and-report).",
        "",
        "## Coverage (headline)",
        "",
        f"- Drugs total (mapped summary): {cov.get('n_drugs_total')}",
        f"- With ≥1 primary-substrate target: {cov.get('n_with_ge1_primary_target')}",
        f"- Gold eligible / misses: {cov.get('n_gold_eligible_for_rediscovery')} / {cov.get('n_gold_coverage_misses')}",
        f"- Target-class breakdown (must include ion_channel): {cov.get('target_class_breakdown_primary')}",
        f"- Gold class breakdown: {cov.get('target_class_breakdown_gold_primary')}",
        "",
        "## Signed-propagation substrate connectivity",
        "",
        json.dumps(conn, indent=2, default=str),
        "",
        "## Signed vs unsigned (co-primary)",
        "",
        json.dumps(comparison, indent=2, default=str),
        "",
        "## Scoring (signed primary)",
        "",
        f"- Scored drugs: {bench.get('n_scored_drugs')}",
        f"- Spearman(score, mean target degree): {bench.get('spearman_score_degree')}",
        f"- BH-FDR < 0.05 survivors: {bench.get('n_bh_fdr_survivors', 0)}",
        "",
        "## Gold rediscovery (signed)",
        "",
        json.dumps({k: v for k, v in redisc.items() if k != "gold_ranks"}, indent=2, default=str)[:4000],
        "",
        "## Gold rediscovery (unsigned)",
        "",
        json.dumps({k: v for k, v in unsigned.items() if k != "gold_ranks"}, indent=2, default=str)[:4000],
        "",
        "## FDR statement",
        "",
        f"**{bench.get('n_bh_fdr_survivors', 0)} drugs** survive BH-FDR < 0.05 on the primary signed-flow score "
        "under the degree-matched null — including if none does.",
        "",
        "## TARGET",
        "",
        json.dumps(result.get("target", {}), indent=2, default=str)[:2000],
        "",
    ]
    path = docs / "EXIT_REPORT.md"
    path.write_text("\n".join(lines) + "\n")
    (stamp_dir / "EXIT_REPORT.md").write_text(path.read_text())
    return str(path)


def run_phase5(
    *,
    phase3_root: str | Path | None = None,
    raw_root: str | Path | None = None,
    skip_acquire: bool = False,
) -> dict[str, object]:
    plan = load_plan()
    verify_benchmark_drugs_sha()
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    out = ensure_directory(PHASE5_OUT / stamp)
    p3 = Path(phase3_root or PHASE3_LATEST).resolve()
    raw = Path(raw_root or PHASE5_RAW)

    substrates = load_phase3_substrates(p3 / "benchmark")
    primary_nodes = substrates["primary"]["nodes"]
    sensitivity_nodes = substrates["sensitivity"]["nodes"]

    drug_data_touched = False
    if skip_acquire and (raw / "pool" / "approved_drug_universe.tsv").exists():
        drug_universe = build_drug_universe(raw, plan=plan)
        drug_data_touched = True
        acquire_info = {
            "raw_root": str(raw),
            "manifests": {},
            "skipped": True,
            "drug_data_touched": True,
            "note": "Acquisition skipped but approved-drug universe + drug tables were read.",
        }
    else:
        # build_drug_universe may create/load pool; acquire fills mechanisms
        if (raw / "pool" / "approved_drug_universe.tsv").exists():
            drug_universe = build_drug_universe(raw, plan=plan)
            drug_data_touched = True
        else:
            drug_universe = load_approved_drug_universe(raw)  # triggers fetch
            drug_data_touched = True
        acquire_info = acquire_drug_resources(drug_universe, raw_root=raw)
        drug_data_touched = True
        drug_universe = build_drug_universe(raw, plan=plan)

    targets = prepare_drug_targets(raw, drug_universe)
    drug_data_touched = True
    targets.to_csv(out / "drug_target_actions.tsv", sep="\t", index=False)

    from .drug_targets import drug_level_target_summary

    drug_summary = drug_level_target_summary(targets, primary_nodes)
    if "is_gold" not in drug_summary.columns:
        drug_summary = drug_summary.merge(
            drug_universe[["chembl_id", "is_gold", "canonical_group", "pref_name"]],
            on="chembl_id",
            how="left",
        )

    cov = write_coverage_report(
        drug_summary,
        targets,
        primary_nodes=primary_nodes,
        sensitivity_nodes=sensitivity_nodes,
        output_dir=out,
        plan=plan,
    )

    scores_root = ensure_directory(out / "scores")
    bench = run_benchmark_evaluation(
        phase3_root=p3,
        output_dir=scores_root,
        drug_universe=drug_universe,
        targets=targets,
        substrates=substrates,
        acquire_info=acquire_info,
    )
    target = run_target_evaluation(
        phase3_root=p3,
        output_dir=scores_root,
        drug_universe=drug_universe,
        targets=targets,
        substrates=substrates,
    )

    pool_path = raw / "pool" / "approved_drug_universe.tsv"
    from . import sha256_path

    result: dict[str, object] = {
        "started_at_utc": stamp,
        "phase3_root": str(p3),
        "output_dir": str(out),
        "drug_data_touched": drug_data_touched,
        "pool_sha256": sha256_path(pool_path) if pool_path.exists() else None,
        "benchmark_drugs_sha": verify_benchmark_drugs_sha(),
        "coverage": cov["report"],
        "benchmark": bench,
        "target": target,
        "acquisition": acquire_info,
        "voided_prior_stamp": plan.get("phase5_voided_rediscovery_stamp"),
    }
    _mark_drug_data_touched(result, resources_read=True)
    assert_drug_data_touched(result, resources_read=True)

    (out / "phase5_result.json").write_text(json.dumps(result, indent=2, default=str) + "\n")
    (out / "substrate_connectivity.json").write_text(
        json.dumps(
            {
                "primary": substrates["primary"].get("connectivity"),
                "sensitivity": substrates["sensitivity"].get("connectivity"),
            },
            indent=2,
            default=str,
        )
        + "\n"
    )
    latest = PHASE5_OUT / "latest"
    if latest.exists() or latest.is_symlink():
        latest.unlink()
    latest.symlink_to(out.name)
    result["exit_report"] = _write_exit_report(out, result)
    return result
