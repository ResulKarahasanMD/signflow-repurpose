"""Pre-score drug-target coverage report (primary + sensitivity substrates)."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import yaml

from . import FROZEN_PLAN, StopAndReportError, ensure_directory


def load_plan() -> dict[str, object]:
    return yaml.safe_load(FROZEN_PLAN.read_text())["plan"]


_TAXONOMY = (
    "gpcr",
    "ion_channel",
    "kinase",
    "nuclear_receptor",
    "transporter",
    "other_signalling_enzyme",
)


def _target_class(gene: str) -> str:
    """Heuristic target-class taxonomy; ion_channel is mandatory in the frozen taxonomy."""

    g = str(gene).upper()
    if g.startswith(("GPR", "DRD", "HTR", "ADRB", "CHRM", "OPR", "HRH")):
        return "gpcr"
    # Ion channels: voltage-gated, ligand-gated, TRP, etc.
    if g.startswith(("KCN", "SCN", "CACNA", "CACN", "GRIN", "GRIA", "GABR", "CHRN", "TRP", "P2RX", "CLCN", "HCN")):
        return "ion_channel"
    if g.startswith(("SLC", "ABC")):
        return "transporter"
    if g.endswith("R") and len(g) <= 6:
        return "nuclear_receptor"
    if g.endswith("K"):
        return "kinase"
    return "other_signalling_enzyme"


def _class_breakdown(genes: list[str]) -> dict[str, int]:
    counts = {k: 0 for k in _TAXONOMY}
    for gene in genes:
        counts[_target_class(gene)] = counts.get(_target_class(gene), 0) + 1
    # Unique genes per class
    by_class: dict[str, set[str]] = {k: set() for k in _TAXONOMY}
    for gene in genes:
        by_class[_target_class(gene)].add(str(gene))
    return {k: int(len(v)) for k, v in by_class.items()}


def build_coverage_tables(
    drug_summary: pd.DataFrame,
    targets: pd.DataFrame,
    *,
    primary_nodes: set[str],
    sensitivity_nodes: set[str],
    plan: dict[str, object] | None = None,
) -> dict[str, object]:
    plan = plan or load_plan()
    min_targets = int(plan["phase5_min_drug_targets_in_substrate"])

    summary = drug_summary.copy()
    summary["eligible_primary_score"] = summary["n_targets_in_substrate"] >= min_targets

    if len(targets):
        sens_rows = []
        for chembl_id, group in targets.groupby("chembl_id"):
            genes = set(group["gene_symbol"].astype(str))
            sens_rows.append(
                {"chembl_id": chembl_id, "n_targets_sensitivity_substrate": len(genes & sensitivity_nodes)}
            )
        summary = summary.merge(pd.DataFrame(sens_rows), on="chembl_id", how="left")
    else:
        summary["n_targets_sensitivity_substrate"] = 0

    gold_col = "is_gold"
    if "is_gold" not in summary.columns:
        if "is_gold_y" in summary.columns:
            gold_col = "is_gold_y"
        elif "is_gold_x" in summary.columns:
            gold_col = "is_gold_x"
        else:
            gold_col = None
    gold = summary.loc[summary[gold_col].astype(bool)].copy() if gold_col else summary.iloc[0:0].copy()
    gold_miss = gold.loc[~gold["eligible_primary_score"]]

    tgt = targets.copy()
    if len(tgt):
        tgt["in_primary"] = tgt["gene_symbol"].astype(str).isin(primary_nodes)
        primary_genes = tgt.loc[tgt["in_primary"], "gene_symbol"].astype(str).tolist()
        class_breakdown = _class_breakdown(primary_genes)
        gold_ids = set(gold["chembl_id"].astype(str)) if len(gold) else set()
        gold_genes = (
            tgt.loc[tgt["in_primary"] & tgt["chembl_id"].astype(str).isin(gold_ids), "gene_symbol"]
            .astype(str)
            .tolist()
        )
        gold_class_breakdown = _class_breakdown(gold_genes)
    else:
        class_breakdown = {k: 0 for k in _TAXONOMY}
        gold_class_breakdown = {k: 0 for k in _TAXONOMY}

    # Taxonomy must always expose ion_channel (coverage gap measurement)
    for key in _TAXONOMY:
        class_breakdown.setdefault(key, 0)
        gold_class_breakdown.setdefault(key, 0)

    hist = summary["n_targets_in_substrate"].value_counts().sort_index().to_dict()
    n_gold = int(len(gold))
    n_total = int(len(summary))
    report = {
        "substrate_primary": plan["phase5_proximity_substrate_primary"],
        "substrate_sensitivity": plan["phase5_proximity_substrate_sensitivity"],
        "min_drug_targets_in_substrate": min_targets,
        "poor_coverage_policy": plan["phase5_poor_coverage_policy"],
        "n_drugs_total": n_total,
        "n_with_ge1_primary_target": int(summary["eligible_primary_score"].sum()),
        "n_coverage_misses": int((~summary["eligible_primary_score"]).sum()),
        "n_gold_total": n_gold,
        "gold_prevalence_in_mapped_summary": float(n_gold / max(n_total, 1)),
        "n_gold_eligible_for_rediscovery": int(gold["eligible_primary_score"].sum()) if n_gold else 0,
        "n_gold_coverage_misses": int(len(gold_miss)),
        "gold_out_of_substrate_policy": plan.get(
            "benchmark_gold_out_of_substrate_policy", "exclude_from_denominator_report_as_coverage_miss"
        ),
        "targets_per_drug_histogram_primary": {str(k): int(v) for k, v in hist.items()},
        "target_class_taxonomy": list(_TAXONOMY),
        "target_class_breakdown_primary": class_breakdown,
        "target_class_breakdown_gold_primary": gold_class_breakdown,
        "coverage_miss_chembl_ids": summary.loc[~summary["eligible_primary_score"], "chembl_id"].astype(str).tolist(),
        "gold_coverage_miss_chembl_ids": gold_miss["chembl_id"].astype(str).tolist(),
        "gold_coverage_miss_groups": gold_miss.get("canonical_group", pd.Series(dtype=str)).astype(str).tolist(),
    }
    return {"summary": summary, "report": report}


def write_coverage_report(
    drug_summary: pd.DataFrame,
    targets: pd.DataFrame,
    *,
    primary_nodes: set[str],
    sensitivity_nodes: set[str],
    output_dir: str | Path,
    plan: dict[str, object] | None = None,
) -> dict[str, object]:
    output_dir = ensure_directory(output_dir)
    built = build_coverage_tables(
        drug_summary, targets, primary_nodes=primary_nodes, sensitivity_nodes=sensitivity_nodes, plan=plan
    )
    path = Path(output_dir) / "coverage_report.json"
    path.write_text(json.dumps(built["report"], indent=2, default=str) + "\n")
    built["report"]["path"] = str(path)
    built["summary"].to_csv(Path(output_dir) / "coverage_drug_summary.tsv", sep="\t", index=False)
    return built
