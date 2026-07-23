"""Stratified bootstrap helpers for Phase 2 signatures.

Bootstrap propagates the standardized-effect (Hedges' g) path used by the
primary random-effects meta-analysis. Full limma eBayes is not re-run inside
each resample (prohibitive at n_boot>=1000 via rpy2); limma tables remain the
observed-sample primary DE artifacts. Logged in ANALYSIS_PLAN_AMENDMENTS.md.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from . import StopAndReportError, ensure_directory


def stratified_bootstrap_samples(
    phenotype_table: pd.DataFrame,
    *,
    n_bootstrap: int = 1000,
    seed: int = 0,
) -> dict[str, object]:
    """Sample within case/control strata while preserving group sizes."""

    case_ids = phenotype_table.loc[phenotype_table["group"] == "case", "sample_id"].tolist()
    control_ids = phenotype_table.loc[phenotype_table["group"] == "control", "sample_id"].tolist()
    if not case_ids or not control_ids:
        raise StopAndReportError(
            "Stratified bootstrap requires both case and control samples",
            reason="ambiguous-phenotype",
        )

    rng = np.random.default_rng(seed)
    resamples: list[list[str]] = []
    degenerate = 0
    for _ in range(n_bootstrap):
        boot_case = rng.choice(case_ids, size=len(case_ids), replace=True).tolist()
        boot_control = rng.choice(control_ids, size=len(control_ids), replace=True).tolist()
        if len(set(boot_case)) == 1 or len(set(boot_control)) == 1:
            degenerate += 1
        resamples.append(boot_case + boot_control)
    return {
        "resamples": resamples,
        "degenerate_resamples": degenerate,
        "method": "stratified_preserving_group_sizes",
    }


def _effects_for_resample(
    expression: pd.DataFrame,
    phenotype_table: pd.DataFrame,
    sample_ids: list[str],
    *,
    study: str,
) -> pd.DataFrame:
    from .de_limma import standardized_effects

    groups = []
    selected_data = {}
    for i, sid in enumerate(sample_ids):
        if sid not in expression.columns:
            continue
        group = phenotype_table.loc[phenotype_table["sample_id"] == sid, "group"].iloc[0]
        col = f"{sid}__{i}"
        groups.append({"sample_id": col, "group": group})
        selected_data[col] = expression[sid].to_numpy()
    selected = pd.DataFrame(selected_data, index=expression.index)
    pheno = pd.DataFrame(groups)
    effects = standardized_effects(selected, pheno)
    effects["study"] = study
    return effects


def _dl_two_study_pooled(e1: np.ndarray, s1: np.ndarray, e2: np.ndarray, s2: np.ndarray) -> np.ndarray:
    v1 = s1 ** 2
    v2 = s2 ** 2
    w1 = 1.0 / v1
    w2 = 1.0 / v2
    wsum = w1 + w2
    fixed = (w1 * e1 + w2 * e2) / wsum
    q = w1 * (e1 - fixed) ** 2 + w2 * (e2 - fixed) ** 2
    c = wsum - (w1**2 + w2**2) / wsum
    tau2 = np.maximum((q - 1.0) / np.where(c > 0, c, np.nan), 0.0)
    tau2 = np.nan_to_num(tau2, nan=0.0)
    rw1 = 1.0 / (v1 + tau2)
    rw2 = 1.0 / (v2 + tau2)
    return (rw1 * e1 + rw2 * e2) / (rw1 + rw2)


def bootstrap_meta_signatures(
    study_expression: dict[str, pd.DataFrame],
    study_phenotypes: dict[str, pd.DataFrame],
    *,
    n_bootstrap: int = 1000,
    seed: int = 0,
    output_dir: str | Path,
) -> dict[str, object]:
    """Propagate stratified bootstraps through Hedges' g + DL meta (not limma each time)."""

    output_dir = ensure_directory(output_dir)
    studies = list(study_expression.keys())
    if len(studies) != 2:
        raise StopAndReportError("bootstrap_meta_signatures currently requires exactly two studies", reason="unavailable-files")

    per_study_resamples = {}
    degenerate_total = 0
    for study, pheno in study_phenotypes.items():
        boot = stratified_bootstrap_samples(pheno, n_bootstrap=n_bootstrap, seed=seed)
        per_study_resamples[study] = boot["resamples"]
        degenerate_total += int(boot["degenerate_resamples"])

    shared = sorted(set.intersection(*(set(expr.index) for expr in study_expression.values())))
    study_expression = {k: v.loc[shared] for k, v in study_expression.items()}
    pooled = np.zeros((n_bootstrap, len(shared)), dtype=float)

    for b in range(n_bootstrap):
        tables = []
        for study in studies:
            tables.append(
                _effects_for_resample(
                    study_expression[study],
                    study_phenotypes[study],
                    per_study_resamples[study][b],
                    study=study,
                ).set_index("gene_symbol").loc[shared]
            )
        pooled[b, :] = _dl_two_study_pooled(
            tables[0]["effect"].to_numpy(dtype=float),
            tables[0]["se"].to_numpy(dtype=float),
            tables[1]["effect"].to_numpy(dtype=float),
            tables[1]["se"].to_numpy(dtype=float),
        )

    summary = pd.DataFrame(
        {
            "gene_symbol": shared,
            "n_boot": n_bootstrap,
            "effect_mean": np.nanmean(pooled, axis=0),
            "effect_sd": np.nanstd(pooled, axis=0, ddof=1),
            "effect_q025": np.nanquantile(pooled, 0.025, axis=0),
            "effect_q975": np.nanquantile(pooled, 0.975, axis=0),
        }
    ).sort_values("gene_symbol")
    summary_path = Path(output_dir) / "bootstrap_effect_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)

    index_path = Path(output_dir) / "bootstrap_resample_index.json"
    index_path.write_text(
        json.dumps(
            {
                "n_bootstrap": n_bootstrap,
                "seed": seed,
                "method": "stratified_preserving_group_sizes",
                "propagated_statistic": "hedges_g_then_dersimonian_laird_meta",
                "degenerate_resamples_sum_over_studies": degenerate_total,
                "n_shared_genes": len(shared),
                "resamples": per_study_resamples,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    return {
        "summary_path": str(summary_path),
        "index_path": str(index_path),
        "degenerate_resamples_sum_over_studies": degenerate_total,
        "n_bootstrap": n_bootstrap,
    }
