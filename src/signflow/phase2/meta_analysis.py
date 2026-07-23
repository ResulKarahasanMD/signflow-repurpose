"""Random-effects meta-analysis on standardized study effects."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def random_effects_meta(effects: np.ndarray, ses: np.ndarray) -> dict[str, float]:
    """DerSimonian-Laird random-effects meta-analysis."""

    variances = ses**2
    fixed_weights = 1.0 / variances
    fixed_mean = float(np.sum(fixed_weights * effects) / np.sum(fixed_weights))
    q = float(np.sum(fixed_weights * (effects - fixed_mean) ** 2))
    df = max(len(effects) - 1, 1)
    c = float(np.sum(fixed_weights) - (np.sum(fixed_weights**2) / np.sum(fixed_weights)))
    tau2 = max((q - df) / c, 0.0) if c > 0 else 0.0
    random_weights = 1.0 / (variances + tau2)
    pooled_effect = float(np.sum(random_weights * effects) / np.sum(random_weights))
    pooled_se = float(math.sqrt(1.0 / np.sum(random_weights)))
    z_value = pooled_effect / pooled_se if pooled_se else float("nan")
    p_value = float(2.0 * stats.norm.sf(abs(z_value))) if np.isfinite(z_value) else float("nan")
    q_p = float(stats.chi2.sf(q, df))
    i2 = max(((q - df) / q) * 100.0, 0.0) if q > 0 else 0.0
    return {
        "pooled_effect": pooled_effect,
        "pooled_se": pooled_se,
        "z": z_value,
        "p": p_value,
        "q": q,
        "q_p": q_p,
        "i2": i2,
        "tau2": tau2,
        "n_studies": int(len(effects)),
    }


def meta_analyze_signatures(study_results: list[pd.DataFrame]) -> pd.DataFrame:
    """Meta-analyze the intersection of genes across studies."""

    prepared = []
    for result in study_results:
        subset = result[["gene_symbol", "effect", "se", "study"]].copy()
        prepared.append(subset.dropna(subset=["gene_symbol", "effect", "se"]))

    shared_genes = set(prepared[0]["gene_symbol"])
    for result in prepared[1:]:
        shared_genes &= set(result["gene_symbol"])
    shared_genes = sorted(shared_genes)

    rows = []
    for gene_symbol in shared_genes:
        gene_rows = pd.concat([result[result["gene_symbol"] == gene_symbol] for result in prepared], ignore_index=True)
        meta = random_effects_meta(
            gene_rows["effect"].to_numpy(dtype=float),
            gene_rows["se"].to_numpy(dtype=float),
        )
        rows.append({"gene_symbol": gene_symbol, **meta})

    output = pd.DataFrame(rows).sort_values("p").reset_index(drop=True)
    output["bh_fdr"] = multipletests(output["p"], method="fdr_bh")[1]
    return output
