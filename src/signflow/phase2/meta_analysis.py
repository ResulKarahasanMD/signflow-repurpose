"""Random-effects meta-analysis on moderated study effects."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def random_effects_meta(effects: np.ndarray, ses: np.ndarray) -> dict[str, float]:
    """DerSimonian-Laird random-effects meta-analysis (scalar)."""

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


def _dl_two_study_vectorized(e1: np.ndarray, s1: np.ndarray, e2: np.ndarray, s2: np.ndarray) -> pd.DataFrame:
    v1 = s1**2
    v2 = s2**2
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
    rwsum = rw1 + rw2
    pooled = (rw1 * e1 + rw2 * e2) / rwsum
    pooled_se = np.sqrt(1.0 / rwsum)
    z = pooled / pooled_se
    p = 2.0 * stats.norm.sf(np.abs(z))
    q_p = stats.chi2.sf(q, 1)
    i2 = np.maximum(((q - 1.0) / np.where(q > 0, q, np.nan)) * 100.0, 0.0)
    i2 = np.nan_to_num(i2, nan=0.0)
    return pd.DataFrame(
        {
            "pooled_effect": pooled,
            "pooled_se": pooled_se,
            "z": z,
            "p": p,
            "q": q,
            "q_p": q_p,
            "i2": i2,
            "tau2": tau2,
            "n_studies": 2,
        }
    )


def meta_analyze_signatures(
    study_results: list[pd.DataFrame],
    *,
    effect_col: str = "effect",
    se_col: str = "se",
) -> pd.DataFrame:
    """Meta-analyze the intersection of genes across studies (vectorized for two studies)."""

    prepared = []
    for result in study_results:
        subset = result[["gene_symbol", effect_col, se_col]].copy()
        subset = subset.rename(columns={effect_col: "effect", se_col: "se"})
        subset = subset.replace([np.inf, -np.inf], np.nan).dropna(subset=["gene_symbol", "effect", "se"])
        subset = subset[subset["se"] > 0]
        prepared.append(subset)

    shared = set(prepared[0]["gene_symbol"])
    for result in prepared[1:]:
        shared &= set(result["gene_symbol"])
    shared_genes = sorted(shared)
    if not shared_genes:
        return pd.DataFrame(
            columns=["gene_symbol", "pooled_effect", "pooled_se", "z", "p", "q", "q_p", "i2", "tau2", "n_studies", "bh_fdr"]
        )

    if len(prepared) == 2:
        left = prepared[0].set_index("gene_symbol").loc[shared_genes]
        right = prepared[1].set_index("gene_symbol").loc[shared_genes]
        out = _dl_two_study_vectorized(
            left["effect"].to_numpy(dtype=float),
            left["se"].to_numpy(dtype=float),
            right["effect"].to_numpy(dtype=float),
            right["se"].to_numpy(dtype=float),
        )
        out.insert(0, "gene_symbol", shared_genes)
    else:
        rows = []
        for gene in shared_genes:
            effects = []
            ses = []
            for table in prepared:
                row = table.loc[table["gene_symbol"] == gene].iloc[0]
                effects.append(float(row["effect"]))
                ses.append(float(row["se"]))
            rows.append({"gene_symbol": gene, **random_effects_meta(np.asarray(effects), np.asarray(ses))})
        out = pd.DataFrame(rows)

    out = out.sort_values("p").reset_index(drop=True)
    out["bh_fdr"] = multipletests(out["p"], method="fdr_bh")[1]
    return out
