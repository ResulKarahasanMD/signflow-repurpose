"""Cross-study concordance summaries for Phase 2 DE results."""

from __future__ import annotations

import pandas as pd


def _corr_with_ci(x: pd.Series, y: pd.Series, *, method: str) -> tuple[float, float, float]:
    import math

    import numpy as np
    from scipy import stats

    if method == "pearson":
        r, _ = stats.pearsonr(x, y)
    else:
        r, _ = stats.spearmanr(x, y)
    n = len(x)
    if n < 4 or abs(r) >= 1:
        return float(r), float("nan"), float("nan")
    z = np.arctanh(r)
    se = 1.0 / math.sqrt(n - 3)
    lo, hi = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
    return float(r), float(lo), float(hi)


def summarize_concordance(
    left: pd.DataFrame,
    right: pd.DataFrame,
    *,
    top_ks: tuple[int, ...] = (50, 100, 250),
) -> dict[str, object]:
    """Summarize concordance between two DE result tables (pre-registered result)."""

    left_cols = ["gene_symbol", "effect"] + (["P.Value"] if "P.Value" in left.columns else [])
    right_cols = ["gene_symbol", "effect"] + (["P.Value"] if "P.Value" in right.columns else [])
    merged = left[left_cols].merge(right[right_cols], on="gene_symbol", suffixes=("_left", "_right"), how="inner")
    merged = merged.dropna(subset=["effect_left", "effect_right"])
    if merged.empty:
        return {
            "n_shared_genes": 0,
            "pearson_r": float("nan"),
            "pearson_ci": [float("nan"), float("nan")],
            "spearman_r": float("nan"),
            "spearman_ci": [float("nan"), float("nan")],
            "directional_agreement_nominal": float("nan"),
            "top_k_overlap": {},
            "plain_statement": "No shared genes; concordance undefined.",
        }

    pearson_r, pearson_lo, pearson_hi = _corr_with_ci(merged["effect_left"], merged["effect_right"], method="pearson")
    spearman_r, spearman_lo, spearman_hi = _corr_with_ci(merged["effect_left"], merged["effect_right"], method="spearman")

    if "P.Value_left" in merged.columns and "P.Value_right" in merged.columns:
        nominal = merged[(merged["P.Value_left"] < 0.05) | (merged["P.Value_right"] < 0.05)]
        directional = float((nominal["effect_left"] * nominal["effect_right"] > 0).mean()) if not nominal.empty else float("nan")
    else:
        directional = float((merged["effect_left"] * merged["effect_right"] > 0).mean())

    overlaps = {}
    for k in top_ks:
        top_left = set(left.reindex(left["effect"].abs().sort_values(ascending=False).index).head(k)["gene_symbol"])
        top_right = set(right.reindex(right["effect"].abs().sort_values(ascending=False).index).head(k)["gene_symbol"])
        overlaps[str(k)] = int(len(top_left & top_right))

    if abs(spearman_r) < 0.2:
        plain = (
            "Low concordance between GSE11783 and GSE57560 effect sizes; "
            "the cohorts do not appear to measure the same signature biology at the gene level. "
            "Phenotype difference (Hunner-selected vs capacity-stratified) is a candidate explanation, not established."
        )
    elif abs(spearman_r) < 0.5:
        plain = (
            "Moderate concordance between cohorts; shared signal is partial. "
            "Phenotype mix (Hunner-selected vs capacity-stratified) remains a candidate explanation, not established."
        )
    else:
        plain = (
            "Higher concordance between cohorts suggests overlapping biology at the gene-effect level; "
            "phenotype differences may still contribute residual discordance."
        )

    return {
        "n_shared_genes": int(len(merged)),
        "pearson_r": pearson_r,
        "pearson_ci": [pearson_lo, pearson_hi],
        "spearman_r": spearman_r,
        "spearman_ci": [spearman_lo, spearman_hi],
        "directional_agreement_nominal": directional,
        "top_k_overlap": overlaps,
        "plain_statement": plain,
    }
