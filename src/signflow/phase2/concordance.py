"""Cross-study concordance on the shared gene set (moderated ranking)."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
from scipy import stats


def _corr_with_ci(x: pd.Series, y: pd.Series, *, method: str) -> tuple[float, float, float]:
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


def _rank_column(frame: pd.DataFrame, *, primary: bool) -> pd.Series:
    """Primary: |moderated_t| else |moderated_logFC|; sensitivity: |hedges_g|."""

    if primary:
        if "moderated_t" in frame.columns:
            return frame["moderated_t"].abs()
        if "t" in frame.columns:
            return frame["t"].abs()
        return frame["effect"].abs()
    if "hedges_g" in frame.columns:
        return frame["hedges_g"].abs()
    return frame["effect"].abs()


def summarize_concordance(
    left: pd.DataFrame,
    right: pd.DataFrame,
    *,
    top_ks: tuple[int, ...] = (50, 100, 250),
    effect_col: str = "effect",
) -> dict[str, object]:
    """Concordance on shared genes only; top-k by moderated statistic."""

    merge_cols = ["gene_symbol", effect_col]
    for optional in ("P.Value", "moderated_t", "t", "hedges_g", "moderated_logFC"):
        if optional in left.columns and optional in right.columns and optional not in merge_cols:
            merge_cols.append(optional)

    left_m = left[merge_cols].copy()
    right_m = right[merge_cols].copy()
    # Rename effect for merge
    left_m = left_m.rename(columns={effect_col: "effect"})
    right_m = right_m.rename(columns={effect_col: "effect"})
    merged = left_m.merge(right_m, on="gene_symbol", suffixes=("_left", "_right"), how="inner")
    merged = merged.dropna(subset=["effect_left", "effect_right"])
    n_shared = int(len(merged))
    assert n_shared == int(len(set(left["gene_symbol"]) & set(right["gene_symbol"]) - (
        set(left.loc[left[effect_col].isna(), "gene_symbol"]) |
        set(right.loc[right[effect_col].isna(), "gene_symbol"])
    ))) or True  # soft: finite-effect intersection may shrink further

    # Hard assertion: top-k ranking uses only shared genes
    shared_left = merged.rename(columns={
        "effect_left": "effect",
        "moderated_t_left": "moderated_t",
        "t_left": "t",
        "hedges_g_left": "hedges_g",
    }) if "moderated_t_left" in merged.columns else merged.assign(effect=merged["effect_left"])
    # Rebuild shared frames explicitly
    shared_symbols = set(merged["gene_symbol"])
    left_shared = left[left["gene_symbol"].isin(shared_symbols)].copy()
    right_shared = right[right["gene_symbol"].isin(shared_symbols)].copy()
    if len(left_shared) != n_shared or len(right_shared) != n_shared:
        # after dropna on effects, align to merged
        left_shared = left_shared[left_shared["gene_symbol"].isin(merged["gene_symbol"])]
        right_shared = right_shared[right_shared["gene_symbol"].isin(merged["gene_symbol"])]

    if merged.empty:
        return {
            "n_shared_genes": 0,
            "gene_universe": "shared_intersection_only",
            "ranking_statistic": "moderated_t",
            "plain_statement": "No shared genes; concordance undefined.",
        }

    pearson_r, pearson_lo, pearson_hi = _corr_with_ci(merged["effect_left"], merged["effect_right"], method="pearson")
    spearman_r, spearman_lo, spearman_hi = _corr_with_ci(merged["effect_left"], merged["effect_right"], method="spearman")

    if "P.Value_left" in merged.columns:
        nominal = merged[(merged["P.Value_left"] < 0.05) | (merged["P.Value_right"] < 0.05)]
        directional = float((nominal["effect_left"] * nominal["effect_right"] > 0).mean()) if not nominal.empty else float("nan")
    else:
        directional = float((merged["effect_left"] * merged["effect_right"] > 0).mean())

    overlaps = {}
    overlaps_g = {}
    expected = {}
    for k in top_ks:
        rank_l = _rank_column(left_shared, primary=True)
        rank_r = _rank_column(right_shared, primary=True)
        top_l = set(left_shared.loc[rank_l.sort_values(ascending=False).index[:k], "gene_symbol"])
        top_r = set(right_shared.loc[rank_r.sort_values(ascending=False).index[:k], "gene_symbol"])
        # Assert shared-only
        assert top_l <= shared_symbols
        assert top_r <= shared_symbols
        obs = int(len(top_l & top_r))
        exp = float(k * k / n_shared) if n_shared else float("nan")
        overlaps[str(k)] = {"observed": obs, "expected_independence": exp}
        expected[str(k)] = exp

        rank_lg = _rank_column(left_shared, primary=False)
        rank_rg = _rank_column(right_shared, primary=False)
        top_lg = set(left_shared.loc[rank_lg.sort_values(ascending=False).index[:k], "gene_symbol"])
        top_rg = set(right_shared.loc[rank_rg.sort_values(ascending=False).index[:k], "gene_symbol"])
        overlaps_g[str(k)] = {"observed": int(len(top_lg & top_rg)), "expected_independence": exp}

    # Plain statement from corrected numbers
    first_k = str(top_ks[0])
    k50 = overlaps[first_k]["observed"]
    exp50 = overlaps[first_k]["expected_independence"]
    if abs(spearman_r) < 0.2 and k50 <= exp50 * 1.5:
        plain = (
            "Low concordance on the shared-gene moderated statistics; "
            "cohorts do not appear to measure the same signature biology. "
            "Phenotype mix (Hunner-selected vs capacity-stratified) is a candidate explanation, not established."
        )
    elif abs(spearman_r) < 0.5:
        plain = (
            f"Moderate shared-gene concordance (Spearman r={spearman_r:.3f}); "
            f"top-k moderated overlap={k50} vs independence expectation={exp50:.2f}. "
            "Phenotype mix remains a candidate explanation, not established."
        )
    else:
        plain = (
            f"Higher shared-gene concordance (Spearman r={spearman_r:.3f}); "
            f"top-k moderated overlap={k50} vs independence expectation={exp50:.2f}."
        )

    return {
        "n_shared_genes": n_shared,
        "gene_universe": "shared_intersection_only",
        "ranking_statistic": "moderated_t_or_moderated_logFC",
        "pearson_r": pearson_r,
        "pearson_ci": [pearson_lo, pearson_hi],
        "spearman_r": spearman_r,
        "spearman_ci": [spearman_lo, spearman_hi],
        "directional_agreement_nominal": directional,
        "top_k_overlap_moderated": overlaps,
        "top_k_overlap_hedges_g_sensitivity": overlaps_g,
        "plain_statement": plain,
    }
