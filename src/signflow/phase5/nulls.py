"""Degree-matched permutation nulls, z-scores, BH-FDR, Spearman diagnostics."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order]
    q = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = ranked[i] * n / rank
        prev = min(prev, val)
        q[i] = prev
    out = np.empty(n, dtype=float)
    out[order] = np.clip(q, 0.0, 1.0)
    return out


def empirical_p_value(observed: float, null_draws: np.ndarray) -> float:
    null_draws = np.asarray(null_draws, dtype=float)
    if null_draws.size == 0:
        return float("nan")
    return float((1 + np.sum(null_draws >= observed)) / (1 + null_draws.size))


def z_from_null(observed: float, null_draws: np.ndarray) -> float:
    null_draws = np.asarray(null_draws, dtype=float)
    if null_draws.size < 2:
        return float("nan")
    sd = float(np.std(null_draws, ddof=1))
    if sd == 0:
        return float("nan")
    return float((observed - float(np.mean(null_draws))) / sd)


def _sample_degree_matched_genes(
    n_targets: int,
    target_mean_degree: float,
    pool_genes: np.ndarray,
    pool_degrees: np.ndarray,
    *,
    rng: np.random.Generator,
    tolerance: float,
    max_attempts: int,
) -> list[str] | None:
    if n_targets <= 0:
        return []
    if n_targets > len(pool_genes):
        return None
    lo = target_mean_degree * (1.0 - tolerance)
    hi = target_mean_degree * (1.0 + tolerance)
    for _ in range(max_attempts):
        idx = rng.choice(len(pool_genes), size=n_targets, replace=False)
        mean_d = float(pool_degrees[idx].mean())
        if target_mean_degree != target_mean_degree:  # nan mean -> any sample ok
            return pool_genes[idx].astype(str).tolist()
        if lo <= mean_d <= hi:
            return pool_genes[idx].astype(str).tolist()
    return None


def degree_matched_null(
    observed_score: float,
    target_genes: list[str],
    degrees: pd.Series,
    score_fn,
    *,
    n_perm: int,
    tolerance: float,
    seed: int,
) -> dict[str, object]:
    """Permute in-substrate targets: match count and mean degree within tolerance."""

    rng = np.random.default_rng(seed)
    genes = [g for g in target_genes if g in degrees.index]
    if not genes:
        return {
            "observed": observed_score,
            "z": float("nan"),
            "p_empirical": float("nan"),
            "null_mean": float("nan"),
            "null_sd": float("nan"),
            "n_perm_requested": n_perm,
            "n_perm_used": 0,
        }
    template_deg = degrees.loc[genes].astype(float).to_numpy()
    target_mean = float(np.mean(template_deg))
    pool_genes = degrees.index.astype(str).to_numpy()
    pool_deg = degrees.astype(float).to_numpy()
    n_targets = len(genes)

    null_scores: list[float] = []
    attempts = 0
    max_attempts = max(n_perm * 50, 5000)
    while len(null_scores) < n_perm and attempts < max_attempts:
        attempts += 1
        matched = _sample_degree_matched_genes(
            n_targets, target_mean, pool_genes, pool_deg, rng=rng, tolerance=tolerance, max_attempts=1
        )
        if matched is None:
            continue
        null_scores.append(float(score_fn(matched)))
    null_arr = np.asarray(null_scores, dtype=float)
    return {
        "observed": float(observed_score),
        "z": z_from_null(observed_score, null_arr),
        "p_empirical": empirical_p_value(observed_score, null_arr),
        "null_mean": float(np.mean(null_arr)) if null_arr.size else float("nan"),
        "null_sd": float(np.std(null_arr, ddof=1)) if null_arr.size > 1 else float("nan"),
        "n_perm_requested": int(n_perm),
        "n_perm_used": int(null_arr.size),
    }


def attach_nulls_and_fdr(
    ranking: pd.DataFrame,
    *,
    degrees: pd.Series,
    score_fn_for_genes,
    n_perm: int,
    tolerance: float,
    seed: int,
) -> pd.DataFrame:
    out_rows = []
    n_rows = len(ranking)
    for i, (_, row) in enumerate(ranking.iterrows(), start=1):
        if i == 1 or i % 100 == 0 or i == n_rows:
            print(f"[phase5] nulls {i}/{n_rows}", flush=True)
        genes = row.get("targets_in_substrate") or []
        if isinstance(genes, str):
            genes = [g for g in genes.split(",") if g]
        null = degree_matched_null(
            float(row["score"]),
            list(genes),
            degrees,
            score_fn_for_genes,
            n_perm=n_perm,
            tolerance=tolerance,
            seed=seed + int(i),
        )
        out_rows.append(null)
    null_df = pd.DataFrame(out_rows)
    merged = pd.concat([ranking.reset_index(drop=True), null_df], axis=1)
    p = merged["p_empirical"].astype(float).to_numpy()
    merged["bh_fdr"] = benjamini_hochberg(p)
    return merged


def spearman_score_degree(scores: np.ndarray, mean_degrees: np.ndarray) -> dict[str, float]:
    rho, pval = spearmanr(scores, mean_degrees, nan_policy="omit")
    return {"rho": float(rho), "p": float(pval)}
