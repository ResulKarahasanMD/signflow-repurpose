"""Phase 2.6 TARGET meta replication-honesty diagnostics."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from .meta_analysis import meta_analyze_signatures


def _study_tables(de_a: pd.DataFrame, de_b: pd.DataFrame, study_a: str, study_b: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    a = de_a.copy()
    b = de_b.copy()
    a["study"] = study_a
    b["study"] = study_b
    return a, b


def random_effects_weights(effect_a: float, se_a: float, effect_b: float, se_b: float) -> tuple[float, float, float]:
    """Return (w_a_share, w_b_share, tau2) for a two-study DL meta."""

    v1, v2 = se_a**2, se_b**2
    w1, w2 = 1.0 / v1, 1.0 / v2
    wsum = w1 + w2
    fixed = (w1 * effect_a + w2 * effect_b) / wsum
    q = w1 * (effect_a - fixed) ** 2 + w2 * (effect_b - fixed) ** 2
    c = wsum - (w1**2 + w2**2) / wsum
    tau2 = max((q - 1.0) / c, 0.0) if c > 0 else 0.0
    rw1, rw2 = 1.0 / (v1 + tau2), 1.0 / (v2 + tau2)
    rsum = rw1 + rw2
    return float(rw1 / rsum), float(rw2 / rsum), float(tau2)


def characterize_meta_significant(
    meta: pd.DataFrame,
    de_a: pd.DataFrame,
    de_b: pd.DataFrame,
    *,
    study_a: str = "GSE11783",
    study_b: str = "GSE57560",
    fdr_col: str = "bh_fdr",
) -> dict[str, object]:
    """Of meta FDR<0.05 genes: counts vs study A FDR, B nominal, directional share."""

    sig = meta.loc[meta[fdr_col] < 0.05, "gene_symbol"].astype(str)
    a = de_a.set_index("gene_symbol")
    b = de_b.set_index("gene_symbol")
    shared = [g for g in sig if g in a.index and g in b.index]
    n = len(shared)
    a_fdr05 = sum(float(a.loc[g, "adj.P.Val"]) < 0.05 for g in shared)
    b_nom05 = sum(float(b.loc[g, "P.Value"]) < 0.05 for g in shared)
    same_dir = sum(float(a.loc[g, "effect"]) * float(b.loc[g, "effect"]) > 0 for g in shared)
    return {
        "n_meta_fdr_lt_0.05": int((meta[fdr_col] < 0.05).sum()),
        "n_annotated_in_both": n,
        "n_also_fdr_lt_0.05_in_study_a": int(a_fdr05),
        "n_nominal_p_lt_0.05_in_study_b": int(b_nom05),
        "n_same_direction_in_study_b": int(same_dir),
        "frac_also_fdr_study_a": float(a_fdr05 / n) if n else float("nan"),
        "frac_nominal_study_b": float(b_nom05 / n) if n else float("nan"),
        "frac_same_direction_study_b": float(same_dir / n) if n else float("nan"),
        "study_a": study_a,
        "study_b": study_b,
    }


def leave_one_study_out_counts(de_a: pd.DataFrame, de_b: pd.DataFrame) -> dict[str, object]:
    """DL 'meta' on each cohort alone = just that study's BH-FDR counts (identity)."""

    def _sig(de: pd.DataFrame) -> dict[str, int]:
        return {
            "n_genes": int(len(de)),
            "bh_fdr_lt_0.05": int((de["adj.P.Val"] < 0.05).sum()),
            "bh_fdr_lt_0.10": int((de["adj.P.Val"] < 0.10).sum()),
        }

    return {
        "GSE11783_alone": _sig(de_a),
        "GSE57560_alone": _sig(de_b),
        "note": (
            "Single-study 'meta' is the study's own limma BH-FDR table. "
            "If two-study meta FDR count < GSE11783-alone, meta reduced rather than added evidence."
        ),
    }


def weight_share_summary(
    meta: pd.DataFrame,
    de_a: pd.DataFrame,
    de_b: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, float]]:
    """Per-gene RE weight shares; mean/median genome-wide and among meta FDR<0.05."""

    a = de_a.set_index("gene_symbol")
    b = de_b.set_index("gene_symbol")
    rows = []
    for _, row in meta.iterrows():
        g = str(row["gene_symbol"])
        if g not in a.index or g not in b.index:
            continue
        wa, wb, tau2 = random_effects_weights(
            float(a.loc[g, "effect"]),
            float(a.loc[g, "se"]),
            float(b.loc[g, "effect"]),
            float(b.loc[g, "se"]),
        )
        rows.append(
            {
                "gene_symbol": g,
                "weight_share_GSE11783": wa,
                "weight_share_GSE57560": wb,
                "tau2": tau2,
                "i2": float(row.get("i2", np.nan)),
                "bh_fdr": float(row.get("bh_fdr", np.nan)),
            }
        )
    frame = pd.DataFrame(rows)
    sig = frame.loc[frame["bh_fdr"] < 0.05]
    summary = {
        "mean_weight_share_GSE11783": float(frame["weight_share_GSE11783"].mean()),
        "median_weight_share_GSE11783": float(frame["weight_share_GSE11783"].median()),
        "mean_weight_share_GSE57560": float(frame["weight_share_GSE57560"].mean()),
        "median_weight_share_GSE57560": float(frame["weight_share_GSE57560"].median()),
        "n_genes": int(len(frame)),
        "meta_sig_mean_weight_share_GSE11783": float(sig["weight_share_GSE11783"].mean()) if len(sig) else float("nan"),
        "meta_sig_median_weight_share_GSE11783": float(sig["weight_share_GSE11783"].median()) if len(sig) else float("nan"),
        "meta_sig_mean_weight_share_GSE57560": float(sig["weight_share_GSE57560"].mean()) if len(sig) else float("nan"),
        "meta_sig_median_weight_share_GSE57560": float(sig["weight_share_GSE57560"].median()) if len(sig) else float("nan"),
        "meta_sig_n_genes": int(len(sig)),
    }
    return frame, summary


def abundance_residual_correlation(
    de_a: pd.DataFrame,
    de_b: pd.DataFrame,
) -> dict[str, float]:
    """Spearman of moderated effects before/after regressing out AveExpr."""

    merged = de_a[["gene_symbol", "effect", "AveExpr"]].merge(
        de_b[["gene_symbol", "effect", "AveExpr"]],
        on="gene_symbol",
        suffixes=("_a", "_b"),
    ).dropna()
    if merged.empty:
        return {"n_shared": 0, "spearman_before": float("nan"), "spearman_after": float("nan")}

    r_before, _ = stats.spearmanr(merged["effect_a"], merged["effect_b"])

    def _resid(y: np.ndarray, x: np.ndarray) -> np.ndarray:
        x1 = np.column_stack([np.ones(len(x)), x])
        beta, _, _, _ = np.linalg.lstsq(x1, y, rcond=None)
        return y - x1 @ beta

    ra = _resid(merged["effect_a"].to_numpy(dtype=float), merged["AveExpr_a"].to_numpy(dtype=float))
    rb = _resid(merged["effect_b"].to_numpy(dtype=float), merged["AveExpr_b"].to_numpy(dtype=float))
    r_after, _ = stats.spearmanr(ra, rb)
    return {
        "n_shared": int(len(merged)),
        "spearman_before": float(r_before),
        "spearman_after_aveexpr_residual": float(r_after),
        "delta": float(r_before - r_after),
    }


def run_target_diagnostics(
    meta_path: str | Path,
    de_a_path: str | Path,
    de_b_path: str | Path,
    output_dir: str | Path,
) -> dict[str, object]:
    """Emit TARGET Phase 2.6 diagnostic tables and summary JSON fields."""

    from . import ensure_directory
    import json

    output_dir = ensure_directory(output_dir)
    meta = pd.read_csv(meta_path, sep="\t")
    de_a = pd.read_csv(de_a_path, sep="\t")
    de_b = pd.read_csv(de_b_path, sep="\t")

    char = characterize_meta_significant(meta, de_a, de_b)
    loso = leave_one_study_out_counts(de_a, de_b)
    weights, weight_summary = weight_share_summary(meta, de_a, de_b)
    abund = abundance_residual_correlation(de_a, de_b)

    two_study_fdr = int((meta["bh_fdr"] < 0.05).sum())
    gse11783_alone = loso["GSE11783_alone"]["bh_fdr_lt_0.05"]
    plain = (
        f"TARGET meta-analytic signature is dominated by a single cohort. "
        f"GSE11783-alone BH-FDR<0.05={gse11783_alone}; two-study meta={two_study_fdr}. "
        f"Genome-wide mean RE weight share GSE11783={weight_summary['mean_weight_share_GSE11783']:.3f}, "
        f"GSE57560={weight_summary['mean_weight_share_GSE57560']:.3f}; "
        f"among meta FDR<0.05 genes mean share GSE11783="
        f"{weight_summary['meta_sig_mean_weight_share_GSE11783']:.3f}, "
        f"GSE57560={weight_summary['meta_sig_mean_weight_share_GSE57560']:.3f}. "
        f"Of {char['n_meta_fdr_lt_0.05']} meta-significant genes: "
        f"{char['n_also_fdr_lt_0.05_in_study_a']} also FDR<0.05 in GSE11783; "
        f"{char['n_nominal_p_lt_0.05_in_study_b']} nominal p<0.05 in GSE57560; "
        f"{char['n_same_direction_in_study_b']} same direction in GSE57560. "
        f"The second cohort is too imprecise to confirm or refute the first."
    )
    if two_study_fdr < gse11783_alone:
        plain += " Meta-analysis reduced rather than added evidence relative to GSE11783 alone."

    abund_plain = (
        f"Cross-cohort Spearman before AveExpr residualization={abund['spearman_before']:.4f}; "
        f"after={abund['spearman_after_aveexpr_residual']:.4f} (delta={abund['delta']:.4f})."
    )
    if abund["spearman_after_aveexpr_residual"] < 0.5 * abs(abund["spearman_before"]):
        abund_plain += (
            " Shared correlation largely vanishes after removing mean-expression association — "
            "consistent with a shared abundance-related technical gradient rather than a disease signature."
        )
    else:
        abund_plain += (
            " Residual correlation remains material after AveExpr adjustment; "
            "a shared low-amplitude component is not fully explained by abundance alone."
        )

    weights.to_csv(Path(output_dir) / "target_re_weight_shares.tsv", sep="\t", index=False)
    summary = {
        "meta_significance_breakdown": char,
        "leave_one_study_out": loso,
        "weight_share_summary": weight_summary,
        "abundance_residual_correlation": abund,
        "plain_statement": plain,
        "abundance_plain_statement": abund_plain,
    }
    Path(output_dir) / "target_meta_diagnostics.json"
    (Path(output_dir) / "target_meta_diagnostics.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary
