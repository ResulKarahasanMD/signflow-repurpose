"""Node features from Phase 2.5 artifacts (unscaled)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from signflow.phase2 import PHASE2_META_ROOT, PHASE2_OUTPUT_ROOT
from signflow.phase2.target_diagnostics import random_effects_weights

from . import StopAndReportError


FEATURE_COLUMNS = [
    "pooled_effect",
    "pooled_se",
    "z",
    "i2",
    "effect_study_a",
    "effect_study_b",
    "weight_share_study_a",
    "weight_share_study_b",
    "mean_expression",
    "expression_variance",
    "collectri_activity",
    "progeny_activity",
    "log_degree",
]


def _gene_mean_var(expression_path: Path) -> pd.DataFrame:
    expr = pd.read_csv(expression_path, sep="\t", index_col=0)
    # genes x samples or samples x genes — Phase 2 writes gene_symbol as first col / index
    if "gene_symbol" in expr.columns:
        expr = expr.set_index("gene_symbol")
    numeric = expr.select_dtypes(include=[np.number])
    return pd.DataFrame(
        {
            "gene_symbol": numeric.index.astype(str),
            "mean_expression": numeric.mean(axis=1).to_numpy(),
            "expression_variance": numeric.var(axis=1).to_numpy(),
        }
    )


def _activity_lookup(path: Path, value_col_candidates: tuple[str, ...]) -> pd.Series:
    """Return a gene-indexed series.

    CollecTRI/PROGENy meta tables are wide (entities as columns). Values map onto
    matching gene symbols; non-matching genes get 0 at assemble time.
    """

    if not path.exists():
        return pd.Series(dtype=float)
    frame = pd.read_csv(path, sep="\t", comment="#")
    if "gene_symbol" in frame.columns:
        for col in value_col_candidates:
            if col in frame.columns:
                return frame.set_index("gene_symbol")[col]
        num = frame.set_index("gene_symbol").select_dtypes(include=[np.number])
        return num.mean(axis=1)
    # Wide TF/pathway x sample: take first numeric row as meta scores keyed by column name
    num = frame.select_dtypes(include=[np.number])
    if num.empty:
        return pd.Series(dtype=float)
    return num.iloc[0].astype(float)


def build_feature_table(
    nodes: list[str],
    *,
    track: str,
    degrees: pd.Series,
    meta_path: Path | None = None,
    study_a: str = "GSE11783",
    study_b: str = "GSE57560",
) -> tuple[pd.DataFrame, dict[str, float]]:
    """Assemble unscaled features for the node universe."""

    if track == "target":
        meta_path = meta_path or (PHASE2_META_ROOT / "signature_meta.tsv")
        de_a = pd.read_csv(PHASE2_OUTPUT_ROOT / study_a / "de_limma.tsv", sep="\t")
        de_b = pd.read_csv(PHASE2_OUTPUT_ROOT / study_b / "de_limma.tsv", sep="\t")
        expr_paths = [
            PHASE2_OUTPUT_ROOT / study_a / "expression_gene_level.tsv",
            PHASE2_OUTPUT_ROOT / study_b / "expression_gene_level.tsv",
        ]
        activity_dir = PHASE2_META_ROOT / "activity_meta"
    elif track == "benchmark":
        meta_path = meta_path or (PHASE2_META_ROOT / "signature_meta_benchmark.tsv")
        de_a = pd.read_csv(PHASE2_OUTPUT_ROOT / "GSE75214" / "de_limma.tsv", sep="\t")
        de_b = pd.read_csv(PHASE2_OUTPUT_ROOT / "GSE87466" / "de_limma.tsv", sep="\t")
        expr_paths = [
            PHASE2_OUTPUT_ROOT / "GSE75214" / "expression_gene_level.tsv",
            PHASE2_OUTPUT_ROOT / "GSE87466" / "expression_gene_level.tsv",
        ]
        activity_dir = PHASE2_META_ROOT / "activity_meta"  # may be target-only; missingness reported
        study_a, study_b = "GSE75214", "GSE87466"
    else:
        raise ValueError(track)

    meta = pd.read_csv(meta_path, sep="\t").set_index("gene_symbol")
    a = de_a.set_index("gene_symbol")
    b = de_b.set_index("gene_symbol")

    mean_var_parts = []
    for path in expr_paths:
        if path.exists():
            mean_var_parts.append(_gene_mean_var(path).set_index("gene_symbol"))
    if mean_var_parts:
        mv = sum(mean_var_parts) / len(mean_var_parts)
    else:
        mv = pd.DataFrame(columns=["mean_expression", "expression_variance"])

    collectri = _activity_lookup(activity_dir / "collectri_activities.tsv", ("activity", "score", "estimate"))
    progeny = _activity_lookup(activity_dir / "progeny_activities.tsv", ("activity", "score", "estimate"))

    rows = []
    for gene in nodes:
        if gene not in meta.index:
            # still emit row with NaNs for missingness report — then fail if any NaN in required?
            # Plan: assert no NaN/inf after fill of structural missing with 0 for activity only? 
            # Better: require meta membership for node universe
            raise StopAndReportError(
                f"Node {gene} missing from meta signature",
                reason="feature-meta-missing",
            )
        m = meta.loc[gene]
        ea = float(a.loc[gene, "effect"]) if gene in a.index else float("nan")
        eb = float(b.loc[gene, "effect"]) if gene in b.index else float("nan")
        if gene in a.index and gene in b.index:
            wa, wb, _ = random_effects_weights(
                float(a.loc[gene, "effect"]),
                float(a.loc[gene, "se"]),
                float(b.loc[gene, "effect"]),
                float(b.loc[gene, "se"]),
            )
        else:
            wa, wb = float("nan"), float("nan")
        deg = float(degrees.get(gene, 0.0))
        rows.append(
            {
                "gene_symbol": gene,
                "pooled_effect": float(m["pooled_effect"]),
                "pooled_se": float(m["pooled_se"]),
                "z": float(m["z"]),
                "i2": float(m["i2"]),
                "effect_study_a": ea,
                "effect_study_b": eb,
                "weight_share_study_a": wa,
                "weight_share_study_b": wb,
                "mean_expression": float(mv.loc[gene, "mean_expression"]) if gene in mv.index else float("nan"),
                "expression_variance": float(mv.loc[gene, "expression_variance"]) if gene in mv.index else float("nan"),
                "collectri_activity": float(collectri.get(gene, 0.0)) if len(collectri) else 0.0,
                "progeny_activity": float(progeny.get(gene, 0.0)) if len(progeny) else 0.0,
                "log_degree": float(np.log1p(deg)),
                "study_a": study_a,
                "study_b": study_b,
            }
        )
    features = pd.DataFrame(rows)
    missingness = {}
    for col in FEATURE_COLUMNS:
        missingness[col] = float(features[col].isna().mean())
    # Assert finite for core columns after dropping rows with NaN in required
    required = ["pooled_effect", "pooled_se", "z", "i2", "log_degree"]
    bad = features[required].isna().any(axis=1) | ~np.isfinite(features[required]).all(axis=1)
    if bad.any():
        features = features.loc[~bad].copy()
    # Fill remaining optional NaNs with 0 and report
    for col in FEATURE_COLUMNS:
        if features[col].isna().any():
            features[col] = features[col].fillna(0.0)
    if not np.isfinite(features[FEATURE_COLUMNS]).all().all():
        raise StopAndReportError("Feature matrix contains non-finite values", reason="feature-nan-inf")
    return features, missingness
