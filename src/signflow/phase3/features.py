"""Node features from Phase 2.5 artifacts (unscaled) with honest activity density QC."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from signflow.phase2 import PHASE2_META_ROOT, PHASE2_OUTPUT_ROOT, PROJECT_ROOT
from signflow.phase2.target_diagnostics import random_effects_weights

from . import FROZEN_PLAN, StopAndReportError

# Core expression/meta features (always present for meta∩network nodes)
CORE_FEATURE_COLUMNS = [
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
    "log_degree",
]

# Activity-derived (with coverage indicators — never silent zero-fill)
ACTIVITY_FEATURE_COLUMNS = [
    "collectri_tf_activity",
    "is_tf",
    "regulon_target_activity",
    "regulon_target_covered",
    "collectri_covered",
    "progeny_score",
    "progeny_covered",
]

FEATURE_COLUMNS = CORE_FEATURE_COLUMNS + ACTIVITY_FEATURE_COLUMNS

PROGENY_PATHWAYS = [
    "Androgen",
    "EGFR",
    "Estrogen",
    "Hypoxia",
    "JAK-STAT",
    "MAPK",
    "NFkB",
    "PI3K",
    "TGFb",
    "TNFa",
    "Trail",
    "VEGF",
    "WNT",
    "p53",
]


def load_min_variance() -> float:
    return float(yaml.safe_load(FROZEN_PLAN.read_text())["plan"]["feature_min_variance"])


def _gene_mean_var(expression_path: Path) -> pd.DataFrame:
    expr = pd.read_csv(expression_path, sep="\t", index_col=0)
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


def _load_pathway_activity(path: Path) -> pd.Series:
    if not path.exists():
        return pd.Series(dtype=float)
    frame = pd.read_csv(path, sep="\t", comment="#")
    num = frame.select_dtypes(include=[np.number])
    if num.empty:
        return pd.Series(dtype=float)
    return num.iloc[0].astype(float)


def _load_tf_activity(path: Path) -> pd.Series:
    if not path.exists():
        return pd.Series(dtype=float)
    frame = pd.read_csv(path, sep="\t", comment="#")
    num = frame.select_dtypes(include=[np.number])
    if num.empty:
        return pd.Series(dtype=float)
    return num.iloc[0].astype(float)


def _progeny_gene_features(
    nodes: list[str],
    pathway_activity: pd.Series,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Per-gene PROGENy: sum_p (gene_weight_p × pathway_activity_p) + coverage indicator."""

    import decoupler as dc

    net = dc.op.progeny(organism="human", top=100, license="academic", verbose=False)
    # net: source=pathway, target=gene, weight=footprint
    act = pathway_activity.reindex(PROGENY_PATHWAYS).fillna(0.0)
    net = net[net["source"].isin(PROGENY_PATHWAYS)].copy()
    net["pathway_activity"] = net["source"].map(act)
    net["contrib"] = net["weight"].astype(float) * net["pathway_activity"].astype(float)
    score = net.groupby("target")["contrib"].sum()
    covered = set(net["target"].astype(str))

    rows = []
    for gene in nodes:
        if gene in covered:
            rows.append({"gene_symbol": gene, "progeny_score": float(score.get(gene, 0.0)), "progeny_covered": 1})
        else:
            # Structural absence — score left as NaN until indicator-aware fill
            rows.append({"gene_symbol": gene, "progeny_score": float("nan"), "progeny_covered": 0})
    frame = pd.DataFrame(rows)
    meta = {
        "formulation": "sum_p(gene_weight_p * pathway_activity_p); progeny_covered indicator",
        "n_progeny_genes_in_model": int(len(covered)),
        "n_nodes_covered": int(frame["progeny_covered"].sum()),
        "pathways": PROGENY_PATHWAYS,
    }
    return frame, meta


def _collectri_gene_features(
    nodes: list[str],
    tf_activity: pd.Series,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """TF activity for TF nodes + is_tf + regulon_target_activity for targets."""

    import decoupler as dc

    net = dc.op.collectri(organism="human", remove_complexes=False, license="academic", verbose=False)
    # source=TF, target=gene, weight=sign
    tfs = set(net["source"].astype(str))
    tf_act = tf_activity.astype(float)

    # Regulon target activity: sum over TFs regulating gene of (TF_activity * weight)
    net = net.copy()
    net["tf_act"] = net["source"].map(tf_act)
    # TFs without measured activity are structural for that edge — exclude from sum, track coverage
    net_has = net.dropna(subset=["tf_act"])
    regulon = net_has.assign(contrib=net_has["tf_act"] * net_has["weight"].astype(float)).groupby("target")[
        "contrib"
    ].sum()
    targets_covered = set(net["target"].astype(str))

    rows = []
    for gene in nodes:
        is_tf = 1 if gene in tfs else 0
        if is_tf and gene in tf_act.index and np.isfinite(tf_act.loc[gene]):
            tf_val = float(tf_act.loc[gene])
            tf_structural = 0
        elif is_tf:
            tf_val = float("nan")
            tf_structural = 1
        else:
            tf_val = float("nan")
            tf_structural = 1  # not a TF → structural absence of TF-activity feature

        if gene in targets_covered and gene in regulon.index:
            reg_val = float(regulon.loc[gene])
            reg_structural = 0
        elif gene in targets_covered:
            # in regulon as target but no TF with activity — structural
            reg_val = float("nan")
            reg_structural = 1
        else:
            reg_val = float("nan")
            reg_structural = 1

        covered = 1 if (is_tf and tf_structural == 0) or reg_structural == 0 else 0
        rows.append(
            {
                "gene_symbol": gene,
                "collectri_tf_activity": tf_val,
                "is_tf": is_tf,
                "regulon_target_activity": reg_val,
                "regulon_target_covered": 0 if reg_structural else 1,
                "collectri_covered": covered,
            }
        )
    frame = pd.DataFrame(rows)
    meta = {
        "formulation": "TF activity for TF nodes; is_tf; regulon_target_activity = sum_tf(act_tf * weight)",
        "n_tfs_in_model": int(len(tfs)),
        "n_nodes_is_tf": int(frame["is_tf"].sum()),
        "n_nodes_collectri_covered": int(frame["collectri_covered"].sum()),
    }
    return frame, meta


def feature_density_qc(features: pd.DataFrame, columns: list[str] | None = None) -> pd.DataFrame:
    """Per-feature density: nonzero, structural absent, true missing, variance."""

    columns = columns or FEATURE_COLUMNS
    rows = []
    for col in columns:
        if col not in features.columns:
            continue
        s = features[col]
        # Companion indicator columns mark structural absence for activity scores
        structural = pd.Series(False, index=features.index)
        if col == "progeny_score" and "progeny_covered" in features.columns:
            structural = features["progeny_covered"] == 0
        elif col == "collectri_tf_activity" and "is_tf" in features.columns:
            structural = features["is_tf"] == 0
        elif col == "regulon_target_activity" and "regulon_target_covered" in features.columns:
            structural = features["regulon_target_covered"] == 0

        true_missing = s.isna() & ~structural
        observed = s[~structural & s.notna()]
        n_nonzero = int((observed != 0).sum()) if len(observed) else 0
        var = float(observed.var()) if len(observed) > 1 else (0.0 if len(observed) == 1 else float("nan"))
        rows.append(
            {
                "feature": col,
                "n": int(len(s)),
                "n_nonzero": n_nonzero,
                "fraction_nonzero": float(n_nonzero / len(s)) if len(s) else float("nan"),
                "n_structural_absent": int(structural.sum()),
                "n_true_missing": int(true_missing.sum()),
                "min": float(observed.min()) if len(observed) else float("nan"),
                "max": float(observed.max()) if len(observed) else float("nan"),
                "variance": var,
                "near_constant_flag": bool(np.isfinite(var) and var < load_min_variance()),
            }
        )
    return pd.DataFrame(rows)


def assert_feature_variance(features: pd.DataFrame, columns: list[str] | None = None) -> None:
    """Fail if any training feature has variance below the frozen threshold (among applicable rows)."""

    thr = load_min_variance()
    qc = feature_density_qc(features, columns)
    bad = qc.loc[qc["near_constant_flag"] & ~qc["feature"].isin(["is_tf", "progeny_covered", "collectri_covered"])]
    # Indicators may be imbalanced but binary with both classes is fine; near-constant continuous scores fail
    continuous = [
        c
        for c in (columns or FEATURE_COLUMNS)
        if c
        not in {
            "is_tf",
            "progeny_covered",
            "collectri_covered",
            "regulon_target_covered",
        }
    ]
    bad = qc.loc[qc["feature"].isin(continuous) & qc["near_constant_flag"]]
    if len(bad):
        raise StopAndReportError(
            f"Near-constant feature columns below variance threshold {thr}: {bad['feature'].tolist()}",
            reason="feature-near-constant",
            details=bad.to_dict(orient="records"),
        )


def build_feature_table(
    nodes: list[str],
    *,
    track: str,
    degrees: pd.Series,
    meta_path: Path | None = None,
    study_a: str = "GSE11783",
    study_b: str = "GSE57560",
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    """Assemble unscaled features. Returns (features, density_qc, activity_meta)."""

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
        # Prefer benchmark activity if present; else TARGET meta activity with explicit note
        bench_act = PHASE2_META_ROOT / "activity_meta_benchmark"
        activity_dir = bench_act if bench_act.exists() else PHASE2_META_ROOT / "activity_meta"
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

    tf_act = _load_tf_activity(activity_dir / "collectri_activities.tsv")
    path_act = _load_pathway_activity(activity_dir / "progeny_activities.tsv")
    collectri_df, collectri_meta = _collectri_gene_features(nodes, tf_act)
    progeny_df, progeny_meta = _progeny_gene_features(nodes, path_act)

    rows = []
    for gene in nodes:
        if gene not in meta.index:
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
                "log_degree": float(np.log1p(deg)),
                "study_a": study_a,
                "study_b": study_b,
            }
        )
    features = pd.DataFrame(rows).merge(collectri_df, on="gene_symbol", how="left").merge(
        progeny_df, on="gene_symbol", how="left"
    )

    # Fill: structural absence → 0 WITH indicator already present; true missing on core → drop row
    required = ["pooled_effect", "pooled_se", "z", "i2", "log_degree", "mean_expression", "expression_variance"]
    bad = features[required].isna().any(axis=1) | ~np.isfinite(features[required]).all(axis=1)
    if bad.any():
        features = features.loc[~bad].copy()

    # Activity scores: NaN where structural → 0 (indicator distinguishes); never claim missingness=0 without indicators
    for col, covered_col in (
        ("progeny_score", "progeny_covered"),
        ("collectri_tf_activity", "is_tf"),
        ("regulon_target_activity", "collectri_covered"),
    ):
        features[col] = features[col].fillna(0.0)

    for col in ("is_tf", "progeny_covered", "collectri_covered"):
        features[col] = features[col].fillna(0).astype(int)

    # Drop helper columns
    features = features.drop(columns=[c for c in features.columns if c.startswith("_")], errors="ignore")

    if not np.isfinite(features[FEATURE_COLUMNS]).all().all():
        raise StopAndReportError("Feature matrix contains non-finite values", reason="feature-nan-inf")

    # Density after indicator-aware fill (structural absence distinguished via companion columns)
    density = feature_density_qc(features)
    assert_feature_variance(features)
    activity_meta = {
        "collectri": collectri_meta,
        "progeny": progeny_meta,
        "activity_dir": str(activity_dir),
        "dropped_dead_columns": ["progeny_activity", "collectri_activity"],
        "note": (
            "Prior Phase 3 zero-fill without indicators made structural absence look like data. "
            "Replaced with coverage indicators + correct per-gene formulations."
        ),
    }
    return features, density, activity_meta
