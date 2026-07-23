"""TARGET seed selection with mandatory replication annotations (Phase 2.6)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml

from . import PROJECT_ROOT
from .target_diagnostics import random_effects_weights

FROZEN_PLAN = PROJECT_ROOT / "docs" / "phase0" / "FROZEN_PLAN.yaml"


def load_seed_constants() -> dict:
    return yaml.safe_load(FROZEN_PLAN.read_text())["plan"]


def _annotate(
    meta: pd.DataFrame,
    de_a: pd.DataFrame,
    de_b: pd.DataFrame,
    genes: list[str],
    *,
    study_a: str = "GSE11783",
    study_b: str = "GSE57560",
) -> pd.DataFrame:
    a = de_a.set_index("gene_symbol")
    b = de_b.set_index("gene_symbol")
    m = meta.set_index("gene_symbol")
    rows = []
    for g in genes:
        wa, wb, _ = random_effects_weights(
            float(a.loc[g, "effect"]),
            float(a.loc[g, "se"]),
            float(b.loc[g, "effect"]),
            float(b.loc[g, "se"]),
        )
        ea, eb = float(a.loc[g, "effect"]), float(b.loc[g, "effect"])
        rows.append(
            {
                "gene_symbol": g,
                "pooled_effect": float(m.loc[g, "pooled_effect"]),
                "pooled_se": float(m.loc[g, "pooled_se"]),
                "i2": float(m.loc[g, "i2"]),
                "bh_fdr_meta": float(m.loc[g, "bh_fdr"]),
                f"effect_{study_a}": ea,
                f"se_{study_a}": float(a.loc[g, "se"]),
                f"p_{study_a}": float(a.loc[g, "P.Value"]),
                f"bh_fdr_{study_a}": float(a.loc[g, "adj.P.Val"]),
                f"effect_{study_b}": eb,
                f"se_{study_b}": float(b.loc[g, "se"]),
                f"p_{study_b}": float(b.loc[g, "P.Value"]),
                f"bh_fdr_{study_b}": float(b.loc[g, "adj.P.Val"]),
                "direction_agree": bool(ea * eb > 0),
                f"weight_share_{study_a}": wa,
                f"weight_share_{study_b}": wb,
            }
        )
    return pd.DataFrame(rows)


def select_seeds(
    meta: pd.DataFrame,
    de_a: pd.DataFrame,
    de_b: pd.DataFrame,
    *,
    require_concordant_direction: bool,
    n: int,
    i2_max: float,
) -> pd.DataFrame:
    """Primary or sensitivity seed selection on the frozen ranking."""

    a = de_a.set_index("gene_symbol")
    b = de_b.set_index("gene_symbol")
    pool = meta.copy()
    pool = pool[pool["i2"] <= i2_max]
    pool = pool[pool["gene_symbol"].isin(a.index) & pool["gene_symbol"].isin(b.index)]
    if require_concordant_direction:
        dirs = []
        for g in pool["gene_symbol"]:
            dirs.append(float(a.loc[g, "effect"]) * float(b.loc[g, "effect"]) > 0)
        pool = pool.loc[dirs]
    pool = pool.assign(_abs_eff=pool["pooled_effect"].abs()).sort_values("_abs_eff", ascending=False)
    genes = pool["gene_symbol"].astype(str).head(n).tolist()
    return _annotate(meta, de_a, de_b, genes)


def emit_target_seed_artifacts(
    meta_path: str | Path,
    de_a_path: str | Path,
    de_b_path: str | Path,
    output_dir: str | Path,
) -> dict[str, object]:
    from . import ensure_directory

    cfg = load_seed_constants()
    output_dir = ensure_directory(output_dir)
    meta = pd.read_csv(meta_path, sep="\t")
    de_a = pd.read_csv(de_a_path, sep="\t")
    de_b = pd.read_csv(de_b_path, sep="\t")

    primary = select_seeds(
        meta,
        de_a,
        de_b,
        require_concordant_direction=False,
        n=int(cfg["target_seed_n"]),
        i2_max=float(cfg["target_seed_i2_max"]),
    )
    primary["seed_set"] = "primary"
    primary["rule"] = cfg["target_seed_rule"]

    sens = select_seeds(
        meta,
        de_a,
        de_b,
        require_concordant_direction=True,
        n=int(cfg["target_seed_sensitivity_n"]),
        i2_max=float(cfg["target_seed_sensitivity_i2_max"]),
    )
    sens["seed_set"] = "sensitivity_concordant_direction"
    sens["rule"] = cfg["target_seed_sensitivity_rule"]

    primary_path = Path(output_dir) / "target_seeds_primary.tsv"
    sens_path = Path(output_dir) / "target_seeds_sensitivity_concordant.tsv"
    primary.to_csv(primary_path, sep="\t", index=False)
    sens.to_csv(sens_path, sep="\t", index=False)

    overlap = set(primary["gene_symbol"]) & set(sens["gene_symbol"])
    return {
        "primary_path": str(primary_path),
        "sensitivity_path": str(sens_path),
        "n_primary": int(len(primary)),
        "n_sensitivity": int(len(sens)),
        "n_overlap": int(len(overlap)),
        "overlap_genes": sorted(overlap),
    }
