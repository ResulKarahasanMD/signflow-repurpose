"""Negative and graph control drug sets."""

from __future__ import annotations

import yaml

import numpy as np
import pandas as pd

from . import FROZEN_PLAN


def load_plan() -> dict[str, object]:
    return yaml.safe_load(FROZEN_PLAN.read_text())["plan"]


def select_unrelated_negative_drugs(
    drug_universe: pd.DataFrame,
    *,
    plan: dict[str, object] | None = None,
) -> pd.DataFrame:
    """Sample unrelated non-gold drugs with frozen seed."""

    plan = plan or load_plan()
    seed = int(plan["phase5_negative_control_seed"])
    n = int(plan["phase5_n_unrelated_negative_drugs"])
    rng = np.random.default_rng(seed)
    pool = drug_universe.loc[~drug_universe.get("is_gold", False).astype(bool)].copy()
    pool = pool.drop_duplicates(subset=["chembl_id"])
    if len(pool) <= n:
        out = pool.copy()
    else:
        out = pool.iloc[rng.choice(len(pool), size=n, replace=False)].copy()
    out["control_kind"] = "unrelated_negative"
    return out.reset_index(drop=True)


def tag_control_graph(name: str) -> str:
    return name
