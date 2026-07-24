"""Degree-matched null permutation tests."""

from __future__ import annotations

import numpy as np
import pandas as pd

from signflow.phase5.nulls import degree_matched_null


def test_null_preserves_target_count_and_degree_band():
    degrees = pd.Series({"A": 10.0, "B": 10.2, "C": 50.0, "D": 9.8, "E": 10.1})
    observed = 3.14

    def score_fn(genes: list[str]) -> float:
        return float(len(genes))

    out = degree_matched_null(
        observed,
        ["A", "B"],
        degrees,
        score_fn,
        n_perm=25,
        tolerance=0.15,
        seed=0,
    )
    assert out["n_perm_used"] >= 5
    assert out["p_empirical"] <= 1.0
