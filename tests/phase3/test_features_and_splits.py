"""Feature scaling policy and split disjointness."""

from __future__ import annotations

import numpy as np
import pandas as pd

from signflow.phase3.features import FEATURE_COLUMNS
from signflow.phase3.splits import assert_disjoint, stratified_cv_splits


def test_feature_columns_include_weight_shares_and_log_degree():
    assert "weight_share_study_a" in FEATURE_COLUMNS
    assert "weight_share_study_b" in FEATURE_COLUMNS
    assert "log_degree" in FEATURE_COLUMNS
    assert "pooled_effect" in FEATURE_COLUMNS


def test_stratified_splits_disjoint_and_report_positives():
    y = pd.Series([1, 1, 1, 1, 0, 0, 0, 0, 1, 0])
    folds = stratified_cv_splits(y, n_folds=2, seed=0)
    for fold in folds:
        assert_disjoint(fold)
        assert len(fold["train"]) + len(fold["val"]) + len(fold["test"]) == len(y)
