"""Feature scaling policy, density QC, and split disjointness."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from signflow.phase2 import StopAndReportError
from signflow.phase3.features import FEATURE_COLUMNS, assert_feature_variance, feature_density_qc
from signflow.phase3.splits import assert_disjoint, stratified_cv_splits


def test_feature_columns_include_weight_shares_activity_indicators_and_log_degree():
    assert "weight_share_study_a" in FEATURE_COLUMNS
    assert "log_degree" in FEATURE_COLUMNS
    assert "progeny_score" in FEATURE_COLUMNS
    assert "progeny_covered" in FEATURE_COLUMNS
    assert "is_tf" in FEATURE_COLUMNS
    assert "regulon_target_activity" in FEATURE_COLUMNS
    assert "progeny_activity" not in FEATURE_COLUMNS
    assert "collectri_activity" not in FEATURE_COLUMNS


def test_density_qc_distinguishes_structural_absence():
    features = pd.DataFrame(
        {
            "progeny_score": [1.0, 0.0, 2.0],
            "progeny_covered": [1, 0, 1],
            "pooled_effect": [0.1, 0.2, 0.3],
            "is_tf": [1, 0, 0],
            "collectri_tf_activity": [0.5, 0.0, 0.0],
            "regulon_target_activity": [0.0, 1.0, 0.0],
            "regulon_target_covered": [0, 1, 0],
            "collectri_covered": [1, 1, 0],
            "log_degree": [1.0, 2.0, 1.5],
            "pooled_se": [0.1, 0.1, 0.1],
            "z": [1.0, 2.0, 3.0],
            "i2": [0.0, 10.0, 5.0],
            "effect_study_a": [0.1, 0.2, 0.3],
            "effect_study_b": [0.1, 0.2, 0.3],
            "weight_share_study_a": [0.6, 0.5, 0.7],
            "weight_share_study_b": [0.4, 0.5, 0.3],
            "mean_expression": [1.0, 2.0, 3.0],
            "expression_variance": [0.1, 0.2, 0.3],
        }
    )
    qc = feature_density_qc(features)
    prog = qc.loc[qc["feature"] == "progeny_score"].iloc[0]
    assert prog["n_structural_absent"] == 1
    assert prog["n_nonzero"] == 2


def test_assert_feature_variance_fails_on_near_constant():
    features = pd.DataFrame(
        {
            "pooled_effect": [1.0, 2.0, 3.0],
            "dead_col": [0.0, 0.0, 0.0],
            "is_tf": [0, 1, 0],
            "progeny_covered": [1, 1, 1],
            "collectri_covered": [1, 1, 1],
            "regulon_target_covered": [1, 1, 1],
        }
    )
    with pytest.raises(StopAndReportError, match="Near-constant"):
        assert_feature_variance(features, columns=["pooled_effect", "dead_col"])


def test_stratified_splits_disjoint_and_report_positives():
    y = pd.Series([1, 1, 1, 1, 0, 0, 0, 0, 1, 0])
    folds = stratified_cv_splits(y, n_folds=2, seed=0)
    for fold in folds:
        assert_disjoint(fold)
        assert len(fold["train"]) + len(fold["val"]) + len(fold["test"]) == len(y)
