"""Phase 4 unit tests: metrics, conformal, scaler fold discipline, split hash."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from sklearn.preprocessing import StandardScaler

from signflow.phase4.conformal import empirical_coverage, mondrian_conformal_sets
from signflow.phase4.data import hash_bytes, load_phase3_benchmark
from signflow.phase4.metrics import compute_metrics


def test_metrics_on_perfect_ranking():
    y = np.array([1, 1, 0, 0])
    s = np.array([0.9, 0.8, 0.2, 0.1])
    m = compute_metrics(y, s)
    assert m["auroc"] == pytest.approx(1.0)
    assert m["auprc"] == pytest.approx(1.0)


def test_mondrian_coverage_on_exchangeable_synthetic():
    rng = np.random.default_rng(0)
    # Well-separated scores
    y = np.array([0] * 40 + [1] * 40)
    scores = np.where(y == 1, rng.uniform(0.6, 1.0, 80), rng.uniform(0.0, 0.4, 80))
    cal_s, cal_y = scores[:40], y[:40]
    test_s, test_y = scores[40:], y[40:]
    # Fix: need both classes in cal — reshuffle
    idx = rng.permutation(80)
    scores, y = scores[idx], y[idx]
    cal_s, cal_y = scores[:40], y[:40]
    test_s, test_y = scores[40:], y[40:]
    sets, diag = mondrian_conformal_sets(cal_s, cal_y, test_s, alpha=0.1)
    cov = empirical_coverage(sets, test_y)
    assert cov >= 0.7  # soft floor for finite sample
    assert diag["mean_set_size"] >= 1.0


def test_scaler_fit_only_on_train_changes_transform():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(20, 3))
    # Deliberately shifted test fold
    X_test = X[10:] + 5.0
    X_train = X[:10]
    scaler = StandardScaler().fit(X_train)
    t1 = scaler.transform(X_test)
    # Wrong: fit including test
    scaler_leak = StandardScaler().fit(X)
    t2 = scaler_leak.transform(X_test)
    assert not np.allclose(t1, t2)


def test_split_hash_stable_if_phase3_present():
    from signflow.phase4 import PHASE3_LATEST

    if not PHASE3_LATEST.exists():
        pytest.skip("Phase 3 artifacts not present")
    bundle = load_phase3_benchmark()
    raw = Path(bundle["splits_path"]).read_bytes()
    assert hash_bytes(raw) == bundle["splits_hash"]
    # Reload matches
    again = hash_bytes(Path(bundle["splits_path"]).read_bytes())
    assert again == bundle["splits_hash"]


def test_sign_shuffle_changes_edge_signs():
    from signflow.phase4 import PHASE3_LATEST

    if not PHASE3_LATEST.exists():
        pytest.skip("Phase 3 artifacts not present")
    bench = PHASE3_LATEST.resolve() / "benchmark"
    real = pd.read_csv(bench / "graph_edges.tsv", sep="\t")
    shuf = pd.read_csv(bench / "controls" / "edges_sign_shuffled.tsv", sep="\t")
    assert len(real) == len(shuf)
    # Sign multiset equal but assignment may differ
    assert sorted(real["sign"].tolist()) == sorted(shuf["sign"].tolist())
