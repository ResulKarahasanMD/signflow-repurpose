"""Stop-and-report guards for disease-agnostic pool floor/ceiling."""

from __future__ import annotations

import pandas as pd
import pytest

from signflow.phase5 import StopAndReportError
from signflow.phase5.acquire import assert_pool_guards
from signflow.phase5.score import build_signed_adjacency, signed_flow
from signflow.phase5.substrate import connectivity_report, signed_causal_edges


def test_gold_prevalence_ceiling_hard_fails():
    # 29 gold / 30 total = 0.967 > 0.10 — the voided failure mode
    uni = pd.DataFrame(
        {
            "chembl_id": [f"CHEMBL{i}" for i in range(30)],
            "is_gold": [True] * 29 + [False],
        }
    )
    with pytest.raises(StopAndReportError) as exc:
        assert_pool_guards(uni, n_in_substrate=30)
    assert exc.value.reason == "gold-prevalence-ceiling-exceeded"


def test_pool_floor_hard_fails():
    uni = pd.DataFrame({"chembl_id": [f"CHEMBL{i}" for i in range(100)], "is_gold": [False] * 100})
    with pytest.raises(StopAndReportError) as exc:
        assert_pool_guards(uni, n_in_substrate=10)
    assert exc.value.reason == "candidate-pool-below-floor"


def test_healthy_pool_passes():
    uni = pd.DataFrame(
        {
            "chembl_id": [f"CHEMBL{i}" for i in range(1000)],
            "is_gold": [True] * 20 + [False] * 980,
        }
    )
    report = assert_pool_guards(uni, n_in_substrate=800)
    assert report["gold_prevalence"] == pytest.approx(0.02)


def test_unsigned_adjacency_ignores_sign():
    edges = pd.DataFrame(
        {
            "source": ["A", "B"],
            "target": ["B", "C"],
            "sign": [1, -1],
            "relation_type": ["stimulation", "inhibition"],
        }
    )
    nodes = ["A", "B", "C"]
    a_signed = build_signed_adjacency(edges, nodes, ignore_sign=False)
    a_unsigned = build_signed_adjacency(edges, nodes, ignore_sign=True)
    assert float(a_signed[1, 2]) < 0
    assert float(a_unsigned[1, 2]) > 0


def test_mixed_injection_not_zeroed_by_sum():
    edges = pd.DataFrame(
        {"source": ["A"], "target": ["B"], "sign": [1], "relation_type": ["stimulation"]}
    )
    flow = signed_flow(["A", "B"], edges, {"A": 0.5, "B": -0.5}, alpha=0.5, max_iter=10)
    assert flow.any()


def test_connectivity_report_counts_components():
    edges = pd.DataFrame(
        {
            "source": ["A", "C"],
            "target": ["B", "D"],
            "sign": [1, -1],
            "relation_type": ["stimulation", "inhibition"],
        }
    )
    report = connectivity_report(edges, all_nodes={"A", "B", "C", "D", "E"}, track="test")
    assert report["n_components"] == 3  # A-B, C-D, E isolate
    assert report["n_isolated_nodes"] == 1
    assert report["largest_component_size"] == 2
