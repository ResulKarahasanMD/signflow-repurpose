"""Tests for Phase 5 §A completions (three-arm, coverage, soft positives, H7 helpers)."""

from __future__ import annotations

import pandas as pd
import pytest

from signflow.phase5.evaluate import (
    acquired_pool_sensitivity,
    compare_connectivity_cost,
    differential_coverage_rates,
    report_soft_positive_ranks,
)
from signflow.phase5.h7 import pairwise_jaccard
from signflow.phase5.score import build_signed_adjacency
from signflow.phase5.substrate import signed_causal_edges, unsigned_full_causal_edges


def test_unsigned_full_keeps_unsigned_directed_excludes_ambiguous():
    edges = pd.DataFrame(
        {
            "source": ["A", "B", "C", "D"],
            "target": ["B", "C", "D", "E"],
            "relation_type": ["stimulation", "unsigned_directed", "ambiguous", "inhibition"],
            "sign": [1, 0, 0, -1],
        }
    )
    signed = signed_causal_edges(edges)
    full = unsigned_full_causal_edges(edges)
    assert set(signed["relation_type"]) == {"stimulation", "inhibition"}
    assert set(full["relation_type"]) == {"stimulation", "inhibition", "unsigned_directed"}
    assert len(full) == 3


def test_unsigned_adjacency_accepts_zero_sign_when_ignore_sign():
    edges = pd.DataFrame(
        {
            "source": ["A", "B"],
            "target": ["B", "C"],
            "sign": [1, 0],
            "relation_type": ["stimulation", "unsigned_directed"],
        }
    )
    nodes = ["A", "B", "C"]
    a = build_signed_adjacency(edges, nodes, ignore_sign=True)
    assert float(a[0, 1]) > 0
    assert float(a[1, 2]) > 0


def test_differential_coverage_and_acquired_sensitivity():
    summary = pd.DataFrame(
        {
            "chembl_id": ["g1", "g2", "n1", "n2", "n3"],
            "is_gold": [True, True, False, False, False],
            "n_targets_in_substrate": [1, 0, 1, 1, 0],
        }
    )
    # Acquired universe larger than mapped summary (A3 denominator)
    diff = differential_coverage_rates({"n_drugs_total": 5}, summary, n_acquired_universe=10)
    assert diff["gold_coverage_rate"] == 0.5
    assert diff["non_gold_coverage_rate"] == pytest.approx(2 / 3)
    assert diff["n_acquired_universe"] == 10
    assert diff["gold_prevalence_acquired"] == pytest.approx(0.2)
    # scored = {g1,n1,n2}; gold_in_scored=1 → prevalence 1/3; inflation (1/3)/0.2 - 1
    assert diff["relative_inflation"] == pytest.approx((1 / 3) / 0.2 - 1)

    ranking = pd.DataFrame(
        {
            "chembl_id": ["g1", "n1", "n2"],
            "is_gold": [True, False, False],
            "canonical_group": ["g", "n", "n"],
            "pref_name": ["G1", "N1", "N2"],
            "score": [0.9, 0.5, 0.1],
            "gold_eligible": [True, False, False],
        }
    )
    universe = pd.DataFrame(
        {
            "chembl_id": ["g1", "g2", "n1", "n2", "n3"],
            "is_gold": [True, True, False, False, False],
            "canonical_group": ["g", "g", "n", "n", "n"],
            "pref_name": ["G1", "G2", "N1", "N2", "N3"],
        }
    )
    sens = acquired_pool_sensitivity(ranking, universe, k_values=[2])
    assert sens["n_acquired_universe"] == 5
    assert sens["n_uncovered"] == 2
    assert sens["n_uncovered_gold"] == 1


def test_soft_positive_ranks_never_forced():
    ranking = pd.DataFrame(
        {
            "chembl_id": ["CHEMBL160", "CHEMBL999"],
            "pref_name": ["CYCLOSPORINE", "OTHER"],
            "score": [0.8, 0.1],
        }
    )
    soft = [
        {"name": "cyclosporine", "chembl_ids": ["CHEMBL160"], "note": "soft"},
        {"name": "pentosan", "chembl_ids": ["CHEMBL1201516"], "note": "soft"},
    ]
    rows = report_soft_positive_ranks(ranking, soft)
    assert rows[0]["rank"] == 1
    assert rows[1]["status"] == "not_in_scored_pool"


def test_pairwise_jaccard_and_connectivity_compare():
    j = pairwise_jaccard([{"a", "b"}, {"a", "c"}, {"a", "b"}])
    assert 0.0 <= j["mean"] <= 1.0
    a = pd.DataFrame({"chembl_id": ["1", "2", "3"], "is_gold": [True, False, False], "score": [0.9, 0.5, 0.1]})
    b = pd.DataFrame({"chembl_id": ["1", "2", "3"], "is_gold": [True, False, False], "score": [0.8, 0.6, 0.05]})
    out = compare_connectivity_cost(a, b, n_boot=20)
    assert "auprc_unsigned_matched" in out


def test_sensitivity_substrate_includes_string_edges():
    causal = pd.DataFrame(
        {
            "source": ["A"],
            "target": ["B"],
            "relation_type": ["stimulation"],
            "sign": [1],
        }
    )
    string = pd.DataFrame(
        {
            "source": ["B"],
            "target": ["C"],
            "sign": [0],
            "relation_type": ["string_physical"],
            "undirected": [True],
        }
    )
    from signflow.phase5.substrate import build_substrates

    built = build_substrates(causal, string)
    sens = built["sensitivity"]["signed_edges"]
    assert "string_unsigned" in set(sens["relation_type"])
    assert {"A", "B", "C"} <= set(sens["source"].astype(str)) | set(sens["target"].astype(str))
    assert built["sensitivity"]["n_sensitivity_edges"] > built["primary"]["n_signed_edges"]
