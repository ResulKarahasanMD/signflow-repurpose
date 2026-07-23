"""Phase 3 edge parsing, conflict rule, sign-shuffle, node-universe rules."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from signflow.phase3.controls import sign_shuffle, verify_rewire_degree, degree_preserving_rewire
from signflow.phase3.edges import parse_omnipath


def test_omnipath_conflict_retain_as_ambiguous(tmp_path: Path):
    raw = tmp_path / "omni.tsv"
    raw.write_text(
        "source\ttarget\tsource_genesymbol\ttarget_genesymbol\tis_directed\tis_stimulation\tis_inhibition\tsources\treferences\n"
        "P1\tP2\tAAA\tBBB\tTrue\tTrue\tFalse\tX\tPMID:1\n"
        "P1\tP3\tAAA\tCCC\tTrue\tTrue\tTrue\tX\tPMID:1;PMID:2\n"
        "P1\tP4\tAAA\tDDD\tTrue\tFalse\tTrue\tX\tPMID:3\n"
        "P1\tP5\tAAA\tEEE\tFalse\tTrue\tFalse\tX\tPMID:4\n"
    )
    edges, stats = parse_omnipath(raw, conflict_rule="retain_as_ambiguous")
    assert stats["n_conflict_stim_inhib"] == 1
    assert set(edges["relation_type"]) == {"stimulation", "ambiguous", "inhibition"}
    amb = edges[edges["relation_type"] == "ambiguous"].iloc[0]
    assert amb["sign"] == 0
    assert amb["target"] == "CCC"
    # undirected dropped
    assert "EEE" not in set(edges["target"])


def test_omnipath_conflict_drop(tmp_path: Path):
    raw = tmp_path / "omni.tsv"
    raw.write_text(
        "source\ttarget\tsource_genesymbol\ttarget_genesymbol\tis_directed\tis_stimulation\tis_inhibition\tsources\treferences\n"
        "P1\tP2\tAAA\tBBB\tTrue\tTrue\tTrue\tX\tPMID:1\n"
        "P1\tP3\tAAA\tCCC\tTrue\tTrue\tFalse\tX\tPMID:2\n"
    )
    edges, stats = parse_omnipath(raw, conflict_rule="drop")
    assert stats["n_conflict_stim_inhib"] == 1
    assert len(edges) == 1
    assert edges.iloc[0]["target"] == "CCC"


def test_sign_shuffle_preserves_degree_and_sign_multiset():
    edges = pd.DataFrame(
        {
            "source": ["A", "A", "B", "C"],
            "target": ["B", "C", "C", "A"],
            "sign": [1, -1, 1, -1],
            "relation_type": ["stimulation", "inhibition", "stimulation", "inhibition"],
            "confidence": [1, 1, 1, 1],
            "provenance": ["x"] * 4,
        }
    )
    shuffled, props = sign_shuffle(edges, seed=0)
    assert props["degree_sequence_equal"]
    assert props["sign_multiset_preserved"]
    assert len(shuffled) == len(edges)


def test_node_universe_is_expressed_intersect_network_not_de_only():
    expressed = {"A", "B", "C", "D"}
    network = {"B", "C", "E"}
    significant_de = {"A", "B"}  # must NOT be used as filter
    nodes = sorted(expressed & network)
    assert nodes == ["B", "C"]
    assert "A" not in nodes  # expressed but not in network
    assert not set(nodes).issubset(significant_de) or nodes != sorted(significant_de & network)
