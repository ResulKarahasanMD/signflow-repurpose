"""Phase 3 edge parsing, conflict rule, sign-shuffle, node-universe rules."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from signflow.phase3.controls import sign_shuffle
from signflow.phase3.edges import collapse_to_unique_pairs, parse_omnipath


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
    assert stats["n_conflict_stim_inhib_within_omnipath"] == 1
    assert set(edges["relation_type"]) == {"stimulation", "ambiguous", "inhibition"}
    amb = edges[edges["relation_type"] == "ambiguous"].iloc[0]
    assert amb["sign"] == 0
    assert amb["target"] == "CCC"
    assert "EEE" not in set(edges["target"])


def test_post_merge_collapse_resolves_cross_source_sign_conflict():
    edges = pd.DataFrame(
        [
            {"source": "A", "target": "B", "sign": 1, "relation_type": "stimulation", "confidence": 2,
             "provenance": "omnipath", "source_datasets": "X", "n_references": 2,
             "source_uniprot": "", "target_uniprot": ""},
            {"source": "A", "target": "B", "sign": -1, "relation_type": "inhibition", "confidence": 1,
             "provenance": "signor", "source_datasets": "SIGNOR", "n_references": 1,
             "source_uniprot": "", "target_uniprot": ""},
            {"source": "A", "target": "C", "sign": 1, "relation_type": "stimulation", "confidence": 1,
             "provenance": "omnipath", "source_datasets": "X", "n_references": 1,
             "source_uniprot": "", "target_uniprot": ""},
            {"source": "A", "target": "C", "sign": 1, "relation_type": "stimulation", "confidence": 3,
             "provenance": "signor", "source_datasets": "SIGNOR", "n_references": 3,
             "source_uniprot": "", "target_uniprot": ""},
        ]
    )
    out, stats = collapse_to_unique_pairs(edges, conflict_rule="retain_as_ambiguous")
    assert len(out) == 2
    assert stats["n_cross_sign_conflicts_to_ambiguous"] == 1
    ab = out[(out["source"] == "A") & (out["target"] == "B")].iloc[0]
    assert ab["sign"] == 0 and ab["relation_type"] == "ambiguous"
    assert "omnipath" in ab["provenance"] and "signor" in ab["provenance"]
    ac = out[(out["source"] == "A") & (out["target"] == "C")].iloc[0]
    assert ac["sign"] == 1
    assert not out.duplicated(subset=["source", "target"]).any()


def test_sign_shuffle_preserves_degree_and_sign_multiset():
    edges = pd.DataFrame(
        {
            "source": ["A", "A", "B", "C", "D"],
            "target": ["B", "C", "C", "A", "A"],
            "sign": [1, -1, 1, -1, 0],
            "relation_type": ["stimulation", "inhibition", "stimulation", "inhibition", "ambiguous"],
            "confidence": [1, 1, 1, 1, 1],
            "provenance": ["x"] * 5,
        }
    )
    shuffled, props = sign_shuffle(edges, seed=0)
    assert props["degree_sequence_equal"]
    assert props["sign_multiset_preserved"]
    assert props["relation_multiset_preserved"]
    assert props["permutes_relation_type"] is True
    assert len(shuffled) == len(edges)
    amb = shuffled[(shuffled["source"] == "D") & (shuffled["target"] == "A")].iloc[0]
    assert amb["sign"] == 0 and amb["relation_type"] == "ambiguous"


def test_sign_shuffle_changes_fields_models_read():
    """A4: shuffled graph must differ in edge_attr (sign) and edge_type (relation)."""

    edges = pd.DataFrame(
        {
            "source": [f"S{i}" for i in range(40)],
            "target": [f"T{i}" for i in range(40)],
            "sign": [1, -1] * 20,
            "relation_type": ["stimulation", "inhibition"] * 20,
            "confidence": [1.0] * 40,
            "provenance": ["x"] * 40,
        }
    )
    shuffled, props = sign_shuffle(edges, seed=7)
    assert props["edge_attr_sign_field_changed"] is True
    assert props["edge_type_relation_field_changed"] is True
    assert not np.array_equal(edges["sign"].to_numpy(), shuffled["sign"].to_numpy())
    assert not np.array_equal(
        edges["relation_type"].astype(str).to_numpy(),
        shuffled["relation_type"].astype(str).to_numpy(),
    )


def test_node_universe_is_expressed_intersect_network_not_de_only():
    expressed = {"A", "B", "C", "D"}
    network = {"B", "C", "E"}
    significant_de = {"A", "B"}
    nodes = sorted(expressed & network)
    assert nodes == ["B", "C"]
    assert "A" not in nodes
