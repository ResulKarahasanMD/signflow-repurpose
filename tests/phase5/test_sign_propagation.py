"""Signed substrate edge policy tests."""

from __future__ import annotations

import pandas as pd

from signflow.phase5.substrate import signed_causal_edges


def test_ambiguous_and_unsigned_excluded_from_signed_edges():
    edges = pd.DataFrame(
        {
            "source": ["A", "B", "C"],
            "target": ["B", "C", "D"],
            "relation_type": ["stimulation", "ambiguous", "inhibition"],
            "sign": [1, 0, -1],
        }
    )
    signed = signed_causal_edges(edges)
    assert len(signed) == 2
    assert set(signed["relation_type"]) == {"stimulation", "inhibition"}
