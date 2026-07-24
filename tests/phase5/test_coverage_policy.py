"""Coverage and rediscovery denominator policy."""

from __future__ import annotations

import pandas as pd

from signflow.phase5.coverage import build_coverage_tables


def test_zero_in_substrate_gold_is_coverage_miss():
    drug_summary = pd.DataFrame(
        [
            {"chembl_id": "G1", "n_targets_in_substrate": 0, "is_gold": True, "canonical_group": "adalimumab"},
            {"chembl_id": "N1", "n_targets_in_substrate": 2, "is_gold": False, "canonical_group": "other"},
        ]
    )
    targets = pd.DataFrame({"chembl_id": ["N1", "N1"], "gene_symbol": ["TP53", "BRAF"]})
    built = build_coverage_tables(
        drug_summary,
        targets,
        primary_nodes={"TP53", "BRAF"},
        sensitivity_nodes={"TP53", "BRAF", "MYC"},
    )
    rep = built["report"]
    assert rep["n_gold_coverage_misses"] == 1
    assert rep["n_gold_eligible_for_rediscovery"] == 0
    assert "G1" in rep["gold_coverage_miss_chembl_ids"]
