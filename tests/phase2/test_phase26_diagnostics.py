"""Unit tests for Phase 2.6 TARGET diagnostics and seed annotations."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from signflow.phase2.seeds import select_seeds
from signflow.phase2.target_diagnostics import (
    abundance_residual_correlation,
    characterize_meta_significant,
    leave_one_study_out_counts,
    random_effects_weights,
    run_target_diagnostics,
    weight_share_summary,
)


def _de(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(rows)


def test_random_effects_weights_sum_to_one():
    wa, wb, _ = random_effects_weights(1.0, 0.1, 1.0, 1.0)
    assert wa + wb == pytest.approx(1.0)
    assert wa > wb


def test_characterize_meta_significant_counts():
    meta = pd.DataFrame(
        {
            "gene_symbol": ["A", "B", "C"],
            "bh_fdr": [0.01, 0.01, 0.2],
            "pooled_effect": [1.0, -1.0, 0.1],
            "i2": [0.0, 10.0, 0.0],
        }
    )
    de_a = _de(
        [
            {"gene_symbol": "A", "effect": 1.0, "se": 0.1, "P.Value": 0.001, "adj.P.Val": 0.01, "AveExpr": 5.0},
            {"gene_symbol": "B", "effect": -1.0, "se": 0.1, "P.Value": 0.001, "adj.P.Val": 0.01, "AveExpr": 5.0},
            {"gene_symbol": "C", "effect": 0.1, "se": 0.5, "P.Value": 0.4, "adj.P.Val": 0.5, "AveExpr": 5.0},
        ]
    )
    de_b = _de(
        [
            {"gene_symbol": "A", "effect": 0.2, "se": 1.0, "P.Value": 0.4, "adj.P.Val": 0.8, "AveExpr": 5.0},
            {"gene_symbol": "B", "effect": -0.5, "se": 0.2, "P.Value": 0.01, "adj.P.Val": 0.2, "AveExpr": 5.0},
            {"gene_symbol": "C", "effect": -0.1, "se": 0.5, "P.Value": 0.5, "adj.P.Val": 0.6, "AveExpr": 5.0},
        ]
    )
    out = characterize_meta_significant(meta, de_a, de_b)
    assert out["n_meta_fdr_lt_0.05"] == 2
    assert out["n_also_fdr_lt_0.05_in_study_a"] == 2
    assert out["n_nominal_p_lt_0.05_in_study_b"] == 1
    assert out["n_same_direction_in_study_b"] == 2


def test_loso_and_weight_summary(tmp_path: Path):
    meta = pd.DataFrame(
        {
            "gene_symbol": ["A", "B"],
            "bh_fdr": [0.01, 0.5],
            "pooled_effect": [1.0, 0.1],
            "pooled_se": [0.2, 0.5],
            "i2": [0.0, 0.0],
        }
    )
    de_a = _de(
        [
            {"gene_symbol": "A", "effect": 1.0, "se": 0.1, "P.Value": 0.001, "adj.P.Val": 0.01, "AveExpr": 8.0},
            {"gene_symbol": "B", "effect": 0.1, "se": 0.4, "P.Value": 0.5, "adj.P.Val": 0.8, "AveExpr": 2.0},
        ]
    )
    de_b = _de(
        [
            {"gene_symbol": "A", "effect": 0.5, "se": 0.8, "P.Value": 0.4, "adj.P.Val": 0.9, "AveExpr": 7.0},
            {"gene_symbol": "B", "effect": 0.05, "se": 0.4, "P.Value": 0.6, "adj.P.Val": 0.9, "AveExpr": 2.0},
        ]
    )
    loso = leave_one_study_out_counts(de_a, de_b)
    assert loso["GSE11783_alone"]["bh_fdr_lt_0.05"] == 1
    assert loso["GSE57560_alone"]["bh_fdr_lt_0.05"] == 0
    _, summary = weight_share_summary(meta, de_a, de_b)
    assert summary["n_genes"] == 2
    assert summary["meta_sig_n_genes"] == 1
    assert 0 < summary["mean_weight_share_GSE11783"] < 1


def test_abundance_residual_and_seeds():
    de_a = _de(
        [
            {"gene_symbol": "A", "effect": 2.0, "se": 0.1, "P.Value": 1e-6, "adj.P.Val": 1e-5, "AveExpr": 10.0},
            {"gene_symbol": "B", "effect": -2.0, "se": 0.1, "P.Value": 1e-6, "adj.P.Val": 1e-5, "AveExpr": 1.0},
            {"gene_symbol": "C", "effect": 0.5, "se": 0.2, "P.Value": 0.01, "adj.P.Val": 0.05, "AveExpr": 5.0},
        ]
    )
    de_b = _de(
        [
            {"gene_symbol": "A", "effect": 1.5, "se": 0.3, "P.Value": 0.01, "adj.P.Val": 0.1, "AveExpr": 9.0},
            {"gene_symbol": "B", "effect": -1.5, "se": 0.3, "P.Value": 0.01, "adj.P.Val": 0.1, "AveExpr": 1.5},
            {"gene_symbol": "C", "effect": 0.4, "se": 0.3, "P.Value": 0.1, "adj.P.Val": 0.3, "AveExpr": 5.0},
        ]
    )
    abund = abundance_residual_correlation(de_a, de_b)
    assert abund["n_shared"] == 3
    assert "spearman_before" in abund and "spearman_after_aveexpr_residual" in abund

    meta = pd.DataFrame(
        {
            "gene_symbol": ["A", "B", "C"],
            "pooled_effect": [1.8, -1.8, 0.45],
            "pooled_se": [0.1, 0.1, 0.2],
            "i2": [0.0, 0.0, 0.0],
            "bh_fdr": [0.001, 0.001, 0.05],
        }
    )
    primary = select_seeds(meta, de_a, de_b, require_concordant_direction=False, n=2, i2_max=50.0)
    assert len(primary) == 2
    for col in (
        "effect_GSE11783",
        "se_GSE57560",
        "p_GSE11783",
        "bh_fdr_GSE57560",
        "direction_agree",
        "i2",
        "weight_share_GSE11783",
        "weight_share_GSE57560",
    ):
        assert col in primary.columns
    sens = select_seeds(meta, de_a, de_b, require_concordant_direction=True, n=2, i2_max=50.0)
    assert set(sens["gene_symbol"]).issubset(set(primary["gene_symbol"]))


def test_run_target_diagnostics_writes_artifacts(tmp_path: Path):
    meta = pd.DataFrame(
        {
            "gene_symbol": ["A", "B"],
            "bh_fdr": [0.01, 0.5],
            "pooled_effect": [1.0, 0.1],
            "pooled_se": [0.2, 0.5],
            "i2": [0.0, 0.0],
        }
    )
    de_a = _de(
        [
            {"gene_symbol": "A", "effect": 1.0, "se": 0.1, "P.Value": 0.001, "adj.P.Val": 0.01, "AveExpr": 8.0},
            {"gene_symbol": "B", "effect": 0.1, "se": 0.4, "P.Value": 0.5, "adj.P.Val": 0.8, "AveExpr": 2.0},
        ]
    )
    de_b = _de(
        [
            {"gene_symbol": "A", "effect": 0.5, "se": 0.8, "P.Value": 0.4, "adj.P.Val": 0.9, "AveExpr": 7.0},
            {"gene_symbol": "B", "effect": 0.05, "se": 0.4, "P.Value": 0.6, "adj.P.Val": 0.9, "AveExpr": 2.0},
        ]
    )
    meta_path = tmp_path / "meta.tsv"
    a_path = tmp_path / "a.tsv"
    b_path = tmp_path / "b.tsv"
    meta.to_csv(meta_path, sep="\t", index=False)
    de_a.to_csv(a_path, sep="\t", index=False)
    de_b.to_csv(b_path, sep="\t", index=False)
    out = tmp_path / "out"
    summary = run_target_diagnostics(meta_path, a_path, b_path, out)
    assert (out / "target_meta_diagnostics.json").exists()
    assert (out / "target_re_weight_shares.tsv").exists()
    loaded = json.loads((out / "target_meta_diagnostics.json").read_text())
    assert "plain_statement" in loaded
    assert "single cohort" in loaded["plain_statement"].lower() or "dominated" in loaded["plain_statement"].lower()
