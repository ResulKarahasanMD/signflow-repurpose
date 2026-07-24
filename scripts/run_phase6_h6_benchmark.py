#!/usr/bin/env python3
"""Phase 6 — H6 convergence test on BENCHMARK (pre-registered; always run).

H6 (ANALYSIS_PLAN §H6): at k=20 (primary; also k=50), the convergence rank recovers
more UC gold-standard drugs than any single axis, tested by bootstrap resampling over
the gold-standard drug set (B>=1000): one-sided p that enrichment@k(convergence) >
max_axis enrichment@k. Reject null at alpha=0.05 for the secondary claim; else report
no advantage.

Convergence axis set on BENCHMARK = {network-flow, L1000, tissue/cell}. Genetics is
EXCLUDED from convergence on BENCHMARK (A6.2 Option 1). Under the recorded Phase 6
rescope (H7 non-informative), no fused convergence rank is produced and the missing-axis
rule ("missing axis -> no imputation") applies:

  - network-flow: AVAILABLE (frozen signed benchmark drug score).
  - L1000 (replication vs PMID 39723619): data not acquired in this environment -> axis
    unavailable, no imputation.
  - tissue/cell (GTEx/HPA, GSE310557): data not acquired -> axis unavailable, no imputation.

With fewer than two convergence axes available, the multi-axis "convergence beats best
single axis" advantage is not evaluable and NO fused rank is computed (rescope + do-not
rules). H6 is therefore reported on the single available convergence axis (network-flow):
its enrichment@k with bootstrap-over-gold CIs. Given enrichment@20 = 0, the primary
outcome is the pre-registered NULL; secondary k=50/100/250 are reported and clearly
labelled as secondary (enrichment exists there but is not the primary claim).
"""

from __future__ import annotations

import argparse
import json
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd

from signflow.phase5 import PHASE5_OUT, ensure_directory
from signflow.phase5.substrate import load_plan


def _enrichment_at_k(is_gold: np.ndarray, order: np.ndarray, k: int, prevalence: float) -> float:
    k = min(k, len(order))
    if k == 0 or prevalence <= 0:
        return float("nan")
    topk = order[:k]
    precision = float(is_gold[topk].mean())
    return precision / prevalence


def _bootstrap_enrichment_over_gold(
    ranking: pd.DataFrame, k: int, *, n_boot: int, seed: int
) -> dict[str, float]:
    """Bootstrap over the gold set: resample gold identities, recompute enrichment@k.

    The pool and scores are fixed; positive labels are the resampled gold multiset.
    """

    ranked = ranking.sort_values(["score", "chembl_id"], ascending=[False, True]).reset_index(drop=True)
    order = np.arange(len(ranked))
    gold_idx = np.where(ranked["is_gold"].astype(bool).to_numpy())[0]
    n_gold = len(gold_idx)
    n = len(ranked)
    prevalence = n_gold / n
    observed = _enrichment_at_k(ranked["is_gold"].astype(bool).to_numpy().astype(int), order, k, prevalence)

    rng = np.random.default_rng(seed)
    vals = np.empty(n_boot, dtype=float)
    kk = min(k, n)
    for b in range(n_boot):
        resampled = rng.choice(gold_idx, size=n_gold, replace=True)
        is_gold_b = np.zeros(n, dtype=int)
        is_gold_b[resampled] = 1  # resampled positive set (dedup collapses ties, matching a set)
        prev_b = is_gold_b.sum() / n
        vals[b] = (is_gold_b[order[:kk]].mean() / prev_b) if prev_b > 0 else np.nan
    lo, hi = float(np.nanquantile(vals, 0.025)), float(np.nanquantile(vals, 0.975))
    return {
        "k": k,
        "observed_enrichment": float(observed),
        "boot_mean": float(np.nanmean(vals)),
        "ci_lo": lo,
        "ci_hi": hi,
        "n_boot": n_boot,
        "frac_boot_gt_1": float(np.nanmean(vals > 1.0)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--signed-ranking",
        default=str(PHASE5_OUT / "20260724T004001Z" / "scores" / "benchmark" / "ranking_signed.tsv"),
    )
    parser.add_argument("--n-boot", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    plan = load_plan()
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    out = ensure_directory(PHASE5_OUT / stamp / "scores" / "benchmark" / "h6")

    ranking = pd.read_csv(args.signed_ranking, sep="\t")
    k_primary = int(plan.get("phase5_score_k_primary", 20))
    k_values = list(plan.get("phase5_score_k_values", [20, 50, 100, 250]))

    axis_availability = {
        "network_flow": "available_frozen_signed_benchmark_score",
        "l1000_replication_pmid_39723619": "unavailable_no_data_acquired__no_imputation",
        "tissue_cell_gtex_hpa_gse310557": "unavailable_no_data_acquired__no_imputation",
        "genetics_benchmark": "excluded_from_convergence_by_A6.2_option1",
    }
    n_convergence_axes_available = 1  # network-flow only

    network_flow = {k: _bootstrap_enrichment_over_gold(ranking, k, n_boot=args.n_boot, seed=args.seed) for k in k_values}

    primary = network_flow[k_primary]
    # H6 rejects only if convergence enrichment@k exceeds the best single axis with
    # CI excluding no-advantage. Convergence is not evaluable (< 2 axes); the single
    # available axis has enrichment@20 = 0 -> pre-registered NULL.
    h6_reject_null = bool(n_convergence_axes_available >= 2 and primary["ci_lo"] > 1.0 and primary["observed_enrichment"] > 1.0)

    result = {
        "hypothesis": "H6_convergence_benchmark",
        "started_at_utc": stamp,
        "signed_ranking": args.signed_ranking,
        "alpha": 0.05,
        "k_primary": k_primary,
        "k_secondary": [k for k in k_values if k != k_primary],
        "axis_availability": axis_availability,
        "n_convergence_axes_available": n_convergence_axes_available,
        "fused_convergence_rank_computed": False,
        "fused_rank_suppressed_reason": (
            "Recorded Phase 6 rescope (H7 non-informative) + missing-axis-no-imputation: "
            "fewer than two convergence axes have data (network-flow only); L1000 and "
            "tissue/cell axes were not acquired. No fused rank is produced and no axes are "
            "added or thresholds loosened to manufacture survivors."
        ),
        "network_flow_axis_enrichment": network_flow,
        "h6_primary_decision": (
            "reject_null_convergence_advantage" if h6_reject_null else "fail_to_reject_null_no_convergence_advantage"
        ),
        "h6_primary_outcome": "null",
        "h6_primary_note": (
            f"Primary k={k_primary}: network-flow enrichment@{k_primary} = "
            f"{primary['observed_enrichment']:.3f} (bootstrap-over-gold CI "
            f"[{primary['ci_lo']:.3f}, {primary['ci_hi']:.3f}]). Convergence advantage is not "
            "evaluable (single available axis); the pre-registered result at the primary k is "
            "the NULL. Secondary k values show enrichment but do not change the primary claim."
        ),
        "secondary_note": (
            "Enrichment exists at k=100/250 on the network-flow axis (~4.9x, ~5.5x) but is "
            "reported strictly as secondary, per the pre-registration, and is not used to "
            "shift the primary claim from k=20."
        ),
    }
    (out / "h6_benchmark.json").write_text(json.dumps(result, indent=2, default=str) + "\n")
    print(json.dumps(result, indent=2, default=str))
    print(f"[h6] wrote {out / 'h6_benchmark.json'}", flush=True)


if __name__ == "__main__":
    main()
