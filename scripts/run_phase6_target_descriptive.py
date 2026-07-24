#!/usr/bin/env python3
"""Phase 6 — TARGET descriptive characterisation of the frozen seed module.

Under the recorded Phase 6 rescope (H7 non-informative), the TARGET axes are a
DESCRIPTIVE characterisation of the frozen seed module only. No fused convergence
rank is produced; per-axis results and disagreements are reported, and every table
carries the H7 non-informative status so no reader mistakes a descriptive list for a
candidate list (§B1).

This script consolidates the already-computed TARGET artifacts (H7 futility, control
graphs, axis disagreement, soft-positive ranks from the Phase 5 completions stamp; the
GSE238208 inflammation-specificity defence) plus the axis-availability notes (§B4), and
ENFORCES the do-not rules in code: it hard-fails if any fused-rank artifact is present or
if H7 is not non-informative while a ranked list is being presented as a finding.
"""

from __future__ import annotations

import argparse
import json
from datetime import UTC, datetime
from pathlib import Path

from signflow.phase5 import PHASE5_OUT, ensure_directory


class Phase6GuardError(RuntimeError):
    pass


def _latest_json(glob_pat: str) -> Path | None:
    hits = sorted(PHASE5_OUT.glob(glob_pat))
    return hits[-1] if hits else None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-summary", default=str(PHASE5_OUT / "20260724T075949Z" / "scores" / "target" / "target_eval_summary.json"))
    parser.add_argument("--gse238208", default=None)
    args = parser.parse_args()

    target = json.loads(Path(args.target_summary).read_text())
    h7 = target.get("h7_futility") or {}
    h7_decision = h7.get("decision")

    gse_path = Path(args.gse238208) if args.gse238208 else _latest_json("*/scores/target/gse238208_specificity/gse238208_specificity.json")
    gse = json.loads(gse_path.read_text()) if gse_path and gse_path.exists() else {"status": "missing"}

    # ---- Enforce do-not rules in code (not just prose) ----
    if h7_decision != "non_informative":
        raise Phase6GuardError(f"Expected H7 non_informative under recorded rescope; got {h7_decision!r}")
    if bool(target.get("present_ranked_list_as_finding", False)):
        raise Phase6GuardError("H7 non-informative but a ranked list is flagged as a finding — do-not violation")
    # No fused convergence rank artifact may exist for TARGET.
    fused_hits = list(PHASE5_OUT.glob("*/scores/target/*convergence*")) + list(PHASE5_OUT.glob("*/scores/target/*fused*"))
    if fused_hits:
        raise Phase6GuardError(f"Fused/convergence TARGET rank artifact present — do-not violation: {fused_hits}")

    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    out = ensure_directory(PHASE5_OUT / stamp / "scores" / "target" / "phase6_descriptive")

    seed_primary = (target.get("seed_annotation") or {}).get("primary", {})
    result = {
        "phase": "phase6_target_descriptive",
        "started_at_utc": stamp,
        "mode": "descriptive_characterisation_of_frozen_seed_module_only",
        "h7_status": {
            "decision": h7_decision,
            "mean_pairwise_top20_jaccard": h7.get("mean_pairwise_top_k_jaccard"),
            "threshold_j": h7.get("threshold_j"),
            "banner": "H7 NON-INFORMATIVE — descriptive characterisation only; NOT a candidate list.",
        },
        "fused_convergence_rank_computed": False,
        "seed_module": {
            "n_seeds": seed_primary.get("n_seeds"),
            "n_direction_agree": seed_primary.get("n_direction_agree"),
            "median_i2": seed_primary.get("median_i2"),
            "median_weight_share_GSE11783": seed_primary.get("median_weight_share_GSE11783"),
            "note": "Single-cohort-dominated (GSE11783); GSE57560 near-zero RE weight (P2.6.1).",
        },
        "axes": {
            "inflammation_specificity_GSE238208": {
                "status": "executed",
                "result": gse.get("verdict") if isinstance(gse, dict) else None,
                "primary_auroc": (((gse.get("seed_sets") or {}).get("primary") or {}).get("module_score_auroc_ic_vs_bcg") if isinstance(gse, dict) else None),
                "empirical_p_vs_random": (((gse.get("seed_sets") or {}).get("primary") or {}).get("empirical_p_auroc_ge_random") if isinstance(gse, dict) else None),
                "role": "primary biological check on the module (no TARGET labels).",
            },
            "network_flow_stability_H7": {
                "status": "executed",
                "h7_decision": h7_decision,
                "axis_disagreement": target.get("axis_disagreement"),
                "control_graph_spearman_vs_real": {
                    k: v.get("spearman_vs_real") for k, v in (target.get("control_graphs") or {}).items()
                },
                "note": "Control-graph top-20 Jaccard near 1.0 reflects injection dominance / score ties; interpret via Spearman.",
            },
            "l1000_replication_pmid_39723619": {
                "status": "not_acquired",
                "framing": "Replication of PMID 39723619 (signature reversal); agreement/disagreement table and cell-context mismatch to be stated. No data acquired in this environment -> no imputation.",
            },
            "tissue_cell_appropriateness": {
                "status": "not_acquired",
                "framing": "GTEx/HPA bladder expression + GSE310557 cell-type support. No data acquired -> no imputation.",
            },
            "genetics": {
                "status": "excluded_on_benchmark_A6.2; TARGET has no genetics channel",
            },
            "depmap": {
                "status": "exploratory_only_out_of_domain_disclaimer",
            },
        },
        "soft_positives_at_ranks": (target.get("soft_positives") or {}).get("primary__primary_substrate"),
        "do_not_enforced": [
            "no_fused_convergence_rank",
            "no_added_axes",
            "no_loosened_thresholds",
            "no_ranked_list_presented_as_finding_under_non_informative_H7",
        ],
    }
    (out / "phase6_target_descriptive.json").write_text(json.dumps(result, indent=2, default=str) + "\n")
    print(json.dumps(result, indent=2, default=str))
    print(f"[phase6-target] guards passed; wrote {out / 'phase6_target_descriptive.json'}", flush=True)


if __name__ == "__main__":
    main()
