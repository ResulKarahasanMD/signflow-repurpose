# Phase 5 Exit Report — VOID / REBUILD IN PROGRESS

**Status:** Prior rediscovery results are **VOID — candidate-pool artifact**. Rebuilt disease-agnostic run pending / in progress.

## Voided stamp (exact signature)

| Field | Value |
|---|---|
| Stamp | `20260723T235004Z` |
| `n_drugs_total` | 30 |
| `n_gold_total` | 29 |
| `precision_at_20` | 1.0 |
| `enrichment_at_20` / `@50` | 1.0 / 1.0 |
| `canonical_group_recall_at_50` | 1.0 |
| `n_coverage_misses` / gold misses | 0 / 0 |
| `n_bh_fdr_survivors` | 10 |
| Spearman(score, degree) | ρ=−0.137, p=0.47 |
| `acquisition` | `{manifests: {}, skipped: true}` |
| `drug_data_touched` | `false` (provenance defect; drug tables were read) |

**Why void:** With 29/30 of the pool gold, precision@20 is arithmetic (≥0.95 for any top-20), enrichment=1.0 is the null, top-50 is the entire pool, and coverage cannot be assessed. Target-class breakdown lacked `ion_channel`. No rediscovery claim from this stamp enters the manuscript.

**Also superseded:** `20260723T235350Z` (OT UC `drugAndClinicalCandidates` — still disease-conditioned).

See `VOID.md` under each run directory and `ANALYSIS_PLAN_AMENDMENTS.md` (P5.3 / P5.4).

## Rebuild specification (frozen before acquisition)

- Pool: ChEMBL `max_phase=4` disease-agnostic ∪ frozen UC gold; disease conditioning **forbidden**.
- Floor: ≥500 drugs with ≥1 in-substrate target; gold-prevalence ceiling 10% (hard-fail).
- Scorers: `phase5_score_variant_primary: signed`, `phase5_score_variant_coprimary: unsigned`.
- Metrics: enrichment/recall@k for k∈{20,50,100,250}; AUROC/AUPRC n_boot=1000; chance baselines; signed∥unsigned.
- Connectivity after `exclude_from_signed_propagation` reported for both tracks.

## Gate

Do not begin Phase 6 until §A–§C are committed and this report is replaced by a completed rebuilt exit report.
