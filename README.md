# SignFlow-Repurpose

Controlled signed-network prioritization and multi-evidence triangulation for **IC/BPS** (understudied application) with **ulcerative colitis** (methods benchmark).

> **Hypothesis generation only.** Scores are not binding, PK/PD, causality, or clinical efficacy. **No treatment recommendations.**

---

## Phase status

| Phase | Status |
|---|---|
| **0 — Gap, data gate, plan** | Complete |
| **0.5 — Mandatory amendments** | **Complete** (2026-07-23) |
| 1 — Environment & provenance harness | Not started (blocked until 0.5 committed) |

**Data gate:** PASS (conditional) — GSE11783 + GSE57560.  
**Label gate:** UC genetics labels *n*=396; IC filtered labels *n*=38 → TARGET **unsupervised primary**.

---

## Phase 0 / 0.5 documents

| Doc | Role |
|---|---|
| [`docs/phase0/GAP_STATEMENT.md`](docs/phase0/GAP_STATEMENT.md) | Re-scoped novelty (not TxGNN/CARNIVAL priority) |
| [`docs/phase0/METHODS_PRIOR_ART.md`](docs/phase0/METHODS_PRIOR_ART.md) | Per-contribution differentiation |
| [`docs/phase0/LABEL_PROVENANCE.md`](docs/phase0/LABEL_PROVENANCE.md) | Frozen OT allowlists + thr≥0.30 |
| [`docs/phase0/DATA_AVAILABILITY.md`](docs/phase0/DATA_AVAILABILITY.md) | Cohort gate + planned control roles |
| [`docs/phase0/ANALYSIS_PLAN.md`](docs/phase0/ANALYSIS_PLAN.md) | Pre-registered plan v0.2.0 |
| [`docs/phase0/ANALYSIS_PLAN_AMENDMENTS.md`](docs/phase0/ANALYSIS_PLAN_AMENDMENTS.md) | Change log |
| [`docs/phase0/LITERATURE_REVIEW.md`](docs/phase0/LITERATURE_REVIEW.md) | Disease-specific lit |
| [`docs/phase0/PREREGISTRATION.md`](docs/phase0/PREREGISTRATION.md) | Timestamp / git tag |

---

## Defaults

| Track | ID | Primary analysis |
|---|---|---|
| BENCHMARK | `MONDO_0005101` | Supervised H1 (filtered genetics labels) |
| TARGET | `MONDO_0018301` | Unsupervised comparator (H7 descriptive) |

Primary TARGET bulk cohorts: `GSE11783`, `GSE57560`. Planned: `GSE238208` (inflammation specificity), `GSE310557` (cell type).
