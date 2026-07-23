# Analysis-plan amendments log

**Rule:** Log every change with date, reason, and whether outcome-dependent. Outcome-dependent amendments after results are discouraged.

---

## 2026-07-23 — Phase 0.5 mandatory amendment (pre-results)

**Outcome-dependent?** No (no pipeline results exist).

| ID | Change | Reason |
|---|---|---|
| A1 | Added `METHODS_PRIOR_ART.md`; rewrote novelty in `GAP_STATEMENT.md` | Methods claims overlapped TxGNN, CARNIVAL, CP literature |
| A2 | Added `LABEL_PROVENANCE.md`; froze OT allowlists + thr≥0.30; TARGET supervised demoted | Drug/literature circularity; IC *n*=38 animal_model only |
| A3 | Mondrian CP; exchangeability = genes/drugs; coverage + mean set size | Unit omitted; patient generalization must not be implied |
| A4 | Per-study DE + random-effects meta-analysis; ComBat sensitivity-only | Study×condition confounding |
| B1–B6 / C1–C7 | Triangulation, L1000 replication, planned cohorts, activity core, stability, hierarchy, ethics | See Phase 0.5 text |

**Claims weakened:** C1/C2/C4/C5 demoted; DepMap demoted; TARGET supervised demoted to exploratory.

---

## 2026-07-23 — Phase 0.6 closure (pre-results)

**Outcome-dependent?** No.

| ID | Change | Reason |
|---|---|---|
| A6.1 | `LABEL_ADJUDICATION.md`; N=15; plausible=6 → **drop animal_model**; UOX out of nodes; seed ban | IMPC labels urinary-tract/renal/stone biased |
| A6.2 | BENCHMARK genetics **excluded from convergence** (Option 1) | Same channel family as training labels |
| A6.3 | External gist receipt; PREREGISTRATION dual-hash clarification | Local tag alone not externally verifiable |
| A6.4 | UC cohorts locked: **GSE75214 + GSE87466** | H1 must not retain cohort degrees of freedom |
| A6.5 | PROGENy/DoRothEA domain-mismatch + activity dispersion QC | Symmetric with DepMap demotion |
| A6.6 | GSE11783 primary = **non-ulcer**; H7 futility **J=0.20**; GSE310557 re-verified | Comparability; descriptive ≠ unbounded |

**Further weakened:** TARGET supervised path **cancelled**. Publishable finding: no usable curated IC/BPS labels under non-drug/non-literature evidence.

---

## 2026-07-23 — Phase 1.5 (pre-Phase 2)

**Outcome-dependent?** No.

| ID | Change | Reason |
|---|---|---|
| 0.1 | Phase 2 download gate cleared after re-verifying public gist receipts | External timestamp required before analysis data download |
| 0.2 | H7 unsupervised-only; mandatory GSE238208 + permuted-signature defences | No TARGET labels after A6.1 |
| A1 | `FROZEN_PLAN.yaml` single source of truth; value-parsing plan-consistency | Substring tests were not real invariants |
| A2–A4 | Append-only manifests; no hardcoded allow-dirty; hard-fail git; plan_hash; device warnings; real GTEx/HPA/SIGNOR checks; seeding; tmp_path tests | Phase 1 defects |
