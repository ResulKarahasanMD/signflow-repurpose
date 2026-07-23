# Analysis-plan amendments log

**Rule:** Log every change with date, reason, and whether outcome-dependent. Outcome-dependent amendments after results are discouraged.

---

## 2026-07-23 — Phase 0.5 mandatory amendment (pre-results)

**Outcome-dependent?** No (no pipeline results exist).

| ID | Change | Reason |
|---|---|---|
| A1 | Added `METHODS_PRIOR_ART.md`; rewrote novelty in `GAP_STATEMENT.md` | Methods claims overlapped TxGNN, CARNIVAL, CP literature; disease-only lit review was insufficient |
| A2 | Added `LABEL_PROVENANCE.md`; froze OT datatype allowlists + thr≥0.30; TARGET supervised demoted | Drug/literature circularity in unfiltered OT scores; IC filtered to n=38 animal_model only |
| A3 | Mondrian CP; exchangeability unit = genes/drugs; coverage + mean set size; patient bootstrap separate | Prior plan omitted unit; patient generalization must not be implied |
| A4 | Primary path = per-study DE + random-effects meta-analysis; ComBat demoted to sensitivity | Study×condition confounding; 2 studies make fold-wise ComBat non-meaningful |
| B1 | Triangulation: replace DepMap primary axis with GTEx/HPA tissue expression + scRNA cell-type | DepMap is cancer dependency; invalid for IC/BPS and non-cancer UC indication scoring |
| B2 | L1000 framed as replication vs PMID 39723619 | Not novel for IC/BPS; cell-context mismatch stated |
| B3 | GSE238208 → inflammation-specificity control; GSE310557 → cell-type attribution | Pre-specify roles reviewer will demand |
| B4 | DoRothEA/collectTRI + PROGENy via decoupleR are core features/scoring inputs | W4 mRNA≠activity was prior project failure |
| B5 | Signature bootstrap + permuted-signature null are first-class | Small n instability otherwise invisible |
| B6 | Unsupervised comparator shares signature dependence — not orthogonal validation | Prevent overclaim |
| C1 | Disease module = conformal α=0.1 set; sensitivity {50,100,250} | Was illegally deferred |
| C2 | H6: fixed k in {20,50}, primary k=20; bootstrap enrichment test | No test/k before |
| C3 | Hypothesis hierarchy: H1 primary (benchmark); H7 descriptive (target); rest secondary | Avoid fake multi-primary |
| C4 | Label thr pinned (see A2) | — |
| C5 | Power / descriptive TARGET paragraph added | Honest n |
| C6 | Pre-registration via git tag `prereg-phase0.5` + commit hash file | External timestamp before Phase 2 download |
| C7 | Ethics one-liner added | Public de-identified data |

**Claims weakened (intentional):** C1, C2 (foundation novelty), C4 (priority vs CARNIVAL), C5 (technique novelty) demoted to method/rigor components. DepMap demoted. TARGET supervised track demoted to exploratory.
