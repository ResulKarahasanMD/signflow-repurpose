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


## 2026-07-23 — Phase 2 build resolutions (pre-results for DE/meta; runtime phenotype lock)

**Outcome-dependent?** No for method choices; phenotype-map strings locked from series-matrix inspection before DE.

| ID | Change | Reason |
|---|---|---|
| P2.1 | GSE238208 phenotype map locked to `disease: hunner_type_interstitial_cystitis` (case) and `disease: mycobacterium bovis bacillus calmette…` (control/BCG) | Prior draft keys (`diagnosis: hunner`) did not match series characteristics; stop-and-report fired correctly |
| P2.2 | Affymetrix/Agilent processing uses author-normalised series matrix; CEL-level RMA not primary in this build | CEL RAW present for GSE11783 but series matrix is the validated author-normalised table used for phenotype-aligned DE |
| P2.3 | Signature bootstrap (≥1000) propagates stratified Hedges' g → DerSimonian–Laird meta; does **not** re-run limma eBayes per resample | Computational resolution; limma remains the observed-sample primary DE artifact; meta identity statistic is Hedges' g |
| P2.4 | Probe collapse rule = max-mean; non-protein-coding rejection via `validate_gene_record` when biotype present | Matches §B2; GPL570 often omits biotype → rejection count may be 0 with note |

| P2.5 | GSE75214 frozen contrast = **UC active colon vs control colon**; CD/Crohn/inactive UC/ileum excluded (109 samples) | Runtime series characteristics; involved-vs-control chosen over involved-vs-uninvolved |


## 2026-07-23 — Phase 2.5 corrections (blocking; pre-Phase 3)

**Outcome-dependent?** Method corrections before accepting Phase 2 headlines; re-freeze after clean-tree re-run.

| ID | Change | Reason |
|---|---|---|
| P2.5.1 | Primary meta = DerSimonian–Laird on **limma moderated logFC + moderated SE**; unmoderated Hedges' g demoted to sensitivity | Prior primary discarded eBayes shrinkage required for n=3 controls |
| P2.5.2 | Concordance restricted to **shared-gene universe**; top-k ranked by moderated t; report k²/N independence expectation; Hedges'g ranking = sensitivity | Prior top-k=0 compared different gene universes / noise-ranked g |
| P2.5.3 | Probe→gene biotype via **HGNC complete set** (cached); never silently keep when GEO lacks biotype | Both platforms lacked biotype → rejected_non_protein_coding=0 was false |
| P2.5.4 | BENCHMARK full DE→meta→concordance→bootstrap→activity; emit `signature_meta_benchmark.tsv` | H1 confirmatory track was phenotype-only |
| P2.5.5 | Record log2 transform fire + min/max; refuse negative clip; vectorize two-study meta; supersede prior `signature_meta.tsv` under `meta/superseded/` | Processing-record gaps |

| P2.5.6 | TARGET unsupervised seed rule frozen **before** seed inspection: top-N=**50** by \|pooled moderated effect\| among genes with I² ≤ **50%** | FDR-based seeding unavailable for TARGET (GSE57560 nearly null at FDR); must not use post-hoc gene list |


## 2026-07-23 — Phase 2.6 (methodological; not outcome-dependent)

**Outcome-dependent?** No. Diagnostics characterise known RE-weighting/I² properties; sensitivity seed freezes a pre-specified concordance guard without changing the primary rule.

| ID | Change | Reason |
|---|---|---|
| P2.6.1 | TARGET meta described as single-cohort–dominated; LOSO + weight-share diagnostics mandatory | GSE57560 n_control=3 → near-zero RE weight; pooled FDR≠replication |
| P2.6.2 | Primary seed rule **unchanged**; sensitivity seed = same ranking + concordant direction in both cohorts (frozen before inspection) | I²→0 when one study is uninformative; primary left intact to avoid outcome-dependent swap |
| P2.6.3 | Seed artifacts must carry per-study effect/SE/p/FDR, direction agreement, I², weight shares | Replication honesty visible at every seed presentation |
| P2.6.4 | Abundance-regressed cross-cohort correlation diagnostic | Distinguish shared technical gradient from disease signature |
| P2.6.5 | Benchmark concordance recorded as positive-control validation of the concordance metric | Makes TARGET near-chance top-k interpretable |
| P3.0 | Edge sign-conflict rule frozen: **retain_as_ambiguous** | Pre-registration before inspecting conflict counts |

| P3.1 | STRING auxiliary score floor frozen at **400** (physical links layer separate) | Pre-inspection; auxiliary undirected layer only |
| P3.2 | Phase 3 CV folds=5; control-graph seed=0; degree-match tolerance=0.15 | Pre-inspection split/control hyperparameters |
