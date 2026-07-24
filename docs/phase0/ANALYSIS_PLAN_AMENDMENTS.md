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


## 2026-07-23 — Phase 3.5 corrections

**Outcome-dependent?** No for A2 (defect fix applying frozen conflict rule to the merged edge set where it was always meant to apply). No for A4 (substrate frozen blind to drug coverage). A1 is methodological (structural absence must not be zero-filled).

| ID | Change | Reason |
|---|---|---|
| P3.5.1 | Apply `retain_as_ambiguous` after OmniPath∪SIGNOR merge on ordered (source,target); assert unique pairs | Defect fix — rule was applied only to OmniPath-internal conflicts |
| P3.5.2 | Redefine CollecTRI/PROGENy gene features with coverage indicators; density QC; variance gate (`feature_min_variance`) | Dead zero-filled columns falsely implied W4 activity fix |
| P3.5.3 | Freeze Phase 5 proximity substrate: primary=`signed_causal_only`, sensitivity=`causal_plus_string_unsigned`; min targets=1; poor coverage → report limitation, do not switch | Blind to drug data (A4) |


## 2026-07-23 — Phase 4

**Outcome-dependent?** No frozen constants changed. H1 W1 rule applied as pre-registered (including negative verdict).

| ID | Change | Reason |
|---|---|---|
| P4.0 | Phase 4 build: baselines + R-GCN/GATv2 + SSL + Mondrian CP + TARGET seed propagation; splits loaded by hash | Execute frozen Phase 4 protocol on corrected Phase 3.5 artifacts |


## 2026-07-24 — Phase 4.5 (deviation-from-protocol repairs; not tuning)

**Outcome-dependent?** No. Undertraining (15 epochs) and single-fold evaluation departed from the pre-registered plan; corrections restore it. Protocol frozen **before** re-run. No hyperparameter search. Original and corrected runs reported side by side. Prior H1 verdict voided as uninformative (rule could not have fired positive under the defective protocol).

| ID | Change | Reason |
|---|---|---|
| P4.5.1 | Freeze GNN training: max_epochs=500, min_epochs=20, patience=50 on val AUPRC, lr=1e-3, wd=1e-4, hidden=32; SSL pretrain_steps=100, finetune matched to scratch | Restore intended training; prior run used 15 full-batch steps |
| P4.5.2 | Evaluate all 5 stratified folds; primary metrics = pooled OOF; also report hub-stratified and LOSO | Plan specifies 5-fold CV; fold 0 only was a protocol deviation |
| P4.5.3 | `phase4_bootstrap_n=1000` on pooled OOF; W1 uses seed-aggregated intervals | Prior n_boot=100 made W1 CI-overlap untrustworthy |
| P4.5.4 | Sign-shuffle permutes `sign` **and** stimulation↔inhibition `relation_type` | Prior shuffle left R-GCN `edge_type` unchanged (inert control) |
| P4.5.5 | Drop conformal `useful` boolean; report mean set size alone | Threshold was analysis-time hardcoded |
| P4.5.6 | Report minimum detectable AUPRC difference beside H1 | Negative from underpowered rule ≠ evidence of absence |


## 2026-07-24 — Phase 5 precondition (§B): UC gold-standard freeze

**Outcome-dependent?** No. Frozen from Open Targets APPROVAL-stage query for `MONDO_0005101` before any drug–target acquisition or Phase 5 scoring. Inclusion is a stage rule plus two documented mechanistic exclusions (non-therapeutics), not a hand-picked rediscovery list.

| ID | Change | Reason |
|---|---|---|
| P5.0.1 | Create `configs/benchmark_drugs.yaml`; hash into `FROZEN_PLAN` (`benchmark_gold_drugs_sha256`) | Close Phase 0 documentation addendum §1 before drug data load |
| P5.0.2 | Out-of-substrate targets: exclude from rediscovery denominator; report as coverage misses | Aligns with `phase5_min_drug_targets_in_substrate=1` and poor-coverage policy |
| P5.0.3 | Eligibility: topical/systemic 5-ASA, corticosteroids, biologics, advanced small molecules included when APPROVAL | Pre-specify class eligibility; not post-hoc after rankings |


## 2026-07-24 — Phase 5.1 freeze (sign-flow policies; before scoring)

**Outcome-dependent?** No. Edge/action/flow constants frozen before drug–target scoring or rankings.

| ID | Change | Reason |
|---|---|---|
| P5.1.1 | `phase5_ambiguous_edge_policy` / `phase5_unsigned_directed_edge_policy` = exclude_from_signed_propagation | Sign-consistent flow cannot invent a sign for ambiguous/unsigned edges |
| P5.1.2 | `phase5_string_unsigned_edge_policy` = unsigned_connectivity_only | Sensitivity substrate adds adjacency without signed path products |
| P5.1.3 | Action sign map + `unknown_action` exclude-from-signed-score | No guessing; retain unknown in coverage |
| P5.1.4 | Flow damping=0.85, max_iter=100, tol=1e-8; aggregation=sum_over_module_opposition | Pre-specify continuous signed diffusion (≠ CARNIVAL ILP) |


## 2026-07-24 — Phase 5.2 ranking universe (pre-rerun) — SUPERSEDED by P5.4

**Outcome-dependent?** No at the time (metric-definition repair). **Superseded:** the OT `drugAndClinicalCandidates` universe remains disease-conditioned on UC and is not a valid rediscovery background under P5.4.

| ID | Change | Reason |
|---|---|---|
| P5.2.1 | Ranking universe = all OT drugAndClinicalCandidates ∪ frozen gold | Gold-only background cannot adjudicate rediscovery |
| P5.2.2 | Rediscovery precision/enrichment computed on full scored universe; report gold ranks explicitly | Fix metric definition defect in first pass |


## 2026-07-24 — Phase 4.6 corrections (blocking; skipped items from prior prompt)

**Outcome-dependent?** No for chance-baseline / sign-finding / conformal-protocol / PU-CI reporting (pre-registered results and protocol repairs). Paired-difference bootstrap is pre-committed **exploratory**; W1 remains binding. Hub/LOSO extension is protocol completeness, not tuning.

| ID | Change | Reason |
|---|---|---|
| P4.6.1 | Chance-baseline narrative: models with CI lower bound >0.5 listed; RWR at chance | Exit report incorrectly claimed no graph method exceeded chance |
| P4.6.2 | Promote sign finding: real≈sign_shuffled but real>rewired (gatv2, ssl_gatv2); MDD≈0.128 caveat into GAP_STATEMENT C1/C4 + Limitations | Topology used; sign not detectably used |
| P4.6.3 | Re-run Mondrian conformal on pooled 5-fold val-calibration protocol | Prior 1.808/0.939 were fold-0-only |
| P4.6.4 | Report PU AUPRC CIs (n_boot=1000); narrate R-GCN PU inversion and ssl_gatv2 PU lead | Point estimates alone hide uncertainty / regime reversal |
| P4.6.5 | Paired-difference bootstrap (gatv2−rwr) as exploratory; W1 binding | Pre-committed exploratory contrast |
| P4.6.6 | Hub-stratified and LOSO extended to rwr and gatv2 (not logistic only) | Protocol completeness |


## 2026-07-24 — Phase 5 VOID (candidate-pool artifact)

**Outcome-dependent?** No. Voiding is a defect recognition about **pool composition**, not a response to an unfavorable rediscovery score. The reported metrics are arithmetically forced by a 29/30 gold pool.

**What was reported (stamp `20260723T235004Z`):** `n_drugs_total=30`, `n_gold_total=29`, `precision_at_20=1.0`, `enrichment_at_20/50=1.0`, `canonical_group_recall_at_50=1.0`, `n_coverage_misses=0`, `n_bh_fdr_survivors=10`, Spearman ρ=−0.137 (p=0.47); `acquisition.skipped=true`; `drug_data_touched=false` while drug tables were read.

**Why void:** With 29/30 gold, any top-20 has precision ≥0.95; enrichment=1.0 is the null; top-50 is the entire pool; coverage cannot be assessed; no ion-channel class in the voided breakdown; W3 degree confounding unanswerable at n=30.

**Manuscript rule:** No rediscovery claim from the voided stamp (or from the disease-conditioned OT-candidate rerun) enters the manuscript.

| ID | Change | Reason |
|---|---|---|
| P5.3.0 | Mark Phase 5 run dirs VOID/SUPERSEDED; retain for audit | Candidate-pool artifact |
| P5.3.1 | Fix `drug_data_touched`: must be true whenever any drug resource is read; stop-and-report + test | Provenance flag inverted reality |


## 2026-07-24 — Phase 5.3 signed-vs-unsigned co-primary freeze

**Outcome-dependent?** Partially exposed, honestly recorded. Motivated by the Phase 4 sign control (pre-registered), **not** by any interpretable drug ranking. The only ranking produced so far is void by pool composition. **Exposure:** 30 drug scores from the voided stamp were nonetheless visible before this freeze. This freeze is therefore not a clean pre-data freeze; the exposure is stated rather than denied.

| ID | Change | Reason |
|---|---|---|
| P5.3.2 | `phase5_score_variant_primary: signed`, `phase5_score_variant_coprimary: unsigned` | Identical directed substrate; unsigned ignores edge sign |
| P5.3.3 | Pre-specify enrichment@20/50 + pool-wide AUROC/AUPRC with paired bootstrap; if signed ≯ unsigned beyond detectable threshold → report no measurable sign benefit | Aligns with Phase 4 control |
| P5.3.4 | Do not drop signed if it underperforms — it remains the pre-registered primary | Finding = the comparison result |
| P5.3.5 | Report connectivity after `exclude_from_signed_propagation` (components, largest component, isolates) for both tracks | ~35% of causal edges removed; fragmentation is a structural finding |


## 2026-07-24 — Phase 5.4 disease-agnostic candidate pool (rebuild)

**Outcome-dependent?** No. Pool inclusion rule frozen **before** acquisition of the approved-drug universe. Motivated by the voided pool-composition defect, not by drug ranks.

| ID | Change | Reason |
|---|---|---|
| P5.4.1 | Pool = ChEMBL `max_phase=4` (disease-agnostic) ∪ frozen UC gold; OT disease candidates annotation-only | Gold must be a small minority (~1–2%) for enrichment to be interpretable |
| P5.4.2 | Freeze pool path + SHA-256; floor `phase5_candidate_pool_min_in_substrate=500`; gold-prevalence ceiling 0.10 (hard-fail) | Prevent re-reporting the voided failure mode |
| P5.4.3 | Metrics: enrichment/recall@k for k∈{20,50,100,250}; AUROC/AUPRC n_boot=1000; chance baselines; signed∥unsigned | Appropriate to corrected prevalence |
| P5.4.4 | Coverage taxonomy includes `ion_channel`; report gold vs full-pool class breakdowns separately | Gap the voided run could not measure |
| P5.4.5 | Stop-and-report: pool below floor; gold prevalence > ceiling; `drug_data_touched` false while resources read; requests to revise frozen constants after rankings | Hard gates |


## 2026-07-24 — Phase 5.4 rebuild completion

**Outcome-dependent?** Results reported as observed under pre-frozen rules; no constants revised after rankings.

| Item | Value |
|---|---|
| Stamp | `20260724T004001Z` |
| Pool | n=4225; gold=29 (0.686%); in-substrate=2833 |
| Signed vs unsigned decision | `sign_contributed_no_measurable_benefit_to_drug_scoring` |
| BH-FDR survivors | **0** |
| Connectivity | 34.3% edges excluded; 1313 components; 1232 isolates |
