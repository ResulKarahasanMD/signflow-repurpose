# Pre-registered analysis plan — SignFlow-Repurpose

**Version:** 0.2.0-phase0.5  
**Date locked:** 2026-07-23  
**Supersedes:** 0.1.0-phase0  
**Amendments:** [`ANALYSIS_PLAN_AMENDMENTS.md`](ANALYSIS_PLAN_AMENDMENTS.md)  
**Label rules:** [`LABEL_PROVENANCE.md`](LABEL_PROVENANCE.md)  
**Methods prior art:** [`METHODS_PRIOR_ART.md`](METHODS_PRIOR_ART.md)

**Amendment rule:** Changes after Phase 2 data freeze must be logged with date, reason, and outcome-dependence (discouraged).

This plan is pre-specified before Phase 2 data download and model fitting.

---

## 0. Tracks, hierarchy, and configuration

| Config key | Default | Notes |
|---|---|---|
| `BENCHMARK_DISEASE_ID` | `MONDO_0005101` | Ulcerative colitis |
| `TARGET_DISEASE_ID` | `MONDO_0018301` | IC/BPS |
| `TARGET_PRIMARY_COHORTS` | `GSE11783`, `GSE57560` | Bulk DE + meta-analysis |
| `TARGET_INFLAMMATION_CONTROL` | `GSE238208` | HIC vs BCG — not healthy DE |
| `TARGET_SCRNA` | `GSE310557` | Cell-type attribution only |
| `LABEL_SCORE_THRESHOLD` | `0.30` | See LABEL_PROVENANCE |
| `RANDOM_SEEDS` | `{0,1,2,3,4}` | Multi-seed mean ± sd |
| `CONFORMAL_ALPHA` | `0.1` | Sensitivities `{0.05, 0.2}` |
| `NULL_PERMUTATIONS` | `≥1000` | Degree-matched drug nulls |
| `SIGNATURE_BOOTSTRAPS` | `≥1000` when feasible | Else max leave-one-out style with documented *B* |
| `TRIANGULATION_K_PRIMARY` | `20` | Also report `k=50` |
| `PYTHON` | `3.12` via `uv` | — |

Device: `cuda → mps → cpu`; `PYTORCH_ENABLE_MPS_FALLBACK=1`; pure-Python PyG convs.

### Hypothesis hierarchy (C3)

| Rank | Hypothesis | Track | Role |
|---|---|---|---|
| **Primary** | **H1** | BENCHMARK | Confirmatory methods endpoint |
| **Descriptive primary (application)** | **H7** | TARGET | Hypothesis-generating; not confirmatory |
| Secondary | H2–H6 | as applicable | Hypothesis-generating; **not** multiplicity-corrected across hypotheses |

Secondary results are reported as exploratory.

---

## 1. Hypotheses

### H1 — PRIMARY (BENCHMARK)
On UC, with **filtered genetics labels** (§2.3), a relational/signed GNN (± SSL pretrain) **improves** gene prioritization vs RWR on the pre-specified primary split (AUPRC primary, AUROC co-reported), **or** we report that it does not.

### H2 — secondary (transfer)
SSL pretrain + few-shot fine-tune improves vs from-scratch under matched budgets (BENCHMARK; TARGET exploratory only).

### H3 — secondary (confounding)
After degree-matched nulls and degree covariates, top drug ranks are not sole degree artifacts; report Spearman ρ(score, degree) and partial correlations.

### H4 — secondary (directionality)
Sign-consistent signaling-flow (± activity features) enriches UC gold-standard drugs vs unsigned/undirected ablations.

### H5 — secondary (conformal)
Mondrian conformal gene/drug sets achieve empirical coverage ≥ 1−α **and** we report mean set size; uselessly large sets = negative result.

### H6 — secondary (triangulation)
At **k=20** (also k=50), convergence rank recovers more UC gold-standard drugs than any single axis, tested by **bootstrap resampling over the gold-standard drug set** (*B*≥1000): one-sided *p* that enrichment@k(convergence) > max_axis enrichment@k. Decision: reject null at α=0.05 for secondary claim; if not, report no advantage.

### H7 — DESCRIPTIVE (TARGET)
IC/BPS yields a hypothesis-generating conformal (or unsupervised-module) drug set with transparent axis disagreements; soft positives reported where they land — not forced.

---

## 2. Data

### 2.1 Expression cohorts and roles

| Accession | Role |
|---|---|
| GSE11783, GSE57560 | **Primary** TARGET bulk DE (independent per study) + meta-analysis |
| GSE238208 | **Inflammation-specificity control:** does the prioritized module/signature separate Hunner-type IC from BCG cystitis? **Forbidden** as IC-vs-healthy DE |
| GSE310557 | **Cell-type attribution** for module genes / target expression axis; **not** merged into bulk DE |
| GSE28242 | Optional urine orthogonal signature only |
| UC ≥2 mucosal cohorts | Locked in Phase 2 after phenotype QC (candidates verified exist: GSE75214, G87466, GSE59071, GSE107499, GSE128682) |

### 2.2 Disease signature — primary path (A4)

1. **DE within each study independently** (limma / limma-voom for microarray; PyDESeq2 only for raw counts).  
2. GSE11783: **one sample per patient** (primary = ulcer site; sensitivity = non-ulcer; patient-mean).  
3. Combine studies by **random-effects meta-analysis of standardised effect sizes**; report pooled effect, SE, *p*, BH-FDR, **I²** and Cochran’s **Q**.  
4. **Concordance analysis** between studies: effect-size correlation, directional agreement, top-*k* overlap. Low concordance is a publishable finding — do not average it away.  
5. **Sensitivity:** rank-based per-study standardisation; pooled ComBat **only** as clearly labelled sensitivity (never primary).  
6. PCA/silhouette = descriptive QC only (“not used for model selection”).

Signature = **feature / scoring target**, never labels.

### 2.3 Labels (A2) — see LABEL_PROVENANCE.md

- **BENCHMARK primary labels:** Open Targets `genetic_association` datatype score ≥ **0.30** (`label_score = max` allowlisted datatype scores). *n*=396 (2026-07-23).  
- **TARGET:** filtered `animal_model` ≥0.30 yields *n*=38 and **no genetics channel** → **supervised TARGET is secondary/exploratory**; **primary TARGET analysis = unsupervised comparator**.  
- Sensitivities: thr ∈ {0.20, 0.40}; UC ± `somatic_mutation`.  
- **Forbidden:** DE-as-labels; `clinical`/`clinical_precedence`/known-drug evidence when drug rediscovery is an endpoint; primary use of `literature`/`europepmc`/`genetic_literature`; `rna_expression`/`expression_atlas` as labels when expression is a feature/signature.

### 2.4 Networks and **core** activity inference (B4)

- OmniPath (`is_directed`, `is_stimulation`, `is_inhibition`) + SIGNOR; STRING v12 auxiliary.  
- **Core (not optional):** TF activity via **DoRothEA / CollecTRI** and pathway activity via **PROGENy**, both through **decoupleR** (versions recorded in run manifest). Cite CARNIVAL/DoRothEA/PROGENy as established footprint approaches; our flow score differs by continuous signed propagation + degree-matched nulls + conformal wrapping (not ILP causal network inference).  
- If a resource is unavailable/unlicensed → **stop that component and report**; do not silently fall back to raw mRNA-only direction scoring.

### 2.5 Drugs and triangulation inputs

- DGIdb GraphQL; ChEMBL; Open Targets drug/clinical candidates (for soft/hard positives and drug–target actions — **not** for gene labels).  
- DrugBank only if user-supplied XML.  
- L1000 / Connectivity Map: **orthogonal replication** of PMID 39723619; compare hit agreement/disagreement; state cancer-line context mismatch (B2).  
- Tissue appropriateness: **GTEx and/or Human Protein Atlas** bladder (TARGET) / colon (BENCHMARK) expression specificity (B1).  
- Cell-type support: GSE310557.  
- DepMap: **exploratory toxicity-proxy only**, labelled out-of-domain for non-cancer indications — not a triangulation axis.  
- Genetics triangulation axis: filtered OT genetics evidence on **drug targets** (independent of gene-label allowlist used for training where applicable).

### 2.6 Provenance and ethics

- Each run: `run_manifest.json` (config hash, package versions, seeds, input SHA-256, source versions, access dates).  
- **Ethics:** All analysis data are public, de-identified, previously published; no IRB approval required for this secondary analysis; no patient-level raw data are committed to the repository.

### 2.7 Power statement (C5)

Effective TARGET bulk sizes are small (GSE11783: 5 IC patients vs 6 controls after one-sample-per-patient; GSE57560: 13 IC vs 3 controls). Under a two-sample *t*-test approximation with α=0.05 two-sided, detecting a large standardized effect (Cohen’s *d*≈1.5) at 80% power already requires on the order of ~8–10 per group; with 3 controls, only very large effects are detectable, and meta-analytic precision remains limited. **TARGET results are descriptive / hypothesis-generating.** BENCHMARK UC cohorts (hundreds of arrays in candidate series) support confirmatory H1. Cohort size **precludes a state-of-the-art methods claim**.

---

## 3. Graph construction

Nodes: protein-coding genes in cohort features ∩ network.  
Edges: directed; signed; relation-typed; confidence retained.  
Mapping QC tables mandatory.

---

## 4. Models

### 4.1 Baselines
Random / degree / PageRank; logistic regression; HistGradientBoosting; feature-only MLP; label propagation; **RWR**; controls: degree-preserving rewire, edge shuffle, **sign shuffle**, features-only.

### 4.2 Proposed model
Relational/signed MP (R-GCN-style and/or GATv2 with edge-type/sign); SSL pretrain (DGI and/or contrastive and/or masked attr + edge-sign pred); few-shot fine-tune; PU framing.

### 4.3 Splits
Stratified CV; degree-aware split; BENCHMARK leave-one-study-out inductive. TARGET supervised exploratory only if labels used.

### 4.4 Metrics
AUROC, AUPRC, MCC, balanced accuracy, ECE, precision@k, recall@k, enrichment@k; bootstrap CIs; multi-seed mean ± sd.  
**H1 primary metric:** AUPRC (PU-aware) + AUROC.  
**W1 honesty:** if RWR AUPRC ≥ GNN within bootstrap CI overlap on primary split → conclude no GNN advantage.

---

## 5. Conformal prediction (A3)

| Stage | Exchangeability unit | Method |
|---|---|---|
| Gene prioritization | **Genes** (labelled calibration genes) | **Mondrian / class-conditional CP** (separate calibration by class; PU-aware protocol documented) |
| Drug ranking | **Drugs** | Mondrian CP on drug scores |

**Non-claim (verbatim for Methods):**  
*Conformal coverage is over the gene/drug population under the labelled-data distribution. It does not quantify sampling uncertainty over patients, and does not imply the disease signature generalises to new IC/BPS cohorts.*

Also report: empirical coverage **and mean set size**.  
Sensitivity: weighted conformal under assumed covariate shift.  
**Patient-level uncertainty:** bootstrap over samples for signature/drug-rank stability (B5) — reported separately from conformal coverage.  
Citations: PMID 33075380; Mondrian PMID 28628322; Norinder PMID 24797111.

α = 0.1 primary; sensitivities 0.05 and 0.2.

---

## 6. Drug scoring (C4) and disease module (C1)

1. **Primary disease module** = Mondrian conformal gene set at α=0.1 (BENCHMARK supervised path). For TARGET primary unsupervised path: module = top fixed-size seeds from signature extremes after unsupervised ranking, with sizes in sensitivity set.  
2. **Sensitivity module sizes:** `{50, 100, 250}` genes.  
3. Signed drug action → signaling-flow / sign-consistent proximity opposing disease signature/activity footprint.  
4. Degree-matched null ≥1000; z, empirical *p*, BH-FDR; score–degree correlation.  
5. Soft/hard positives reported at ranks (not forced).

Differentiate from CARNIVAL: continuous flow score + nulls + conformal sets + triangulation; not integer linear causal subnetwork inference.

---

## 7. Triangulation (B1–B2) — axis table

| Axis | Role | Definition |
|---|---|---|
| Network-flow | Primary scoring | Sign-consistent drug score (C4) |
| L1000 | **Replication** of PMID 39723619 | Signature reversal; agreement table vs prior compounds; cell-context limitation stated |
| Tissue/cell appropriateness | Replaces DepMap | GTEx/HPA expression in bladder or colon; GSE310557 cell-type support |
| Genetics | Orthogonal | Filtered OT genetics support for drug targets |
| DepMap | Exploratory only | Cancer dependency; out-of-domain disclaimer |

**Convergence rank:** reciprocal rank fusion over axes with non-missing scores. Disagreements listed. Missing axis → no imputation.

---

## 8. Unsupervised comparator (B6; TARGET primary)

SSL embeddings + propagation from **signature-extreme seeds** (seeds ≠ curated labels).  

**Shared-dependence non-claim:** This comparator and the drug score both use the disease signature; agreement with the supervised track is **not** orthogonal validation.

---

## 9. Validation and stability (B3, B5)

- Positives / negatives / ablations as in Phase 0 plan (signed vs unsigned; SSL vs scratch; etc.).  
- **GSE238208:** inflammation-specificity of module/signature (IC vs BCG).  
- **GSE310557:** cell-type attribution.  
- **Signature bootstrap (≥1000 when feasible):** recompute drug ranking; top-*k* Jaccard; rank correlation; per-drug rank IQR.  
- **Permuted-signature null:** shuffle case/control labels, rerun drug scoring; top hits must collapse — else pipeline reports topology, not disease biology.  
- Leave-one-dataset-out; top-*k* stability.

---

## 10. Reporting rules

1. Real numbers only, including negatives and uselessly large conformal sets.  
2. No clinical recommendations.  
3. Stop-and-report if a core resource is missing.  
4. No SOTA methods claim given TARGET *n*.  
5. Novelty as combination + controlled benchmark + application — not priority over TxGNN/CARNIVAL/CP.

---

## 11. Phase exit criteria

| Phase | Exit when |
|---|---|
| 0.5 | Amendments + label provenance + methods prior art + prereg timestamp committed |
| 1 | Env + config + API harness smoke tests green |
| 2 | Per-study DE + meta-analysis + concordance artifacts |
| 3–7 | As previously scoped with amended designs |

---

## 12. Commitment

Follow this plan unless amended in `ANALYSIS_PLAN_AMENDMENTS.md`. RWR may win. Label filters will not be relaxed post-hoc to manufacture TARGET supervised power.
