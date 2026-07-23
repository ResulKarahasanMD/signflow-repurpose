# Pre-registered analysis plan — SignFlow-Repurpose

**Version:** 0.3.0-phase0.6  
**Date locked:** 2026-07-23  
**Supersedes:** 0.2.0-phase0.5  
**Amendments:** [`ANALYSIS_PLAN_AMENDMENTS.md`](ANALYSIS_PLAN_AMENDMENTS.md)  
**Label rules:** [`LABEL_PROVENANCE.md`](LABEL_PROVENANCE.md) · [`LABEL_ADJUDICATION.md`](LABEL_ADJUDICATION.md)  
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
| `LABEL_ADJUDICATION_N` | `15` | Min `plausible-ICBPS`; see LABEL_ADJUDICATION |
| `H7_FUTILITY_J` | `0.20` | Signature-bootstrap top-20 Jaccard futility |
| `BENCHMARK_COHORTS` | `GSE75214`, `GSE87466` | Locked by §A6.4 rule (below) |
| `GSE11783_PRIMARY_SITE` | `non-ulcer` | Ulcer = sensitivity (A6.6) |
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
SSL pretrain + few-shot fine-tune improves vs from-scratch under matched budgets (**BENCHMARK only**; TARGET supervised path cancelled after A6.1).

### H3 — secondary (confounding)
After degree-matched nulls and degree covariates, top drug ranks are not sole degree artifacts; report Spearman ρ(score, degree) and partial correlations.

### H4 — secondary (directionality)
Sign-consistent signaling-flow (± activity features) enriches UC gold-standard drugs vs unsigned/undirected ablations.

### H5 — secondary (conformal)
Mondrian conformal gene/drug sets achieve empirical coverage ≥ 1−α **and** we report mean set size; uselessly large sets = negative result.

### H6 — secondary (triangulation)
At **k=20** (also k=50), convergence rank recovers more UC gold-standard drugs than any single axis, tested by **bootstrap resampling over the gold-standard drug set** (*B*≥1000): one-sided *p* that enrichment@k(convergence) > max_axis enrichment@k. Decision: reject null at α=0.05 for secondary claim; if not, report no advantage.

### H7 — DESCRIPTIVE (TARGET) — unsupervised only
IC/BPS has **no usable curated disease-gene labels** after A6.1 (`plausible-ICBPS`=6 < N=15; `animal_model` dropped). H7 is therefore **unsupervised-only**: signature-extreme seeds → module → drug scoring → triangulation, with soft positives reported where they land — not forced.  
**Futility (A6.6):** if signature-bootstrap top-20 drug Jaccard < **J = 0.20**, report TARGET as **non-informative** (stability stats + disagreements only; **no ranked candidate list as a finding**).  
**Mandatory defences (Phase 1.5):** with no supervised anchor, **GSE238208 inflammation-specificity** and the **permuted-signature null** are the primary checks against reporting graph topology as disease biology. Failure of either is a **reportable negative result** — not optional.

---

## 2. Data

### 2.1 Expression cohorts and roles

| Accession | Role |
|---|---|
| GSE11783, GSE57560 | **Primary** TARGET bulk DE (independent per study) + meta-analysis |
| GSE238208 | **Inflammation-specificity control** (HIC vs BCG); **forbidden** as IC-vs-healthy DE |
| GSE310557 | **Cell-type attribution** (re-verified: *n*=17 = 10 HIC + 7 PSUI; GPL34284; PMID 41457077); not merged into bulk DE |
| GSE28242 | Optional urine orthogonal signature only |
| **GSE75214, GSE87466** | **BENCHMARK UC cohorts (locked by A6.4 rule)** |

**UC lock rule (A6.4, frozen):** Among `{GSE75214, GSE87466, GSE59071, GSE107499, GSE128682}`, keep series with series-matrix HEAD `Content-Length > 100000` and ≥1 control/healthy/normal sample; pick the two with largest control counts; ties → ascending GSE. **Expected selection (2026-07-23):** GSE75214 (22 controls) + GSE87466 (21). GSE107499 has 0 controls; GSE128682 matrix ~4 KB stub → ineligible.

### 2.2 Disease signature — primary path (A4)

1. **DE within each study independently** (limma / limma-voom for microarray; PyDESeq2 only for raw counts).  
2. GSE11783: **one sample per patient**. **Primary site = non-ulcer (“ni”)** for cross-cohort comparability with non-Hunner-selected GSE57560; **ulcer site = sensitivity**; patient-mean = additional sensitivity. Rationale (A6.6): ulcer-vs-normal maximises granulation/immune-infiltrate divergence from GSE57560 and would bias the meta-analysis toward lesion necrosis rather than shared IC/BPS biology; non-ulcer primary reduces that divergence. This choice **increases** the importance of GSE238208 inflammation-specificity control.  
3. Combine studies by **random-effects meta-analysis of standardised effect sizes**; report pooled effect, SE, *p*, BH-FDR, **I²** and Cochran’s **Q**.  
4. **Concordance analysis** between studies: effect-size correlation, directional agreement, top-*k* overlap. Low concordance is a publishable finding — do not average it away.  
5. **Sensitivity:** rank-based per-study standardisation; pooled ComBat **only** as clearly labelled sensitivity (never primary).  
6. PCA/silhouette = descriptive QC only (“not used for model selection”).

Signature = **feature / scoring target**, never labels.

### 2.3 Labels (A2) — see LABEL_PROVENANCE.md

- **BENCHMARK primary labels:** Open Targets `genetic_association` datatype score ≥ **0.30**. *n*=396 (2026-07-23). Sensitivities: thr ∈ {0.20, 0.40}; ± `somatic_mutation`.  
- **TARGET labels:** After IMPC adjudication ([`LABEL_ADJUDICATION.md`](LABEL_ADJUDICATION.md)), `plausible-ICBPS` = 6 < **N=15** → **`animal_model` channel DROPPED**. Finding: *IC/BPS has no usable curated disease-gene label set under non-drug, non-literature evidence.* TARGET supervised path **cancelled**.  
- **Forbidden:** DE-as-labels; drug/clinical evidence labels when rediscovery is an endpoint; literature labels in primary definition; expression_atlas labels when expression is a feature; non-`protein_coding` genes (e.g. **UOX**) as labels or nodes; adjudicated animal_model genes as unsupervised seeds.

### 2.4 Networks and **core** activity inference (B4)

- OmniPath (`is_directed`, `is_stimulation`, `is_inhibition`) + SIGNOR; STRING v12 auxiliary.  
- **Core (not optional):** TF activity via **DoRothEA / CollecTRI** and pathway activity via **PROGENy**, both through **decoupleR** (versions recorded in run manifest). Cite CARNIVAL/DoRothEA/PROGENy as established footprint approaches; our flow score differs by continuous signed propagation + degree-matched nulls + conformal wrapping (not ILP causal network inference).  
- **Domain-mismatch limitation (A6.5, symmetric with DepMap):** PROGENy footprints and DoRothEA/CollecTRI regulons are substantially informed by cancer/cell-line perturbation biology. They remain the best open activity-inference tools and are retained, but are **out-of-domain relative to benign inflammatory bladder/colon disease** — stated in Limitations.  
- **Precision limitation:** TF/pathway activity from *n*≈16 microarray cohorts is low-precision; report **per-sample activity score dispersion** as QC.  
- If a resource is unavailable/unlicensed → **stop that component and report**; do not silently fall back to raw mRNA-only direction scoring.

### 2.5 Drugs and triangulation inputs

- DGIdb GraphQL; ChEMBL; Open Targets drug/clinical candidates (for soft/hard positives and drug–target actions — **not** for gene labels).  
- DrugBank only if user-supplied XML.  
- L1000 / Connectivity Map: **orthogonal replication** of PMID 39723619; compare hit agreement/disagreement; state cancer-line context mismatch (B2).  
- Tissue appropriateness: **GTEx and/or Human Protein Atlas** bladder (TARGET) / colon (BENCHMARK) expression specificity (B1).  
- Cell-type support: GSE310557.  
- DepMap: **exploratory toxicity-proxy only**, labelled out-of-domain for non-cancer indications — not a triangulation axis.

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
Relational/signed MP (R-GCN-style and/or GATv2 with edge-type/sign); SSL pretrain (DGI and/or contrastive and/or masked attr + edge-sign pred); few-shot fine-tune on **BENCHMARK** labels only; PU framing. TARGET has no curated labels after A6.1.

### 4.3 Splits
Stratified CV; degree-aware split; BENCHMARK leave-one-study-out inductive.

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
| Network-flow | In convergence | Sign-consistent drug score (C4) |
| L1000 | In convergence; **replication** of PMID 39723619 | Signature reversal; agreement/disagreement table; cell-context limitation |
| Tissue/cell appropriateness | In convergence | GTEx/HPA bladder or colon expression; GSE310557 cell-type support |
| Genetics (BENCHMARK) | **Descriptive only — EXCLUDED from convergence (A6.2 Option 1)** | Filtered OT genetics on drug targets; same channel family as training labels → would make H6 partly self-fulfilling. Methods note: genetic support predicts drug-development success by construction, not by error. TARGET unaffected (no genetics channel). |
| DepMap | Exploratory only | Cancer dependency; out-of-domain disclaimer |

**Convergence rank (BENCHMARK):** reciprocal rank fusion over {network-flow, L1000, tissue/cell} only. Genetics reported separately, never folded into H6. Disagreements listed. Missing axis → no imputation.

---

## 8. Unsupervised comparator (B6; TARGET primary — only TARGET path)

SSL embeddings + propagation from **signature-extreme seeds only**.  

**Hard constraint (A6.1):** Open Targets / IMPC / adjudicated animal_model genes — including the former 38-gene list — **must never be used as seeds** for this path.  

**Shared-dependence non-claim:** This comparator and the drug score both use the disease signature; agreement with any supervised BENCHMARK analysis is **not** orthogonal validation of TARGET.

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
