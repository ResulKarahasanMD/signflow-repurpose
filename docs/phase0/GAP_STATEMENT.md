# Gap statement — SignFlow-Repurpose (Phase 0.5)

**Project working title:** *SignFlow-Repurpose*: controlled signed-network prioritization and multi-evidence triangulation for IC/BPS (application) with ulcerative colitis (methods benchmark).

**Date:** 2026-07-24 (UTC)  
**Status:** Phase 5 completions (§A) done; TARGET H7 **non-informative**; Phase 6 rescope recorded; paper identity reframed (P5.5.C).  
**See:** [`DATA_AVAILABILITY.md`](DATA_AVAILABILITY.md) · [`METHODS_PRIOR_ART.md`](METHODS_PRIOR_ART.md) · [`LABEL_PROVENANCE.md`](LABEL_PROVENANCE.md) · [`ANALYSIS_PLAN.md`](ANALYSIS_PLAN.md) · [`PHASE6_RESCOPE_DECISION.json`](PHASE6_RESCOPE_DECISION.json)

---

## 1. Clinical / therapeutic gap (TARGET: IC/BPS)

Interstitial cystitis / bladder pain syndrome (IC/BPS; Open Targets `MONDO_0018301`) lacks a consensus disease-modifying therapy.

**Verified evidence:**
- NMA of 70 RCTs (*n*=3,651): no treatment improved **all** outcomes (PMID [38658949](https://pubmed.ncbi.nlm.nih.gov/38658949/); DOI `10.1186/s12894-024-01485-w`).
- Pharmacotherapy NMA: symptomatic/partial signals without DMT standard (PMID [33638677](https://pubmed.ncbi.nlm.nih.gov/33638677/)).
- Open Targets (2026-07-23): 49 drug/clinical candidates; sole APPROVAL-stage = pentosan polysulfate sodium (`CHEMBL1201516`).
- ClinicalTrials.gov: 238 studies for interstitial cystitis (API v2, 2026-07-23).

**Interpretation:** Appropriate understudied application for **hypothesis generation only** — not a clinical recommendation engine.

**Power limitation (plain statement):** Primary TARGET bulk cohorts are tiny (after one-sample-per-patient: ~5 vs 6 in GSE11783; 13 vs 3 in GSE57560). This **precludes a state-of-the-art methods claim**. TARGET outputs are descriptive / hypothesis-generating.

---

## 2. Novelty — methods prior art + IC/BPS application (re-scoped)

### 2.1 Disease-specific computational priors (still true)

Prior IC/BPS *in silico* work (e.g. PMID [39723619](https://pubmed.ncbi.nlm.nih.gov/39723619/) L1000/cMAP reversal; PMID [36619718](https://pubmed.ncbi.nlm.nih.gov/36619718/); PMID [40777016](https://pubmed.ncbi.nlm.nih.gov/40777016/)) does not implement our **controlled** stack (filtered labels, meta-analytic dual-cohort signature, Mondrian CP with set size, inflammation-specificity + scRNA attribution, honest GNN-vs-RWR). L1000 here is **replication**, not novelty.

### 2.2 Methods priors that block priority claims (mandatory)

| Prior | PMID / DOI | Blocks claiming… |
|---|---|---|
| **TxGNN** | [39322717](https://pubmed.ncbi.nlm.nih.gov/39322717/); `10.1038/s41591-024-03233-x` | Priority on GNN transfer / zero-shot repurposing for diseases with limited treatments |
| **CARNIVAL** | [31728204](https://pubmed.ncbi.nlm.nih.gov/31728204/); `10.1038/s41540-019-0118-z` | Priority on signed directed OmniPath causal reasoning with signed drug effects |
| **CP in drug discovery** | [33075380](https://pubmed.ncbi.nlm.nih.gov/33075380/) | Novelty of conformal prediction as a technique |
| DoRothEA / PROGENy / decoupleR | [31340985](https://pubmed.ncbi.nlm.nih.gov/31340985/), [29295995](https://pubmed.ncbi.nlm.nih.gov/29295995/), [36699385](https://pubmed.ncbi.nlm.nih.gov/36699385/) | Novelty of footprint activity inference |
| Network proximity | [26831545](https://pubmed.ncbi.nlm.nih.gov/26831545/), [30002366](https://pubmed.ncbi.nlm.nih.gov/30002366/) | Novelty of proximity-based repurposing |
| PrimeKG | [36732524](https://pubmed.ncbi.nlm.nih.gov/36732524/) | Novelty of multi-relational biomedical KGs |

Full differentiation table: [`METHODS_PRIOR_ART.md`](METHODS_PRIOR_ART.md).

### 2.3 Authoritative novelty claim (combination + application + controlled benchmark)

**We claim:** a **pre-registered, leakage-aware combination** of established components (signed OmniPath/SIGNOR graph; footprint activities; relational GNN; RWR baseline; Mondrian CP; multi-axis triangulation with DepMap demoted) applied to IC/BPS with a UC rediscovery benchmark — with **honest negative-result capability** (RWR may win; conformal sets may be uselessly large; inter-study concordance may be low).

**We do not claim:** priority over TxGNN, CARNIVAL, or conformal-prediction methodology; nor SOTA methods performance given TARGET *n*.

**Demotions (success):** C1, C2-as-foundation, C4-as-first-OmniPath-drug-score, and C5-as-technique → **method/rigor components**. Residual publishable value = controlled evaluation + application + leakage/label design.

### Phase 4.5/4.6 control finding (binding for C1/C4)

| Claim component | Phase 4.5 result (pooled 5-fold OOF; repaired sign-shuffle) |
|---|---|
| **C1 (signed/directed graph as information)** | For valid signed models: **real ≈ sign_shuffled** but **real > rewired** (gatv2 0.590 / 0.588 / 0.516; ssl_gatv2 0.579 / 0.589 / 0.503). **Topology is used; sign is not detectably used.** |
| **C4 (sign-consistent drug scoring)** | Phase 5 three-arm: signed ≈ unsigned_matched on the **same reduced** edge set (paired AUPRC Δ ≈ 0; binding decision = no measurable sign benefit). unsigned_full (connectivity diagnostic) does not rescue sign. Consistent with Phase 4 control. |
| **Power caveat** | W1 MDD ≈ 0.128; observed GATv2−RWR gap ≈ 0.066 < MDD. Negative / null sign findings are limited by power, not proof of absence. |

---

## 3. Prior LUAD weaknesses → design responses (updated)

| ID | Weakness | Response after 0.5 |
|---|---|---|
| W1 | RWR beat GNN | Fair comparison retained; honesty rule retained; RWR itself at chance on corrected run |
| W2 | Label scarcity | SSL transfer retained as component; TARGET primary = **unsupervised** (*n*=38 filtered labels) |
| W3 | Degree confounding | Degree-matched nulls retained |
| W4 | mRNA ≠ activity | DoRothEA/PROGENy **core** (not optional) |
| W5 | Overclaimed CIs | Mondrian CP + mean set size; patient bootstrap separate |
| W6 | Single cohort / weak validation | Per-study DE + meta-analysis; GSE238208 specificity; GSE310557 cell type; L1000 replication; tissue expression axis |

---

## 4. Scientific question (verbatim)

**Primary objective:** Prioritize mechanism-consistent repurposing *hypotheses* for an understudied disease without disease-modifying therapy, by (1) learning or propagating disease-relevant genes on a signed directed multi-relational network with activity-aware features (supervised on UC; unsupervised-primary on IC/BPS), (2) scoring drugs by signed signaling opposition to the disease signature/activity footprint, and (3) retaining candidates with multi-axis support and Mondrian conformal sets — without claiming clinical efficacy.

Task definitions (nodes/edges/features/labels/outputs/falsifiers) remain as in Phase 0 README, with labels = **filtered Open Targets** only ([`LABEL_PROVENANCE.md`](LABEL_PROVENANCE.md)) and DepMap removed from primary triangulation.

---

## 5. Scope configuration

| Variable | Default |
|---|---|
| BENCHMARK | `MONDO_0005101` (UC) — confirmatory H1 |
| TARGET | `MONDO_0018301` (IC/BPS) — H7 **unsupervised-only**; **no curated labels** after adjudication (6&lt;N=15) |

**Reported evidence-base finding:** *The absence of usable curated IC/BPS disease-gene labels under non-drug, non-literature Open Targets evidence is itself a result.* With no supervised TARGET anchor, GSE238208 inflammation-specificity and the permuted-signature null are **mandatory** defences against topology-as-biology.
| Primary bulk | GSE11783 + GSE57560 |
| Planned controls | GSE238208 (inflammation specificity); GSE310557 (cell type) |

---

## 6. Venue suggestions (re-scoped — Phase 5.5 / P5.5.C)

**Reframing (results, not preference):** With zero FDR survivors, zero enrichment at the primary k values (@20/@50) on the disease-agnostic BENCHMARK, and TARGET failing H7 futility (mean top-20 Jaccard ≈ 0.051 ≪ 0.20), this is **not** an IC/BPS repurposing application paper with candidates. The honest contribution is an empirical demonstration, under pre-registration and leakage control, that **curated signed causal networks are too sparse and too incompletely signed to support signed-flow drug prioritization**, together with a negative-result-capable benchmark design that made that conclusion reportable rather than hidden.

| Finding | Status |
|---|---|
| No usable curated IC/BPS disease-gene labels under non-drug/non-literature evidence | Reported negative |
| TARGET cohorts do not replicate at the extremes; meta is single-cohort dominated | Reported negative |
| GNN shows no advantage over RWR (MDD 0.128); RWR itself at chance | Reported negative |
| Features + `log_degree` is the strongest model | Positive but unflattering to the design |
| Edge sign contributes nothing detectable (Phase 4 control + Phase 5 signed vs unsigned_matched) | Reported negative |
| Conformal sets cover but are uninformative (mean size ~1.79 of 2) | Reported negative |
| Drug scoring: weak global signal (~3.4× chance AUPRC), no top-k@20/50, no FDR survivors | Mixed → practically null |

**Recurring explanation:** substrate sparsity — 7 773 nodes, median degree 3, 274 components (1 313 after signed exclusion), receptor-class coverage dominated by residual “other signalling enzyme.”

- **Preferred now:** bioinformatics **methods** journals, or dedicated **negative-results / registered-report** tracks emphasizing pre-registered null capability and signed-network sparsity limits.  
- **Demoted:** translational/urology application venues (*Neurourology and Urodynamics*, *BMC Urology*) that fit a paper with surviving candidates — there are none.  
- **Not appropriate:** pitching as TxGNN competitor or *Nat Med*-class methods priority.

Logged in `ANALYSIS_PLAN_AMENDMENTS.md` (P5.5.C). Phase 6 mode = descriptive characterisation of the frozen seed module only (`docs/phase0/PHASE6_RESCOPE_DECISION.json`); do not manufacture survivors.

---

## 7. Non-claims / ethics

- No treatment recommendations; no efficacy claims from scores.  
- Conformal coverage ≠ patient generalization.  
- Unsupervised–supervised agreement ≠ orthogonal validation (shared signature).  
- Ethics: public, de-identified, previously published data; no IRB; no patient-level raw data in git.
