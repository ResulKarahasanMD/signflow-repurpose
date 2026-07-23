# Methods prior art and honest novelty re-scoping

**Access date (UTC):** 2026-07-23  
**Purpose:** Phase 0.5 amendment A1 — methods literature review for claims C1–C6. Disease-specific IC/BPS *in silico* work remains in [`LITERATURE_REVIEW.md`](LITERATURE_REVIEW.md); this file covers **methodological** prior art.

All PMIDs/DOIs below were verified via NCBI E-utilities (or resource APIs) on 2026-07-23 unless noted.

---

## 1. Mandatory priors (verified)

| Work | Verified citation | Overlap with this project |
|---|---|---|
| **TxGNN** | Huang et al., *Nat Med* 2024; PMID [39322717](https://pubmed.ncbi.nlm.nih.gov/39322717/); DOI `10.1038/s41591-024-03233-x` | Foundation GNN + metric-learning transfer for **zero-shot repurposing in diseases with limited/no treatments** across ~17k diseases. Directly overlaps **C2, C3**, and the “understudied / no DMT” framing. |
| **CARNIVAL** | Liu et al., *npj Syst Biol Appl* 2019; PMID [31728204](https://pubmed.ncbi.nlm.nih.gov/31728204/); DOI `10.1038/s41540-019-0118-z` | Causal reasoning on a **signed, directed OmniPath** prior; uses perturbation/drug targets with **signed effects**. Directly overlaps **C1, C4**. |
| **Conformal prediction in drug discovery** | Alvarsson et al., *J Pharm Sci* 2021; PMID [33075380](https://pubmed.ncbi.nlm.nih.gov/33075380/); DOI `10.1016/j.xphs.2020.09.055` | CP as established uncertainty tool in drug discovery. **C5 is technique adoption, not a novel contribution.** |
| **Mondrian / class-conditional CP (imbalanced)** | Sun et al., *J Chem Inf Model* 2017; PMID [28628322](https://pubmed.ncbi.nlm.nih.gov/28628322/); DOI `10.1021/acs.jcim.7b00159`. Also Norinder & Boyer lineage: PMID [24797111](https://pubmed.ncbi.nlm.nih.gov/24797111/); DOI `10.1021/ci5001168` | Class-conditional calibration for imbalanced bioactivity — required for our PU/popularity-biased labels (amendment A3). |

---

## 2. Additional methods searches (runtime; recorded)

| Query theme | Example verified hits |
|---|---|
| Footprint activity (DoRothEA) | Garcia-Alonso et al., *Genome Res* 2019; PMID [31340985](https://pubmed.ncbi.nlm.nih.gov/31340985/); DOI `10.1101/gr.240663.118` |
| Pathway footprints (PROGENy) | Schubert et al., *Nat Commun* 2018; PMID [29295995](https://pubmed.ncbi.nlm.nih.gov/29295995/); DOI `10.1038/s41467-017-02391-6` |
| TF activity ↔ drug sensitivity | Garcia-Alonso et al., *Cancer Res* 2018; PMID [29229604](https://pubmed.ncbi.nlm.nih.gov/29229604/) |
| decoupleR | Badia-i-Mompel et al., *Bioinform Adv* 2022; PMID [36699385](https://pubmed.ncbi.nlm.nih.gov/36699385/); DOI `10.1093/bioadv/vbac016` |
| Network proximity repurposing | Guney et al., *Nat Commun* 2016; PMID [26831545](https://pubmed.ncbi.nlm.nih.gov/26831545/); Cheng et al., *Nat Commun* 2018; PMID [30002366](https://pubmed.ncbi.nlm.nih.gov/30002366/) |
| PrimeKG | Chandak et al., *Sci Data* 2023; PMID [36732524](https://pubmed.ncbi.nlm.nih.gov/36732524/); DOI `10.1038/s41597-023-01960-3` |
| GNN disease–gene prioritization | e.g. XGDAG, PMID [37531293](https://pubmed.ncbi.nlm.nih.gov/37531293/); additional GCN prioritization methods exist (PMID [36241966](https://pubmed.ncbi.nlm.nih.gov/36241966/)) |
| IC/BPS L1000 reversal (application prior) | PMID [39723619](https://pubmed.ncbi.nlm.nih.gov/39723619/) — frames our L1000 axis as **replication**, not novelty |

**PubMed queries used (abbreviated):** `TxGNN[Title]`; `CARNIVAL[Title] AND Saez-Rodriguez`; conformal drug-discovery title search; Mondrian conformal + Norinder/Sun; DoRothEA / PROGENy footprint titles; `Building a knowledge graph to enable precision medicine[Title]`; network proximity Guney/Cheng titles; `graph neural network disease gene prioritization[Title/Abstract]`.

**Note:** Classic R-GCN (Schlichtkrull et al.) and SGCN (Derr et al.) are foundational ML citations often indexed outside PubMed; we treat them as established architectural baselines, not novelty claims.

---

## 3. Per-contribution differentiation table

| Contribution | Closest prior art (verified) | What is genuinely different here | Residual novelty (honest) |
|---|---|---|---|
| **C1 — Signed/directed multi-relational mechanistic graph** | CARNIVAL/OmniPath (PMID 31728204); SIGNOR; PrimeKG multi-relational KG (PMID 36732524) | We reuse OmniPath/SIGNOR as **input**, not invent the signed directed prior. Difference is packaging for **relational GNN + controlled ablations** (sign/direction/type shuffle) under a fixed analysis plan. | **Low as a standalone claim.** Demote C1 from “novel resource” to **method component**: *use of* an established signed directed prior with explicit graph-information controls. |
| **C2 — SSL pretraining + few-shot transfer** | TxGNN zero-shot / transfer across diseases with limited treatments (PMID 39322717); broader SSL on bio graphs | TxGNN already targets scarce-indication transfer at KG scale. Our difference is transfer on a **signed causal molecular graph** with **expression/activity node features** and an honest **from-scratch vs pretrained** table on a small-n disease — not a new foundation model. | **Low–moderate.** Residual = *controlled transfer ablation under leakage-safe omics features*, not priority over TxGNN. Demote from “foundation-model novelty” to **method component + evaluation**. |
| **C3 — Relational/signed GNN that can beat propagation** | TxGNN; many biomedical GNN prioritizers (e.g. PMID 37531293); RWR/proximity baselines (Guney/Cheng) | Fair **head-to-head vs RWR** with degree/sign/edge controls and pre-registered honesty rule if RWR wins (W1). | **Moderate as evaluation science**, not as architecture invention. Architecture itself is standard; the **negative-result-capable benchmark** is the publishable piece. |
| **C4 — Signaling-flow / sign-consistent drug scoring** | CARNIVAL signed perturbation reasoning (PMID 31728204); network proximity (PMID 26831545, 30002366); footprint activity (DoRothEA/PROGENy) | CARNIVAL already does signed causal reasoning with drug targets on OmniPath. Residual difference: **score definition + degree-matched nulls + conformal wrapping + triangulation**, with footprint activities as **core features** (not optional). | **Low–moderate.** Demote from “first signed drug scoring on OmniPath” to **method component**: *propagation/flow scoring with degree-matched nulls and activity features*, explicitly differentiated from ILP/CARNIVAL. |
| **C5 — Conformal candidate sets** | Alvarsson et al. (PMID 33075380); Mondrian CP (PMID 28628322) | Application of **Mondrian CP** to gene/drug prioritization sets with **coverage + mean set size** reported; explicit non-claim on patient exchangeability. | **None as technique novelty.** Demote C5 to **rigor/reporting component**. |
| **C6 — Multi-evidence triangulation** | Standard in translational bioinformatics; IC/BPS already used L1000 (PMID 39723619) | Axis set revised (tissue expression / scRNA instead of DepMap as primary); L1000 framed as **replication**; genetics filtered to avoid drug-label circularity. | **Low as novelty; moderate as disciplined reporting.** Residual = transparent disagreement tables + pre-specified convergence test on UC gold standards. |

---

## 4. Rewritten novelty claim (authoritative)

**Permitted claim:** This project offers a **combination + application + rigorously controlled benchmark** for mechanism-aware drug *hypothesis generation* on an understudied pain/inflammatory indication (IC/BPS), with a well-powered inflammatory benchmark (UC). It does **not** claim priority over TxGNN for transfer/zero-shot repurposing, over CARNIVAL/DoRothEA/PROGENy for signed OmniPath causal reasoning or footprint activity, or over the conformal-prediction literature for uncertainty quantification.

**What remains worth publishing if results are honest:**
1. A **pre-registered, leakage-aware** pipeline (filtered Open Targets labels; per-study DE + meta-analysis; Mondrian CP with set-size reporting).  
2. A **fair GNN-vs-RWR (+ graph information controls)** evaluation that is allowed to conclude RWR wins.  
3. An **IC/BPS application** with inflammation-specificity (GSE238208) and cell-type attribution (GSE310557), plus L1000 agreement/disagreement vs PMID 39723619.  
4. Explicit **power limitations**: cohort sizes preclude a state-of-the-art methods claim.

**Demotions (success signal):** C1, C2 (as foundation novelty), C4 (as first signed OmniPath drug score), and C5 are **method components**, not priority claims. C3/C6 retain residual novelty only as **controlled evaluation / triangulation design**.

---

## 5. Venue implications

Given residual novelty and *n*:

- Prefer **application + transparent methods** venues (*Neurourol Urodyn*, *BMC Urology*, computational biology application tracks) or **methods-honesty / negative-result** framing in bioinformatics journals.  
- Do **not** pitch as competing with TxGNN (*Nat Med*-class foundation-model priority).  
- A methods journal is appropriate only if the **controlled benchmark + leakage/conformal design** is the lead story — not “a new GNN.”
