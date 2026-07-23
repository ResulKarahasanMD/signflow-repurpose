# Phase 0 literature review (verified identifiers only)

**Access date (UTC):** 2026-07-23  
**Scope:** Confirm IC/BPS therapeutic gap; map prior computational repurposing/bioinformatics; support novelty differentiation.  
**Method:** NCBI E-utilities (PubMed), ClinicalTrials.gov API v2, Open Targets GraphQL, Exa/web for discovery then PMID verification.

Unverified DOIs/PMIDs are omitted.

---

## 1. Therapy landscape

| Citation | Key finding used here |
|---|---|
| PMID [38658949](https://pubmed.ncbi.nlm.nih.gov/38658949/) — BMC Urol 2024; DOI `10.1186/s12894-024-01485-w` | NMA of 70 RCTs; no treatment improved **all** outcomes in NMA |
| PMID [33638677](https://pubmed.ncbi.nlm.nih.gov/33638677/) — Int Urogynecol J 2021; DOI `10.1007/s00192-020-04659-w` | Pharmacotherapy NMA; cyclosporine A / amitriptyline signals without disease-modifying standard |
| PMID [21497847](https://pubmed.ncbi.nlm.nih.gov/21497847/) — AUA guideline 2011 | Foundational IC/BPS diagnosis/treatment guideline |
| PMID [25623737](https://pubmed.ncbi.nlm.nih.gov/25623737/) — AUA guideline amendment 2015 | Updated multimodal symptomatic framework |
| Open Targets `MONDO_0018301` drug candidates (API 2026-07-23) | 49 candidates; sole APPROVAL-stage: pentosan polysulfate sodium (`CHEMBL1201516`) |

**ClinicalTrials.gov:** condition query `Interstitial Cystitis` → `totalCount=238` (2026-07-23). Example interventional drug studies: NCT05744908 (pentosan), NCT01990898 (cyclosporine), NCT00124306 (amitriptyline), NCT02787083 (mirabegron).

---

## 2. Primary expression biology papers tied to our cohorts

| PMID | Cohort / note |
|---|---|
| [19400928](https://pubmed.ncbi.nlm.nih.gov/19400928/) | GSE11783 ulcerative IC bladder tissue microarray |
| [24840534](https://pubmed.ncbi.nlm.nih.gov/24840534/) | GSE57560 capacity-stratified IC/BPS bladder biopsies |
| [38026177](https://pubmed.ncbi.nlm.nih.gov/38026177/) | GSE238208 HIC vs BCG (secondary only) |
| [22177197](https://pubmed.ncbi.nlm.nih.gov/22177197/) | GSE28242 urine sediment (orthogonal biospecimen) |

---

## 3. Prior computational / ML / repurposing for IC/BPS

| PMID | Approach | Differentiation vs this project |
|---|---|---|
| [39723619](https://pubmed.ncbi.nlm.nih.gov/39723619/) | cMAP / L1000CDS2 signature reversal | No signed multi-relational GNN, SSL transfer, conformal, DepMap/genetics triangulation |
| [36619718](https://pubmed.ncbi.nlm.nih.gov/36619718/) | LASSO/SVM-RFE/RF + CMap + docking | Feature ML; not relational signed message passing |
| [40777016](https://pubmed.ncbi.nlm.nih.gov/40777016/) | Multi-algorithm ML + docking + mouse | Diagnostic/therapeutic gene focus; not causal-flow GNN stack |
| [35333927](https://pubmed.ncbi.nlm.nih.gov/35333927/) | WGCNA + docking | Coexpression, not signed GNN |
| [35401223](https://pubmed.ncbi.nlm.nih.gov/35401223/) | Adhesion-molecule target framing | Narrative/target biology; not our methods stack |
| [39917301](https://pubmed.ncbi.nlm.nih.gov/39917301/), [38338847](https://pubmed.ncbi.nlm.nih.gov/38338847/) | Biomarker / transcriptome ML | Diagnosis-oriented |

PubMed search used (runtime):

```text
("interstitial cystitis"[Title] OR "bladder pain syndrome"[Title])
AND (repurpos*[Title/Abstract] OR "connectivity map"[Title/Abstract]
     OR bioinformatic*[Title/Abstract] OR "machine learning"[Title/Abstract])
```

**Result:** `count=50` (2026-07-23); top hits reviewed via esummary; table above retains those relevant to computational drug/target prioritization.

---

## 4. Resource endpoints verified for later phases

| Resource | Verification |
|---|---|
| OmniPath interactions API | Directed + stimulation/inhibition fields present |
| STRING | version `12.0` |
| SIGNOR | download portal live |
| DGIdb | GraphQL at `https://dgidb.org/api/graphql` live |
| Open Targets | GraphQL API live; disease IDs above |

---

## 5. Literature gaps remaining (not claimed as searched exhaustively)

- No PubMed/Consensus MCP server was available in this environment; search used NCBI E-utilities + Exa/web.  
- Phase 7 will re-run the PubMed query and bioRxiv/medRxiv title search before manuscript submission.  
- CP/CPPS as alternate TARGET was not data-gated in Phase 0; requires its own ≥2-cohort gate before switching.
