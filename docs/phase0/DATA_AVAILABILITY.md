# Phase 0 — Data-availability gate

**Access date (UTC):** 2026-07-23  
**Verdict:** **PASS (conditional)** — ≥2 independent open human expression cohorts exist for IC/BPS bladder tissue with case–control labels. Proceed with TARGET = IC/BPS, with **honest power and phenotype caveats** below. Benchmark track (ulcerative colitis) is independently well-powered.

This document records **runtime-verified** queries only. No accessions were taken from memory without NCBI/Open Targets confirmation.

---

## 1. Gate criteria (pre-specified)

| Criterion | Threshold | Result |
|---|---|---|
| Independent open expression datasets for TARGET disease | ≥2 | **Met** |
| Human, disease-relevant tissue (not adjacent disease) | bladder ± documented alternate biospecimen | **Met** (bladder biopsies) |
| Case vs non-disease control contrast available | yes | **Met** (with small control *n*) |
| Do not substitute adjacent disease (e.g. BCG cystitis ≠ healthy) | yes | **Enforced** (see exclusions) |
| Benchmark disease open data + known drugs | yes | **Met** (UC) |

---

## 2. Queries executed (reproduce)

### 2.1 GEO / NCBI Gene Expression Omnibus (GDS)

```text
db=gds
term=(interstitial cystitis[Title]) AND (Expression profiling by array[DataSet Type]
     OR Expression profiling by high throughput sequencing[DataSet Type])
```

**Result (2026-07-23):** `count=8` series UIDs  
`200328464, 200310557, 200238208, 200057560, 200028242, 200011783, 200011839, 200000621`

Mapped accessions (esummary):  
`GSE328464, GSE310557, GSE238208, GSE57560, GSE28242, GSE11783, GSE11839, GSE621`

### 2.2 ClinicalTrials.gov API v2

```text
GET https://clinicaltrials.gov/api/v2/studies?query.cond=Interstitial%20Cystitis&pageSize=10&countTotal=true
```

**Result:** `totalCount=238` studies matching the condition query (mixed observational/interventional; relevance varies).

Drug-focused follow-up (examples returned): NCT05744908 (pentosan), NCT01990898 (cyclosporine), NCT00124306 (amitriptyline), NCT02787083 (mirabegron).

### 2.3 Open Targets Platform GraphQL

```graphql
{ search(queryString: "interstitial cystitis", entityNames: ["disease"], page: {index: 0, size: 10}) {
    hits { id name } } }
```

**Primary disease ID used:** `MONDO_0018301` (interstitial cystitis).  
Related: `EFO_1000869` (chronic interstitial cystitis).

```graphql
{ disease(efoId: "MONDO_0018301") {
    associatedTargets(page: {index: 0, size: 15}) { count rows { score target { approvedSymbol } } }
    drugAndClinicalCandidates { count rows { maxClinicalStage drug { id name } } } } }
```

**Result:** 514 associated targets; 49 drug/clinical candidates.  
Only **APPROVAL**-stage agent in the returned list for this ID: **pentosan polysulfate sodium** (`CHEMBL1201516`). Remaining entries are Phase 1–3 / unknown (symptomatic or investigational).

Benchmark disease: `MONDO_0005101` (ulcerative colitis) — 7161 associated targets; 245 drug/clinical candidates including multiple APPROVAL-stage agents (e.g. adalimumab, golimumab, budesonide, olsalazine).

### 2.4 Network / drug-target resource liveness

| Resource | Check | Result |
|---|---|---|
| OmniPath | `GET https://omnipathdb.org/interactions?...&format=json` | HTTP 200; fields `is_directed`, `is_stimulation`, `is_inhibition` present |
| STRING | `GET https://string-db.org/api/json/version` | `string_version: "12.0"` |
| SIGNOR | `https://signor.uniroma2.it/downloads.php` | Site live (download portal) |
| DGIdb | `POST https://dgidb.org/api/graphql` | HTTP 200 (GraphQL); legacy `/api/v2/` returns SPA HTML — **use GraphQL** |

---

## 3. TARGET disease cohorts — inclusion decision

### 3.1 Primary analysis cohorts (bulk bladder tissue, case vs control)

| Accession | Platform | *n* (verified) | Design (from GEO SOFT / series matrix) | Role |
|---|---|---|---|
| **GSE11783** | GPL570 (Affymetrix U133 Plus 2.0) | 16 arrays | Ulcerative IC bladder biopsies (ulcer + non-ulcer sites from same patients) vs healthy controls; PMID [19400928](https://pubmed.ncbi.nlm.nih.gov/19400928/) | **Primary cohort A** |
| **GSE57560** | GPL16699 (Agilent) | 16 arrays | IC/BPS bladder biopsy: low capacity / normal capacity / control; disease-state counts: IC=13, control=3; PMID [24840534](https://pubmed.ncbi.nlm.nih.gov/24840534/) | **Primary cohort B** |

**GSE11783 phenotype counts (series matrix, runtime):**  
`disease state: no` = 6 (controls); `IC, non-ulcer tissue` = 5; `IC, ulcer` = 5.  
Note: ulcer and non-ulcer samples from the **same patients** are **not** independent biological replicates — analysis plan must account for patient-level correlation.

**GSE57560 phenotype counts (series matrix, runtime):**  
`Interstitial cystitis` = 13; `control` = 3. Capacity phenotype: low=4, normal=9, control=3.

### 3.2 Planned non-primary cohorts (Phase 0.5 roles; still not counted toward the ≥2 healthy case–control gate)

| Accession | Planned role (pre-specified) | Restrictions |
|---|---|---|
| **GSE238208** | **Inflammation-specificity control:** does the prioritized module/signature separate Hunner-type IC from BCG cystitis (generic bladder inflammation)? *n*=63. PMID [38026177](https://pubmed.ncbi.nlm.nih.gov/38026177/). | **Forbidden** as IC-vs-healthy DE |
| **GSE310557** | **Cell-type attribution.** **Re-verified 2026-07-23:** *n*=17 (10 HIC + 7 PSUI); GPL34284; scRNA-seq; PMID [41457077](https://pubmed.ncbi.nlm.nih.gov/41457077/) DOI `10.1038/s41467-025-68060-1`. | **Not** merged into bulk DE |
| **GSE28242** | Optional urine orthogonal signature only. PMID [22177197](https://pubmed.ncbi.nlm.nih.gov/22177197/). | Different biospecimen |
| **GSE11839** | Not planned for primary analyses | Cultured cells |
| **GSE621** | Not planned for primary analyses | Historical explants |
| **GSE328464** | Sensitivity only after QC (*n*=6) | Underpowered |

### 3.3 Honest power / phenotype caveats (required in Methods/Limitations)

1. **Small *n***: both primary cohorts have 16 arrays; GSE57560 has only **3** controls → wide uncertainty; conformal sets may be large (report mean set size).
2. **Platform heterogeneity**: Affymetrix vs Agilent → **primary path is per-study DE + random-effects meta-analysis** (Phase 0.5 A4); pooled ComBat is sensitivity-only, not primary.
3. **Phenotype mix**: Hunner/ulcerative (GSE11783) vs capacity-stratified IC/BPS (GSE57560) — report concordance; low concordance is a real finding.
4. **Paired biopsies** in GSE11783 inflate apparent sample size if treated as independent (one sample per patient primary).
5. **No fabricated cohorts**: we do **not** treat BCG cystitis, overactive bladder, or bacterial cystitis as IC/BPS healthy controls or as substitutes for the target disease.

**Gate decision:** PASS for TARGET = IC/BPS using **GSE11783 + GSE57560** as the minimum dual-cohort set. Label all IC/BPS application results as **data-limited / hypothesis-generating**.

---

## 4. BENCHMARK disease (ulcerative colitis) — positive control track

**Disease ID:** `MONDO_0005101`.

**GEO search (runtime):** ulcerative colitis + human + colon + expression → **count=106** series.

**Verified example series with usable *n* (esummary):**

| Accession | *n* | Type |
|---|---|---|
| GSE75214 | 194 | array (IBD mucosal) |
| GSE87466 | 108 | array (UC mucosal) |
| GSE59071 | 116 | array (IBD mucosal) |
| GSE107499 | 119 | array (UC) |
| GSE128682 | 44 | RNA-seq (UC remission signatures) |

Final UC cohort pair will be locked in Phase 2 after downloading sample phenotype tables and confirming case/control or inflamed/uninflamed contrasts (no memory-hardcoding of group labels).

**Gold-standard rediscovery set (Open Targets APPROVAL-stage examples for UC, runtime):** adalimumab, golimumab, budesonide, methylprednisolone, prednisolone, olsalazine sodium, among others. Exact rediscovery list will be frozen in `config/benchmark_drugs.yaml` in Phase 1 from Open Targets + curated guideline agents — not from DE.

---

## 5. What would have failed the gate

- Fewer than 2 independent **IC/BPS** expression series with open download → STOP TARGET application; run BENCHMARK-only methods paper.
- Only urine / cell-line / BCG-contrast data available → STOP TARGET bulk DE; report missing bladder case–control data.
- Inventing a “proxy cohort” (e.g. general cystitis, BPS animal model only) → forbidden.

**None of these failure modes apply to the dual bladder-biopsy set above.**

---

## 6. Next phase gate

Phase 1 may proceed. Phase 2 must re-verify each accession’s supplementary files are downloadable and record SHA-256 in `run_manifest.json` before DE.
