# IC/BPS `animal_model` label adjudication (Phase 0.6 A6.1)

**Access date (UTC):** 2026-07-23  
**Outcome-dependent?** No (pre-results; no model output consulted).  
**Disease:** `MONDO_0018301`  
**Channel:** Open Targets datasource `impc` / datatype `animal_model`  
**Input set:** 38 genes with `animal_model` datatype score ≥ 0.30 (from [`LABEL_PROVENANCE.md`](LABEL_PROVENANCE.md))

---

## 0. Pre-committed decision rule (written BEFORE adjudication)

| Parameter | Frozen value |
|---|---|
| **N (minimum `plausible-ICBPS`)** | **15** |
| Rule | If fewer than **N=15** genes are classified `plausible-ICBPS`, the `animal_model` channel is **dropped from TARGET labels entirely** |
| Reported finding if dropped | *IC/BPS has no usable curated disease-gene label set under non-drug, non-literature evidence.* |
| Hard constraint | These genes (whether kept or dropped as labels) **must never seed** the unsupervised TARGET path |

This section was authored before the classification table below was filled.

---

## 1. Adjudication rubric (written BEFORE classification)

Apply using **IMPC mouse phenotype terms returned by Open Targets**, not gene-name intuition alone.

| Class | Definition |
|---|---|
| **`plausible-ICBPS`** | Evidence includes **urinary bladder inflammation** and/or **bladder-related pain/nociception** phenotypes (or closely equivalent bladder inflammatory/fibrotic urothelial injury) **without** the dominant phenotype being Mendelian stone disease, primary nephropathy, congenital urinary-tract malformation, or multi-organ lethality where bladder is incidental. |
| **`urinary-tract-nonspecific`** | Bladder morphology/physiology, incontinence, frequency, detrusor contractility, or urothelial structure **without** IC/BPS-like inflammatory/pain signature; or bladder phenotype is secondary to obstruction/reflux/developmental defects. |
| **`renal-or-stone-disease`** | Dominant phenotypes are nephrolithiasis, cystinuria, crystalluria, glomerulopathy, tubulointerstitial nephropathy, polycystic kidney, renal failure, or analogous stone/renal Mendelian disease. |
| **`implausible`** | Non-protein-coding human biotype; wrong organ system dominance; reproductive-only; or no bladder-relevant phenotype. |

Biotype rule: any human gene that is **not** `protein_coding` is **`implausible`** for the label set **and** must be excluded from the node universe (mapping QC regression).

---

## 2. Evidence pull (runtime)

**Query pattern (verified schema requires `ensemblIds`):**

```graphql
{
  disease(efoId: "MONDO_0018301") {
    evidences(ensemblIds: ["<ENSG...>"], datasourceIds: ["impc"]) {
      count
      rows {
        score resourceScore
        biologicalModelAllelicComposition biologicalModelId biologicalModelGeneticBackground
        diseaseModelAssociatedModelPhenotypes { id label }
        target { approvedSymbol id }
      }
    }
  }
}
```

All 38 genes returned ≥1 IMPC evidence row (2026-07-23).

**Biotype sources:** Open Targets `target.biotype` for `UOX`; mygene.info `type_of_gene` for all 38 (Ensembl REST rate-limited/500 during pull).

| Symbol | Human biotype (runtime) | Action |
|---|---|---|
| **UOX** | `transcribed_unitary_pseudogene` (OT `ENSG00000240520`; mygene: urate oxidase pseudogene) | **Remove from labels and node universe** |
| All others in the 38 | `protein-coding` (mygene) | Eligible for rubric classification |

---

## 3. Per-gene adjudication table

| Gene | am_score | Key IMPC phenotypes (runtime) | Class |
|---|---|---|---|
| SLC3A1 | 0.503 | cystolithiasis, nephrolithiasis, cystinuria, detrusor changes | renal-or-stone-disease |
| HTR3A | 0.480 | urinary incontinence, pyelonephritis, prostate inflammation | urinary-tract-nonspecific |
| UPK3A | 0.442 | abnormal urothelium, increased urination frequency, VUR | urinary-tract-nonspecific |
| KCNMA1 | 0.415 | increased urination frequency, urinary incontinence | urinary-tract-nonspecific |
| APRT | 0.408 | hydronephrosis, renal tubule dilation, colorless urine | renal-or-stone-disease |
| CALD1 | 0.406 | abnormal detrusor morphology/physiology | urinary-tract-nonspecific |
| GATA5 | 0.398 | female reproductive morphology; vaginal inflammation | implausible |
| UMOD | 0.385 | increased creatinine, abnormal urine chemistries | renal-or-stone-disease |
| CHRM3 | 0.373 | hydronephrosis, distended bladder, impaired detrusor contractility | urinary-tract-nonspecific |
| HOXA13 | 0.373 | absent urinary bladder; abnormal renal/urinary morphology | urinary-tract-nonspecific |
| AHR | 0.369 | urinary bladder inflammation, fibrosis, abnormal urothelium | **plausible-ICBPS** |
| SCARA5 | 0.366 | bladder inflammation amid multi-organ inflammation | urinary-tract-nonspecific |
| UOX | 0.357 | crystalluria, renal fibrosis (human pseudogene) | **implausible** |
| CHRNB4 | 0.354 | abnormal bladder physiology/morphology | urinary-tract-nonspecific |
| PAPOLB | 0.350 | abnormal bladder morphology; infertility/uterus | urinary-tract-nonspecific |
| MEP1A | 0.341 | urinary bladder inflammation; kidney brush-border loss | **plausible-ICBPS** |
| TAC1 | 0.339 | abnormal pain threshold; bladder physiology; incontinence | **plausible-ICBPS** |
| HPSE2 | 0.338 | incontinence, frequency, bladder fibrosis (urofacial spectrum) | urinary-tract-nonspecific |
| CHRNA3 | 0.335 | bladder inflammation + incontinence; kidney cortex cyst | urinary-tract-nonspecific |
| PAX2 | 0.335 | VUR, abnormal ureter/kidney developmental phenotypes | urinary-tract-nonspecific |
| DEFB1 | 0.334 | abnormal urinary bladder physiology only | urinary-tract-nonspecific |
| TRPV4 | 0.331 | abnormal pain threshold; nociception after inflammation; bladder | **plausible-ICBPS** |
| MAZ | 0.329 | bladder hypoplasia; abnormal bladder urine volume | urinary-tract-nonspecific |
| UPK2 | 0.329 | abnormal urine electrolytes/creatinine | urinary-tract-nonspecific |
| GSTK1 | 0.326 | bladder inflammation; kidney weight/glomerulus changes | urinary-tract-nonspecific |
| CD209 | 0.326 | urinary bladder inflammation; macrophage/infection phenotypes | **plausible-ICBPS** |
| CLEC4M | 0.326 | urinary bladder inflammation; macrophage/infection phenotypes | **plausible-ICBPS** |
| MYOCD | 0.322 | distended bladder; abnormal detrusor morphology | urinary-tract-nonspecific |
| PI4K2A | 0.319 | urinary incontinence; systemic phenotypes | urinary-tract-nonspecific |
| UPK1A | 0.319 | bladder fibrosis; abnormal urothelium; multi-organ | urinary-tract-nonspecific |
| KRT7 | 0.316 | abnormal urinary bladder physiology | urinary-tract-nonspecific |
| GEN1 | 0.315 | duplex kidney, VUR, ureteric bud defects | urinary-tract-nonspecific |
| ARHGDIA | 0.313 | renal cast, kidney failure | renal-or-stone-disease |
| ACTG2 | 0.311 | distended bladder; oliguria; gut smooth-muscle phenotypes | urinary-tract-nonspecific |
| DST | 0.309 | urinary incontinence; neuropathy phenotypes | urinary-tract-nonspecific |
| COL4A3 | 0.308 | glomerulosclerosis, renal tubule atrophy, albuminuria | renal-or-stone-disease |
| MYH11 | 0.307 | distended bladder; vascular/cardiac phenotypes | urinary-tract-nonspecific |
| ARIH2 | 0.301 | polycystic kidney; cystolithiasis; abnormal bladder morphology | renal-or-stone-disease |

---

## 4. Counts and decision

| Class | *n* |
|---|---|
| plausible-ICBPS | **6** (AHR, MEP1A, TAC1, TRPV4, CD209, CLEC4M) |
| urinary-tract-nonspecific | 24 |
| renal-or-stone-disease | 6 |
| implausible | 2 (GATA5, UOX) |

**N = 15; plausible-ICBPS = 6 < 15 → DROP `animal_model` from TARGET labels entirely.**

### Frozen consequence

1. TARGET has **no usable curated disease-gene label set** under non-drug, non-literature evidence (publishable negative finding).  
2. TARGET **supervised / few-shot GNN path is cancelled** (not merely demoted).  
3. TARGET **primary path remains unsupervised** (signature-extreme seeds only).  
4. **UOX** excluded from node universe; mapping QC must reject non-`protein_coding` biotypes.  
5. The 38-gene list is **never** used as unsupervised seeds.
