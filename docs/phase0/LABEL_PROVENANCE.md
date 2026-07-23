# Label provenance — Open Targets evidence filters (frozen)

**Access date (UTC):** 2026-07-23  
**Purpose:** Phase 0.5 amendment A2 — eliminate drug-channel and popularity (literature) circularity in disease–gene labels.  
**Outcome-dependent?** No (no model results exist).

---

## 1. Problem

Unfiltered Open Targets `associatedTargets.score` aggregates evidence channels. Two circularities are fatal for this project:

1. **Drug/clinical circularity:** `clinical` / `clinical_precedence` elevates genes because drugs already target them — then “rediscovering” those drugs is grading the model on leaked information.  
2. **Popularity circularity:** `literature` / `europepmc` elevates well-studied genes — exactly the bias an understudied-disease project claims to overcome.

Additionally, **`rna_expression` / `expression_atlas`** must not define labels when expression DE is used as node features / signature (feature–label circularity).

---

## 2. Runtime GraphQL (exact)

Endpoint: `https://api.platform.opentargets.org/api/v4/graphql`

```graphql
{
  disease(efoId: "<DISEASE_ID>") {
    id
    name
    associatedTargets(page: {index: 0, size: 500}) {
      count
      rows {
        score
        datatypeScores { id score }
        datasourceScores { id score }
        target { id approvedSymbol }
      }
    }
  }
}
```

Pagination: increment `page.index` until all `count` rows retrieved.  
Diseases queried: `MONDO_0018301` (IC/BPS), `MONDO_0005101` (UC).

Machine-readable dump of this pull is **not** committed (large); counts below are frozen from the 2026-07-23 query recorded in `provenance/phase0_query_log.json` (amended) and this file.

---

## 3. Datatype / datasource breakdown (presence counts)

### 3.1 IC/BPS — `MONDO_0018301` (*n* associations = 514)

| datatypeId | Targets with datatype | Primary allowlist? |
|---|---|---|
| `literature` | 341 | **EXCLUDE** |
| `animal_model` | 139 | **ALLOW** (only surviving non-drug/non-literature channel) |
| `clinical` | 61 | **EXCLUDE** |

| datasourceId | Targets with datasource |
|---|---|
| `europepmc` | 341 |
| `impc` | 139 |
| `clinical_precedence` | 61 |

**No `genetic_association` channel is present for this disease ID in the 2026-07-23 pull.**

### 3.2 UC — `MONDO_0005101` (*n* associations = 7161)

| datatypeId | Targets with datatype | Primary allowlist? |
|---|---|---|
| `rna_expression` | 4914 | **EXCLUDE** (feature circularity with DE) |
| `literature` | 2588 | **EXCLUDE** |
| `genetic_association` | 1050 | **ALLOW** |
| `clinical` | 239 | **EXCLUDE** |
| `somatic_mutation` | 71 | **ALLOW (sensitivity union only)** |
| `genetic_literature` | 3 | **EXCLUDE** (literature-adjacent) |

| datasourceId | Targets with datasource |
|---|---|
| `expression_atlas` | 4914 |
| `europepmc` | 2588 |
| `gwas_credible_sets` | 1049 |
| `clinical_precedence` | 239 |
| `cancer_gene_census` | 71 |
| `genomics_england` | 3 |
| `gene_burden` | 1 |

---

## 4. Frozen primary label definition

### 4.1 Score used for thresholding

For each target–disease row, compute:

```text
label_score = max({ datatypeScores[d] | d ∈ ALLOWLIST and d present })
```

If no allowlisted datatype is present → target is **unlabeled** (not a negative).

### 4.2 Allowlists

| Track | Disease ID | Primary ALLOWLIST | Explicit EXCLUDE |
|---|---|---|---|
| BENCHMARK | `MONDO_0005101` | `{genetic_association}` | `clinical`, `literature`, `rna_expression`, `genetic_literature` |
| TARGET | `MONDO_0018301` | `{animal_model}` | `clinical`, `literature` |

**Sensitivity (pre-specified, not primary):**
- UC: ALLOWLIST `{genetic_association, somatic_mutation}` at same thresholds.  
- Both: thresholds in §4.3.  
- **Never** re-introduce `clinical` or `literature` into labels used when drug rediscovery is an endpoint.

### 4.3 Thresholds (pinned now)

| Role | Threshold on `label_score` |
|---|---|
| **Primary** | **≥ 0.30** |
| Sensitivity low | ≥ 0.20 |
| Sensitivity high | ≥ 0.40 |

No post-hoc threshold selection after seeing GNN/RWR metrics.

---

## 5. Positive counts after filtering (2026-07-23)

| Disease | Definition | thr≥0.20 | thr≥0.30 (primary) | thr≥0.40 |
|---|---|---|---|---|
| UC | `genetic_association` only | 502 | **396** | 284 |
| UC | genetics + somatic (sensitivity) | 568 | 463 | 317 |
| IC/BPS | `animal_model` only | 139 | **38** | 6 |

IC/BPS positives at primary threshold (symbols, runtime):  
`SLC3A1, HTR3A, UPK3A, KCNMA1, APRT, CALD1, GATA5, UMOD, CHRM3, HOXA13, AHR, SCARA5, UOX, CHRNB4, PAPOLB, MEP1A, TAC1, HPSE2, CHRNA3, PAX2, DEFB1, TRPV4, MAZ, UPK2, GSTK1, CD209, CLEC4M, MYOCD, PI4K2A, UPK1A, KRT7, GEN1, ARHGDIA, ACTG2, DST, COL4A3, MYH11, ARIH2` (*n*=38; none co-carry `clinical` at this threshold).

Unfiltered overall OT score for comparison (not used): IC thr≥0.3 → *n*=3; UC thr≥0.3 → *n*=174 — illustrating how drug/literature channels dominate naive scores.

---

## 6. Track decision under label scarcity (A2.4)

- **BENCHMARK (UC):** *n*=396 genetics positives at thr≥0.30 → **supervised prioritization remains primary** for H1.  
- **TARGET (IC/BPS):** *n*=38 animal-model positives at thr≥0.30 with **no human genetics channel**. This is too scarce for a credible supervised primary claim under Mondrian CP / few-shot evaluation.

**Decision (frozen):** For the TARGET track, the **primary analysis is the unsupervised / label-free comparator** (signature-seeded propagation + SSL embeddings). Supervised few-shot fine-tuning on the 38 animal-model labels is **secondary / exploratory only**. We do **not** relax filters to manufacture more labels.

---

## 7. Forbidden label sources (analysis plan §2.3)

1. DE genes / signature genes as ground-truth disease labels.  
2. Open Targets `clinical` / `clinical_precedence` / ChEMBL-derived known-drug evidence when drug rediscovery is an evaluation endpoint.  
3. Open Targets `literature` / `europepmc` / `genetic_literature` in the **primary** label definition.  
4. Open Targets `rna_expression` / `expression_atlas` as labels when cohort expression is used as features or signature.
