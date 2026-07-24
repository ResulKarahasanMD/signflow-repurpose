# Phase 5 Exit Report (rebuilt)

**Run stamp (UTC):** 20260724T004001Z  
**Output:** `/Users/rak/bio-research/data/processed/phase5/20260724T004001Z`  
**Phase 3 root:** `/Users/rak/bio-research/data/processed/phase3/20260723T203525Z`  
**Drug data touched:** `True`  
**Pool SHA-256:** `f748e4b93a24bd491d6040b0797ed081adc9b010484d19532ff3b4e4261d294c`

## Status

Prior rediscovery stamps are **VOID/SUPERSEDED** (candidate-pool artifacts):

| Stamp | Status | Why |
|---|---|---|
| `20260723T235004Z` | VOID | 30 drugs / 29 gold; precision@20=1.0 is arithmetic |
| `20260723T235350Z` | SUPERSEDED | OT UC disease-conditioned candidates, not disease-agnostic |

No rediscovery claim from those stamps enters the manuscript. This stamp is the rebuilt result under the frozen disease-agnostic pool.

## Pool construction

| Item | Value |
|---|---|
| Rule | ChEMBL `max_phase=4`, disease-agnostic ∪ frozen UC gold |
| Disease conditioning | **forbidden** |
| Universe size | **4225** |
| Gold in universe | **29** (prevalence **0.686%**) |
| With ≥1 in-substrate target | **2833** (floor 500 ✓) |
| Gold-prevalence ceiling | 10% (observed ≪ ceiling ✓) |
| Access date | 2026-07-24 |

## Coverage (headline)

- Mapped drug–target summary rows: **3197**; ≥1 primary-substrate target: **2833**; coverage misses: **364**
- Gold eligible / misses: **29 / 0**
- Primary-substrate target-class counts (taxonomy includes `ion_channel`):

| Class | Full pool | Gold |
|---|---:|---:|
| gpcr | 35 | 3 |
| ion_channel | 83 | 1 |
| kinase | 30 | 4 |
| nuclear_receptor | 76 | 10 |
| transporter | 39 | 7 |
| other_signalling_enzyme | 1649 | 181 |

Receptor-class coverage is limited relative to other signalling enzymes — structural limitation of the curated causal substrate; primary substrate unchanged (`phase5_poor_coverage_policy`).

## Signed-propagation substrate connectivity

After `exclude_from_signed_propagation` (ambiguous + unsigned_directed removed):

| Metric | Primary signed substrate |
|---|---:|
| Causal edges before exclusion | 33 797 |
| Edges excluded | 11 599 (**34.3%**) |
| Signed edges retained | 22 198 |
| Nodes (pre-exclusion causal set) | 7 511 |
| Components | **1 313** |
| Largest component | **6 089** |
| Isolated nodes | **1 232** |

Fragmentation after sign-consistent edge exclusion is a structural finding about signed causal networks and is reported, not worked around.

## Signed vs unsigned (co-primary)

Pre-registered primary = **signed**; co-primary = **unsigned** (same directed edges; edge sign ignored). Signed is retained even if it underperforms.

| Metric | Signed | Unsigned | Chance |
|---|---:|---:|---:|
| enrichment@20 | 0.0 | 0.0 | — |
| enrichment@50 | 0.0 | 0.0 | — |
| AUROC | 0.665 | 0.592 | 0.5 |
| AUPRC | 0.0350 | 0.0351 | 0.0102 |
| Paired AUPRC Δ (signed−unsigned) | −0.0002 [−0.0068, 0.0053] | | |

**Decision:** `sign_contributed_no_measurable_benefit_to_drug_scoring` — consistent with the Phase 4 sign control (topology used; sign not detectably used). Signed remains the reported primary.

**Exposure note (freeze):** this co-primary was frozen after Phase 4 sign controls and after the voided 30-drug ranking was visible; the freeze is motivated by the Phase 4 control, not by interpretable drug ranks. Exposure is recorded in `ANALYSIS_PLAN_AMENDMENTS.md` (P5.3).

## Scoring (signed primary)

- Scored drugs: **2833**
- Null: 1000 degree-matched permutations; BH-FDR across scored drugs
- Spearman(score, mean target degree): ρ=**0.047**, p=**0.012** (weak positive association on the full pool; W3 not dominated by degree alone)

### Rediscovery metrics (chance baselines beside figures)

| k | precision@k | enrichment@k | recall@k | E[gold in top-k] at chance |
|---:|---:|---:|---:|---:|
| 20 | 0.000 | 0.000 | 0.000 | 0.20 |
| 50 | 0.000 | 0.000 | 0.000 | 0.51 |
| 100 | 0.050 | **4.88** | 0.172 | 1.02 |
| 250 | 0.056 | **5.47** | 0.483 | 2.56 |

Pool-wide AUROC = **0.665** [0.537, 0.785]; AUPRC = **0.035** [0.017, 0.069] (chance AUPRC = 0.010). Enrichment@20/@50 is 0 — expected to be noisy at ~1% prevalence (chance gold in top-20 ≈ 0.2).

### Gold ranks (never forced; first 10)

| Rank | Drug | Group | Score | BH-FDR |
|---:|---|---|---:|---:|
| 89 | BALSALAZIDE DISODIUM | balsalazide | 0.241 | 0.153 |
| 90 | OLSALAZINE SODIUM | olsalazine | 0.241 | 0.155 |
| 91 | MESALAMINE | mesalamine | 0.241 | 0.074 |
| 97 | UPADACITINIB | upadacitinib | 0.236 | 0.153 |
| 100 | UPADACITINIB HEMIHYDRATE | upadacitinib | 0.236 | 0.210 |
| 129 | GUSELKUMAB | guselkumab | 0.203 | 0.463 |
| 130 | MIRIKIZUMAB | mirikizumab | 0.203 | 0.296 |
| 131 | RISANKIZUMAB | risankizumab | 0.203 | 0.498 |
| 161 | USTEKINUMAB | ustekinumab | 0.185 | 0.074 |
| 163 | ADALIMUMAB | adalimumab | 0.182 | 0.074 |

## FDR statement

**0 drugs** survive BH-FDR < 0.05 on the primary signed-flow score under the degree-matched null. Computational enrichment only — not a clinical recommendation.

## TARGET

Primary seed module scored on the same pool (n_scored=2833). TARGET remains unsupervised / hypothesis-generating; no clinical claim. Seed provenance / replication annotations carried from Phase 3.

## Controls

- Unrelated negative-control drug list under `scores/benchmark/controls/`
- Degree-matched null ≥1000; BH-FDR as above
- Signed substrate connectivity reported above
- Control-graph full re-scoring (rewired / edge-shuffled / repaired sign-shuffled) remains available via Phase 3 control graphs; this stamp’s primary defence is the degree-matched null + signed-vs-unsigned co-primary on the real substrate

## Gate status

- §A void + `drug_data_touched` fix: done
- §B Phase 4.6 corrections: done
- §C signed/unsigned freeze + exposure note: done
- §D rebuilt exit report: **this document**
- **Do not begin Phase 6 until this report is committed**
