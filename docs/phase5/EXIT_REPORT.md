# Phase 5 Exit Report (completions — §A)

**Run stamp (UTC):** 20260724T075949Z
**Output:** `data/processed/phase5/20260724T075949Z`
**Prior rebuild reused:** `data/processed/phase5/20260724T004001Z`
**Phase 3 root:** `data/processed/phase3/20260723T203525Z`
**Drug data touched:** `True`

## Status

Prior rediscovery stamps with gold-dominated / disease-conditioned pools remain **VOID/SUPERSEDED**.
This stamp completes the rebuilt Phase 5 report under the frozen disease-agnostic pool (§A1–A5).

## A1 — Both halves of the result (required pair)

Chance AUPRC on the scored pool = 29/2833 = **0.0102**.

1. The signed-flow score achieves **weak but statistically distinguishable global discrimination** of UC therapeutics against a disease-agnostic approved universe (~3.4× chance AUPRC; AUROC CI lower bound 0.537).
   Observed AUPRC **0.035** [0.017, 0.069] (chance 0.0102); AUROC **0.665** [0.537, 0.785] (chance 0.5).
2. It produces **no top-end enrichment** (expected gold in top-20 under chance = 0.20, observed 0) and **no drug survives BH-FDR** (n_survivors=0), so it has **no practical screening utility**.

Do not report either half alone. Chance baselines sit beside every metric below.

## A2 — Sign vs connectivity (three-arm)

| Arm | Edge set | Sign used | AUPRC (chance) | AUROC | enrichment@20 |
|---|---|---|---:|---:|---:|
| signed | signed_only_reduced | True | 0.03497828632090649 (0.010236498411577833) | 0.6646109006837522 | 0.0 |
| unsigned_matched | signed_only_reduced | False | 0.0350709014423084 (0.010236498411577833) | 0.592029858822372 | 0.0 |
| unsigned_full | full_causal_directed_excluding_ambiguous | False | 0.03344539021361963 (0.010006901311249138) | 0.6030997223591061 | 0.0 |

- Binding sign decision (signed vs unsigned_matched): **sign_contributed_no_measurable_benefit_to_drug_scoring**
- Connectivity contrast: `{"k_primary": 20, "enrichment_unsigned_matched": 0.0, "enrichment_unsigned_full": 0.0, "auroc_unsigned_matched": 0.592029858822372, "auroc_unsigned_full": 0.6021385705150277, "auprc_unsigned_matched": 0.0350709014423084, "auprc_unsigned_full": 0.033576896543434424, "paired_auprc_matched_minus_full": {"mean": 0.0016982769337518875, "lo": -0.0016200784705884199, "hi": 0.005783721992750885, "n_boot": 1000, "width": 0.007403800463339305, "p_le_zero": 0.161, "status": "exploratory_paired_difference", "note": "Exploratory; W1 CI-overlap verdict remains binding."}, "interpretation": "unsigned_matched vs unsigned_full isolates the connectivity cost of excluding ambiguous + unsigned_directed edges from the signed substrate."}`
- Edge counts: signed/matched n_edges=22198; unsigned_full n_edges=33007

signed vs unsigned_matched isolates the sign contribution; unsigned_matched vs unsigned_full quantifies the connectivity cost of `exclude_from_signed_propagation`. Logged as a **defect-driven completion** of the frozen comparison — decision rule unchanged.

## A3 — Differential gold coverage

- Gold coverage rate (mapped summary): **100%** (29/29)
- Non-gold coverage rate (mapped summary): **88.5%** (2804/3168)
- Gold prevalence **acquired universe** 29/4225 = **0.686%** → **scored pool** 29/2833 = **1.024%** (**+49%** relative inflation) — substrate coverage correlates with approval status
- Acquired-pool sensitivity (uncovered drugs counted as misses): AUPRC=0.035 (chance 0.0069); enrichment@20=0.0 (E[gold@20]=0.14)

## A4 — All frozen k values

| k | precision@k | enrichment@k | recall@k | E[gold in top-k] at chance |
|---:|---:|---:|---:|---:|
| 20 | 0.0 | 0.0 | 0.0 | 0.20472996823155665 |
| 50 | 0.0 | 0.0 | 0.0 | 0.5118249205788916 |
| 100 | 0.05 | 4.88448275862069 | 0.1724137931034483 | 1.0236498411577832 |
| 250 | 0.056 | 5.470620689655172 | 0.4827586206896552 | 2.559124602894458 |

## A5 — TARGET track + H7

- H7 mean pairwise top-20 Jaccard: **0.05058371762948879** (threshold 0.2)
- H7 decision: **non_informative** — non-informative — stability statistics and axis disagreements only; no ranked candidate list presented as a finding
- Present ranked list as finding: `False`
- Soft positives (primary seeds × primary substrate): [{"name": "pentosan_polysulfate_sodium", "chembl_ids": ["CHEMBL1201516"], "status": "not_in_scored_pool", "rank": null, "score": null, "note": "Sole OT APPROVAL-stage agent for MONDO_0018301"}, {"name": "amitriptyline", "chembl_id": "CHEMBL1200964", "chembl_ids": ["CHEMBL629", "CHEMBL1200964"], "status": "scored", "rank": 502, "score": 0.00017048280872290666, "pref_name": "AMITRIPTYLINE HYDROCHLORIDE", "note": "Soft positive; report whichever ID is scored (prefer hydrochloride if present)"}, {"name": "cyclosporine", "chembl_id": "CHEMBL160", "chembl_ids": ["CHEMBL160"], "status": "scored", "rank": 20, "score": 0.05209333048283155, "pref_name": "CYCLOSPORINE", "note": "Soft positive (pharmacotherapy NMA signal)"}, {"name": "mirabegron", "chembl_id": "CHEMBL2095212", "chembl_ids": ["CHEMBL2095212"], "status": "scored", "rank": 2432, "score": -0.0001085771566279755, "pref_name": "MIRABEGRON", "note": "Soft positive (symptomatic / investigational context)"}]
- Axis disagreement: primary vs sensitivity-seed top-20 Jaccard=**1.0** (same substrate); primary vs sensitivity-substrate (STRING wired) top-20 Jaccard=**0.905**. Control-graph top-20 Jaccard near 1.0 can reflect score ties / injection dominance — interpret with Spearman vs real (rewired ρ≈0.33; sign-shuffled ρ≈−0.27).
- Control-graph top-20 Jaccard vs real: {"rewired": 1.0, "edge_shuffled": 0.905, "sign_shuffled": 1.0}
- Outstanding mandatory Phase 1.5 defences (not executed here): ['GSE238208_inflammation_specificity', 'permuted_signature_null']

## Coverage headline

- Mapped summary / in-substrate / misses: 3197 / 2833 / 364

## FDR

**0 drugs** survive BH-FDR < 0.05 on the primary signed-flow score.

## Gate

- §A completions: **this document**
- Phase 6 must not begin until H7 decision + §C reframing are committed

