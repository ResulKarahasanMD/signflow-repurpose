# Phase 6 Exit Report — build under the recorded rescope

**Rescope artifact:** `docs/phase0/PHASE6_RESCOPE_DECISION.json`
(`phase6_mode = descriptive_characterisation_of_frozen_seed_module_only`;
`h7_decision = non_informative`; `benchmark_fdr_survivors = 0`).

Phase 6 is built **to the recorded decision**: TARGET axes are a descriptive
characterisation of the frozen seed module (no fused rank), H6 runs on BENCHMARK as
pre-registered, and the two outstanding mandatory TARGET defences are executed. The
"do-not" rules (no added axes, no loosened thresholds, no fused rank to manufacture
survivors) are enforced **in code** (`scripts/run_phase6_target_descriptive.py`,
`scripts/run_phase6_h6_benchmark.py`), not just in prose.

Reproduced by:
- `scripts/run_phase6_h6_benchmark.py` — H6 on BENCHMARK
- `scripts/run_phase6_permuted_signature_null.py` — permuted-signature null defence
- `scripts/run_phase6_gse238208_specificity.py` — GSE238208 inflammation-specificity defence
- `scripts/run_phase6_target_descriptive.py` — TARGET descriptive consolidation + guards

---

## 1. TARGET axes — descriptive characterisation only (§B1)

**H7 NON-INFORMATIVE** (mean pairwise top-20 Jaccard = **0.051**, threshold 0.20).
Every TARGET table below carries this banner. **No fused convergence rank is produced**
(`fused_convergence_rank_computed = false`); `present_ranked_list_as_finding = false`.
The code guard hard-fails if any fused/convergence TARGET artifact exists or if a ranked
list is flagged as a finding under non-informative H7.

**Frozen seed module (primary = sensitivity_concordant; Jaccard 1.0, n=50).** All 50
seeds direction-agree across both cohorts; median I² = 0; the meta is single-cohort
dominated (median GSE11783 weight share ≈ 0.79; GSE57560 near-zero RE weight, P2.6.1).
Top seeds: CHI3L1, CXCL13, WIF1, CCL18, S100A8, MMP9, CXCL1, S100A9, LCN2, … —
an inflammatory / neutrophil / chemokine module.

**Per-axis results and disagreements (descriptive; no fusion):**

| Axis | Status | Result |
|---|---|---|
| Inflammation-specificity (GSE238208) | **executed** | seed module separates Hunner-type IC from BCG cystitis, AUROC **0.94** (see §3a) |
| Network-flow stability (H7) | **executed** | H7 non-informative (0.051); control-graph Spearman vs real: rewired 0.33, edge_shuffled 0.18, sign_shuffled −0.27 |
| Substrate/seed disagreement | **executed** | primary vs sensitivity-seed top-20 Jaccard 1.0; primary vs STRING-wired substrate 0.905 |
| L1000 replication (PMID 39723619) | **not acquired** | framed as replication with agreement/disagreement table + cell-context mismatch; **no data acquired → no imputation** |
| Tissue/cell (GTEx/HPA, GSE310557) | **not acquired** | bladder/colon expression + cell-type support; **no data → no imputation** |
| Genetics | excluded on BENCHMARK (A6.2); TARGET has no genetics channel |
| DepMap | exploratory only, out-of-domain disclaimer |

**Soft positives at their ranks** (primary seeds × primary substrate; never forced):
cyclosporine rank 20, amitriptyline rank 502, mirabegron rank 2432; pentosan polysulfate
`not_in_scored_pool` (no protein-target mechanism — see Phase 5 §A4).

---

## 2. H6 on BENCHMARK — pre-registered (§B2)

Convergence axis set = {network-flow, L1000, tissue/cell}; genetics excluded (A6.2).
Only **network-flow** has data (L1000 and tissue/cell not acquired → no imputation), so
**fewer than two convergence axes are available**: no fused convergence rank is computed,
and the "convergence beats best single axis" advantage is **not evaluable**. H6 is
reported on the single available axis with bootstrap resampling over the gold set.

| k | network-flow enrichment@k [95% CI over gold] | role |
|---:|---|---|
| **20** | **0.0 [0.0, 0.0]** | **primary → NULL** |
| 50 | 0.0 [0.0, 0.0] | secondary |
| 100 | 4.88 [1.67, 7.87] | secondary |
| 250 | 5.47 [3.97, 7.08] | secondary |

**Primary decision: fail to reject the null — no convergence advantage.** Enrichment@20 =
0, exactly the pre-registered likely outcome. Enrichment exists at k=100/250 but is
reported **strictly as secondary** and does **not** shift the primary claim from k=20.

---

## 3. Outstanding mandatory TARGET defences — executed

### 3a. GSE238208 inflammation-specificity (primary biological check)

Contrast: Hunner-type IC (n=50) vs BCG cystitis (n=13), an out-of-sample test against a
different inflammatory bladder condition. Directional seed-module score (signature-signed
mean of z-scored log2 FPKM):

| Seed set | AUROC IC vs BCG | Mann-Whitney p | seed-gene direction agreement | random-set null AUROC (mean / q95) | empirical p | verdict |
|---|---:|---:|---:|---|---:|---|
| primary | **0.937** | 1.5e-6 | 100% | 0.65 / 0.84 | **0.006** | **Hunner-specific** |
| sensitivity_concordant | 0.937 | 1.5e-6 | 100% | 0.65 / 0.85 | 0.001 | Hunner-specific |

The module **separates Hunner-type IC from BCG cystitis** far beyond size-matched random
gene sets: it captures **Hunner-type-IC-specific inflammation**, not generic bladder
inflammation. This is a genuinely positive biological check on the frozen module.
(Descriptive only; H7 non-informative; no candidate list derived.)

### 3b. Permuted-signature null (BENCHMARK)

Signature label-permutation (n_perm=1000; flow computed once per drug, re-dotted against
permuted signatures — R-free and tractable). Confirms whether the ranking collapses when
the disease case/control signal is shuffled.

| Quantity | Value |
|---|---|
| observed AUPRC | 0.033 (chance 0.0102) |
| permutation-null AUPRC (mean / q95) | 0.017 / 0.039 |
| empirical p (observed ≥ null) | 0.073 |
| mean \|Spearman\| permuted-vs-observed | **0.15** |
| ranking order collapses | **yes** |
| danger flag (pipeline reports non-disease-biology) | **NOT raised** |

**PASS (§B3).** Under signature permutation the ranking **order collapses** (scores
decorrelate, mean |Spearman| ≈ 0.15) and AUPRC halves toward chance. The pipeline is
therefore **not** reporting a signature-independent artifact. The residual discrimination
is only marginally above the magnitude-matched null (p ≈ 0.073) — consistent with the
weak-signal / no-screening-utility headline and the §A2 signature-dominance finding.
Given the small gold set (n=29) the null AUPRC has a heavy right tail.

---

## 4. Axis notes (§B4)

- **L1000** — framed as **replication** of PMID 39723619 (signature reversal), with an
  agreement/disagreement table and the cancer-cell-line context mismatch stated. Data not
  acquired in this environment; no imputation.
- **Tissue/cell appropriateness** — GTEx/HPA bladder (TARGET) / colon (BENCHMARK)
  expression plus GSE310557 cell-type support. Data not acquired; no imputation.
- **Genetics** — excluded from convergence on BENCHMARK (A6.2 Option 1; same channel
  family as training labels); TARGET has no genetics channel.
- **DepMap** — exploratory only, out-of-domain disclaimer (non-cancer indications).

---

## 5. Do-not, enforced in code (§B5)

`no_fused_convergence_rank`, `no_added_axes`, `no_loosened_thresholds`,
`no_ranked_list_presented_as_finding_under_non_informative_H7`. The TARGET consolidation
script raises `Phase6GuardError` on any violation; the H6 script suppresses the fused rank
and records the missing-axis / no-imputation reason.

---

## Exit criteria (§C) status

| Item | Status |
|---|---|
| A1 control-graph table + pre-committed interpretation | ✅ `docs/phase5/FINAL_DIAGNOSTICS.md` §A1 |
| A2 invariance/instability panel + Spearman quantification | ✅ §A2 |
| A3 top-50 profile | ✅ §A3 |
| A4 Limitations text | ✅ §A4 |
| Phase 6 built to recorded rescope | ✅ this report |
| Both TARGET defences executed | ✅ §3a GSE238208, §3b permuted-signature null |
| H6 reported whatever it shows | ✅ §2 (null at k=20) |

Analytical work complete. Remaining task: the manuscript.
