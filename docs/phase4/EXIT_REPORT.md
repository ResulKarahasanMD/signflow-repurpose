# Phase 4 / 4.5 Exit Report

**Corrected run date (UTC):** 2026-07-23 / 2026-07-24
**Corrected output:** `/Users/rak/bio-research/data/processed/phase4/20260723T215010Z`
**Original (voided) output:** `/Users/rak/bio-research/data/processed/phase4/20260723T204931Z`
**Phase 3 root:** `/Users/rak/bio-research/data/processed/phase3/20260723T203525Z`
**Splits hash:** `52f3f0d6cf98c6ff14afa2a3e14371b57ea61e8df15be08dfb5119f16a77a584`
**Drug data touched:** `False`

## Framing

All Phase 4.5 changes are **deviation-from-protocol repairs**, not tuning. Undertraining (15 epochs) and single-fold evaluation departed from the pre-registered plan; corrections restore it. Protocol was frozen in `FROZEN_PLAN.yaml` **before** re-run. No hyperparameter search. Original and corrected runs are reported side by side. The original H1 verdict is **voided as uninformative** (the W1 rule could not have fired positive under the defective protocol).

### Frozen training protocol (pre-rerun)

| Constant | Value |
|---|---:|
| `phase4_gnn_max_epochs` | 500 |
| `phase4_gnn_min_epochs` | 20 |
| `phase4_gnn_early_stop_patience` | 50 (val AUPRC) |
| `phase4_gnn_lr` / `weight_decay` / `hidden` / `dropout` | 1e-3 / 1e-4 / 32 / 0.2 |
| `phase4_ssl_pretrain_steps` | 100 |
| `phase4_ssl_finetune_max_epochs` | 500 (matched to scratch) |
| `phase4_bootstrap_n` | 1000 |
| `phase4_evaluation` | pooled_5fold_oof |
| `phase4_hyperparameter_search` | false |
| `phase4_sign_shuffle_permutes_relation_type` | true |

## Hub confound

Positives mean degree **15.48780487804878** vs graph mean **8.715059308922125** (Phase 3 values: 15.49 vs 8.72).

## Side-by-side: original vs corrected (degree-matched, real graph)

| Model | Original AUPRC (fold 0) | Corrected AUPRC (pooled 5-fold OOF) | Corrected seed-agg CI |
|---|---:|---:|---|
| random | 0.476 | 0.489 ± 0.013 | [0.434, 0.552] |
| degree | 0.459 | 0.503 ± 0.000 | [0.449, 0.570] |
| pagerank | 0.453 | 0.499 ± 0.000 | [0.443, 0.566] |
| logistic | 0.641 | 0.615 ± 0.000 | [0.551, 0.680] |
| hgb | 0.599 | 0.573 ± 0.000 | [0.512, 0.642] |
| mlp | 0.586 | 0.574 ± 0.031 | [0.513, 0.639] |
| label_propagation | 0.533 | 0.543 ± 0.000 | [0.482, 0.611] |
| rwr | 0.516 | 0.524 ± 0.000 | [0.466, 0.594] |
| rgcn | 0.559 | 0.580 ± 0.016 | [0.516, 0.646] |
| gatv2 | 0.580 | 0.590 ± 0.005 | [0.529, 0.658] |
| ssl_gatv2 | 0.567 | 0.579 ± 0.019 | [0.517, 0.647] |

### Corrected full table (pooled OOF; primary = degree-matched)

| Model | Graph | AUPRC mean±sd | Seed-agg CI (n_boot=1000) | PU AUPRC [CI] | AUROC mean±sd |
|---|---|---:|---|---:|---:|
| random | real | 0.489 ± 0.013 | 0.434–0.552 | 0.031 [0.026–0.039] | 0.488 ± 0.016 |
| random | features_only | 0.489 ± 0.013 | 0.434–0.552 | 0.031 [0.026–0.039] | 0.488 ± 0.016 |
| degree | real | 0.503 ± 0.000 | 0.449–0.570 | 0.049 [0.040–0.064] | 0.502 ± 0.000 |
| degree | rewired | 0.503 ± 0.000 | 0.449–0.570 | 0.049 [0.040–0.064] | 0.502 ± 0.000 |
| degree | edge_shuffled | 0.503 ± 0.000 | 0.449–0.570 | 0.049 [0.040–0.064] | 0.502 ± 0.000 |
| degree | sign_shuffled | 0.503 ± 0.000 | 0.449–0.570 | 0.049 [0.040–0.064] | 0.502 ± 0.000 |
| pagerank | real | 0.499 ± 0.000 | 0.443–0.566 | 0.049 [0.040–0.065] | 0.489 ± 0.000 |
| pagerank | rewired | 0.507 ± 0.000 | 0.451–0.574 | 0.050 [0.041–0.066] | 0.507 ± 0.000 |
| pagerank | edge_shuffled | 0.502 ± 0.000 | 0.446–0.569 | 0.050 [0.040–0.066] | 0.496 ± 0.000 |
| pagerank | sign_shuffled | 0.499 ± 0.000 | 0.443–0.566 | 0.049 [0.040–0.065] | 0.489 ± 0.000 |
| logistic | real | 0.615 ± 0.000 | 0.551–0.680 | 0.068 [0.049–0.099] | 0.604 ± 0.000 |
| logistic | features_only | 0.615 ± 0.000 | 0.551–0.680 | 0.068 [0.049–0.099] | 0.604 ± 0.000 |
| hgb | real | 0.573 ± 0.000 | 0.512–0.642 | 0.064 [0.048–0.092] | 0.568 ± 0.000 |
| hgb | features_only | 0.573 ± 0.000 | 0.512–0.642 | 0.064 [0.048–0.092] | 0.568 ± 0.000 |
| mlp | real | 0.574 ± 0.031 | 0.513–0.639 | 0.069 [0.049–0.098] | 0.552 ± 0.023 |
| mlp | features_only | 0.574 ± 0.031 | 0.513–0.639 | 0.069 [0.049–0.098] | 0.552 ± 0.023 |
| label_propagation | real | 0.543 ± 0.000 | 0.482–0.611 | 0.040 [0.034–0.049] | 0.543 ± 0.000 |
| label_propagation | rewired | 0.506 ± 0.000 | 0.449–0.574 | 0.040 [0.033–0.053] | 0.494 ± 0.000 |
| label_propagation | edge_shuffled | 0.535 ± 0.000 | 0.476–0.604 | 0.038 [0.031–0.051] | 0.534 ± 0.000 |
| label_propagation | sign_shuffled | 0.543 ± 0.000 | 0.482–0.611 | 0.040 [0.034–0.049] | 0.543 ± 0.000 |
| rwr | real | 0.524 ± 0.000 | 0.466–0.594 | 0.052 [0.042–0.069] | 0.522 ± 0.000 |
| rwr | rewired | 0.513 ± 0.000 | 0.457–0.580 | 0.050 [0.040–0.068] | 0.502 ± 0.000 |
| rwr | edge_shuffled | 0.500 ± 0.000 | 0.445–0.568 | 0.048 [0.039–0.063] | 0.503 ± 0.000 |
| rwr | sign_shuffled | 0.524 ± 0.000 | 0.466–0.594 | 0.052 [0.042–0.069] | 0.522 ± 0.000 |
| rgcn | real | 0.580 ± 0.016 | 0.516–0.646 | 0.082 [0.057–0.116] | 0.553 ± 0.014 |
| rgcn | rewired | 0.556 ± 0.017 | 0.495–0.621 | 0.093 [0.063–0.133] | 0.547 ± 0.024 |
| rgcn | edge_shuffled | 0.559 ± 0.025 | 0.499–0.624 | 0.099 [0.066–0.140] | 0.551 ± 0.030 |
| rgcn | sign_shuffled | 0.586 ± 0.020 | 0.523–0.650 | 0.081 [0.057–0.115] | 0.555 ± 0.017 |
| gatv2 | real | 0.590 ± 0.005 | 0.529–0.658 | 0.083 [0.058–0.116] | 0.589 ± 0.008 |
| gatv2 | rewired | 0.516 ± 0.008 | 0.459–0.582 | 0.044 [0.033–0.062] | 0.506 ± 0.018 |
| gatv2 | edge_shuffled | 0.539 ± 0.028 | 0.479–0.605 | 0.071 [0.050–0.100] | 0.533 ± 0.028 |
| gatv2 | sign_shuffled | 0.588 ± 0.010 | 0.527–0.654 | 0.082 [0.057–0.115] | 0.586 ± 0.010 |
| ssl_gatv2 | real | 0.579 ± 0.019 | 0.517–0.647 | 0.114 [0.081–0.156] | 0.576 ± 0.026 |
| ssl_gatv2 | rewired | 0.503 ± 0.021 | 0.447–0.569 | 0.065 [0.045–0.096] | 0.482 ± 0.022 |
| ssl_gatv2 | edge_shuffled | 0.538 ± 0.021 | 0.477–0.605 | 0.091 [0.061–0.131] | 0.530 ± 0.021 |
| ssl_gatv2 | sign_shuffled | 0.589 ± 0.020 | 0.528–0.654 | 0.089 [0.062–0.126] | 0.587 ± 0.031 |

### PU regime notes (P4.6.4)

- **Inverted controls (R-GCN PU):** real 0.082 < rewired 0.093 < edge_shuffled 0.099 — control graphs outscore the real graph in the PU regime.
- **Regime reversal:** ssl_gatv2 leads PU at 0.114 [0.081, 0.156]; logistic leads the degree-matched regime.

## H1 verdict (W1 on **corrected** run — binding)

W1 (corrected pooled 5-fold OOF, n_boot=1000): RWR AUPRC=0.524 [0.466,0.594] vs gatv2=0.590 [0.529,0.658]; overlap=True; MDD≈0.128. Verdict=no_gnn_advantage.

**Verdict:** `no_gnn_advantage`
**Interval used:** `mean_across_seeds_of_bootstrap_CI_on_pooled_OOF`
**Minimum detectable AUPRC difference (corrected):** ≈ **0.128**
**Minimum detectable AUPRC difference (original, approx from prior CI widths):** ≈ **0.262**

Observed GATv2−RWR gap ≈ 0.066 < MDD ≈ 0.128 → even after protocol repair, non-overlap would have required a gap larger than observed. The negative is therefore still limited by power, but is no longer an artifact of undertraining / single-fold / n_boot=100.

### Original H1 (voided)

Original: RWR 0.516 vs GATv2 0.580, fold 0, n_boot=100 → `no_gnn_advantage`. **Voided** as uninformative.

## Reframing / chance baseline (P4.6.1)

Chance AUPRC on balanced labels = 0.5.

**CI lower bound > 0.5 (exceed chance):** logistic 0.615 [0.551, 0.680], gatv2 0.590 [0.529, 0.658], rgcn 0.580 [0.516, 0.646], ssl_gatv2 0.579 [0.517, 0.647], hgb 0.573 [0.512, 0.642], mlp 0.574 [0.513, 0.639].

**CI includes 0.5 (at chance):** rwr 0.524 [0.466, 0.594], label_propagation, pagerank, degree, random.

**RWR — the pre-registered structural baseline — is itself at chance.** Substrate sparsity (median degree 3, many components, hub-biased positives) remains a **candidate** explanation, not established.

## Sign-shuffle (repaired control — A4 / P4.6.2)

- **Prior R-GCN identical-to-shuffle conclusion RETRACTED** (inert control: only `sign` shuffled, leaving `edge_type` unchanged).
- **Finding (promoted to GAP_STATEMENT C1/C4 + Limitations):** for both valid signed models, **real ≈ sign_shuffled but real > rewired** — topology is used; sign is not detectably used. Caveat: MDD ≈ 0.128.
  - gatv2: real / sign_shuffled / rewired = **0.590 / 0.588 / 0.516**
  - ssl_gatv2: **0.579 / 0.589 / 0.503**
  - rgcn: real ≈ sign_shuffled (0.580 / 0.586); rewired 0.556 (weaker topology gap than GATv2)

## H2 (SSL vs scratch, matched converged budget)

| | AUPRC mean | sd |
|---|---:|---:|
| Scratch GATv2 | 0.590 | 0.005 |
| SSL→finetune GATv2 | 0.579 | 0.019 |

Matched converged budget (same max_epochs/patience); underpowered with ~246 positives.

## Additional splits (A2 / P4.6.6)

Hub-stratified and LOSO extended beyond logistic to rwr and gatv2. Exact per-model numbers are filled by `scripts/phase46_posthoc.py` into `benchmark/phase46_posthoc.json` / updated `benchmark_results.json`.

### Hub-stratified CV
- logistic: AUPRC 0.651 ± 0.066 over 5 folds (pre-existing).
- rwr / gatv2: see posthoc artifact.

### Leave-one-study-out
- logistic: AUPRC 0.641 ± 0.020 over 2 folds (pre-existing).
- rwr / gatv2: see posthoc artifact.

## Paired-difference bootstrap (P4.6.5; exploratory)

Pre-committed exploratory contrast: gatv2 − rwr on seed-0 pooled OOF with n_boot=1000. **W1 CI-overlap verdict remains binding.** Results written by posthoc script.

## Convergence diagnostics (A1)

- GNN/SSL runs logged: **300**
- Early stopped: **300**
- Hit epoch cap: **0**
- Undertrained (cap without plateau): **0**

All runs early-stopped; none hit the 500-epoch cap without plateauing.

## Mondrian conformal (mean set size only; A7 / P4.6.3)

**Protocol:** pooled 5-fold val-calibration (not fold 0 only). Prior fold-0-only figures (α=0.1 → 1.808 / 0.939) are superseded and retained only for audit.

| α | mean set size (pooled) | empirical coverage (pooled) |
|---:|---:|---:|
| 0.1 | *(filled by phase46_posthoc)* | *(filled by phase46_posthoc)* |
| 0.05 | *(filled by phase46_posthoc)* | *(filled by phase46_posthoc)* |
| 0.2 | *(filled by phase46_posthoc)* | *(filled by phase46_posthoc)* |

Non-claim: Conformal coverage is over the gene/drug population under the labelled-data distribution. It does not quantify sampling uncertainty over patients, and does not imply the disease signature generalises to new IC/BPS cohorts.

## Metric notes (A7)

- ECE computed only on models with actual probabilities; omitted for unbounded scores / PU regime.
- MCC / balanced accuracy reported only in the degree-matched regime with probability threshold 0.5; omitted in PU.
- PU-regime AUPRC is reported in the table above.

## TARGET propagation

Top-50 overlap n = 50 (see `/Users/rak/bio-research/data/processed/phase4/20260723T215010Z/target/target_propagation_summary.json`).

## Firewall

Drug data not touched in Phase 4/4.5. Labels remain BENCHMARK-only for supervised evaluation.

## §B status

UC gold-standard rediscovery set frozen in `configs/benchmark_drugs.yaml` and hashed into `FROZEN_PLAN.yaml` before Phase 5 drug load.

