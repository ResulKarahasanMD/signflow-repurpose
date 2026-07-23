# Phase 4 Exit Report

**Date (UTC):** 2026-07-23  
**Output:** `/Users/rak/bio-research/data/processed/phase4/20260723T204931Z`  
**Phase 3 root:** `/Users/rak/bio-research/data/processed/phase3/20260723T203525Z`  
**Splits hash:** `52f3f0d6cf98c6ff14afa2a3e14371b57ea61e8df15be08dfb5119f16a77a584` (loaded from Phase 3; not re-split)  
**Drug data touched:** `False`  
**Seeds:** {0,1,2,3,4} reported as mean ± sd on primary stratified fold (degree-matched test labels).

## Hub confound (visible in results)

Positives mean degree **15.49** vs graph mean **8.72**. Degree-matched negatives control hub bias in regime 1; PU regime does not.

## Model × control-graph AUPRC (degree-matched primary split)

| Model | Graph | AUPRC mean±sd | AUPRC bootstrap CI (seed 0) | AUROC mean±sd |
|---|---|---:|---|---:|
| random | real | 0.476 ± 0.016 | 0.373–0.600 | 0.483 ± 0.030 |
| random | features_only | 0.476 ± 0.016 | 0.373–0.600 | 0.483 ± 0.030 |
| degree | real | 0.459 ± 0.000 | 0.380–0.628 | 0.424 ± 0.000 |
| degree | rewired | 0.459 ± 0.000 | 0.380–0.628 | 0.424 ± 0.000 |
| degree | edge_shuffled | 0.459 ± 0.000 | 0.380–0.628 | 0.424 ± 0.000 |
| degree | sign_shuffled | 0.459 ± 0.000 | 0.380–0.628 | 0.424 ± 0.000 |
| pagerank | real | 0.453 ± 0.000 | 0.380–0.611 | 0.419 ± 0.000 |
| pagerank | rewired | 0.461 ± 0.000 | 0.385–0.633 | 0.423 ± 0.000 |
| pagerank | edge_shuffled | 0.452 ± 0.000 | 0.371–0.617 | 0.419 ± 0.000 |
| pagerank | sign_shuffled | 0.453 ± 0.000 | 0.380–0.611 | 0.419 ± 0.000 |
| logistic | real | 0.641 ± 0.000 | 0.515–0.767 | 0.669 ± 0.000 |
| logistic | features_only | 0.641 ± 0.000 | 0.515–0.767 | 0.669 ± 0.000 |
| hgb | real | 0.599 ± 0.000 | 0.485–0.734 | 0.579 ± 0.000 |
| hgb | features_only | 0.599 ± 0.000 | 0.485–0.734 | 0.579 ± 0.000 |
| mlp | real | 0.586 ± 0.064 | 0.440–0.687 | 0.571 ± 0.071 |
| mlp | features_only | 0.586 ± 0.064 | 0.440–0.687 | 0.571 ± 0.071 |
| label_propagation | real | 0.533 ± 0.000 | 0.447–0.678 | 0.518 ± 0.000 |
| label_propagation | rewired | 0.543 ± 0.000 | 0.418–0.714 | 0.528 ± 0.000 |
| label_propagation | edge_shuffled | 0.497 ± 0.000 | 0.407–0.642 | 0.463 ± 0.000 |
| label_propagation | sign_shuffled | 0.533 ± 0.000 | 0.447–0.678 | 0.518 ± 0.000 |
| rwr | real | 0.516 ± 0.000 | 0.415–0.663 | 0.467 ± 0.000 |
| rwr | rewired | 0.495 ± 0.000 | 0.392–0.658 | 0.439 ± 0.000 |
| rwr | edge_shuffled | 0.456 ± 0.000 | 0.369–0.614 | 0.412 ± 0.000 |
| rwr | sign_shuffled | 0.516 ± 0.000 | 0.415–0.663 | 0.467 ± 0.000 |
| rgcn | real | 0.559 ± 0.054 | 0.430–0.673 | 0.543 ± 0.047 |
| rgcn | rewired | 0.559 ± 0.052 | 0.416–0.642 | 0.508 ± 0.047 |
| rgcn | edge_shuffled | 0.527 ± 0.041 | 0.399–0.664 | 0.492 ± 0.059 |
| rgcn | sign_shuffled | 0.559 ± 0.054 | 0.430–0.673 | 0.543 ± 0.047 |
| gatv2 | real | 0.580 ± 0.042 | 0.462–0.738 | 0.576 ± 0.048 |
| gatv2 | rewired | 0.537 ± 0.028 | 0.401–0.655 | 0.531 ± 0.049 |
| gatv2 | edge_shuffled | 0.523 ± 0.056 | 0.360–0.565 | 0.509 ± 0.058 |
| gatv2 | sign_shuffled | 0.574 ± 0.047 | 0.451–0.737 | 0.568 ± 0.055 |
| ssl_gatv2 | real | 0.567 ± 0.039 | 0.455–0.722 | 0.541 ± 0.048 |
| ssl_gatv2 | rewired | 0.529 ± 0.022 | 0.400–0.652 | 0.504 ± 0.025 |
| ssl_gatv2 | edge_shuffled | 0.546 ± 0.070 | 0.373–0.598 | 0.523 ± 0.070 |
| ssl_gatv2 | sign_shuffled | 0.568 ± 0.039 | 0.446–0.711 | 0.545 ± 0.048 |

### Features-only including `log_degree` (honest comparator)

**Logistic regression on unscaled Phase 3.5 features (includes `log_degree`) AUPRC = 0.641 (CI 0.515–0.767).** This exceeds RWR (0.516) and GATv2 (0.580) on the real graph under the same primary split. The graph did not beat node-level statistics plus connectivity under this comparator.

## H1 verdict (W1 honesty rule — binding)

W1: RWR AUPRC CI overlaps GNN/gatv2 CI on primary split (RWR 0.516 vs GNN 0.580) — conclude no demonstrated GNN advantage under the pre-registered rule.

**Verdict:** `no_gnn_advantage`  
Relational GNN point estimate (GATv2 0.580) exceeded RWR (0.516), but bootstrap CIs overlap → **no demonstrated GNN advantage**. No further tuning was performed.

## Sign-shuffle

- GATv2: Signed model performs indistinguishably on sign-shuffled graph — not using sign information in a detectable way.
- R-GCN: Signed model performs indistinguishably on sign-shuffled graph — not using sign information in a detectable way.

## H2 (SSL vs from-scratch, matched budget)

| | AUPRC mean | sd |
|---|---:|---:|
| Scratch GATv2 | 0.580 | 0.042 |
| SSL→finetune GATv2 | 0.567 | 0.039 |

n_in_graph_positives = 246. With ~246 in-graph positives this comparison is underpowered; report CI/seed-spread width rather than over-interpreting the point estimate. No clear SSL transfer gain within seed spread.

## Mondrian conformal (gene stage; logistic scores; cal=val, test=test)

Non-claim (verbatim):  
*Conformal coverage is over the gene/drug population under the labelled-data distribution. It does not quantify sampling uncertainty over patients, and does not imply the disease signature generalises to new IC/BPS cohorts.*

| α | Empirical coverage | Mean set size | Useful (size < 1.8)? |
|---|---:|---:|---|
| 0.10 | 0.939 | 1.808 | False |
| 0.05 | 0.970 | 1.889 | False |
| 0.20 | 0.869 | 1.515 | True |

At α=0.1 coverage is high but mean set size ≈1.81 — a covering but largely uninformative set (**negative utility result** at the primary α).

## TARGET unsupervised propagation

Primary and sensitivity seed RWR rankings emitted under `target/` for real / rewired / edge-shuffled / sign-shuffled graphs. Top-50 overlap primary∩sensitivity on real graph: **50**. Seed replication annotations carried into ranking tables.

## Plain statement

**The relational GNN did not demonstrate an advantage over RWR under the pre-registered W1 rule (CI overlap). Features-only logistic regression including `log_degree` was the strongest degree-matched primary-split model in this run.**

## Phase 5 gate

Do not begin Phase 5 until this report is committed and H1 verdict recorded (done here, including the negative W1 conclusion). Proximity substrate already frozen in Phase 3.5 A4.
