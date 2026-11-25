# SignFlow-Repurpose

A pre-registered, controlled study of signed-network drug prioritization, with
**interstitial cystitis / bladder pain syndrome (IC/BPS)** as the understudied target
application and **ulcerative colitis (UC)** as the methods benchmark.

> **Hypothesis generation only.** Scores are not binding, PK/PD, causality, or clinical efficacy. **No treatment recommendations.**

**The main result is negative.** On the UC benchmark the signed-flow score separates known
UC therapeutics from a disease-agnostic approved-drug universe only weakly, puts none of
them in the top 20, and no drug survives FDR control. On IC/BPS the pre-registered
robustness check (H7) is non-informative, so no ranked candidate list is reported.
The repository is published so that the plan, the code and the null results can be checked.

---

## Results at a glance

All numbers are copied from the phase reports in `docs/`; each report names the script that produced it.

| Question | Result | Source |
|---|---|---|
| UC benchmark, signed-flow score vs. 29 gold drugs in a pool of 2,833 | AUPRC 0.035 [0.017, 0.069] (chance 0.0102); AUROC 0.665 [0.537, 0.785] | `docs/phase5/EXIT_REPORT.md` §A1 |
| Top-end enrichment | 0 gold drugs in the top 20 (0.20 expected by chance); 0 survivors after BH-FDR | `docs/phase5/EXIT_REPORT.md` §A1 |
| Does the real network matter? | Degree-preserving rewired, edge-shuffled and sign-shuffled graphs give AUPRC 0.033, 0.027 and 0.027 with overlapping intervals | `docs/phase5/FINAL_DIAGNOSTICS.md` §A1 |
| GNN vs. random walk with restart for disease-gene prioritization (H1) | RWR AUPRC 0.524 [0.466, 0.594] vs. GATv2 0.590 [0.529, 0.658]; verdict `no_gnn_advantage` | `docs/phase4/EXIT_REPORT.md` |
| IC/BPS candidate-list stability (H7) | mean pairwise top-20 Jaccard 0.051, threshold 0.20: non-informative | `docs/phase0/PHASE6_RESCOPE_DECISION.md` |
| Convergence of evidence axes on the benchmark (H6) | enrichment@20 = 0.0 [0.0, 0.0]; null not rejected. Only the network-flow axis had data | `docs/phase6/EXIT_REPORT.md` §2 |

Two earlier runs were voided and are kept in the reports next to their corrections:
the first Phase 4 run (undertrained, single fold) and the first Phase 5 pools
(dominated by gold drugs). See `docs/phase4/EXIT_REPORT.md` and `docs/phase5/EXIT_REPORT.md`.

**Where the numbers come from.** Every figure above was produced by the runs of July 2026 and is quoted
from the phase reports written at the time. Nothing was re-run for this publication; only the offline
test suite was run on the published tree (see *Running*).

## Limitations

- The IC/BPS cohorts are very small (5 vs. 6 and 13 vs. 3 samples after one sample per patient). This rules out any state-of-the-art methods claim.
- The gold set on the benchmark is small (29 drugs), so intervals are wide.
- The L1000 replication axis and the tissue/cell-type axis of Phase 6 were **not acquired**; nothing was imputed in their place.
- The results were produced in July 2026 with the environment described in the phase reports. They have not been re-run for this release; only the offline test suite was run (see below).

## Phase status

| Phase | Content | Status | Report |
|---|---|---|---|
| 0 / 0.5 / 0.6 | Gap statement, data gate, analysis plan, pre-registration | Complete (2026-07-23) | `docs/phase0/` |
| 1 | Environment, frozen config, provenance harness | Complete | `docs/phase1/EXIT_REPORT.md` |
| 2 | Expression acquisition, differential expression, meta-analysis | Complete | `docs/phase2/EXIT_REPORT.md` |
| 3 | Signed network construction, features, splits, control graphs | Complete | `docs/phase3/EXIT_REPORT.md` |
| 4 / 4.5 | Disease-gene prioritization models (first run voided and corrected) | Complete | `docs/phase4/EXIT_REPORT.md` |
| 5 | Drug scoring against a frozen disease-agnostic pool, diagnostics | Complete | `docs/phase5/` |
| 6 | Descriptive characterisation under the recorded rescope | Complete, two axes not acquired | `docs/phase6/EXIT_REPORT.md` |

## Pre-registration

The analysis plan was frozen before any Phase 2 data download.
[`docs/phase0/PREREGISTRATION.md`](docs/phase0/PREREGISTRATION.md) lists the frozen commits, the annotated tag
`prereg-phase0.5` and two public gists that serve as the external timestamp (2026-07-23).
Deviations are logged in [`docs/phase0/ANALYSIS_PLAN_AMENDMENTS.md`](docs/phase0/ANALYSIS_PLAN_AMENDMENTS.md);
[`docs/phase0/FROZEN_PLAN.yaml`](docs/phase0/FROZEN_PLAN.yaml) is the machine-readable source of every frozen constant.

**Former name.** This repository was developed privately under the name `bio-research` and is published
as `signflow-repurpose` with its original commit history. The two gists and the commit hashes in the
pre-registration record refer to that repository; because the history is unchanged, those hashes resolve here:

```bash
git rev-parse 2b1a39af2a2cb6ee9595bc1b8f251d4210a389ab   # plan-content freeze commit
git rev-parse d2e4bef71b7f924db761bbfb9c10603df82b17a2   # Phase 0.6 commit
git rev-parse prereg-phase0.5^{}                          # d9751f90965b17d1ffe7c1138975ddc2afcfc5f9
```

| Doc | Role |
|---|---|
| [`docs/phase0/GAP_STATEMENT.md`](docs/phase0/GAP_STATEMENT.md) | Re-scoped novelty (not TxGNN/CARNIVAL priority) |
| [`docs/phase0/METHODS_PRIOR_ART.md`](docs/phase0/METHODS_PRIOR_ART.md) | Per-contribution differentiation |
| [`docs/phase0/LABEL_PROVENANCE.md`](docs/phase0/LABEL_PROVENANCE.md) | Frozen Open Targets allowlists, threshold ≥ 0.30 |
| [`docs/phase0/LABEL_ADJUDICATION.md`](docs/phase0/LABEL_ADJUDICATION.md) | Label adjudication record |
| [`docs/phase0/DATA_AVAILABILITY.md`](docs/phase0/DATA_AVAILABILITY.md) | Cohort gate and planned control roles |
| [`docs/phase0/ANALYSIS_PLAN.md`](docs/phase0/ANALYSIS_PLAN.md) | Pre-registered plan |
| [`docs/phase0/LITERATURE_REVIEW.md`](docs/phase0/LITERATURE_REVIEW.md) | Disease-specific literature |

## Running

Python ≥ 3.12 and [uv](https://docs.astral.sh/uv/). Phase 2 calls limma through `rpy2`, so it also needs R.

```bash
uv sync
uv run pytest -m 'not live'        # offline tests
uv run python scripts/run_smoke.py
```

`scripts/run_phase*.py` are the phase runners; they must be run in order, and each phase has its own gate.
Run manifests are written under `runs/`.

## Data

No data is stored in this repository; `data/` holds only placeholders. The runners download from the
original providers, and each provider's own terms apply:

- NCBI GEO: `GSE11783`, `GSE57560`, `GSE238208` (IC/BPS); `GSE75214`, `GSE87466` (UC)
- Open Targets Platform (disease-gene labels, drug candidates)
- OmniPath, SIGNOR and STRING v12 (networks)
- ChEMBL and DGIdb (drug-target links, from Phase 5 on)

Access dates and query counts of the Phase 0 data gate are in `provenance/phase0_query_log.json`.

## Defaults

| Track | ID | Primary analysis |
|---|---|---|
| BENCHMARK | `MONDO_0005101` (ulcerative colitis) | Supervised H1 (filtered genetics labels) |
| TARGET | `MONDO_0018301` (IC/BPS) | Unsupervised comparator (H7 descriptive) |

## License

Code and documents: [MIT](LICENSE). Third-party data are not redistributed and keep their own licenses.
