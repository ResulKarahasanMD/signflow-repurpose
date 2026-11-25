# CLAUDE.md — bio-research (SignFlow-Repurpose)

## Read this first
In this repository **`README.md` and the pre-registered plan under `docs/phase0/` are binding.** The phase status, data gate, pre-registration and amendments log are there; I do not repeat them here. Do not run any analysis outside the plan — `docs/phase0/ANALYSIS_PLAN.md` is pre-registered, and deviations are logged in `ANALYSIS_PLAN_AMENDMENTS.md`.

## What this is
Controlled signed-network prioritization + multi-evidence triangulation for IC/BPS (the main application); ulcerative colitis is the methods benchmark arm. The package name is `signflow-repurpose`, and the CLI entry point is `signflow` (`src/signflow/cli.py`).

Python ≥3.12, uv (`uv.lock`, `.venv`). pydeseq2 / decoupler / gseapy / torch / torch-geometric / rpy2.

## Commands
```bash
uv sync
uv run pytest                      # pythonpath = ["src"]
uv run python scripts/run_smoke.py
uv run python scripts/run_phase5.py
```

`scripts/run_phase*.py` are the phase runners; do not run them out of order — each phase passes through its own gate.

## Caution
- **Hypothesis generation, not treatment recommendations.** Do not present scores as PK/PD, causality or clinical efficacy.
- Labels have been frozen using the Open Targets allowlist and a threshold of ≥0.30 (`docs/phase0/LABEL_PROVENANCE.md`). Do not use DE results as labels — that would be circular.
- Label gate: IC n=38 → the primary analysis is **unsupervised**. Do not turn it into a supervised one.
- `provenance/`, `runs/`, `logs/` are generated records — do not edit them by hand, and ask before deleting.
- The data is not in the repository (`data/` holds only `.gitkeep`); the phase runners download it. Write new large outputs to a dated directory under `runs/`.
- No direct work on the main branch, and no commit/push without permission.
