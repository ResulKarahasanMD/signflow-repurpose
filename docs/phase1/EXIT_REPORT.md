# Phase 1 Exit Report

## Scope

Built Phase 1 infrastructure only for `SignFlow-Repurpose`:

- Environment and packaging via `uv`
- Frozen YAML config encoding Phase 0.6 constants
- Provenance manifest writing
- Resource verification harness
- CLI and smoke command
- Offline unit tests plus one live verification pass

Explicitly not implemented in Phase 1:

- Differential expression
- Graph construction
- Models
- Bulk data downloads

## Files Created

- `src/signflow/__init__.py`
- `src/signflow/config.py`
- `src/signflow/provenance.py`
- `src/signflow/resources.py`
- `src/signflow/cli.py`
- `src/signflow/pipeline.py`
- `src/signflow/mapping.py`
- `configs/default.yaml`
- `configs/benchmark_uc.yaml`
- `configs/target_icbps.yaml`
- `scripts/run_verify_resources.py`
- `scripts/run_smoke.py`
- `tests/test_config.py`
- `tests/test_plan_consistency.py`
- `tests/test_provenance.py`
- `tests/test_resources.py`
- `tests/test_mapping.py`
- `tests/test_seeding.py`
- `tests/test_cli.py`
- `environment.yml`
- `uv.lock`

## Commands Run

```bash
uv venv --python 3.12
uv lock
uv sync
uv run pytest -m 'not live'
uv run python scripts/run_verify_resources.py
```

## Test Results

- Offline pytest: `15 passed, 1 deselected, 0 failed`
- Live verification command: `0` exit status after verifier fix
- IDE diagnostics: no linter errors found in edited files

## Live Resource Verification Summary

- Open Targets: OK
  - UC genetics-only counts matched frozen values: `thr_0.20=502`, `thr_0.30=396`, `thr_0.40=284`
  - UC genetics + somatic sensitivity at `thr_0.30`: `463`
  - IC/BPS animal-model historical counts matched frozen values: `thr_0.20=139`, `thr_0.30=38`, `thr_0.40=6`
  - `animal_model` remains recorded as dropped for TARGET labels
- GEO E-utils: OK (`count=8`)
- OmniPath: OK
- SIGNOR: OK
- STRING: OK (`12.0`)
- DGIdb GraphQL: OK
- ChEMBL: OK
- ClinicalTrials.gov v2: OK (`totalCount=238`)
- decoupler import/liveness: OK
- GTEx best-effort: OK
- HPA best-effort: OK

## Provenance and Data Handling

- `run_manifest.json` is written by CLI commands with config hash, resolved config, package versions, seeds, input SHA-256, resource records, and git commit/dirty state.
- Dirty repositories are refused by default unless `--allow-dirty` is supplied.
- No GEO expression matrices were downloaded.
- No analysis data were downloaded into `data/raw` or committed to the repository.
