# Phase 2 Exit Report

**Date (UTC):** 2026-07-23  
**Tracks run:** TARGET (expression → DE → meta → concordance → bootstrap → activity), BENCHMARK (acquisition + phenotype only)  
**R / limma:** 4.6.1 / 3.68.4  
**plan_hash:** `see run manifest`  
**Firewall:** No disease-gene labels, drugs, or networks were loaded in Phase 2 code paths.

## Commands

```bash
export KMP_DUPLICATE_LIB_OK=TRUE OMP_NUM_THREADS=1
export R_HOME=/opt/homebrew/Cellar/r/4.6.1/lib/R
uv run python scripts/run_phase2.py --track target --allow-dirty --n-bootstrap 1000
uv run python scripts/run_phase2.py --track benchmark --allow-dirty
uv run pytest -m 'not live'
```

`--allow-dirty` was required because Phase 2 sources were uncommitted during the live run. Manifests remain append-only under `runs/phase2/`.

## Per-accession sample counts and phenotype mapping

| Accession | Role | Derived counts | Mapping notes |
|---|---|---|---|
| GSE11783 | TARGET primary A | case 5 / control 6 (one sample/patient, **non-ulcer**) | Raw: IC non-ulcer=5, IC ulcer=5, control=6; primary keeps non-ulcer + controls |
| GSE57560 | TARGET primary B | case 13 / control 3 | Agilent **1-colour** (runtime); disease state IC vs control |
| GSE238208 | Inflammation-specificity control | case 50 (HIC) / control 13 (BCG cystitis) | Phenotype only; **forbidden** as IC-vs-healthy DE |
| GSE75214 | BENCHMARK | case 74 / control 11 after exclusion | Frozen contrast: **UC active colon vs control colon**; excluded 109 (CD/Crohn/inactive UC/ileum) |
| GSE87466 | BENCHMARK | case 87 / control 21 | Phenotype from series characteristics |

**Drift vs DATA_AVAILABILITY.md:** GSE11783 and GSE57560 primary counts match frozen expectations (6/5 and 3/13). No count drift stop fired.

## Probe → gene mapping

| Study | probes_in | genes_out | unmapped | rejected non-coding | many-to-one |
|---|---:|---:|---:|---:|---:|
| GSE11783 (GPL570) | 54675 | 21755 | 9557 | 0 | 23363 |
| GSE57560 (GPL16699) | 62976 | 32080 | 8681 | 0 | 22215 |

Collapse rule: **max-mean probe per gene**. Expression source: **author-normalised series matrix** (not CEL RMA); recorded in mapping QC.

## DE and meta-analysis

- Per-study limma `~ group`, moderated t, BH-FDR; standardised effect = Hedges' g with SE.
- GSE11783 DE genes: 21755; GSE57560 DE genes: 32080.
- Primary meta (`signature_meta.tsv`): **DerSimonian–Laird** on Hedges' g; intersection **17529** genes; lost GSE11783-only=4226; lost GSE57560-only=14551.
- Genes with BH-FDR < 0.05 in primary meta: **0** (reportable underpowered dual-cohort result; not adapted post-hoc).

## Concordance (pre-registered result — stated plainly)

Moderate concordance between cohorts; shared signal is partial. Phenotype mix (Hunner-selected vs capacity-stratified) remains a candidate explanation, not established.

| Metric | Value |
|---|---|
| Shared genes | 17529 |
| Pearson r (95% CI) | 0.2570 (0.2431, 0.2707) |
| Spearman r (95% CI) | 0.2772 (0.2634, 0.2908) |
| Directional agreement (nominal p<0.05 in either) | 0.7441 |
| Top-k overlap k=50/100/250 | 0 / 0 / 0 |
| I² median / mean / q75 | 0.00 / 18.85 / 40.31 |

**Top-k overlap of zero** among the largest |effect| genes is itself a publishable finding of discordance at the extremes, despite moderate genome-wide rank correlation.

## Sensitivity (labelled, never primary)

- GSE11783 ulcer-site and patient-mean variants carried through meta.
- Rank-based combination and ComBat-like study mean-centering emitted under `meta/sensitivity/` with confounding headers (6:5 vs 3:13).
- Top-50 absolute-effect overlap vs primary: ulcer ∩ primary = 5; patient-mean ∩ primary = 13.

## Bootstrap

- n=1000 stratified resamples preserving group sizes.
- Propagates Hedges' g → DL meta (not limma per resample; amendment P2.3).
- Degenerate resamples (sum over studies): 116 / 2000 study-level draws.
- Artifact: `meta/bootstrap/bootstrap_effect_summary.tsv` + resample index JSON.

## Activity inference QC (A6.5)

- CollecTRI + PROGENy via decoupler `ulm`; dispersion tables written per study and on the meta signature vector.
- GSE11783 activity score std median ≈ 0.343 across 776 sources.
- **Domain-mismatch limitation (symmetric with DepMap):** PROGENy footprints and CollecTRI regulons are substantially cancer/cell-line derived; activity inference from n≈16 microarray cohorts is low-precision and must not overrule expression-derived results.

## Descriptive QC

Boxplots and PCA under each study `qc/` directory carry the mandated header: descriptive QC — not used for model selection.

## Stop-and-report events during build

1. Initial GSE238208 map keys mismatched series text → fixed to runtime `disease:` strings before DE (P2.1).
2. Initial GSE75214 map missed `disease: CD` → rules-based map + frozen UC-active-colon contrast (P2.5).
3. decoupler≥2 returns `(estimate, pvalue)` tuple → handled; activity no longer aborts DE.

Final TARGET and BENCHMARK runs: **`stop_conditions: []`**.

## Phase 3 gate

Do **not** begin Phase 3 until this exit report is reviewed/committed and Phase 2 offline tests remain green.


---

## Phase 2.5 Addendum (2026-07-23)

**Status:** Corrections A1–A5 implemented and re-run from a **clean tree without `--allow-dirty`**.  
**Commits:** `5d7d2fa` (Phase 2), `cd967df` / `bec55f0` (Phase 2.5).  
**plan_hash:** `f09ab91730402e037888d498b4fc67076741f95a9b91b8ac4648b431fc9c272c`

### Phase 1.5 §A precondition checklist

| Item | Status |
|---|---|
| A1 `FROZEN_PLAN.yaml` single source of truth + plan-consistency value parsing | **Done** (commit `cfd6704`) |
| A2 append-only manifests + hard-failing `git_state` | **Done** |
| A2.4 `plan_hash` (semantic) vs `resolved_config_hash` (full) split | **Done** |

### Numbers retracted as artifacts

| Prior Phase 2 claim | Why artifactual | Replacement |
|---|---|---|
| Primary meta **0 genes** at BH-FDR < 0.05 | Meta pooled unmoderated Hedges' g, discarding limma eBayes | Moderated logFC+SE DL meta: **640** genes at BH-FDR < 0.05 (intersection n=16983) |
| Top-k overlap **0** for k=50/100/250 presented as dramatic discordance | Top-k from full per-study tables (different universes) ranked by noisy abs(g) | Shared-gene only; moderated-t ranking: k50=0 (E=0.15); k100=2 (E=0.59); k250=4 (E=3.68) |
| `rejected_non_protein_coding: 0` | GEO platforms lack biotype; validate path skipped | HGNC: GSE11783 rejected=3743, unmatched=2910, genes_out=17221; GSE57560 rejected=2753, unmatched=16698, genes_out=18402 (was 32080) |

Prior `signature_meta.tsv` archived under `data/processed/phase2/meta/superseded/`.

### Corrected TARGET headline results

**Per-study limma (BH-FDR):**
- GSE11783: 775 at <0.05; 1698 at <0.10 (n=17221)
- GSE57560: 0 at <0.05; 1 at <0.10 (n=18402) — control-limited (n_control=3)

**Primary meta (moderated logFC + moderated SE, DL):** BH-FDR<0.05 = **640** / 16983 shared genes.

**Concordance (shared genes only, moderated ranking):**  
Moderate shared-gene concordance (Spearman r=0.356); top-k moderated overlap=0 vs independence expectation=0.15. Phenotype mix remains a candidate explanation, not established.  
Spearman r=0.3557; directional agreement (nominal)=0.7465.  
Top-k moderated overlaps remain near independence expectation at k=50 and modestly above at k=100/250 — **not** the previously implied absolute discordance from mismatched universes.

**Hedges' g sensitivity:** top-50 overlap with moderated primary = 3 (labelled sensitivity only).

### BENCHMARK

Full DE→meta→concordance→bootstrap path implemented; clean-tree results under `runs/phase2/benchmark/` and `data/processed/phase2/meta/signature_meta_benchmark.tsv`.
