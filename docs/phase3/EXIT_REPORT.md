# Phase 3 Exit Report

**Date (UTC):** 2026-07-23  
**Output root:** `/Users/rak/bio-research/data/processed/phase3/20260723T200209Z` (symlink `data/processed/phase3/latest`)  
**Conflict rule (frozen pre-inspection):** `retain_as_ambiguous`  
**Drug data touched:** `False`  
**Firewall:** Disease-gene labels allowed on **BENCHMARK** only; drug/DGIdb/ChEMBL resources not imported in Phase 3.

## Commands

```bash
uv run python scripts/run_phase26.py   # §A precondition
uv run python scripts/run_phase3.py --track both
uv run pytest -m 'not live'
```

## B1 — Network construction

### Resources (manifest)

Versions/access recorded in `network_manifest.json` under the run directory.

| Resource | Role |
|---|---|
| OmniPath `datasets=omnipath` (human 9606) | Signed directed causal core |
| SIGNOR latest release (`getLatestRelease.php`), TAX_ID=9606 protein–protein | Second causal source |
| STRING v12 physical links (score ≥ 400) | Auxiliary **undirected** layer — separate edge set, never merged into signed causal |

### Conflict rule impact

Rule: **retain_as_ambiguous** (edges both stimulatory and inhibitory kept as `relation_type=ambiguous`, `sign=0`).

| Metric | Value |
|---|---:|
| OmniPath directed inputs | 86189 |
| OmniPath stim∩inhib conflicts | 1418 |
| SIGNOR human protein edges | 20886 |
| Combined causal edges (deduped) | 101372 |
| Sign balance (+1 / 0 / −1) | 0 / 0 / 0 |
| Relation types | {'stimulation': 73927, 'unsigned_directed': 13744, 'inhibition': 11978, 'ambiguous': 1723} |

## B2 — Identifier mapping (HGNC only)

Network symbol QC (shared): matched=8116, alias=37, unmatched=1060, non-coding rejected=72, one-to-many UniProt=25.

Phase 2.5 unmatched asymmetry propagates into the measurable universe before ∩ network: GSE57560 unmatched=16 698 vs GSE11783=2 910 (meta intersection n=16983).

## B3 — Node universe and features

| Track | Nodes (meta∩network) | Causal edges in subgraph | STRING aux edges | Largest CC |
|---|---:|---:|---:|---:|
| TARGET | 7773 | 34322 | 178546 | 7438 |
| BENCHMARK | 7756 | 33797 | 178278 | 7397 |

Node universe = **meta-intersection genes ∩ causal network** (not DE-significance filtered).  
Features unscaled; columns include pooled effect/SE/z/I², per-study effects and **weight shares**, mean expression/variance, CollecTRI TF activity (0 if non-TF), PROGENy pathway activity when gene name matches pathway label else 0, **log-degree**. Missingness after assemble: all reported 0.0 after fill policy for activity.

## B4 — Labels (BENCHMARK) and seeds (TARGET)

- UC genetics positives via Phase 1 path (`genetic_association`, thr≥0.30): **n=396** (asserted vs frozen 396; drift would hard-fail).
- In-graph positives: **246**; degree-matched controls: **246**.
- Degree QC: pos mean=20.35, neg mean=19.88.
- A6.2 carried in `label_meta.json`: genetics channel used for labels excluded from later convergence rank.
- PU framing: unlabeled ≠ negative.
- TARGET: **no labels**. Primary and sensitivity seed sets copied with replication annotation columns; overlap **50/50**.

## B5 — Splits

BENCHMARK: stratified CV (5-fold), hub-stratified CV, leave-one-study-out inductive indices in `benchmark/splits/splits.json` with positive counts per fold. Train/val/test asserted disjoint.

TARGET: supervised splits not emitted (label-less).

## B6 — Graph-information controls

| Control | TARGET verification |
|---|---|
| Degree-preserving rewire | preserves degree exactly = `True` |
| Sign-shuffle | degree equal = `True`; sign multiset preserved = `True` |
| Edge-shuffle | n=45917 |
| Features-only | config JSON recorded |

Control seed = 0 (frozen).

## B7 — PyG

`Data(x, edge_index, edge_attr, edge_type, y)` saved per track under `pyg/`; causal edges directed (not bidirected); STRING auxiliary bidirected; GraphML causal export written; index↔symbol↔UniProt↔Ensembl map alongside.

| Track | PyG nodes | PyG edges (incl. STRING bi-dir) |
|---|---:|---:|
| TARGET | 7773 | 403009 |
| BENCHMARK | 7756 | 401739 |

## Confirmation

- No drug data were touched (`drug_data_touched=false`).
- Phase 4 must not begin until this report is committed and Phase 3 tests are green.

## §A precondition (Phase 2.6)

Completed and committed (`08b77a0`): TARGET meta dominated by GSE11783; sensitivity seed frozen; abundance residual Spearman 0.356→0.353; BENCHMARK concordance recorded as positive-control validation of the metric.


---

## Phase 3.5 Addendum (2026-07-23) — corrections before Phase 4

**Status:** Activity features redefined; post-merge unique-pair conflict rule applied; control graphs regenerated; Phase 5 proximity substrate frozen blind to drug data.  
**Commits:** this addendum; run artifacts under `/Users/rak/bio-research/data/processed/phase3/20260723T203525Z`.  
**Outcome-dependence:** A2 = defect fix (not a rule change). A4 frozen pre-inspection. A1 methodological.

### Numbers that changed vs prior Phase 3 exit report

| Quantity | Prior (incorrect/incomplete) | Corrected |
|---|---|---|
| Causal edge records | 101 372 (multi-record pairs allowed) | **89306 unique ordered pairs** |
| Multi-record pairs collapsed | (unreported; implicit duplicate weighting) | **12352 multi-record → 17769 records removed** |
| Cross-source ±1 conflicts → ambiguous | 0 (rule not applied post-merge) | **162** |
| TARGET track causal edges | 45 917 records (uncollapsed) | **34322 unique pairs** |
| TARGET sign balance (+ / 0 / −) | exit report printed 0/0/0 (bug) | **15498 / 11860 / 6964** |
| Full causal sign balance (+ / 0 / −) | misprinted 0/0/0 | **66617 / 14420 / 8269** |
| `progeny_activity` / `collectri_activity` | dead zero-filled columns (1 and 647 nonzero; missingness reported 0.0) | **Dropped**; replaced by `progeny_score`+`progeny_covered`, `collectri_tf_activity`+`is_tf`+`regulon_target_activity`+coverage indicators |

### A1 — Activity feature density (BENCHMARK n=7756)

See `feature_density_qc.tsv`. Headline nonzero fractions among nodes:
- `progeny_score`: covered **721** nodes (9.3%); structural absent flagged via `progeny_covered`
- `collectri_tf_activity`: **647** nonzero among TFs; `is_tf`=854 (11.0%)
- `regulon_target_activity`: **4499** nonzero (58.0%) — non-TF nodes now carry activity-derived signal
- Near-constant gate (`feature_min_variance=1e-8`): **passed**

Formulations: PROGENy = Σ_p (gene_weight_p × pathway_activity_p); CollecTRI = TF activity for TF nodes + regulon-weighted target activity. Structural absence never silently equals observed zero without an indicator.

### A2 — Post-merge conflict rule (defect fix)

Frozen rule `retain_as_ambiguous` now applied on ordered `(source, target)` **after** OmniPath∪SIGNOR merge. Each pair appears exactly once. The prior 11 884+ multi-record pairs constituted **undocumented edge multiplicity weighting** with no ablation — removed by collapse. Controls regenerated from the unique-pair graph (TARGET rewire degree-exact; sign-shuffle multiset preserved; n_edges=34322).

### A3 — Node-universe provenance

measurable genes (16983) → causal-network symbols (9285, of which 8116 matched + 37 alias) → node universe (7773); binding constraint is the curated network, not expression.

Binding constraint: **curated causal network**, not expression breadth.

### A4 — Phase 5 proximity substrate (frozen blind)

| Key | Value |
|---|---|
| Primary | `signed_causal_only` |
| Sensitivity | `causal_plus_string_unsigned` |
| Min drug targets in substrate | 1 |
| Poor coverage policy | `report_structural_limitation_do_not_switch_substrate` |

### Phase 4 gate

§A corrections committed. Proceed to Phase 4. Do not begin Phase 5 until Phase 4 tests are green, the Phase 4 exit report is committed, and the H1 verdict (including if negative) is recorded.
