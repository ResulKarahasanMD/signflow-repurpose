# Phase 6 rescope decision (binding)

**Decided:** 2026-07-24T075949Z  
**Artifact:** `docs/phase0/PHASE6_RESCOPE_DECISION.json`  
**Amendment:** P5.5.B in `ANALYSIS_PLAN_AMENDMENTS.md`

## Inputs (pre-registered facts)

| Input | Value |
|---|---|
| BENCHMARK BH-FDR survivors | **0** |
| BENCHMARK enrichment@20 / @50 | **0** / **0** |
| TARGET H7 mean top-20 Jaccard | **0.051** (threshold **0.20**) |
| TARGET H7 decision | **non_informative** |

## Decision

1. **TARGET failed H7** → do **not** run convergence ranking on a non-informative list. Run Phase 6 axes as **descriptive characterisation of the frozen seed module only** (tissue/cell appropriateness via GTEx/HPA and GSE310557; L1000 as replication versus PMID 39723619). Report axis disagreement **without a fused rank**.
2. **Always run H6 on BENCHMARK** as pre-registered (k=20 primary, bootstrap over the gold set) — answerable even with weak signal; a null is itself the result.
3. **Do not** add axes or loosen thresholds to manufacture survivors.
4. Carry forward outstanding mandatory TARGET defences: **GSE238208 inflammation-specificity** and the **permuted-signature null** (Phase 1.5 §0.2).

## Mode

`descriptive_characterisation_of_frozen_seed_module_only`

**Do not begin Phase 6 implementation until this decision and the §C GAP reframing are committed.**
