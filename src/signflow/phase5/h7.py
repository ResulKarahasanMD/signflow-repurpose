"""H7 signature-bootstrap top-k drug Jaccard futility for the TARGET track."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from signflow.phase2 import PROJECT_ROOT
from signflow.phase2.bootstrap_signature import _dl_two_study_pooled, _effects_for_resample
from signflow.phase2.seeds import load_seed_constants

from . import FROZEN_PLAN, ensure_directory
from .score import CachedFlowScorer, drug_injection_vector


def load_plan() -> dict[str, object]:
    return yaml.safe_load(FROZEN_PLAN.read_text())["plan"]


def _i2_two_study(e1: np.ndarray, s1: np.ndarray, e2: np.ndarray, s2: np.ndarray) -> np.ndarray:
    v1 = s1 ** 2
    v2 = s2 ** 2
    w1 = 1.0 / v1
    w2 = 1.0 / v2
    wsum = w1 + w2
    fixed = (w1 * e1 + w2 * e2) / wsum
    q = w1 * (e1 - fixed) ** 2 + w2 * (e2 - fixed) ** 2
    # I² = max(0, (Q - df) / Q) with df=1 for two studies
    with np.errstate(divide="ignore", invalid="ignore"):
        i2 = np.maximum((q - 1.0) / q, 0.0) * 100.0
    i2 = np.where(np.isfinite(i2), i2, 0.0)
    i2 = np.where(q <= 1.0, 0.0, i2)
    return i2


def pairwise_jaccard(sets: list[set[str]]) -> dict[str, float]:
    if len(sets) < 2:
        return {"mean": float("nan"), "median": float("nan"), "n_pairs": 0}
    vals = []
    for i in range(len(sets)):
        for j in range(i + 1, len(sets)):
            a, b = sets[i], sets[j]
            union = len(a | b)
            vals.append((len(a & b) / union) if union else 1.0)
    arr = np.asarray(vals, dtype=float)
    return {
        "mean": float(np.mean(arr)),
        "median": float(np.median(arr)),
        "lo": float(np.quantile(arr, 0.025)),
        "hi": float(np.quantile(arr, 0.975)),
        "n_pairs": int(len(arr)),
    }


def _select_bootstrap_seeds(
    genes: list[str],
    pooled: np.ndarray,
    i2: np.ndarray,
    *,
    n: int,
    i2_max: float,
) -> list[str]:
    mask = i2 <= i2_max
    order = np.argsort(-np.abs(pooled))
    picked = []
    for idx in order:
        if not mask[idx]:
            continue
        picked.append(genes[idx])
        if len(picked) >= n:
            break
    return picked


def precompute_drug_flows(
    *,
    nodes: list[str],
    signed_edges: pd.DataFrame,
    targets: pd.DataFrame,
    drug_ids: list[str],
    plan: dict[str, object] | None = None,
) -> tuple[np.ndarray, list[str]]:
    """Return (n_drugs × n_nodes) flow matrix for fixed injections on a fixed graph."""

    plan = plan or load_plan()
    scorer = CachedFlowScorer(
        nodes,
        signed_edges,
        module_genes=set(),
        disease_direction=pd.Series(dtype=float),
        alpha=float(plan["phase5_flow_damping"]),
        max_iter=int(plan["phase5_flow_max_iter"]),
        tol=float(plan["phase5_flow_tol"]),
        ignore_sign=False,
        variant="h7_flow_cache",
    )
    flows = np.zeros((len(drug_ids), len(nodes)), dtype=float)
    id_to_group = {cid: g for cid, g in targets.groupby("chembl_id")}
    kept = []
    rows = []
    for cid in drug_ids:
        group = id_to_group.get(cid)
        if group is None:
            continue
        injection, _ = drug_injection_vector(
            group, unknown_policy=str(plan["phase5_unknown_action_policy"])
        )
        if not injection:
            # Still keep a zero flow row so ranks stay aligned when injection empty.
            kept.append(cid)
            rows.append(np.zeros(len(nodes), dtype=float))
            continue
        from .score import signed_flow_adjacency

        flow = signed_flow_adjacency(
            scorer.adjacency_T,
            injection,
            scorer.node_index,
            scorer.n,
            alpha=scorer.alpha,
            max_iter=scorer.max_iter,
            tol=scorer.tol,
        )
        kept.append(cid)
        rows.append(flow)
    if not rows:
        return np.zeros((0, len(nodes)), dtype=float), []
    return np.vstack(rows), kept


def score_drugs_for_module(
    flows: np.ndarray,
    drug_ids: list[str],
    nodes: list[str],
    module_genes: list[str],
    disease_direction: pd.Series,
) -> pd.DataFrame:
    index = {g: i for i, g in enumerate(nodes)}
    idx = []
    ys = []
    for g in module_genes:
        if g not in index or g not in disease_direction.index:
            continue
        y = float(disease_direction.loc[g])
        if not np.isfinite(y):
            continue
        idx.append(index[g])
        ys.append(y)
    if not idx:
        scores = np.zeros(len(drug_ids), dtype=float)
    else:
        module_idx = np.asarray(idx, dtype=int)
        module_y = np.asarray(ys, dtype=float)
        scores = -(flows[:, module_idx] @ module_y)
    return (
        pd.DataFrame({"chembl_id": drug_ids, "score": scores})
        .sort_values(["score", "chembl_id"], ascending=[False, True])
        .reset_index(drop=True)
    )


def run_h7_futility(
    *,
    nodes: list[str],
    signed_edges: pd.DataFrame,
    targets: pd.DataFrame,
    scored_drug_ids: list[str],
    output_dir: str | Path,
    plan: dict[str, object] | None = None,
    n_bootstrap: int | None = None,
) -> dict[str, object]:
    """Apply frozen H7: mean pairwise top-k drug Jaccard across signature bootstraps."""

    plan = plan or load_plan()
    seed_cfg = load_seed_constants()
    output_dir = ensure_directory(output_dir)
    k = int(plan.get("phase5_h7_top_k", 20))
    n_boot = int(n_bootstrap if n_bootstrap is not None else plan.get("phase5_h7_n_bootstrap", 1000))
    threshold = float(plan["h7_futility_j"])

    print(f"[phase5] H7 precompute flows for {len(scored_drug_ids)} drugs", flush=True)
    flows, kept_ids = precompute_drug_flows(
        nodes=nodes,
        signed_edges=signed_edges,
        targets=targets,
        drug_ids=scored_drug_ids,
        plan=plan,
    )

    boot_idx_path = PROJECT_ROOT / "data/processed/phase2/meta/bootstrap/bootstrap_resample_index.json"
    boot_idx = json.loads(boot_idx_path.read_text())
    resamples = boot_idx["resamples"]
    expr = {
        "GSE11783": pd.read_csv(
            PROJECT_ROOT / "data/processed/phase2/GSE11783/expression_gene_level.tsv",
            sep="\t",
            index_col=0,
        ),
        "GSE57560": pd.read_csv(
            PROJECT_ROOT / "data/processed/phase2/GSE57560/expression_gene_level.tsv",
            sep="\t",
            index_col=0,
        ),
    }
    pheno = {
        "GSE11783": pd.read_csv(PROJECT_ROOT / "data/processed/phase2/GSE11783/phenotype_table.tsv", sep="\t"),
        "GSE57560": pd.read_csv(PROJECT_ROOT / "data/processed/phase2/GSE57560/phenotype_table.tsv", sep="\t"),
    }
    shared = sorted(set(expr["GSE11783"].index) & set(expr["GSE57560"].index))
    expr = {k: v.loc[shared] for k, v in expr.items()}
    n_boot = min(n_boot, int(boot_idx["n_bootstrap"]))

    top_sets: list[set[str]] = []
    seed_jaccards: list[set[str]] = []
    for b in range(n_boot):
        if b == 0 or (b + 1) % 100 == 0 or b + 1 == n_boot:
            print(f"[phase5] H7 bootstrap {b + 1}/{n_boot}", flush=True)
        tables = []
        for study in ("GSE11783", "GSE57560"):
            tables.append(
                _effects_for_resample(
                    expr[study],
                    pheno[study],
                    resamples[study][b],
                    study=study,
                )
                .set_index("gene_symbol")
                .loc[shared]
            )
        e1 = tables[0]["effect"].to_numpy(dtype=float)
        s1 = tables[0]["se"].to_numpy(dtype=float)
        e2 = tables[1]["effect"].to_numpy(dtype=float)
        s2 = tables[1]["se"].to_numpy(dtype=float)
        pooled = _dl_two_study_pooled(e1, s1, e2, s2)
        i2 = _i2_two_study(e1, s1, e2, s2)
        seeds = _select_bootstrap_seeds(
            shared,
            pooled,
            i2,
            n=int(seed_cfg["target_seed_n"]),
            i2_max=float(seed_cfg["target_seed_i2_max"]),
        )
        seed_jaccards.append(set(seeds))
        direction = pd.Series(pooled, index=shared)
        ranking = score_drugs_for_module(flows, kept_ids, nodes, seeds, direction)
        top_sets.append(set(ranking.head(k)["chembl_id"].astype(str)))

    drug_j = pairwise_jaccard(top_sets)
    seed_j = pairwise_jaccard(seed_jaccards)
    mean_j = float(drug_j["mean"]) if drug_j["n_pairs"] else float("nan")
    futile = bool(mean_j < threshold) if np.isfinite(mean_j) else True
    decision = "non_informative" if futile else "informative_enough_for_triangulation"
    result = {
        "rule": "H7_signature_bootstrap_top_k_drug_jaccard",
        "k": k,
        "n_bootstrap": n_boot,
        "threshold_j": threshold,
        "mean_pairwise_top_k_jaccard": mean_j,
        "jaccard_summary": drug_j,
        "seed_set_jaccard_summary": seed_j,
        "decision": decision,
        "futile": futile,
        "reporting_rule": (
            "non-informative — stability statistics and axis disagreements only; "
            "no ranked candidate list presented as a finding"
            if futile
            else "H7 passed — triangulation may present ranked candidates with caveats"
        ),
        "n_drugs_scored": len(kept_ids),
        "seed_rule": seed_cfg["target_seed_rule"],
        "seed_n": int(seed_cfg["target_seed_n"]),
        "seed_i2_max": float(seed_cfg["target_seed_i2_max"]),
        "propagated_statistic": boot_idx.get("propagated_statistic"),
        "note": (
            "Flows cached once on the signed primary substrate; bootstrap reseeds the "
            "module only. Binding futility uses drug top-k Jaccard, not seed Jaccard."
        ),
    }
    path = Path(output_dir) / "h7_futility.json"
    path.write_text(json.dumps(result, indent=2, default=str) + "\n")
    result["path"] = str(path)
    return result
