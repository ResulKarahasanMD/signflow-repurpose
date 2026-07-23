"""Run all models × control graphs under both evaluation regimes."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from signflow.phase3.features import FEATURE_COLUMNS

from . import FROZEN_PLAN, ensure_directory
from .baselines import (
    feature_matrix,
    fit_predict_tabular,
    label_propagation_scores,
    rwr_scores,
    score_degree,
    score_pagerank,
    score_random,
)
from .conformal import NONCLAIM, empirical_coverage, mondrian_conformal_sets
from .data import labelled_frame, load_frozen_seeds
from .metrics import bootstrap_ci, ci_overlap, compute_metrics
from .models_gnn import train_gnn_scores
from .ssl import pretrain_and_finetune


def _feature_cols_with_degree() -> list[str]:
    return list(FEATURE_COLUMNS)


def _feature_cols_without_graph_free() -> list[str]:
    # still includes log_degree for the honest comparator when requested
    return list(FEATURE_COLUMNS)


def evaluate_scores(
    nodes: list[str],
    y: pd.Series,
    scores: np.ndarray,
    *,
    regime: str,
    seed: int,
    with_ci: bool = True,
) -> dict[str, object]:
    score_map = dict(zip(nodes, scores, strict=False))
    if regime == "degree_matched":
        mask = y >= 0
    elif regime == "pu":
        # positives vs all non-positives
        mask = y.notna()
        y_eval = y.copy()
        y_bin = (y_eval == 1).astype(int)
        genes = y_eval.index[mask].tolist()
        s = np.array([score_map[g] for g in genes], dtype=float)
        yt = y_bin.loc[genes].to_numpy()
        metrics = compute_metrics(yt, s)
        if with_ci:
            metrics["auprc_ci"] = bootstrap_ci(yt, s, metric="auprc", seed=seed)
            metrics["auroc_ci"] = bootstrap_ci(yt, s, metric="auroc", seed=seed)
        return metrics
    else:
        raise ValueError(regime)

    genes = y.index[mask].tolist()
    yt = y.loc[genes].to_numpy().astype(int)
    s = np.array([score_map[g] for g in genes], dtype=float)
    metrics = compute_metrics(yt, s)
    if with_ci:
        metrics["auprc_ci"] = bootstrap_ci(yt, s, metric="auprc", seed=seed)
        metrics["auroc_ci"] = bootstrap_ci(yt, s, metric="auroc", seed=seed)
    return metrics


def run_benchmark_suite(bundle: dict[str, object], output_dir: str | Path) -> dict[str, object]:
    output_dir = ensure_directory(output_dir)
    plan = yaml.safe_load(FROZEN_PLAN.read_text())["plan"]
    seeds = load_frozen_seeds()
    nodes = bundle["nodes"]
    y = bundle["y"]
    features = bundle["features"]
    splits = bundle["splits"]
    # Primary split = stratified fold 0
    primary = splits["stratified_cv"][0]
    labelled = labelled_frame(bundle)
    labelled_genes = labelled["gene_symbol"].tolist()

    def genes_from_idx(idx_list: list[int]) -> list[str]:
        return [labelled_genes[i] for i in idx_list]

    train_genes = genes_from_idx(primary["train"])
    val_genes = genes_from_idx(primary["val"])
    test_genes = genes_from_idx(primary["test"])
    # Assert disjointness already in Phase 3; re-check
    assert set(train_genes).isdisjoint(val_genes)
    assert set(train_genes).isdisjoint(test_genes)

    X_all = feature_matrix(features.set_index("gene_symbol").loc[nodes].reset_index(), FEATURE_COLUMNS)
    feat_by_gene = features.set_index("gene_symbol")

    results: dict[str, object] = {
        "primary_split": {"train": train_genes, "val": val_genes, "test": test_genes},
        "splits_hash": bundle["splits_hash"],
        "pos_degree_mean": bundle["pos_degree_mean"],
        "graph_degree_mean": bundle["graph_degree_mean"],
        "models": {},
    }

    graph_sets = {
        "real": bundle["controls"]["real"],
        "rewired": bundle["controls"]["rewired"],
        "edge_shuffled": bundle["controls"]["edge_shuffled"],
        "sign_shuffled": bundle["controls"]["sign_shuffled"],
        "features_only": None,
    }

    model_names = [
        "random",
        "degree",
        "pagerank",
        "logistic",  # features-only including log_degree (honest graph-vs-stats comparator)
        "hgb",
        "mlp",
        "label_propagation",
        "rwr",
        "rgcn",
        "gatv2",
        "ssl_gatv2",
    ]

    for model in model_names:
        results["models"][model] = {}
        seed_metrics = {regime: {g: [] for g in graph_sets} for regime in ("degree_matched", "pu")}
        print(f"[phase4] model={model}", flush=True)
        for seed in seeds:
            for gname, edges in graph_sets.items():
                if model in {"logistic", "hgb", "mlp", "random"}:
                    # Tabular/random: features_only is the primary setting; also record under real for table completeness
                    if gname not in {"real", "features_only"}:
                        continue
                elif gname == "features_only":
                    continue

                print(f"[phase4]   seed={seed} graph={gname}", flush=True)
                if model == "random":
                    scores = score_random(len(nodes), seed)
                elif model == "degree":
                    scores = score_degree(nodes, edges)
                elif model == "pagerank":
                    scores = score_pagerank(nodes, edges)
                elif model == "rwr":
                    seeds_pos = [g for g in train_genes if y.loc[g] == 1]
                    scores = rwr_scores(nodes, edges, seeds_pos)
                elif model == "label_propagation":
                    scores = label_propagation_scores(nodes, edges, y, train_genes, seed=seed)
                elif model in {"logistic", "hgb", "mlp"}:
                    cols = FEATURE_COLUMNS
                    train_X = feat_by_gene.loc[train_genes, cols].to_numpy(dtype=float)
                    train_y = y.loc[train_genes].to_numpy().astype(int)
                    score_X = feat_by_gene.loc[nodes, cols].to_numpy(dtype=float)
                    scores = fit_predict_tabular(model, train_X, train_y, score_X, seed=seed)
                elif model in {"rgcn", "gatv2"}:
                    scores = train_gnn_scores(
                        model,
                        nodes,
                        features,
                        edges,
                        y,
                        train_genes,
                        val_genes,
                        FEATURE_COLUMNS,
                        seed=seed,
                        epochs=15,
                    )
                elif model == "ssl_gatv2":
                    scores = pretrain_and_finetune(
                        nodes, features, edges, y, train_genes, FEATURE_COLUMNS, seed=seed,
                        pretrain_steps=15, finetune_epochs=15,
                    )
                else:
                    raise ValueError(model)

                for regime in ("degree_matched", "pu"):
                    # Primary metrics on test genes for degree_matched; full for PU ranking
                    if regime == "degree_matched":
                        # Restrict evaluation to test labelled genes
                        score_map = dict(zip(nodes, scores, strict=False))
                        yt = y.loc[test_genes].to_numpy().astype(int)
                        st = np.array([score_map[g] for g in test_genes], dtype=float)
                        m = compute_metrics(yt, st)
                        if seed == seeds[0]:
                            m["auprc_ci"] = bootstrap_ci(yt, st, metric="auprc", seed=seed)
                            m["auroc_ci"] = bootstrap_ci(yt, st, metric="auroc", seed=seed)
                    else:
                        m = evaluate_scores(nodes, y, scores, regime="pu", seed=seed, with_ci=(seed == seeds[0]))
                        if seed != seeds[0]:
                            m.pop("auprc_ci", None)
                            m.pop("auroc_ci", None)
                    seed_metrics[regime][gname].append(m)

        # Aggregate mean±sd across seeds for each graph/regime
        agg = {}
        for regime, gdict in seed_metrics.items():
            agg[regime] = {}
            for gname, lst in gdict.items():
                if not lst:
                    continue
                auprcs = [m["auprc"] for m in lst]
                aurocs = [m["auroc"] for m in lst]
                # Use first seed's CI as representative bootstrap; also report seed spread
                agg[regime][gname] = {
                    "auprc_mean": float(np.nanmean(auprcs)),
                    "auprc_sd": float(np.nanstd(auprcs)),
                    "auroc_mean": float(np.nanmean(aurocs)),
                    "auroc_sd": float(np.nanstd(aurocs)),
                    "auprc_ci_seed0": lst[0].get("auprc_ci"),
                    "n_seeds": len(lst),
                    "per_seed": lst,
                }
        results["models"][model] = agg

    # W1 honesty on primary split, degree_matched, real graph
    rwr = results["models"]["rwr"]["degree_matched"].get("real", {})
    gnn = results["models"]["gatv2"]["degree_matched"].get("real", {})
    rgcn = results["models"]["rgcn"]["degree_matched"].get("real", {})
    # Prefer better of gatv2/rgcn as "GNN"
    def _pick(a, b):
        if not a:
            return b
        if not b:
            return a
        return a if a.get("auprc_mean", -1) >= b.get("auprc_mean", -1) else b

    gnn_best = _pick(gnn, rgcn)
    gnn_name = "gatv2" if gnn_best is gnn else "rgcn"
    w1 = {"gnn_model": gnn_name}
    if rwr and gnn_best and rwr.get("auprc_ci_seed0") and gnn_best.get("auprc_ci_seed0"):
        overlap = ci_overlap(rwr["auprc_ci_seed0"], gnn_best["auprc_ci_seed0"])
        gnn_ge = gnn_best["auprc_mean"] > rwr["auprc_mean"] and not overlap
        # Binding rule: if RWR >= GNN within CI overlap → no GNN advantage
        if rwr["auprc_mean"] >= gnn_best["auprc_mean"] and overlap:
            verdict = "no_gnn_advantage"
            plain = (
                f"W1: RWR AUPRC ({rwr['auprc_mean']:.3f}) ≥ GNN/{gnn_name} AUPRC "
                f"({gnn_best['auprc_mean']:.3f}) with overlapping bootstrap CIs on the primary split — "
                "conclude the GNN demonstrated no advantage."
            )
        elif overlap and abs(rwr["auprc_mean"] - gnn_best["auprc_mean"]) < 1e-6:
            verdict = "no_gnn_advantage"
            plain = "W1: RWR and GNN statistically indistinguishable (CI overlap, equal mean) — no GNN advantage."
        elif rwr["auprc_mean"] >= gnn_best["auprc_mean"] or (overlap and rwr["auprc_ci_seed0"]["lo"] <= gnn_best["auprc_ci_seed0"]["hi"]):
            # CI overlap case when RWR is not lower
            if overlap:
                verdict = "no_gnn_advantage"
                plain = (
                    f"W1: RWR AUPRC CI overlaps GNN/{gnn_name} CI on primary split "
                    f"(RWR {rwr['auprc_mean']:.3f} vs GNN {gnn_best['auprc_mean']:.3f}) — "
                    "conclude no demonstrated GNN advantage under the pre-registered rule."
                )
            else:
                verdict = "rwr_superior"
                plain = f"W1: RWR superior to GNN/{gnn_name} without CI overlap."
        else:
            verdict = "gnn_advantage"
            plain = (
                f"W1: GNN/{gnn_name} AUPRC ({gnn_best['auprc_mean']:.3f}) exceeds RWR "
                f"({rwr['auprc_mean']:.3f}) without CI overlap on the primary split."
            )
        w1.update(
            {
                "verdict": verdict,
                "plain": plain,
                "rwr_auprc_mean": rwr["auprc_mean"],
                "gnn_auprc_mean": gnn_best["auprc_mean"],
                "ci_overlap": overlap,
            }
        )
    else:
        w1.update({"verdict": "incomplete", "plain": "W1 could not be evaluated (missing metrics)."})
    results["h1_w1"] = w1

    # Sign-shuffle conclusion for signed models
    for signed_model in ("rgcn", "gatv2"):
        real_m = results["models"][signed_model]["degree_matched"].get("real", {})
        shuf_m = results["models"][signed_model]["degree_matched"].get("sign_shuffled", {})
        if real_m and shuf_m:
            if abs(real_m["auprc_mean"] - shuf_m["auprc_mean"]) < 0.01 or (
                real_m.get("auprc_ci_seed0")
                and shuf_m.get("auprc_ci_seed0")
                and ci_overlap(real_m["auprc_ci_seed0"], shuf_m["auprc_ci_seed0"])
            ):
                results["models"][signed_model]["sign_shuffle_conclusion"] = (
                    "Signed model performs indistinguishably on sign-shuffled graph — "
                    "not using sign information in a detectable way."
                )
            else:
                results["models"][signed_model]["sign_shuffle_conclusion"] = (
                    "Performance drops under sign-shuffle relative to real graph — "
                    "consistent with use of sign information."
                )

    # H2: SSL vs from-scratch matched budget
    ssl_m = results["models"]["ssl_gatv2"]["degree_matched"].get("real", {})
    scratch = results["models"]["gatv2"]["degree_matched"].get("real", {})
    results["h2"] = {
        "ssl_auprc_mean": ssl_m.get("auprc_mean"),
        "ssl_auprc_sd": ssl_m.get("auprc_sd"),
        "scratch_auprc_mean": scratch.get("auprc_mean"),
        "scratch_auprc_sd": scratch.get("auprc_sd"),
        "n_in_graph_positives": int((y == 1).sum()),
        "note": (
            "With ~246 in-graph positives this comparison is underpowered; "
            "report CI/seed-spread width rather than over-interpreting the point estimate."
        ),
    }

    # Mondrian conformal on primary split using logistic scores (cal=val, test=test)
    train_X = feat_by_gene.loc[train_genes, FEATURE_COLUMNS].to_numpy(dtype=float)
    train_y = y.loc[train_genes].to_numpy().astype(int)
    val_X = feat_by_gene.loc[val_genes, FEATURE_COLUMNS].to_numpy(dtype=float)
    test_X = feat_by_gene.loc[test_genes, FEATURE_COLUMNS].to_numpy(dtype=float)
    # score val+test from model fit on train
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    scaler = StandardScaler().fit(train_X)
    clf = LogisticRegression(max_iter=2000, class_weight="balanced", random_state=0)
    clf.fit(scaler.transform(train_X), train_y)
    cal_scores = clf.predict_proba(scaler.transform(val_X))[:, 1]
    test_scores = clf.predict_proba(scaler.transform(test_X))[:, 1]
    cal_y = y.loc[val_genes].to_numpy().astype(int)
    test_y = y.loc[test_genes].to_numpy().astype(int)
    conf = {}
    for alpha in [plan["conformal_alpha"], *plan["conformal_alpha_sensitivities"]]:
        sets, diag = mondrian_conformal_sets(cal_scores, cal_y, test_scores, alpha=float(alpha))
        cov = empirical_coverage(sets, test_y)
        conf[str(alpha)] = {
            **diag,
            "empirical_coverage": cov,
            "useful": diag["mean_set_size"] < 1.8,
        }
    results["conformal"] = {
        "method": "Mondrian class-conditional",
        "nonclaim": NONCLAIM,
        "by_alpha": conf,
    }

    (Path(output_dir) / "benchmark_results.json").write_text(json.dumps(results, indent=2, default=str) + "\n")
    return results
