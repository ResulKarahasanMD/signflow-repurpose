#!/usr/bin/env python3
"""Phase 4.6 post-hoc corrections on the corrected Phase 4.5 stamp.

Recomputes pooled-5-fold conformal, hub/LOSO for logistic+rwr+gatv2, and
exploratory paired-difference bootstrap without re-running the full suite.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import yaml
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from signflow.phase3.features import FEATURE_COLUMNS
from signflow.phase4.baselines import fit_predict_tabular, rwr_scores
from signflow.phase4.conformal import NONCLAIM, empirical_coverage, mondrian_conformal_sets
from signflow.phase4.data import labelled_frame, load_phase3_benchmark
from signflow.phase4.metrics import (
    compute_metrics,
    load_bootstrap_n,
    paired_difference_bootstrap,
)
from signflow.phase4.models_gnn import train_gnn_scores

PHASE4_STAMP = Path("data/processed/phase4/20260723T215010Z")
FROZEN = Path("docs/phase0/FROZEN_PLAN.yaml")


def _genes_from_fold(labelled_genes: list[str], idx_list: list[int]) -> list[str]:
    return [labelled_genes[i] for i in idx_list]


def main() -> None:
    plan = yaml.safe_load(FROZEN.read_text())["plan"]
    n_boot = load_bootstrap_n()
    bundle = load_phase3_benchmark(PHASE4_STAMP.parent.parent / "phase3" / "20260723T203525Z")
    # Prefer stamp's phase3 root if recorded
    stamp_meta = json.loads((PHASE4_STAMP / "phase4_result.json").read_text())
    if stamp_meta.get("phase3_root"):
        bundle = load_phase3_benchmark(stamp_meta["phase3_root"])

    nodes = bundle["nodes"]
    y = bundle["y"]
    features = bundle["features"]
    splits = bundle["splits"]
    labelled = labelled_frame(bundle)
    labelled_genes = labelled["gene_symbol"].tolist()
    feat_by_gene = features.set_index("gene_symbol")
    node_index = {g: i for i, g in enumerate(nodes)}
    real_edges = bundle["controls"]["real"]
    stratified = splits["stratified_cv"]
    n_lab = len(labelled_genes)

    results_path = PHASE4_STAMP / "benchmark" / "benchmark_results.json"
    results = json.loads(results_path.read_text())

    # --- pooled 5-fold conformal ---
    conf_by_alpha: dict[str, dict[str, object]] = {}
    alphas = [plan["conformal_alpha"], *plan["conformal_alpha_sensitivities"]]
    for alpha in alphas:
        set_sizes: list[float] = []
        coverages: list[float] = []
        for fold in stratified:
            tr = _genes_from_fold(labelled_genes, fold["train"])
            va = _genes_from_fold(labelled_genes, fold["val"])
            te = _genes_from_fold(labelled_genes, fold["test"])
            scaler = StandardScaler().fit(feat_by_gene.loc[tr, FEATURE_COLUMNS])
            clf = LogisticRegression(max_iter=2000, class_weight="balanced", random_state=0)
            clf.fit(
                scaler.transform(feat_by_gene.loc[tr, FEATURE_COLUMNS]),
                y.loc[tr].to_numpy().astype(int),
            )
            cal_scores = clf.predict_proba(scaler.transform(feat_by_gene.loc[va, FEATURE_COLUMNS]))[:, 1]
            test_scores = clf.predict_proba(scaler.transform(feat_by_gene.loc[te, FEATURE_COLUMNS]))[:, 1]
            sets, diag = mondrian_conformal_sets(
                cal_scores, y.loc[va].to_numpy().astype(int), test_scores, alpha=float(alpha)
            )
            set_sizes.append(float(diag["mean_set_size"]))
            coverages.append(float(empirical_coverage(sets, y.loc[te].to_numpy().astype(int))))
        conf_by_alpha[str(alpha)] = {
            "alpha": float(alpha),
            "mean_set_size": float(np.mean(set_sizes)),
            "empirical_coverage": float(np.mean(coverages)),
            "per_fold_mean_set_size": set_sizes,
            "per_fold_empirical_coverage": coverages,
            "protocol": "pooled_5fold_val_calibration",
        }
    results["conformal"] = {
        "method": "Mondrian class-conditional",
        "protocol": "pooled_5fold_val_calibration",
        "nonclaim": NONCLAIM,
        "by_alpha": conf_by_alpha,
        "fold0_only_superseded": {
            "0.1": {"mean_set_size": 1.808, "empirical_coverage": 0.939},
            "note": "Prior fold-0-only figures retained for audit; superseded by pooled_5fold protocol.",
        },
    }

    def score_model(model: str, train_genes: list[str], val_genes: list[str]):
        if model == "logistic":
            train_X = feat_by_gene.loc[train_genes, FEATURE_COLUMNS].to_numpy(dtype=float)
            train_y = y.loc[train_genes].to_numpy().astype(int)
            score_X = feat_by_gene.loc[nodes, FEATURE_COLUMNS].to_numpy(dtype=float)
            scores = fit_predict_tabular("logistic", train_X, train_y, score_X, seed=0)
            return scores, scores
        if model == "rwr":
            seeds_pos = [g for g in train_genes if y.loc[g] == 1]
            scores = rwr_scores(nodes, real_edges, seeds_pos)
            return scores, None
        if model == "gatv2":
            scores, _ = train_gnn_scores(
                "gatv2", nodes, features, real_edges, y, train_genes, val_genes, FEATURE_COLUMNS, seed=0
            )
            return scores, scores
        raise ValueError(model)

    # --- hub / LOSO for logistic, rwr, gatv2 ---
    extra: dict[str, object] = {}
    for split_name, fold_list in (
        ("hub_stratified_cv", splits.get("hub_stratified_cv", [])),
        ("leave_one_study_out", splits.get("leave_one_study_out", [])),
    ):
        by_model: dict[str, object] = {}
        for model in ("logistic", "rwr", "gatv2"):
            print(f"[phase4.6] {split_name} model={model}", flush=True)
            fold_auprcs = []
            for fold in fold_list:
                if split_name == "hub_stratified_cv":
                    tr = _genes_from_fold(labelled_genes, fold["train"])
                    te = _genes_from_fold(labelled_genes, fold["test"])
                    va = _genes_from_fold(labelled_genes, fold["val"])
                else:
                    tr = [nodes[i] for i in fold["train"] if y.loc[nodes[i]] >= 0]
                    te = [nodes[i] for i in fold["test"] if y.loc[nodes[i]] >= 0]
                    va = [nodes[i] for i in fold["val"] if y.loc[nodes[i]] >= 0]
                    if len(te) < 5 or len(tr) < 10:
                        continue
                scores_nodes, probs = score_model(model, tr, va if va else tr[: max(1, len(tr) // 5)])
                st = np.array([scores_nodes[node_index[g]] for g in te])
                yt = y.loc[te].to_numpy().astype(int)
                if len(np.unique(yt)) < 2:
                    continue
                fold_auprcs.append(
                    compute_metrics(
                        yt,
                        st,
                        probabilities=st if probs is not None else None,
                        regime="degree_matched",
                    )["auprc"]
                )
            by_model[model] = {
                "model": model,
                "n_folds_evaluated": len(fold_auprcs),
                "auprc_mean": float(np.mean(fold_auprcs)) if fold_auprcs else float("nan"),
                "auprc_sd": float(np.std(fold_auprcs)) if fold_auprcs else float("nan"),
                "per_fold_auprc": fold_auprcs,
            }
        extra[split_name] = by_model
    results["additional_splits"] = extra

    # --- exploratory paired bootstrap (recompute seed-0 OOF for rwr/gatv2) ---
    print("[phase4.6] paired-difference OOF for rwr/gatv2", flush=True)
    oof_rwr = np.full(n_lab, np.nan)
    oof_gat = np.full(n_lab, np.nan)
    oof_y = np.array([int(y.loc[g]) for g in labelled_genes], dtype=int)
    for fold in stratified:
        tr = _genes_from_fold(labelled_genes, fold["train"])
        va = _genes_from_fold(labelled_genes, fold["val"])
        te = _genes_from_fold(labelled_genes, fold["test"])
        s_rwr, _ = score_model("rwr", tr, va)
        s_gat, _ = score_model("gatv2", tr, va)
        for g in te:
            li = labelled_genes.index(g)
            oof_rwr[li] = s_rwr[node_index[g]]
            oof_gat[li] = s_gat[node_index[g]]
    paired = paired_difference_bootstrap(oof_y, oof_gat, oof_rwr, metric="auprc", seed=0, n_boot=n_boot)
    paired["contrast"] = "gatv2_minus_rwr_seed0_pooled_OOF"
    paired["w1_binding"] = True
    results.setdefault("h1_w1", {})
    results["h1_w1"]["paired_difference_bootstrap_exploratory"] = paired

    # Chance / sign narrative fields for consumers
    results["h1_w1"]["reframing"] = (
        "Chance AUPRC on balanced labels = 0.5. logistic 0.615 [0.551, 0.680], "
        "gatv2 0.590 [0.529, 0.658], rgcn 0.580 [0.516, 0.646], ssl_gatv2 0.579 [0.517, 0.647], "
        "hgb and mlp all have CI lower bounds above 0.5; rwr 0.524 [0.466, 0.594], "
        "label_propagation, pagerank, degree and random all include 0.5. "
        "RWR — the pre-registered structural baseline — is itself at chance."
    )
    for signed_model, real_v, shuf_v, rew_v in (
        ("gatv2", 0.590, 0.588, 0.516),
        ("ssl_gatv2", 0.579, 0.589, 0.503),
        ("rgcn", 0.580, 0.586, 0.556),
    ):
        results["models"][signed_model]["sign_shuffle_conclusion"] = (
            f"real≈sign_shuffled ({real_v:.3f} / {shuf_v:.3f}) but real>rewired ({rew_v:.3f}): "
            "topology is used; sign is not detectably used."
            + (
                " PRIOR R-GCN identical-to-shuffle conclusion RETRACTED (inert control artifact)."
                if signed_model == "rgcn"
                else ""
            )
        )

    results_path.write_text(json.dumps(results, indent=2, default=str) + "\n")
    out = {
        "conformal": results["conformal"],
        "additional_splits": results["additional_splits"],
        "paired_difference_bootstrap_exploratory": paired,
    }
    (PHASE4_STAMP / "benchmark" / "phase46_posthoc.json").write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(json.dumps(out, indent=2, default=str))


if __name__ == "__main__":
    main()
