"""BENCHMARK conformal / top-k modules and TARGET seed modules."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import yaml
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from signflow.phase3.features import FEATURE_COLUMNS
from signflow.phase4.conformal import mondrian_conformal_sets
from signflow.phase4.data import load_phase3_benchmark

from . import FROZEN_PLAN, SIGNATURE_BENCHMARK, SIGNATURE_TARGET, ensure_directory


def load_plan() -> dict[str, object]:
    return yaml.safe_load(FROZEN_PLAN.read_text())["plan"]


def load_signature_direction(meta_path: str | Path) -> pd.Series:
    meta = pd.read_csv(meta_path, sep="\t")
    meta = meta.dropna(subset=["gene_symbol", "pooled_effect"])
    return meta.set_index("gene_symbol")["pooled_effect"].astype(float)


def _genes_from_fold(labelled_genes: list[str], idx_list: list[int]) -> list[str]:
    return [labelled_genes[i] for i in idx_list]


def build_benchmark_modules(phase3_root: str | Path | None = None, *, plan: dict[str, object] | None = None) -> dict:
    plan = plan or load_plan()
    bundle = load_phase3_benchmark(phase3_root)
    features = bundle["features"]
    y = bundle["y"]
    splits = bundle["splits"]["stratified_cv"]
    nodes = bundle["nodes"]
    labelled_genes = features.loc[features["gene_symbol"].isin(y.index[y >= 0]), "gene_symbol"].astype(str).tolist()
    feat_by_gene = features.set_index("gene_symbol")

    fold0 = splits[0]
    tr = _genes_from_fold(labelled_genes, fold0["train"])
    va = _genes_from_fold(labelled_genes, fold0["val"])

    scaler = StandardScaler().fit(feat_by_gene.loc[tr, FEATURE_COLUMNS])
    clf = LogisticRegression(max_iter=2000, class_weight="balanced", random_state=0)
    clf.fit(scaler.transform(feat_by_gene.loc[tr, FEATURE_COLUMNS]), y.loc[tr].to_numpy().astype(int))
    cal_scores = clf.predict_proba(scaler.transform(feat_by_gene.loc[va, FEATURE_COLUMNS]))[:, 1]
    all_scores = clf.predict_proba(scaler.transform(feat_by_gene.loc[nodes, FEATURE_COLUMNS]))[:, 1]
    score_series = pd.Series(all_scores, index=nodes, name="logistic_score")

    alpha = float(plan["conformal_alpha"])
    _, diag = mondrian_conformal_sets(cal_scores, y.loc[va].to_numpy().astype(int), cal_scores, alpha=alpha)
    thr1 = float(diag["thresholds"].get(1, 1.0))
    conformal_genes = [g for g, s in score_series.items() if (1.0 - float(s)) <= thr1]

    modules = {"conformal_a0.1": conformal_genes}
    for k in plan.get("module_sizes", [50, 100, 250]):
        modules[f"size_{k}"] = score_series.sort_values(ascending=False).head(int(k)).index.astype(str).tolist()

    disease_direction = load_signature_direction(SIGNATURE_BENCHMARK)
    gene_table = pd.DataFrame(
        {
            "gene_symbol": nodes,
            "logistic_score": score_series.values,
            "in_conformal_module": [g in set(conformal_genes) for g in nodes],
            "pooled_effect": disease_direction.reindex(nodes).values,
        }
    )
    return {
        "track": "benchmark",
        "estimator": plan.get("phase5_conformal_module_estimator", "logistic_features_only"),
        "conformal_alpha": alpha,
        "conformal_threshold_class1": thr1,
        "modules": modules,
        "gene_scores": gene_table,
        "disease_direction": disease_direction,
        "splits_hash": bundle["splits_hash"],
    }


def load_target_modules(phase3_root: str | Path | None = None) -> dict:
    from signflow.phase4 import PHASE3_LATEST

    root = Path(phase3_root or PHASE3_LATEST).resolve()
    tdir = root / "target" / "seeds"
    primary = pd.read_csv(tdir / "target_seeds_primary.tsv", sep="\t")
    sens = pd.read_csv(tdir / "target_seeds_sensitivity_concordant.tsv", sep="\t")
    disease_direction = load_signature_direction(SIGNATURE_TARGET)
    return {
        "track": "target",
        "modules": {
            "primary": primary["gene_symbol"].astype(str).tolist(),
            "sensitivity_concordant": sens["gene_symbol"].astype(str).tolist(),
        },
        "seed_frames": {"primary": primary, "sensitivity_concordant": sens},
        "disease_direction": disease_direction,
    }


def write_module_artifacts(module_bundle: dict, output_dir: str | Path) -> dict[str, str]:
    output_dir = ensure_directory(output_dir)
    paths: dict[str, str] = {}
    for name, genes in module_bundle["modules"].items():
        path = Path(output_dir) / f"{name}.tsv"
        pd.DataFrame({"gene_symbol": genes}).to_csv(path, sep="\t", index=False)
        paths[name] = str(path)
    if "gene_scores" in module_bundle:
        gs_path = Path(output_dir) / "gene_scores.tsv"
        module_bundle["gene_scores"].to_csv(gs_path, sep="\t", index=False)
        paths["gene_scores"] = str(gs_path)
    meta_path = Path(output_dir) / "module_meta.json"
    meta = {k: v for k, v in module_bundle.items() if k not in {"gene_scores", "seed_frames", "disease_direction"}}
    meta["n_modules"] = {k: len(v) for k, v in module_bundle["modules"].items()}
    meta_path.write_text(json.dumps(meta, indent=2, default=str) + "\n")
    paths["module_meta"] = str(meta_path)
    return paths
