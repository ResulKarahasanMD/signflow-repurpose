"""Load Phase 3 artifacts; assert split hashes; never re-split."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import yaml

from . import FROZEN_PLAN, PHASE3_LATEST, StopAndReportError, sha256_path


def hash_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def load_frozen_seeds() -> list[int]:
    return list(yaml.safe_load(FROZEN_PLAN.read_text())["plan"]["random_seeds"])


def load_phase3_benchmark(phase3_root: Path | None = None) -> dict[str, object]:
    root = Path(phase3_root or PHASE3_LATEST)
    if root.is_symlink():
        root = root.resolve()
    bench = root / "benchmark"
    if not bench.exists():
        raise StopAndReportError("Phase 3 benchmark artifacts missing", reason="phase3-missing")

    features = pd.read_csv(bench / "node_features.tsv", sep="\t")
    edges = pd.read_csv(bench / "graph_edges.tsv", sep="\t")
    nodes = features["gene_symbol"].astype(str).tolist()
    splits_path = bench / "splits" / "splits.json"
    splits_raw = splits_path.read_bytes()
    splits_hash = hash_bytes(splits_raw)
    splits = json.loads(splits_raw.decode())

    pos = pd.read_csv(bench / "labels" / "uc_positives_in_graph.tsv", sep="\t")
    neg = pd.read_csv(bench / "labels" / "degree_matched_controls.tsv", sep="\t")
    y = pd.Series(-1, index=nodes, dtype=int)
    y.loc[y.index.isin(pos["gene_symbol"])] = 1
    y.loc[y.index.isin(neg["gene_symbol"])] = 0

    # Control graphs
    controls = {
        "real": edges,
        "rewired": pd.read_csv(bench / "controls" / "edges_rewired.tsv", sep="\t"),
        "edge_shuffled": pd.read_csv(bench / "controls" / "edges_edge_shuffled.tsv", sep="\t"),
        "sign_shuffled": pd.read_csv(bench / "controls" / "edges_sign_shuffled.tsv", sep="\t"),
    }

    density = pd.read_csv(bench / "feature_density_qc.tsv", sep="\t")
    # Fail if near-constant continuous features survived
    bad = density.loc[density["near_constant_flag"] & ~density["feature"].isin(
        ["is_tf", "progeny_covered", "collectri_covered", "regulon_target_covered"]
    )]
    if len(bad):
        raise StopAndReportError(
            f"Near-constant features in Phase 3 artifacts: {bad['feature'].tolist()}",
            reason="feature-near-constant",
        )

    # Assert unique pairs
    if edges.duplicated(subset=["source", "target"]).any():
        raise StopAndReportError("Duplicate causal pairs in Phase 3 edges", reason="duplicate-causal-pairs")

    pyg_path = bench / "pyg" / "pyg_data.pt"
    data = torch.load(pyg_path, weights_only=False) if pyg_path.exists() else None

    return {
        "root": str(root),
        "benchmark_dir": str(bench),
        "nodes": nodes,
        "features": features,
        "edges": edges,
        "controls": controls,
        "y": y,
        "splits": splits,
        "splits_path": str(splits_path),
        "splits_hash": splits_hash,
        "density": density,
        "pyg": data,
        "pos_degree_mean": float(pd.read_csv(bench / "labels" / "degree_match_qc.tsv", sep="\t").iloc[0]["pos_degree_mean"]),
        "graph_degree_mean": float(json.loads((bench / "network_summary.json").read_text())["degree_mean"]),
    }


def assert_splits_hash(bundle: dict[str, object], expected_hash: str | None = None) -> str:
    path = Path(bundle["splits_path"])
    observed = hash_bytes(path.read_bytes())
    if observed != bundle["splits_hash"]:
        raise StopAndReportError("Split hash mismatch vs loaded bundle", reason="split-hash-mismatch")
    if expected_hash is not None and observed != expected_hash:
        raise StopAndReportError(
            f"Split hash mismatch vs recorded Phase 3 hash: {observed} != {expected_hash}",
            reason="split-hash-mismatch",
        )
    return observed


def labelled_frame(bundle: dict[str, object]) -> pd.DataFrame:
    """Rows with y in {{0,1}} in labelled-subset order matching stratified split indices."""

    y = bundle["y"]
    labelled = y[y >= 0]
    feats = bundle["features"].set_index("gene_symbol").loc[labelled.index]
    out = feats.copy()
    out["label"] = labelled.to_numpy()
    out["gene_symbol"] = labelled.index
    return out.reset_index(drop=True)
