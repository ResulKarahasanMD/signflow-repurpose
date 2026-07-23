"""TARGET unsupervised propagation from frozen seed sets."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from . import PHASE3_LATEST, ensure_directory
from .baselines import rwr_scores


def run_target_propagation(output_dir: str | Path, phase3_root: Path | None = None) -> dict[str, object]:
    root = Path(phase3_root or PHASE3_LATEST).resolve()
    tdir = root / "target"
    output_dir = ensure_directory(output_dir)
    features = pd.read_csv(tdir / "node_features.tsv", sep="\t")
    nodes = features["gene_symbol"].astype(str).tolist()
    edges = pd.read_csv(tdir / "graph_edges.tsv", sep="\t")
    controls = {
        "real": edges,
        "rewired": pd.read_csv(tdir / "controls" / "edges_rewired.tsv", sep="\t"),
        "edge_shuffled": pd.read_csv(tdir / "controls" / "edges_edge_shuffled.tsv", sep="\t"),
        "sign_shuffled": pd.read_csv(tdir / "controls" / "edges_sign_shuffled.tsv", sep="\t"),
    }
    primary = pd.read_csv(tdir / "seeds" / "target_seeds_primary.tsv", sep="\t")
    sens = pd.read_csv(tdir / "seeds" / "target_seeds_sensitivity_concordant.tsv", sep="\t")
    results = {}
    for seed_name, seed_df in (("primary", primary), ("sensitivity", sens)):
        seed_genes = seed_df["gene_symbol"].astype(str).tolist()
        rankings = {}
        for gname, gedges in controls.items():
            scores = rwr_scores(nodes, gedges, seed_genes)
            rank = pd.DataFrame({"gene_symbol": nodes, "score": scores}).sort_values("score", ascending=False)
            # Carry seed replication annotations for seed genes
            rank = rank.merge(seed_df, on="gene_symbol", how="left", suffixes=("", "_seed"))
            path = Path(output_dir) / f"target_rank_{seed_name}_{gname}.tsv"
            rank.to_csv(path, sep="\t", index=False)
            rankings[gname] = str(path)
        results[seed_name] = {
            "n_seeds": len(seed_genes),
            "rankings": rankings,
            "top50": pd.read_csv(rankings["real"], sep="\t").head(50)["gene_symbol"].tolist(),
        }
    overlap = set(results["primary"]["top50"]) & set(results["sensitivity"]["top50"])
    summary = {
        "primary_top50_n": 50,
        "sensitivity_top50_n": 50,
        "top50_overlap": sorted(overlap),
        "n_overlap": len(overlap),
        "note": "Sign-shuffle/rewire rankings are the primary topology vs biology controls for label-less TARGET.",
        "paths": results,
    }
    (Path(output_dir) / "target_propagation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary
