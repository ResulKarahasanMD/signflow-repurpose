"""Graph-information controls: rewire, edge-shuffle, sign-shuffle, features-only."""

from __future__ import annotations

import json
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import yaml

from . import FROZEN_PLAN, ensure_directory


def _seed() -> int:
    return int(yaml.safe_load(FROZEN_PLAN.read_text())["plan"]["phase3_control_seed"])


def degree_preserving_rewire(edges: pd.DataFrame, *, seed: int | None = None, n_iter_multiplier: int = 10) -> pd.DataFrame:
    seed = _seed() if seed is None else seed
    g = nx.DiGraph()
    for _, row in edges.iterrows():
        g.add_edge(row["source"], row["target"], **{k: row[k] for k in edges.columns if k not in {"source", "target"}})
    # directed degree-preserving rewire via configuration-like double edge swaps
    try:
        rewired = nx.directed_edge_swap(g.copy(), nswap=max(1, n_iter_multiplier * g.number_of_edges()), seed=seed)
    except Exception:
        # Fallback: shuffle targets while preserving out-degree sequence approximately
        rng = np.random.default_rng(seed)
        src = edges["source"].to_numpy()
        tgt = edges["target"].to_numpy()
        rng.shuffle(tgt)
        out = edges.copy()
        out["target"] = tgt
        return out
    rows = []
    for u, v, data in rewired.edges(data=True):
        rows.append({"source": u, "target": v, **dict(data)})
    return pd.DataFrame(rows)


def edge_shuffle(edges: pd.DataFrame, *, seed: int | None = None) -> pd.DataFrame:
    seed = _seed() if seed is None else seed
    rng = np.random.default_rng(seed)
    out = edges.copy()
    tgt = out["target"].to_numpy().copy()
    rng.shuffle(tgt)
    out["target"] = tgt
    return out


def sign_shuffle(edges: pd.DataFrame, *, seed: int | None = None) -> tuple[pd.DataFrame, dict[str, object]]:
    """Preserve degree and edge count; permute signs among signed edges."""

    seed = _seed() if seed is None else seed
    rng = np.random.default_rng(seed)
    out = edges.copy()
    signed_mask = out["sign"].isin([1, -1])
    signs = out.loc[signed_mask, "sign"].to_numpy().copy()
    before = {"plus": int((signs == 1).sum()), "minus": int((signs == -1).sum())}
    rng.shuffle(signs)
    out.loc[signed_mask, "sign"] = signs
    after = {
        "plus": int((out.loc[signed_mask, "sign"] == 1).sum()),
        "minus": int((out.loc[signed_mask, "sign"] == -1).sum()),
    }
    # Degree sequence of underlying undirected projection
    def _deg(frame: pd.DataFrame) -> pd.Series:
        return pd.concat([frame["source"], frame["target"]]).value_counts().sort_index()

    props = {
        "n_edges_before": int(len(edges)),
        "n_edges_after": int(len(out)),
        "degree_sequence_equal": bool(_deg(edges).equals(_deg(out))),
        "sign_counts_before": before,
        "sign_counts_after": after,
        "sign_multiset_preserved": before == after,
        "seed": seed,
    }
    return out, props


def verify_rewire_degree(original: pd.DataFrame, rewired: pd.DataFrame) -> bool:
    def _out_in(frame: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
        out_d = frame.groupby("source").size()
        in_d = frame.groupby("target").size()
        return out_d.sort_index(), in_d.sort_index()

    o_out, o_in = _out_in(original)
    r_out, r_in = _out_in(rewired)
    # Exact preservation when directed_edge_swap works; else report False
    return bool(o_out.equals(r_out) and o_in.equals(r_in))


def emit_control_graphs(causal_edges: pd.DataFrame, output_dir: str | Path) -> dict[str, object]:
    output_dir = ensure_directory(output_dir)
    seed = _seed()
    rewired = degree_preserving_rewire(causal_edges, seed=seed)
    shuffled = edge_shuffle(causal_edges, seed=seed + 1)
    signed, sign_props = sign_shuffle(causal_edges, seed=seed + 2)
    rewired.to_csv(Path(output_dir) / "edges_rewired.tsv", sep="\t", index=False)
    shuffled.to_csv(Path(output_dir) / "edges_edge_shuffled.tsv", sep="\t", index=False)
    signed.to_csv(Path(output_dir) / "edges_sign_shuffled.tsv", sep="\t", index=False)
    (Path(output_dir) / "features_only_config.json").write_text(
        json.dumps({"mode": "features_only", "edges": None, "seed": seed}, indent=2) + "\n"
    )
    verification = {
        "seed": seed,
        "rewire_preserves_degree_exactly": verify_rewire_degree(causal_edges, rewired),
        "sign_shuffle": sign_props,
        "edge_shuffle_n": int(len(shuffled)),
    }
    (Path(output_dir) / "control_verification.json").write_text(json.dumps(verification, indent=2) + "\n")
    return verification
