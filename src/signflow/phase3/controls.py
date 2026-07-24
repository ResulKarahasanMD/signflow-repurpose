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
    try:
        rewired = nx.directed_edge_swap(g.copy(), nswap=max(1, n_iter_multiplier * g.number_of_edges()), seed=seed)
    except Exception:
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
    """Permute signs among signed edges AND swap stimulation↔inhibition relation_type.

    Leaves unsigned_directed and ambiguous untouched. Preserves degree and the
    sign / signed-relation multiset. Required so R-GCN (edge_type) and GATv2
    (edge_attr sign) both see a changed input (Phase 4.5 A4).
    """

    seed = _seed() if seed is None else seed
    rng = np.random.default_rng(seed)
    out = edges.copy()
    signed_mask = out["sign"].isin([1, -1])
    # Only permute among stimulation/inhibition rows (signed causal)
    stim_inhib = out["relation_type"].isin(["stimulation", "inhibition"])
    mask = signed_mask & stim_inhib

    signs = out.loc[mask, "sign"].to_numpy().copy()
    rels = out.loc[mask, "relation_type"].astype(str).to_numpy().copy()
    before = {
        "plus": int((signs == 1).sum()),
        "minus": int((signs == -1).sum()),
        "stimulation": int((rels == "stimulation").sum()),
        "inhibition": int((rels == "inhibition").sum()),
    }
    # Joint shuffle of (sign, relation) pairs to preserve consistency: +1↔stimulation, -1↔inhibition
    pairs = list(zip(signs.tolist(), rels.tolist(), strict=False))
    rng.shuffle(pairs)
    new_signs = np.array([p[0] for p in pairs], dtype=int)
    new_rels = np.array([p[1] for p in pairs], dtype=object)
    out.loc[mask, "sign"] = new_signs
    out.loc[mask, "relation_type"] = new_rels

    after = {
        "plus": int((out.loc[mask, "sign"] == 1).sum()),
        "minus": int((out.loc[mask, "sign"] == -1).sum()),
        "stimulation": int((out.loc[mask, "relation_type"] == "stimulation").sum()),
        "inhibition": int((out.loc[mask, "relation_type"] == "inhibition").sum()),
    }

    def _deg(frame: pd.DataFrame) -> pd.Series:
        return pd.concat([frame["source"], frame["target"]]).value_counts().sort_index()

    # Fields models read
    edge_attr_changed = not np.array_equal(
        edges.loc[mask, "sign"].to_numpy(), out.loc[mask, "sign"].to_numpy()
    ) if mask.any() else False
    edge_type_changed = not np.array_equal(
        edges.loc[mask, "relation_type"].astype(str).to_numpy(),
        out.loc[mask, "relation_type"].astype(str).to_numpy(),
    ) if mask.any() else False

    props = {
        "n_edges_before": int(len(edges)),
        "n_edges_after": int(len(out)),
        "degree_sequence_equal": bool(_deg(edges).equals(_deg(out))),
        "sign_counts_before": {k: before[k] for k in ("plus", "minus")},
        "sign_counts_after": {k: after[k] for k in ("plus", "minus")},
        "relation_counts_before": {k: before[k] for k in ("stimulation", "inhibition")},
        "relation_counts_after": {k: after[k] for k in ("stimulation", "inhibition")},
        "sign_multiset_preserved": before["plus"] == after["plus"] and before["minus"] == after["minus"],
        "relation_multiset_preserved": before["stimulation"] == after["stimulation"]
        and before["inhibition"] == after["inhibition"],
        "edge_attr_sign_field_changed": bool(edge_attr_changed),
        "edge_type_relation_field_changed": bool(edge_type_changed),
        "permutes_relation_type": True,
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
