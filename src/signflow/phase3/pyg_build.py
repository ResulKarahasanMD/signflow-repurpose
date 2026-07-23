"""Assemble PyTorch Geometric Data objects with mappings and GraphML export."""

from __future__ import annotations

import json
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import torch
from torch_geometric.data import Data

from . import ensure_directory
from .features import FEATURE_COLUMNS


def assemble_pyg(
    nodes: list[str],
    features: pd.DataFrame,
    causal_edges: pd.DataFrame,
    string_edges: pd.DataFrame | None,
    labels: pd.Series | None,
    masks: dict[str, np.ndarray] | None,
    id_map: pd.DataFrame,
    output_dir: str | Path,
    *,
    include_string_bidirectional: bool = True,
) -> dict[str, object]:
    output_dir = ensure_directory(output_dir)
    node_index = {g: i for i, g in enumerate(nodes)}
    feat = features.set_index("gene_symbol").loc[nodes, FEATURE_COLUMNS]
    x = torch.tensor(feat.to_numpy(dtype=np.float32), dtype=torch.float32)
    assert torch.isfinite(x).all()

    # Causal directed edges — do NOT bidirect
    src, tgt, attr, etype = [], [], [], []
    relation_types = sorted(causal_edges["relation_type"].unique().tolist())
    rel_to_id = {r: i for i, r in enumerate(relation_types)}
    for _, row in causal_edges.iterrows():
        if row["source"] not in node_index or row["target"] not in node_index:
            continue
        src.append(node_index[row["source"]])
        tgt.append(node_index[row["target"]])
        attr.append([float(row["sign"]), float(row["confidence"])])
        etype.append(rel_to_id[row["relation_type"]])

    string_offset = len(relation_types)
    if string_edges is not None and include_string_bidirectional:
        relation_types.append("string_physical_undirected")
        for _, row in string_edges.iterrows():
            if row["source"] not in node_index or row["target"] not in node_index:
                continue
            i, j = node_index[row["source"]], node_index[row["target"]]
            for a, b in ((i, j), (j, i)):  # bidirectional for undirected auxiliary only
                src.append(a)
                tgt.append(b)
                attr.append([0.0, float(row["confidence"])])
                etype.append(string_offset)

    if not src:
        edge_index = torch.zeros((2, 0), dtype=torch.long)
        edge_attr = torch.zeros((0, 2), dtype=torch.float32)
        edge_type = torch.zeros((0,), dtype=torch.long)
    else:
        edge_index = torch.tensor([src, tgt], dtype=torch.long)
        edge_attr = torch.tensor(attr, dtype=torch.float32)
        edge_type = torch.tensor(etype, dtype=torch.long)
        n = len(nodes)
        assert int(edge_index.min()) >= 0
        assert int(edge_index.max()) < n

    if labels is None:
        y = torch.full((len(nodes),), -1, dtype=torch.long)
    else:
        y = torch.tensor(labels.reindex(nodes).fillna(-1).to_numpy(dtype=np.int64), dtype=torch.long)

    data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr, edge_type=edge_type, y=y)
    if masks:
        for name, arr in masks.items():
            setattr(data, name, torch.tensor(arr, dtype=torch.bool))
            # disjointness checked by caller for train/val/test

    torch.save(data, Path(output_dir) / "pyg_data.pt")
    mapping = {
        "index_to_symbol": {str(i): g for g, i in node_index.items()},
        "symbol_to_index": node_index,
        "relation_types": relation_types,
    }
    # attach uniprot/ensembl
    id_lookup = id_map.drop_duplicates("approved_symbol").set_index("approved_symbol")
    map_rows = []
    for g, i in node_index.items():
        row = {"index": i, "symbol": g}
        if g in id_lookup.index:
            row["ensembl_gene_id"] = id_lookup.loc[g, "ensembl_gene_id"]
            row["uniprot_ids"] = id_lookup.loc[g, "uniprot_ids"]
        map_rows.append(row)
    pd.DataFrame(map_rows).to_csv(Path(output_dir) / "node_index_map.tsv", sep="\t", index=False)
    (Path(output_dir) / "pyg_mapping.json").write_text(json.dumps(mapping, indent=2) + "\n")

    # GraphML of causal layer for inspection
    g = nx.DiGraph()
    for gene in nodes:
        g.add_node(gene)
    for _, row in causal_edges.iterrows():
        if row["source"] in node_index and row["target"] in node_index:
            g.add_edge(row["source"], row["target"], sign=int(row["sign"]), relation=row["relation_type"])
    nx.write_graphml(g, Path(output_dir) / "causal_graph.graphml")
    return {"n_nodes": len(nodes), "n_edges": int(edge_index.shape[1]), "relation_types": relation_types}
