"""Relational/signed message-passing models (R-GCN-style and GATv2)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.preprocessing import StandardScaler
from torch_geometric.data import Data
from torch_geometric.nn import GATv2Conv, RGCNConv


class SignedRGCN(nn.Module):
    """2-layer R-GCN on relation types; edge signs concatenated into node init via mean neighbor sign (simple)."""

    def __init__(self, in_dim: int, hidden: int, n_relations: int, dropout: float = 0.2):
        super().__init__()
        self.conv1 = RGCNConv(in_dim, hidden, num_relations=max(n_relations, 1))
        self.conv2 = RGCNConv(hidden, hidden, num_relations=max(n_relations, 1))
        self.dropout = dropout
        self.lin = nn.Linear(hidden, 1)

    def forward(self, x, edge_index, edge_type):
        h = self.conv1(x, edge_index, edge_type)
        h = F.relu(h)
        h = F.dropout(h, p=self.dropout, training=self.training)
        h = self.conv2(h, edge_index, edge_type) + h  # residual
        h = F.relu(h)
        return self.lin(h).squeeze(-1)


class SignedGATv2(nn.Module):
    """GATv2 with edge attributes [sign, confidence]."""

    def __init__(self, in_dim: int, hidden: int, dropout: float = 0.2):
        super().__init__()
        self.conv1 = GATv2Conv(in_dim, hidden, edge_dim=2, heads=2, concat=False, dropout=dropout)
        self.conv2 = GATv2Conv(hidden, hidden, edge_dim=2, heads=2, concat=False, dropout=dropout)
        self.dropout = dropout
        self.lin = nn.Linear(hidden, 1)

    def forward(self, x, edge_index, edge_attr):
        h = self.conv1(x, edge_index, edge_attr)
        h = F.relu(h)
        h = F.dropout(h, p=self.dropout, training=self.training)
        h = self.conv2(h, edge_index, edge_attr) + h
        h = F.relu(h)
        return self.lin(h).squeeze(-1)


def build_pyg_from_edges(
    nodes: list[str],
    features: pd.DataFrame,
    edges: pd.DataFrame,
    feature_cols: list[str],
    *,
    scaler: StandardScaler | None = None,
    fit_scaler: bool = False,
    train_idx: np.ndarray | None = None,
) -> tuple[Data, StandardScaler]:
    index = {g: i for i, g in enumerate(nodes)}
    X = features.set_index("gene_symbol").loc[nodes, feature_cols].to_numpy(dtype=np.float32)
    if scaler is None:
        scaler = StandardScaler()
    if fit_scaler:
        if train_idx is None:
            raise ValueError("train_idx required to fit scaler")
        scaler.fit(X[train_idx])
    X = scaler.transform(X).astype(np.float32)

    rels = sorted(edges["relation_type"].astype(str).unique().tolist())
    rel_id = {r: i for i, r in enumerate(rels)}
    src_s = edges["source"].astype(str).to_numpy()
    tgt_s = edges["target"].astype(str).to_numpy()
    keep = np.array([s in index and t in index for s, t in zip(src_s, tgt_s)], dtype=bool)
    if keep.any():
        src = [index[s] for s in src_s[keep]]
        tgt = [index[t] for t in tgt_s[keep]]
        etype = [rel_id[str(r)] for r in edges.loc[keep, "relation_type"].astype(str)]
        signs = edges.loc[keep, "sign"].astype(float).to_numpy()
        conf = edges.loc[keep, "confidence"].astype(float).to_numpy() if "confidence" in edges.columns else np.ones(keep.sum())
        edge_index = torch.tensor([src, tgt], dtype=torch.long)
        edge_type = torch.tensor(etype, dtype=torch.long)
        edge_attr = torch.tensor(np.column_stack([signs, conf]), dtype=torch.float32)
    else:
        edge_index = torch.zeros((2, 0), dtype=torch.long)
        edge_type = torch.zeros((0,), dtype=torch.long)
        edge_attr = torch.zeros((0, 2), dtype=torch.float32)
    data = Data(
        x=torch.tensor(X, dtype=torch.float32),
        edge_index=edge_index,
        edge_type=edge_type,
        edge_attr=edge_attr,
        n_relations=len(rels),
    )
    return data, scaler


def train_gnn_scores(
    model_name: str,
    nodes: list[str],
    features: pd.DataFrame,
    edges: pd.DataFrame,
    y: pd.Series,
    train_genes: list[str],
    val_genes: list[str],
    feature_cols: list[str],
    *,
    seed: int,
    epochs: int = 40,
    hidden: int = 32,
) -> np.ndarray:
    torch.manual_seed(seed)
    np.random.seed(seed)
    index = {g: i for i, g in enumerate(nodes)}
    train_idx = np.array([index[g] for g in train_genes if g in index and y.loc[g] >= 0], dtype=int)
    val_idx = np.array([index[g] for g in val_genes if g in index and y.loc[g] >= 0], dtype=int)
    if len(train_idx) < 4:
        return np.zeros(len(nodes))

    data, _ = build_pyg_from_edges(
        nodes, features, edges, feature_cols, fit_scaler=True, train_idx=train_idx
    )
    y_full = torch.tensor(y.reindex(nodes).fillna(-1).to_numpy(), dtype=torch.float32)

    if model_name == "rgcn":
        model = SignedRGCN(data.x.size(1), hidden, int(data.n_relations))
    elif model_name == "gatv2":
        model = SignedGATv2(data.x.size(1), hidden)
    else:
        raise ValueError(model_name)

    opt = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-4)
    best_state = None
    best_val = -1.0
    for _ in range(epochs):
        model.train()
        opt.zero_grad()
        if model_name == "rgcn":
            logits = model(data.x, data.edge_index, data.edge_type)
        else:
            logits = model(data.x, data.edge_index, data.edge_attr)
        loss = F.binary_cross_entropy_with_logits(logits[train_idx], y_full[train_idx])
        loss.backward()
        opt.step()
        model.eval()
        with torch.no_grad():
            if model_name == "rgcn":
                logits = model(data.x, data.edge_index, data.edge_type)
            else:
                logits = model(data.x, data.edge_index, data.edge_attr)
            if len(val_idx):
                val_prob = torch.sigmoid(logits[val_idx]).numpy()
                val_y = y_full[val_idx].numpy()
                if len(np.unique(val_y)) > 1:
                    from sklearn.metrics import average_precision_score

                    ap = average_precision_score(val_y, val_prob)
                    if ap > best_val:
                        best_val = ap
                        best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        if model_name == "rgcn":
            logits = model(data.x, data.edge_index, data.edge_type)
        else:
            logits = model(data.x, data.edge_index, data.edge_attr)
        scores = torch.sigmoid(logits).cpu().numpy()
    return scores
