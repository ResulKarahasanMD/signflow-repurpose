"""Relational/signed message-passing models (R-GCN-style and GATv2)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
import yaml
from sklearn.metrics import average_precision_score
from sklearn.preprocessing import StandardScaler
from torch_geometric.data import Data
from torch_geometric.nn import GATv2Conv, RGCNConv

from . import FROZEN_PLAN


_GNN_HPARAMS: dict[str, Any] | None = None


def load_gnn_hparams() -> dict[str, Any]:
    global _GNN_HPARAMS
    if _GNN_HPARAMS is not None:
        return _GNN_HPARAMS
    plan = yaml.safe_load(FROZEN_PLAN.read_text())["plan"]
    _GNN_HPARAMS = {
        "max_epochs": int(plan["phase4_gnn_max_epochs"]),
        "min_epochs": int(plan["phase4_gnn_min_epochs"]),
        "patience": int(plan["phase4_gnn_early_stop_patience"]),
        "lr": float(plan["phase4_gnn_lr"]),
        "weight_decay": float(plan["phase4_gnn_weight_decay"]),
        "hidden": int(plan["phase4_gnn_hidden"]),
        "dropout": float(plan["phase4_gnn_dropout"]),
    }
    return _GNN_HPARAMS


class SignedRGCN(nn.Module):
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
        h = self.conv2(h, edge_index, edge_type) + h
        h = F.relu(h)
        return self.lin(h).squeeze(-1)


class SignedGATv2(nn.Module):
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

    # Fixed global relation vocabulary so edge_type IDs are stable across shuffles
    rels = ["ambiguous", "inhibition", "stimulation", "unsigned_directed"]
    present = sorted(edges["relation_type"].astype(str).unique().tolist())
    for r in present:
        if r not in rels:
            rels.append(r)
    rel_id = {r: i for i, r in enumerate(rels)}
    src_s = edges["source"].astype(str).to_numpy()
    tgt_s = edges["target"].astype(str).to_numpy()
    keep = np.fromiter((s in index and t in index for s, t in zip(src_s, tgt_s)), dtype=bool, count=len(src_s))
    if keep.any():
        src = [index[s] for s in src_s[keep]]
        tgt = [index[t] for t in tgt_s[keep]]
        etype = [rel_id[str(r)] for r in edges.loc[keep, "relation_type"].astype(str)]
        signs = edges.loc[keep, "sign"].astype(float).to_numpy()
        conf = (
            edges.loc[keep, "confidence"].astype(float).to_numpy()
            if "confidence" in edges.columns
            else np.ones(int(keep.sum()))
        )
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
        relation_types=rels,
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
) -> tuple[np.ndarray, dict[str, Any]]:
    """Train to frozen protocol; return scores + convergence diagnostic."""

    hp = load_gnn_hparams()
    torch.manual_seed(seed)
    np.random.seed(seed)
    index = {g: i for i, g in enumerate(nodes)}
    train_idx = np.array([index[g] for g in train_genes if g in index and y.loc[g] >= 0], dtype=int)
    val_idx = np.array([index[g] for g in val_genes if g in index and y.loc[g] >= 0], dtype=int)
    diag: dict[str, Any] = {
        "model": model_name,
        "seed": seed,
        "max_epochs": hp["max_epochs"],
        "min_epochs": hp["min_epochs"],
        "patience": hp["patience"],
        "train_loss_curve": [],
        "val_auprc_curve": [],
        "best_epoch": None,
        "early_stopped": False,
        "hit_epoch_cap": False,
        "undertrained_no_plateau": False,
    }
    if len(train_idx) < 4:
        return np.zeros(len(nodes)), {**diag, "status": "insufficient_train"}

    data, _ = build_pyg_from_edges(
        nodes, features, edges, feature_cols, fit_scaler=True, train_idx=train_idx
    )
    y_full = torch.tensor(y.reindex(nodes).fillna(-1).to_numpy(), dtype=torch.float32)

    if model_name == "rgcn":
        model = SignedRGCN(data.x.size(1), hp["hidden"], int(data.n_relations), dropout=hp["dropout"])
    elif model_name == "gatv2":
        model = SignedGATv2(data.x.size(1), hp["hidden"], dropout=hp["dropout"])
    else:
        raise ValueError(model_name)

    opt = torch.optim.Adam(model.parameters(), lr=hp["lr"], weight_decay=hp["weight_decay"])
    best_state = None
    best_val = -1.0
    best_epoch = 0
    stall = 0

    for epoch in range(1, hp["max_epochs"] + 1):
        model.train()
        opt.zero_grad()
        if model_name == "rgcn":
            logits = model(data.x, data.edge_index, data.edge_type)
        else:
            logits = model(data.x, data.edge_index, data.edge_attr)
        loss = F.binary_cross_entropy_with_logits(logits[train_idx], y_full[train_idx])
        loss.backward()
        opt.step()
        diag["train_loss_curve"].append(float(loss.detach().cpu()))

        model.eval()
        with torch.no_grad():
            if model_name == "rgcn":
                logits = model(data.x, data.edge_index, data.edge_type)
            else:
                logits = model(data.x, data.edge_index, data.edge_attr)
            ap = float("nan")
            if len(val_idx) and len(np.unique(y_full[val_idx].numpy())) > 1:
                ap = float(average_precision_score(y_full[val_idx].numpy(), torch.sigmoid(logits[val_idx]).numpy()))
            diag["val_auprc_curve"].append(ap)
            improved = np.isfinite(ap) and ap > best_val + 1e-6
            if improved:
                best_val = ap
                best_epoch = epoch
                best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
                stall = 0
            else:
                stall += 1

        if epoch >= hp["min_epochs"] and stall >= hp["patience"]:
            diag["early_stopped"] = True
            break
    else:
        diag["hit_epoch_cap"] = True

    diag["best_epoch"] = best_epoch
    diag["best_val_auprc"] = best_val
    # Cap hit without recent improvement near the end → still undertrained signal
    if diag["hit_epoch_cap"] and stall < hp["patience"]:
        diag["undertrained_no_plateau"] = True

    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        if model_name == "rgcn":
            logits = model(data.x, data.edge_index, data.edge_type)
        else:
            logits = model(data.x, data.edge_index, data.edge_attr)
        scores = torch.sigmoid(logits).cpu().numpy()
    return scores, diag
