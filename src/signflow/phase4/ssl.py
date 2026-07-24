"""SSL pretraining: edge-sign prediction, then few-shot fine-tune (matched converged budget)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
import yaml
from sklearn.metrics import average_precision_score

from . import FROZEN_PLAN
from .models_gnn import SignedGATv2, build_pyg_from_edges, load_gnn_hparams


def pretrain_and_finetune(
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
    plan = yaml.safe_load(FROZEN_PLAN.read_text())["plan"]
    hp = load_gnn_hparams()
    pretrain_steps = int(plan["phase4_ssl_pretrain_steps"])
    # Matched converged budget with scratch
    max_epochs = int(plan["phase4_ssl_finetune_max_epochs"])
    min_epochs = hp["min_epochs"]
    patience = hp["patience"]

    torch.manual_seed(seed)
    np.random.seed(seed)
    index = {g: i for i, g in enumerate(nodes)}
    train_idx = np.array([index[g] for g in train_genes if g in index and y.loc[g] >= 0], dtype=int)
    val_idx = np.array([index[g] for g in val_genes if g in index and y.loc[g] >= 0], dtype=int)
    diag: dict[str, Any] = {
        "model": "ssl_gatv2",
        "seed": seed,
        "pretrain_steps": pretrain_steps,
        "max_epochs": max_epochs,
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
    encoder = SignedGATv2(data.x.size(1), hp["hidden"], dropout=hp["dropout"])
    encoder.lin = nn.Identity()
    sign_head = nn.Linear(2 * hp["hidden"], 3)
    opt = torch.optim.Adam(
        list(encoder.parameters()) + list(sign_head.parameters()),
        lr=hp["lr"],
        weight_decay=hp["weight_decay"],
    )
    if data.edge_index.numel():
        for _ in range(pretrain_steps):
            encoder.train()
            opt.zero_grad()
            h = encoder(data.x, data.edge_index, data.edge_attr)
            if h.ndim == 1:
                break
            src, tgt = data.edge_index
            logits = sign_head(torch.cat([h[src], h[tgt]], dim=-1))
            signs = (data.edge_attr[:, 0].long() + 1).clamp(0, 2)
            loss = F.cross_entropy(logits, signs)
            loss.backward()
            opt.step()

    encoder.lin = nn.Linear(hp["hidden"], 1)
    opt2 = torch.optim.Adam(encoder.parameters(), lr=hp["lr"], weight_decay=hp["weight_decay"])
    y_full = torch.tensor(y.reindex(nodes).fillna(-1).to_numpy(), dtype=torch.float32)
    best_state = None
    best_val = -1.0
    best_epoch = 0
    stall = 0
    for epoch in range(1, max_epochs + 1):
        encoder.train()
        opt2.zero_grad()
        logits = encoder(data.x, data.edge_index, data.edge_attr).squeeze(-1)
        loss = F.binary_cross_entropy_with_logits(logits[train_idx], y_full[train_idx])
        loss.backward()
        opt2.step()
        diag["train_loss_curve"].append(float(loss.detach().cpu()))
        encoder.eval()
        with torch.no_grad():
            logits = encoder(data.x, data.edge_index, data.edge_attr).squeeze(-1)
            ap = float("nan")
            if len(val_idx) and len(np.unique(y_full[val_idx].numpy())) > 1:
                ap = float(average_precision_score(y_full[val_idx].numpy(), torch.sigmoid(logits[val_idx]).numpy()))
            diag["val_auprc_curve"].append(ap)
            if np.isfinite(ap) and ap > best_val + 1e-6:
                best_val = ap
                best_epoch = epoch
                best_state = {k: v.detach().cpu().clone() for k, v in encoder.state_dict().items()}
                stall = 0
            else:
                stall += 1
        if epoch >= min_epochs and stall >= patience:
            diag["early_stopped"] = True
            break
    else:
        diag["hit_epoch_cap"] = True
    diag["best_epoch"] = best_epoch
    diag["best_val_auprc"] = best_val
    if diag["hit_epoch_cap"] and stall < patience:
        diag["undertrained_no_plateau"] = True
    if best_state is not None:
        encoder.load_state_dict(best_state)
    encoder.eval()
    with torch.no_grad():
        logits = encoder(data.x, data.edge_index, data.edge_attr).squeeze(-1)
        return torch.sigmoid(logits).cpu().numpy(), diag
