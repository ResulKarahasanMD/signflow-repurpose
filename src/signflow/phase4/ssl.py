"""SSL pretraining: edge-sign prediction, then few-shot fine-tune (matched-budget vs scratch)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F

from .models_gnn import SignedGATv2, build_pyg_from_edges


def pretrain_and_finetune(
    nodes: list[str],
    features: pd.DataFrame,
    edges: pd.DataFrame,
    y: pd.Series,
    train_genes: list[str],
    feature_cols: list[str],
    *,
    seed: int,
    pretrain_steps: int = 20,
    finetune_epochs: int = 20,
    hidden: int = 32,
) -> np.ndarray:
    torch.manual_seed(seed)
    np.random.seed(seed)
    index = {g: i for i, g in enumerate(nodes)}
    train_idx = np.array([index[g] for g in train_genes if g in index and y.loc[g] >= 0], dtype=int)
    data, _ = build_pyg_from_edges(
        nodes,
        features,
        edges,
        feature_cols,
        fit_scaler=True,
        train_idx=train_idx if len(train_idx) else np.arange(len(nodes)),
    )
    encoder = SignedGATv2(data.x.size(1), hidden)
    # Temporarily use hidden via replacing lin during pretrain
    encoder.lin = nn.Identity()
    sign_head = nn.Linear(2 * hidden, 3)
    opt = torch.optim.Adam(list(encoder.parameters()) + list(sign_head.parameters()), lr=1e-3)
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
    # Fine-tune with classification head
    encoder.lin = nn.Linear(hidden, 1)
    opt2 = torch.optim.Adam(encoder.parameters(), lr=1e-3)
    y_full = torch.tensor(y.reindex(nodes).fillna(-1).to_numpy(), dtype=torch.float32)
    for _ in range(finetune_epochs):
        if len(train_idx) < 2:
            break
        encoder.train()
        opt2.zero_grad()
        logits = encoder(data.x, data.edge_index, data.edge_attr).squeeze(-1)
        loss = F.binary_cross_entropy_with_logits(logits[train_idx], y_full[train_idx])
        loss.backward()
        opt2.step()
    encoder.eval()
    with torch.no_grad():
        logits = encoder(data.x, data.edge_index, data.edge_attr).squeeze(-1)
        return torch.sigmoid(logits).cpu().numpy()
