"""Stratified, degree-aware, and leave-one-study-out splits (indices as artifacts)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.model_selection import StratifiedKFold

from . import FROZEN_PLAN, ensure_directory


def _load_folds() -> int:
    return int(yaml.safe_load(FROZEN_PLAN.read_text())["plan"]["phase3_cv_folds"])


def stratified_cv_splits(
    labels: pd.Series,
    *,
    n_folds: int | None = None,
    seed: int = 0,
) -> list[dict[str, list[int]]]:
    n_folds = n_folds or _load_folds()
    y = labels.to_numpy()
    idx = np.arange(len(y))
    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=seed)
    folds = []
    for fold_i, (train, test) in enumerate(skf.split(idx, y)):
        # hold out 20% of train as val
        rng = np.random.default_rng(seed + fold_i)
        train = train.copy()
        rng.shuffle(train)
        n_val = max(1, int(0.2 * len(train)))
        val, tr = train[:n_val], train[n_val:]
        folds.append(
            {
                "fold": fold_i,
                "train": tr.tolist(),
                "val": val.tolist(),
                "test": test.tolist(),
            }
        )
    return folds


def hub_stratified_splits(
    labels: pd.Series,
    degrees: pd.Series,
    *,
    n_folds: int | None = None,
    seed: int = 0,
) -> list[dict[str, list[int]]]:
    """Stratify on label × degree tertile."""

    n_folds = n_folds or _load_folds()
    deg = degrees.reindex(labels.index).fillna(0).to_numpy(dtype=float)
    tertile = pd.qcut(deg, q=3, labels=False, duplicates="drop")
    strata = labels.to_numpy().astype(int) * 10 + tertile.astype(int)
    # Fallback if a stratum too small
    try:
        return stratified_cv_splits(pd.Series(strata, index=labels.index), n_folds=n_folds, seed=seed)
    except ValueError:
        return stratified_cv_splits(labels, n_folds=n_folds, seed=seed)


def leave_one_study_out_inductive(
    node_symbols: list[str],
    study_membership: dict[str, set[str]],
) -> list[dict[str, object]]:
    """BENCHMARK LOSO: train on genes measurable in study A, test on study-B-only relative holdout.

    Practical inductive split: for each held-out study, test indices are nodes whose
    expression evidence comes primarily from that study's DE table membership flag.
    """

    folds = []
    symbols = list(node_symbols)
    index = {g: i for i, g in enumerate(symbols)}
    studies = list(study_membership.keys())
    for held in studies:
        test_genes = [g for g in symbols if g in study_membership[held]]
        train_genes = [g for g in symbols if g not in study_membership[held] or any(
            g in study_membership[s] for s in studies if s != held
        )]
        # Inductive: test = in held study; train = in other study (intersection genes go to train)
        other = [g for g in symbols if any(g in study_membership[s] for s in studies if s != held)]
        test = [g for g in symbols if g in study_membership[held] and g not in set(other)]
        if len(test) < 10:
            # If almost all genes shared, use random half of held-study genes as inductive proxy
            shared = [g for g in symbols if g in study_membership[held] and g in set(other)]
            rng = np.random.default_rng(abs(hash(held)) % (2**31))
            rng.shuffle(shared)
            test = shared[: max(10, len(shared) // 5)]
            train_genes = [g for g in symbols if g not in set(test)]
        else:
            train_genes = [g for g in symbols if g not in set(test)]
        train_idx = [index[g] for g in train_genes if g in index]
        test_idx = [index[g] for g in test if g in index]
        n_val = max(1, int(0.2 * len(train_idx)))
        val_idx, train_idx = train_idx[:n_val], train_idx[n_val:]
        folds.append(
            {
                "held_out_study": held,
                "train": train_idx,
                "val": val_idx,
                "test": test_idx,
            }
        )
    return folds


def assert_disjoint(fold: dict[str, list[int]]) -> None:
    sets = {k: set(fold[k]) for k in ("train", "val", "test")}
    assert sets["train"].isdisjoint(sets["val"])
    assert sets["train"].isdisjoint(sets["test"])
    assert sets["val"].isdisjoint(sets["test"])


def emit_splits(
    output_dir: str | Path,
    labels: pd.Series,
    degrees: pd.Series,
    node_symbols: list[str],
    study_membership: dict[str, set[str]],
    *,
    seed: int = 0,
) -> dict[str, object]:
    output_dir = ensure_directory(output_dir)
    # Stratified/hub splits operate on the labelled subset index space (0..n_labelled-1)
    strat = stratified_cv_splits(labels, seed=seed)
    hub = hub_stratified_splits(labels, degrees.reindex(labels.index).fillna(0), seed=seed + 1)
    # LOSO operates on full node_symbols index space
    loso = leave_one_study_out_inductive(node_symbols, study_membership)
    for fold in strat + hub:
        assert_disjoint(fold)
    for fold in loso:
        assert_disjoint(fold)

    def _pos_counts_labelled(folds: list[dict], y: np.ndarray) -> list[dict]:
        out = []
        for f in folds:
            rec = {"fold": f.get("fold", f.get("held_out_study"))}
            for split in ("train", "val", "test"):
                idx = f[split]
                rec[f"n_{split}"] = len(idx)
                rec[f"n_pos_{split}"] = int(y[idx].sum()) if len(idx) else 0
            out.append(rec)
        return out

    def _pos_counts_loso(folds: list[dict], node_symbols: list[str], label_by_symbol: pd.Series) -> list[dict]:
        out = []
        for f in folds:
            rec = {"fold": f.get("held_out_study")}
            for split in ("train", "val", "test"):
                idx = f[split]
                rec[f"n_{split}"] = len(idx)
                genes = [node_symbols[i] for i in idx]
                vals = label_by_symbol.reindex(genes).fillna(-1)
                rec[f"n_pos_{split}"] = int((vals == 1).sum())
            out.append(rec)
        return out

    y = labels.to_numpy()
    # Full-graph label series for LOSO reporting (unlabelled = -1)
    full_y = pd.Series(-1, index=node_symbols, dtype=int)
    full_y.loc[full_y.index.isin(labels.index[labels == 1])] = 1
    full_y.loc[full_y.index.isin(labels.index[labels == 0])] = 0
    payload = {
        "stratified_cv": strat,
        "hub_stratified_cv": hub,
        "leave_one_study_out": loso,
        "positive_counts_stratified": _pos_counts_labelled(strat, y),
        "positive_counts_hub": _pos_counts_labelled(hub, y),
        "positive_counts_loso": _pos_counts_loso(loso, node_symbols, full_y),
        "note": "stratified/hub indices refer to labelled subset order; LOSO indices refer to full node_symbols order",
    }
    (Path(output_dir) / "splits.json").write_text(json.dumps(payload, indent=2) + "\n")
    return payload
