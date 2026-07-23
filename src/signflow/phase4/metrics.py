"""Classification and ranking metrics with bootstrap CIs."""

from __future__ import annotations

import numpy as np
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    matthews_corrcoef,
    precision_score,
    recall_score,
    roc_auc_score,
)


def expected_calibration_error(y_true: np.ndarray, y_prob: np.ndarray, n_bins: int = 10) -> float:
    bins = np.linspace(0.0, 1.0, n_bins + 1)
    ece = 0.0
    for i in range(n_bins):
        mask = (y_prob >= bins[i]) & (y_prob < bins[i + 1] if i < n_bins - 1 else y_prob <= bins[i + 1])
        if not mask.any():
            continue
        ece += mask.mean() * abs(y_true[mask].mean() - y_prob[mask].mean())
    return float(ece)


def precision_at_k(y_true: np.ndarray, scores: np.ndarray, k: int) -> float:
    k = min(k, len(scores))
    if k == 0:
        return float("nan")
    top = np.argsort(scores)[::-1][:k]
    return float(y_true[top].mean())


def enrichment_at_k(y_true: np.ndarray, scores: np.ndarray, k: int) -> float:
    base = float(y_true.mean()) if len(y_true) else float("nan")
    if not base:
        return float("nan")
    return precision_at_k(y_true, scores, k) / base


def compute_metrics(y_true: np.ndarray, scores: np.ndarray, *, k: int = 50) -> dict[str, float]:
    y_true = np.asarray(y_true).astype(int)
    scores = np.asarray(scores).astype(float)
    # Map scores to [0,1] probs via rank for ECE when unbounded
    ranks = scores.argsort().argsort().astype(float)
    probs = (ranks + 1) / (len(scores) + 1)
    pred = (scores >= np.median(scores)).astype(int)
    out = {
        "auroc": float(roc_auc_score(y_true, scores)) if len(np.unique(y_true)) > 1 else float("nan"),
        "auprc": float(average_precision_score(y_true, scores)) if len(np.unique(y_true)) > 1 else float("nan"),
        "mcc": float(matthews_corrcoef(y_true, pred)) if len(np.unique(y_true)) > 1 else float("nan"),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, pred)),
        "ece": expected_calibration_error(y_true, probs),
        "precision_at_k": precision_at_k(y_true, scores, k),
        "recall_at_k": float(y_true[np.argsort(scores)[::-1][:k]].sum() / max(y_true.sum(), 1)),
        "enrichment_at_k": enrichment_at_k(y_true, scores, k),
        "n_pos": int(y_true.sum()),
        "n": int(len(y_true)),
    }
    return out


def bootstrap_ci(
    y_true: np.ndarray,
    scores: np.ndarray,
    *,
    metric: str = "auprc",
    n_boot: int = 100,
    seed: int = 0,
    max_n: int = 2000,
) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true)
    scores = np.asarray(scores)
    if len(y_true) > max_n:
        base = rng.choice(len(y_true), size=max_n, replace=False)
        y_true = y_true[base]
        scores = scores[base]
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(y_true), size=len(y_true))
        if len(np.unique(y_true[idx])) < 2:
            continue
        m = compute_metrics(y_true[idx], scores[idx])
        vals.append(m[metric])
    if not vals:
        return {"mean": float("nan"), "lo": float("nan"), "hi": float("nan")}
    arr = np.asarray(vals, dtype=float)
    return {
        "mean": float(np.mean(arr)),
        "lo": float(np.quantile(arr, 0.025)),
        "hi": float(np.quantile(arr, 0.975)),
    }


def ci_overlap(a: dict[str, float], b: dict[str, float]) -> bool:
    return not (a["hi"] < b["lo"] or b["hi"] < a["lo"])
