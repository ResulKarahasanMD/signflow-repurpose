"""Mondrian (class-conditional) conformal prediction for gene prioritization."""

from __future__ import annotations

import numpy as np


def mondrian_conformal_sets(
    cal_scores: np.ndarray,
    cal_labels: np.ndarray,
    test_scores: np.ndarray,
    *,
    alpha: float = 0.1,
) -> tuple[list[set[int]], dict[str, float]]:
    """Return prediction sets over classes {{0,1}} for each test point + diagnostics.

    Nonconformity = 1 - p_y for the true/candidate class (score treated as P(y=1)).
    """

    cal_scores = np.asarray(cal_scores, dtype=float)
    cal_labels = np.asarray(cal_labels, dtype=int)
    test_scores = np.asarray(test_scores, dtype=float)

    def _nc(scores: np.ndarray, labels: np.ndarray) -> np.ndarray:
        # For class 1: 1 - score; for class 0: score
        return np.where(labels == 1, 1.0 - scores, scores)

    sets: list[set[int]] = []
    thresholds = {}
    for cls in (0, 1):
        mask = cal_labels == cls
        if not mask.any():
            thresholds[cls] = 1.0
            continue
        nc = _nc(cal_scores[mask], np.full(mask.sum(), cls))
        # quantile level
        q_level = min(1.0, np.ceil((1 - alpha) * (len(nc) + 1)) / len(nc))
        thresholds[cls] = float(np.quantile(nc, q_level, method="higher")) if len(nc) else 1.0

    for s in test_scores:
        pred_set = set()
        for cls in (0, 1):
            nc = (1.0 - s) if cls == 1 else s
            if nc <= thresholds[cls]:
                pred_set.add(cls)
        if not pred_set:
            # fallback: include argmax
            pred_set.add(1 if s >= 0.5 else 0)
        sets.append(pred_set)

    mean_size = float(np.mean([len(s) for s in sets]))
    return sets, {"alpha": alpha, "mean_set_size": mean_size, "thresholds": thresholds}


def empirical_coverage(sets: list[set[int]], y_true: np.ndarray) -> float:
    y_true = np.asarray(y_true, dtype=int)
    hits = [int(y_true[i]) in sets[i] for i in range(len(y_true))]
    return float(np.mean(hits)) if hits else float("nan")


NONCLAIM = (
    "Conformal coverage is over the gene/drug population under the labelled-data distribution. "
    "It does not quantify sampling uncertainty over patients, and does not imply the disease "
    "signature generalises to new IC/BPS cohorts."
)
