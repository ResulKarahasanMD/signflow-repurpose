"""Classification and ranking metrics with bootstrap CIs."""

from __future__ import annotations

import numpy as np
import yaml
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    matthews_corrcoef,
    roc_auc_score,
)

from . import FROZEN_PLAN


_BOOTSTRAP_N: int | None = None


def load_bootstrap_n() -> int:
    global _BOOTSTRAP_N
    if _BOOTSTRAP_N is None:
        _BOOTSTRAP_N = int(yaml.safe_load(FROZEN_PLAN.read_text())["plan"]["phase4_bootstrap_n"])
    return _BOOTSTRAP_N


def expected_calibration_error(y_true: np.ndarray, y_prob: np.ndarray, n_bins: int = 10) -> float:
    """ECE on actual predicted probabilities in [0, 1]."""

    y_prob = np.asarray(y_prob, dtype=float)
    y_true = np.asarray(y_true, dtype=float)
    if np.any((y_prob < 0) | (y_prob > 1)):
        raise ValueError("ECE requires probabilities in [0, 1]; omit ECE for unbounded scores")
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


def compute_metrics(
    y_true: np.ndarray,
    scores: np.ndarray,
    *,
    k: int = 50,
    probabilities: np.ndarray | None = None,
    regime: str = "degree_matched",
) -> dict[str, float]:
    y_true = np.asarray(y_true).astype(int)
    scores = np.asarray(scores).astype(float)
    out: dict[str, float] = {
        "auroc": float(roc_auc_score(y_true, scores)) if len(np.unique(y_true)) > 1 else float("nan"),
        "auprc": float(average_precision_score(y_true, scores)) if len(np.unique(y_true)) > 1 else float("nan"),
        "precision_at_k": precision_at_k(y_true, scores, k),
        "recall_at_k": float(y_true[np.argsort(scores)[::-1][:k]].sum() / max(y_true.sum(), 1)),
        "enrichment_at_k": enrichment_at_k(y_true, scores, k),
        "n_pos": int(y_true.sum()),
        "n": int(len(y_true)),
        "chance_auprc": float(y_true.mean()) if len(y_true) else float("nan"),
    }
    # MCC / balanced accuracy only in balanced degree-matched regime with a principled threshold
    if regime == "degree_matched" and probabilities is not None:
        pred = (np.asarray(probabilities) >= 0.5).astype(int)
        out["mcc"] = float(matthews_corrcoef(y_true, pred)) if len(np.unique(y_true)) > 1 else float("nan")
        out["balanced_accuracy"] = float(balanced_accuracy_score(y_true, pred))
        try:
            out["ece"] = expected_calibration_error(y_true, np.asarray(probabilities))
        except ValueError:
            out["ece"] = float("nan")
            out["ece_omitted_reason"] = "probabilities_out_of_range"  # type: ignore[assignment]
    else:
        out["mcc"] = float("nan")
        out["balanced_accuracy"] = float("nan")
        out["ece"] = float("nan")
        out["ece_omitted_reason"] = "unbounded_scores_or_pu_regime"  # type: ignore[assignment]
    return out


def bootstrap_ci(
    y_true: np.ndarray,
    scores: np.ndarray,
    *,
    metric: str = "auprc",
    n_boot: int | None = None,
    seed: int = 0,
) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true)
    scores = np.asarray(scores)
    n_boot = int(n_boot if n_boot is not None else load_bootstrap_n())
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(y_true), size=len(y_true))
        if len(np.unique(y_true[idx])) < 2:
            continue
        m = compute_metrics(y_true[idx], scores[idx], regime="degree_matched")
        vals.append(m[metric])
    if not vals:
        return {"mean": float("nan"), "lo": float("nan"), "hi": float("nan"), "n_boot": 0, "width": float("nan")}
    arr = np.asarray(vals, dtype=float)
    lo, hi = float(np.quantile(arr, 0.025)), float(np.quantile(arr, 0.975))
    return {
        "mean": float(np.mean(arr)),
        "lo": lo,
        "hi": hi,
        "n_boot": int(len(vals)),
        "width": hi - lo,
    }


def ci_overlap(a: dict[str, float], b: dict[str, float]) -> bool:
    return not (a["hi"] < b["lo"] or b["hi"] < a["lo"])


def paired_difference_bootstrap(
    y_true: np.ndarray,
    scores_a: np.ndarray,
    scores_b: np.ndarray,
    *,
    metric: str = "auprc",
    n_boot: int | None = None,
    seed: int = 0,
) -> dict[str, float]:
    """Paired bootstrap of metric(A) − metric(B) over the same index resamples.

    Pre-committed exploratory contrast; W1 CI-overlap verdict remains binding.
    """

    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true)
    scores_a = np.asarray(scores_a, dtype=float)
    scores_b = np.asarray(scores_b, dtype=float)
    n_boot = int(n_boot if n_boot is not None else load_bootstrap_n())
    diffs: list[float] = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(y_true), size=len(y_true))
        if len(np.unique(y_true[idx])) < 2:
            continue
        ma = compute_metrics(y_true[idx], scores_a[idx], regime="degree_matched")[metric]
        mb = compute_metrics(y_true[idx], scores_b[idx], regime="degree_matched")[metric]
        diffs.append(float(ma - mb))
    if not diffs:
        return {
            "mean": float("nan"),
            "lo": float("nan"),
            "hi": float("nan"),
            "n_boot": 0,
            "width": float("nan"),
            "p_le_zero": float("nan"),
            "status": "exploratory_paired_difference",
        }
    arr = np.asarray(diffs, dtype=float)
    lo, hi = float(np.quantile(arr, 0.025)), float(np.quantile(arr, 0.975))
    return {
        "mean": float(np.mean(arr)),
        "lo": lo,
        "hi": hi,
        "n_boot": int(len(arr)),
        "width": hi - lo,
        "p_le_zero": float(np.mean(arr <= 0.0)),
        "status": "exploratory_paired_difference",
        "note": "Exploratory; W1 CI-overlap verdict remains binding.",
    }


def minimum_detectable_difference(ci_width_a: float, ci_width_b: float) -> float:
    """Approximate gap needed for non-overlap if intervals are centered on equal means.

    For two symmetric CIs of half-widths ha, hb, non-overlap requires |Δ| > ha + hb.
    Using full widths: MDD ≈ (width_a + width_b) / 2.
    """

    return float((ci_width_a + ci_width_b) / 2.0)
