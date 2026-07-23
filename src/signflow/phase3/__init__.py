"""Phase 3: signed directed graph, features, labels, splits (no drug data)."""

from __future__ import annotations

from pathlib import Path

from signflow.phase2 import PROJECT_ROOT, StopAndReportError, ensure_directory, sha256_path

PHASE3_RAW = PROJECT_ROOT / "data" / "raw" / "networks"
PHASE3_OUT = PROJECT_ROOT / "data" / "processed" / "phase3"
FROZEN_PLAN = PROJECT_ROOT / "docs" / "phase0" / "FROZEN_PLAN.yaml"

__all__ = [
    "PROJECT_ROOT",
    "PHASE3_RAW",
    "PHASE3_OUT",
    "FROZEN_PLAN",
    "StopAndReportError",
    "ensure_directory",
    "sha256_path",
]
