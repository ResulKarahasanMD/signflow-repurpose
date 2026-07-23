"""Phase 4: baselines, relational GNN, graph controls, Mondrian conformal (no drug data)."""

from __future__ import annotations

from pathlib import Path

from signflow.phase2 import PROJECT_ROOT, StopAndReportError, ensure_directory, sha256_path

PHASE3_LATEST = PROJECT_ROOT / "data" / "processed" / "phase3" / "latest"
PHASE4_OUT = PROJECT_ROOT / "data" / "processed" / "phase4"
FROZEN_PLAN = PROJECT_ROOT / "docs" / "phase0" / "FROZEN_PLAN.yaml"

__all__ = [
    "PROJECT_ROOT",
    "PHASE3_LATEST",
    "PHASE4_OUT",
    "FROZEN_PLAN",
    "StopAndReportError",
    "ensure_directory",
    "sha256_path",
]
