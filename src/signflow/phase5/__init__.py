"""Phase 5: signed signaling-flow drug scoring (BENCHMARK + TARGET)."""

from __future__ import annotations

from pathlib import Path

from signflow.phase2 import PROJECT_ROOT, StopAndReportError, ensure_directory, sha256_path

PHASE3_LATEST = PROJECT_ROOT / "data" / "processed" / "phase3" / "latest"
PHASE5_OUT = PROJECT_ROOT / "data" / "processed" / "phase5"
PHASE5_RAW = PROJECT_ROOT / "data" / "raw" / "drugs"
FROZEN_PLAN = PROJECT_ROOT / "docs" / "phase0" / "FROZEN_PLAN.yaml"
BENCHMARK_DRUGS = PROJECT_ROOT / "configs" / "benchmark_drugs.yaml"
SIGNATURE_BENCHMARK = PROJECT_ROOT / "data" / "processed" / "phase2" / "meta" / "signature_meta_benchmark.tsv"
SIGNATURE_TARGET = PROJECT_ROOT / "data" / "processed" / "phase2" / "meta" / "signature_meta.tsv"

__all__ = [
    "PROJECT_ROOT",
    "PHASE3_LATEST",
    "PHASE5_OUT",
    "PHASE5_RAW",
    "FROZEN_PLAN",
    "BENCHMARK_DRUGS",
    "SIGNATURE_BENCHMARK",
    "SIGNATURE_TARGET",
    "StopAndReportError",
    "ensure_directory",
    "sha256_path",
]
