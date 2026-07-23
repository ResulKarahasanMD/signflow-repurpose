"""Phase 4 orchestrator."""

from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

from . import PHASE4_OUT, ensure_directory
from .data import assert_splits_hash, load_phase3_benchmark
from .evaluate import run_benchmark_suite
from .target_prop import run_target_propagation


def run_phase4(phase3_root: str | Path | None = None) -> dict[str, object]:
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    out = ensure_directory(PHASE4_OUT / stamp)
    bundle = load_phase3_benchmark(phase3_root)
    splits_hash = assert_splits_hash(bundle)
    (out / "splits_hash.txt").write_text(splits_hash + "\n")

    bench = run_benchmark_suite(bundle, out / "benchmark")
    target = run_target_propagation(out / "target", Path(bundle["root"]))

    result = {
        "started_at_utc": stamp,
        "phase3_root": bundle["root"],
        "splits_hash": splits_hash,
        "drug_data_touched": False,
        "h1_w1": bench.get("h1_w1"),
        "h2": bench.get("h2"),
        "conformal": bench.get("conformal"),
        "target_propagation": {
            "n_overlap_top50": target.get("n_overlap"),
            "summary_path": str(Path(out) / "target" / "target_propagation_summary.json"),
        },
        "output_dir": str(out),
    }
    (out / "phase4_result.json").write_text(json.dumps(result, indent=2, default=str) + "\n")
    latest = PHASE4_OUT / "latest"
    if latest.exists() or latest.is_symlink():
        latest.unlink()
    latest.symlink_to(out.name)
    return result
