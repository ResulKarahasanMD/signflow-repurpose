#!/usr/bin/env python3
"""Phase 2.6 TARGET replication-honesty diagnostics + annotated seed artifacts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from signflow.phase2 import PHASE2_META_ROOT, PHASE2_OUTPUT_ROOT, ensure_directory
from signflow.phase2.seeds import emit_target_seed_artifacts
from signflow.phase2.target_diagnostics import run_target_diagnostics


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PHASE2_META_ROOT / "diagnostics_phase26",
    )
    args = parser.parse_args()
    out = ensure_directory(args.output_dir)
    summary = run_target_diagnostics(
        PHASE2_META_ROOT / "signature_meta.tsv",
        PHASE2_OUTPUT_ROOT / "GSE11783" / "de_limma.tsv",
        PHASE2_OUTPUT_ROOT / "GSE57560" / "de_limma.tsv",
        out,
    )
    seeds = emit_target_seed_artifacts(
        PHASE2_META_ROOT / "signature_meta.tsv",
        PHASE2_OUTPUT_ROOT / "GSE11783" / "de_limma.tsv",
        PHASE2_OUTPUT_ROOT / "GSE57560" / "de_limma.tsv",
        out,
    )
    payload = {"diagnostics": summary, "seeds": seeds}
    (out / "phase26_summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
