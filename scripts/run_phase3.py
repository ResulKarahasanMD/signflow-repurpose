#!/usr/bin/env python3
"""Run Phase 3 signed-graph / features / labels / splits build."""

from __future__ import annotations

import argparse
import json

from signflow.phase3.run import run_phase3


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--track", choices=["target", "benchmark", "both"], default="both")
    args = parser.parse_args()
    result = run_phase3(args.track)
    print(json.dumps(result, indent=2, default=str))


if __name__ == "__main__":
    main()
