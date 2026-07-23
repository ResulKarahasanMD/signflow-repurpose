#!/usr/bin/env python3
"""Run Phase 4 baselines / GNN / controls / conformal / TARGET propagation."""

from __future__ import annotations

import argparse
import json

from signflow.phase4.run import run_phase4


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase3-root", default=None)
    args = parser.parse_args()
    print(json.dumps(run_phase4(args.phase3_root), indent=2, default=str))


if __name__ == "__main__":
    main()
