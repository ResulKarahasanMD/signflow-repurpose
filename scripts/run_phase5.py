#!/usr/bin/env python3
"""Run Phase 5 drug scoring pipeline."""

from __future__ import annotations

import argparse
import json

from signflow.phase5.run import run_phase5


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase3-root", default=None)
    parser.add_argument("--raw-root", default=None)
    parser.add_argument("--skip-acquire", action="store_true")
    args = parser.parse_args()
    print(
        json.dumps(
            run_phase5(phase3_root=args.phase3_root, raw_root=args.raw_root, skip_acquire=args.skip_acquire),
            indent=2,
            default=str,
        )
    )


if __name__ == "__main__":
    main()
