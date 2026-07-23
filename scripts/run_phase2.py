from __future__ import annotations

import argparse
import json
import os

# Homebrew R + OpenMP coexistence on macOS
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
os.environ.setdefault("OMP_NUM_THREADS", "1")
if "R_HOME" not in os.environ:
    candidate = "/opt/homebrew/Cellar/r/4.6.1/lib/R"
    if os.path.isdir(candidate):
        os.environ["R_HOME"] = candidate

from signflow.phase2.run import run_phase2


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--track", choices=["target", "benchmark"], required=True)
    parser.add_argument("--allow-dirty", action="store_true")
    parser.add_argument("--n-bootstrap", type=int, default=1000)
    args = parser.parse_args()

    result = run_phase2(track=args.track, allow_dirty=args.allow_dirty, n_bootstrap=args.n_bootstrap)
    print(json.dumps(result.to_dict(), indent=2, sort_keys=True))
    return 1 if result.stop_conditions else 0


if __name__ == "__main__":
    raise SystemExit(main())
