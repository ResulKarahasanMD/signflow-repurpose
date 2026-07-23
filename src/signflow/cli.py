"""Command line interface for the Phase 1 harness."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .config import ConfigError, load_config
from .pipeline import StageNotImplementedError, run_smoke, run_stage
from .provenance import DirtyRepositoryError, write_run_manifest
from .resources import PACKAGE_VERSION_NAMES, verify_resources

DEFAULT_CONFIG = "configs/default.yaml"


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="signflow")
    parser.add_argument("-c", "--config", default=DEFAULT_CONFIG)
    parser.add_argument("--set", action="append", default=[], dest="overrides")
    parser.add_argument("--allow-dirty", action="store_true")

    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("verify-resources")
    subparsers.add_parser("smoke")
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--stage", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)

    try:
        config = load_config(args.config, args.overrides)
    except ConfigError as exc:
        print(f"Config error: {exc}", file=sys.stderr)
        return 2

    try:
        if args.command == "verify-resources":
            records = verify_resources(config)
            write_run_manifest(
                Path("runs/verify-resources/run_manifest.json"),
                config=config,
                input_paths=[args.config],
                resource_records=records,
                package_names=PACKAGE_VERSION_NAMES,
                repo_root=".",
                allow_dirty=args.allow_dirty,
            )
            print(json.dumps(records, indent=2, sort_keys=True))
            has_failure = any((not record["ok"]) or record["status"] in {"drift", "error"} for record in records)
            return 1 if has_failure else 0

        if args.command == "smoke":
            result = run_smoke(config)
            write_run_manifest(
                Path("runs/smoke/run_manifest.json"),
                config=config,
                input_paths=[args.config],
                resource_records=[],
                package_names=PACKAGE_VERSION_NAMES,
                repo_root=".",
                allow_dirty=args.allow_dirty,
            )
            print(json.dumps(result, indent=2, sort_keys=True))
            return 0

        if args.command == "run":
            run_stage(args.stage, config)
            return 0

        parser.error(f"Unknown command: {args.command}")
        return 2
    except ConfigError as exc:
        print(f"Config error: {exc}", file=sys.stderr)
        return 2
    except DirtyRepositoryError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    except StageNotImplementedError as exc:
        print(str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
