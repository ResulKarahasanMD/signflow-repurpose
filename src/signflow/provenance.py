"""Run-manifest creation and provenance helpers."""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import UTC, datetime
from importlib import metadata
from pathlib import Path
from typing import Any

from .config import get_seed_list, plan_hash, resolved_config_hash


class DirtyRepositoryError(RuntimeError):
    """Raised when writing a manifest from a dirty repository."""


class GitStateError(RuntimeError):
    """Raised when git metadata cannot be resolved for a manifest write."""


def utc_now_iso() -> str:
    return datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def utc_now_path_component() -> str:
    return datetime.now(UTC).strftime("%Y%m%dT%H%M%S%fZ")


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def package_versions(package_names: list[str]) -> dict[str, str]:
    versions: dict[str, str] = {}
    for name in package_names:
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = "not-installed"
    return versions


def git_state(repo_root: str | Path = ".") -> dict[str, Any]:
    repo_root = str(repo_root)
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo_root,
        capture_output=True,
        check=False,
        text=True,
    )
    if commit.returncode != 0 or not commit.stdout.strip():
        detail = commit.stderr.strip() or commit.stdout.strip() or "git rev-parse HEAD failed"
        raise GitStateError(f"Unable to resolve git commit in {repo_root}: {detail}")

    status = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=repo_root,
        capture_output=True,
        check=False,
        text=True,
    )
    if status.returncode != 0:
        detail = status.stderr.strip() or status.stdout.strip() or "git status --porcelain failed"
        raise GitStateError(f"Unable to resolve git status in {repo_root}: {detail}")

    return {
        "commit": commit.stdout.strip(),
        "dirty": bool(status.stdout.strip()),
        "status_lines": [line for line in status.stdout.splitlines() if line.strip()],
    }


def build_manifest(
    *,
    config: dict[str, Any],
    input_paths: list[str | Path],
    resource_records: list[dict[str, Any]],
    package_names: list[str],
    git: dict[str, Any],
) -> dict[str, Any]:
    return {
        "created_at_utc": utc_now_iso(),
        "plan_hash": plan_hash(config),
        "resolved_config_hash": resolved_config_hash(config),
        "resolved_config": config,
        "runtime_environment": config.get("runtime", {}),
        "package_versions": package_versions(package_names),
        "seeds": get_seed_list(config),
        "input_sha256": {str(path): sha256_file(path) for path in input_paths},
        "resource_records": resource_records,
        "git": {"commit": git["commit"], "dirty": git["dirty"]},
    }


def write_run_manifest(
    output_root: str | Path,
    *,
    config: dict[str, Any],
    input_paths: list[str | Path],
    resource_records: list[dict[str, Any]],
    package_names: list[str],
    repo_root: str | Path = ".",
    allow_dirty: bool = False,
) -> dict[str, Any]:
    git = git_state(repo_root)
    if git["dirty"] and not allow_dirty:
        raise DirtyRepositoryError("Repository is dirty; re-run with allow_dirty=True to continue.")

    manifest = build_manifest(
        config=config,
        input_paths=input_paths,
        resource_records=resource_records,
        package_names=package_names,
        git=git,
    )
    output = Path(output_root) / f"{utc_now_path_component()}_{manifest['plan_hash'][:8]}" / "run_manifest.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest
