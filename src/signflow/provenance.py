"""Run-manifest creation and provenance helpers."""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import UTC, datetime
from importlib import metadata
from pathlib import Path
from typing import Any

from .config import config_hash, get_seed_list


class DirtyRepositoryError(RuntimeError):
    """Raised when writing a manifest from a dirty repository."""


def utc_now_iso() -> str:
    return datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


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
    status = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=repo_root,
        capture_output=True,
        check=False,
        text=True,
    )
    return {
        "commit": commit.stdout.strip() or "unknown",
        "dirty": bool(status.stdout.strip()),
        "status_lines": [line for line in status.stdout.splitlines() if line.strip()],
    }


def build_manifest(
    *,
    config: dict[str, Any],
    input_paths: list[str | Path],
    resource_records: list[dict[str, Any]],
    package_names: list[str],
    repo_root: str | Path = ".",
) -> dict[str, Any]:
    git = git_state(repo_root)
    return {
        "created_at_utc": utc_now_iso(),
        "config_hash": config_hash(config),
        "resolved_config": config,
        "package_versions": package_versions(package_names),
        "seeds": get_seed_list(config),
        "input_sha256": {str(path): sha256_file(path) for path in input_paths},
        "resource_records": resource_records,
        "git": {"commit": git["commit"], "dirty": git["dirty"]},
    }


def write_run_manifest(
    output_path: str | Path,
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
        repo_root=repo_root,
    )
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest
