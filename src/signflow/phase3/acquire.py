"""Acquire OmniPath, SIGNOR, and STRING network resources (human only)."""

from __future__ import annotations

import gzip
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd
import requests

from . import PHASE3_RAW, StopAndReportError, ensure_directory, sha256_path

REQUEST_TIMEOUT = 600

OMNIPATH_URL = (
    "https://omnipathdb.org/interactions"
    "?datasets=omnipath&genesymbols=1&organisms=9606&fields=sources,references"
)
SIGNOR_URL = "https://signor.uniroma2.it/releases/getLatestRelease.php"
STRING_PHYSICAL_URL = (
    "https://stringdb-downloads.org/download/protein.physical.links.v12.0/"
    "9606.protein.physical.links.v12.0.txt.gz"
)
STRING_INFO_URL = (
    "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz"
)


def _download(url: str, destination: Path, *, binary: bool = True) -> dict[str, str]:
    ensure_directory(destination.parent)
    accessed = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    if destination.exists() and destination.stat().st_size > 0:
        return {
            "url": url,
            "path": str(destination),
            "accessed_at_utc": accessed,
            "sha256": sha256_path(destination),
            "status": "reused",
            "size_bytes": str(destination.stat().st_size),
        }
    response = requests.get(url, timeout=REQUEST_TIMEOUT, stream=True)
    if response.status_code != 200:
        raise StopAndReportError(
            f"Network resource unavailable: {url} HTTP {response.status_code}",
            reason="network-resource-unavailable",
            details={"url": url, "status_code": response.status_code},
        )
    with destination.open("wb") as handle:
        for chunk in response.iter_content(chunk_size=1 << 20):
            if chunk:
                handle.write(chunk)
    return {
        "url": url,
        "path": str(destination),
        "accessed_at_utc": accessed,
        "sha256": sha256_path(destination),
        "status": "downloaded",
        "size_bytes": str(destination.stat().st_size),
    }


def acquire_networks(raw_root: str | Path = PHASE3_RAW) -> dict[str, object]:
    """Download/cached OmniPath, SIGNOR, STRING physical links + protein info."""

    raw_root = ensure_directory(raw_root)
    omnipath_path = Path(raw_root) / "omnipath" / "interactions_omnipath_9606.tsv"
    signor_path = Path(raw_root) / "signor" / "signor_latest.txt"
    string_links = Path(raw_root) / "string" / "9606.protein.physical.links.v12.0.txt.gz"
    string_info = Path(raw_root) / "string" / "9606.protein.info.v12.0.txt.gz"

    records = {
        "omnipath": _download(OMNIPATH_URL, omnipath_path),
        "signor": _download(SIGNOR_URL, signor_path),
        "string_physical": _download(STRING_PHYSICAL_URL, string_links),
        "string_info": _download(STRING_INFO_URL, string_info),
    }

    # Validate OmniPath has required sign/direction fields
    omni = pd.read_csv(omnipath_path, sep="\t", nrows=5)
    required = {"source_genesymbol", "target_genesymbol", "is_directed", "is_stimulation", "is_inhibition"}
    missing = required - set(omni.columns)
    if missing:
        raise StopAndReportError(
            f"OmniPath missing sign/direction fields: {sorted(missing)}",
            reason="omnipath-fields-absent",
            details={"columns": omni.columns.tolist()},
        )

    signor = pd.read_csv(signor_path, sep="\t", nrows=5)
    if "EFFECT" not in signor.columns or "IDA" not in signor.columns:
        raise StopAndReportError(
            "SIGNOR missing EFFECT/IDA fields",
            reason="signor-fields-absent",
            details={"columns": signor.columns.tolist()},
        )

    return {
        "manifest": records,
        "paths": {
            "omnipath": str(omnipath_path),
            "signor": str(signor_path),
            "string_physical": str(string_links),
            "string_info": str(string_info),
        },
        "versions": {
            "omnipath_dataset": "omnipath",
            "signor": "getLatestRelease.php",
            "string": "v12.0 physical links",
            "taxon": 9606,
        },
    }
