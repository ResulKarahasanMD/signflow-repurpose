"""GEO acquisition helpers for Phase 2 expression data."""

from __future__ import annotations

import re
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path
from urllib.parse import urljoin

import requests

from . import DownloadRecord, RAW_GEO_ROOT, StopAndReportError, ensure_directory, sha256_bytes, sha256_path

REQUEST_TIMEOUT_SECONDS = 120
ALLOWED_SUPPL_SUFFIXES = (".gz", ".tar", ".txt", ".csv", ".tsv", ".soft", ".cel", ".zip")


def utc_now_iso() -> str:
    return datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def geo_series_prefix(accession: str) -> str:
    """Return the GEO series path prefix, e.g. GSE11783 -> GSE11nnn."""

    return f"{accession[:-3]}nnn"


def ftp_to_https(url: str) -> str:
    """Normalize GEO FTP links to HTTPS."""

    return url.replace("ftp://ftp.ncbi.nlm.nih.gov/", "https://ftp.ncbi.nlm.nih.gov/")


def series_matrix_url(accession: str) -> str:
    """Build the canonical series matrix URL for a GEO series."""

    prefix = geo_series_prefix(accession)
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{accession}/matrix/{accession}_series_matrix.txt.gz"


def supplementary_index_url(accession: str) -> str:
    """Return the supplementary directory listing URL for a GEO series."""

    prefix = geo_series_prefix(accession)
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{accession}/suppl/"


def download_remote_file(url: str, destination: str | Path, *, accession: str) -> DownloadRecord:
    """Download a remote GEO artifact and return its hashed record."""

    accessed_at = utc_now_iso()
    response = requests.get(url, timeout=REQUEST_TIMEOUT_SECONDS)
    if response.status_code != 200:
        raise StopAndReportError(
            f"Required Phase 2 file is unavailable: {url}",
            accession=accession,
            reason="unavailable-files",
            details={"url": url, "status_code": response.status_code},
        )

    destination = Path(destination)
    ensure_directory(destination.parent)
    payload = response.content
    destination.write_bytes(payload)
    return DownloadRecord(
        accession=accession,
        url=url,
        path=str(destination),
        sha256=sha256_bytes(payload),
        size_bytes=len(payload),
        status="downloaded",
        accessed_at_utc=accessed_at,
    )


def reuse_or_download(url: str, destination: str | Path, *, accession: str) -> DownloadRecord:
    """Reuse a cached file if present, otherwise download it."""

    destination = Path(destination)
    accessed_at = utc_now_iso()
    if destination.exists():
        return DownloadRecord(
            accession=accession,
            url=url,
            path=str(destination),
            sha256=sha256_path(destination),
            size_bytes=destination.stat().st_size,
            status="reused",
            accessed_at_utc=accessed_at,
        )
    return download_remote_file(url, destination, accession=accession)


def _is_geo_supplementary_href(href: str, accession: str) -> bool:
    lowered = href.lower()
    if href in {"../", "/"} or href.startswith("?") or href.startswith("#"):
        return False
    if "hhs.gov" in lowered or "vulnerability" in lowered:
        return False
    if not any(lowered.endswith(suffix) for suffix in ALLOWED_SUPPL_SUFFIXES):
        return False
    # Prefer files that mention the accession or look like GEO artifacts.
    name = Path(href.split("?")[0]).name.lower()
    return accession.lower() in name or name.startswith("gsm") or name in {"filelist.txt"}


def discover_supplementary_urls(accession: str) -> list[str]:
    """Return supplementary file URLs from the GEO directory index."""

    index_url = supplementary_index_url(accession)
    response = requests.get(index_url, timeout=REQUEST_TIMEOUT_SECONDS)
    if response.status_code != 200:
        raise StopAndReportError(
            f"Unable to list supplementary files for {accession}",
            accession=accession,
            reason="unavailable-files",
            details={"url": index_url, "status_code": response.status_code},
        )

    hrefs = re.findall(r'href="([^"]+)"', response.text)
    urls = []
    for href in hrefs:
        if not _is_geo_supplementary_href(href, accession):
            continue
        urls.append(urljoin(index_url, href))

    if not urls:
        raise StopAndReportError(
            f"No supplementary files were discoverable for {accession}",
            accession=accession,
            reason="unavailable-files",
            details={"url": index_url},
        )
    return sorted(dict.fromkeys(urls))


def acquire_geo_series(accession: str, raw_root: str | Path = RAW_GEO_ROOT) -> dict[str, object]:
    """Download the series matrix and all supplementary GEO files."""

    series_root = ensure_directory(Path(raw_root) / accession)
    matrix_url = series_matrix_url(accession)
    matrix_path = series_root / "matrix" / Path(matrix_url).name
    records = [reuse_or_download(matrix_url, matrix_path, accession=accession)]

    supplementary_paths: list[str] = []
    for url in discover_supplementary_urls(accession):
        normalized_url = ftp_to_https(url)
        destination = series_root / "suppl" / Path(normalized_url).name
        record = reuse_or_download(normalized_url, destination, accession=accession)
        records.append(record)
        supplementary_paths.append(record.path)

    return {
        "accession": accession,
        "series_matrix_path": str(matrix_path),
        "supplementary_paths": supplementary_paths,
        "download_records": [asdict(record) for record in records],
    }
