"""Authoritative HGNC biotype/alias reference for Phase 2 gene collapse."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd
import requests

from . import PROJECT_ROOT, StopAndReportError, ensure_directory, sha256_path

HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
CACHE_DIR = PROJECT_ROOT / "data" / "raw" / "reference" / "hgnc"
REQUEST_TIMEOUT = 180


@dataclass(frozen=True)
class BiotypeResolution:
    symbol: str
    status: str  # matched | alias_resolved | unmatched | non_coding_rejected | biotype_unknown
    approved_symbol: str | None
    locus_group: str | None


class HgncReference:
    def __init__(self, frame: pd.DataFrame, *, source_url: str, path: str, accessed_at_utc: str, sha256: str):
        self.frame = frame
        self.source_url = source_url
        self.path = path
        self.accessed_at_utc = accessed_at_utc
        self.sha256 = sha256
        self._by_symbol = {}
        self._alias_to_symbol = {}
        for _, row in frame.iterrows():
            approved = str(row["symbol"]).strip()
            locus_group = str(row.get("locus_group") or "").strip().lower()
            self._by_symbol[approved.upper()] = (approved, locus_group)
            aliases = str(row.get("alias_symbol") or "")
            prev = str(row.get("prev_symbol") or "")
            for chunk in (aliases + "|" + prev).replace(",", "|").split("|"):
                alias = chunk.strip()
                if alias:
                    self._alias_to_symbol[alias.upper()] = approved

    def manifest_record(self) -> dict[str, str]:
        return {
            "resource": "HGNC complete set",
            "url": self.source_url,
            "path": self.path,
            "accessed_at_utc": self.accessed_at_utc,
            "sha256": self.sha256,
            "n_rows": str(len(self.frame)),
        }

    def resolve(self, symbol: str) -> BiotypeResolution:
        key = symbol.strip().upper()
        if not key or key in {"---", "NA"}:
            return BiotypeResolution(symbol, "unmatched", None, None)
        if key in self._by_symbol:
            approved, locus = self._by_symbol[key]
            if locus == "protein-coding gene":
                return BiotypeResolution(symbol, "matched", approved, locus)
            if not locus:
                return BiotypeResolution(symbol, "biotype_unknown", approved, locus or None)
            return BiotypeResolution(symbol, "non_coding_rejected", approved, locus)
        if key in self._alias_to_symbol:
            approved = self._alias_to_symbol[key]
            locus = self._by_symbol.get(approved.upper(), (approved, ""))[1]
            if locus == "protein-coding gene":
                return BiotypeResolution(symbol, "alias_resolved", approved, locus)
            if not locus:
                return BiotypeResolution(symbol, "biotype_unknown", approved, locus or None)
            return BiotypeResolution(symbol, "non_coding_rejected", approved, locus)
        return BiotypeResolution(symbol, "unmatched", None, None)


def load_hgnc_reference(cache_dir: str | Path = CACHE_DIR) -> HgncReference:
    cache_dir = ensure_directory(cache_dir)
    destination = Path(cache_dir) / "hgnc_complete_set.txt"
    accessed = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    if not destination.exists():
        response = requests.get(HGNC_URL, timeout=REQUEST_TIMEOUT)
        if response.status_code != 200:
            raise StopAndReportError(
                f"Authoritative HGNC biotype reference unavailable: HTTP {response.status_code}",
                reason="biotype-reference-unavailable",
                details={"url": HGNC_URL, "status_code": response.status_code},
            )
        destination.write_bytes(response.content)
        accessed = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    frame = pd.read_csv(destination, sep="\t", dtype=str, low_memory=False)
    required = {"symbol", "locus_group"}
    if not required.issubset(frame.columns):
        raise StopAndReportError(
            "HGNC reference missing required columns",
            reason="biotype-reference-unavailable",
            details={"columns": frame.columns.tolist()},
        )
    return HgncReference(
        frame,
        source_url=HGNC_URL,
        path=str(destination),
        accessed_at_utc=accessed,
        sha256=sha256_path(destination),
    )
