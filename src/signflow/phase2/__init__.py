"""Phase 2 expression-only workflow for SignFlow-Repurpose."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[3]
CONFIG_ROOT = PROJECT_ROOT / "configs"
RAW_GEO_ROOT = PROJECT_ROOT / "data" / "raw" / "geo"
PHASE2_OUTPUT_ROOT = PROJECT_ROOT / "data" / "processed" / "phase2"
PHASE2_META_ROOT = PHASE2_OUTPUT_ROOT / "meta"
PHENOTYPE_MAP_PATH = CONFIG_ROOT / "phase2_phenotype_maps.yaml"

MANDATED_QC_HEADER = (
    "DESCRIPTIVE QC ONLY. Phase 2 uses expression data only and does not import disease-gene labels, "
    "drug priors, or network priors. Activity outputs from CollecTRI/PROGENy are retained with an "
    "explicit domain-mismatch limitation and are not used to overrule expression-derived results."
)


class StopAndReportError(RuntimeError):
    """Raised when a preregistered Phase 2 guardrail requires a hard stop."""

    def __init__(
        self,
        message: str,
        *,
        accession: str | None = None,
        reason: str | None = None,
        details: dict[str, Any] | None = None,
    ) -> None:
        super().__init__(message)
        self.accession = accession
        self.reason = reason or "stop-and-report"
        self.details = details or {}

    def to_record(self) -> dict[str, Any]:
        return {
            "status": "stopped",
            "reason": self.reason,
            "accession": self.accession,
            "message": str(self),
            "details": self.details,
        }


@dataclass(frozen=True)
class DownloadRecord:
    """One downloaded or reused remote artifact."""

    accession: str
    url: str
    path: str
    sha256: str
    size_bytes: int
    status: str = "downloaded"
    accessed_at_utc: str = ""


@dataclass
class Phase2Result:
    """Serializable summary returned by the Phase 2 runner."""

    track: str
    files_created: list[str] = field(default_factory=list)
    per_accession_counts: dict[str, dict[str, int]] = field(default_factory=dict)
    concordance: dict[str, Any] = field(default_factory=dict)
    stop_conditions: list[dict[str, Any]] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "track": self.track,
            "files_created": self.files_created,
            "per_accession_counts": self.per_accession_counts,
            "concordance": self.concordance,
            "stop_conditions": self.stop_conditions,
            "notes": self.notes,
        }


def ensure_directory(path: str | Path) -> Path:
    """Create a directory tree if needed and return the path."""

    output = Path(path)
    output.mkdir(parents=True, exist_ok=True)
    return output


def sha256_bytes(payload: bytes) -> str:
    """Return the SHA256 of in-memory bytes."""

    return hashlib.sha256(payload).hexdigest()


def sha256_path(path: str | Path) -> str:
    """Return the SHA256 of a local file."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()
