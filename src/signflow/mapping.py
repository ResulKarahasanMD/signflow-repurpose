"""Gene identifier and biotype mapping safeguards."""

from __future__ import annotations

from dataclasses import dataclass


class MappingError(ValueError):
    """Raised when a gene mapping violates frozen plan rules."""


def strip_ensembl_version(ensembl_id: str) -> str:
    return ensembl_id.split(".", 1)[0]


def normalize_biotype(biotype: str | None) -> str:
    return (biotype or "").strip().lower().replace("-", "_")


def is_protein_coding(biotype: str | None) -> bool:
    return normalize_biotype(biotype) == "protein_coding"


@dataclass(frozen=True)
class GeneRecord:
    symbol: str
    ensembl_id: str
    biotype: str


def validate_gene_record(record: GeneRecord) -> GeneRecord:
    clean_id = strip_ensembl_version(record.ensembl_id)
    if not is_protein_coding(record.biotype):
        raise MappingError(
            f"Rejected non-protein-coding gene {record.symbol} ({record.biotype}); "
            "frozen plan forbids labels or nodes such as UOX."
        )
    return GeneRecord(symbol=record.symbol, ensembl_id=clean_id, biotype=normalize_biotype(record.biotype))
