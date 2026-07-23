"""Expression matrix parsing and probe-to-gene processing for Phase 2."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from signflow.mapping import GeneRecord, MappingError, validate_gene_record

from . import PROJECT_ROOT, StopAndReportError, ensure_directory, sha256_bytes
from .phenotype import parse_series_matrix_samples

REQUEST_TIMEOUT_SECONDS = 120
PLATFORM_CACHE_ROOT = PROJECT_ROOT / "data" / "raw" / "geo" / "platforms"
AUTHOR_NORMALIZED_NOTE = (
    "Processing uses the author-normalised GEO series matrix values "
    "(not CEL-level RMA). Choice recorded for provenance."
)


def _read_text(path: str | Path) -> str:
    payload = Path(path).read_bytes()
    if str(path).endswith(".gz"):
        return gzip.decompress(payload).decode("utf-8", errors="replace")
    return payload.decode("utf-8", errors="replace")


def detect_agilent_channel_count(series_matrix_path: str | Path) -> int:
    """Detect Agilent channel count from sample metadata."""

    samples = parse_series_matrix_samples(series_matrix_path)
    channel_values = set()
    if "Sample_channel_count" in samples:
        for values in samples["Sample_channel_count"]:
            if isinstance(values, list):
                channel_values.update(str(value).strip() for value in values if str(value).strip())
    if not channel_values:
        raise StopAndReportError(
            "Could not determine Agilent channel count from GEO metadata",
            accession="GSE57560",
            reason="undeterminable-agilent-channel",
        )
    if len(channel_values) != 1:
        raise StopAndReportError(
            "Agilent channel count metadata disagrees across samples",
            accession="GSE57560",
            reason="undeterminable-agilent-channel",
            details={"channel_values": sorted(channel_values)},
        )
    return int(next(iter(channel_values)))


def parse_series_matrix_expression(series_matrix_path: str | Path) -> pd.DataFrame:
    """Load the expression table from a GEO series matrix."""

    text = _read_text(series_matrix_path)
    start_marker = "!series_matrix_table_begin"
    end_marker = "!series_matrix_table_end"
    if start_marker not in text or end_marker not in text:
        raise StopAndReportError(
            "Series matrix is missing an expression table block",
            reason="unavailable-files",
            details={"path": str(series_matrix_path)},
        )
    body = text.split(start_marker, 1)[1].split(end_marker, 1)[0].strip().splitlines()
    reader = csv.reader(body, delimiter="\t", quotechar='"')
    header = next(reader)
    if header[0] != "ID_REF":
        raise StopAndReportError(
            "Expected ID_REF as the first expression-table column",
            reason="unavailable-files",
            details={"header": header[:5]},
        )
    frame = pd.DataFrame(reader, columns=header)
    frame = frame.rename(columns={"ID_REF": "probe_id"})
    for column in frame.columns[1:]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return frame


def _download_text(url: str, destination: Path) -> None:
    response = requests.get(url, timeout=REQUEST_TIMEOUT_SECONDS)
    if response.status_code != 200:
        raise requests.HTTPError(f"{response.status_code} {url}")
    ensure_directory(destination.parent)
    destination.write_bytes(response.content)
    _ = sha256_bytes(response.content)


def _geo_platform_prefix(gpl: str) -> str:
    return f"{gpl[:-3]}nnn"


def load_platform_annotation(gpl: str, cache_root: str | Path = PLATFORM_CACHE_ROOT) -> pd.DataFrame:
    """Download and parse a GEO platform annotation table."""

    cache_root = Path(cache_root)
    platform_root = ensure_directory(cache_root / gpl)
    sources = [
        (
            f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{_geo_platform_prefix(gpl)}/{gpl}/annot/{gpl}.annot.gz",
            platform_root / f"{gpl}.annot.gz",
            True,
        ),
        (
            f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gpl}&targ=self&form=text&view=full",
            platform_root / f"{gpl}.txt",
            False,
        ),
    ]

    text = None
    for url, path, is_gzip in sources:
        try:
            if not path.exists():
                _download_text(url, path)
            payload = path.read_bytes()
            text = gzip.decompress(payload).decode("utf-8", errors="replace") if is_gzip else payload.decode("utf-8", errors="replace")
            break
        except Exception:
            continue

    if text is None:
        raise StopAndReportError(
            f"Unable to fetch platform annotation for {gpl}",
            reason="unavailable-files",
            details={"platform": gpl},
        )

    start_marker = "!platform_table_begin"
    end_marker = "!platform_table_end"
    if start_marker not in text or end_marker not in text:
        raise StopAndReportError(
            f"Platform annotation table not found for {gpl}",
            reason="unavailable-files",
            details={"platform": gpl},
        )

    body = text.split(start_marker, 1)[1].split(end_marker, 1)[0].strip().splitlines()
    reader = csv.reader(body, delimiter="\t", quotechar='"')
    header = next(reader)
    return pd.DataFrame(reader, columns=header)


def _candidate_column(columns: list[str], names: list[str]) -> str | None:
    lowered = {column.lower(): column for column in columns}
    for name in names:
        if name.lower() in lowered:
            return lowered[name.lower()]
    return None


def _normalize_gene_type(value: str | None) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text or text in {"---", "NA", "na", "None"}:
        return None
    return text.lower().replace("-", "_")


def collapse_probes_to_genes(
    expression: pd.DataFrame,
    annotation: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Collapse probes to genes by max-mean; reject non-protein-coding via validate_gene_record."""

    probe_column = _candidate_column(annotation.columns.tolist(), ["ID", "ID_REF", "probe_id"])
    symbol_column = _candidate_column(
        annotation.columns.tolist(),
        ["Gene symbol", "GENE_SYMBOL", "gene_symbol", "Symbol"],
    )
    type_column = _candidate_column(
        annotation.columns.tolist(),
        ["gene_type", "Gene type", "gene_biotype", "Biotype"],
    )
    ensembl_column = _candidate_column(
        annotation.columns.tolist(),
        ["ENSEMBL_ID", "Ensembl ID", "ensembl_id", "Gene ID"],
    )
    if not probe_column or not symbol_column:
        raise StopAndReportError(
            "Platform annotation is missing probe or gene symbol columns",
            reason="unavailable-files",
            details={"columns": annotation.columns.tolist()},
        )

    keep_columns = [probe_column, symbol_column]
    if type_column:
        keep_columns.append(type_column)
    if ensembl_column:
        keep_columns.append(ensembl_column)
    merged = expression.merge(annotation[keep_columns], left_on="probe_id", right_on=probe_column, how="left")
    merged = merged.rename(columns={symbol_column: "gene_symbol"})
    if type_column:
        merged = merged.rename(columns={type_column: "gene_type"})
    else:
        merged["gene_type"] = None
    if ensembl_column:
        merged = merged.rename(columns={ensembl_column: "ensembl_id"})
    else:
        merged["ensembl_id"] = ""

    probes_in = int(len(expression))
    merged["gene_symbol"] = merged["gene_symbol"].fillna("").astype(str).str.split("///").str[0].str.strip()
    merged["gene_type"] = merged["gene_type"].map(_normalize_gene_type)
    unmapped = int(((merged["gene_symbol"] == "") | (merged["gene_symbol"] == "---")).sum())
    merged = merged[(merged["gene_symbol"] != "") & (merged["gene_symbol"] != "---")]

    kept_rows = []
    rejected_noncoding = 0
    for _, row in merged.iterrows():
        biotype = row["gene_type"]
        if biotype is None:
            kept_rows.append(row)
            continue
        try:
            validate_gene_record(
                GeneRecord(
                    symbol=str(row["gene_symbol"]),
                    ensembl_id=str(row.get("ensembl_id") or "UNKNOWN"),
                    biotype=str(biotype),
                )
            )
            kept_rows.append(row)
        except MappingError:
            rejected_noncoding += 1

    if not kept_rows:
        raise StopAndReportError("All probes were filtered during probe-to-gene mapping", reason="unavailable-files")

    filtered = pd.DataFrame(kept_rows)
    value_columns = [column for column in expression.columns if column != "probe_id"]
    filtered["probe_mean"] = filtered[value_columns].mean(axis=1, skipna=True)
    before_collapse = int(len(filtered))
    collapsed = (
        filtered.sort_values(["gene_symbol", "probe_mean"], ascending=[True, False])
        .drop_duplicates(subset=["gene_symbol"], keep="first")
        .set_index("gene_symbol")[value_columns]
        .sort_index()
    )
    many_to_one = before_collapse - int(len(collapsed))
    qc = pd.DataFrame(
        [
            {
                "probes_in": probes_in,
                "genes_out": int(len(collapsed)),
                "unmapped_probes": unmapped,
                "rejected_non_protein_coding": rejected_noncoding,
                "many_to_one_collapses": many_to_one,
                "collapse_rule": "max-mean probe per gene",
                "processing_note": AUTHOR_NORMALIZED_NOTE,
            }
        ]
    )
    return collapsed, qc


def maybe_log2_transform(expression: pd.DataFrame) -> pd.DataFrame:
    """Log2-transform strictly positive matrices that appear unlogged."""

    values = expression.to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size and float(np.nanmax(finite)) > 100:
        return np.log2(expression.clip(lower=1.0))
    return expression


def process_expression_study(
    accession: str,
    platform_id: str,
    series_matrix_path: str | Path,
    selected_sample_ids: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    """Parse, normalize, annotate, and collapse one expression study."""

    processing_meta: dict[str, object] = {
        "accession": accession,
        "platform_id": platform_id,
        "expression_source": "author_normalised_series_matrix",
        "note": AUTHOR_NORMALIZED_NOTE,
    }

    if platform_id == "GPL16699":
        channel_count = detect_agilent_channel_count(series_matrix_path)
        processing_meta["agilent_channel_count"] = channel_count
        if channel_count != 1:
            raise StopAndReportError(
                f"{accession} uses {channel_count}-color Agilent data; only 1-color handling is implemented",
                accession=accession,
                reason="undeterminable-agilent-channel",
                details={"channel_count": channel_count},
            )

    table = parse_series_matrix_expression(series_matrix_path)
    missing_samples = sorted(set(selected_sample_ids) - set(table.columns))
    if missing_samples:
        raise StopAndReportError(
            f"Expression matrix is missing selected samples for {accession}",
            accession=accession,
            reason="unavailable-files",
            details={"missing_samples": missing_samples},
        )
    table = table[["probe_id", *selected_sample_ids]].copy()
    table[selected_sample_ids] = maybe_log2_transform(table[selected_sample_ids])
    annotation = load_platform_annotation(platform_id)
    collapsed, qc = collapse_probes_to_genes(table, annotation)
    return collapsed, qc, processing_meta
