"""HGNC/Ensembl reference mapping for Phase 3 node identifiers."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from signflow.phase2.biotype_reference import HgncReference, load_hgnc_reference

from . import StopAndReportError


class IdentifierMapper:
    """Map symbols ↔ UniProt ↔ Ensembl via HGNC reference only (never resource annotation)."""

    def __init__(self, ref: HgncReference):
        self.ref = ref
        frame = ref.frame
        if "ensembl_gene_id" not in frame.columns or "uniprot_ids" not in frame.columns:
            raise StopAndReportError(
                "HGNC/Ensembl reference unavailable (missing ensembl_gene_id/uniprot_ids)",
                reason="hgnc-reference-unavailable",
                details={"columns": frame.columns.tolist()},
            )
        self._symbol_to_ids: dict[str, dict[str, object]] = {}
        self._uniprot_to_symbol: dict[str, str] = {}
        self._ensembl_to_symbol: dict[str, str] = {}
        for _, row in frame.iterrows():
            symbol = str(row["symbol"]).strip()
            locus = str(row.get("locus_group") or "").strip().lower()
            ensembl = str(row.get("ensembl_gene_id") or "").strip()
            uniprots = [u.strip() for u in str(row.get("uniprot_ids") or "").replace(",", "|").split("|") if u.strip()]
            self._symbol_to_ids[symbol.upper()] = {
                "approved_symbol": symbol,
                "ensembl_gene_id": ensembl or None,
                "uniprot_ids": uniprots,
                "locus_group": locus,
                "n_uniprot": len(uniprots),
            }
            if ensembl:
                self._ensembl_to_symbol[ensembl] = symbol
            for u in uniprots:
                self._uniprot_to_symbol[u] = symbol

    def map_symbols(self, symbols: list[str]) -> tuple[pd.DataFrame, dict[str, int]]:
        rows = []
        counts = {
            "input": 0,
            "matched": 0,
            "alias_resolved": 0,
            "one_to_many_uniprot": 0,
            "unmatched": 0,
            "non_coding_rejected": 0,
        }
        seen = set()
        for raw in symbols:
            if raw is None:
                continue
            symbol = str(raw).strip()
            if not symbol or symbol in seen:
                continue
            seen.add(symbol)
            counts["input"] += 1
            res = self.ref.resolve(symbol)
            if res.status == "unmatched":
                counts["unmatched"] += 1
                rows.append(
                    {
                        "input_symbol": symbol,
                        "status": "unmatched",
                        "approved_symbol": None,
                        "ensembl_gene_id": None,
                        "uniprot_ids": "",
                        "n_uniprot": 0,
                        "locus_group": None,
                    }
                )
                continue
            if res.status == "non_coding_rejected":
                counts["non_coding_rejected"] += 1
                rows.append(
                    {
                        "input_symbol": symbol,
                        "status": "non_coding_rejected",
                        "approved_symbol": res.approved_symbol,
                        "ensembl_gene_id": None,
                        "uniprot_ids": "",
                        "n_uniprot": 0,
                        "locus_group": res.locus_group,
                    }
                )
                continue
            info = self._symbol_to_ids.get((res.approved_symbol or "").upper(), {})
            n_u = int(info.get("n_uniprot") or 0)
            if n_u > 1:
                counts["one_to_many_uniprot"] += 1
            if res.status == "alias_resolved":
                counts["alias_resolved"] += 1
            else:
                counts["matched"] += 1
            rows.append(
                {
                    "input_symbol": symbol,
                    "status": res.status,
                    "approved_symbol": res.approved_symbol,
                    "ensembl_gene_id": info.get("ensembl_gene_id"),
                    "uniprot_ids": "|".join(info.get("uniprot_ids") or []),
                    "n_uniprot": n_u,
                    "locus_group": info.get("locus_group"),
                }
            )
        return pd.DataFrame(rows), counts


def build_mapper() -> IdentifierMapper:
    return IdentifierMapper(load_hgnc_reference())
