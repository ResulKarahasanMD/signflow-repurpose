#!/usr/bin/env python3
"""Phase 5 §A3 — Profile the top-50 non-gold BENCHMARK drugs.

Gold drugs are absent from the top 50 yet enriched ~4.9x at k=100 and ~5.5x at
k=250. Something systematically outranks them. This profiles the top-50 non-gold
drugs by the signed-flow score: number of mapped targets, mean/max target degree,
and ATC (therapeutic-area) distribution, each compared to the gold set and the full
scored pool. Promiscuous / hub-targeting top drugs support the injection-dominance
account and explain the enrichment-curve shape.
"""

from __future__ import annotations

import argparse
import json
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd

from signflow.phase5 import PHASE3_LATEST, PHASE5_OUT, PHASE5_RAW, ensure_directory
from signflow.phase5.evaluate import build_drug_universe, prepare_drug_targets
from signflow.phase5.score import node_degrees
from signflow.phase5.substrate import load_plan, signed_causal_edges


def _load_atc_map(pool_json: Path) -> dict[str, list[str]]:
    payload = json.loads(pool_json.read_text())
    out: dict[str, list[str]] = {}
    for mol in payload.get("molecules", []):
        cid = str(mol.get("molecule_chembl_id") or mol.get("chembl_id") or "")
        if not cid:
            continue
        out[cid] = [str(a) for a in (mol.get("atc_classifications") or [])]
    return out


_ATC_L1 = {
    "A": "Alimentary tract & metabolism",
    "B": "Blood & blood-forming organs",
    "C": "Cardiovascular system",
    "D": "Dermatologicals",
    "G": "Genito-urinary & sex hormones",
    "H": "Systemic hormonal (excl. sex hormones)",
    "J": "Antiinfectives (systemic)",
    "L": "Antineoplastic & immunomodulating",
    "M": "Musculo-skeletal system",
    "N": "Nervous system",
    "P": "Antiparasitic",
    "R": "Respiratory system",
    "S": "Sensory organs",
    "V": "Various",
}


def _atc_l1(codes: list[str]) -> list[str]:
    return sorted({_ATC_L1.get(c[0], f"Unknown ({c[0]})") for c in codes if c})


def _degree_stats(genes: list[str], degrees: pd.Series) -> dict[str, float]:
    present = [g for g in genes if g in degrees.index]
    if not present:
        return {"n_targets_in_substrate": 0, "mean_target_degree": float("nan"), "max_target_degree": float("nan")}
    vals = degrees.loc[present].to_numpy(dtype=float)
    return {
        "n_targets_in_substrate": int(len(present)),
        "mean_target_degree": float(np.mean(vals)),
        "max_target_degree": float(np.max(vals)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase3-root", default=None)
    parser.add_argument(
        "--signed-ranking",
        default=str(PHASE5_OUT / "20260724T004001Z" / "scores" / "benchmark" / "ranking_signed.tsv"),
    )
    parser.add_argument("--top-k", type=int, default=50)
    args = parser.parse_args()

    plan = load_plan()
    p3 = Path(args.phase3_root or PHASE3_LATEST).resolve()
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    out = ensure_directory(PHASE5_OUT / stamp / "scores" / "benchmark" / "top50_profile")

    ranking = pd.read_csv(args.signed_ranking, sep="\t")
    drug_universe = build_drug_universe(PHASE5_RAW, plan=plan)
    targets = prepare_drug_targets(PHASE5_RAW, drug_universe)

    real_raw = pd.read_csv(p3 / "benchmark" / "graph_edges.tsv", sep="\t")
    signed = signed_causal_edges(real_raw, plan)
    nodes = sorted(set(signed["source"].astype(str)) | set(signed["target"].astype(str)))
    node_set = set(nodes)
    degrees = pd.Series(node_degrees(signed, nodes), index=nodes)

    atc_map = _load_atc_map(PHASE5_RAW / "pool" / "chembl_max_phase4_molecules.json")

    # Total mapped targets per drug (any source) and signed in-substrate targets.
    targets_by_drug = targets.groupby("chembl_id")
    n_targets_total = targets_by_drug["gene_symbol"].nunique().to_dict()
    in_sub_by_drug = {
        cid: sorted(set(g["gene_symbol"].astype(str)) & node_set) for cid, g in targets_by_drug
    }

    ranked = ranking.sort_values(["score", "chembl_id"], ascending=[False, True]).reset_index(drop=True)
    ranked["rank"] = np.arange(1, len(ranked) + 1)

    def _row_profile(row: pd.Series) -> dict[str, object]:
        cid = str(row["chembl_id"])
        genes = in_sub_by_drug.get(cid, [])
        stats = _degree_stats(genes, degrees)
        codes = atc_map.get(cid, [])
        return {
            "rank": int(row["rank"]),
            "chembl_id": cid,
            "pref_name": row.get("pref_name", ""),
            "is_gold": bool(row.get("is_gold", False)),
            "score": float(row["score"]),
            "n_targets_total": int(n_targets_total.get(cid, 0)),
            **stats,
            "atc_codes": ";".join(codes),
            "atc_l1": ";".join(_atc_l1(codes)),
            "n_atc": len(codes),
        }

    top = ranked.head(args.top_k)
    top_profiles = pd.DataFrame([_row_profile(r) for _, r in top.iterrows()])
    top_nongold = top_profiles.loc[~top_profiles["is_gold"]].copy()
    top_nongold.to_csv(out / "top50_non_gold_profile.tsv", sep="\t", index=False)

    # Reference distributions: gold set and full scored pool.
    all_profiles = pd.DataFrame([_row_profile(r) for _, r in ranked.iterrows()])
    gold = all_profiles.loc[all_profiles["is_gold"]]

    def _dist(df: pd.DataFrame, col: str) -> dict[str, float]:
        vals = df[col].to_numpy(dtype=float)
        vals = vals[np.isfinite(vals)]
        if not len(vals):
            return {"n": 0}
        return {
            "n": int(len(vals)),
            "median": float(np.median(vals)),
            "mean": float(np.mean(vals)),
            "p90": float(np.quantile(vals, 0.9)),
            "max": float(np.max(vals)),
        }

    def _atc_counts(df: pd.DataFrame) -> dict[str, int]:
        counts: dict[str, int] = {}
        for cell in df["atc_l1"].fillna(""):
            for area in [a for a in str(cell).split(";") if a]:
                counts[area] = counts.get(area, 0) + 1
        return dict(sorted(counts.items(), key=lambda kv: kv[1], reverse=True))

    summary = {
        "diagnostic": "A3_top50_non_gold_profile",
        "started_at_utc": stamp,
        "phase3_root": str(p3),
        "signed_ranking": args.signed_ranking,
        "top_k": args.top_k,
        "n_top_non_gold": int(len(top_nongold)),
        "n_top_gold": int(top_profiles["is_gold"].sum()),
        "distributions": {
            "top50_non_gold": {
                "n_targets_total": _dist(top_nongold, "n_targets_total"),
                "n_targets_in_substrate": _dist(top_nongold, "n_targets_in_substrate"),
                "mean_target_degree": _dist(top_nongold, "mean_target_degree"),
                "max_target_degree": _dist(top_nongold, "max_target_degree"),
            },
            "gold_set": {
                "n_targets_total": _dist(gold, "n_targets_total"),
                "n_targets_in_substrate": _dist(gold, "n_targets_in_substrate"),
                "mean_target_degree": _dist(gold, "mean_target_degree"),
                "max_target_degree": _dist(gold, "max_target_degree"),
            },
            "full_scored_pool": {
                "n_targets_total": _dist(all_profiles, "n_targets_total"),
                "n_targets_in_substrate": _dist(all_profiles, "n_targets_in_substrate"),
                "mean_target_degree": _dist(all_profiles, "mean_target_degree"),
                "max_target_degree": _dist(all_profiles, "max_target_degree"),
            },
        },
        "atc_l1_counts": {
            "top50_non_gold": _atc_counts(top_nongold),
            "gold_set": _atc_counts(gold),
        },
        "interpretation": (
            "If top-50 non-gold drugs carry more targets and/or higher-degree (hub) targets "
            "than the gold set and pool, the ranking is promiscuity/hub-driven via the drug "
            "injection term, consistent with the A1 network-invariance finding and the shape "
            "of the k=100/250 enrichment curve."
        ),
    }
    (out / "a3_top50_profile.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(json.dumps(summary, indent=2, default=str))
    print(f"[a3] wrote {out / 'a3_top50_profile.json'}", flush=True)


if __name__ == "__main__":
    main()
