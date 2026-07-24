#!/usr/bin/env python3
"""Phase 6 mandatory TARGET defence — GSE238208 inflammation-specificity.

Phase 1.5 §0.2 mandates this as the primary biological check on the TARGET seed
module (there are no TARGET disease-gene labels). GSE238208 contrasts Hunner-type
interstitial cystitis (IC) against BCG cystitis — a different inflammatory bladder
condition. Question: does the frozen TARGET seed module separate Hunner-type IC from
BCG cystitis (Hunner-specific inflammation), or does it merely mark generic bladder
inflammation (no separation)?

The seed module was derived from IC-vs-non-ulcer contrasts in GSE11783 + GSE57560,
so GSE238208 (IC vs BCG) is an out-of-sample specificity test. We build a directional
per-sample module score (signature-signed mean of z-scored log2 FPKM over seed genes),
measure how well it separates IC from BCG (AUROC + Mann-Whitney), and calibrate that
AUROC against random gene sets of matched size (permutation null). Both the primary and
sensitivity-concordant seed sets are reported. This is descriptive characterisation of
the frozen module (H7 non-informative); no ranked candidate list is produced.
"""

from __future__ import annotations

import argparse
import gzip
import json
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score

from signflow.phase5 import PHASE3_LATEST, PHASE5_OUT, ensure_directory

FPKM = Path("data/raw/geo/GSE238208/suppl/GSE238208_gene_expression_yoshiyukiakiyama.txt.gz")


def _load_expression() -> tuple[pd.DataFrame, np.ndarray]:
    """Return (log2 FPKM genes x samples, group array) with IC=1, BCG=0.

    Column layout (from the supplementary header): 1-25 Hunner lesion, 26-50
    non-Hunner lesion (both Hunner-type IC = case), 51-63 BCG cystitis (control).
    """

    with gzip.open(FPKM, "rt") as f:
        f.readline()  # group header row
        header = f.readline().rstrip("\n").split("\t")
    fpkm_cols = [c for c in header if c.endswith(".FPKM")]
    df = pd.read_csv(FPKM, sep="\t", skiprows=1, index_col=0, usecols=["Name", *fpkm_cols])
    df = df[fpkm_cols]
    # 63 samples: first 50 = Hunner-type IC (case), last 13 = BCG cystitis (control).
    if df.shape[1] != 63:
        raise SystemExit(f"Unexpected FPKM sample count: {df.shape[1]} (expected 63)")
    group = np.array([1] * 50 + [0] * 13, dtype=int)  # 1 = IC, 0 = BCG
    df = df[~df.index.duplicated(keep="first")]
    logx = np.log2(df.astype(float) + 1.0)
    return logx, group


def _zscore_genes(logx: pd.DataFrame) -> pd.DataFrame:
    mu = logx.mean(axis=1)
    sd = logx.std(axis=1).replace(0, np.nan)
    z = logx.sub(mu, axis=0).div(sd, axis=0)
    return z.dropna(how="all")


def _module_score(z: pd.DataFrame, genes: list[str], directions: dict[str, float]) -> tuple[np.ndarray, list[str]]:
    present = [g for g in genes if g in z.index]
    if not present:
        return np.zeros(z.shape[1]), []
    signs = np.array([np.sign(directions.get(g, 1.0)) or 1.0 for g in present], dtype=float)
    sub = z.loc[present].to_numpy(dtype=float)
    score = np.nanmean(sub * signs[:, None], axis=0)
    return score, present


def _evaluate_seed_set(
    *,
    name: str,
    seed_frame: pd.DataFrame,
    z: pd.DataFrame,
    group: np.ndarray,
    rng: np.random.Generator,
    n_perm: int,
) -> dict[str, object]:
    genes = seed_frame["gene_symbol"].astype(str).tolist()
    directions = dict(zip(seed_frame["gene_symbol"].astype(str), seed_frame["pooled_effect"].astype(float)))
    score, present = _module_score(z, genes, directions)
    n_present = len(present)
    if n_present < 3:
        return {"seed_set": name, "status": "insufficient_gene_coverage", "n_seed_genes": len(genes), "n_present": n_present}

    auroc = float(roc_auc_score(group, score))
    mw = mannwhitneyu(score[group == 1], score[group == 0], alternative="two-sided")

    # Per-seed-gene direction agreement between IC-vs-BCG and the TARGET signature.
    ic = z.loc[present, :].to_numpy(dtype=float)[:, group == 1].mean(axis=1)
    bcg = z.loc[present, :].to_numpy(dtype=float)[:, group == 0].mean(axis=1)
    obs_dir = np.sign(ic - bcg)
    sig_dir = np.array([np.sign(directions.get(g, 1.0)) or 1.0 for g in present])
    frac_agree = float(np.mean(obs_dir == sig_dir))

    # Random-gene-set null: does a size-matched random module separate IC from BCG as well?
    all_genes = z.index.to_numpy()
    null_auroc = np.empty(n_perm, dtype=float)
    for b in range(n_perm):
        pick = rng.choice(len(all_genes), size=n_present, replace=False)
        rand_genes = list(all_genes[pick])
        # Random directions matched to the seed sign multiset (removes directional prior).
        rand_dir = {g: s for g, s in zip(rand_genes, rng.permutation(sig_dir))}
        rscore, _ = _module_score(z, rand_genes, rand_dir)
        # AUROC of |centered| separation is orientation-agnostic; take max(auroc, 1-auroc)
        a = roc_auc_score(group, rscore)
        null_auroc[b] = max(a, 1.0 - a)
    directed_auroc = max(auroc, 1.0 - auroc)
    p_emp = float((np.sum(null_auroc >= directed_auroc) + 1) / (n_perm + 1))

    separates = bool(directed_auroc > np.quantile(null_auroc, 0.95))
    return {
        "seed_set": name,
        "status": "ok",
        "n_seed_genes": len(genes),
        "n_present_in_gse238208": n_present,
        "n_ic": int((group == 1).sum()),
        "n_bcg": int((group == 0).sum()),
        "module_score_auroc_ic_vs_bcg": auroc,
        "orientation_agnostic_auroc": directed_auroc,
        "mannwhitney_p": float(mw.pvalue),
        "seed_gene_direction_agreement_ic_vs_bcg": frac_agree,
        "random_gene_set_null_auroc": {
            "n_perm": n_perm,
            "mean": float(np.mean(null_auroc)),
            "q95": float(np.quantile(null_auroc, 0.95)),
            "max": float(np.max(null_auroc)),
        },
        "empirical_p_auroc_ge_random": p_emp,
        "module_separates_ic_from_bcg": separates,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase3-root", default=None)
    parser.add_argument("--n-perm", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    p3 = Path(args.phase3_root or PHASE3_LATEST).resolve()
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    out = ensure_directory(PHASE5_OUT / stamp / "scores" / "target" / "gse238208_specificity")

    logx, group = _load_expression()
    z = _zscore_genes(logx)
    print(f"[gse238208] expression genes={z.shape[0]} samples={z.shape[1]} IC={int((group==1).sum())} BCG={int((group==0).sum())}", flush=True)

    seeds_dir = p3 / "target" / "seeds"
    primary = pd.read_csv(seeds_dir / "target_seeds_primary.tsv", sep="\t")
    sensitivity = pd.read_csv(seeds_dir / "target_seeds_sensitivity_concordant.tsv", sep="\t")

    rng = np.random.default_rng(args.seed)
    results = [
        _evaluate_seed_set(name="primary", seed_frame=primary, z=z, group=group, rng=rng, n_perm=args.n_perm),
        _evaluate_seed_set(name="sensitivity_concordant", seed_frame=sensitivity, z=z, group=group, rng=rng, n_perm=args.n_perm),
    ]

    def _verdict(r: dict) -> str:
        if r.get("status") != "ok":
            return r.get("status", "unknown")
        return "hunner_specific" if r["module_separates_ic_from_bcg"] else "generic_bladder_inflammation"

    summary = {
        "diagnostic": "gse238208_inflammation_specificity",
        "started_at_utc": stamp,
        "phase3_root": str(p3),
        "contrast": "Hunner-type interstitial cystitis (n=50) vs BCG cystitis (n=13)",
        "h7_status": "non_informative_descriptive_only",
        "seed_sets": {r["seed_set"]: r for r in results},
        "verdict": {r["seed_set"]: _verdict(r) for r in results},
        "interpretation": (
            "The TARGET seed module was derived from IC-vs-non-ulcer contrasts; GSE238208 "
            "is an out-of-sample test against a different inflammatory bladder condition "
            "(BCG cystitis). If the directional module score separates IC from BCG beyond "
            "size-matched random gene sets, the module captures Hunner-type-IC-specific "
            "inflammation rather than generic bladder inflammation. Descriptive only "
            "(H7 non-informative); no ranked candidate list is derived from this test."
        ),
    }
    (out / "gse238208_specificity.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(json.dumps(summary, indent=2, default=str))
    print(f"[gse238208] wrote {out / 'gse238208_specificity.json'}", flush=True)


if __name__ == "__main__":
    main()
