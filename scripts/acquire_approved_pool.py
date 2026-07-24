#!/usr/bin/env python3
"""Acquire disease-agnostic approved-drug universe + ChEMBL/DGIdb targets."""

from __future__ import annotations

import json

from signflow.phase5 import PHASE5_RAW, sha256_path
from signflow.phase5.acquire import (
    assert_pool_guards,
    fetch_chembl_approved_molecules,
    fetch_chembl_mechanisms_for_universe,
    fetch_dgidb_by_names,
    load_approved_drug_universe,
)


def main() -> None:
    raw = PHASE5_RAW
    print("acquiring approved universe...", flush=True)
    rec = fetch_chembl_approved_molecules(raw)
    print(json.dumps(rec, indent=2, default=str), flush=True)
    uni = load_approved_drug_universe(raw)
    print(
        f"universe={len(uni)} gold={int(uni['is_gold'].astype(bool).sum())} "
        f"sha={sha256_path(raw / 'pool' / 'approved_drug_universe.tsv')}",
        flush=True,
    )
    assert_pool_guards(uni)
    print("fetching mechanisms...", flush=True)
    mech = fetch_chembl_mechanisms_for_universe(uni["chembl_id"].astype(str).tolist(), raw)
    print(json.dumps(mech, indent=2, default=str), flush=True)
    print("fetching dgidb...", flush=True)
    dg = fetch_dgidb_by_names(uni["pref_name"].astype(str).tolist(), raw)
    print(json.dumps(dg, indent=2, default=str), flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
