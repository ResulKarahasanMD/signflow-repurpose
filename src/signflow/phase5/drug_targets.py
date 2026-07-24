"""Merge ChEMBL + DGIdb drug-target rows with normalized actions."""
from __future__ import annotations
import pandas as pd
from .actions import UNKNOWN, action_sign, is_signed_action, load_action_sign_map, normalize_action_token

def _attach_actions(frame: pd.DataFrame) -> pd.DataFrame:
    sign_map = load_action_sign_map()
    out = frame.copy()
    out["canonical_action"] = out["raw_action"].map(normalize_action_token)
    out["action_sign"] = out["raw_action"].map(lambda x: action_sign(x, sign_map=sign_map))
    out["signed_action"] = out["canonical_action"].map(lambda a: a if a != UNKNOWN else None)
    out["is_signed"] = out["raw_action"].map(lambda x: is_signed_action(x, sign_map=sign_map))
    return out

def build_drug_target_table(chembl_frame, dgidb_frame, drug_universe):
    uni = drug_universe.copy()
    uni["chembl_id"] = uni["chembl_id"].astype(str)
    uni["pref_name_upper"] = uni.get("pref_name", "").astype(str).str.upper()
    chembl = chembl_frame.merge(uni[["chembl_id","pref_name","pref_name_upper","is_gold","canonical_group"]], on="chembl_id", how="inner") if len(chembl_frame) else pd.DataFrame()
    dgidb = dgidb_frame.copy()
    if len(dgidb):
        dgidb["drug_name_upper"] = dgidb["drug_name"].astype(str).str.upper()
        dgidb = dgidb.merge(uni[["chembl_id","pref_name","pref_name_upper","is_gold","canonical_group"]], left_on="drug_name_upper", right_on="pref_name_upper", how="inner")
    else:
        dgidb = pd.DataFrame()
    parts = [p for p in (chembl, dgidb) if len(p)]
    if not parts:
        return pd.DataFrame(columns=["chembl_id","pref_name","gene_symbol","raw_action","source","canonical_action","action_sign","signed_action","is_signed","is_gold","canonical_group"])
    merged = _attach_actions(pd.concat(parts, ignore_index=True))
    merged["gene_symbol"] = merged["gene_symbol"].astype(str).str.strip()
    merged = merged.loc[merged["gene_symbol"].str.len() > 0]
    return merged.drop_duplicates(subset=["chembl_id","gene_symbol","source","raw_action"], keep="first").reset_index(drop=True)

def drug_level_target_summary(targets, substrate_nodes):
    rows = []
    for chembl_id, group in targets.groupby("chembl_id"):
        genes = set(group["gene_symbol"].astype(str))
        in_sub = genes & substrate_nodes
        signed = group.loc[group["is_signed"]]
        rows.append({"chembl_id": chembl_id, "pref_name": group["pref_name"].iloc[0] if "pref_name" in group else "",
            "canonical_group": group["canonical_group"].iloc[0] if "canonical_group" in group.columns else "",
            "is_gold": bool(group["is_gold"].iloc[0]) if "is_gold" in group.columns else False,
            "n_targets_total": len(genes), "n_targets_in_substrate": len(in_sub),
            "n_targets_signed": int(signed["gene_symbol"].nunique()),
            "n_targets_signed_in_substrate": len(set(signed["gene_symbol"].astype(str)) & substrate_nodes),
            "targets_in_substrate": sorted(in_sub)})
    return pd.DataFrame(rows)
