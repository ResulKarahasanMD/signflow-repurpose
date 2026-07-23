"""Phase 3 orchestrator: networks → mapping → features → labels/seeds → splits → controls → PyG."""

from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import yaml

from signflow.phase2 import PHASE2_META_ROOT, PHASE2_OUTPUT_ROOT

from . import FROZEN_PLAN, PHASE3_OUT, StopAndReportError, ensure_directory
from .acquire import acquire_networks
from .controls import emit_control_graphs
from .edges import build_causal_edge_table, build_string_edge_table
from .features import FEATURE_COLUMNS, build_feature_table
from .id_map import build_mapper
from .labels import write_benchmark_labels
from .pyg_build import assemble_pyg
from .splits import emit_splits


def _network_summary(nodes: list[str], edges: pd.DataFrame) -> dict[str, object]:
    g = nx.DiGraph()
    g.add_nodes_from(nodes)
    for _, row in edges.iterrows():
        if row["source"] in g and row["target"] in g:
            g.add_edge(row["source"], row["target"])
    undirected = g.to_undirected()
    components = list(nx.connected_components(undirected))
    degrees = pd.Series({n: g.degree(n) for n in g.nodes()})
    return {
        "n_nodes": g.number_of_nodes(),
        "n_edges": g.number_of_edges(),
        "sign_balance": edges["sign"].value_counts().to_dict() if len(edges) else {},
        "relation_type_distribution": edges["relation_type"].value_counts().to_dict() if len(edges) else {},
        "degree_mean": float(degrees.mean()) if len(degrees) else 0.0,
        "degree_median": float(degrees.median()) if len(degrees) else 0.0,
        "degree_max": int(degrees.max()) if len(degrees) else 0,
        "n_connected_components": len(components),
        "largest_component_size": max((len(c) for c in components), default=0),
    }


def _expressed_genes(track: str) -> tuple[set[str], dict[str, set[str]], dict[str, object]]:
    """Return genes measurable in the track.

    Primary measurable set = meta intersection (present in both cohorts' processed
    tables / signature_meta), NOT significance-filtered DE. Per-study membership is
    retained for LOSO splits and QC.
    """

    if track == "target":
        studies = {"GSE11783": PHASE2_OUTPUT_ROOT / "GSE11783", "GSE57560": PHASE2_OUTPUT_ROOT / "GSE57560"}
        meta_path = PHASE2_META_ROOT / "signature_meta.tsv"
    else:
        studies = {"GSE75214": PHASE2_OUTPUT_ROOT / "GSE75214", "GSE87466": PHASE2_OUTPUT_ROOT / "GSE87466"}
        meta_path = PHASE2_META_ROOT / "signature_meta_benchmark.tsv"
    membership: dict[str, set[str]] = {}
    mapping_notes: dict[str, object] = {}
    per_study: list[set[str]] = []
    for name, root in studies.items():
        de = pd.read_csv(root / "de_limma.tsv", sep="\t")
        genes = set(de["gene_symbol"].astype(str))
        membership[name] = genes
        per_study.append(genes)
        qc_path = root / "mapping_qc.tsv"
        if qc_path.exists():
            qc = pd.read_csv(qc_path, sep="\t")
            if "status" in qc.columns:
                mapping_notes[name] = qc["status"].value_counts().to_dict()
    meta_genes = set(pd.read_csv(meta_path, sep="\t")["gene_symbol"].astype(str))
    # Measurable for pooled features = meta intersection (equals ∩ of study DE tables)
    measurable = meta_genes
    mapping_notes["n_union_expressed"] = int(len(set.union(*per_study)))
    mapping_notes["n_intersection_expressed"] = int(len(set.intersection(*per_study)))
    mapping_notes["n_meta_genes"] = int(len(meta_genes))
    mapping_notes["phase25_unmatched_asymmetry"] = {
        "GSE57560_unmatched": 16698,
        "GSE11783_unmatched": 2910,
        "note": "Asymmetry propagates into measurable gene universe before network intersection",
    }
    return measurable, membership, mapping_notes


def run_phase3(track: str = "both") -> dict[str, object]:
    plan = yaml.safe_load(FROZEN_PLAN.read_text())["plan"]
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    out_root = ensure_directory(PHASE3_OUT / stamp)
    manifest: dict[str, object] = {"started_at_utc": stamp, "track": track, "drug_data_touched": False}

    # B1 acquire
    acquired = acquire_networks()
    manifest["network_acquire"] = acquired["manifest"]
    (out_root / "network_manifest.json").write_text(json.dumps(acquired, indent=2, default=str) + "\n")

    causal, causal_summary = build_causal_edge_table(acquired["paths"]["omnipath"], acquired["paths"]["signor"])
    string_edges, string_summary = build_string_edge_table(
        acquired["paths"]["string_physical"],
        acquired["paths"]["string_info"],
        score_min=int(plan["string_score_min"]),
    )
    causal.to_csv(out_root / "graph_edges_causal_full.tsv", sep="\t", index=False)
    string_edges.to_csv(out_root / "graph_edges_string_aux_full.tsv", sep="\t", index=False)

    mapper = build_mapper()
    tracks = ["target", "benchmark"] if track == "both" else [track]
    results: dict[str, object] = {"causal_summary": causal_summary, "string_summary": string_summary}

    for tr in tracks:
        tr_dir = ensure_directory(out_root / tr)
        expressed, membership, map_notes = _expressed_genes(tr)
        # Map expressed symbols through HGNC
        map_table, map_counts = mapper.map_symbols(sorted(expressed))
        map_table.to_csv(tr_dir / "mapping_qc.tsv", sep="\t", index=False)
        protein_ok = map_table.loc[map_table["status"].isin(["matched", "alias_resolved"]), "approved_symbol"]
        approved = set(protein_ok.dropna().astype(str))

        # Network node symbols (from causal edges) mapped via HGNC
        net_symbols = sorted(set(causal["source"]) | set(causal["target"]))
        net_map, net_counts = mapper.map_symbols(net_symbols)
        net_map.to_csv(tr_dir / "network_symbol_mapping_qc.tsv", sep="\t", index=False)
        net_approved = set(
            net_map.loc[net_map["status"].isin(["matched", "alias_resolved"]), "approved_symbol"].dropna().astype(str)
        )

        # Remap causal edges to approved symbols
        sym_map = dict(zip(net_map["input_symbol"], net_map["approved_symbol"], strict=False))
        edges = causal.copy()
        edges["source"] = edges["source"].map(lambda s: sym_map.get(s, s))
        edges["target"] = edges["target"].map(lambda s: sym_map.get(s, s))
        edges = edges[edges["source"].isin(net_approved) & edges["target"].isin(net_approved)]

        string_tr = string_edges.copy()
        string_tr["source"] = string_tr["source"].map(lambda s: sym_map.get(s, s) if s in sym_map else s)
        string_tr["target"] = string_tr["target"].map(lambda s: sym_map.get(s, s) if s in sym_map else s)

        # Node universe = expressed ∩ network (NOT DE-filtered)
        nodes = sorted(approved & net_approved)
        if len(nodes) < 500:
            raise StopAndReportError(
                f"Node universe collapsed to {len(nodes)} after correct mapping",
                reason="node-universe-collapse",
                details={"n_nodes": len(nodes), "track": tr, "n_expressed_approved": len(approved), "n_network_approved": len(net_approved)},
            )

        edges = edges[edges["source"].isin(nodes) & edges["target"].isin(nodes)].copy()
        string_tr = string_tr[string_tr["source"].isin(nodes) & string_tr["target"].isin(nodes)].copy()

        edges.to_csv(tr_dir / "graph_edges.tsv", sep="\t", index=False)
        string_tr.to_csv(tr_dir / "graph_edges_string_aux.tsv", sep="\t", index=False)

        # Degrees from causal graph
        deg = pd.concat([edges["source"], edges["target"]]).value_counts()
        degrees = deg.reindex(nodes).fillna(0).astype(float)
        degrees.index = nodes

        nodes_df = pd.DataFrame(
            {
                "gene_symbol": nodes,
                "degree": [float(degrees.loc[g]) for g in nodes],
            }
        )
        # attach ids
        id_lookup = map_table.drop_duplicates("approved_symbol").set_index("approved_symbol")
        nodes_df["ensembl_gene_id"] = [
            id_lookup.loc[g, "ensembl_gene_id"] if g in id_lookup.index else None for g in nodes
        ]
        nodes_df["uniprot_ids"] = [
            id_lookup.loc[g, "uniprot_ids"] if g in id_lookup.index else "" for g in nodes
        ]
        nodes_df.to_csv(tr_dir / "graph_nodes.tsv", sep="\t", index=False)

        summary = _network_summary(nodes, edges)
        summary["mapping_qc_counts"] = map_counts
        summary["network_mapping_counts"] = net_counts
        summary["phase25_mapping_notes"] = map_notes
        summary["n_string_aux_edges"] = int(len(string_tr))
        (tr_dir / "network_summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")

        features, missingness = build_feature_table(nodes, track=tr, degrees=degrees)
        # Align to nodes that survived feature build
        nodes = features["gene_symbol"].tolist()
        degrees = degrees.reindex(nodes).fillna(0)
        edges = edges[edges["source"].isin(nodes) & edges["target"].isin(nodes)]
        string_tr = string_tr[string_tr["source"].isin(nodes) & string_tr["target"].isin(nodes)]
        features.to_csv(tr_dir / "node_features.tsv", sep="\t", index=False)
        pd.DataFrame([missingness]).to_csv(tr_dir / "feature_missingness.tsv", sep="\t", index=False)

        labels_series = None
        if tr == "benchmark":
            label_meta = write_benchmark_labels(
                tr_dir / "labels",
                nodes,
                degrees,
                seed=int(plan["phase3_control_seed"]),
                tolerance=float(plan["degree_match_tolerance"]),
            )
            pos = pd.read_csv(tr_dir / "labels" / "uc_positives_in_graph.tsv", sep="\t")
            neg = pd.read_csv(tr_dir / "labels" / "degree_matched_controls.tsv", sep="\t")
            y = pd.Series(-1, index=nodes, dtype=int)
            y.loc[y.index.isin(pos["gene_symbol"])] = 1
            y.loc[y.index.isin(neg["gene_symbol"])] = 0
            # For splits use labelled subset only
            labelled = y[y >= 0]
            splits = emit_splits(
                tr_dir / "splits",
                labelled,
                degrees.reindex(labelled.index).fillna(0),
                nodes,
                membership,
                seed=int(plan["phase3_control_seed"]),
            )
            labels_series = y
            (tr_dir / "label_meta.json").write_text(json.dumps(label_meta, indent=2) + "\n")
        else:
            # TARGET: seeds only, no labels
            seed_src = PHASE2_META_ROOT / "diagnostics_phase26"
            seed_dir = ensure_directory(tr_dir / "seeds")
            for name in ("target_seeds_primary.tsv", "target_seeds_sensitivity_concordant.tsv"):
                src = seed_src / name
                if src.exists():
                    pd.read_csv(src, sep="	").to_csv(seed_dir / name, sep="	", index=False)
            if (seed_src / "target_seeds_primary.tsv").exists():
                p = set(pd.read_csv(seed_src / "target_seeds_primary.tsv", sep="\t")["gene_symbol"])
                s = set(pd.read_csv(seed_src / "target_seeds_sensitivity_concordant.tsv", sep="\t")["gene_symbol"])
                overlap = {
                    "n_primary": len(p),
                    "n_sensitivity": len(s),
                    "n_overlap": len(p & s),
                    "overlap_genes": sorted(p & s),
                }
                (tr_dir / "seeds" / "seed_overlap.json").write_text(json.dumps(overlap, indent=2) + "\n")
            # dummy labels of -1 for split artifact structure not required; skip supervised splits
            splits = {"note": "TARGET has no labels; supervised splits not emitted"}
            (ensure_directory(tr_dir / "splits") / "splits.json").write_text(json.dumps(splits, indent=2) + "\n")

        controls = emit_control_graphs(edges, tr_dir / "controls")
        pyg_stats = assemble_pyg(
            nodes,
            features,
            edges,
            string_tr,
            labels_series,
            None,
            map_table,
            tr_dir / "pyg",
        )
        results[tr] = {
            "n_nodes": len(nodes),
            "network_summary": summary,
            "feature_missingness": missingness,
            "controls": controls,
            "pyg": pyg_stats,
            "output_dir": str(tr_dir),
        }

    results["drug_data_touched"] = False
    results["conflict_rule"] = plan["edge_sign_conflict_rule"]
    results["output_root"] = str(out_root)
    (out_root / "phase3_result.json").write_text(json.dumps(results, indent=2, default=str) + "\n")
    # symlink/latest pointer
    latest = PHASE3_OUT / "latest"
    if latest.is_symlink() or latest.exists():
        latest.unlink()
    latest.symlink_to(out_root.name)
    return results


if __name__ == "__main__":
    print(json.dumps(run_phase3("both"), indent=2, default=str))
