"""Build frozen primary / sensitivity proximity substrates."""

from __future__ import annotations

from pathlib import Path

import networkx as nx
import pandas as pd
import yaml

from . import FROZEN_PLAN, StopAndReportError


def load_plan() -> dict[str, object]:
    return yaml.safe_load(FROZEN_PLAN.read_text())["plan"]


def signed_causal_edges(edges: pd.DataFrame, plan: dict[str, object] | None = None) -> pd.DataFrame:
    """Keep only stimulation/inhibition edges for signed propagation."""

    plan = plan or load_plan()
    if plan["phase5_ambiguous_edge_policy"] != "exclude_from_signed_propagation":
        raise StopAndReportError("Unexpected ambiguous edge policy", reason="policy-drift")
    if plan["phase5_unsigned_directed_edge_policy"] != "exclude_from_signed_propagation":
        raise StopAndReportError("Unexpected unsigned directed policy", reason="policy-drift")
    keep = edges["relation_type"].isin(["stimulation", "inhibition"]) & edges["sign"].isin([1, -1])
    return edges.loc[keep].copy()


def unsigned_full_causal_edges(edges: pd.DataFrame, plan: dict[str, object] | None = None) -> pd.DataFrame:
    """Directed causal edges for unsigned-full: stim + inhib + unsigned_directed.

    Ambiguous edges remain excluded (conflicting signs ≠ merely unsigned).
    Used only with ignore_sign=True so zero-signed unsigned_directed edges
    contribute unit weight connectivity without inventing a sign.
    """

    plan = plan or load_plan()
    if plan["phase5_ambiguous_edge_policy"] != "exclude_from_signed_propagation":
        raise StopAndReportError("Unexpected ambiguous edge policy", reason="policy-drift")
    keep = edges["relation_type"].isin(["stimulation", "inhibition", "unsigned_directed"])
    return edges.loc[keep].copy()


def substrate_nodes(edges: pd.DataFrame) -> set[str]:
    if edges.empty:
        return set()
    return set(edges["source"].astype(str)) | set(edges["target"].astype(str))


def connectivity_report(
    edges: pd.DataFrame,
    *,
    all_nodes: set[str] | None = None,
    track: str = "signed_propagation",
) -> dict[str, object]:
    """Components / largest component / isolates after exclude_from_signed_propagation."""

    nodes = set(all_nodes) if all_nodes is not None else substrate_nodes(edges)
    g = nx.DiGraph()
    g.add_nodes_from(nodes)
    if len(edges):
        for _, row in edges.iterrows():
            g.add_edge(str(row["source"]), str(row["target"]))
    # Weakly connected components treat direction as adjacency for reachability structure
    undirected = g.to_undirected()
    components = list(nx.connected_components(undirected)) if undirected.number_of_nodes() else []
    sizes = sorted((len(c) for c in components), reverse=True)
    isolates = [n for n in nodes if undirected.degree(n) == 0]
    n_excluded_note = (
        "exclude_from_signed_propagation removes ambiguous + unsigned_directed edges "
        "from the signed-propagation substrate."
    )
    return {
        "track": track,
        "n_nodes": int(len(nodes)),
        "n_edges": int(len(edges)),
        "n_components": int(len(components)),
        "largest_component_size": int(sizes[0]) if sizes else 0,
        "component_size_histogram_top10": sizes[:10],
        "n_isolated_nodes": int(len(isolates)),
        "note": n_excluded_note,
    }


def build_substrates(
    causal_edges: pd.DataFrame,
    string_edges: pd.DataFrame | None,
    *,
    plan: dict[str, object] | None = None,
) -> dict[str, object]:
    plan = plan or load_plan()
    n_causal = int(len(causal_edges))
    n_excluded = int(
        (
            causal_edges["relation_type"].isin(["ambiguous", "unsigned_directed"])
            | ~causal_edges["sign"].isin([1, -1])
        ).sum()
    )
    signed = signed_causal_edges(causal_edges, plan)
    full_unsigned = unsigned_full_causal_edges(causal_edges, plan)
    primary_nodes = substrate_nodes(signed)
    # Full causal node set (pre-exclusion) for isolate accounting relative to the TARGET/BENCHMARK graph
    all_causal_nodes = substrate_nodes(causal_edges)
    conn_signed = connectivity_report(signed, all_nodes=all_causal_nodes, track="signed_primary")
    conn_signed["n_causal_edges_before_exclusion"] = n_causal
    conn_signed["n_edges_excluded"] = n_excluded
    conn_signed["fraction_edges_excluded"] = float(n_excluded / max(n_causal, 1))
    conn_full = connectivity_report(
        full_unsigned, all_nodes=all_causal_nodes, track="unsigned_full_directed"
    )
    conn_full["n_causal_edges_before_exclusion"] = n_causal
    conn_full["n_ambiguous_excluded"] = int((causal_edges["relation_type"] == "ambiguous").sum())
    primary = {
        "name": plan["phase5_proximity_substrate_primary"],
        "signed_edges": signed,
        "full_unsigned_edges": full_unsigned,
        "causal_edges": causal_edges,
        "nodes": primary_nodes,
        "full_unsigned_nodes": substrate_nodes(full_unsigned),
        "n_signed_edges": int(len(signed)),
        "n_full_unsigned_edges": int(len(full_unsigned)),
        "n_nodes": int(len(primary_nodes)),
        "connectivity": conn_signed,
        "connectivity_unsigned_full": conn_full,
    }
    if string_edges is None or string_edges.empty:
        combined_nodes = primary_nodes
        sens_edges = signed
    else:
        string_nodes = substrate_nodes(string_edges)
        combined_nodes = primary_nodes | string_nodes
        sens_edges = signed
    conn_sens = connectivity_report(sens_edges, all_nodes=combined_nodes, track="sensitivity_signed_core")
    sensitivity = {
        "name": plan["phase5_proximity_substrate_sensitivity"],
        "signed_edges": signed,
        "string_edges": string_edges,
        "nodes": combined_nodes,
        "n_signed_edges": int(len(signed)),
        "n_string_edges": int(len(string_edges)) if string_edges is not None else 0,
        "n_nodes": int(len(combined_nodes)),
        "string_policy": plan["phase5_string_unsigned_edge_policy"],
        "connectivity": conn_sens,
    }
    return {"primary": primary, "sensitivity": sensitivity}


def load_phase3_substrates(bench_dir: str | Path) -> dict[str, object]:
    bench_dir = Path(bench_dir)
    causal = pd.read_csv(bench_dir / "graph_edges.tsv", sep="\t")
    string_path = bench_dir / "graph_edges_string_aux.tsv"
    string_edges = pd.read_csv(string_path, sep="\t") if string_path.exists() else pd.DataFrame()
    return build_substrates(causal, string_edges)
