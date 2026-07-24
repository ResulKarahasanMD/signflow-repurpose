"""Signed signaling-flow propagation and module opposition score."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse


def node_degrees(edges: pd.DataFrame, nodes: list[str]) -> np.ndarray:
    """Undirected degree on signed edge projection."""

    index = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    if edges.empty:
        return np.zeros(n, dtype=float)
    src = edges["source"].astype(str).to_numpy()
    tgt = edges["target"].astype(str).to_numpy()
    keep = np.array([s in index and t in index for s, t in zip(src, tgt)], dtype=bool)
    if not keep.any():
        return np.zeros(n, dtype=float)
    ii = np.array([index[s] for s in src[keep]], dtype=int)
    jj = np.array([index[t] for t in tgt[keep]], dtype=int)
    rows = np.concatenate([ii, jj])
    deg = np.bincount(rows, minlength=n).astype(float)
    return deg


def build_signed_adjacency(
    edges: pd.DataFrame,
    nodes: list[str],
    *,
    ignore_sign: bool = False,
) -> sparse.csr_matrix:
    """Build row-normalized adjacency. If ignore_sign, treat every retained edge as +1."""

    index = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    for _, row in edges.iterrows():
        s, t = str(row["source"]), str(row["target"])
        if s not in index or t not in index:
            continue
        sign = int(row["sign"])
        if ignore_sign:
            # Allow unsigned_directed (sign=0) as unit-weight connectivity.
            if sign not in (-1, 0, 1):
                continue
            weight = 1.0
        else:
            if sign not in (1, -1):
                continue
            weight = float(sign)
        rows.append(index[s])
        cols.append(index[t])
        data.append(weight)
    if not rows:
        return sparse.csr_matrix((n, n), dtype=float)
    mat = sparse.csr_matrix((data, (rows, cols)), shape=(n, n), dtype=float)
    # Row-normalize signed out-edges so flow iteration remains stable (|row| sum -> 1).
    abs_mat = mat.copy()
    abs_mat.data = np.abs(abs_mat.data)
    row_sum = np.asarray(abs_mat.sum(axis=1)).ravel()
    row_sum[row_sum == 0] = 1.0
    inv = sparse.diags(1.0 / row_sum)
    return inv @ mat


def signed_flow(
    nodes: list[str],
    signed_edges: pd.DataFrame,
    injection: dict[str, float],
    *,
    alpha: float = 0.85,
    max_iter: int = 100,
    tol: float = 1e-8,
    ignore_sign: bool = False,
) -> np.ndarray:
    """f = (1-alpha)*s + alpha * A^T f."""

    n = len(nodes)
    index = {g: i for i, g in enumerate(nodes)}
    s = np.zeros(n, dtype=float)
    for gene, val in injection.items():
        if gene in index:
            s[index[gene]] = float(val)
    # Non-zero mixed signed injections may sum to 0; do not early-return on sum alone.
    if not np.any(s):
        return s
    A = build_signed_adjacency(signed_edges, nodes, ignore_sign=ignore_sign)
    f = s.copy()
    for _ in range(max_iter):
        f_new = (1.0 - alpha) * s + alpha * (A.T @ f)
        if np.linalg.norm(f_new - f, 1) < tol:
            f = f_new
            break
        f = f_new
    return np.asarray(f).ravel()


def drug_injection_vector(
    targets: pd.DataFrame,
    *,
    unknown_policy: str = "exclude_from_signed_score_retain_in_coverage",
) -> tuple[dict[str, float], dict[str, object]]:
    """Build normalized signed injection; unknown actions excluded from denominator."""

    signed_rows = targets.loc[targets["signed_action"].notna() & (targets["action_sign"].notna())].copy()
    n_signed = int(len(signed_rows))
    diag = {
        "n_targets_total": int(len(targets)),
        "n_targets_signed": n_signed,
        "unknown_policy": unknown_policy,
    }
    if n_signed == 0:
        return {}, diag
    denom = float(n_signed)
    injection: dict[str, float] = {}
    for _, row in signed_rows.iterrows():
        gene = str(row["gene_symbol"])
        sign = int(row["action_sign"])
        injection[gene] = injection.get(gene, 0.0) + sign / denom
    return injection, diag


def module_opposition_score(
    flow: np.ndarray,
    nodes: list[str],
    module_genes: set[str],
    disease_direction: pd.Series,
) -> float:
    """S(d) = -sum_{v in M} f_v * y_v."""

    index = {g: i for i, g in enumerate(nodes)}
    total = 0.0
    for gene in module_genes:
        if gene not in index or gene not in disease_direction.index:
            continue
        y = float(disease_direction.loc[gene])
        if not np.isfinite(y):
            continue
        total += flow[index[gene]] * y
    return float(-total)

def signed_flow_adjacency(
    adjacency_T: sparse.csr_matrix,
    injection: dict[str, float],
    node_index: dict[str, int],
    n: int,
    *,
    alpha: float = 0.85,
    max_iter: int = 100,
    tol: float = 1e-8,
) -> np.ndarray:
    """Flow using pre-built A^T (sparse)."""
    s = np.zeros(n, dtype=float)
    for gene, val in injection.items():
        if gene in node_index:
            s[node_index[gene]] = float(val)
    if not np.any(s):
        return s
    f = s.copy()
    for _ in range(max_iter):
        f_new = (1.0 - alpha) * s + alpha * (adjacency_T @ f)
        if np.linalg.norm(f_new - f, 1) < tol:
            return np.asarray(f_new).ravel()
        f = f_new
    return np.asarray(f).ravel()


class CachedFlowScorer:
    """Reuse signed adjacency and module projection for repeated null draws."""

    def __init__(
        self,
        nodes: list[str],
        signed_edges: pd.DataFrame,
        module_genes: set[str],
        disease_direction: pd.Series,
        *,
        alpha: float,
        max_iter: int,
        tol: float,
        default_action_sign: int = -1,
        ignore_sign: bool = False,
        variant: str = "signed",
    ) -> None:
        self.nodes = nodes
        self.node_index = {g: i for i, g in enumerate(nodes)}
        self.n = len(nodes)
        self.ignore_sign = ignore_sign
        self.variant = variant
        self.adjacency_T = build_signed_adjacency(
            signed_edges, nodes, ignore_sign=ignore_sign
        ).T.tocsr()
        self.alpha = alpha
        self.max_iter = max_iter
        self.tol = tol
        self.default_action_sign = default_action_sign
        # align module_idx with module_y lengths by filtering together
        idx_y = []
        ys = []
        for g in module_genes:
            if g not in self.node_index or g not in disease_direction.index:
                continue
            y = float(disease_direction.loc[g])
            if not np.isfinite(y):
                continue
            idx_y.append(self.node_index[g])
            ys.append(y)
        self.module_idx = np.array(idx_y, dtype=int)
        self.module_y = np.array(ys, dtype=float)

    def score_genes(self, genes: list[str], *, action_sign: int | None = None) -> float:
        sign = self.default_action_sign if action_sign is None else action_sign
        genes = [g for g in genes if g in self.node_index]
        if not genes:
            return 0.0
        denom = float(len(genes))
        injection = {g: sign / denom for g in genes}
        flow = signed_flow_adjacency(
            self.adjacency_T, injection, self.node_index, self.n,
            alpha=self.alpha, max_iter=self.max_iter, tol=self.tol,
        )
        if self.module_idx.size == 0:
            return 0.0
        return float(-np.dot(flow[self.module_idx], self.module_y))

    def score_drug_targets(self, drug_targets: pd.DataFrame, *, unknown_policy: str) -> float:
        injection, _ = drug_injection_vector(drug_targets, unknown_policy=unknown_policy)
        if not injection:
            return 0.0
        flow = signed_flow_adjacency(
            self.adjacency_T, injection, self.node_index, self.n,
            alpha=self.alpha, max_iter=self.max_iter, tol=self.tol,
        )
        if self.module_idx.size == 0:
            return 0.0
        return float(-np.dot(flow[self.module_idx], self.module_y))

