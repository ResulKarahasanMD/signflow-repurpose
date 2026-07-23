"""Baselines: random, degree, PageRank, LR, HGB, MLP, label propagation, RWR."""

from __future__ import annotations

import numpy as np
import pandas as pd
import networkx as nx
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.semi_supervised import LabelPropagation

from signflow.phase3.features import FEATURE_COLUMNS


def _graph_from_edges(edges: pd.DataFrame, nodes: list[str]) -> nx.DiGraph:
    g = nx.DiGraph()
    g.add_nodes_from(nodes)
    for _, row in edges.iterrows():
        if row["source"] in g and row["target"] in g:
            g.add_edge(row["source"], row["target"], sign=int(row["sign"]))
    return g


def score_random(n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.random(n)


def score_degree(nodes: list[str], edges: pd.DataFrame) -> np.ndarray:
    deg = pd.concat([edges["source"], edges["target"]]).value_counts()
    return np.array([float(deg.get(g, 0.0)) for g in nodes], dtype=float)


def score_pagerank(nodes: list[str], edges: pd.DataFrame) -> np.ndarray:
    g = _graph_from_edges(edges, nodes).to_undirected()
    pr = nx.pagerank(g, alpha=0.85) if g.number_of_edges() else {n: 0.0 for n in nodes}
    return np.array([float(pr.get(n, 0.0)) for n in nodes], dtype=float)


def rwr_scores(
    nodes: list[str],
    edges: pd.DataFrame,
    seed_genes: list[str],
    *,
    alpha: float = 0.85,
    tol: float = 1e-8,
    max_iter: int = 100,
) -> np.ndarray:
    """Random-walk-with-restart from seed genes on the undirected projection."""

    from scipy import sparse

    index = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    if n == 0 or not seed_genes:
        return np.zeros(n)
    src = edges["source"].astype(str).to_numpy()
    tgt = edges["target"].astype(str).to_numpy()
    keep = np.array([s in index and t in index for s, t in zip(src, tgt)], dtype=bool)
    if not keep.any():
        return np.zeros(n)
    ii = np.array([index[s] for s in src[keep]], dtype=int)
    jj = np.array([index[t] for t in tgt[keep]], dtype=int)
    rows = np.concatenate([ii, jj])
    cols = np.concatenate([jj, ii])
    data = np.ones(len(rows), dtype=float)
    A = sparse.csr_matrix((data, (rows, cols)), shape=(n, n))
    A.data[:] = 1.0
    A = (A > 0).astype(float)
    row_sum = np.asarray(A.sum(axis=1)).ravel()
    row_sum[row_sum == 0] = 1.0
    P = sparse.diags(1.0 / row_sum) @ A
    e = np.zeros(n)
    seeds = [index[s] for s in seed_genes if s in index]
    if not seeds:
        return np.zeros(n)
    e[seeds] = 1.0 / len(seeds)
    r = e.copy()
    for _ in range(max_iter):
        r_new = alpha * (P.T @ r) + (1 - alpha) * e
        if np.linalg.norm(r_new - r, 1) < tol:
            r = r_new
            break
        r = r_new
    return np.asarray(r).ravel()


def fit_predict_tabular(
    model_name: str,
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_score: np.ndarray,
    *,
    seed: int,
) -> np.ndarray:
    """Fit scaler on train only; return decision scores for X_score."""

    scaler = StandardScaler()
    Xtr = scaler.fit_transform(X_train)
    Xs = scaler.transform(X_score)
    if model_name == "logistic":
        clf = LogisticRegression(max_iter=2000, random_state=seed, class_weight="balanced")
    elif model_name == "hgb":
        clf = HistGradientBoostingClassifier(random_state=seed, max_depth=4)
    elif model_name == "mlp":
        clf = MLPClassifier(
            hidden_layer_sizes=(64, 32),
            max_iter=200,
            random_state=seed,
            early_stopping=True,
        )
    else:
        raise ValueError(model_name)
    clf.fit(Xtr, y_train)
    if hasattr(clf, "predict_proba"):
        return clf.predict_proba(Xs)[:, 1]
    return clf.decision_function(Xs)


def label_propagation_scores(
    nodes: list[str],
    edges: pd.DataFrame,
    y: pd.Series,
    train_genes: list[str],
    *,
    seed: int,
) -> np.ndarray:
    """Clamp train labels and diffuse on sparse undirected adjacency."""

    from scipy import sparse

    index = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    labels = np.full(n, -1, dtype=int)
    for gene in train_genes:
        if gene in index and y.loc[gene] >= 0:
            labels[index[gene]] = int(y.loc[gene])
    if (labels >= 0).sum() < 2 or len(np.unique(labels[labels >= 0])) < 2:
        return np.zeros(n)
    src = edges["source"].astype(str).to_numpy()
    tgt = edges["target"].astype(str).to_numpy()
    keep = np.fromiter((s in index and t in index for s, t in zip(src, tgt)), dtype=bool, count=len(src))
    ii = np.fromiter((index[s] for s in src[keep]), dtype=int, count=int(keep.sum()))
    jj = np.fromiter((index[t] for t in tgt[keep]), dtype=int, count=int(keep.sum()))
    rows = np.concatenate([ii, jj])
    cols = np.concatenate([jj, ii])
    A = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    A = (A > 0).astype(float)
    row_sum = np.asarray(A.sum(axis=1)).ravel()
    row_sum[row_sum == 0] = 1.0
    P = sparse.diags(1.0 / row_sum) @ A
    known = labels >= 0
    probs = np.zeros((n, 2))
    probs[known, labels[known]] = 1.0
    for _ in range(30):
        probs = P @ probs
        probs[known] = 0.0
        probs[known, labels[known]] = 1.0
    return np.asarray(probs[:, 1]).ravel()


def feature_matrix(features: pd.DataFrame, columns: list[str] | None = None) -> np.ndarray:
    cols = columns or FEATURE_COLUMNS
    return features[cols].to_numpy(dtype=float)
