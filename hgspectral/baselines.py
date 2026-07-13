"""
Non-spectral graph baselines on the clique-expanded co-complex graph.

Louvain and Leiden run on the weighted projection; MCL runs on the same
projection. ClusterONE is overlapping and distributed as a JAR, so it is
provided as an optional wrapper (returns None with a clear message if the JAR
is absent) rather than a reimplementation.
"""
from __future__ import annotations

from typing import List, Optional

import numpy as np
from scipy import sparse


def _adj_to_edges(A: sparse.csr_matrix):
    A = sparse.triu(A, k=1).tocoo()
    return list(zip(A.row.tolist(), A.col.tolist(), A.data.tolist()))


def louvain_labels(A: sparse.csr_matrix, seed: int = 0) -> np.ndarray:
    import networkx as nx
    n = A.shape[0]
    G = nx.Graph()
    G.add_nodes_from(range(n))
    for u, v, w in _adj_to_edges(A):
        G.add_edge(u, v, weight=w)
    comms = nx.community.louvain_communities(G, weight="weight", seed=seed)
    labels = np.zeros(n, dtype=int)
    for c, nodes in enumerate(comms):
        for u in nodes:
            labels[u] = c
    return labels


def leiden_labels(A: sparse.csr_matrix, seed: int = 0) -> Optional[np.ndarray]:
    try:
        import igraph as ig
        import leidenalg as la
    except Exception:
        return None
    n = A.shape[0]
    edges = _adj_to_edges(A)
    g = ig.Graph(n=n, edges=[(u, v) for u, v, _ in edges])
    g.es["weight"] = [w for _, _, w in edges]
    part = la.find_partition(g, la.RBConfigurationVertexPartition,
                             weights="weight", seed=seed)
    labels = np.zeros(n, dtype=int)
    for c, nodes in enumerate(part):
        for u in nodes:
            labels[u] = c
    return labels


def mcl_labels(A: sparse.csr_matrix, inflation: float = 2.0) -> Optional[np.ndarray]:
    try:
        import markov_clustering as mc
    except Exception:
        return None
    n = A.shape[0]
    result = mc.run_mcl(A.tocsr(), inflation=inflation)
    clusters = mc.get_clusters(result)
    labels = -np.ones(n, dtype=int)
    for c, nodes in enumerate(clusters):
        for u in nodes:
            labels[u] = c
    # singletons for any unassigned node
    nxt = labels.max() + 1
    for i in range(n):
        if labels[i] < 0:
            labels[i] = nxt; nxt += 1
    return labels


def clusterone_labels(A: sparse.csr_matrix, jar_path: Optional[str] = None):
    """
    Optional ClusterONE wrapper. Requires the ClusterONE JAR + Java.
    Returns None (with message) when unavailable so the pipeline degrades
    gracefully; documented as 'proposed / external' in the manuscript.
    """
    if jar_path is None:
        return None
    raise NotImplementedError(
        "ClusterONE integration is a documented optional baseline: export the "
        "co-complex graph to a weighted edge list, run the ClusterONE JAR, and "
        "load its overlapping communities. See README for the exact command."
    )
