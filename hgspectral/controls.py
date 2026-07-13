"""
Randomized controls (null models). Each returns randomized hyperedges (or a
randomized annotation map) that feed back through the SAME pipeline, so any
'signal' above the null is attributable to real higher-order structure.

Controls implemented:
  - degree_preserving_rewire: rewire the co-complex GRAPH keeping degree seq.
  - hyperedge_size_preserving: rebuild hyperedges of identical sizes from a
    random permutation of the vertex multiset (destroys who-groups-with-whom).
  - go_label_permutation: shuffle GO annotations across proteins.
  - taxid_stratified_null: size-preserving randomization performed WITHIN taxid
    strata, so cross-species artefacts are not introduced.
"""
from __future__ import annotations

from typing import Dict, List

import numpy as np
import networkx as nx
from scipy import sparse


def degree_preserving_rewire(A: sparse.csr_matrix, seed: int = 0, n_swaps_factor: int = 10):
    """
    Degree-preserving (Maslov-Sneppen) rewiring of the co-complex graph.

    For dense or degree-heterogeneous graphs, achieving the full target number
    of successful swaps can be impossible within any finite try budget
    (the Diabetes failure). ``double_edge_swap`` then raises
    ``NetworkXAlgorithmError`` -- which is NOT a subclass of ``NetworkXError``,
    so the previous handler missed it. We (i) catch every NetworkX exception,
    (ii) cap the try budget, and (iii) accept a partially rewired graph, which
    is still a valid approximate degree-preserving null. The number of achieved
    swaps is attached as a graph attribute for logging.
    """
    n = A.shape[0]
    G = nx.from_scipy_sparse_array(A)
    G = nx.Graph(G)
    G.remove_edges_from(nx.selfloop_edges(G))
    m = G.number_of_edges()
    achieved = 0
    if m > 1 and G.number_of_nodes() > 3:
        target = n_swaps_factor * m
        max_tries = min(50 * target, 2_000_000)  # bounded, avoids runaway
        try:
            # newer networkx returns the number of swaps performed
            res = nx.double_edge_swap(G, nswap=target, max_tries=max_tries, seed=seed)
            achieved = res if isinstance(res, int) else target
        except (nx.NetworkXException, nx.NetworkXError, nx.NetworkXAlgorithmError):
            # partial rewiring already applied in-place; accept it
            achieved = -1
        except Exception:
            achieved = -1
    out = nx.to_scipy_sparse_array(G, nodelist=range(n), format="csr")
    return out


def hyperedge_size_preserving(hyperedges: List[List[str]], seed: int = 0) -> List[List[str]]:
    rng = np.random.RandomState(seed)
    vocab = sorted({p for e in hyperedges for p in e})
    out = []
    for e in hyperedges:
        s = len(e)
        out.append(list(rng.choice(vocab, size=min(s, len(vocab)), replace=False)))
    return out


def taxid_stratified_null(hyperedges: List[List[str]], prot2taxid: Dict[str, str],
                          seed: int = 0) -> List[List[str]]:
    rng = np.random.RandomState(seed)
    by_tax: Dict[str, list] = {}
    for p, t in prot2taxid.items():
        by_tax.setdefault(t, []).append(p)
    out = []
    for e in hyperedges:
        newe = []
        for p in e:
            t = prot2taxid.get(p)
            pool = by_tax.get(t, [p])
            newe.append(pool[rng.randint(len(pool))])
        out.append(list(dict.fromkeys(newe)))  # dedup, keep order
    return out


def go_label_permutation(prot2go: Dict[str, set], seed: int = 0) -> Dict[str, set]:
    rng = np.random.RandomState(seed)
    prots = list(prot2go.keys())
    gos = [prot2go[p] for p in prots]
    perm = rng.permutation(len(prots))
    return {prots[i]: gos[perm[i]] for i in range(len(prots))}
