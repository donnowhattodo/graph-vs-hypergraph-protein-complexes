"""
Graph and hypergraph operators (sparse).

Representation-controlled comparison: both operators are built from the SAME
hyperedge set. The graph is the clique/co-complex expansion; the hypergraph is
the Zhou normalized Laplacian, with an optional degree-aware penalty lambda*Dv^-1.
"""
from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np
from scipy import sparse


def index_proteins(proteins: List[str]) -> Dict[str, int]:
    return {p: i for i, p in enumerate(proteins)}


def build_incidence(n: int, comp_members: List[np.ndarray]) -> Tuple[sparse.csr_matrix, np.ndarray]:
    """Sparse incidence H (n x m) and hyperedge sizes."""
    rows, cols = [], []
    sizes = np.zeros(len(comp_members), dtype=float)
    for j, idx in enumerate(comp_members):
        idx = np.asarray(idx, dtype=int)
        sizes[j] = len(idx)
        rows.extend(idx.tolist())
        cols.extend([j] * len(idx))
    data = np.ones(len(rows), dtype=float)
    H = sparse.csr_matrix((data, (rows, cols)), shape=(n, max(1, len(comp_members))))
    return H, sizes


def hypergraph_laplacian(H: sparse.csr_matrix, sizes: np.ndarray,
                         weights: np.ndarray = None, penalty: float = 0.0):
    """
    Zhou normalized hypergraph Laplacian
        L0 = I - Dv^{-1/2} H W De^{-1} H^T Dv^{-1/2}
    plus optional degree-aware penalty  L = L0 + penalty * Dv^{-1}.
    Returns (L, Dv, Dv_inv).
    """
    n, m = H.shape
    W = np.ones(m) if weights is None else np.asarray(weights, float)
    De = np.where(sizes > 0, sizes, 1.0)
    Dv = np.asarray((H @ sparse.diags(W)).sum(axis=1)).ravel()
    Dv_safe = np.where(Dv > 0, Dv, 1.0)
    Dv_isqrt = sparse.diags(1.0 / np.sqrt(Dv_safe))
    HW = H @ sparse.diags(W)
    inner = HW @ sparse.diags(1.0 / De) @ H.T           # n x n
    G = Dv_isqrt @ inner @ Dv_isqrt
    L = sparse.identity(n, format="csr") - G
    Dv_inv = 1.0 / Dv_safe
    if penalty and penalty > 0:
        L = L + penalty * sparse.diags(Dv_inv)
    return L.tocsr(), Dv, Dv_inv


def graph_from_incidence(H: sparse.csr_matrix, sizes: np.ndarray,
                         weight_by_inverse: bool = True):
    """
    Co-complex (clique-expanded) adjacency A and normalized graph Laplacian.
    Each hyperedge of size s contributes weight 1/(s-1) (if weight_by_inverse)
    or 1 to every internal pair. Returns (A, L, deg).
    """
    n = H.shape[0]
    if weight_by_inverse:
        wcol = np.where(sizes > 1, 1.0 / (sizes - 1), 0.0)
    else:
        wcol = np.ones_like(sizes)
    Hw = H @ sparse.diags(wcol)
    A = (Hw @ H.T).tolil()
    A.setdiag(0.0)
    A = A.tocsr()
    A.eliminate_zeros()
    deg = np.asarray(A.sum(axis=1)).ravel()
    deg_safe = np.where(deg > 0, deg, 1.0)
    Disq = sparse.diags(1.0 / np.sqrt(deg_safe))
    L = sparse.identity(n, format="csr") - Disq @ A @ Disq
    return A.tocsr(), L.tocsr(), deg
