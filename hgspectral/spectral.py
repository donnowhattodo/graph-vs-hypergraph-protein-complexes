"""
Spectral embedding + k-means, with explicit matched-k and adaptive-k modes.

adaptive-k criteria: 'eigengap', 'silhouette', 'modularity'. Each method may
pick its own k in adaptive mode; in matched mode both operators receive the
same externally supplied k.
"""
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


def smallest_eigs(L: sparse.spmatrix, k: int):
    """
    k smallest eigenpairs of a symmetric PSD Laplacian L.

    ARPACK's which='SM' (smallest magnitude) frequently fails to converge on
    large, spectrally-clustered normalized Laplacians (the Cancer failure). The
    robust recipe is shift-invert with sigma strictly *below* zero: since L is
    PSD, (L - sigma I) = L + |sigma| I is positive definite and well
    conditioned, and which='LM' on the shift-inverted operator returns the
    eigenvalues of L nearest sigma, i.e. the smallest. We try progressively
    stronger shifts, then fall back to LOBPCG and finally a dense solve.
    """
    n = L.shape[0]
    k = int(min(max(2, k), n - 1))
    L = L.astype(float).tocsc()

    # small problems: dense is fastest and always converges
    if n <= 400:
        vals, vecs = np.linalg.eigh(L.toarray())
        return vals[:k], vecs[:, :k]

    maxiter = max(2000, 40 * n)
    # 1) shift-invert from below zero (primary, most robust)
    for sigma in (-1e-3, -1e-2, -1e-1, -1.0):
        try:
            vals, vecs = eigsh(L, k=k, sigma=sigma, which="LM", maxiter=maxiter, tol=1e-6)
            order = np.argsort(vals)
            return vals[order], vecs[:, order]
        except Exception:
            continue
    # 2) smallest-algebraic with a generous iteration budget
    try:
        vals, vecs = eigsh(L, k=k, which="SA", maxiter=maxiter, tol=1e-5)
        order = np.argsort(vals)
        return vals[order], vecs[:, order]
    except Exception:
        pass
    # 3) LOBPCG with a random orthonormal block (handles the hardest spectra)
    try:
        from scipy.sparse.linalg import lobpcg
        rng = np.random.RandomState(0)
        X = rng.rand(n, k)
        X, _ = np.linalg.qr(X)
        M = sparse.identity(n, format="csc")
        vals, vecs = lobpcg(L, X, M=None, tol=1e-5, maxiter=2000, largest=False)
        order = np.argsort(vals)
        return vals[order], vecs[:, order]
    except Exception:
        pass
    # 4) last resort: dense (only reached if everything above failed)
    vals, vecs = np.linalg.eigh(L.toarray())
    return vals[:k], vecs[:, :k]


def embed(L: sparse.spmatrix, k_embed: int) -> Tuple[np.ndarray, np.ndarray]:
    vals, vecs = smallest_eigs(L, k_embed + 1)
    U = vecs[:, 1:k_embed + 1]  # skip trivial eigenvector
    norms = np.linalg.norm(U, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    return U / norms, vals


def kmeans_labels(Y: np.ndarray, k: int, seed: int) -> np.ndarray:
    k = int(min(max(2, k), Y.shape[0]))
    km = KMeans(n_clusters=k, n_init=10, random_state=seed)
    return km.fit_predict(Y)


def eigengap_k(vals: np.ndarray, kmax: int = 15) -> int:
    ev = np.sort(np.asarray(vals))
    kcap = int(min(kmax, len(ev) - 1))
    if kcap < 2:
        return 2
    gaps = np.diff(ev[:kcap + 1])
    return int(max(2, np.argmax(gaps[1:]) + 2))


def spectral_cluster(L, k: int, seed: int = 0, return_embed: bool = False):
    """Matched-k spectral clustering: embed with k, k-means with k."""
    Y, vals = embed(L, k)
    labels = kmeans_labels(Y, k, seed)
    if return_embed:
        return labels, Y, vals
    return labels


def spectral_cluster_adaptive(L, A_for_mod=None, criterion: str = "eigengap",
                              k_grid=None, kmax: int = 15, seed: int = 0):
    """
    Adaptive-k spectral clustering. Returns (labels, k, embedding, diagnostics).
    - 'eigengap': one-shot k from the spectrum.
    - 'silhouette'/'modularity': sweep k_grid, score, pick best.
    """
    from .metrics import modularity_weighted
    n = L.shape[0]
    if k_grid is None:
        k_grid = list(range(2, min(kmax, n - 1) + 1))
    if criterion == "eigengap":
        Ymax, vals = embed(L, max(k_grid))
        k = eigengap_k(vals, kmax=max(k_grid))
        Y = Ymax[:, :k]
        labels = kmeans_labels(Y, k, seed)
        return labels, k, Y, {"criterion": "eigengap"}

    best = (-np.inf, None, None, None)
    Ymax, _ = embed(L, max(k_grid))
    for k in k_grid:
        Y = Ymax[:, :k]
        labels = kmeans_labels(Y, k, seed)
        if len(set(labels)) < 2:
            continue
        if criterion == "silhouette":
            score = silhouette_score(Y, labels)
        elif criterion == "modularity":
            score = modularity_weighted(A_for_mod, labels)
        else:
            raise ValueError(criterion)
        if score > best[0]:
            best = (score, k, labels, Y)
    _, k, labels, Y = best
    return labels, k, Y, {"criterion": criterion, "score": best[0]}
