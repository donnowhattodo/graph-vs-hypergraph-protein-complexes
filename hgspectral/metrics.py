"""
Evaluation metrics (vectorized / sparse) against overlapping complex gold.

All external metrics take a gold set of complexes (as member-index arrays) that
can be the HELD-OUT test hyperedges, enabling leakage-free evaluation.
"""
from __future__ import annotations

from typing import Dict, List

import numpy as np
from scipy import sparse
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


def _labels_to_indicator(labels: np.ndarray) -> sparse.csr_matrix:
    labels = np.asarray(labels)
    uniq = {l: i for i, l in enumerate(sorted(set(labels.tolist())))}
    rows = np.arange(len(labels))
    cols = np.array([uniq[l] for l in labels])
    data = np.ones(len(labels))
    return sparse.csr_matrix((data, (rows, cols)), shape=(len(labels), len(uniq)))


def gold_incidence(n: int, gold_members: List[np.ndarray]) -> sparse.csr_matrix:
    rows, cols = [], []
    for j, idx in enumerate(gold_members):
        idx = np.asarray(idx, int)
        rows.extend(idx.tolist())
        cols.extend([j] * len(idx))
    return sparse.csr_matrix((np.ones(len(rows)), (rows, cols)),
                             shape=(n, max(1, len(gold_members))))


def gold_adjacency(n: int, gold_members: List[np.ndarray]) -> sparse.csr_matrix:
    Hg = gold_incidence(n, gold_members)
    A = (Hg @ Hg.T)
    A = (A > 0).astype(float).tolil()
    A.setdiag(0.0)
    return A.tocsr()


def pairwise_prf1(labels: np.ndarray, gold_adj: sparse.csr_matrix):
    C = _labels_to_indicator(labels)
    sizes = np.asarray(C.sum(axis=0)).ravel()
    pred_pairs = float(np.sum(sizes * (sizes - 1) / 2.0))
    gold_pairs = float(gold_adj.nnz) / 2.0
    tp = 0.0
    labels = np.asarray(labels)
    for l in set(labels.tolist()):
        idx = np.where(labels == l)[0]
        if len(idx) < 2:
            continue
        sub = gold_adj[np.ix_(idx, idx)]
        tp += sub.nnz / 2.0
    P = tp / pred_pairs if pred_pairs else 0.0
    R = tp / gold_pairs if gold_pairs else 0.0
    F1 = 2 * P * R / (P + R) if (P + R) else 0.0
    return P, R, F1


def bestmatch_f1(labels: np.ndarray, gold_members: List[np.ndarray], n: int):
    """Symmetric best-match F1 (average of pred->gold and gold->pred)."""
    C = _labels_to_indicator(labels)
    Hg = gold_incidence(n, gold_members)
    O = (C.T @ Hg).toarray()                       # k x m overlaps
    csize = np.asarray(C.sum(axis=0)).ravel()[:, None]
    gsize = np.asarray(Hg.sum(axis=0)).ravel()[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        P = np.where(csize > 0, O / csize, 0.0)
        R = np.where(gsize > 0, O / gsize, 0.0)
        F = np.where((P + R) > 0, 2 * P * R / (P + R), 0.0)
    # pred->gold
    csize_f = csize.ravel()
    pg = (F.max(axis=1) * csize_f).sum() / max(csize_f.sum(), 1)
    gsize_f = gsize.ravel()
    gp = (F.max(axis=0) * gsize_f).sum() / max(gsize_f.sum(), 1)
    return float(pg), float(gp), float(0.5 * (pg + gp))


def b3_overlapping(labels: np.ndarray, gold_members: List[np.ndarray], n: int):
    labels = np.asarray(labels)
    prot_gold = [set() for _ in range(n)]
    for idx in gold_members:
        s = set(int(x) for x in idx)
        for i in s:
            prot_gold[i] |= s
    clusters = {l: set(np.where(labels == l)[0].tolist()) for l in set(labels.tolist())}
    precs, recs = [], []
    for i in range(n):
        G = set(prot_gold[i]); G.discard(i)
        C = set(clusters[labels[i]]); C.discard(i)
        inter = len(C & G)
        precs.append(inter / max(len(C), 1))
        recs.append(inter / max(len(G), 1))
    P, R = float(np.mean(precs)), float(np.mean(recs))
    F1 = 2 * P * R / (P + R) if (P + R) else 0.0
    return P, R, F1


def hard_gold_labels(n: int, gold_members: List[np.ndarray]) -> np.ndarray:
    gold = -np.ones(n, dtype=int)
    for cid, idx in enumerate(gold_members):
        for i in idx:
            if gold[i] == -1:
                gold[i] = cid
    return gold


def nmi_ari(labels: np.ndarray, gold_labels: np.ndarray):
    mask = gold_labels >= 0
    if mask.sum() < 2 or len(set(gold_labels[mask].tolist())) < 2 or len(set(labels[mask].tolist())) < 2:
        return float("nan"), float("nan")
    return (float(normalized_mutual_info_score(gold_labels[mask], labels[mask])),
            float(adjusted_rand_score(gold_labels[mask], labels[mask])))


def hyperedge_recovery(labels: np.ndarray, test_members: List[np.ndarray], n: int,
                       f1_threshold: float = 0.25):
    """
    Held-out hyperedge recovery: for each TEST complex, best F1 against any
    predicted cluster restricted to proteins seen in training (index space).
    Returns mean best-F1 and fraction of test complexes recovered above threshold.
    """
    labels = np.asarray(labels)
    clusters = [set(np.where(labels == l)[0].tolist()) for l in set(labels.tolist())]
    bests = []
    for idx in test_members:
        G = set(int(x) for x in idx if 0 <= int(x) < n)
        if len(G) < 2:
            continue
        best = 0.0
        for S in clusters:
            inter = len(S & G)
            if not inter:
                continue
            P = inter / len(S); R = inter / len(G)
            f1 = 2 * P * R / (P + R) if (P + R) else 0.0
            best = max(best, f1)
        bests.append(best)
    if not bests:
        return float("nan"), float("nan")
    bests = np.asarray(bests)
    return float(bests.mean()), float((bests >= f1_threshold).mean())


def modularity_weighted(A: sparse.csr_matrix, labels: np.ndarray) -> float:
    if A is None:
        return float("nan")
    labels = np.asarray(labels)
    deg = np.asarray(A.sum(axis=1)).ravel()
    m2 = float(A.sum())
    if m2 <= 0:
        return 0.0
    Q = 0.0
    for l in set(labels.tolist()):
        idx = np.where(labels == l)[0]
        e_in = float(A[np.ix_(idx, idx)].sum())
        d = float(deg[idx].sum())
        Q += e_in - d * d / m2
    return Q / m2


def evaluate(labels, gold_members, n, A=None, embed=None, gold_labels=None):
    """Bundle of metrics for one clustering vs a (test) gold set."""
    gadj = gold_adjacency(n, gold_members)
    pP, pR, pF = pairwise_prf1(labels, gadj)
    bmg, bmp, bsym = bestmatch_f1(labels, gold_members, n)
    b3P, b3R, b3F = b3_overlapping(labels, gold_members, n)
    if gold_labels is None:
        gold_labels = hard_gold_labels(n, gold_members)
    nmi, ari = nmi_ari(np.asarray(labels), gold_labels)
    rec_f1, rec_frac = hyperedge_recovery(labels, gold_members, n)
    out = dict(
        k=int(len(set(np.asarray(labels).tolist()))),
        pairwise_precision=pP, pairwise_recall=pR, pairwise_f1=pF,
        bestmatch_pred_gold=bmg, bestmatch_gold_pred=bmp, bestmatch_f1_sym=bsym,
        b3_precision=b3P, b3_recall=b3R, b3_f1=b3F,
        NMI=nmi, ARI=ari,
        hyperedge_recovery_f1=rec_f1, hyperedge_recovery_frac=rec_frac,
        modularity=modularity_weighted(A, labels) if A is not None else float("nan"),
    )
    return out
