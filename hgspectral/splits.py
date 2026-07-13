"""
Split hyperedges (experiments) into train/val/test, stratified by hyperedge
size so small and large complexes are represented fairly in every fold.

Train defines the vertex universe (clustering is done on training proteins).
Val is used for nested lambda/k selection. Test is used only for final scoring.
"""
from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np


def size_bin(sz: int) -> str:
    if sz <= 3:
        return "small"
    if sz <= 6:
        return "medium"
    return "large"


def split_hyperedges(hyperedges: List[List[str]],
                     fracs=(0.6, 0.2, 0.2),
                     seed: int = 0,
                     stratify_by_size: bool = True
                     ) -> Tuple[List, List, List, Dict]:
    rng = np.random.RandomState(seed)
    m = len(hyperedges)
    idx = np.arange(m)
    if stratify_by_size:
        bins: Dict[str, list] = {}
        for i, e in enumerate(hyperedges):
            bins.setdefault(size_bin(len(e)), []).append(i)
        tr, va, te = [], [], []
        for b, ids in bins.items():
            ids = np.array(ids); rng.shuffle(ids)
            n = len(ids)
            n_tr = int(round(fracs[0] * n))
            n_va = int(round(fracs[1] * n))
            tr += ids[:n_tr].tolist()
            va += ids[n_tr:n_tr + n_va].tolist()
            te += ids[n_tr + n_va:].tolist()
    else:
        rng.shuffle(idx)
        n_tr = int(round(fracs[0] * m)); n_va = int(round(fracs[1] * m))
        tr, va, te = idx[:n_tr].tolist(), idx[n_tr:n_tr + n_va].tolist(), idx[n_tr + n_va:].tolist()

    g = lambda L: [hyperedges[i] for i in L]
    stats = dict(m_total=m, m_train=len(tr), m_val=len(va), m_test=len(te),
                 stratified=stratify_by_size)
    return g(tr), g(va), g(te), stats


def to_index_members(edges: List[List[str]], pid: Dict[str, int]) -> List[np.ndarray]:
    """Map protein-id hyperedges to index arrays, dropping unseen proteins."""
    out = []
    for e in edges:
        idx = np.array([pid[p] for p in e if p in pid], dtype=int)
        if len(idx) >= 2:
            out.append(idx)
    return out
