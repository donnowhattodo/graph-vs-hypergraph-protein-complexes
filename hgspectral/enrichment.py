"""
GO over-representation analysis (hypergeometric) with Benjamini-Hochberg FDR.

Uses scipy's hypergeometric survival function. Ontology (BP/MF/CC) separation
and semantic-similarity redundancy reduction require a real GO DAG; hooks are
left for goatools but the core enrichment works from the MITAB-parsed GO map.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Dict, List

import numpy as np
import pandas as pd
from scipy.stats import hypergeom


def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, float)
    m = len(p)
    if m == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * m / (np.arange(1, m + 1))
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(m)
    q[order] = np.clip(ranked, 0, 1)
    return q


def enrich_partition(proteins: List[str], labels: np.ndarray,
                     prot2go: Dict[str, set], model: str, tag: str,
                     min_cluster_size: int = 5, min_term_size: int = 5,
                     min_overlap: int = 2) -> pd.DataFrame:
    labels = np.asarray(labels)
    idx_go = [prot2go.get(proteins[i], set()) for i in range(len(proteins))]
    univ = [i for i, g in enumerate(idx_go) if g]
    N = len(univ)
    if N == 0:
        return pd.DataFrame()
    univ_set = set(univ)
    term_counts = defaultdict(int)
    for i in univ:
        for go in idx_go[i]:
            term_counts[go] += 1
    rows = []
    for cid in sorted(set(labels.tolist())):
        members = [i for i in np.where(labels == cid)[0] if i in univ_set]
        nk = len(members)
        if nk < min_cluster_size:
            continue
        cc = defaultdict(int)
        for i in members:
            for go in idx_go[i]:
                cc[go] += 1
        for go, k in cc.items():
            K = term_counts[go]
            if K < min_term_size or k < min_overlap:
                continue
            pval = hypergeom.sf(k - 1, N, K, nk)
            rows.append(dict(model=model, tag=tag, cluster_id=int(cid),
                             cluster_size=nk, go_id=go, K=K, k=k, N=N,
                             p_value=float(pval)))
    df = pd.DataFrame(rows)
    if not df.empty:
        df["q_value"] = benjamini_hochberg(df["p_value"].values)
    return df


def enrichment_summary(proteins, labels, enr_df: pd.DataFrame, q_thr=0.05) -> dict:
    labels = np.asarray(labels)
    total = len(set(labels.tolist()))
    if enr_df is None or enr_df.empty:
        return dict(total_clusters=total, enriched_clusters=0,
                    enriched_fraction_clusters=0.0, enriched_fraction_proteins=0.0)
    sig = enr_df[enr_df["q_value"] <= q_thr]
    eids = set(sig["cluster_id"].unique().tolist())
    frac_c = len(eids) / total if total else 0.0
    frac_p = float(np.isin(labels, list(eids)).mean()) if eids else 0.0
    return dict(total_clusters=total, enriched_clusters=len(eids),
                enriched_fraction_clusters=frac_c, enriched_fraction_proteins=frac_p)
