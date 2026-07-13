"""
Statistical testing across datasets/seeds for method comparison.

- paired_wilcoxon: graph vs hypergraph on matched (dataset, seed) pairs.
- friedman_nemenyi: >2 methods across datasets, with post-hoc Nemenyi.
- bootstrap_ci: CI for a paired performance difference.
- Effect sizes: rank-biserial correlation (Wilcoxon) and Cliff's delta.
Report effect sizes alongside p-values, never p-values alone.
"""
from __future__ import annotations

from typing import Dict, List

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, friedmanchisquare


def rank_biserial(x: np.ndarray, y: np.ndarray) -> float:
    d = np.asarray(x) - np.asarray(y)
    d = d[d != 0]
    if len(d) == 0:
        return 0.0
    ranks = np.argsort(np.argsort(np.abs(d))) + 1
    Rpos = ranks[d > 0].sum()
    Rneg = ranks[d < 0].sum()
    T = Rpos + Rneg
    return float((Rpos - Rneg) / T) if T else 0.0


def cliffs_delta(x, y) -> float:
    x, y = np.asarray(x), np.asarray(y)
    gt = sum((xi > y).sum() for xi in x)
    lt = sum((xi < y).sum() for xi in x)
    n = len(x) * len(y)
    return float((gt - lt) / n) if n else 0.0


def paired_wilcoxon(x, y) -> Dict:
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 3 or np.allclose(x, y):
        return dict(n=len(x), statistic=np.nan, p_value=np.nan,
                    rank_biserial=0.0, median_diff=float(np.median(x - y)))
    stat, p = wilcoxon(x, y)
    return dict(n=len(x), statistic=float(stat), p_value=float(p),
                rank_biserial=rank_biserial(x, y),
                median_diff=float(np.median(x - y)))


def bootstrap_ci(x, y, n_boot: int = 10000, alpha: float = 0.05, seed: int = 0) -> Dict:
    rng = np.random.RandomState(seed)
    d = np.asarray(x, float) - np.asarray(y, float)
    n = len(d)
    if n == 0:
        return dict(mean_diff=np.nan, lo=np.nan, hi=np.nan)
    boots = np.array([d[rng.randint(0, n, n)].mean() for _ in range(n_boot)])
    lo, hi = np.percentile(boots, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return dict(mean_diff=float(d.mean()), lo=float(lo), hi=float(hi))


def friedman_nemenyi(wide: pd.DataFrame) -> Dict:
    """
    wide: rows = blocks (datasets), columns = methods, values = a metric.
    Returns Friedman result, mean ranks, and Nemenyi post-hoc matrix.
    """
    data = wide.dropna(axis=0, how="any")
    methods = list(data.columns)
    if data.shape[0] < 3 or data.shape[1] < 3:
        return dict(friedman_p=np.nan, mean_ranks=None, nemenyi=None,
                    note="need >=3 blocks and >=3 methods")
    arrays = [data[m].values for m in methods]
    stat, p = friedmanchisquare(*arrays)
    ranks = data.rank(axis=1, ascending=False)  # rank 1 = best
    mean_ranks = ranks.mean(axis=0)
    nem = None
    try:
        import scikit_posthocs as sp
        nem = sp.posthoc_nemenyi_friedman(data.values)
        nem.index = methods; nem.columns = methods
    except Exception:
        pass
    return dict(friedman_stat=float(stat), friedman_p=float(p),
                mean_ranks=mean_ranks, nemenyi=nem)
