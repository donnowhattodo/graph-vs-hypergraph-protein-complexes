"""
go_analysis.py -- GO over-representation analysis and REVIGO-style semantic
redundancy reduction, with publication-grade figures.

Design
------
* Enrichment (ORA): hypergeometric test + Benjamini-Hochberg FDR
  (implemented in hgspectral.enrichment; re-exported here for convenience).
* REVIGO-style reduction: reduce a long, redundant list of enriched GO terms to
  a small set of representatives.
    - PRIMARY (DAG-based, publication standard): if a `go-basic.obo` is supplied
      and `goatools` is installed, information content is estimated from term
      annotation frequency and Lin/Resnik semantic similarity is used. This is
      the true REVIGO recipe.
    - FALLBACK (no external files): term-term similarity is estimated from the
      overlap (Jaccard) of the gene sets annotated to each term. Clearly labelled
      as an annotation-overlap approximation.
  Terms are grouped by average-linkage on the similarity, a representative is
  chosen per group (most significant), and 2-D coordinates are produced by
  classical MDS for the canonical REVIGO scatter.
* Figures: REVIGO semantic scatter, GO dot-plot of top terms, and a graph-vs-
  hypergraph enrichment comparison. All saved as PDF + PNG at 300 dpi.

Inputs are the pipeline's per-term table `<dataset>_go_terms.csv` (columns:
model, tag/dataset, cluster_id, go_id, K, k, N, p_value, q_value) and, for the
gene-overlap similarity, a protein->GO map. Everything degrades gracefully when
optional pieces are missing.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt  # type: ignore[import]
from scipy.cluster.hierarchy import linkage, fcluster  # type: ignore[import]
from scipy.spatial.distance import squareform  # type: ignore[import]

from hgspectral.enrichment import benjamini_hochberg, enrich_partition  # noqa: F401

# ----------------------------------------------------------------------
# Semantic similarity
# ----------------------------------------------------------------------
def jaccard_similarity(term_genes: Dict[str, set]) -> (List[str], np.ndarray):
    """Annotation-overlap (Jaccard) similarity between GO terms. No external files."""
    terms = list(term_genes.keys())
    n = len(terms)
    S = np.eye(n)
    for i in range(n):
        gi = term_genes[terms[i]]
        for j in range(i + 1, n):
            gj = term_genes[terms[j]]
            u = len(gi | gj)
            s = len(gi & gj) / u if u else 0.0
            S[i, j] = S[j, i] = s
    return terms, S


def lin_similarity_dag(go_ids: List[str], obo_path: str,
                       prot2go: Dict[str, set]) -> Optional[np.ndarray]:
    """
    DAG-based Lin semantic similarity (primary REVIGO backend).
    Requires goatools + a go-basic.obo file. Returns None if unavailable so the
    caller falls back to the annotation-overlap similarity.
    """
    try:
        from goatools.obo_parser import GODag
        from goatools.semantic import TermCounts, lin_sim
    except Exception:
        return None
    if not obo_path or not Path(obo_path).exists():
        return None
    godag = GODag(obo_path)
    assoc = {}
    for p, gos in prot2go.items():
        assoc.setdefault(p, set()).update(gos)
    tc = TermCounts(godag, assoc)
    n = len(go_ids)
    S = np.eye(n)
    for i in range(n):
        for j in range(i + 1, n):
            try:
                s = lin_sim(go_ids[i], go_ids[j], godag, tc)
            except Exception:
                s = 0.0
            s = 0.0 if s is None or np.isnan(s) else float(s)
            S[i, j] = S[j, i] = s
    return S


# ----------------------------------------------------------------------
# REVIGO-style reduction
# ----------------------------------------------------------------------
def revigo_reduce(enriched: pd.DataFrame,
                  term_genes: Optional[Dict[str, set]] = None,
                  prot2go: Optional[Dict[str, set]] = None,
                  obo_path: Optional[str] = None,
                  similarity_threshold: float = 0.5,
                  q_col: str = "q_value",
                  go_col: str = "go_id") -> pd.DataFrame:
    """
    Reduce redundant enriched GO terms to representatives.

    enriched : one row per (significant) GO term with columns [go_id, q_value,
               k (gene count)]; duplicates across clusters are aggregated to the
               most significant occurrence.
    Returns a dataframe of unique terms with: representative flag, group id,
    dispensability (1 - sim to representative), MDS coords (x, y), and size.
    """
    df = (enriched.sort_values(q_col)
                  .drop_duplicates(subset=[go_col])
                  .reset_index(drop=True)).copy()
    go_ids = df[go_col].tolist()
    n = len(go_ids)
    if n == 0:
        return df.assign(representative=[], group=[], x=[], y=[], dispensability=[])
    if n == 1:
        return df.assign(representative=[True], group=[0], x=[0.0], y=[0.0],
                         dispensability=[0.0], neg_log10_q=[-np.log10(max(df[q_col].iloc[0], 1e-300))])

    # similarity: DAG (primary) or annotation overlap (fallback)
    S = None
    backend = "annotation-overlap (Jaccard)"
    if obo_path and prot2go:
        S = lin_similarity_dag(go_ids, obo_path, prot2go)
        if S is not None:
            backend = "GO-DAG Lin semantic similarity"
    if S is None:
        if term_genes is None:
            raise ValueError("Provide term_genes (protein sets per GO term) or an OBO + prot2go.")
        tg = {g: term_genes.get(g, set()) for g in go_ids}
        _, S = jaccard_similarity(tg)

    # cluster terms on distance = 1 - similarity
    D = 1.0 - S
    np.fill_diagonal(D, 0.0)
    D = np.clip((D + D.T) / 2, 0, 1)
    Z = linkage(squareform(D, checks=False), method="average")
    groups = fcluster(Z, t=1.0 - similarity_threshold, criterion="distance")

    # representative per group = most significant term
    df["group"] = groups
    df["neg_log10_q"] = -np.log10(np.clip(df[q_col].values, 1e-300, 1.0))
    rep_idx = df.groupby("group")["neg_log10_q"].idxmax().values
    df["representative"] = df.index.isin(rep_idx)

    # dispensability = similarity to the group's representative (0 for reps)
    disp = np.zeros(n)
    for g in np.unique(groups):
        members = np.where(groups == g)[0]
        r = df.index.get_loc(df[df.group == g]["neg_log10_q"].idxmax())
        for m in members:
            disp[m] = 1.0 - S[m, r]
    df["dispensability"] = disp

    # 2-D layout via classical MDS on the similarity (REVIGO-style scatter)
    df["x"], df["y"] = _classical_mds(S)
    df.attrs["similarity_backend"] = backend
    return df


def _classical_mds(S: np.ndarray) -> (np.ndarray, np.ndarray):
    D2 = (1.0 - S) ** 2
    n = D2.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ D2 @ J
    w, V = np.linalg.eigh(B)
    idx = np.argsort(w)[::-1][:2]
    L = np.sqrt(np.clip(w[idx], 0, None))
    X = V[:, idx] * L
    return X[:, 0], X[:, 1]


# ----------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------
def plot_revigo_scatter(reduced: pd.DataFrame, out: Path, title: str,
                        name_col: Optional[str] = None, top_labels: int = 12):
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    size = 40 + 260 * (reduced["neg_log10_q"] / max(reduced["neg_log10_q"].max(), 1e-9))
    sc = ax.scatter(reduced["x"], reduced["y"], s=size,
                    c=reduced["neg_log10_q"], cmap="viridis",
                    alpha=0.85, edgecolor="k", linewidth=0.4)
    reps = reduced[reduced["representative"]].sort_values("neg_log10_q", ascending=False).head(top_labels)
    for _, r in reps.iterrows():
        lab = str(r[name_col]) if name_col and name_col in reduced.columns else str(r["go_id"])
        ax.annotate(lab[:28], (r["x"], r["y"]), fontsize=7,
                    xytext=(3, 3), textcoords="offset points")
    cb = fig.colorbar(sc, ax=ax); cb.set_label(r"$-\log_{10} q$")
    ax.set_xlabel("semantic space x"); ax.set_ylabel("semantic space y")
    ax.set_title(title)
    backend = reduced.attrs.get("similarity_backend", "")
    if backend:
        ax.text(0.99, 0.01, backend, transform=ax.transAxes, ha="right",
                va="bottom", fontsize=6, color="0.4")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{out}.{ext}", dpi=300)
    plt.close(fig)


def plot_go_dotplot(enriched: pd.DataFrame, out: Path, title: str,
                    top_n: int = 15, name_col: Optional[str] = None,
                    model_col: str = "model"):
    """Top-N enriched terms per model as a dot plot (x=-log10 q, size=gene count)."""
    df = enriched.copy()
    df["neg_log10_q"] = -np.log10(np.clip(df["q_value"].values, 1e-300, 1.0))
    label = name_col if (name_col and name_col in df.columns) else "go_id"
    top = (df.sort_values("neg_log10_q", ascending=False)
             .drop_duplicates(subset=[label]).head(top_n)
             .sort_values("neg_log10_q"))
    fig, ax = plt.subplots(figsize=(6.6, 0.38 * len(top) + 1.2))
    models = list(df[model_col].unique()) if model_col in df.columns else [None]
    cmap = {m: c for m, c in zip(models, plt.cm.tab10.colors)}
    for _, r in top.iterrows():
        color = cmap.get(r.get(model_col), "steelblue")
        ax.scatter(r["neg_log10_q"], str(r[label])[:40],
                   s=30 + 12 * r.get("k", 5), color=color, edgecolor="k", linewidth=0.4)
    ax.set_xlabel(r"$-\log_{10} q$ (BH-FDR)")
    ax.set_title(title); ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{out}.{ext}", dpi=300)
    plt.close(fig)


def plot_enrichment_comparison(summary: pd.DataFrame, out: Path):
    """Graph vs hypergraph enriched-fraction (clusters & proteins) across datasets."""
    piv_c = summary.pivot_table(index="dataset", columns="method",
                                values="enriched_fraction_clusters")
    piv_p = summary.pivot_table(index="dataset", columns="method",
                                values="enriched_fraction_proteins")
    order = piv_c.index.tolist()
    x = np.arange(len(order)); w = 0.38
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, piv, ttl in ((axes[0], piv_c, "enriched clusters"),
                         (axes[1], piv_p, "enriched proteins")):
        if "graph_spectral" in piv:
            ax.bar(x - w / 2, piv["graph_spectral"].loc[order], w, label="Graph", color="#4C72B0")
        if "hyper_penalized" in piv:
            ax.bar(x + w / 2, piv["hyper_penalized"].loc[order], w, label="Hypergraph", color="#DD8452")
        ax.set_xticks(x); ax.set_xticklabels(order, rotation=35, ha="right")
        ax.set_title(f"Fraction of {ttl} ($q\\leq 0.05$)"); ax.set_ylim(0, 1.05)
        ax.grid(axis="y", alpha=0.3)
    axes[0].set_ylabel("fraction"); axes[0].legend(frameon=False)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{out}.{ext}", dpi=300)
    plt.close(fig)
