"""
Consume ALL_test_metrics.csv and produce: (1) per-(dataset,method) summary with
bootstrap CIs, (2) paired graph-vs-hypergraph Wilcoxon + effect sizes,
(3) Friedman + Nemenyi across methods, (4) publication figures (PDF+PNG).

Usage: python analyze_results.py --results fin_result/ALL_test_metrics.csv --out fin_result/analysis
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from hgspectral import stats

PRIMARY = "bestmatch_f1_sym"
METRICS = ["pairwise_f1", "bestmatch_f1_sym", "b3_f1", "NMI", "ARI", "hyperedge_recovery_f1"]


def summarize(df, out):
    rows = []
    for (ds, method, kmode), g in df.groupby(["dataset", "method", "k_mode"]):
        row = dict(dataset=ds, method=method, k_mode=kmode, n_seeds=len(g))
        for m in METRICS:
            if m in g:
                v = g[m].dropna().values
                row[f"{m}_mean"] = float(np.mean(v)) if len(v) else np.nan
                row[f"{m}_std"] = float(np.std(v)) if len(v) else np.nan
        rows.append(row)
    s = pd.DataFrame(rows)
    s.to_csv(out / "summary_by_method.csv", index=False)
    return s


def paired_tests(df, out):
    """Graph vs hypergraph (matched-k) per dataset, and pooled."""
    rows = []
    md = df[df.k_mode == "matched"]
    for ds, g in md.groupby("dataset"):
        gg = g[g.method == "graph_spectral"].sort_values("seed")
        hh = g[g.method == "hyper_penalized"].sort_values("seed")
        common = sorted(set(gg.seed) & set(hh.seed))
        if len(common) < 3:
            continue
        x = hh.set_index("seed").loc[common, PRIMARY].values
        y = gg.set_index("seed").loc[common, PRIMARY].values
        w = stats.paired_wilcoxon(x, y)
        ci = stats.bootstrap_ci(x, y)
        rows.append(dict(dataset=ds, metric=PRIMARY, hyper_mean=float(np.mean(x)),
                         graph_mean=float(np.mean(y)), **w, **{f"ci_{k}": v for k, v in ci.items()}))
    res = pd.DataFrame(rows)
    res.to_csv(out / "paired_graph_vs_hyper.csv", index=False)
    return res


def friedman_across_methods(df, out):
    md = df[df.k_mode == "matched"]
    wide = md.groupby(["dataset", "method"])[PRIMARY].mean().unstack("method")
    keep = [c for c in ["graph_spectral", "hyper_zhou_lam0", "hyper_penalized", "louvain", "leiden"]
            if c in wide.columns]
    wide = wide[keep]
    res = stats.friedman_nemenyi(wide)
    with open(out / "friedman_nemenyi.txt", "w") as f:
        f.write(f"Friedman p = {res.get('friedman_p')}\n\nMean ranks (1=best):\n{res.get('mean_ranks')}\n")
        if res.get("nemenyi") is not None:
            f.write(f"\nNemenyi post-hoc p-values:\n{res['nemenyi']}\n")
    return res


def fig_graph_vs_hyper(summary, out):
    md = summary[summary.k_mode == "matched"]
    piv = md.pivot_table(index="dataset", columns="method", values=f"{PRIMARY}_mean")
    err = md.pivot_table(index="dataset", columns="method", values=f"{PRIMARY}_std")
    for col in ["graph_spectral", "hyper_penalized"]:
        if col not in piv:
            return
    x = np.arange(len(piv.index)); w = 0.38
    fig, ax = plt.subplots(figsize=(max(6, len(piv) * 1.1), 4.2))
    ax.bar(x - w / 2, piv["graph_spectral"], w, yerr=err.get("graph_spectral"),
           capsize=3, label="Graph spectral")
    ax.bar(x + w / 2, piv["hyper_penalized"], w, yerr=err.get("hyper_penalized"),
           capsize=3, label="Hypergraph (penalized)")
    ax.set_xticks(x); ax.set_xticklabels(piv.index, rotation=30, ha="right")
    ax.set_ylabel("Held-out best-match $F_1$ (matched $k$)")
    ax.set_title("Graph vs hypergraph on held-out test hyperedges (mean ± SD over seeds)")
    ax.legend(frameon=False); ax.grid(axis="y", alpha=0.3); fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(out / f"graph_vs_hyper_bestmatchF1.{ext}", dpi=300)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    df = pd.read_csv(args.results)
    out = Path(args.out or (Path(args.results).parent / "analysis"))
    out.mkdir(parents=True, exist_ok=True)
    s = summarize(df, out)
    paired_tests(df, out)
    friedman_across_methods(df, out)
    fig_graph_vs_hyper(s, out)
    print(f"[ANALYSIS] wrote summaries, tests, and figures to {out}")


if __name__ == "__main__":
    main()
