"""
make_paper_figures.py -- generate publication-grade figures from the completed
pipeline results in fin_result/, plus a REVIGO demonstration.

Figures (real data):
  F1 forest_plot_hyper_vs_graph      -- paired improvement + 95% CI + significance
  F2 matched_k_bestmatch_by_method   -- graph / Zhou(lam=0) / penalized, mean+-SD
  F3 critical_difference             -- Friedman/Nemenyi CD diagram
  F4 go_enrichment_comparison        -- graph vs hyper enriched fractions
  F5 lambda_selected                 -- validation-selected lambda* per dataset
  F6 real_vs_null_controls           -- structure vs degree/size-preserving nulls

REVIGO demo (synthetic SynthA, clearly labelled): runs ORA + revigo_reduce and
saves a semantic scatter + GO dot-plot to prove the go_analysis code end-to-end.
"""
from pathlib import Path
import glob, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = Path("fin_result")
OUT = RES / "figures"
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

SPECTRAL = ["graph_spectral", "hyper_zhou_lam0", "hyper_penalized"]
PRETTY = {"graph_spectral": "Graph spectral", "hyper_zhou_lam0": r"Hypergraph ($\lambda{=}0$)",
          "hyper_penalized": "Hypergraph (penalized)"}
COLORS = {"graph_spectral": "#4C72B0", "hyper_zhou_lam0": "#55A868", "hyper_penalized": "#DD8452"}


def stars(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "n.s."


# ---------- F1: forest plot of paired improvement ----------
def forest_plot():
    p = pd.read_csv(RES / "analysis" / "paired_graph_vs_hyper.csv")
    p = p.sort_values("ci_mean_diff")
    y = np.arange(len(p))
    fig, ax = plt.subplots(figsize=(7, 0.5 * len(p) + 1.4))
    ax.axvline(0, color="0.6", lw=1, ls="--")
    ax.hlines(y, p["ci_lo"], p["ci_hi"], color="#333", lw=2)
    ax.scatter(p["ci_mean_diff"], y, s=60, color="#DD8452", zorder=3, edgecolor="k")
    for yi, (_, r) in zip(y, p.iterrows()):
        ax.text(r["ci_hi"] + 0.001, yi,
                f"{stars(r['p_value'])}  (r={r['rank_biserial']:.2f})",
                va="center", fontsize=8)
    ax.set_yticks(y); ax.set_yticklabels(p["dataset"])
    ax.set_xlabel(r"$\Delta$ best-match $F_1$ (hypergraph $-$ graph), 95% bootstrap CI")
    ax.set_title("Paired improvement on held-out test hyperedges (matched $k$, 25 seeds)")
    fig.tight_layout()
    for e in ("pdf", "png"):
        fig.savefig(OUT / f"F1_forest_plot_hyper_vs_graph.{e}", dpi=300)
    plt.close(fig)


# ---------- F2: matched-k bar by method ----------
def matched_k_bars():
    s = pd.read_csv(RES / "analysis" / "summary_by_method.csv")
    s = s[s.k_mode == "matched"]
    datasets = sorted(s.dataset.unique())
    x = np.arange(len(datasets)); w = 0.26
    fig, ax = plt.subplots(figsize=(10, 4.4))
    for i, mth in enumerate(SPECTRAL):
        sub = s[s.method == mth].set_index("dataset")
        means = [sub.loc[d, "bestmatch_f1_sym_mean"] if d in sub.index else np.nan for d in datasets]
        errs = [sub.loc[d, "bestmatch_f1_sym_std"] if d in sub.index else 0 for d in datasets]
        ax.bar(x + (i - 1) * w, means, w, yerr=errs, capsize=2.5,
               label=PRETTY[mth], color=COLORS[mth])
    ax.set_xticks(x); ax.set_xticklabels(datasets, rotation=30, ha="right")
    ax.set_ylabel("held-out best-match $F_1$"); ax.legend(frameon=False, ncol=3)
    ax.set_title("Matched-$k$ comparison (mean $\\pm$ SD over 25 seeds)")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    for e in ("pdf", "png"):
        fig.savefig(OUT / f"F2_matched_k_bestmatch_by_method.{e}", dpi=300)
    plt.close(fig)


# ---------- F3: critical-difference (Nemenyi) diagram ----------
def critical_difference():
    s = pd.read_csv(RES / "analysis" / "summary_by_method.csv")
    s = s[(s.k_mode == "matched") & (s.method.isin(SPECTRAL))]
    wide = s.pivot_table(index="dataset", columns="method", values="bestmatch_f1_sym_mean")[SPECTRAL].dropna()
    ranks = wide.rank(axis=1, ascending=False).mean(axis=0)
    N, k = wide.shape[0], wide.shape[1]
    q_alpha = 2.343  # Nemenyi, k=3, alpha=0.05
    CD = q_alpha * np.sqrt(k * (k + 1) / (6.0 * N))

    order = ranks.sort_values()
    lo, hi = 1, k
    fig, ax = plt.subplots(figsize=(7, 2.6))
    ax.set_xlim(lo - 0.2, hi + 0.2); ax.set_ylim(0, 1); ax.axis("off")
    ax.hlines(0.8, lo, hi, color="k")
    for r in range(lo, hi + 1):
        ax.vlines(r, 0.78, 0.82, color="k"); ax.text(r, 0.86, str(r), ha="center")
    for i, (mth, rk) in enumerate(order.items()):
        yv = 0.55 - 0.14 * i
        ax.plot([rk, rk], [0.8, yv], color=COLORS[mth])
        ax.plot([rk, lo - 0.15 if i % 2 == 0 else hi + 0.15], [yv, yv], color=COLORS[mth])
        ax.text(lo - 0.18 if i % 2 == 0 else hi + 0.18, yv,
                f"{PRETTY[mth]} ({rk:.2f})", va="center",
                ha="right" if i % 2 == 0 else "left", fontsize=9)
    # CD bar
    ax.hlines(0.95, lo, lo + CD, color="k", lw=3)
    ax.text(lo + CD / 2, 0.98, f"CD = {CD:.2f}", ha="center", fontsize=8)
    ax.set_title(f"Critical-difference diagram (Friedman+Nemenyi, N={N} datasets)", fontsize=10)
    fig.tight_layout()
    for e in ("pdf", "png"):
        fig.savefig(OUT / f"F3_critical_difference.{e}", dpi=300)
    plt.close(fig)


# ---------- F4: GO enrichment comparison ----------
def go_comparison():
    import go_analysis
    frames = [pd.read_csv(f) for f in glob.glob(str(RES / "*" / "*_go_enrichment_summary.csv"))]
    if not frames:
        return
    summ = pd.concat(frames, ignore_index=True)
    go_analysis.plot_enrichment_comparison(summ, OUT / "F4_go_enrichment_comparison")


# ---------- F5: selected lambda* ----------
def lambda_plot():
    frames = [pd.read_csv(f) for f in glob.glob(str(RES / "*" / "*_manifest.csv"))]
    man = pd.concat(frames, ignore_index=True).sort_values("lambda_star")
    fig, ax = plt.subplots(figsize=(7.5, 4))
    lam = man["lambda_star"].replace(0, 1e-5)
    ax.barh(man["dataset"], lam, color="#8172B3")
    ax.set_xscale("log"); ax.set_xlabel(r"validation-selected penalty $\lambda^\star$ (log scale)")
    ax.set_title(r"Selected regularization strength per dataset")
    fig.tight_layout()
    for e in ("pdf", "png"):
        fig.savefig(OUT / f"F5_lambda_selected.{e}", dpi=300)
    plt.close(fig)


# ---------- F6: real vs null controls ----------
def controls_plot():
    rows = []
    for f in glob.glob(str(RES / "*" / "*_controls.csv")):
        d = pd.read_csv(f); ds = Path(f).stem.replace("_controls", "")
        for ctl, g in d.groupby("control"):
            rows.append(dict(dataset=ds, control=ctl, bm=g["bestmatch_f1_sym"].mean()))
    ctl = pd.DataFrame(rows)
    real = (pd.read_csv(RES / "analysis" / "summary_by_method.csv")
              .query("k_mode=='matched' and method=='hyper_penalized'")
              .set_index("dataset")["bestmatch_f1_sym_mean"])
    datasets = sorted(real.index)
    x = np.arange(len(datasets)); w = 0.26
    fig, ax = plt.subplots(figsize=(10, 4.4))
    ax.bar(x - w, [real.get(d, np.nan) for d in datasets], w, label="Real (hypergraph)", color="#DD8452")
    for j, c in enumerate(["degree_preserving_rewire", "hyperedge_size_preserving"]):
        vals = [ctl.query("dataset==@d and control==@c")["bm"].mean() for d in datasets]
        ax.bar(x + j * w, vals, w, label=c.replace("_", " "), color=["#4C72B0", "#55A868"][j])
    ax.set_xticks(x); ax.set_xticklabels(datasets, rotation=30, ha="right")
    ax.set_ylabel("best-match $F_1$"); ax.legend(frameon=False)
    ax.set_title("Real structure vs randomized null models (held-out test)")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    for e in ("pdf", "png"):
        fig.savefig(OUT / f"F6_real_vs_null_controls.{e}", dpi=300)
    plt.close(fig)


# ---------- REVIGO demo on synthetic SynthA ----------
def revigo_demo():
    import go_analysis
    from hgspectral.mitab import load_mitab, mitab_to_complexes, build_protein_go
    from hgspectral.operators import build_incidence, hypergraph_laplacian, index_proteins
    from hgspectral.splits import to_index_members
    from hgspectral.spectral import spectral_cluster
    from hgspectral.enrichment import enrich_partition
    from collections import defaultdict
    synth = Path("SynthA/SynthA.txt")
    if not synth.exists():
        synth = Path("smoke_out/SynthA/SynthA.txt")
    if not synth.exists():
        print("[revigo demo] SynthA not found; skipping"); return
    df, _ = load_mitab(synth, top_k_taxids=3)
    complexes = mitab_to_complexes(df, min_size=2)
    proteins = sorted({p for e in complexes.values() for p in e})
    pid = index_proteins(proteins)
    prot2go, _ = build_protein_go(df, proteins)
    members = to_index_members(list(complexes.values()), pid)
    H, sizes = build_incidence(len(proteins), members)
    L, _, _ = hypergraph_laplacian(H, sizes, penalty=0.1)
    lab = spectral_cluster(L, 8, seed=0)
    enr = enrich_partition(proteins, lab, prot2go, "hyper", "SynthA", min_term_size=3, min_overlap=2)
    if enr.empty:
        print("[revigo demo] no enriched terms; skipping"); return
    sig = enr[enr["q_value"] <= 0.1]
    term_genes = defaultdict(set)
    for p, gos in prot2go.items():
        for g in gos:
            term_genes[g].add(p)
    reduced = go_analysis.revigo_reduce(sig, term_genes=term_genes, similarity_threshold=0.4)
    go_analysis.plot_revigo_scatter(reduced, OUT / "DEMO_revigo_scatter_SYNTHETIC",
                                    "REVIGO semantic scatter (SYNTHETIC demo)")
    go_analysis.plot_go_dotplot(sig, OUT / "DEMO_go_dotplot_SYNTHETIC",
                                "Top enriched GO terms (SYNTHETIC demo)")
    print(f"[revigo demo] {len(sig)} sig terms -> {reduced['representative'].sum()} representatives "
          f"({reduced.attrs.get('similarity_backend')})")


if __name__ == "__main__":
    forest_plot();        print("F1 forest plot")
    matched_k_bars();     print("F2 matched-k bars")
    critical_difference();print("F3 critical-difference")
    go_comparison();      print("F4 GO enrichment comparison")
    lambda_plot();        print("F5 lambda selected")
    controls_plot();      print("F6 real vs null controls")
    revigo_demo()
    print("\nFigures written to", OUT.resolve())
