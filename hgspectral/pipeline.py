"""
End-to-end orchestration for one dataset and the full benchmark.

Design guarantees the audit called for:
  * Vertex universe = TRAIN hyperedges only (no test-protein leakage).
  * lambda AND k are selected on VALIDATION hyperedges, never on test.
  * Final metrics are computed on HELD-OUT TEST hyperedges.
  * Matched-k: every method gets the SAME k. Adaptive-k: each picks its own.
  * Multiple seeds for k-means so results carry dispersion, not a single point.
  * Randomized controls run through the identical evaluation path.
"""
from __future__ import annotations

import time
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
from scipy import sparse

from . import baselines, controls, metrics, spectral
from .config import Config
from .enrichment import enrich_partition, enrichment_summary
from .mitab import build_protein_go, find_mitab_file, load_mitab, mitab_to_complexes
from .operators import build_incidence, graph_from_incidence, hypergraph_laplacian, index_proteins
from .splits import split_hyperedges, to_index_members


def _matched_k(cfg: Config, n_train: int) -> int:
    if cfg.k_matched == "rule":
        return int(max(2, min(n_train // cfg.min_cluster_size, cfg.kmax, n_train - 1)))
    return int(cfg.k_matched)


def select_lambda_k_on_val(H, sizes, A, val_members, n, cfg: Config, seed: int):
    """Nested selection: pick (lambda, k) maximising held-out VAL bestmatch-F1."""
    best = (-np.inf, cfg.lambda_grid[0], _matched_k(cfg, n))
    k_grid = list(range(2, min(cfg.kmax, n - 1) + 1))
    for lam in cfg.lambda_grid:
        L, _, _ = hypergraph_laplacian(H, sizes, penalty=lam)
        Ymax, _ = spectral.embed(L, max(k_grid))
        for k in k_grid:
            labels = spectral.kmeans_labels(Ymax[:, :k], k, seed)
            _, _, bsym = metrics.bestmatch_f1(labels, val_members, n)
            if bsym > best[0]:
                best = (bsym, lam, k)
    return best[1], best[2], best[0]


def run_dataset(name: str, cfg: Config) -> Dict[str, pd.DataFrame]:
    path = find_mitab_file(name, cfg.search_roots)
    if path is None:
        print(f"[SKIP] {name}: MITAB file not found under {cfg.search_roots}")
        return {}
    out_dir = Path(cfg.out_root) / name
    out_dir.mkdir(parents=True, exist_ok=True)

    df, load_stats = load_mitab(path, top_k_taxids=cfg.top_k_taxids)
    complexes = mitab_to_complexes(df, min_size=cfg.min_complex_size)
    if len(complexes) < 6:
        print(f"[SKIP] {name}: too few complexes ({len(complexes)})")
        return {}
    hyperedges = list(complexes.values())

    tr, va, te, split_stats = split_hyperedges(
        hyperedges, fracs=cfg.split_fracs, seed=cfg.split_seed,
        stratify_by_size=cfg.stratify_by_size)

    proteins = sorted({p for e in tr for p in e})
    n = len(proteins)
    if n < 6:
        print(f"[SKIP] {name}: too few train proteins ({n})")
        return {}
    pid = index_proteins(proteins)
    prot2go, _ = build_protein_go(df, proteins)

    train_members = to_index_members(tr, pid)
    val_members = to_index_members(va, pid)
    test_members = to_index_members(te, pid)
    if len(test_members) < 3 or len(val_members) < 3:
        print(f"[SKIP] {name}: too few val/test hyperedges")
        return {}

    H, sizes = build_incidence(n, train_members)
    A, Lg, _ = graph_from_incidence(H, sizes, weight_by_inverse=True)

    k_match = _matched_k(cfg, n)
    lam_star, k_star, val_score = select_lambda_k_on_val(H, sizes, A, val_members, n, cfg, cfg.base_seed)
    print(f"[{name}] n_train={n} k_matched={k_match} selected(lambda*={lam_star}, k*={k_star}) val_bestF1={val_score:.3f}")

    # Precompute hypergraph Laplacians we need
    L_hyp0, _, _ = hypergraph_laplacian(H, sizes, penalty=0.0)
    L_hyp_star, _, _ = hypergraph_laplacian(H, sizes, penalty=lam_star)

    rows = []

    def add(method, seed, labels, kmode, lam):
        m = metrics.evaluate(labels, test_members, n, A=A)
        m.update(dict(dataset=name, method=method, seed=seed, k_mode=kmode, lam=lam))
        rows.append(m)

    for seed in range(cfg.base_seed, cfg.base_seed + cfg.n_seeds):
        # ---- MATCHED-k spectral methods (all share k_match) ----
        add("graph_spectral", seed, spectral.spectral_cluster(Lg, k_match, seed), "matched", 0.0)
        add("hyper_zhou_lam0", seed, spectral.spectral_cluster(L_hyp0, k_match, seed), "matched", 0.0)
        add("hyper_penalized", seed, spectral.spectral_cluster(L_hyp_star, k_match, seed), "matched", lam_star)

        # ---- ADAPTIVE-k spectral (each picks own k via eigengap) ----
        lab_g, kg, _, _ = spectral.spectral_cluster_adaptive(Lg, A, "eigengap", kmax=cfg.kmax, seed=seed)
        add("graph_spectral_adaptive", seed, lab_g, "adaptive", 0.0)
        lab_h, kh, _, _ = spectral.spectral_cluster_adaptive(L_hyp_star, A, "eigengap", kmax=cfg.kmax, seed=seed)
        add("hyper_penalized_adaptive", seed, lab_h, "adaptive", lam_star)

    # ---- Non-spectral baselines (deterministic-ish; run a few seeds) ----
    for seed in range(cfg.base_seed, cfg.base_seed + min(cfg.n_seeds, 5)):
        if cfg.run_louvain:
            add("louvain", seed, baselines.louvain_labels(A, seed), "native", 0.0)
        if cfg.run_leiden:
            lab = baselines.leiden_labels(A, seed)
            if lab is not None:
                add("leiden", seed, lab, "native", 0.0)
        if cfg.run_mcl and seed == cfg.base_seed:
            lab = baselines.mcl_labels(A)
            if lab is not None:
                add("mcl", seed, lab, "native", 0.0)

    results = pd.DataFrame(rows)
    results.to_csv(out_dir / f"{name}_test_metrics.csv", index=False)

    # ---- Randomized controls (matched-k, best hypergraph operator) ----
    control_rows = []
    if cfg.run_controls:
        for cseed in range(cfg.n_control_seeds):
            # hyperedge-size-preserving null: rebuild operators from randomized edges
            rnd_edges = controls.hyperedge_size_preserving(tr, seed=cseed)
            rnd_members = to_index_members(rnd_edges, pid)
            if len(rnd_members) >= 2:
                Hr, szr = build_incidence(n, rnd_members)
                Lr, _, _ = hypergraph_laplacian(Hr, szr, penalty=lam_star)
                lab = spectral.spectral_cluster(Lr, k_match, cfg.base_seed)
                mrow = metrics.evaluate(lab, test_members, n, A=A)
                mrow.update(dict(dataset=name, control="hyperedge_size_preserving", seed=cseed))
                control_rows.append(mrow)
            # degree-preserving graph rewire
            Ar = controls.degree_preserving_rewire(A, seed=cseed)
            degr = np.asarray(Ar.sum(axis=1)).ravel(); degr = np.where(degr > 0, degr, 1.0)
            Disq = sparse.diags(1.0 / np.sqrt(degr))
            Lgr = sparse.identity(n, format="csr") - Disq @ Ar @ Disq
            lab = spectral.spectral_cluster(Lgr.tocsr(), k_match, cfg.base_seed)
            mrow = metrics.evaluate(lab, test_members, n, A=A)
            mrow.update(dict(dataset=name, control="degree_preserving_rewire", seed=cseed))
            control_rows.append(mrow)
    controls_df = pd.DataFrame(control_rows)
    if not controls_df.empty:
        controls_df.to_csv(out_dir / f"{name}_controls.csv", index=False)

    # ---- GO enrichment for the two key matched-k partitions ----
    enr_rows = []
    per_term_frames = []
    if prot2go:
        for method, L in [("graph_spectral", Lg), ("hyper_penalized", L_hyp_star)]:
            lab = spectral.spectral_cluster(L, k_match, cfg.base_seed)
            edf = enrich_partition(proteins, lab, prot2go, method, name,
                                   min_term_size=cfg.min_term_size)
            summ = enrichment_summary(proteins, lab, edf, q_thr=cfg.q_threshold)
            summ.update(dict(dataset=name, method=method))
            enr_rows.append(summ)
            if edf is not None and not edf.empty:
                per_term_frames.append(edf)
            # save cluster membership so REVIGO / case studies can be reproduced
            pd.DataFrame({"protein": proteins, "cluster": lab, "method": method,
                          "dataset": name}).to_csv(
                out_dir / f"{name}_{method}_clusters.csv", index=False)
    enr_df = pd.DataFrame(enr_rows)
    if not enr_df.empty:
        enr_df.to_csv(out_dir / f"{name}_go_enrichment_summary.csv", index=False)
    # NEW: full per-term enrichment table (inputs to GO dot-plots and REVIGO)
    if per_term_frames:
        pd.concat(per_term_frames, ignore_index=True).to_csv(
            out_dir / f"{name}_go_terms.csv", index=False)

    manifest = dict(dataset=name, mitab_path=str(path), **load_stats, **split_stats,
                    n_train_proteins=n, n_train_hyperedges=len(train_members),
                    n_val_hyperedges=len(val_members), n_test_hyperedges=len(test_members),
                    lambda_star=lam_star, k_star=k_star, k_matched=k_match,
                    go_coverage=float(len(prot2go) / n) if n else 0.0)
    pd.DataFrame([manifest]).to_csv(out_dir / f"{name}_manifest.csv", index=False)

    return dict(results=results, controls=controls_df, enrichment=enr_df,
                manifest=pd.DataFrame([manifest]))


def run_all(cfg: Config):
    t0 = time.time()
    all_results, all_manifests = [], []
    for name in cfg.dataset_names:
        try:
            out = run_dataset(name, cfg)
            if out:
                all_results.append(out["results"])
                all_manifests.append(out["manifest"])
        except Exception as e:
            print(f"[ERROR] {name}: {type(e).__name__}: {e}")
    root = Path(cfg.out_root)
    root.mkdir(parents=True, exist_ok=True)
    if all_results:
        pd.concat(all_results, ignore_index=True).to_csv(root / "ALL_test_metrics.csv", index=False)
    if all_manifests:
        pd.concat(all_manifests, ignore_index=True).to_csv(root / "ALL_manifest.csv", index=False)
    print(f"[DONE] {len(all_results)} datasets in {time.time() - t0:.1f}s -> {root}")
