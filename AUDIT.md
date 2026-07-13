# Audit & Q1-Upgrade Plan
### "Uncovering Higher-Order Functional Modules: Graph vs. Hypergraph Clustering"

This document audits the four project files, lists required fixes, gives a
revised experimental plan, and ends with a Q1-readiness checklist. It is
paired with the `hgspectral/` pipeline, which implements the fixes.

> **Reproducibility caveat.** No MITAB `.txt` files or CORUM file were provided,
> so the real numbers could not be regenerated. Per the project rule, no results
> were invented. All result cells in the revised manuscript are placeholders to
> be filled from CSVs produced by the pipeline on your own data.

---

## 1. What is already implemented

**`Final_Mintab_spectral-mod.ipynb`** (the manuscript's apparent source):
MITAB parsing, top-3 taxid filter, pull-down hyperedge construction by
experiment key, dense Zhou hypergraph Laplacian with penalty `λI` / `λDv⁻¹`,
clique-expansion graph Laplacian, dense-eigendecomposition spectral clustering,
metrics (pairwise F1, best-match F1, B³, NMI, ARI, modularity, conductance),
hypergeometric GO enrichment + BH FDR, keyword-based BP/MF/CC split,
Jaccard-based "REVIGO-like" reduction, radar/scree plotting.

**`research_approach.ipynb`** (a more advanced, unreferenced pipeline): sparse
operators + `eigsh`, **hyperedge train/test split**, **stability-based λ
selection on train**, held-out test evaluation, a real GO DAG loader
(`go-basic.obo`) with Lin semantic similarity, and an MGSA MCMC posterior. This
notebook is substantially stronger than what the manuscript describes and is the
right foundation — but it is not the pipeline the paper reports.

---

## 2. Methodological weaknesses (ordered by severity)

1. **Unmatched k presented as matched.** `tab:mitab_best_structural` reports
   graph at k = 2–5 and hypergraph at k = 11–13. The narrative says k is "held
   fixed for both." Extra clusters mechanically raise pairwise-F1/NMI, so the
   "3–4×" gains are confounded by cluster count. **This alone can sink the paper
   in review.**
2. **Results not reproducible from the code.** The committed
   `Final_Mintab` uses `force_k = n // 5` for both operators, which produces
   neither the k=2 nor k=13 values in the table. The reported numbers come from
   an older eigengap-per-model variant that is not in the repo.
3. **Selection on the test set.** "Best hypergraph over λ" is chosen using the
   same structural metric that is then reported. There is no validation split
   for λ (or k) in the manuscript pipeline; this is optimistic bias.
4. **Ground truth = input.** The "structural ground truth" (complexes) is
   identical to the hyperedge set used to build both operators. Recovering your
   own input is not held-out evaluation. No train/test complex split in the
   reported pipeline.
5. **Single seed, no dispersion, no tests.** One `SEED=42`; no repeats, no
   confidence intervals, no significance testing, no effect sizes.
6. **GO leakage / circularity.** GO terms are parsed from the same MITAB xref
   column used to define proteins, then used as "functional validation." At
   minimum this must be stated; ideally validate against an independent GO
   release.
7. **Ontology split is a keyword heuristic.** `infer_go_ontology` guesses
   BP/MF/CC from the term name. This is not the GO DAG and will mislabel terms.
   (The `research_approach` notebook already has the real DAG loader — use it.)
8. **"REVIGO-like" is not semantic.** It uses Jaccard of protein sets, not
   ontology-based semantic similarity. Call it "annotation-overlap reduction,"
   not REVIGO.
9. **Scalability.** `pairwise_F1` enumerates all O(n²) pairs in Python; for
   Cancer (n≈6749) this is ~22M iterations and will effectively hang. Dense
   `eigh` on n×n is also infeasible at that scale.
10. **Graph baseline is weak and singular.** Only clique-expansion spectral is
    compared. No Louvain/Leiden/MCL/ClusterONE, and no unregularized-vs-penalized
    hypergraph ablation isolating the effect of λ.
11. **No randomized controls.** Nothing shows the signal exceeds a null model.

---

## 3. Places where claims are stronger than the evidence

- Abstract: "hypergraph … **consistently** matches or improves." Under matched-k
  and held-out evaluation this is unproven; likely conditional.
- "pairwise and B³ F1 scores increase by **≈2–4×**." Confounded by k mismatch.
- "converting graph clusterings with **no** significant GO enrichment into
  hypergraph clusterings where **most** proteins belong to enriched modules."
  Enrichment fraction depends on k and cluster sizes; not attributable to
  representation alone as stated.
- Conclusion: "hypergraph spectral clustering is a **natural and effective
  extension** … particularly [for] intrinsically higher-order [designs]." Keep
  the mechanism claim; drop the universal-superiority framing.

---

## 4. Required-fixes checklist (maps to the pipeline)

- [x] Matched-k comparison (all methods same k) — `pipeline.run_dataset`, k_mode="matched".
- [x] Adaptive-k comparison (eigengap/silhouette/modularity) — `spectral.spectral_cluster_adaptive`.
- [x] Nested λ,k selection on validation only — `pipeline.select_lambda_k_on_val`.
- [x] ≥20 seeds — `Config.n_seeds` (default 20/25).
- [x] Split by hyperedges, size-stratified — `splits.split_hyperedges`.
- [x] Held-out hyperedge recovery metric — `metrics.hyperedge_recovery`.
- [x] Statistical tests + effect sizes + bootstrap CIs — `stats.py`, `analyze_results.py`.
- [x] Randomized controls (rewire, size-preserving, GO-perm, taxid-null) — `controls.py`.
- [x] Baselines: graph spectral, Zhou λ=0, penalized, Louvain, Leiden, MCL — `baselines.py`, pipeline.
- [x] λ grid {0,1e-4,1e-3,1e-2,1e-1,1,10} + λ=0 ablation — `Config.lambda_grid`.
- [x] Vectorized/sparse metrics for scale — `metrics.py`, sparse `operators.py`.
- [x] BH-FDR GO enrichment via scipy hypergeom — `enrichment.py`.
- [x] Config files, CSV outputs, PDF+PNG figures, data manifest — pipeline + analyze.
- [x] Unit tests — `tests/test_pipeline.py` (9 passing).
- [x] environment.yml / Dockerfile / README / archival plan.
- [ ] ClusterONE overlapping baseline — wrapper provided; run externally (needs JAR).
- [ ] Real GO DAG ontology split + semantic similarity — loader exists in
      `research_approach.ipynb`; wire into `enrichment.py` (marked *proposed*).
- [ ] Hyperedge weighting ablation (uniform / inverse-size / confidence /
      evidence-count) — `operators.hypergraph_laplacian(weights=...)` supports
      it; add the sweep to the runner (*proposed*).
- [ ] Min-complex-size sweep {2,3,4,5} and top-1/top-3/all species sweep —
      one-line loops over `Config` (*proposed*, cheap).
- [ ] Runtime/memory table — time each method in `pipeline` (*proposed*).

---

## 5. Revised experimental design (what to run)

For each dataset, hold `split_seed` fixed and:

1. **Dataset summary** from `*_manifest.csv` (proteins, hyperedges, size
   distribution, GO coverage, top taxids).
2. **Matched-k** (Table 1): graph_spectral vs hyper_zhou_lam0 vs
   hyper_penalized vs Louvain vs Leiden vs MCL, all at the same k, mean ± SD
   over ≥20 seeds, on held-out test hyperedges.
3. **Adaptive-k** (Table 2): each method's best practical k (eigengap; report
   silhouette/modularity variants in supplement).
4. **λ sensitivity** (Figure): held-out best-match F1 vs λ, per dataset, with
   λ chosen on validation marked.
5. **Ablations**: λ=0 vs λ\*; weighting schemes; min-size {2..5}; species scope.
6. **Randomized controls** (Table 3): real vs degree-preserving-rewire vs
   size-preserving null; Δ with bootstrap CI.
7. **GO enrichment** (Table 4): enriched cluster/protein fractions, graph vs
   hypergraph, BH-FDR; graph-only vs hypergraph-only terms.
8. **Statistics**: paired Wilcoxon (graph vs hypergraph) per dataset + pooled,
   with rank-biserial and bootstrap CI; Friedman + Nemenyi across ≥3 methods
   over datasets (needs ≥3 datasets).
9. **Case studies**: 1–2 hypergraph-only enriched modules from a strong dataset
   (e.g. Cancer/Coronavirus), described qualitatively.
10. **Cost**: runtime & peak memory per method per dataset.

---

## 6. Concrete code changes (summary)

- Replace dense `eigh` with sparse `eigsh` (done in `spectral.py`).
- Replace O(n²) `pairwise_F1` loop with cluster-blocked sparse counting
  (`metrics.pairwise_prf1`).
- Move λ/k selection to a validation split (`pipeline.select_lambda_k_on_val`);
  never touch test during selection.
- Force identical k across methods in the matched condition.
- Add Louvain/Leiden/MCL, controls, multi-seed loop, statistics, manifests.
- Keep the authors' MITAB/hyperedge logic verbatim for behavioural parity
  (`mitab.py`).

---

## 7. Q1-readiness checklist

| Requirement | Status | Where |
|---|---|---|
| Fair (matched-k) comparison | Implemented | `pipeline`, Table 1 |
| Adaptive-k "best practical" comparison | Implemented | `spectral_cluster_adaptive` |
| Held-out validation (split by hyperedge) | Implemented | `splits`, `pipeline` |
| Nested λ,k selection (no test leakage) | Implemented | `select_lambda_k_on_val` |
| ≥20 seeds, dispersion reported | Implemented | `Config.n_seeds` |
| Strong baselines (Louvain/Leiden/MCL) | Implemented | `baselines` |
| Unregularized vs penalized ablation | Implemented | λ grid incl. 0 |
| Randomized controls | Implemented | `controls` |
| Statistical tests + effect sizes + CIs | Implemented | `stats`, `analyze_results` |
| Vectorized/sparse for scale | Implemented | `operators`, `metrics` |
| GO enrichment + BH FDR | Implemented | `enrichment` |
| Real GO DAG ontology split + semantic sim | Proposed | reuse `research_approach` loader |
| Overlapping baseline (ClusterONE) | Proposed/external | `baselines.clusterone_labels` |
| Weighting / size / species ablations | Proposed (cheap) | `Config` loops |
| Runtime/memory profile | Proposed | add timers |
| Reproducible pipeline + config + tests | Implemented | repo root |
| Docker/conda + README + archival plan | Implemented | repo root |
| Claim-safe manuscript | Implemented | `manuscript_revisions.tex` |
| Data manifest | Implemented | `*_manifest.csv` |
| Case studies | Pending data | fill from real runs |

**Bottom line:** the machinery for a Q1 submission is in place and tested. The
remaining gap is running it on the real MITAB/CORUM files and filling the
manuscript's result placeholders — plus the small set of items marked
*proposed* above.
