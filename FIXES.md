# Bug fixes for the Cancer and Diabetes failures

## 1. Cancer — `ArpackNoConvergence` in the eigensolver
**Symptom:** `ARPACK error -1: No convergence (54811 iterations, 15/16
eigenvectors converged)` on n_train=5481.

**Cause:** `hgspectral/spectral.py::smallest_eigs` used
`eigsh(L, k, which="SM")`. Smallest-magnitude ARPACK is notoriously slow to
converge on large, spectrally-clustered normalized Laplacians. The fallback
`eigsh(..., sigma=1e-6, which="LM")` sits almost exactly on the Laplacian's zero
eigenvalue, so `(L − σI)` is near-singular and shift-invert is ill-conditioned —
it can fail too.

**Fix:** shift-invert with σ strictly **below** zero. Because L is PSD,
`(L − σI) = L + |σ| I` is positive definite and well conditioned, and
`which="LM"` then returns the eigenvalues of L nearest σ (i.e. the smallest). We
try progressively stronger shifts (−1e-3, −1e-2, −1e-1, −1), then fall back to
`which="SA"` with a large iteration budget, then LOBPCG, then a dense solve.
Verified: the exact n=5481 case now converges in ~0.2 s.

## 2. Diabetes — `NetworkXAlgorithmError` in the rewire null
**Symptom:** `Maximum number of swap attempts (134021000) exceeded before desired
swaps achieved (1340210)` in `degree_preserving_rewire`.

**Cause:** two issues. (a) The dense/degree-heterogeneous co-complex graph makes
achieving 10×m successful double-edge swaps effectively impossible within any
finite budget. (b) The handler caught only `nx.NetworkXError`, but
`double_edge_swap` raises `NetworkXAlgorithmError`, which is a **sibling** of
`NetworkXError` (both subclass `NetworkXException`), so the exception escaped.

**Fix:** catch `nx.NetworkXException` (covers both), bound the try budget
(`max_tries = min(50·target, 2e6)`), and **accept the partially-rewired graph** —
a valid approximate degree-preserving (Maslov–Sneppen) null; double-edge swaps
preserve the degree sequence exactly whether the target count is reached or not.
Verified: near-complete and heterogeneous graphs no longer crash and the degree
sequence is preserved.

## Re-running the two datasets
With the fixed code and the raw MITAB files in place
(`<root>/Cancer/Cancer.txt`, `<root>/Diabetes/Diabetes.txt`):

```bash
PYTHONPATH=. python run_pipeline.py --datasets Cancer Diabetes \
    --search-roots HG_Human --seeds 25 --out fin_result
PYTHONPATH=. python analyze_results.py --results fin_result/ALL_test_metrics.csv
PYTHONPATH=. python make_paper_figures.py
```

The pipeline now also writes `<dataset>_go_terms.csv` and
`<dataset>_<method>_clusters.csv`, which are the inputs to the REVIGO scatter and
GO dot-plot in `go_analysis.py`.
