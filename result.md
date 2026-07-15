# Results

This document summarizes the evaluation of graph and hypergraph spectral clustering across nine IntAct thematic MITAB datasets. Unless otherwise stated, performance is measured using the **symmetric best-match \(F_1\)** score on held-out test hyperedges. Matched-\(k\) and adaptive-\(k\) spectral results are averaged over **25 seeds**. The regularization parameter \(\lambda\) and the matched number of clusters \(k\) were selected using validation data only.

## Main findings

- The **penalized hypergraph** significantly improved held-out recovery on four datasets: **Chromatin, Cyanobacteria, Affinomics, and Crohn's disease**.
- Graph and penalized-hypergraph performance was statistically indistinguishable on **BioCreative, Cancer, and Coronavirus**.
- The penalized hypergraph was significantly lower on **Cardiac and Diabetes**, although the absolute reductions were small: approximately \(0.003\) and \(0.007\) \(F_1\), respectively.
- The unregularized Zhou hypergraph operator was not consistently competitive. The observed hypergraph gains came primarily from the degree-aware penalty rather than from the hypergraph representation alone.
- Hypergraph performance depended strongly on controlling cluster granularity. Under independently selected adaptive \(k\), the graph frequently outperformed the hypergraph because the hypergraph eigengap often selected an unsuitable number of clusters.
- Real data substantially outperformed degree-preserving and hyperedge-size-preserving null models on every dataset except Cancer.
- Louvain and Markov clustering often achieved higher absolute recovery because they produced many more communities. This shows that **cluster granularity affects absolute complex recovery more strongly than the graph-versus-hypergraph operator choice**.
- Both graph and hypergraph modules were strongly GO-enriched on most datasets.

---

## 1. Dataset summary

The processed datasets span almost three orders of magnitude in size. The training splits contain between **71 and 5,481 proteins** and between **87 and 4,171 hyperedges**. Gene Ontology coverage is high throughout, ranging from **0.86 to 0.97**.

| Dataset | Filtered PPI records | Training proteins | Train / validation / test hyperedges | \(\lambda^\star\) | Matched \(k\) | GO coverage |
|---|---:|---:|---:|---:|---:|---:|
| Crohn's disease | 249 | 71 | 93 / 31 / 30 | 1 | 14 | 0.86 |
| Cyanobacteria | 575 | 198 | 87 / 28 / 29 | 0.001 | 15 | 0.97 |
| Affinomics | 2,148 | 326 | 254 / 85 / 85 | 0 | 15 | 0.86 |
| BioCreative | 870 | 327 | 309 / 103 / 103 | 0.001 | 15 | 0.95 |
| Diabetes | 1,053 | 658 | 134 / 45 / 43 | 1 | 15 | 0.97 |
| Cardiac | 4,172 | 1,152 | 889 / 296 / 296 | 0 | 15 | 0.97 |
| Chromatin | 8,546 | 1,467 | 847 / 282 / 283 | 10 | 15 | 0.96 |
| Coronavirus | 13,935 | 3,045 | 1,534 / 511 / 511 | 0 | 15 | 0.94 |
| Cancer | 33,657 | 5,481 | 4,171 / 1,391 / 1,390 | 0 | 15 | 0.91 |

The validation-selected penalty varies from \(0\) to \(10\). Four datasets select no regularization, whereas Chromatin selects the strongest tested penalty, \(\lambda^\star=10\).

---

## 2. Matched-\(k\) graph-versus-hypergraph comparison

The following table reports mean held-out symmetric best-match \(F_1\), with standard deviation in parentheses, over 25 seeds.

| Dataset | Graph spectral | Hypergraph, \(\lambda=0\) | Penalized hypergraph | Penalized minus graph |
|---|---:|---:|---:|---:|
| Affinomics | 0.182 (0.006) | **0.203 (0.007)** | **0.203 (0.007)** | +0.022 |
| BioCreative | 0.192 (0.005) | 0.184 (0.005) | **0.194 (0.005)** | +0.002 |
| Cancer | **0.049 (0.001)** | 0.048 (0.002) | 0.049 (0.001) | -0.000 |
| Cardiac | **0.115 (0.003)** | 0.113 (0.004) | 0.112 (0.003) | -0.003 |
| Chromatin | 0.207 (0.002) | 0.209 (0.002) | **0.244 (0.004)** | +0.037 |
| Coronavirus | **0.242 (0.006)** | 0.241 (0.004) | 0.241 (0.003) | -0.001 |
| Crohn's disease | 0.532 (0.015) | 0.526 (0.011) | **0.544 (0.005)** | +0.012 |
| Cyanobacteria | 0.471 (0.005) | 0.473 (0.001) | **0.494 (0.000)** | +0.023 |
| Diabetes | **0.241 (0.012)** | 0.238 (0.011) | 0.234 (0.006) | -0.007 |

The paired comparison separates the datasets into three regimes:

1. **Significant hypergraph improvements:** Chromatin, Cyanobacteria, Affinomics, and Crohn's disease.
2. **No statistically detectable difference:** BioCreative, Cancer, and Coronavirus.
3. **Significant but small absolute decrements:** Cardiac and Diabetes.

### Paired statistical tests

Two-sided paired Wilcoxon signed-rank tests were applied to the 25 seed-matched scores. Confidence intervals are paired percentile-bootstrap intervals for the mean difference, computed from 20,000 bootstrap samples.

| Dataset | Mean \(\Delta F_1\) | 95% bootstrap CI | Wilcoxon \(p\) | Rank-biserial \(r\) | Interpretation |
|---|---:|---:|---:|---:|---|
| Affinomics | +0.0217 | [0.0185, 0.0251] | \(5.96\times10^{-8}\) | +1.000 | Hypergraph better |
| BioCreative | +0.0019 | [-0.0016, 0.0051] | 0.210 | +0.292 | No clear difference |
| Cancer | -0.0003 | [-0.0008, 0.0003] | 0.411 | -0.194 | No clear difference |
| Cardiac | -0.0025 | [-0.0038, -0.0012] | 0.0023 | -0.674 | Graph better; small absolute gap |
| Chromatin | +0.0370 | [0.0352, 0.0384] | \(5.96\times10^{-8}\) | +1.000 | Hypergraph better |
| Coronavirus | -0.0013 | [-0.0041, 0.0016] | 0.491 | -0.163 | No clear difference |
| Crohn's disease | +0.0116 | [0.0058, 0.0176] | 0.0004 | +0.809 | Hypergraph better |
| Cyanobacteria | +0.0226 | [0.0203, 0.0243] | \(6.25\times10^{-6}\) | +1.000 | Hypergraph better |
| Diabetes | -0.0067 | [-0.0118, -0.0015] | 0.0067 | -0.606 | Graph better; small absolute gap |

[Open the paired-effect forest plot (PDF)](figures/F1_forest_plot_hyper_vs_graph.pdf)

---

## 3. The penalty drives the hypergraph gains

The unregularized Zhou operator,

\[
L_{\mathrm{hyp}}^{(0)},
\]

was never the unique best method for a dataset. The improved model uses the degree-aware penalty

\[
L_{\mathrm{hyp}}(\lambda)
=
L_{\mathrm{hyp}}^{(0)}
+
\lambda D_v^{-1}.
\]

Across the nine dataset-level mean scores, the average ranks were:

| Spectral operator | Mean rank |
|---|---:|
| Penalized hypergraph | **1.72** |
| Graph spectral | 1.89 |
| Unregularized Zhou hypergraph | 2.39 |

Although the penalized hypergraph has the best mean rank, the omnibus Friedman test is not significant:

\[
\chi^2_F=2.23,\qquad p=0.328.
\]

Therefore, the rankings are consistent with a possible advantage for the penalized hypergraph, but they do not establish a universal global ordering. The evidence is strongest at the individual-dataset level.

[Open the critical-difference diagram (PDF)](figures/F3_critical_difference.pdf)

---

## 4. When the hypergraph helps

The largest benefits occur on the smaller, more group-structured datasets—Crohn's disease, Cyanobacteria, and Affinomics—and on Chromatin, where validation selected a strong penalty of \(\lambda^\star=10\).

The two largest datasets show little separation between the operators:

- **Coronavirus:** graph \(F_1=0.242\), penalized hypergraph \(F_1=0.241\).
- **Cancer:** graph \(F_1=0.049\), penalized hypergraph \(F_1=0.049\) after rounding.

At a fixed small \(k\), thousands of proteins are compressed into only 15 clusters. Consequently, neither representation can adequately recover the much finer gold-standard complex structure.

Across the nine datasets, the penalized-hypergraph improvement has:

- a descriptive negative association with training protein count: Spearman \(\rho=-0.37\), \(p=0.332\);
- a descriptive positive association with \(\lambda^\star\): Spearman \(\rho=+0.35\), \(p=0.354\).

Because \(N=9\), neither correlation is statistically significant. They should be interpreted as exploratory evidence rather than as confirmed cross-dataset relationships.

---

## 5. Adaptive-\(k\) results

When graph and hypergraph methods independently select \(k\) using the eigengap heuristic, the graph frequently obtains a higher \(F_1\). The main problem is not necessarily the hypergraph operator itself: its eigengap often selects a very small and unsuitable number of clusters.

| Dataset | Adaptive graph \(F_1\) | Graph \(k\) | Adaptive hypergraph \(F_1\) | Hypergraph \(k\) |
|---|---:|---:|---:|---:|
| Affinomics | 0.099 (0.001) | 6 | 0.076 (0.003) | 4 |
| BioCreative | 0.129 (0.003) | 9 | 0.101 (0.001) | 6 |
| Cancer | 0.039 (0.001) | 9.12 | 0.041 (0.001) | 9.12 |
| Cardiac | 0.096 (0.003) | 11 | 0.095 (0.003) | 11 |
| Chromatin | 0.200 (0.008) | 12.68 | 0.104 (0.000) | 2 |
| Coronavirus | 0.208 (0.004) | 8 | 0.167 (0.025) | 8 |
| Crohn's disease | 0.530 (0.010) | 14 | 0.297 (0.000) | 6 |
| Cyanobacteria | 0.160 (0.000) | 2 | **0.492 (0.002)** | 14 |
| Diabetes | 0.240 (0.015) | 15 | 0.064 (0.000) | 3 |

Cancer and Chromatin occasionally selected different \(k\) values across seeds, so their table entries show the mean selected \(k\). Cyanobacteria is the major exception to the general adaptive-\(k\) pattern: the graph selects \(k=2\), whereas the hypergraph selects \(k=14\) and retains strong performance.

These results show that the observed hypergraph advantage is specific to the **controlled matched-\(k\) comparison**. Adaptive-\(k\) results are a caveat and should not be used as evidence of a general hypergraph advantage.

---

## 6. Randomized controls

The penalized hypergraph on real data substantially exceeds both null models on all datasets except Cancer. Null values are means over 10 randomizations per control type.

| Dataset | Real penalized hypergraph | Degree-preserving rewire null | Hyperedge-size-preserving null |
|---|---:|---:|---:|
| Affinomics | **0.203** | 0.114 | 0.116 |
| BioCreative | **0.194** | 0.099 | 0.108 |
| Cancer | 0.049 | 0.049 | 0.034 |
| Cardiac | **0.112** | 0.047 | 0.047 |
| Chromatin | **0.244** | 0.080 | 0.088 |
| Coronavirus | **0.241** | 0.066 | 0.074 |
| Crohn's disease | **0.544** | 0.305 | 0.315 |
| Cyanobacteria | **0.494** | 0.173 | 0.176 |
| Diabetes | **0.234** | 0.061 | 0.053 |

For example:

- Chromatin improves from \(0.080/0.088\) under the two nulls to \(0.244\) on real data.
- Cyanobacteria improves from \(0.173/0.176\) to \(0.494\).
- Coronavirus improves from \(0.066/0.074\) to \(0.241\).

Cancer is the exception. Its real score of approximately \(0.049\) is essentially equal to the degree-preserving null. This indicates a **floor effect** at \(k=15\), not a tie between two high-quality solutions.

[Open the real-versus-null control figure (PDF)](figures/F6_real_vs_null_controls.pdf)

---

## 7. Baselines with independently selected granularity

Louvain and Markov clustering (MCL) produce their own numbers of communities rather than using the matched spectral \(k\). Their scores are therefore useful as absolute-recovery references, but they are not controlled operator comparisons.

Louvain results are averaged over five runs. MCL has one recorded run per dataset.

| Dataset | Penalized hypergraph, matched \(k\) | Louvain \(F_1\) | Louvain communities | MCL \(F_1\) | MCL communities |
|---|---:|---:|---:|---:|---:|
| Affinomics | 0.203 | 0.320 (0.006) | 54.2 | **0.410** | 116 |
| BioCreative | 0.194 | 0.484 (0.000) | 94 | **0.523** | 107 |
| Cancer | 0.049 | 0.083 (0.002) | 137.6 | **0.237** | 3,426 |
| Cardiac | 0.112 | 0.215 (0.003) | 63 | **0.319** | 574 |
| Chromatin | 0.244 | 0.297 (0.019) | 46.8 | **0.371** | 742 |
| Coronavirus | 0.241 | **0.260 (0.005)** | 32 | 0.188 | 2,164 |
| Crohn's disease | 0.544 | 0.528 (0.000) | 13 | **0.552** | 16 |
| Cyanobacteria | 0.494 | **0.506 (0.000)** | 18 | 0.365 | 110 |
| Diabetes | 0.234 | 0.241 (0.000) | 20 | **0.284** | 581 |

Best-match \(F_1\) rewards partitions whose granularity resembles the often-small gold complexes. The spectral methods use only 14 or 15 clusters, whereas MCL can create hundreds or thousands. The higher baseline scores therefore emphasize that absolute recovery depends heavily on the number and sizes of predicted modules.

The contribution of this experiment is a **controlled within-spectral comparison** between graph and hypergraph operators, not a claim that fixed-\(k\) spectral clustering is the strongest available protein-complex detection method.

---

## 8. GO enrichment

Both representations produce strongly enriched modules. The table reports the fraction of clusters and proteins associated with at least one significant GO enrichment result at \(q\leq0.05\), using a hypergeometric test with Benjamini-Hochberg correction.

| Dataset | Graph enriched clusters | Hypergraph enriched clusters | Graph enriched proteins | Hypergraph enriched proteins |
|---|---:|---:|---:|---:|
| Affinomics | 1.000 | 0.933 | 1.000 | 0.979 |
| BioCreative | 1.000 | 1.000 | 1.000 | 1.000 |
| Cancer | 1.000 | 1.000 | 1.000 | 1.000 |
| Cardiac | 1.000 | 1.000 | 1.000 | 1.000 |
| Chromatin | 1.000 | 1.000 | 1.000 | 1.000 |
| Coronavirus | 1.000 | 1.000 | 1.000 | 1.000 |
| Crohn's disease | 0.357 | 0.429 | 0.535 | 0.577 |
| Cyanobacteria | 0.667 | 0.533 | 0.934 | 0.909 |
| Diabetes | 0.800 | 0.667 | 0.979 | 0.956 |

Six datasets have at least 90% enriched clusters for both graph and hypergraph partitions. Crohn's disease, Cyanobacteria, and Diabetes have lower or mixed cluster-level enrichment, although protein-level enrichment remains high for Cyanobacteria and Diabetes.

Functional coherence therefore does not clearly separate the two operators. The main effect of the operator appears to concern **which proteins are grouped together**, rather than whether the resulting modules contain significant biological annotations.

[Open the graph-versus-hypergraph GO enrichment figure (PDF)](figures/F4_go_enrichment_comparison.pdf)

---

## 9. Semantic content of recovered modules

A REVIGO-style semantic-space reduction was used to summarize enriched GO terms. The available result files contain GO identifiers but not the GO directed acyclic graph. Therefore, term similarity in the reported analysis is based on **module co-occurrence**: terms are considered similar when they are enriched in the same recovered modules.

For Chromatin, penalized-hypergraph modules form coherent groups represented by terms associated with:

- the nucleoplasm;
- chromatin remodeling;
- ribosomal and protein-synthesis-related assemblies.

These categories are biologically compatible with chromatin-associated complexes. A GO-DAG-based Lin or Resnik analysis can be run when `go-basic.obo` is supplied. Namespace-resolved Biological Process, Molecular Function, and Cellular Component summaries should likewise be generated from the GO ontology rather than from keyword heuristics.

[Open the Chromatin semantic scatter plot (PDF)](figures/R1_revigo_Chromatin_hyper_penalized.pdf)

---

## 10. Statistical interpretation

The complete evidence supports a **conditional** conclusion:

- Four datasets show reproducible improvements from the penalized hypergraph.
- Three datasets provide no statistically detectable graph-versus-hypergraph difference.
- Two datasets favor the graph by small absolute \(F_1\) margins.
- The cross-dataset Friedman test does not establish a universal method ordering.
- The unregularized hypergraph alone does not explain the gains.
- The matched-\(k\) advantage does not reliably transfer to adaptive-\(k\) selection.
- Null-model and GO analyses confirm that most recovered modules contain non-random biological structure, even when graph and hypergraph recovery scores are similar.

Accordingly, the results support the claim that **degree-aware hypergraph regularization can improve complex recovery for particular interaction structures**, not that hypergraph spectral clustering is universally superior.

---

## 11. Run completeness and source files

All nine dataset runs completed successfully. The recorded end-to-end runtime was approximately **2,821 seconds (47.0 minutes)**.
