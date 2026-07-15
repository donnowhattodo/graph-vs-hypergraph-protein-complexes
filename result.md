# Results

This study evaluated graph and hypergraph spectral clustering on nine thematic IntAct MITAB datasets. Unless otherwise stated, performance is reported as the **symmetric best-match $F_1$ score** on held-out test hyperedges. Results for the matched-$k$ and adaptive-$k$ spectral methods are averaged over **25 random seeds**. The regularization parameter $\lambda$ and the matched number of clusters $k$ were selected using validation data only.

## 1. Overview of the findings

The penalized hypergraph method produced significant improvements on **Chromatin, Cyanobacteria, Affinomics, and Crohn's disease**. It was statistically indistinguishable from graph spectral clustering on **BioCreative, Cancer, and Coronavirus**. Graph spectral clustering performed significantly better on **Cardiac and Diabetes**, although the absolute differences were small.

The results also show that the improvement was not caused by the hypergraph representation alone. The unregularized Zhou hypergraph operator was not consistently competitive, whereas the degree-aware penalty produced the strongest hypergraph results. In addition, the hypergraph advantage was observed primarily when graph and hypergraph methods were evaluated using the same number of clusters. When each method selected its own $k$ using the eigengap heuristic, the graph method frequently performed better because the hypergraph method often selected an unsuitable number of clusters.

Real-data performance exceeded both randomized controls on almost every dataset, and both graph and hypergraph solutions produced biologically enriched modules. Overall, the evidence supports a **conditional advantage of degree-aware hypergraph regularization**, rather than a universal advantage of hypergraph spectral clustering.

---

## 2. Dataset summary

The processed datasets vary substantially in scale. The training splits contain between **71 and 5,481 proteins** and between **87 and 4,171 hyperedges**. Gene Ontology (GO) coverage is consistently high, ranging from **0.86 to 0.97**.

| Dataset | Filtered PPI records | Training proteins | Train / validation / test hyperedges | $\lambda^\star$ | Matched $k$ | GO coverage |
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

The validation-selected penalty varies from $0$ to $10$. No regularization was selected for Affinomics, Cardiac, Coronavirus, or Cancer. Chromatin selected the strongest tested penalty, $\lambda^\star=10$, indicating that the degree-aware regularization term was particularly important for this dataset.

---

## 3. Matched-$k$ comparison of graph and hypergraph methods

The matched-$k$ experiment provides the main controlled comparison because graph and hypergraph methods are evaluated using the same number of clusters. Table values are the mean symmetric best-match $F_1$ scores over 25 seeds, with standard deviations shown in parentheses.

| Dataset | Graph spectral | Hypergraph, $\lambda=0$ | Penalized hypergraph | Penalized hypergraph $-$ graph |
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

The penalized hypergraph equals or exceeds the graph score on six datasets after rounding and is lower on three. Statistical testing gives a more precise interpretation:

- **Significant improvement:** Affinomics, Chromatin, Crohn's disease, and Cyanobacteria.
- **No statistically detectable difference:** BioCreative, Cancer, and Coronavirus.
- **Significant but small decrease:** Cardiac and Diabetes.

### 3.1 Paired statistical analysis

A two-sided paired Wilcoxon signed-rank test was applied to the 25 seed-matched graph and penalized-hypergraph scores for each dataset. The confidence intervals are paired percentile-bootstrap intervals for the mean difference, based on 20,000 bootstrap samples.

| Dataset | Mean $\Delta F_1$ | 95% bootstrap CI | Wilcoxon $p$ | Rank-biserial $r$ | Interpretation |
|---|---:|---:|---:|---:|---|
| Affinomics | +0.0217 | [0.0185, 0.0251] | $5.96\times10^{-8}$ | +1.000 | Hypergraph significantly better |
| BioCreative | +0.0019 | [-0.0016, 0.0051] | 0.210 | +0.292 | No significant difference |
| Cancer | -0.0003 | [-0.0008, 0.0003] | 0.411 | -0.194 | No significant difference |
| Cardiac | -0.0025 | [-0.0038, -0.0012] | 0.0023 | -0.674 | Graph significantly better; small absolute difference |
| Chromatin | +0.0370 | [0.0352, 0.0384] | $5.96\times10^{-8}$ | +1.000 | Hypergraph significantly better |
| Coronavirus | -0.0013 | [-0.0041, 0.0016] | 0.491 | -0.163 | No significant difference |
| Crohn's disease | +0.0116 | [0.0058, 0.0176] | 0.0004 | +0.809 | Hypergraph significantly better |
| Cyanobacteria | +0.0226 | [0.0203, 0.0243] | $6.25\times10^{-6}$ | +1.000 | Hypergraph significantly better |
| Diabetes | -0.0067 | [-0.0118, -0.0015] | 0.0067 | -0.606 | Graph significantly better; small absolute difference |

Chromatin shows the largest improvement, with an average increase of $0.0370$ in $F_1$. Affinomics and Cyanobacteria also show strong and consistent improvements, with rank-biserial effect sizes of $1.0$. Crohn's disease shows a smaller but still statistically significant gain. Cardiac and Diabetes favor the graph method, but their mean absolute differences remain below $0.01$.

---

## 4. Contribution of the degree-aware penalty

The unregularized Zhou hypergraph operator was never the unique best spectral method. The penalized operator is defined as

$$
L_{\mathrm{hyp}}(\lambda)
=
L_{\mathrm{hyp}}^{(0)}
+
\lambda D_v^{-1}.
$$

Across the nine dataset-level mean scores, the average method ranks were:

| Spectral method | Mean rank |
|---|---:|
| Penalized hypergraph | **1.72** |
| Graph spectral | 1.89 |
| Unregularized Zhou hypergraph | 2.39 |

The penalized hypergraph has the best average rank, while the unregularized hypergraph has the worst. However, the Friedman test across the three spectral methods is not statistically significant:

$$
\chi_F^2=2.23, \qquad p=0.328.
$$

Therefore, the cross-dataset ranking does not establish a universal ordering of the methods. Instead, the strongest evidence comes from the paired, dataset-specific comparisons. The results indicate that the successful hypergraph improvements are primarily attributable to the degree-aware penalty rather than to replacing the graph with an unregularized hypergraph.

---

## 5. Conditions under which the hypergraph helps

The hypergraph improvements are concentrated on the smaller or more strongly group-structured datasets: Crohn's disease, Cyanobacteria, and Affinomics. A substantial improvement is also observed for Chromatin, where validation selected the strongest regularization value, $\lambda^\star=10$.

The two largest datasets show almost no difference between graph and hypergraph methods:

- Coronavirus: graph $F_1=0.242$ and penalized-hypergraph $F_1=0.241$.
- Cancer: graph $F_1=0.049$ and penalized-hypergraph $F_1=0.049$ after rounding.

For these large datasets, a fixed value of $k=15$ compresses thousands of proteins into a very small number of clusters. This severe under-segmentation prevents either spectral representation from recovering the much finer gold-standard complexes. Consequently, the graph-versus-hypergraph choice has little effect on absolute performance.

Across the nine datasets, the hypergraph improvement has a descriptive negative association with the number of training proteins and a descriptive positive association with the selected penalty:

| Association | Spearman $\rho$ | $p$-value |
|---|---:|---:|
| $\Delta F_1$ versus training protein count | -0.37 | 0.332 |
| $\Delta F_1$ versus $\lambda^\star$ | +0.35 | 0.354 |

Neither relationship is statistically significant because only nine datasets are available. These correlations should therefore be interpreted as exploratory patterns rather than confirmed general relationships.

---

## 6. Adaptive-$k$ comparison

In the adaptive-$k$ experiment, each spectral method independently selected its number of clusters using the eigengap heuristic. Under this setting, graph spectral clustering frequently outperformed the hypergraph method.

| Dataset | Adaptive graph $F_1$ | Graph $k$ | Adaptive hypergraph $F_1$ | Hypergraph $k$ |
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

Cancer and Chromatin selected different values of $k$ across some seeds; therefore, their reported $k$ values are averages. The largest adaptive-$k$ failures occur when the hypergraph eigengap selects very few clusters. For example, it selects $k=2$ for Chromatin, $k=6$ for Crohn's disease, and $k=3$ for Diabetes. These values are too small to represent the underlying complex structure effectively.

Cyanobacteria is the major exception. The graph method selects only two clusters and obtains an $F_1$ score of $0.160$, whereas the hypergraph selects $k=14$ and achieves $0.492$.

These results demonstrate that the hypergraph advantage is specific to the **controlled matched-$k$ setting**. They do not support a general claim that the hypergraph method remains superior when cluster granularity is selected independently.

---

## 7. Randomized-control analysis

The penalized hypergraph was compared with two null models:

1. a degree-preserving interaction-rewiring control; and
2. a hyperedge-size-preserving control.

The null results are averaged over 10 randomizations per control type.

| Dataset | Real penalized hypergraph | Degree-preserving null | Hyperedge-size-preserving null |
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

Real-data performance exceeds both controls on eight of the nine datasets. The contrast is particularly strong for Chromatin, Cyanobacteria, Coronavirus, and Diabetes. For example, Chromatin achieves $F_1=0.244$ on the real data, compared with $0.080$ and $0.088$ under the two null models.

Cancer is the only exception. Its real-data score of approximately $0.049$ is equal to the degree-preserving null score and only slightly above the size-preserving null score. This result supports the interpretation that Cancer is at a performance floor under $k=15$. The apparent graph-hypergraph tie is therefore a tie between two poorly resolved solutions rather than between two high-quality partitions.

---

## 8. Comparison with Louvain and Markov clustering

Louvain and Markov clustering (MCL) determine their own numbers of communities. Their results provide useful reference points for absolute complex recovery, but they are not directly controlled comparisons with the matched-$k$ spectral methods.

Louvain results are averaged over five runs. One MCL result is available for each dataset.

| Dataset | Penalized hypergraph, matched $k$ | Louvain $F_1$ | Louvain communities | MCL $F_1$ | MCL communities |
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

Louvain or MCL obtains a higher absolute best-match $F_1$ score than the matched-$k$ spectral methods on most datasets. This is expected because best-match $F_1$ favors partitions whose granularity resembles the often-small gold-standard complexes. The spectral methods produce only 14 or 15 clusters, whereas MCL produces hundreds or thousands of communities on several datasets.

These results show that cluster granularity has a larger effect on absolute complex recovery than the choice between the graph and hypergraph spectral operators. The main contribution of the spectral analysis is therefore the controlled comparison between graph and hypergraph representations, not a claim that fixed-$k$ spectral clustering is the best overall complex-detection approach.

---

## 9. GO enrichment of recovered modules

Both graph and penalized-hypergraph partitions produce strongly GO-enriched modules. Enrichment was assessed using a hypergeometric test followed by Benjamini-Hochberg correction, with significance defined as $q\leq0.05$.

| Dataset | Graph: enriched clusters | Hypergraph: enriched clusters | Graph: enriched proteins | Hypergraph: enriched proteins |
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

Six datasets have at least 90% enriched clusters for both methods. Crohn's disease, Cyanobacteria, and Diabetes show lower or mixed cluster-level enrichment, although protein-level enrichment remains high for Cyanobacteria and Diabetes.

GO enrichment therefore does not clearly distinguish the graph and hypergraph operators. Both methods generally recover biologically coherent modules. The operator choice appears to affect **which proteins are grouped together** more than it affects whether the resulting modules are functionally enriched.

For Chromatin, the enriched terms recovered from penalized-hypergraph modules include nucleoplasm, chromatin-remodeling, and ribosomal categories. These terms are consistent with the expected biological composition of chromatin-associated protein assemblies. Because the available results contain GO identifiers but not the complete GO directed acyclic graph, ontology-based Lin or Resnik semantic-similarity analysis requires an additional `go-basic.obo` file.

---

## 10. Overall interpretation

The combined results support the following conclusions:

1. Degree-aware hypergraph regularization provides significant improvements on four of the nine datasets.
2. The improvement is dataset-dependent and is strongest for smaller, group-structured datasets and for Chromatin, where a strong regularization value was selected.
3. The unregularized hypergraph operator does not consistently outperform graph spectral clustering.
4. The matched-$k$ hypergraph advantage does not reliably persist when each method selects $k$ independently.
5. Real interaction structure produces substantially better recovery than randomized controls on all datasets except Cancer.
6. Graph and hypergraph partitions are both strongly GO-enriched on most datasets.
7. Community granularity has a larger influence on absolute recovery than the choice of spectral operator.

Accordingly, the evidence supports the conclusion that **degree-aware hypergraph regularization can improve protein-complex recovery under particular structural and granularity conditions**. It does not support the stronger claim that hypergraph spectral clustering is universally superior to graph spectral clustering.

---

## 11. Run completeness

All nine dataset analyses completed successfully. The recorded end-to-end runtime was approximately **2,821 seconds**, or **47.0 minutes**.

The numerical results summarized here were obtained from:

- `ALL_manifest.csv`
- `ALL_test_metrics.csv`
- `ALL_controls.csv`
- `ALL_go_enrichment_summary.csv`
- `dataset_progress.csv`
