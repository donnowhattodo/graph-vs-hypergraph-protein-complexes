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

## 9. Complete evaluation metric matrices

The preceding sections emphasize symmetric best-match $F_1$ because it is the primary held-out recovery measure. However, `ALL_test_metrics.csv` also contains pairwise, best-match, B$^3$, information-theoretic, hyperedge-recovery, and modularity metrics. This section reports all of those values.

For methods evaluated over multiple seeds or runs, each entry is shown as **mean (standard deviation)**. Matched-$k$ and adaptive-$k$ spectral methods use 25 seeds per dataset; Louvain uses five runs per dataset; MCL has one run per dataset and therefore has no standard deviation. The $k$ column is also averaged when the selected number of clusters varies across runs. `NA` indicates that the metric was not defined in the source output. In particular, NMI and ARI are unavailable for adaptive graph clustering on Cyanobacteria.

### 9.1 Metric interpretation

- **Pairwise precision, recall, and $F_1$** evaluate whether pairs of proteins are placed together consistently with the gold complexes.
- **Best-match P$\to$G** averages the best gold-complex match for each predicted cluster, whereas **Best-match G$\to$P** averages the best predicted-cluster match for each gold complex. Their harmonic mean is the **symmetric best-match $F_1$**.
- **B$^3$ precision, recall, and $F_1$** evaluate clustering quality at the individual-protein level.
- **NMI** and **ARI** measure agreement between the predicted and reference partitions, with adjustment for chance in ARI.
- **Hyperedge recovery $F_1$** and **hyperedge recovery fraction** quantify direct recovery of held-out complexes.
- **Modularity** measures the strength of within-module connectivity under the corresponding representation.

### 9.2 Matched-$k$ spectral metrics

All three spectral operators use the same validation-selected $k$ for each dataset, making these rows the controlled graph-versus-hypergraph comparison.

| `dataset` | `method` | `k` | `lam` | `pairwise_precision` | `pairwise_recall` | `pairwise_f1` | `bestmatch_pred_gold` | `bestmatch_gold_pred` | `bestmatch_f1_sym` | `b3_precision` | `b3_recall` | `b3_f1` | `NMI` | `ARI` | `hyperedge_recovery_f1` | `hyperedge_recovery_frac` | `modularity` |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Affinomics | Graph spectral | 15 | 0 | 0.0199 (0.0006) | 0.8381 (0.0163) | 0.0390 (0.0012) | 0.1891 (0.0057) | 0.1739 (0.0093) | 0.1815 (0.0060) | 0.0210 (0.0012) | 0.2549 (0.0034) | 0.0388 (0.0020) | 0.7329 (0.0074) | 0.1727 (0.0058) | 0.1662 (0.0084) | 0.1771 (0.0575) | 0.7990 (0.0084) |
| Affinomics | Hypergraph ($\lambda=0$) | 15 | 0 | 0.0209 (0.0022) | 0.9434 (0.0000) | 0.0408 (0.0042) | 0.2047 (0.0087) | 0.2017 (0.0074) | 0.2032 (0.0068) | 0.0288 (0.0011) | 0.2797 (0.0000) | 0.0522 (0.0018) | 0.7643 (0.0093) | 0.2067 (0.0097) | 0.1891 (0.0089) | 0.2337 (0.0351) | 0.8964 (0.0052) |
| Affinomics | Penalized hypergraph | 15 | 0 | 0.0209 (0.0022) | 0.9434 (0.0000) | 0.0408 (0.0042) | 0.2047 (0.0087) | 0.2017 (0.0074) | 0.2032 (0.0068) | 0.0288 (0.0011) | 0.2797 (0.0000) | 0.0522 (0.0018) | 0.7643 (0.0093) | 0.2067 (0.0097) | 0.1891 (0.0089) | 0.2337 (0.0351) | 0.8964 (0.0052) |
| BioCreative | Graph spectral | 15 | 0 | 0.0202 (0.0006) | 0.9462 (0.0060) | 0.0395 (0.0012) | 0.2035 (0.0065) | 0.1806 (0.0074) | 0.1921 (0.0052) | 0.0227 (0.0012) | 0.3203 (0.0011) | 0.0425 (0.0020) | 0.7652 (0.0091) | 0.1941 (0.0186) | 0.1732 (0.0067) | 0.1463 (0.0520) | 0.8927 (0.0041) |
| BioCreative | Hypergraph ($\lambda=0$) | 15 | 0 | 0.0204 (0.0004) | 0.9383 (0.0000) | 0.0399 (0.0009) | 0.1951 (0.0049) | 0.1735 (0.0054) | 0.1843 (0.0046) | 0.0212 (0.0008) | 0.3188 (0.0000) | 0.0398 (0.0014) | 0.7658 (0.0064) | 0.2076 (0.0091) | 0.1676 (0.0050) | 0.0767 (0.0446) | 0.9233 (0.0033) |
| BioCreative | Penalized hypergraph | 15 | 0.0010 (0.0000) | 0.0172 (0.0010) | 0.9383 (0.0000) | 0.0337 (0.0019) | 0.2025 (0.0061) | 0.1854 (0.0053) | 0.1940 (0.0050) | 0.0252 (0.0010) | 0.3188 (0.0000) | 0.0467 (0.0017) | 0.7691 (0.0081) | 0.2014 (0.0121) | 0.1743 (0.0056) | 0.1556 (0.0331) | 0.9038 (0.0061) |
| Cancer | Graph spectral | 15 | 0 | 0.0055 (0.0006) | 0.6928 (0.0300) | 0.0109 (0.0012) | 0.0736 (0.0009) | 0.0240 (0.0012) | 0.0488 (0.0006) | 0.0055 (0.0002) | 0.2409 (0.0114) | 0.0107 (0.0003) | 0.3719 (0.0220) | 0.0175 (0.0040) | 0.0148 (0.0007) | 0.0008 (0.0000) | 0.4765 (0.0204) |
| Cancer | Hypergraph ($\lambda=0$) | 15 | 0 | 0.0053 (0.0005) | 0.7730 (0.0405) | 0.0104 (0.0009) | 0.0735 (0.0019) | 0.0221 (0.0014) | 0.0478 (0.0016) | 0.0049 (0.0002) | 0.2572 (0.0147) | 0.0096 (0.0005) | 0.3350 (0.0301) | 0.0144 (0.0035) | 0.0140 (0.0009) | 0.0001 (0.0003) | 0.4518 (0.0482) |
| Cancer | Penalized hypergraph | 15 | 0 | 0.0054 (0.0003) | 0.7558 (0.0280) | 0.0108 (0.0006) | 0.0743 (0.0015) | 0.0228 (0.0012) | 0.0485 (0.0013) | 0.0050 (0.0002) | 0.2518 (0.0102) | 0.0098 (0.0003) | 0.3470 (0.0186) | 0.0157 (0.0025) | 0.0144 (0.0009) | 0.0001 (0.0003) | 0.4699 (0.0123) |
| Cardiac | Graph spectral | 15 | 0 | 0.0085 (0.0002) | 0.8225 (0.0038) | 0.0169 (0.0004) | 0.1249 (0.0019) | 0.1049 (0.0049) | 0.1149 (0.0026) | 0.0155 (0.0013) | 0.2284 (0.0022) | 0.0290 (0.0022) | 0.6030 (0.0054) | 0.0625 (0.0019) | 0.0752 (0.0026) | 0.0556 (0.0155) | 0.8202 (0.0020) |
| Cardiac | Hypergraph ($\lambda=0$) | 15 | 0 | 0.0085 (0.0002) | 0.8017 (0.0026) | 0.0168 (0.0004) | 0.1237 (0.0014) | 0.1018 (0.0074) | 0.1128 (0.0036) | 0.0150 (0.0022) | 0.2229 (0.0009) | 0.0282 (0.0038) | 0.6046 (0.0060) | 0.0599 (0.0022) | 0.0742 (0.0030) | 0.0452 (0.0154) | 0.8236 (0.0025) |
| Cardiac | Penalized hypergraph | 15 | 0 | 0.0085 (0.0001) | 0.8026 (0.0025) | 0.0167 (0.0001) | 0.1233 (0.0015) | 0.1014 (0.0058) | 0.1124 (0.0028) | 0.0148 (0.0017) | 0.2232 (0.0019) | 0.0278 (0.0029) | 0.6041 (0.0057) | 0.0598 (0.0024) | 0.0741 (0.0027) | 0.0452 (0.0156) | 0.8234 (0.0025) |
| Chromatin | Graph spectral | 15 | 0 | 0.0732 (0.0018) | 0.9800 (0.0032) | 0.1363 (0.0031) | 0.2502 (0.0038) | 0.1633 (0.0061) | 0.2068 (0.0021) | 0.0579 (0.0024) | 0.4185 (0.0022) | 0.1017 (0.0037) | 0.5649 (0.0078) | 0.1382 (0.0016) | 0.0917 (0.0035) | 0.0652 (0.0061) | 0.6514 (0.0026) |
| Chromatin | Hypergraph ($\lambda=0$) | 15 | 0 | 0.0749 (0.0016) | 0.9781 (0.0019) | 0.1391 (0.0027) | 0.2518 (0.0033) | 0.1667 (0.0012) | 0.2093 (0.0021) | 0.0593 (0.0005) | 0.4145 (0.0020) | 0.1037 (0.0007) | 0.5686 (0.0045) | 0.1383 (0.0014) | 0.0932 (0.0023) | 0.0653 (0.0074) | 0.6502 (0.0024) |
| Chromatin | Penalized hypergraph | 15 | 10 | 0.0936 (0.0025) | 0.7028 (0.0174) | 0.1651 (0.0039) | 0.2762 (0.0031) | 0.2113 (0.0056) | 0.2438 (0.0042) | 0.0859 (0.0045) | 0.3434 (0.0039) | 0.1374 (0.0059) | 0.6148 (0.0053) | 0.2423 (0.0103) | 0.0856 (0.0009) | 0.0802 (0.0053) | 0.6445 (0.0028) |
| Coronavirus | Graph spectral | 15 | 0 | 0.0863 (0.0039) | 0.7376 (0.0337) | 0.1545 (0.0062) | 0.2658 (0.0067) | 0.2179 (0.0086) | 0.2418 (0.0060) | 0.0671 (0.0035) | 0.2226 (0.0067) | 0.1031 (0.0043) | 0.5374 (0.0091) | 0.4772 (0.0209) | 0.0320 (0.0010) | 0.0095 (0.0020) | 0.6182 (0.0050) |
| Coronavirus | Hypergraph ($\lambda=0$) | 15 | 0 | 0.0749 (0.0019) | 0.7598 (0.0076) | 0.1364 (0.0031) | 0.2730 (0.0043) | 0.2082 (0.0051) | 0.2406 (0.0038) | 0.0620 (0.0021) | 0.2247 (0.0029) | 0.0972 (0.0027) | 0.5357 (0.0074) | 0.4596 (0.0118) | 0.0377 (0.0014) | 0.0110 (0.0029) | 0.6073 (0.0060) |
| Coronavirus | Penalized hypergraph | 15 | 0 | 0.0752 (0.0023) | 0.7610 (0.0102) | 0.1369 (0.0037) | 0.2724 (0.0038) | 0.2087 (0.0053) | 0.2406 (0.0034) | 0.0625 (0.0023) | 0.2247 (0.0032) | 0.0978 (0.0029) | 0.5359 (0.0063) | 0.4611 (0.0119) | 0.0371 (0.0015) | 0.0107 (0.0031) | 0.6061 (0.0056) |
| Crohn's disease | Graph spectral | 14 | 0 | 0.2325 (0.0040) | 0.9175 (0.0042) | 0.3710 (0.0052) | 0.5052 (0.0151) | 0.5587 (0.0158) | 0.5320 (0.0150) | 0.2328 (0.0131) | 0.4737 (0.0035) | 0.3120 (0.0121) | 0.8383 (0.0128) | 0.4648 (0.0246) | 0.5249 (0.0185) | 1.0000 (0.0000) | 0.8065 (0.0014) |
| Crohn's disease | Hypergraph ($\lambda=0$) | 14 | 0 | 0.2350 (0.0081) | 0.9167 (0.0000) | 0.3740 (0.0104) | 0.5000 (0.0082) | 0.5527 (0.0133) | 0.5264 (0.0108) | 0.2240 (0.0082) | 0.4730 (0.0000) | 0.3039 (0.0076) | 0.8296 (0.0101) | 0.4488 (0.0354) | 0.5194 (0.0123) | 1.0000 (0.0000) | 0.8078 (0.0018) |
| Crohn's disease | Penalized hypergraph | 14 | 1 | 0.2178 (0.0078) | 0.7875 (0.0361) | 0.3412 (0.0129) | 0.5094 (0.0038) | 0.5778 (0.0086) | 0.5436 (0.0046) | 0.2376 (0.0069) | 0.4050 (0.0172) | 0.2995 (0.0096) | 0.8571 (0.0102) | 0.5446 (0.0208) | 0.5467 (0.0101) | 1.0000 (0.0000) | 0.7552 (0.0146) |
| Cyanobacteria | Graph spectral | 15 | 0 | 0.1170 (0.0006) | 0.9699 (0.0000) | 0.2089 (0.0009) | 0.4511 (0.0011) | 0.4913 (0.0094) | 0.4712 (0.0052) | 0.1407 (0.0034) | 0.3667 (0.0000) | 0.2033 (0.0035) | 0.8885 (0.0034) | 0.8617 (0.0017) | 0.3946 (0.0178) | 0.7500 (0.0000) | 0.8637 (0.0004) |
| Cyanobacteria | Hypergraph ($\lambda=0$) | 15 | 0 | 0.1169 (0.0000) | 0.9676 (0.0000) | 0.2086 (0.0001) | 0.4522 (0.0003) | 0.4937 (0.0017) | 0.4729 (0.0007) | 0.1485 (0.0000) | 0.3566 (0.0000) | 0.2097 (0.0000) | 0.8925 (0.0000) | 0.8665 (0.0000) | 0.3949 (0.0035) | 0.7500 (0.0000) | 0.8533 (0.0000) |
| Cyanobacteria | Penalized hypergraph | 15 | 0.0010 (0.0000) | 0.1039 (0.0000) | 0.9699 (0.0000) | 0.1877 (0.0000) | 0.4617 (0.0000) | 0.5258 (0.0000) | 0.4938 (0.0000) | 0.1580 (0.0000) | 0.3667 (0.0000) | 0.2209 (0.0000) | 0.9052 (0.0000) | 0.8702 (0.0000) | 0.4606 (0.0000) | 0.7500 (0.0000) | 0.8252 (0.0000) |
| Diabetes | Graph spectral | 15 | 0 | 0.0005 (0.0000) | 0.9434 (0.0079) | 0.0010 (0.0000) | 0.0605 (0.0027) | 0.4208 (0.0239) | 0.2406 (0.0122) | 0.0206 (0.0019) | 0.0615 (0.0008) | 0.0308 (0.0022) | 0.7887 (0.0125) | 0.2689 (0.0125) | 0.3964 (0.0269) | 0.7394 (0.0537) | 0.5753 (0.0193) |
| Diabetes | Hypergraph ($\lambda=0$) | 15 | 0 | 0.0004 (0.0000) | 0.8131 (0.0450) | 0.0008 (0.0001) | 0.0600 (0.0035) | 0.4161 (0.0221) | 0.2381 (0.0115) | 0.0195 (0.0021) | 0.0538 (0.0021) | 0.0285 (0.0023) | 0.7973 (0.0183) | 0.2936 (0.0226) | 0.3858 (0.0246) | 0.8549 (0.0723) | 0.5650 (0.0108) |
| Diabetes | Penalized hypergraph | 15 | 1 | 0.0004 (0.0000) | 0.7055 (0.0048) | 0.0007 (0.0000) | 0.0653 (0.0022) | 0.4025 (0.0119) | 0.2339 (0.0065) | 0.0158 (0.0010) | 0.0502 (0.0008) | 0.0240 (0.0011) | 0.8429 (0.0114) | 0.3966 (0.0279) | 0.3807 (0.0123) | 0.7646 (0.0322) | 0.5147 (0.0021) |

### 9.3 Adaptive-$k$ spectral metrics

In this experiment, graph and penalized-hypergraph methods independently select $k$ using the eigengap heuristic.

| `dataset` | `method` | `k` | `lam` | `pairwise_precision` | `pairwise_recall` | `pairwise_f1` | `bestmatch_pred_gold` | `bestmatch_gold_pred` | `bestmatch_f1_sym` | `b3_precision` | `b3_recall` | `b3_f1` | `NMI` | `ARI` | `hyperedge_recovery_f1` | `hyperedge_recovery_frac` | `modularity` |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Affinomics | Graph spectral | 6 | 0 | 0.0095 (0.0003) | 0.8600 (0.0148) | 0.0187 (0.0006) | 0.1141 (0.0003) | 0.0840 (0.0022) | 0.0990 (0.0011) | 0.0102 (0.0003) | 0.2615 (0.0041) | 0.0196 (0.0006) | 0.5847 (0.0055) | 0.0900 (0.0046) | 0.0752 (0.0021) | 0.0032 (0.0065) | 0.7077 (0.0182) |
| Affinomics | Penalized hypergraph | 4 | 0 | 0.0066 (0.0003) | 0.9921 (0.0035) | 0.0132 (0.0007) | 0.0908 (0.0040) | 0.0619 (0.0022) | 0.0763 (0.0029) | 0.0077 (0.0004) | 0.2863 (0.0023) | 0.0151 (0.0007) | 0.4835 (0.0151) | 0.0545 (0.0056) | 0.0561 (0.0017) | 0.0000 (0.0000) | 0.7182 (0.0202) |
| BioCreative | Graph spectral | 9 | 0 | 0.0115 (0.0007) | 0.9383 (0.0000) | 0.0227 (0.0013) | 0.1457 (0.0028) | 0.1130 (0.0038) | 0.1294 (0.0026) | 0.0137 (0.0006) | 0.3188 (0.0000) | 0.0263 (0.0011) | 0.6758 (0.0061) | 0.1240 (0.0074) | 0.1074 (0.0034) | 0.0296 (0.0204) | 0.8440 (0.0094) |
| BioCreative | Penalized hypergraph | 6 | 0.0010 (0.0000) | 0.0030 (0.0001) | 0.9383 (0.0000) | 0.0060 (0.0002) | 0.1061 (0.0005) | 0.0963 (0.0022) | 0.1012 (0.0009) | 0.0127 (0.0003) | 0.3188 (0.0000) | 0.0244 (0.0006) | 0.4428 (0.0078) | 0.0223 (0.0013) | 0.0868 (0.0021) | 0.0970 (0.0176) | 0.6115 (0.0145) |
| Cancer | Graph spectral | 9.12 (0.60) | 0 | 0.0039 (0.0005) | 0.8944 (0.0572) | 0.0078 (0.0010) | 0.0641 (0.0012) | 0.0146 (0.0010) | 0.0393 (0.0010) | 0.0038 (0.0002) | 0.2982 (0.0211) | 0.0076 (0.0004) | 0.2191 (0.0318) | 0.0066 (0.0023) | 0.0077 (0.0006) | 0.0004 (0.0010) | 0.2617 (0.0711) |
| Cancer | Penalized hypergraph | 9.12 (0.33) | 0 | 0.0041 (0.0006) | 0.8808 (0.0617) | 0.0081 (0.0012) | 0.0659 (0.0011) | 0.0153 (0.0010) | 0.0406 (0.0010) | 0.0040 (0.0002) | 0.2920 (0.0232) | 0.0078 (0.0005) | 0.2315 (0.0385) | 0.0081 (0.0029) | 0.0082 (0.0004) | 0.0002 (0.0004) | 0.3004 (0.0941) |
| Cardiac | Graph spectral | 11 | 0 | 0.0063 (0.0001) | 0.8614 (0.0021) | 0.0125 (0.0001) | 0.1087 (0.0012) | 0.0830 (0.0052) | 0.0959 (0.0027) | 0.0126 (0.0010) | 0.2350 (0.0016) | 0.0239 (0.0019) | 0.5518 (0.0072) | 0.0478 (0.0019) | 0.0543 (0.0022) | 0.0165 (0.0040) | 0.7881 (0.0032) |
| Cardiac | Penalized hypergraph | 11 | 0 | 0.0053 (0.0009) | 0.8821 (0.0373) | 0.0105 (0.0018) | 0.1068 (0.0058) | 0.0838 (0.0048) | 0.0953 (0.0032) | 0.0126 (0.0009) | 0.2403 (0.0055) | 0.0240 (0.0017) | 0.5551 (0.0112) | 0.0459 (0.0044) | 0.0566 (0.0026) | 0.0248 (0.0107) | 0.7713 (0.0247) |
| Chromatin | Graph spectral | 12.68 (0.90) | 0 | 0.0683 (0.0051) | 0.9858 (0.0019) | 0.1278 (0.0090) | 0.2438 (0.0072) | 0.1556 (0.0099) | 0.1997 (0.0083) | 0.0553 (0.0035) | 0.4244 (0.0044) | 0.0978 (0.0055) | 0.5352 (0.0294) | 0.1332 (0.0084) | 0.0779 (0.0089) | 0.0591 (0.0070) | 0.6445 (0.0097) |
| Chromatin | Penalized hypergraph | 2 | 10 | 0.0198 (0.0000) | 0.9999 (0.0000) | 0.0388 (0.0000) | 0.1524 (0.0000) | 0.0553 (0.0000) | 0.1038 (0.0000) | 0.0194 (0.0000) | 0.4566 (0.0000) | 0.0372 (0.0000) | 0.0702 (0.0000) | 0.0049 (0.0000) | 0.0127 (0.0000) | 0.0000 (0.0000) | 0.0903 (0.0000) |
| Coronavirus | Graph spectral | 8 | 0 | 0.0543 (0.0046) | 0.7690 (0.0225) | 0.1014 (0.0080) | 0.2479 (0.0045) | 0.1673 (0.0096) | 0.2076 (0.0043) | 0.0492 (0.0031) | 0.2374 (0.0114) | 0.0814 (0.0038) | 0.4245 (0.0220) | 0.3319 (0.0295) | 0.0176 (0.0015) | 0.0096 (0.0029) | 0.5527 (0.0054) |
| Coronavirus | Penalized hypergraph | 8 | 0 | 0.0327 (0.0023) | 0.6829 (0.0625) | 0.0623 (0.0040) | 0.2074 (0.0253) | 0.1264 (0.0243) | 0.1669 (0.0247) | 0.0332 (0.0068) | 0.2149 (0.0173) | 0.0570 (0.0099) | 0.3325 (0.0268) | 0.1834 (0.0306) | 0.0220 (0.0015) | 0.0262 (0.0029) | 0.3958 (0.0126) |
| Crohn's disease | Graph spectral | 14 | 0 | 0.2398 (0.0034) | 0.9167 (0.0000) | 0.3801 (0.0044) | 0.5054 (0.0124) | 0.5552 (0.0077) | 0.5303 (0.0096) | 0.2228 (0.0091) | 0.4730 (0.0000) | 0.3028 (0.0085) | 0.8381 (0.0089) | 0.4702 (0.0175) | 0.5201 (0.0081) | 1.0000 (0.0000) | 0.8084 (0.0019) |
| Crohn's disease | Penalized hypergraph | 6 | 1 | 0.0544 (0.0000) | 0.9167 (0.0000) | 0.1027 (0.0000) | 0.3311 (0.0000) | 0.2619 (0.0000) | 0.2965 (0.0000) | 0.0740 (0.0000) | 0.4730 (0.0000) | 0.1279 (0.0000) | 0.5326 (0.0000) | 0.0918 (0.0000) | 0.2598 (0.0000) | 0.4400 (0.0000) | 0.6343 (0.0000) |
| Cyanobacteria | Graph spectral | 2 | 0 | 0.0238 (0.0000) | 1.0000 (0.0000) | 0.0465 (0.0000) | 0.2233 (0.0000) | 0.0962 (0.0000) | 0.1598 (0.0000) | 0.0230 (0.0000) | 0.3939 (0.0000) | 0.0434 (0.0000) | NA | NA | 0.0407 (0.0000) | 0.0000 (0.0000) | 0.1360 (0.0000) |
| Cyanobacteria | Penalized hypergraph | 14 | 0.0010 (0.0000) | 0.1038 (0.0000) | 0.9699 (0.0000) | 0.1875 (0.0000) | 0.4641 (0.0015) | 0.5193 (0.0050) | 0.4917 (0.0019) | 0.1532 (0.0031) | 0.3667 (0.0000) | 0.2161 (0.0031) | 0.9052 (0.0000) | 0.8702 (0.0000) | 0.4472 (0.0102) | 0.7500 (0.0000) | 0.8251 (0.0000) |
| Diabetes | Graph spectral | 15 | 0 | 0.0005 (0.0000) | 0.9441 (0.0075) | 0.0010 (0.0000) | 0.0600 (0.0040) | 0.4205 (0.0270) | 0.2403 (0.0148) | 0.0206 (0.0019) | 0.0616 (0.0008) | 0.0308 (0.0021) | 0.7842 (0.0233) | 0.2639 (0.0234) | 0.3958 (0.0317) | 0.7269 (0.0584) | 0.5777 (0.0198) |
| Diabetes | Penalized hypergraph | 3 | 1 | 0.0003 (0.0000) | 0.9138 (0.0000) | 0.0005 (0.0000) | 0.0294 (0.0000) | 0.0985 (0.0000) | 0.0640 (0.0000) | 0.0028 (0.0000) | 0.0577 (0.0000) | 0.0053 (0.0000) | 0.3330 (0.0000) | 0.0236 (0.0000) | 0.1025 (0.0000) | 0.2286 (0.0000) | 0.1349 (0.0000) |

### 9.4 Native-granularity baseline metrics

Louvain and MCL determine their own numbers of communities. Their values are therefore useful for absolute-performance comparison but are not matched-granularity comparisons with the spectral methods.

| `dataset` | `method` | `k` | `lam` | `pairwise_precision` | `pairwise_recall` | `pairwise_f1` | `bestmatch_pred_gold` | `bestmatch_gold_pred` | `bestmatch_f1_sym` | `b3_precision` | `b3_recall` | `b3_f1` | `NMI` | `ARI` | `hyperedge_recovery_f1` | `hyperedge_recovery_frac` | `modularity` |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Affinomics | Louvain | 54.20 (0.45) | 0 | 0.0439 (0.0012) | 0.9075 (0.0103) | 0.0837 (0.0021) | 0.2603 (0.0030) | 0.3802 (0.0101) | 0.3203 (0.0064) | 0.0629 (0.0017) | 0.2676 (0.0024) | 0.1018 (0.0022) | 0.8439 (0.0062) | 0.2951 (0.0225) | 0.3893 (0.0094) | 0.6635 (0.0519) | 0.9334 (0.0021) |
| Affinomics | MCL | 116 | 0 | 0.0856 | 0.5943 | 0.1496 | 0.2882 | 0.5326 | 0.4104 | 0.0973 | 0.2205 | 0.1350 | 0.8964 | 0.4078 | 0.5441 | 0.9524 | 0.8195 |
| BioCreative | Louvain | 94 | 0 | 0.1186 (0.0000) | 0.9383 (0.0000) | 0.2105 (0.0000) | 0.3333 (0.0000) | 0.6343 (0.0000) | 0.4838 (0.0000) | 0.1578 (0.0000) | 0.3188 (0.0000) | 0.2111 (0.0000) | 0.9328 (0.0000) | 0.5581 (0.0000) | 0.6279 (0.0000) | 0.9041 (0.0000) | 0.9737 (0.0000) |
| BioCreative | MCL | 107 | 0 | 0.1447 | 0.8148 | 0.2458 | 0.3543 | 0.6926 | 0.5234 | 0.1698 | 0.2793 | 0.2112 | 0.9570 | 0.7000 | 0.6863 | 1.0000 | 0.9164 |
| Cancer | Louvain | 137.60 (2.07) | 0 | 0.0090 (0.0003) | 0.4716 (0.0102) | 0.0177 (0.0006) | 0.0997 (0.0021) | 0.0661 (0.0026) | 0.0829 (0.0023) | 0.0175 (0.0005) | 0.2075 (0.0022) | 0.0322 (0.0008) | 0.5504 (0.0033) | 0.0477 (0.0006) | 0.0694 (0.0027) | 0.0574 (0.0040) | 0.6548 (0.0007) |
| Cancer | MCL | 3426 | 0 | 0.0538 | 0.0564 | 0.0551 | 0.1676 | 0.3072 | 0.2374 | 0.0600 | 0.0984 | 0.0746 | 0.8047 | 0.0475 | 0.3917 | 0.6858 | 0.3335 |
| Cardiac | Louvain | 63.00 (0.71) | 0 | 0.0121 (0.0006) | 0.7724 (0.0122) | 0.0238 (0.0012) | 0.1730 (0.0026) | 0.2578 (0.0048) | 0.2154 (0.0034) | 0.0556 (0.0007) | 0.2198 (0.0034) | 0.0887 (0.0008) | 0.7384 (0.0077) | 0.1162 (0.0100) | 0.2340 (0.0059) | 0.2909 (0.0145) | 0.8545 (0.0004) |
| Cardiac | MCL | 574 | 0 | 0.0407 | 0.4535 | 0.0747 | 0.1963 | 0.4415 | 0.3189 | 0.0792 | 0.1533 | 0.1044 | 0.8424 | 0.1731 | 0.4407 | 0.7273 | 0.5789 |
| Chromatin | Louvain | 46.80 (1.10) | 0 | 0.1366 (0.0044) | 0.8191 (0.0162) | 0.2340 (0.0063) | 0.3188 (0.0241) | 0.2748 (0.0137) | 0.2968 (0.0187) | 0.1056 (0.0057) | 0.3753 (0.0030) | 0.1647 (0.0069) | 0.6777 (0.0199) | 0.2744 (0.0384) | 0.1815 (0.0096) | 0.2255 (0.0140) | 0.7377 (0.0015) |
| Chromatin | MCL | 742 | 0 | 0.4118 | 0.4119 | 0.4119 | 0.2928 | 0.4499 | 0.3714 | 0.1521 | 0.2271 | 0.1822 | 0.7623 | 0.1995 | 0.4520 | 0.8008 | 0.5350 |
| Coronavirus | Louvain | 32.00 (1.00) | 0 | 0.1037 (0.0046) | 0.6878 (0.0116) | 0.1802 (0.0072) | 0.2736 (0.0024) | 0.2465 (0.0127) | 0.2601 (0.0054) | 0.0756 (0.0059) | 0.2199 (0.0023) | 0.1124 (0.0067) | 0.5818 (0.0067) | 0.4833 (0.0134) | 0.0691 (0.0018) | 0.0619 (0.0052) | 0.7030 (0.0014) |
| Coronavirus | MCL | 2164 | 0 | 0.1939 | 0.0794 | 0.1126 | 0.1158 | 0.2593 | 0.1876 | 0.0457 | 0.0524 | 0.0488 | 0.6699 | 0.1203 | 0.2537 | 0.4260 | 0.2351 |
| Crohn's disease | Louvain | 13 | 0 | 0.2035 (0.0000) | 0.9792 (0.0000) | 0.3369 (0.0000) | 0.5131 (0.0000) | 0.5430 (0.0000) | 0.5281 (0.0000) | 0.2262 (0.0000) | 0.5035 (0.0000) | 0.3122 (0.0000) | 0.8517 (0.0000) | 0.5027 (0.0000) | 0.5062 (0.0000) | 1.0000 (0.0000) | 0.8497 (0.0000) |
| Crohn's disease | MCL | 16 | 0 | 0.1958 | 0.7708 | 0.3122 | 0.5062 | 0.5974 | 0.5518 | 0.2461 | 0.4099 | 0.3076 | 0.8763 | 0.5869 | 0.5756 | 1.0000 | 0.7658 |
| Cyanobacteria | Louvain | 18 | 0 | 0.1428 (0.0003) | 0.9699 (0.0000) | 0.2490 (0.0004) | 0.4708 (0.0003) | 0.5409 (0.0006) | 0.5059 (0.0002) | 0.1639 (0.0002) | 0.3667 (0.0000) | 0.2265 (0.0002) | 0.9245 (0.0000) | 0.8949 (0.0000) | 0.4773 (0.0003) | 0.7500 (0.0000) | 0.8901 (0.0000) |
| Cyanobacteria | MCL | 110 | 0 | 0.2534 | 0.2153 | 0.2328 | 0.2724 | 0.4577 | 0.3650 | 0.1398 | 0.1950 | 0.1629 | 0.7546 | 0.2819 | 0.5172 | 0.9167 | 0.6304 |
| Diabetes | Louvain | 20 | 0 | 0.0005 (0.0000) | 0.9483 (0.0000) | 0.0010 (0.0000) | 0.0578 (0.0000) | 0.4244 (0.0000) | 0.2411 (0.0000) | 0.0247 (0.0000) | 0.0620 (0.0000) | 0.0353 (0.0000) | 0.7837 (0.0000) | 0.2627 (0.0000) | 0.4164 (0.0000) | 0.6857 (0.0000) | 0.6055 (0.0000) |
| Diabetes | MCL | 581 | 0 | 0.1640 | 0.8966 | 0.2773 | 0.0603 | 0.5077 | 0.2840 | 0.0297 | 0.0549 | 0.0386 | 0.8134 | 0.2800 | 0.4883 | 0.8000 | 0.2947 |

---

## 10. GO enrichment of recovered modules

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

## 11. Overall interpretation

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

## 12. Run completeness

All nine dataset analyses completed successfully. The recorded end-to-end runtime was approximately **2,821 seconds**, or **47.0 minutes**.
