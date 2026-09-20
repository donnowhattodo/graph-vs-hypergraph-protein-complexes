# Graph vs Hypergraph Protein-Complex Recovery

[![Preprint](https://img.shields.io/badge/preprint-10.20944%2Fpreprints202608.0397.v1-blue)](https://doi.org/10.20944/preprints202608.0397.v1)
[![License: CC BY 4.0](https://img.shields.io/badge/paper%20license-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

Official code and reproducibility repository for:

> **When Does Higher-Order Representation Improve Protein-Complex Recovery? A Leakage-Controlled, Matched-*k* Comparison of Graph and Penalized Hypergraph Spectral Clustering**  
> Kazi Hafiz Md Asad, Rafi Majid, Md Tanjeelur Rahman Labib, and Ahsanur Rahman  
> *Preprints.org*, posted 6 August 2026  
> DOI: [10.20944/preprints202608.0397.v1](https://doi.org/10.20944/preprints202608.0397.v1)

## Overview

Protein complexes are higher-order molecular assemblies, but interaction data are commonly reduced to pairwise protein-protein interaction graphs before clustering. That projection can discard experiment membership and allow large groups to dominate the representation.

This repository provides a reproducible, leakage-controlled benchmark that isolates the contribution of representation. It compares three spectral operators constructed from the **same training hyperedges**:

1. inverse-size-weighted clique graph spectral clustering;
2. the normalized Zhou hypergraph operator; and
3. a degree-aware penalized hypergraph operator.

The primary question is:

> When the input groups, protein universe, downstream clustering procedure, and number of output clusters are controlled, does retaining higher-order incidence structure improve recovery of held-out protein groups?

## Main findings

Across nine archived IntAct-derived thematic PSI-MITAB collections and 25 paired initializations per dataset:

- the penalized hypergraph achieved a higher mean held-out symmetric best-match F1 than graph spectral clustering on seven of nine datasets;
- six positive differences remained significant after within-dataset Benjamini-Hochberg correction;
- two datasets were statistically indistinguishable;
- Cancer showed a statistically significant but very small decrease of 0.0015 at a low recovery floor;
- positive mean differences ranged from 0.0098 to 0.0356;
- the penalized operator had the best mean rank, but the across-dataset Friedman test was not significant (`chi-square = 5.20`, `p = 0.074`), and the graph-versus-penalized-hypergraph Nemenyi comparison was also non-significant (`p = 0.111`);
- observed recovery exceeded both structural null models on eight of nine datasets; and
- graph and hypergraph partitions had high, substantially overlapping GO coherence, with pooled term Jaccard overlap of 0.78-0.86 and no aspect-level difference surviving multiplicity correction.

The evidence therefore supports a **conditional**, not universal, advantage for higher-order modeling: retain experiment-level group membership when available, select regularization on validation data, control output granularity, and report neutral and negative datasets alongside improvements.

## Matched-*k* results

The table below reports mean held-out symmetric best-match F1 over 25 external seeds. Delta is penalized hypergraph minus graph spectral clustering.

| Dataset | Graph spectral | Hypergraph (`lambda = 0`) | Penalized hypergraph | Delta | BH-adjusted *p* |
|---|---:|---:|---:|---:|---:|
| Affinomics | 0.182 | 0.202 | **0.202** | +0.0208 | <0.001 |
| BioCreative | 0.192 | 0.184 | **0.194** | +0.0019 | 0.210 |
| Cancer | **0.049** | 0.047 | 0.047 | -0.0015 | <0.001 |
| Cardiac | **0.114** | 0.113 | 0.113 | -0.0010 | 0.210 |
| Chromatin | 0.207 | 0.209 | **0.243** | +0.0356 | <0.001 |
| Coronavirus | 0.230 | 0.234 | **0.245** | +0.0152 | <0.001 |
| Crohn's disease | 0.532 | 0.526 | **0.544** | +0.0116 | <0.001 |
| Cyanobacteria | 0.471 | 0.473 | **0.494** | +0.0226 | <0.001 |
| Diabetes | 0.226 | **0.240** | 0.236 | +0.0098 | 0.0029 |

Standard deviations, paired-bootstrap confidence intervals, effect sizes, secondary metrics, adaptive-*k* results, null-model results, and GO analyses are provided in the paper and generated result files.

## Experimental design

### Leakage control

- Hyperedges are stratified by size (`2-3`, `4-6`, and `>6`) and split into nominal 60% training, 20% validation, and 20% test partitions.
- Only training hyperedges define the protein universe and graph/hypergraph operators.
- Validation and test groups are restricted to proteins observed during training.
- Validation data select the penalty and validation-optimal cluster count.
- Test hyperedges remain untouched until final evaluation.

### Controlled operator comparison

- All spectral representations use the same training hyperedges.
- The clique graph uses inverse-size weighting, `1 / (|e| - 1)`, to limit domination by large groups.
- The penalized hypergraph adds a low-degree-sensitive diagonal term to the normalized Zhou operator.
- The validation grid is `lambda in {0, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10}` and `k in {2, ..., 15}`.
- The primary comparison uses a shared cluster count, `k_matched = min(floor(n / 5), 15, n - 1)`.
- Spectral embeddings use identical downstream row normalization and k-means clustering.
- Twenty-five external seeds (`42-66`) are paired across spectral methods, with 10 internal k-means starts.

### Evaluation

The primary endpoint is symmetric best-match F1. Secondary analyses include:

- pairwise precision, recall, and F1;
- overlapping B3;
- normalized mutual information;
- adjusted Rand index;
- mean held-out hyperedge-recovery F1;
- stability-medoid partitions;
- degree-preserving graph nulls;
- hyperedge-size-preserving nulls;
- Louvain, Leiden, and MCL context baselines; and
- Gene Ontology enrichment for BP, CC, and MF with multiple-testing correction and semantic reduction.

Statistical analysis includes paired Wilcoxon signed-rank tests, rank-biserial effects, paired bootstrap confidence intervals, Benjamini-Hochberg correction, an across-dataset Friedman test, and Nemenyi post-hoc comparisons.

## Datasets

The benchmark uses nine archived IntAct-derived thematic PSI-MITAB 2.7 collections:

- Affinomics
- BioCreative
- Cancer
- Cardiac
- Chromatin
- Coronavirus
- Crohn's disease
- Cyanobacteria
- Diabetes

Because IntAct is continuously curated, archived input snapshots, checksums, retrieval metadata, processed manifests, and the frozen Gene Ontology release metadata are essential parts of the reproducible benchmark.

## Repository layout

```text
.
├── hgspectral/                 # Parsing, operators, clustering, metrics, controls, GO, and statistics
├── configs/                    # Pipeline configuration files
├── scripts/                    # Data-generation and workflow helpers
├── tests/                      # Regression tests
├── fin_result/                 # Completed benchmark outputs
├── figures/                    # Generated paper figures
├── run_pipeline.py             # Main benchmark entry point
├── analyze_results.py          # Statistical analysis
├── make_paper_figures.py       # Publication figure generation
├── environment.yml             # Conda environment
└── requirements.txt            # Pip dependencies
```

## Installation

### Conda

```bash
conda env create -f environment.yml
conda activate hgspectral
```

### Pip

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

On Windows, activate the environment with `.venv\\Scripts\\activate`.

## Reproducing the analysis

### 1. Run a smoke test

```bash
PYTHONPATH=. python run_pipeline.py \
  --datasets SynthA \
  --search-roots . \
  --seeds 6 \
  --out smoke_out
```

### 2. Prepare real datasets

Place each PSI-MITAB dataset in the expected layout:

```text
<data-root>/<Dataset>/<Dataset>.txt
```

### 3. Run the full benchmark

```bash
PYTHONPATH=. python run_pipeline.py \
  --config configs/full_mitab.json \
  --search-roots <data-root>
```

### 4. Analyze results and regenerate figures

```bash
PYTHONPATH=. python analyze_results.py \
  --results fin_result/ALL_test_metrics.csv

PYTHONPATH=. python make_paper_figures.py
```

### 5. Run regression tests

```bash
PYTHONPATH=. python tests/test_pipeline.py
```

For Windows PowerShell, set the module path first with `$env:PYTHONPATH = "."`, then run the corresponding Python command.

## Outputs

Depending on the configuration, the workflow produces:

- per-seed held-out test metrics;
- validation-selected hyperparameters;
- graph and hypergraph cluster assignments;
- paired statistical comparisons and bootstrap intervals;
- adaptive-*k* and native graph-baseline results;
- degree-preserving and size-preserving null-model results;
- GO enrichment and semantic-reduction tables;
- dataset and run manifests;
- input checksums and frozen GO-release metadata; and
- publication-ready tables and figures.

## Interpretation notes

- **Matched granularity is central.** Native Louvain, Leiden, and MCL outputs can have many more communities than the matched spectral solutions, so they provide practical context rather than a controlled representation comparison.
- **Penalty selection matters.** Validation selected a non-zero penalty in six datasets; `lambda = 0` remains an explicit ablation.
- **Eigengap selection is not uniformly reliable.** The complete adaptive-*k* hypergraph pipeline selected overly coarse partitions in several collections even when the penalized operator helped under matched *k*.
- **GO enrichment is supportive, not independent ground truth.** Interaction and ontology curation can share provenance.
- **Repeated seeds are algorithmic replicates.** They measure sensitivity to initialization; they do not turn one archived biological collection into 25 independent datasets.

## Citation

If you use this repository, please cite the preprint:

```bibtex
@article{asad2026higherorder,
  title   = {When Does Higher-Order Representation Improve Protein-Complex Recovery? A Leakage-Controlled, Matched-k Comparison of Graph and Penalized Hypergraph Spectral Clustering},
  author  = {Asad, Kazi Hafiz Md and Majid, Rafi and Labib, Md Tanjeelur Rahman and Rahman, Ahsanur},
  journal = {Preprints},
  year    = {2026},
  doi     = {10.20944/preprints202608.0397.v1},
  url     = {https://doi.org/10.20944/preprints202608.0397.v1}
}
```

## Data and code availability

This repository contains the computational workflow, final configuration, processed tables, figures, per-seed metrics, structural controls, run manifests, input checksums, and frozen GO-release metadata associated with the study. IntAct source records remain available through the IntAct/IMEx resources.

## License

The preprint is distributed under the [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).

Before publishing or redistributing the software, add an explicit repository-level software license (for example, MIT, BSD-3-Clause, or Apache-2.0) in a `LICENSE` file. The paper's CC BY 4.0 license does not automatically define the software license.

## Contact

- **Kazi Hafiz Md Asad** - [@donnowhattodo](https://github.com/donnowhattodo)
- **Corresponding author:** Ahsanur Rahman - `ahsanur.rahman@northsouth.edu`

Department of Electrical & Computer Engineering, North South University, Dhaka, Bangladesh.

---

Copyright (c) 2026 Kazi Hafiz Md Asad and contributors.
