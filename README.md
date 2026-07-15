# Graph vs Hypergraph Protein Complexes

This repository implements a reproducible benchmark for comparing graph-based and hypergraph-based spectral clustering on protein-complex recovery from MITAB-derived interaction data. The workflow is designed to be fair, leakage-free, and statistically grounded, with a particular focus on evaluating whether a degree-aware penalized hypergraph Laplacian improves held-out recovery of protein modules compared with clique-expanded graph spectral clustering.

## Project summary

Protein complexes are naturally higher-order interactions involving multiple proteins at once. Standard graph methods convert these relationships into pairwise edges, which can distort the original biological structure. Hypergraph methods preserve the multi-protein relationships explicitly and can recover more meaningful modules.

This project provides a complete analysis pipeline that:

- parses MITAB interaction files,
- constructs hyperedges from protein-complex data,
- builds matched graph and hypergraph operators from the same underlying hyperedge set,
- evaluates methods on held-out test hyperedges,
- compares clustering quality across many seeds and null models,
- performs GO enrichment analysis and statistical testing.

## Why this repository exists

The main scientific goal is to answer a simple but important question:

> Does a penalized hypergraph representation recover protein complexes better than a clique-expanded graph representation when both are compared under the same experimental conditions?

The pipeline is structured to avoid common pitfalls in this comparison, including:

- unequal cluster counts across methods,
- parameter selection using test data,
- leakage between training and evaluation,
- weak or insufficient statistical analysis.

## Main methodological features

The project includes the following design choices:

- representation-controlled comparison: graph and hypergraph are built from the same hyperedge set,
- held-out evaluation: hyperedges are split into train/validation/test partitions,
- matched-k and adaptive-k conditions,
- multi-seed stability analysis,
- randomized controls such as degree-preserving rewiring and size-preserving nulls,
- GO enrichment with multiple-testing correction,
- statistical tests including Wilcoxon, Friedman, and Nemenyi analyses.

## Repository layout

- hgspectral/: core Python package implementing parsing, spectral clustering, metrics, controls, enrichment, and statistics.
- configs/: configuration files for the pipeline.
- scripts/: helper scripts for data generation and workflow setup.
- tests/: regression tests for pipeline behavior.
- fin_result/: output files from completed benchmark runs.
- figures/: generated figures.

## Installation

### Conda

```bash
conda env create -f environment.yml
conda activate hgspectral
```

### Pip

```bash
pip install -r requirements.txt
```

## Quick start

### Smoke test

```bash
PYTHONPATH=. python run_pipeline.py --datasets SynthA --search-roots . --seeds 6 --out smoke_out
```

### Full run on real data

Place each dataset in the expected layout:

```text
<root>/<Dataset>/<Dataset>.txt
```

Then run:

```bash
PYTHONPATH=. python run_pipeline.py --config configs/full_mitab.json --search-roots <data-root>
```

### Analysis and figure generation

```bash
PYTHONPATH=. python analyze_results.py --results fin_result/ALL_test_metrics.csv
PYTHONPATH=. python make_paper_figures.py
```

## Outputs

The pipeline writes per-dataset result files such as:

- test metrics,
- control-model metrics,
- GO enrichment summaries,
- dataset manifests,
- figures and analysis tables.

## Tests

```bash
PYTHONPATH=. python tests/test_pipeline.py
```

## License
This repository is provided as a research codebase. If you intend to publish or redistribute it, please review and adjust the licensing terms according to your institutional or publication requirements and cite properly.

Copyright © Kazi Hafiz Md Asad (GitHub: [@donnowhattodo](https://github.com/donnowhattodo))
