# Graph vs Hypergraph Protein Complexes

## Overview

This repository contains a reproducible computational workflow for comparing graph-based and hypergraph-based spectral clustering for protein-complex recovery from MITAB interaction data. The project evaluates whether a degree-aware penalized hypergraph Laplacian improves held-out recovery of protein modules compared with clique-expanded graph spectral clustering.

## Why this project exists

Protein complexes are often represented as hyperedges connecting multiple proteins. Traditional graph methods convert these higher-order interactions into pairwise edges, which can distort the original structure. Hypergraph methods retain the higher-order relationships directly and can recover more meaningful modules.

This repository provides a benchmark pipeline that:

- parses MITAB interaction files,
- constructs hyperedges from protein complexes,
- builds graph and hypergraph operators from the same source data,
- evaluates methods on held-out test hyperedges,
- compares clustering quality with statistical tests and controls,
- analyzes GO enrichment for biological interpretation.

## Main goals

1. Compare graph spectral clustering and hypergraph spectral clustering fairly.
2. Test whether penalized hypergraph regularization improves recovery.
3. Provide a reproducible reference pipeline with figures, analysis outputs, and tests.
4. Enable downstream paper-style benchmarking and validation.

## Repository structure

- hgspectral/: core Python package implementing parsing, clustering, metrics, controls, enrichment, and statistics.
- configs/: experiment configuration files.
- scripts/: helper scripts for synthetic data generation and pipeline setup.
- tests/: regression tests for pipeline behavior.
- fin_result/: outputs from completed benchmark runs.
- figures/: generated figures.
- datasets/: place your MITAB datasets here in the expected folder structure.

## Scientific framing

The key comparison is between:

- graph spectral clustering on a clique-expanded graph,
- unregularized hypergraph spectral clustering,
- penalized hypergraph spectral clustering with a degree-aware penalty.

The pipeline is designed to avoid common evaluation pitfalls such as leakage from test data, inconsistent cluster counts across methods, and unfair selection of parameters.

## Core workflow

1. Collect MITAB datasets and place them under a root directory.
2. Run the pipeline with a chosen configuration.
3. Evaluate held-out complex recovery with multiple seeds.
4. Compare methods using paired statistical tests and null models.
5. Generate summary tables and plots for downstream analysis.

## Installation

Use either conda or pip:

```bash
conda env create -f environment.yml
conda activate hgspectral
```

or:

```bash
pip install -r requirements.txt
```

## Quick start

```bash
PYTHONPATH=. python run_pipeline.py --datasets SynthA --search-roots . --seeds 6 --out smoke_out
```

For real datasets:

```bash
PYTHONPATH=. python run_pipeline.py --config configs/full_mitab.json --search-roots <data-root>
```

## Outputs

The pipeline writes per-dataset result files such as:

- test metrics,
- controls,
- GO enrichment summaries,
- manifests describing dataset characteristics,
- figures and analysis tables.

## Notes

This project is intended as a benchmark and reproducibility package. The repository includes both a smoke-test path and a full pipeline for larger runs.
