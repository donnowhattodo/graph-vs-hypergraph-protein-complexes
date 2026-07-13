"""Configuration for the hypergraph-vs-graph benchmark pipeline."""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import List, Tuple
import json


@dataclass
class Config:
    dataset_names: List[str] = field(default_factory=lambda: ["Cardiac", "BioCreative"])
    search_roots: List[str] = field(default_factory=lambda: [".", "HG_Human1"])
    top_k_taxids: int = 3
    min_complex_size: int = 2

    # split
    split_fracs: Tuple[float, float, float] = (0.6, 0.2, 0.2)
    stratify_by_size: bool = True
    split_seed: int = 12345

    # model grids
    lambda_grid: Tuple[float, ...] = (0.0, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0)
    k_matched: str = "rule"          # 'rule' => n_train // min_cluster_size, or int
    min_cluster_size: int = 5
    adaptive_criteria: Tuple[str, ...] = ("eigengap", "silhouette", "modularity")
    kmax: int = 15

    # seeds
    n_seeds: int = 20
    base_seed: int = 42

    # controls
    n_control_seeds: int = 10
    run_controls: bool = True

    # baselines
    run_louvain: bool = True
    run_leiden: bool = True
    run_mcl: bool = True

    # enrichment
    q_threshold: float = 0.05
    min_term_size: int = 5

    out_root: str = "fin_result"

    def save(self, path):
        with open(path, "w") as f:
            json.dump(asdict(self), f, indent=2)

    @classmethod
    def load(cls, path):
        with open(path) as f:
            d = json.load(f)
        d["split_fracs"] = tuple(d["split_fracs"])
        d["lambda_grid"] = tuple(d["lambda_grid"])
        d["adaptive_criteria"] = tuple(d["adaptive_criteria"])
        return cls(**d)
