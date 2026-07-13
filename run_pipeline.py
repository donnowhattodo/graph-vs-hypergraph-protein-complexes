"""
Reproducible pipeline runner.

Usage:
  python run_pipeline.py --config configs/default.json
  python run_pipeline.py --datasets Cardiac BioCreative --seeds 20

Point --search-roots at the directory that holds <Dataset>/<Dataset>.txt.
"""
import argparse
from hgspectral.config import Config
from hgspectral import pipeline


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    ap.add_argument("--datasets", nargs="*", default=None)
    ap.add_argument("--search-roots", nargs="*", default=None)
    ap.add_argument("--seeds", type=int, default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--no-controls", action="store_true")
    args = ap.parse_args()

    cfg = Config.load(args.config) if args.config else Config()
    if args.datasets:
        cfg.dataset_names = args.datasets
    if args.search_roots:
        cfg.search_roots = args.search_roots
    if args.seeds:
        cfg.n_seeds = args.seeds
    if args.out:
        cfg.out_root = args.out
    if args.no_controls:
        cfg.run_controls = False

    from pathlib import Path
    Path(cfg.out_root).mkdir(parents=True, exist_ok=True)
    cfg.save(Path(cfg.out_root) / "config_used.json")
    pipeline.run_all(cfg)


if __name__ == "__main__":
    main()
