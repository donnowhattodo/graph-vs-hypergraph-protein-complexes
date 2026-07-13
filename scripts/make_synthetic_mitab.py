"""
Generate a SYNTHETIC MITAB 2.7 file with planted higher-order (complex)
structure, so the pipeline can be smoke-tested without real data.

*** The numbers this produces are NOT real results and must never appear in
the manuscript. It exists only to verify the code runs end-to-end. ***
"""
import random
import sys
from pathlib import Path

TAXIDS = ["9606", "10090", "10116"]
GO_POOL = [f"GO:{7000000 + i:07d}" for i in range(40)]
GO_NAMES = {g: f"synthetic process {i}" for i, g in enumerate(GO_POOL)}


def make(path, n_proteins=180, n_complexes=120, seed=0):
    rng = random.Random(seed)
    proteins = [f"P{idx:05d}" for idx in range(n_proteins)]
    # give each protein a taxid and a small GO signature tied to a latent module
    n_modules = 8
    module_of = {p: rng.randrange(n_modules) for p in proteins}
    tax_of = {p: rng.choice(TAXIDS) for p in proteins}
    go_of = {}
    for p in proteins:
        base = module_of[p] * 3
        go_of[p] = [GO_POOL[(base + j) % len(GO_POOL)] for j in range(3)]

    def xref(p):
        return "|".join(f'go:"{g}"({GO_NAMES[g]})' for g in go_of[p])

    lines = []
    for c in range(n_complexes):
        module = rng.randrange(n_modules)
        members = [p for p in proteins if module_of[p] == module]
        if len(members) < 3:
            members = rng.sample(proteins, 4)
        size = rng.randint(3, min(8, len(members)))
        chosen = rng.sample(members, size)
        bait = chosen[0]
        pub = f"pubmed:{1000 + c}"
        for prey in chosen[1:]:
            cols = ["-"] * 32
            cols[0] = f"uniprotkb:{bait}"
            cols[1] = f"uniprotkb:{prey}"
            cols[8] = pub
            cols[9] = f"taxid:{tax_of[bait]}(name)"
            cols[10] = f"taxid:{tax_of[prey]}(name)"
            cols[11] = 'psi-mi:"MI:0915"(physical association)'
            cols[15] = 'psi-mi:"MI:1060"(spoke expansion)'
            cols[18] = 'psi-mi:"MI:0496"(bait)'
            cols[19] = 'psi-mi:"MI:0498"(prey)'
            cols[22] = xref(bait)
            cols[23] = xref(prey)
            cols[27] = "taxid:9606(host)"
            lines.append("\t".join(cols))
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("\n".join(lines) + "\n")
    print(f"[SYNTHETIC] wrote {len(lines)} PPI rows -> {path}")


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "SynthA/SynthA.txt"
    make(out, seed=0)
    make("SynthB/SynthB.txt", n_proteins=140, n_complexes=90, seed=1)
