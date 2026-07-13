"""
MITAB parsing and hyperedge (pull-down complex) construction.

The parsing / hyperedge logic is intentionally kept close to the authors'
original notebooks so behaviour on real files is unchanged; the functions are
just factored out, documented, and made unit-testable.
"""
from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

# MITAB 2.7 column indices used by the authors' files.
COLS = {
    "A": 0, "B": 1, "pub_ids": 8, "taxA": 9, "taxB": 10,
    "interaction_type": 11, "source_db": 12, "int_id": 13,
    "confidence": 14, "expansion": 15,
    "bio_role_A": 16, "bio_role_B": 17, "exp_role_A": 18, "exp_role_B": 19,
    "xrefA": 22, "xrefB": 23, "host": 27, "creation_date": 30, "update_date": 31,
}
USECOLS = sorted(set(COLS.values()))
NAMES = {v: k for k, v in COLS.items() if v in USECOLS}

TAXID_PATTERN = re.compile(r"taxid:(-?\d+)")
GO_PATTERN = re.compile(r'go:"(GO:\d+)"\(([^)]*)\)', re.IGNORECASE)


def clean_uniprot(idstr: str, collapse_isoforms: bool = True) -> str:
    if not isinstance(idstr, str) or idstr == "-":
        return "-"
    first = idstr.split("|")[0]
    m = re.match(r"^([^:]+):([A-Za-z0-9]+)(-\d+)?$", first)
    if not m:
        return first
    db, acc, iso = m.groups()
    return f"{db}:{acc}" if collapse_isoforms else f"{db}:{acc}{iso or ''}"


def extract_primary_taxid(cell) -> Optional[str]:
    m = TAXID_PATTERN.search(str(cell))
    return m.group(1) if m else None


def parse_roles(cell) -> set:
    s = str(cell).lower()
    roles = {tag for tag in ("bait", "prey") if f"({tag})" in s}
    if not roles and "unspecified role" in s:
        roles.add("unspecified")
    return roles or {"unspecified"}


def is_pull_down_like(row: dict) -> bool:
    det = str(row.get("interaction_type", "-")).lower()
    meth = str(row.get("expansion", "-")).lower()
    rolesA = parse_roles(row.get("exp_role_A", "-"))
    rolesB = parse_roles(row.get("exp_role_B", "-"))
    baitlike = ("bait" in rolesA) or ("bait" in rolesB)
    complexlike = ("spoke" in meth) or ("matrix" in meth)
    phys_assoc = ("physical association" in det)
    return baitlike or complexlike or phys_assoc


def experiment_key(row: dict, bait: str) -> tuple:
    return (row.get("pub_ids", "-"), row.get("host", "-"),
            row.get("interaction_type", "-"), row.get("expansion", "-"), bait)


def extract_go_terms(cell) -> List[Tuple[str, str]]:
    out = []
    for m in GO_PATTERN.finditer(str(cell)):
        out.append((m.group(1).upper(), m.group(2).strip()))
    return out


def find_mitab_file(dataset_name: str, search_roots: Optional[List[str]] = None) -> Optional[Path]:
    roots = [Path(r) for r in (search_roots or [".", "HG_Human"])]
    cands = []
    for r in roots:
        cands += [r / dataset_name / f"{dataset_name}.txt", r / f"{dataset_name}.txt"]
    cands.append(Path(dataset_name))  # allow passing a direct path
    for p in cands:
        if p.exists():
            return p
    return None


def load_mitab(path: Path, top_k_taxids: int = 3) -> Tuple[pd.DataFrame, dict]:
    """Load a MITAB file, apply a top-k taxid filter, return (df, stats)."""
    df = pd.read_csv(path, sep="\t", header=None, comment="#",
                     low_memory=False, usecols=USECOLS, dtype=str)
    df = df.rename(columns=NAMES).fillna("-")
    df["A"] = df["A"].map(clean_uniprot)
    df["B"] = df["B"].map(clean_uniprot)
    df["taxidA"] = df["taxA"].map(extract_primary_taxid)
    df["taxidB"] = df["taxB"].map(extract_primary_taxid)

    stats = {"n_ppi_raw": int(len(df))}
    tax = pd.concat([df["taxidA"], df["taxidB"]], ignore_index=True).dropna()
    counts = tax.value_counts()
    if not counts.empty and top_k_taxids:
        top = list(counts.head(top_k_taxids).index)
        mask = df["taxidA"].isin(top) | df["taxidB"].isin(top)
        df = df[mask].copy()
        stats["top_taxids"] = top
    stats["n_ppi_filtered"] = int(len(df))
    df = df[(df["A"] != "-") & (df["B"] != "-")]
    return df, stats


def mitab_to_complexes(df: pd.DataFrame, min_size: int = 2) -> Dict[str, List[str]]:
    """Group pull-down-like rows into hyperedges keyed by experiment."""
    complexes = defaultdict(set)
    for row in df.itertuples(index=False):
        r = row._asdict()
        if not is_pull_down_like(r):
            continue
        idA, idB = r["A"], r["B"]
        rolesA, rolesB = parse_roles(r.get("exp_role_A", "-")), parse_roles(r.get("exp_role_B", "-"))
        if "bait" in rolesA and "prey" in rolesB:
            bait, prey = idA, idB
        elif "bait" in rolesB and "prey" in rolesA:
            bait, prey = idB, idA
        elif "bait" in rolesA and "bait" not in rolesB:
            bait, prey = idA, idB
        elif "bait" in rolesB and "bait" not in rolesA:
            bait, prey = idB, idA
        else:
            bait, prey = (idA, idB) if idA <= idB else (idB, idA)
        key = experiment_key(r, bait)
        complexes[key].update((bait, prey))
    out = {}
    for i, (_, members) in enumerate(complexes.items()):
        if len(members) >= min_size:
            out[f"exp_{i}"] = sorted(members)
    return out


def build_protein_go(df: pd.DataFrame, proteins=None):
    prot2go = defaultdict(set)
    go_id_to_name = {}
    for row in df.itertuples(index=False):
        r = row._asdict()
        for side, xref in (("A", "xrefA"), ("B", "xrefB")):
            pid = r.get(side, "-")
            if pid in (None, "-"):
                continue
            for go_id, go_name in extract_go_terms(r.get(xref, "-")):
                prot2go[pid].add(go_id)
                if go_id not in go_id_to_name and go_name:
                    go_id_to_name[go_id] = go_name
    if proteins is not None:
        ps = set(proteins)
        prot2go = {p: g for p, g in prot2go.items() if p in ps and g}
    return dict(prot2go), go_id_to_name
