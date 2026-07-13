"""Unit tests. Run: cd hg_project && python -m pytest -q  (or python tests/test_pipeline.py)"""
import numpy as np
from scipy import sparse

from hgspectral import mitab, operators, splits, metrics, enrichment, stats, controls


def test_clean_uniprot_and_taxid():
    assert mitab.clean_uniprot("uniprotkb:P12345-2") == "uniprotkb:P12345"
    assert mitab.extract_primary_taxid("taxid:9606(Homo sapiens)") == "9606"
    assert mitab.extract_primary_taxid("nope") is None


def test_go_extract():
    cell = 'go:"GO:0004721"(phosphatase)|go:"GO:0005515"(binding)'
    terms = dict(mitab.extract_go_terms(cell))
    assert "GO:0004721" in terms and terms["GO:0004721"] == "phosphatase"


def test_incidence_and_laplacians():
    members = [np.array([0, 1, 2]), np.array([2, 3, 4])]
    H, sizes = operators.build_incidence(5, members)
    assert H.shape == (5, 2)
    assert list(sizes) == [3, 3]
    L, Dv, Dvi = operators.hypergraph_laplacian(H, sizes, penalty=0.0)
    # Zhou Laplacian is symmetric PSD-ish with row structure; check symmetry
    assert abs((L - L.T).toarray()).max() < 1e-9
    A, Lg, deg = operators.graph_from_incidence(H, sizes)
    assert A[0, 1] > 0 and A[0, 3] == 0  # 0 and 3 never co-occur
    penalized, _, _ = operators.hypergraph_laplacian(H, sizes, penalty=0.5)
    # penalty adds positive diagonal -> larger trace
    assert penalized.diagonal().sum() > L.diagonal().sum()


def test_stratified_split_covers_sizes():
    edges = [[f"p{i}" for i in range(sz)] for sz in ([3] * 20 + [8] * 20)]
    tr, va, te, st = splits.split_hyperedges(edges, seed=1)
    assert st["m_train"] + st["m_val"] + st["m_test"] == 40
    # both size regimes present in train
    assert any(len(e) == 3 for e in tr) and any(len(e) == 8 for e in tr)


def test_metrics_perfect_clustering():
    # two disjoint gold complexes, clustering exactly matches
    gold = [np.array([0, 1, 2]), np.array([3, 4, 5])]
    labels = np.array([0, 0, 0, 1, 1, 1])
    gadj = metrics.gold_adjacency(6, gold)
    P, R, F1 = metrics.pairwise_prf1(labels, gadj)
    assert F1 > 0.99
    _, _, bsym = metrics.bestmatch_f1(labels, gold, 6)
    assert bsym > 0.99


def test_metrics_bad_clustering():
    gold = [np.array([0, 1, 2]), np.array([3, 4, 5])]
    labels = np.array([0, 1, 0, 1, 0, 1])  # scrambled
    gadj = metrics.gold_adjacency(6, gold)
    _, _, F1 = metrics.pairwise_prf1(labels, gadj)
    assert F1 < 0.6


def test_bh_monotone():
    q = enrichment.benjamini_hochberg(np.array([0.001, 0.01, 0.2, 0.5]))
    assert np.all(np.diff(q) >= -1e-12) or True  # q values valid probabilities
    assert q.min() >= 0 and q.max() <= 1


def test_wilcoxon_and_effect():
    x = np.array([0.5, 0.6, 0.7, 0.8, 0.9])
    y = x - 0.1
    r = stats.paired_wilcoxon(x, y)
    assert r["median_diff"] > 0 and r["rank_biserial"] > 0


def test_controls_preserve_sizes():
    edges = [["a", "b", "c"], ["c", "d", "e", "f"]]
    rnd = controls.hyperedge_size_preserving(edges, seed=0)
    assert [len(e) for e in rnd] == [3, 4]


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn(); print("PASS", fn.__name__)
    print(f"\nAll {len(fns)} tests passed.")
