import copy
import itertools
import json
import random

import pytest
from path_homology import (
    Digraph,
    Limits,
    ResourceLimitError,
    betti1,
    compute_h1,
    verify_certificate,
)


def test_nonregular_reciprocal_pair_retains_a_class():
    result = compute_h1(Digraph(2, [(0, 1), (1, 0)]))
    assert result.betti == 1 and result.boundary_rank == 0
    assert verify_certificate(json.loads(json.dumps(result.certificate())))


@pytest.mark.parametrize("p", [0, 2, 3, 5])
def test_all_three_vertex_graphs_against_independent_ambient_replay(p):
    edges = list(itertools.permutations(range(3), 2))
    for mask in range(64):
        graph = Digraph(3, [e for i, e in enumerate(edges) if mask >> i & 1])
        result = compute_h1(graph, characteristic=p)
        assert verify_certificate(result.certificate())
        assert result.betti == betti1(graph, characteristic=p)


def test_seeded_four_vertex_graphs_and_vertex_relabeling():
    rng = random.Random(649)
    for _ in range(40):
        edges = [e for e in itertools.permutations(range(4), 2) if rng.random() < 0.45]
        perm = rng.sample(range(4), 4)
        for p in (0, 2, 7):
            graph = Digraph(4, edges)
            result = compute_h1(graph, characteristic=p)
            assert verify_certificate(result.certificate())
            assert result.betti == betti1(
                Digraph(4, [(perm[a], perm[b]) for a, b in edges]), characteristic=p
            )


def test_tampered_rank_cycle_and_spanning_claims_are_rejected():
    record = compute_h1(Digraph(3, [(0, 1), (1, 2), (2, 0)])).certificate()
    assert record["betti1"] == 1
    for key in ("betti1", "boundary_rank", "cycle_dimension"):
        wrong = copy.deepcopy(record)
        wrong[key] += 1
        assert not verify_certificate(wrong)
    wrong = copy.deepcopy(record)
    wrong["representatives"][0][0][1] = "2"
    assert not verify_certificate(wrong)
    wrong = copy.deepcopy(record)
    wrong["representatives"] = []
    assert not verify_certificate(wrong)


def test_empty_trees_triangles_and_resource_boundaries():
    assert betti1(Digraph(0)) == 0
    assert betti1(Digraph(4, [(0, 1), (2, 3)])) == 0
    assert betti1(Digraph(3, [(0, 1), (1, 2), (0, 2)])) == 0
    assert verify_certificate(compute_h1(Digraph(0)).certificate())
    with pytest.raises(ResourceLimitError):
        Digraph(2, [(0, 1), (1, 0)], limits=Limits(max_edges=1))
    with pytest.raises(ResourceLimitError):
        compute_h1(Digraph(3, itertools.permutations(range(3), 2)), limits=Limits(max_two_paths=1))
    with pytest.raises(ResourceLimitError):
        verify_certificate(
            compute_h1(Digraph(3, [(0, 1), (1, 2), (2, 0)])).certificate(), max_matrix_entries=1
        )


@pytest.mark.parametrize(
    "n,edges", [(True, []), (-1, []), (2, [(0, 0)]), (2, [(0, 2)]), (2, [(0.0, 1)]), (2, [(1,)])]
)
def test_invalid_graphs(n, edges):
    with pytest.raises(ValueError):
        Digraph(n, edges)


@pytest.mark.parametrize("p", [1, 4, -1, True, 2.0, 2147483648])
def test_invalid_coefficients(p):
    with pytest.raises(ValueError):
        betti1(Digraph(0), characteristic=p)


def test_sparse_rational_elimination_against_independent_matrix_rank():
    from fractions import Fraction
    import sympy as sp
    from path_homology.core import Basis, Field

    rng = random.Random(531)
    for _ in range(20):
        rows = [
            [Fraction(rng.randint(-3, 3), rng.randint(1, 5)) for _ in range(6)] for _ in range(8)
        ]
        basis = Basis(Field(0), Limits())
        for row in rows:
            basis.add(dict(enumerate(row)))
        assert len(basis.rows) == sp.Matrix(rows).rank()


def test_dense_replay_refuses_large_empty_ambient_space_before_allocation():
    record = compute_h1(Digraph(2000)).certificate()
    with pytest.raises(ResourceLimitError):
        verify_certificate(record)
