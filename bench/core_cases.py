"""Same exact rational Betti dimension, two independent algorithms."""

import itertools
import random
from path_homology import Digraph, betti1
import sympy as sp


def ambient_betti(graph):
    edges = graph.edges
    ei = {e: j for j, e in enumerate(edges)}
    paths = [(a, b, d) for a, b in edges for c, d in edges if b == c]
    illegal_pairs = [e for e in itertools.product(range(graph.n), repeat=2) if e not in ei]
    fi = {e: j for j, e in enumerate(illegal_pairs)}
    illegal, allowed = sp.zeros(len(fi), len(paths)), sp.zeros(len(edges), len(paths))
    for j, (a, b, c) in enumerate(paths):
        for edge, sign in [((b, c), 1), ((a, c), -1), ((a, b), 1)]:
            if edge in ei:
                allowed[ei[edge], j] += sign
            else:
                illegal[fi[edge], j] += sign
    kernel = illegal.nullspace()
    image = allowed * sp.Matrix.hstack(*kernel) if kernel else sp.zeros(len(edges), 0)
    incidence = sp.zeros(graph.n, len(edges))
    for j, (a, b) in enumerate(edges):
        incidence[a, j] -= 1
        incidence[b, j] += 1
    return len(edges) - incidence.rank() - image.rank()


def make(case, variant):
    n = {"small": 5, "medium": 8, "large": 12}[case]
    rng = random.Random(649 + n)
    edges = [e for e in itertools.permutations(range(n), 2) if rng.random() < 0.28]
    graph = Digraph(n, edges)
    reference = ambient_betti(graph)

    def call():
        return ambient_betti(graph) if variant == "reference" else betti1(graph)

    def validate(result):
        if result != reference:
            raise ArithmeticError("Independent ambient Q homology disagrees")
        return {"passed": True, "betti1": result, "coefficient_field": "Q"}

    return (
        call,
        validate,
        {
            "n": n,
            "edges": len(edges),
            "seed": 649 + n,
            "requested_output": "Betti dimension only; no representatives",
            "variant": variant,
        },
    )
