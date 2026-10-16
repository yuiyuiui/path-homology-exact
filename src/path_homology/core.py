"""Sparse short boundaries and exact field elimination, without floating rank tests."""

from collections import defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral


class ResourceLimitError(RuntimeError):
    """A configured work/storage limit was reached; no partial answer is returned."""


def integer(value, name, low=0):
    if isinstance(value, bool) or not isinstance(value, Integral) or value < low:
        raise ValueError(f"{name} must be an integer >= {low}")
    return int(value)


@dataclass(frozen=True)
class Limits:
    max_edges: int = 20000
    max_two_paths: int = 200000
    max_stored_entries: int = 1000000
    max_reduction_steps: int = 2000000

    def __post_init__(self):
        for key, value in vars(self).items():
            integer(value, key, 1)


class Field:
    def __init__(self, characteristic):
        p = integer(characteristic, "characteristic")
        if p and (p < 2 or p > 2147483647 or any(p % d == 0 for d in range(2, int(p**0.5) + 1))):
            raise ValueError("Use characteristic 0 (Q) or a prime <= 2^31-1")
        self.p = p

    def scalar(self, x):
        if isinstance(x, int):
            return x % self.p if self.p else x
        q = x if isinstance(x, Fraction) else Fraction(x)
        if self.p:
            if q.denominator % self.p == 0:
                raise ValueError("Denominator is zero in the coefficient field")
            return q.numerator * pow(q.denominator, -1, self.p) % self.p
        return q.numerator if q.denominator == 1 else q

    def inverse(self, x):
        return pow(x, -1, self.p) if self.p else Fraction(1, x)

    def clean(self, x):
        return (
            x % self.p
            if self.p
            else (x.numerator if isinstance(x, Fraction) and x.denominator == 1 else x)
        )


@dataclass(frozen=True)
class Digraph:
    """Vertices are 0..n-1. Duplicate edges are coalesced; self-loops are rejected."""

    n: int
    edges: tuple

    def __init__(self, n, edges=(), *, limits=Limits()):
        n = integer(n, "n")
        if n > 10000:
            raise ResourceLimitError("At most 10000 vertices are admitted")
        checked = set()
        for count, edge in enumerate(edges, 1):
            if count > limits.max_edges:
                raise ResourceLimitError("Input edge count exceeds max_edges")
            try:
                a, b = edge
            except (TypeError, ValueError) as exc:
                raise ValueError("Each edge must contain two vertex indices") from exc
            a, b = integer(a, "vertex"), integer(b, "vertex")
            if a >= n or b >= n or a == b:
                raise ValueError("Edges must join distinct vertices in 0..n-1")
            checked.add((a, b))
        object.__setattr__(self, "n", n)
        object.__setattr__(self, "edges", tuple(sorted(checked)))


def boundary_columns(graph, limits):
    """Generate a spanning image of Omega_2, retaining forbidden diagonal pairs."""
    index = {e: j for j, e in enumerate(graph.edges)}
    outgoing = [[] for _ in range(graph.n)]
    for a, b in graph.edges:
        outgoing[a].append(b)
    groups = defaultdict(list)
    total = 0
    for a, b in graph.edges:
        for c in outgoing[b]:
            total += 1
            if total > limits.max_two_paths:
                raise ResourceLimitError("Allowed two-path count exceeds max_two_paths")
            groups[a, c].append(b)
    for (a, c), middles in sorted(groups.items()):

        def row(b):
            result = defaultdict(int)
            result[index[a, b]] += 1
            result[index[b, c]] += 1
            if (a, c) in index:
                result[index[a, c]] -= 1
            return dict(result)

        if (a, c) in index:
            for b in middles:
                yield row(b)
        else:
            first = row(middles[0])
            for b in middles[1:]:
                value = row(b)
                for j, coefficient in first.items():
                    value[j] = value.get(j, 0) - coefficient
                yield {j: v for j, v in value.items() if v}


class Basis:
    def __init__(self, field, limits):
        self.field, self.limits = field, limits
        self.rows = {}
        self.entries = self.steps = 0

    def add(self, column):
        f = self.field
        value = {j: c for j, v in column.items() if (c := f.scalar(v))}
        while value:
            pivot = max(value)
            if pivot not in self.rows:
                scale = f.inverse(value[pivot])
                if scale != 1:
                    value = {j: f.clean(v * scale) for j, v in value.items()}
                if self.entries + len(value) > self.limits.max_stored_entries:
                    raise ResourceLimitError("Elimination fill exceeds max_stored_entries")
                self.rows[pivot] = value
                self.entries += len(value)
                return True
            factor = value[pivot]
            for j, v in self.rows[pivot].items():
                self.steps += 1
                if self.steps > self.limits.max_reduction_steps:
                    raise ResourceLimitError("Elimination exceeds max_reduction_steps")
                product = v if factor == 1 else -v if factor == -1 else factor * v
                updated = f.clean(value.get(j, 0) - product)
                if updated:
                    value[j] = updated
                else:
                    value.pop(j, None)
        return False


def cycles(graph):
    """Integer fundamental cycles of the underlying edge multigraph."""
    parent = list(range(graph.n))

    def find(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v

    forest = [[] for _ in range(graph.n)]
    chords = []
    for j, (a, b) in enumerate(graph.edges):
        x, y = find(a), find(b)
        if x == y:
            chords.append((j, a, b))
        else:
            parent[x] = y
            forest[a].append((b, j, 1))
            forest[b].append((a, j, -1))
    components = len({find(v) for v in range(graph.n)})

    def iterator():
        for j, a, b in chords:
            todo = deque([b])
            previous = {b: None}
            while todo and a not in previous:
                v = todo.popleft()
                for w, edge, sign in forest[v]:
                    if w not in previous:
                        previous[w] = v, edge, sign
                        todo.append(w)
            value = {j: 1}
            v = a
            while v != b:
                old, edge, sign = previous[v]
                value[edge] = sign
                v = old
            yield value

    return components, iterator()


@dataclass(frozen=True)
class Homology:
    graph: Digraph
    characteristic: int
    betti: int
    boundary_rank: int
    cycle_dimension: int
    representatives: tuple

    def certificate(self):
        return {
            "format_version": 1,
            "convention": "non-regular",
            "n": self.graph.n,
            "edges": [list(e) for e in self.graph.edges],
            "characteristic": self.characteristic,
            "betti1": self.betti,
            "boundary_rank": self.boundary_rank,
            "cycle_dimension": self.cycle_dimension,
            "representatives": [[[j, str(c)] for j, c in chain] for chain in self.representatives],
        }


def compute_h1(graph, *, characteristic=0, limits=Limits(), representatives=True):
    """Exact H1 over Q (0) or F_p; columns use the sorted edge order in graph."""
    if not isinstance(graph, Digraph):
        raise TypeError("graph must be a Digraph")
    if len(graph.edges) > limits.max_edges:
        raise ResourceLimitError("Graph exceeds max_edges")
    if not isinstance(representatives, bool):
        raise ValueError("representatives must be boolean")
    field = Field(characteristic)
    basis = Basis(field, limits)
    for column in boundary_columns(graph, limits):
        basis.add(column)
    rank = len(basis.rows)
    components, candidates = cycles(graph)
    dim = len(graph.edges) - graph.n + components
    answer = []
    if representatives:
        for cycle in candidates:
            if basis.add(cycle):
                answer.append(tuple(sorted((j, field.scalar(v)) for j, v in cycle.items())))
        if len(answer) != dim - rank:
            raise ArithmeticError("Cycle/image rank invariant failed")
    return Homology(graph, field.p, dim - rank, rank, dim, tuple(answer))


def betti1(graph, *, characteristic=0, limits=Limits()):
    """Compute only the dimension, avoiding representative-cycle elimination."""
    return compute_h1(
        graph, characteristic=characteristic, limits=limits, representatives=False
    ).betti
