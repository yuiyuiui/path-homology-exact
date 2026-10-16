"""Independent dense ambient-boundary replay of a small finite certificate.

This path does not use the production short-boundary generator or sparse basis.
It intentionally pays more time/memory and has its own explicit admission cap.
"""

from .core import Digraph, Field, ResourceLimitError


def rref(rows, width, field):
    rows = [list(map(field.scalar, row)) for row in rows]
    pivots = []
    for j in range(width):
        pivot = next((i for i in range(len(pivots), len(rows)) if rows[i][j]), None)
        if pivot is None:
            continue
        k = len(pivots)
        rows[k], rows[pivot] = rows[pivot], rows[k]
        inverse = field.inverse(rows[k][j])
        rows[k] = [field.clean(x * inverse) for x in rows[k]]
        for i in range(len(rows)):
            if i != k and rows[i][j]:
                factor = rows[i][j]
                rows[i] = [field.clean(a - factor * b) for a, b in zip(rows[i], rows[k])]
        pivots.append(j)
    return rows, pivots


def verify_certificate(record, *, max_matrix_entries=2000000):
    """Return False for a malformed/false certificate; raise on the size cap.

    Recompute Omega_2 as a nullspace in all ordered pairs, including (a,a).
    A certificate must contain a complete independent representative basis.
    """
    try:
        if record["format_version"] != 1 or record["convention"] != "non-regular":
            return False
        graph = Digraph(record["n"], record["edges"])
        field = Field(record["characteristic"])
        n, edges = graph.n, graph.edges
        ei = {e: j for j, e in enumerate(edges)}
        if max(n * n, len(edges) ** 2) > max_matrix_entries:
            raise ResourceLimitError("Independent dense replay exceeds max_matrix_entries")
        paths = [(a, b, d) for a, b in edges for c, d in edges if b == c]
        k, m = len(paths), len(edges)
        if max(n * n * k, m * k, m * m, k * k) > max_matrix_entries:
            raise ResourceLimitError("Independent dense replay exceeds max_matrix_entries")
        illegal_pairs = [(a, b) for a in range(n) for b in range(n) if (a, b) not in ei]
        fi = {pair: j for j, pair in enumerate(illegal_pairs)}
        illegal = [[0] * k for _ in fi]
        allowed = [[0] * k for _ in edges]
        incidence = [[0] * m for _ in range(n)]
        for j, (a, b) in enumerate(edges):
            incidence[a][j] -= 1
            incidence[b][j] += 1
        for j, (a, b, c) in enumerate(paths):
            for pair, sign in [((b, c), 1), ((a, c), -1), ((a, b), 1)]:
                if pair in ei:
                    allowed[ei[pair]][j] += sign
                else:
                    illegal[fi[pair]][j] += sign
        rows, pivots = rref(illegal, k, field)
        image = []
        for free in set(range(k)) - set(pivots):
            vector = [field.scalar(0)] * k
            vector[free] = field.scalar(1)
            for i, pivot in enumerate(pivots):
                vector[pivot] = field.clean(-rows[i][free])
            image.append([field.clean(sum(x * y for x, y in zip(row, vector))) for row in allowed])
        rank = len(rref(image, m, field)[1])
        dim = m - len(rref(incidence, m, field)[1])
        beta = dim - rank
        if (record["betti1"], record["boundary_rank"], record["cycle_dimension"]) != (
            beta,
            rank,
            dim,
        ):
            return False
        representatives = record["representatives"]
        if len(representatives) != beta:
            return False
        for terms in representatives:
            value = [field.scalar(0)] * m
            seen = set()
            for j, c in terms:
                if not isinstance(j, int) or isinstance(j, bool) or not 0 <= j < m or j in seen:
                    return False
                seen.add(j)
                value[j] = field.scalar(c)
            if any(field.clean(sum(x * y for x, y in zip(row, value))) for row in incidence):
                return False
            image.append(value)
        return len(rref(image, m, field)[1]) == rank + beta
    except (KeyError, TypeError, ValueError, ZeroDivisionError, IndexError):
        return False
