# Contract and API

- `Digraph(n, edges, limits=Limits())`: vertices 0..n−1, n>=0; integral indices,
  self-loops rejected, duplicate edges coalesced, sorted immutable edge tuple.
  Input storage is copied; no caller graph is mutated. Isolated vertices count.
- `compute_h1(graph, characteristic=0, limits=Limits(), representatives=True)`:
  exact arithmetic over Q or F_p for a prime p<=2^31−1. No float rank threshold,
  randomized algorithm, global thread/precision change or numerical tolerance.
  Returns a frozen `Homology` record with `betti`, `boundary_rank`,
  `cycle_dimension`, original graph and immutable representative tuples.
- A chain is a tuple of `(edge_index, coefficient)` pairs. Representatives are
  cycles whose cosets form a basis modulo the two-boundaries; they are not a
  preferred geometric or shortest basis. Rational coefficients use `Fraction`;
  finite-field coefficients are residues 0..p−1.
- `betti1(...)` requests dimension only, avoiding representative elimination.
  `compute_h1(..., representatives=False)` has no complete certificate when the
  Betti number is positive; its dimension remains exact.
- `result.certificate()` returns detached JSON-compatible data with rational
  coefficients as strings. `verify_certificate(record, max_matrix_entries=2_000_000)`
  rebuilds the full ambient boundary/nullspace, checks the ranks, divergence,
  independence and spanning of the supplied representatives. It does not trust
  the claimed ranks or use the production short-generator routine. False or
  malformed records return False; excessive replay size raises `ResourceLimitError`.

`Limits` defaults: 20000 input edges, 200000 allowed two-paths, 1000000 stored
elimination entries and 2000000 scalar reduction updates. Graph admission also
limits n<=10000. Limits can be reduced or explicitly increased where available;
no partial result is returned after exhaustion. Exact rational elimination can
have substantial fill and coefficient growth. These are admission/work guards,
not a hard OS memory sandbox or a universal time bound. The dense verifier has
quadratic storage and is intended for smaller graphs than the production path.
Invalid vertices, composite characteristics and invalid limits raise ValueError.
No graph-theoretic asymptotic guarantee is inferred from a finite calculation.

Run `python examples/cycle.py` after installation. For development install
`.[dev]`, run `python -m pytest`, and build distributions with `python -m build`.
