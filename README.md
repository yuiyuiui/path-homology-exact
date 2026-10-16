# Path Homology Exact

Compute **non-regular first path homology** of finite loopless directed graphs
over the rationals or a prime field. Return a Betti dimension, representative
cycles in the original directed-edge coordinates, and a JSON certificate that
can be replayed by a separate dense ambient-boundary algorithm.

```sh
python -m pip install 'git+https://github.com/yuiyuiui/path-homology-exact.git'
```

```python
from path_homology import Digraph, compute_h1, verify_certificate

graph = Digraph(3, [(0, 1), (1, 2), (2, 0)])
result = compute_h1(graph, characteristic=0)  # 0 means Q; a prime means F_p.
assert result.betti == 1
assert verify_certificate(result.certificate())
print(result.representatives)  # Edge indices refer to graph.edges (sorted).
```

Use `betti1(graph, characteristic=2)` when only the dimension is needed. Python
>=3.11 is required; the library and certificate verifier have **no runtime
dependencies**. SymPy is only a development/benchmark reference.

## When this package helps

Researchers comparing directed-network invariants need precise coefficient and
regularity conventions, signed representatives, and independently reproducible
small calculations. This package provides that focused interface in pure Python.
It does not provide persistence, higher-dimensional homology, weighted shortest-
path filtrations, or a large-network backend. For grounded persistent homology,
consider [GrPPHATI](https://github.com/tomchaplin/grpphati). Non-regular homology
already has implementations, including
[PerformantPathHomology](https://github.com/SteveHuntsmanBAESystems/PerformantPathHomology).
This package makes no first-implementation or state-of-the-art speed claim.

**Non-regular matters:** the digraph with both 0→1 and 1→0 has Betti number 1
here. Its forbidden diagonal pairs are retained in the ambient boundary test;
a regular convention can give a different answer. A rank deficit over F_2 does
not establish a rational class. Always record the coefficient field.

See [contract and API](docs/API.md), [algorithm and provenance](docs/METHOD.md),
[benchmarks](bench/README.md), [readiness](docs/RELEASE_READINESS.md),
[contributing](CONTRIBUTING.md) and [citation](CITATION.cff).

The implementation, tests and documentation were developed with AI assistance
under Kaiwen Jin's direction. Software tests do not imply human peer review,
formal verification or a resolution of an asymptotic random-graph conjecture.
