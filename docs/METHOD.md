# Method, conventions and provenance

For an allowed path (a,b,c), use the ordinary boundary
(b,c) − (a,c) + (a,b) in the space of **all ordered vertex pairs**, including
(a,a). Let A2 be the span of allowed two-paths and Omega2 those whose boundary
has no forbidden pairs. H1 = ker(partial1) / partial2(Omega2).

Each forbidden endpoint pair (a,c) imposes one equation: its two-path coefficients
sum to zero. Differences from one fixed middle vertex therefore span that
constraint's kernel. When (a,c) is an edge, each two-path contributes its triangle
boundary separately. These groups have disjoint path coordinates, so the generated
short boundaries span the complete two-boundary image. Sparse field elimination
gives its rank. A forest of the underlying edge multigraph supplies integer
fundamental cycles; extending the boundary basis by these cycles selects
representatives. The Betti dimension is m−n+components−rank(partial2).

The separate verifier instead constructs both allowed and forbidden boundary
matrices, solves the forbidden nullspace by dense elimination, and reconstructs
the image. Tests exhaust every three-vertex loopless digraph over Q, F2, F3 and
F5, then check seeded four-vertex graphs, relabeling, false certificates, exact
rational rank controls and resource/error paths. This finite coverage is not a
proof about arbitrary inputs or a formal verification of Python.

Path homology originates in Grigor'yan, Lin, Muranov and Yau,
[Homologies of path complexes and digraphs](https://arxiv.org/abs/1207.2834).
Existing computational work includes Chowdhury, Huntsman and Yutin,
[Path homologies of motifs and temporal network representations](https://doi.org/10.1007/s41109-021-00441-z).
The package's short-boundary implementation was adapted from the user's AIM-649
and AIM-650 computational seeds; no third-party repository code was copied.
Exact signed elimination principles also informed the AIM-531 computational work.
The software contribution is the bounded reusable interface, representative and
replay contract, validation and packaging; no new mathematical priority is claimed.

## Direct-use result

For a finite loopless digraph on n labelled vertices with m sorted directed
edges, and coefficients Q or a prime field F_p, `compute_h1` returns the exact
finite Betti dimension m−n+c−r and, when requested, cycle representatives for
that quotient. Here c is the number of weak components, including isolated
vertices, and r is the rank of the complete non-regular two-boundary image.
No diagonal ordered pair is discarded. The algorithm returns no partial answer
after a configured resource limit; its certificate replay has a separate smaller
dense-storage limit. Exact arithmetic here is a computational contract, not a
formal verification or a statement about random-graph asymptotics.
