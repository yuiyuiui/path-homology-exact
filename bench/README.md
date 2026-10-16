# Reproducible core benchmarks

[Raw results](results/core_v0.1.0.json) include five warm samples per variant,
first-call and setup times, source hashes, profiles, validation, settings and
dependency versions. macOS arm64, CPython 3.13.3, one numerical worker and
BLAS/OpenMP thread. Variants alternate execution order between cases and run in
separate fresh processes. No random input seed is hidden.

| Case | Reference median s | Production median s | Ratio | Python peak B (ref/prod) | Process RSS B (ref/prod) |
|---|---:|---:|---:|---:|---:|
| small | 0.000148125 | 5.12495e-06 | 28.90x | 12622 / 2376 | 73056256 / 72728576 |
| medium | 0.00306513 | 5.8667e-05 | 52.25x | 140004 / 10032 | 73433088 / 72826880 |
| large | 0.0187015 | 0.000231791 | 80.68x | 395252 / 22432 | 73629696 / 73416704 |

Both paths request only the exact rational Betti dimension, with the same graph and field. The reference constructs the full ambient forbidden-pair nullspace and uses SymPy exact ranks; it is deliberately simple, not an optimized competing topology library. The production path uses short boundaries and sparse elimination. Profiling then identified avoidable Fraction creation and unit-pivot arithmetic: retaining integral intermediates reduced that work. The final profile and raw source hash are recorded below. No superiority over GrPPHATI or other optimized topology software is claimed. Representative generation and certificate replay are tested separately, not included in this dimension-only timing.

Validation is excluded from timed calls, repeated after every sample, and must
agree with the independent reference. Separate profiling and tracemalloc calls
do not contribute to timing samples. Setup includes imports, case construction
and the independent reference; it is not a pure production startup measurement.
Process RSS includes these reference/validation dependencies. Python-traced peak
measures one separate call and excludes untraced native memory. Small RSS changes
are noise; no universal process-memory reduction is claimed. These finite cases
do not establish asymptotic speedups or performance on arbitrary hardware.

```sh
python -m pip install '.[dev]'
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python bench/run.py \
  --repeats 5 --output /tmp/core-comparison.json
```

The output path must be new. The driver runs bounded small/medium/large cases;
performance thresholds are intentionally absent from shared CI. No network,
paid backend, GPU or hidden input data is needed. Each JSON `settings` object
specifies its requested output and cost boundary. Keep accuracy settings fixed
when comparing variants.
