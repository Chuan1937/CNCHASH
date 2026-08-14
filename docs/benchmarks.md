# Benchmarks

## Accuracy (synthetic recovery)

![Accuracy Verification](images/accuracy_verification.png)

Synthetic events with uniform full-focal-sphere coverage: without
polarity errors CNCHASH recovers 10/10 known mechanisms (median Kagan
rotation ~8 deg, the dang=5 grid discretization floor); with ~10%
flipped polarities 8/10 remain within 25 deg. These checks are part of
the test suite (tests/test_accuracy.py).

**Note:** Strike differences (40-80 deg) are normal - focal mechanisms
have two orthogonal nodal planes that both satisfy polarity data.

Historical notebook results (HASH_Tests.ipynb): CNCHASH synthetic
trials 300/300 successful solutions; HASH v1.2 60/60 (direct-runs).

## Native backend measurements

Run locally with:

```bash
python benchmarks/benchmark_backends.py --events 1000 --seconds 2
```

Measured on a development machine (x86_64, gfortran). Values are for
30-station events with 30 MC trials (`dang=5`).

## Thread scaling (portable build, 30 stations)

| Threads | ms/event | Speedup | Efficiency |
|---------|----------|---------|------------|
| 1       | 49.5     | 1.00x   | 100%       |
| 2       | 33.0     | 1.50x   | 75%        |
| 4       | 19.4     | 2.56x   | 64%        |
| 8       | 15.5     | 3.19x   | 40%        |
| 16      | 14.0     | 3.53x   | 22%        |

The native build (`-DCNCHASH_NATIVE_MARCH=ON`, AVX2) is roughly 2x
faster: ~23 ms/event at 1 thread and ~4 ms/event at 16 threads.

## Config comparison (portable build, ms/event)

| Config         | 1 thread | 4 threads |
|----------------|----------|-----------|
| small (12 st)  | 39       | 25        |
| medium (30 st) | 49       | 20        |
| large (80 st)  | 119      | 41        |
| dense (dang=2) | 49       | 19        |
| nmc=100        | 156      | 63        |

Weakly constrained events (few stations, multiple solutions) are
dominated by the MECH_PROB clustering, which now matches the original
HASH working-set behavior; the parallel path still helps there.

## Comparison with the original HASH v1.2 (same machine)

| Implementation                | ms/event |
|-------------------------------|----------|
| Original HASH v1.2 (1 thread) | ~145     |
| CNCHASH portable (1 thread)   | ~49      |
| CNCHASH portable (4 threads)  | ~20      |
| CNCHASH native AVX2 (4 thr.)  | ~8       |

## Batch mode (medium config, portable build)

| Threads | Batch ev/s | Scalar ev/s | Batch speedup |
|---------|------------|-------------|---------------|
| 1       | 18.6       | 18.7        | 1.00x         |
| 4       | 56.9       | 42.6        | 1.34x         |

Numbers are machine-specific. Always re-run the benchmark on the
target hardware before publishing performance claims.
