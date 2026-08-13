# Benchmarks

## Speed Comparison (CNCHASH vs original HASH)

![Speed Comparison](images/speed_comparison.png)

| Metric | CNCHASH | HASH v1.2 (Fortran) | Speedup |
|--------|---------|---------------------|---------|
| 24 events | 0.068s | 0.473s | **6.9x** |
| Per event | 2.85ms | 19.7ms | **6.9x** |
| 5000 events | 13.0s | 98.5s | **7.6x** |
| 10000 events | 26.0s | 197.0s | **7.6x** |

## Accuracy

![Accuracy Verification](images/accuracy_verification.png)

| Metric | CNCHASH | HASH v1.2 |
|--------|---------|-----------|
| Dip error median | 5.6° | 11.1° |
| Rake error median | 24.1° | 26.6° |

**Note:** Strike differences (40-80°) are normal - focal mechanisms have
two orthogonal nodal planes that both satisfy polarity data.

Evaluation status (latest notebook run):
- CNCHASH synthetic trials: 300/300 successful solutions
- HASH v1.2 synthetic trials: 60/60 successful solutions (direct-runs)

## Native backend measurements

Run locally with:

```bash
python benchmarks/benchmark_backends.py --events 1000 --seconds 2
```

Measured on a development machine (x86_64, gfortran, portable SSE2
build, no `-march=native`). Values are events/second for 30-station
events with 30 MC trials (`dang=5`).

## Fortran backend thread scaling (native build)

| Threads | ms/event | Speedup | Efficiency |
|---------|----------|---------|------------|
| 1       | 18.6     | 1.00x   | 100%       |
| 2       | 10.3     | 1.81x   | 90%        |
| 4       | 6.1      | 3.06x   | 77%        |
| 8       | 5.6      | 3.30x   | 41%        |
| 16      | 4.0      | 4.66x   | 29%        |

The native build (`-DCNCHASH_NATIVE_MARCH=ON`, AVX2) is roughly 2.5x
faster than the portable build at one thread; the portable build at
4 threads matches the native build at 1 thread.

## Config comparison (portable build, ms/event)

| Config         | 1 thread | 4 threads |
|----------------|----------|-----------|
| small (12 st)  | 81       | 52        |
| medium (30 st) | 53       | 18        |
| large (80 st)  | 127      | 36        |
| dense (dang=2) | 52       | 20        |
| nmc=100        | 166      | 54        |

Weakly constrained events (few stations, multiple solutions) are
dominated by the MECH_PROB clustering; the parallel path helps ~1.5x
there, while well-constrained events scale ~3x at 4 threads.

## Batch mode (medium config, portable build)

| Threads | Batch ev/s | Scalar ev/s | Batch speedup |
|---------|------------|-------------|---------------|
| 1       | 15.5       | 15.4        | 1.00x         |
| 4       | 49.7       | 36.2        | 1.37x         |

Numbers are machine-specific. Always re-run the benchmark on the
target hardware before publishing performance claims.
