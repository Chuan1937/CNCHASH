# Benchmarks

## Speed Comparison (Numba vs original Fortran)

![Speed Comparison](images/speed_comparison.png)

| Metric | Python+Numba | Fortran | Speedup |
|--------|-------------|---------|---------|
| 24 events | 0.068s | 0.473s | **6.9x** |
| Per event | 2.85ms | 19.7ms | **6.9x** |
| 5000 events | 13.0s | 98.5s | **7.6x** |
| 10000 events | 26.0s | 197.0s | **7.6x** |

**Optimization:** Numba JIT compilation + NumPy vectorization

## Accuracy

![Accuracy Verification](images/accuracy_verification.png)

| Metric | Python | Fortran |
|--------|--------|---------|
| Dip error median | 5.6° | 11.1° |
| Rake error median | 24.1° | 26.6° |

**Note:** Strike differences (40-80°) are normal - focal mechanisms have
two orthogonal nodal planes that both satisfy polarity data.

Evaluation status (latest notebook run):
- Python synthetic trials: 300/300 successful solutions
- Fortran synthetic trials: 60/60 successful solutions (direct-runs)

## Native backend measurements

Run locally with:

```bash
python benchmarks/benchmark_backends.py --events 1000 --seconds 2
```

Measured on a development machine (x86_64, gfortran, portable build,
no `-march=native`). Values are events/second for 30-station events
with 30 MC trials (`dang=5`).

## Fortran backend thread scaling

| Threads | ev/s | Speedup | Efficiency |
|---------|------|---------|------------|
| 1       | 13.6 | 1.00x   | 100%       |
| 2       | 24.7 | 1.82x   | 91%        |
| 4       | 41.1 | 3.03x   | 76%        |
| 8       | 55.6 | 4.10x   | 51%        |
| 16      | 82.8 | 6.10x   | 38%        |

## Backend comparison (1 thread)

| Config        | Numba (ms/ev) | Fortran (ms/ev) |
|---------------|---------------|-----------------|
| small (12 st) | 306           | 82              |
| medium (30 st)| 36            | 73              |
| large (80 st) | 44            | 205             |
| dense (dang=2)| 36            | 73              |
| nmc=100       | 80            | 232             |

## Batch mode (medium config)

| Backend | Threads | Batch ev/s | Scalar ev/s | Batch speedup |
|---------|---------|------------|-------------|---------------|
| Fortran | 1       | 12.8       | 11.9        | 1.08x         |
| Fortran | 4       | 49.3       | 28.4        | 1.74x         |
| Numba   | 1       | 12.0       | 12.0        | 1.00x         |

Numbers are machine-specific. Always re-run the benchmark on the
target hardware before publishing performance claims.
