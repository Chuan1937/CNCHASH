# Benchmarks

## Speed Comparison

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

**Note:** Strike differences (40-80°) are normal - focal mechanisms have two orthogonal nodal planes that both satisfy polarity data.

Evaluation status (latest notebook run):
- Python synthetic trials: 300/300 successful solutions
- Fortran synthetic trials: 60/60 successful solutions (direct-runs)
