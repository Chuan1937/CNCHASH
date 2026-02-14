# CNCHASH

Python implementation of HASH for earthquake focal mechanism determination from P-wave polarities.

![Python](https://img.shields.io/badge/python-3.10+-orange.svg)
![License](https://img.shields.io/badge/license-BSD%203--blue.svg)
![Numba](https://img.shields.io/badge/numba-0.53+-red.svg)
![Numpy](https://img.shields.io/badge/numpy-1.19+-yellow.svg)
![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)

Python uses Numba JIT compilation optimization and vectorization, achieving speed improvements while maintaining complete consistency with the core Fortran algorithm.

![Speed Comparison](docs/images/speed_comparison.png)

## Performance

| Metric | Python+Numba | Fortran | Speedup |
|--------|-------------|---------|---------|
| 24 events | 0.068s | 0.473s | **6.9x** |
| Per event | 2.85ms | 19.7ms | **6.9x** |
| 1000 events | 2.6s | 19.7s | **7.5x** |

### For comprehensive benchmarks and analysis:

![Comprehensive Comparison](docs/images/comprehensive_comparison.png)

## Accuracy Verification

![Accuracy Verification](docs/images/accuracy_verification.png)

**Key Results:**
- Dip error median: < 10° 
- Rake error median: < 15° 
- Core algorithm matches Fortran exactly 

**Note:** Strike differences (40-80°) are normal - focal mechanisms have two orthogonal nodal planes that both satisfy polarity data.

## Quick Start

```bash
pip install cnchash
```

### Basic Usage (P-wave polarities only)

```python
from cnchash import run_hash
import numpy as np

# Azimuths, takeoff angles, polarities, quality
p_azi = np.array([45.0, 135.0, 225.0, 315.0])
p_the = np.array([30.0, 45.0, 60.0, 75.0])
p_pol = np.array([1, -1, 1, -1])  # 1=up, -1=down
p_qual = np.array([0, 0, 0, 0])

result = run_hash(p_azi, p_the, p_pol, p_qual)

print(f"Strike: {result['strike_avg']:.1f}")
print(f"Dip: {result['dip_avg']:.1f}")
print(f"Rake: {result['rake_avg']:.1f}")
print(f"Quality: {result['quality']}")
```

### With S/P Amplitude Ratio

```python
from cnchash import run_hash_with_amp
import numpy as np

# Same inputs as above, plus S/P amplitude ratios (log10 scale)
# sp_amp = 0.0 means no amplitude data for that station
sp_amp = np.array([0.3, -0.2, 0.5, 0.0])  # log10(S/P), 0.0 = no data

result = run_hash_with_amp(p_azi, p_the, p_pol, sp_amp)

print(f"Strike: {result['strike_avg']:.1f}")
print(f"Dip: {result['dip_avg']:.1f}")
print(f"Rake: {result['rake_avg']:.1f}")
print(f"Quality: {result['quality']}")
print(f"Polarity misfit: {result['mfrac']*100:.1f}%")
print(f"Amplitude misfit: {result['mavg']:.2f}")
```

## Features

- Grid search for focal mechanism determination
- Monte Carlo uncertainty analysis
- S/P amplitude ratio constraint
- Quality rating (A-D, E, F)
- Multiple phase file formats
- Core algorithm matches Fortran exactly

## Documentation

See [docs/README.md](docs/README.md) for full documentation including:
- API reference
- Algorithm details
- File format specifications
- Performance optimization

## Run Tests

```bash
jupyter notebook HASH_Tests.ipynb
```

## Project Structure

```
cnchash/
├── core.py        # Grid search algorithm (focalmc)
├── amp_subs.py    # S/P amplitude ratio (focalamp_mc)
├── uncertainty.py # Uncertainty analysis (mech_prob)
├── driver.py      # Main driver (run_hash, run_hash_with_amp)
├── io.py          # File I/O
└── utils.py       # Utilities

HASH_complete/     # Complete Fortran code with examples
```

## Author

He XingChen

## License

BSD 3-Clause

## References

Hardebeck, Jeanne L. and Peter M. Shearer, A new method for determining first-motion
focal mechanisms, Bulletin of the Seismological Society of America, 92,
2264-2276, 2002.

Hardebeck, Jeanne L. and Peter M. Shearer, Using S/P Amplitude Ratios to
Constrain the Focal Mechanisms of Small Earthquakes, Bulletin of the
Seismological Society of America, 93, 2434-2444, 2003.
