# Installation

## Requirements

- Python >= 3.10
- NumPy >= 1.20.0
- A Fortran compiler with OpenMP support (gfortran >= 9) to build the
  native backend

## Install

```bash
pip install cnchash
```

After installing, build the native backend (or let the package build
helper do it):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

The package auto-detects `libhashhp` (also via `CNCHASH_HASHHP_LIB`).

## Development Setup

```bash
pip install -e ".[dev]"
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

## Verify Installation

```python
from cnchash import run_hash, get_backend_info
print(get_backend_info())
```

`get_backend_info()` reports the backend, OpenMP thread count, and
library path. Without the compiled library a clear build error is
raised.
