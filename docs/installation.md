# Installation

## Requirements

- Python >= 3.10
- NumPy >= 1.20.0
- Numba >= 0.56.0
- SciPy >= 1.7.0

## Install

```bash
pip install cnchash
```

## Install from Source

```bash
pip install -r requirements.txt
pip install .
```

## Development Setup

```bash
pip install -r requirements.txt
pip install -e ".[dev]"
```

## Verify Installation

```python
from nchash import run_hash
print("CNCHASH installed successfully!")
```
