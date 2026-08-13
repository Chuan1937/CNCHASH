# Native Fortran backend (HASH-HP)

CNCHASH is built around a Modern Fortran + OpenMP backend. All numerical
work (grid search, S/P amplitudes, uncertainty, velocity tables) runs in
`libhashhp`; the Python layer provides the API, file handling, and batch
dispatch. Without the compiled library the package raises a clear build
error.

## Architecture

```
HASH_complete/            # original Fortran HASH v1.2 (immutable upstream)
cnchash/                  # Python package
  ├── backend/
  │   ├── base.py         # HashBackend contract
  │   └── fortran_backend.py  # ctypes binding to libhashhp
  ├── driver.py           # run_hash / run_hash_batch
  ├── io.py               # phase/station file handling
  └── utils.py            # pure-Python helpers
src/hash_hp/              # Modern Fortran core
  ├── hash_kinds.f90      # real64/int32 kinds
  ├── hash_geometry.f90   # TO_CAR/FPCOOR/cross
  ├── hash_rotation.f90   # cached rotation grid
  ├── hash_runtime.f90    # grid cache + OpenMP thread control
  ├── hash_focalmc.f90    # FOCALMC, OpenMP over rotations
  ├── hash_misfit.f90     # GET_GAP/GET_MISF
  ├── hash_uncertainty.f90# MECH_ROT/MECH_AVG/MECH_PROB
  ├── hash_amplitude.f90  # FOCALAMP_MC/GET_MISF_AMP
  ├── hash_velocity.f90   # LAYERTRACE/MK_TABLE/GET_TTS
  ├── hash_batch.f90      # full event pipeline + event-parallel batch
  └── hash_c_api.f90      # stable ISO_C_BINDING ABI (ctypes, no f2py)
```

## Design rules

- **No I/O in the kernel.** `src/hash_hp` receives arrays and returns
  arrays; Python handles files, formats, and exceptions.
- **No allocations in hot loops.** Per-event workspaces are reused;
  the rotation grid is built once per `dang` and shared read-only.
- **No nested OpenMP.** Batch mode parallelizes over events *or* over
  rotations, never both.
- **Deterministic by default.** Acceptable-mechanism selection takes
  the first `maxout` in grid order (`selection=0`). Pass
  `selection=1` to reproduce the original HASH random selection.
- **No `-ffast-math` by default.** Scientific correctness first; set
  `CNCHASH_FAST_MATH=ON` explicitly.
- **`HASH_complete/` is untouched.** It is the golden reference.

## Building

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

Portable binaries: add `-DCNCHASH_NATIVE_MARCH=OFF`. Without a Fortran
compiler, CNCHASH installs and runs on the Numba backend.

## Backend selection

```python
from cnchash import run_hash, available_backends, get_backend_info

run_hash(az, the, pol, qual, backend="auto")     # default (fortran)
run_hash(az, the, pol, qual, backend="fortran")  # explicit
print(available_backends())
print(get_backend_info())
```

Environment variables: `CNCHASH_BACKEND`, `CNCHASH_NUM_THREADS`,
`CNCHASH_HASHHP_LIB` (explicit library path).

## Batch mode

`run_hash_batch` sends many events to the Fortran core in one call
(CSR-style flat arrays), with event-level OpenMP parallelism for large
catalogs:

```python
from cnchash import run_hash_batch

results = run_hash_batch(events, nmc=30, backend="fortran", num_threads=16)
```

## Deliberate differences from HASH v1.2

| Item | HASH v1.2 | HASH-HP |
|------|-----------|---------|
| GET_GAP gaps | INTEGER (implicit typing truncates) | REAL, identical to Numba |
| GET_MISF fault vectors | radians misread as degrees (unit bug) | fixed: degrees passed to FPCOOR |
| GET_MISF_AMP stdr | weights summed for amplitude stations | identical |
| selection when nf > maxout | `rand(0)` random | deterministic default |
| velocity GET_TTS depth check | `d(nd0)` (uninitialized) | `d(ndep)` |
| velocity ray buffers | fixed 10001 (can overflow) | 20001 heap-allocated |

All deviations are bug fixes that make the backend internally
consistent; parity tests pin the Fortran and Numba backends to the same
behavior.
