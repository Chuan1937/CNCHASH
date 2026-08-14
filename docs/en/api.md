# API Reference

## Main Functions

### run_hash()

Determine focal mechanism from P-wave polarities.

```python
result = run_hash(p_azi, p_the, p_pol, p_qual, **kwargs)
```

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `p_azi` | ndarray | required | Azimuths (degrees), shape (nsta,) or (nsta, nmc) |
| `p_the` | ndarray | required | Takeoff angles (degrees) |
| `p_pol` | ndarray | required | Polarities: 1=up, -1=down |
| `p_qual` | ndarray | required | Quality: 0=impulsive, 1=emergent |
| `dang` | float | 5.0 | Grid angle increment |
| `nmc` | int | 30 | Monte Carlo trials |
| `maxout` | int | 500 | Max output mechanisms |
| `badfrac` | float | 0.1 | Allowed bad polarity fraction |
| `npolmin` | int | 8 | Minimum polarities |
| `max_agap` | float | 90.0 | Max azimuth gap |
| `max_pgap` | float | 60.0 | Max plunge gap |

**Returns:**

| Key | Type | Description |
|-----|------|-------------|
| `success` | bool | Solution found |
| `strike_avg` | float | Strike (degrees) |
| `dip_avg` | float | Dip (degrees) |
| `rake_avg` | float | Rake (degrees) |
| `quality` | str | A, B, C, D, E, or F |
| `mfrac` | float | Misfit fraction (0-1) |
| `prob` | float | Solution probability |
| `stdr` | float | Station distribution ratio |
| `nout2` | int | Number of solutions |

---

### run_hash_with_amp()

Determine mechanism using polarities + S/P amplitude ratios.

```python
result = run_hash_with_amp(p_azi, p_the, p_pol, sp_amp, **kwargs)
```

**Additional Parameter:**

| Name | Type | Description |
|------|------|-------------|
| `sp_amp` | ndarray | S/P ratios (log10), 0.0 = no data |

**Additional Returns:**

| Key | Type | Description |
|-----|------|-------------|
| `mavg` | float | Amplitude misfit (log10) |
| `npol` | int | Polarity count |
| `nspr` | int | S/P ratio count |

---

### run_hash_batch() / run_hash_batch_with_amp()

Process many events in one native call (event-level OpenMP).

```python
results = run_hash_batch(events, nmc=30, backend="fortran", num_threads=16)
```

### run_hash_from_file()

Process events from HASH input file.

```python
results = run_hash_from_file("example.inp")
```

---

## Quality Rating

| Grade | Criteria |
|-------|----------|
| **A** | prob > 0.8, var ≤ 25°, misfit ≤ 15%, stdr ≥ 0.5 |
| **B** | prob > 0.6, var ≤ 35°, misfit ≤ 20%, stdr ≥ 0.4 |
| **C** | prob > 0.5, var ≤ 45°, misfit ≤ 30%, stdr ≥ 0.3 |
| **D** | Solution found but below C criteria |
| **E** | Azimuth or plunge gap too large |
| **F** | No acceptable mechanism found |

---

## Modules

| Module | Functions | Description |
|--------|-----------|-------------|
| `driver.py` | `run_hash`, `run_hash_with_amp`, `run_hash_batch`, `run_hash_from_file` | Main entry points |
| `backend/fortran_backend.py` | `run_event`, `run_event_amp`, `run_batch`, `build_velocity_table`, `get_tts` | Native backend binding |
| `io.py` | `read_phase_file`, `write_mechanism_output` | File I/O |
| `utils.py` | `fp_coord_angles_to_vectors`, `fp_coord_vectors_to_angles` | Coordinate conversions |

---

## Low-Level Backend Methods

These live on the Fortran backend instance (via `backend.get_backend("fortran")`)
and are used by the test suite to pin results to the original HASH:

### run_event()

Full polarity-only pipeline for one event (equivalent to `run_hash` without
input normalization).

### run_event_amp()

Full polarity + S/P amplitude pipeline for one event.

### run_batch() / run_batch_amp()

Batch pipelines (CSR-style flat arrays) used by `run_hash_batch` and
`run_hash_batch_with_amp`.

### mech_prob()

Preferred mechanism(s) and probability from a set of acceptable mechanisms.

```python
result = backend.mech_prob(nf, faults, slips, cangle=45.0, prob_max=0.1)
```

### get_misfit() / get_misfit_amp() / get_gap()

Misfit fractions, station distribution ratio, and azimuthal/takeoff gaps for
a given mechanism.

### build_velocity_table() / get_tts()

Takeoff-angle table construction and interpolation for 1D velocity models.

### get_rotation_grid()

The cached rotation grid for a given `dang` (b1/b2/b3 direction cosines).
