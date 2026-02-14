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
| `driver.py` | `run_hash`, `run_hash_with_amp`, `run_hash_from_file` | Main entry points |
| `core.py` | `focalmc`, `get_misfit`, `get_gap` | Grid search algorithm |
| `amp_subs.py` | `focalamp_mc`, `get_misf_amp` | S/P amplitude ratio |
| `uncertainty.py` | `mech_prob`, `mech_rot`, `mech_avg` | Uncertainty analysis |
| `io.py` | `read_phase_file`, `write_mechanism_output` | File I/O |
| `utils.py` | `fp_coord_angles_to_vectors`, `fp_coord_vectors_to_angles` | Coordinate conversions |
| `velocity.py` | `make_table`, `get_angle` | Velocity model tables |

---

## Low-Level Functions

### focalmc()

Core grid search algorithm.

```python
result = focalmc(p_azi_mc, p_the_mc, p_pol, p_qual, npol, nmc, dang, maxout, nextra, ntotal)
```

Returns all acceptable mechanisms with fault normals and slip vectors.

### get_misfit()

Calculate polarity misfit.

```python
mfrac, stdr = get_misfit(npol, p_azi, p_the, p_pol, p_qual, strike, dip, rake)
```

### get_gap()

Calculate azimuthal and plunge gaps.

```python
magap, mpgap = get_gap(npol, p_azi, p_the)
```

### mech_prob()

Calculate mechanism probability and average.

```python
result = mech_prob(nf, norm1, norm2, cangle=45.0, prob_max=0.1)
```
