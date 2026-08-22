# Usage

## Basic Usage

Determine focal mechanism from P-wave polarities:

```python
from cnchash import run_hash
import numpy as np

# Station data
p_azi = np.array([45, 135, 225, 315, 0, 90, 180, 270])    # Azimuths (degrees)
p_the = np.array([30, 45, 60, 75, 40, 50, 55, 65])         # Takeoff angles (degrees)
p_pol = np.array([1, -1, 1, -1, 1, -1, 1, -1])             # Polarities: 1=up, -1=down
p_qual = np.zeros(8)                                        # Quality: 0=impulsive, 1=emergent

# Run HASH
result = run_hash(p_azi, p_the, p_pol, p_qual)

# Output
print(f"Strike: {result['strike_avg']:.1f}°")
print(f"Dip: {result['dip_avg']:.1f}°")
print(f"Rake: {result['rake_avg']:.1f}°")
print(f"Quality: {result['quality']}")
```

## With S/P Amplitude

Add S/P amplitude ratios for better constraint:

```python
from cnchash import run_hash_with_amp

# S/P ratios in log10 scale
# 0.0 = no data, typical range: -1.0 to 2.0
sp_amp = np.array([0.3, -0.2, 0.5, 0.0, 0.4, -0.1, 0.6, 0.2])

result = run_hash_with_amp(p_azi, p_the, p_pol, sp_amp)

print(f"Polarity misfit: {result['mfrac']*100:.1f}%")
print(f"Amplitude misfit: {result['mavg']:.2f}")
```

## From Input File

Process multiple events from HASH input file:

```python
from cnchash import run_hash_from_file

results = run_hash_from_file("example.inp")

for result in results:
    if result['success']:
        print(f"S={result['strike_avg']:.0f}° D={result['dip_avg']:.0f}° "
              f"R={result['rake_avg']:.0f}° Q={result['quality']}")
```

Velocity model file names can be listed at the end of the input file
(`vz.*` lines, like the original HASH); `run_hash_from_file` loads
them and rotates them per Monte Carlo trial (hash_driver2.f logic).
Models can also be passed explicitly:

```python
from cnchash import run_hash_from_file

depth = [0.0, 5.0, 20.0, 35.0]
velocity = [5.5, 5.8, 6.2, 6.8]

results = run_hash_from_file("example.inp", velocity_models=[(depth, velocity)])
```

## Takeoff Angles from Velocity Models

`TakeoffTable` wraps the native takeoff-angle tables (the original
HASH MK_TABLE + GET_TTS) with scalar and vectorized lookups:

```python
from cnchash import TakeoffTable, compute_takeoff_angles
import numpy as np

depth = np.array([0.0, 5.0, 20.0, 35.0])   # km
velocity = np.array([5.5, 5.8, 6.2, 6.8])  # km/s

table = TakeoffTable(depth, velocity)
angle = table.takeoff(distance=30.0, source_depth=8.0)   # degrees
angles = table.takeoff_batch([10.0, 30.0, 60.0], [8.0, 12.0, 15.0])

# Convenience: one table build plus a batch lookup
angles = compute_takeoff_angles(depth, velocity, [10.0, 30.0], 8.0)
```

Convention (same as HASH): degrees from the vertical, 0 = straight up,
90 = horizontal, 180 = straight down. Lookups outside the table depth
range raise `TakeoffRangeError` (batch lookups return NaN). Table
ranges and steps are configurable through `params`, see
`DEFAULT_TABLE_PARAMS`.

### HASH Input File (Original Fortran Format)

The original Fortran HASH program uses a fixed-format input file (`example.inp`) with one value per line. No comments are supported.

**Format 1 Structure:**

| Line | Content | Description |
|------|---------|-------------|
| 1 | polfile | Polarity reversal file |
| 2 | phasefile | Phase data file |
| 3 | outfile1 | Output file (mechanisms) |
| 4 | outfile2 | Output file (acceptable planes) |
| 5-14 | parameters | Algorithm parameters (see below) |

**Algorithm Parameters (Lines 5-14):**

| Line | Parameter | Description |
|------|-----------|-------------|
| 5 | npolmin | Minimum polarities required |
| 6 | max_agap | Max azimuth gap (degrees) |
| 7 | max_pgap | Max plunge gap (degrees) |
| 8 | dang | Grid angle increment (degrees) |
| 9 | nmc | Monte Carlo trials |
| 10 | maxout | Max output mechanisms |
| 11 | badfrac | Allowed bad polarity fraction |
| 12 | delmax | Max distance (km) |
| 13 | cangle | Clustering angle (degrees) |
| 14 | prob_max | Probability threshold |

**Other Format Variants:**

| Format | File Structure | Handling |
|--------|----------------|----------|
| 1 | polfile, phasefile, out1, out2, params | File azimuths/takeoffs + per-station uncertainties (driver1) |
| 2/4 | stationfile, polfile, phasefile, out1, out2, params + vz lines | Velocity takeoff tables + depth perturbation (driver2) |
| 3 | stationfile, polfile, statcor, amp, phasefile, out1, params | Amplitude pipeline: SNR filtering + station corrections + FOCALAMP (driver3) |
| 5 | polfile, simulfile, phasefile, out1, out2, params | SIMULPS takeoff file + polarity matching (driver5) |

All pipelines mirror the original HASH drivers:

- **Amplitudes**: `read_amp_file` reads `ICUSP NIN` + station records;
  `read_statcor_file` provides GET_COR station corrections
  (`log10(S/P) - qcor`); `ratmin` (default 3.0) filters on P/S
  signal-to-noise, and only impulsive (I) picks are used.
- **SIMULPS takeoffs**: `read_simul_takeoff_file` matches polarities by
  event id and station name (first three characters), with per-station
  uncertainties (sazi/sthe, default 0.5).
- **File geometry**: format-1 phase files carry distance/azimuth/
  takeoff/uncertainties which are used directly for the Monte Carlo
  trials (driver1).

### Supported File Formats

#### Phase Files (`read_phase_file`)

Supports multiple formats, auto-detected:

| Format | Example | Description |
|--------|---------|-------------|
| 1 | north1.phase | 2-digit year, compressed |
| 2 | north2.phase | 4-digit year, separate stations |
| 3 | north4.phase | 8-digit date format |
| 4 | north5.simul | SIMUL2000 format |

**Format 1** (2-digit year):
```
YY MDD hhmmss lat lon depth mag ... event_id
STATION polcode ...
...
event_id
```

**Format 2** (4-digit year):
```
YYYY MDDHHMM lat lon depth mag ... event_id
STATION NETWORK COMPONENT ONSET POLARITY
```

#### Station Files (`read_station_file`)

```
STATION COMPONENT lat lon elevation
```

Example:
```
PAS   HHZ   34.1512  -118.1567   0.405
CAL   HHZ   34.1296  -117.9258   0.757
```

#### Polarity Reversal Files (`read_polarity_reversal_file`)

Records station polarity reversals over time:

```
STATION start_date end_date
```

Dates in YYYYMMDD format. `end_date = 0` means ongoing.

Example:
```
PAS  20100101  20121231
CAL  20150601  0
```

#### Velocity Model Files (`read_velocity_model`)

```
depth velocity
```

Example:
```
0.0   5.5
5.0   5.8
10.0  6.2
20.0  6.8
```

#### Input File Variants

Different HASH formats supported:

| Format | Files |
|--------|-------|
| 1 | polfile, phasefile, out1, out2, params |
| 2/4 | stationfile, polfile, phasefile, out1, out2, params |
| 3 | stationfile, polfile, statcor, amp, phasefile, out1, params |
| 5 | polfile, simulfile, phasefile, out1, out2, params |

### Reading Files Directly

```python
from cnchash.io import (
    read_phase_file,
    read_station_file,
    read_polarity_reversal_file,
    read_velocity_model,
    read_hash_input_file
)

# Read phase data
events = read_phase_file("north1.phase")

# Read station coordinates
stations = read_station_file("scsn.stations")

# Read polarity reversals
reversals = read_polarity_reversal_file("scsn.pol")

# Read velocity model
depths, velocities = read_velocity_model("vz.layer")

# Parse input file
params = read_hash_input_file("example.inp")
```

## Parameters

```python
result = run_hash(p_azi, p_the, p_pol, p_qual,
                  dang=5.0,      # Grid angle increment (degrees)
                  nmc=30,        # Monte Carlo trials
                  maxout=500,    # Max output mechanisms
                  badfrac=0.1,   # Allowed bad polarity fraction
                  npolmin=8,     # Minimum polarities required
                  max_agap=90.0, # Max azimuth gap (degrees)
                  max_pgap=60.0) # Max plunge gap (degrees)
```

## Tips

- Minimum 8 polarities recommended
- Azimuth gap < 90° for well-constrained solutions
- Use S/P ratios when available
- `nmc=30` is usually sufficient
- `dang=5°` balances speed and accuracy
