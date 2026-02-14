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

### Input File Format

HASH input file (`example.inp`) specifies files and parameters:

```
scsn.pol           # Polarity reversal file
north1.phase       # Phase data file
output.out         # Output file 1 (mechanisms)
out2.out           # Output file 2 (acceptable planes)
8                  # npolmin: minimum polarities
90                 # max_agap: max azimuth gap (degrees)
60                 # max_pgap: max plunge gap (degrees)
5                  # dang: grid angle increment
30                 # nmc: Monte Carlo trials
500                # maxout: max output mechanisms
0.1                # badfrac: allowed bad fraction
120                # delmax: max distance (km)
45                 # cangle: clustering angle
0.1                # prob_max: probability threshold
```

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
