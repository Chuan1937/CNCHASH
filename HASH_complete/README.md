# HASH v1.2 - Complete Distribution

HASH (Hardebeck and Shearer, 2002, 2003) is a Fortran package for computing earthquake focal mechanisms from first-motion polarity data, with optional S/P amplitude ratio constraints.

## Directory Structure

```
HASH_complete/
├── src/                    # Source code
│   ├── drivers/           # Main driver programs
│   │   ├── hash_driver1.f - Basic driver (polarity only)
│   │   ├── hash_driver2.f - Driver with velocity model
│   │   ├── hash_driver3.f - Driver with S/P amplitude ratios
│   │   ├── hash_driver4.f - Driver with 5-char station names
│   │   └── hash_driver5.f - Driver with simulated data
│   ├── subs/              # Subroutines
│   │   ├── fmamp_subs.f   - S/P amplitude ratio calculations (FOCALAMP_MC, GET_MISF_AMP)
│   │   ├── fmech_subs.f   - Core focal mechanism grid search (FOCALMC)
│   │   ├── pol_subs.f     - Polarity data processing (GET_MISF, GET_GAP)
│   │   ├── station_subs.f - Station information processing
│   │   ├── uncert_subs.f  - Uncertainty analysis (MECH_PROB)
│   │   ├── util_subs.f    - Utility functions
│   │   └── vel_subs.f     - Velocity model processing
│   └── include/           # Header files
│       ├── param.inc      - Parameter definitions
│       ├── rot.inc        - Rotation grid parameters
│       └── vel.inc        - Velocity model parameters
├── examples/              # Example data
│   ├── input/            # Input files
│   │   ├── example*.inp  - Example input files
│   │   ├── north*.phase  - Phase data files
│   │   ├── north3.amp    - S/P amplitude ratio data
│   │   └── north3.statcor - Station corrections
│   └── output/           - Expected output files
│       └── example*.out* - Example outputs
├── data/                  # Data files
│   ├── stations/         - Station information
│   │   ├── scsn.stations - SCSN station list
│   │   └── scsn.reverse  - Polarity reversal file
│   └── velocity/         - Velocity models
│       ├── vz.socal      - Southern California
│       ├── vz.north      - Northern California
│       └── vz.*          - Other velocity models
├── doc/                   # Documentation
│   └── HASH_manual_v1.2.pdf
├── Makefile              - Build configuration
└── README.md             # This file
```

## Driver Programs

1. **hash_driver1.f** - Basic driver using only P-wave polarity data
2. **hash_driver2.f** - Driver with 1D velocity model for takeoff angle calculation
3. **hash_driver3.f** - Driver using both P-wave polarity AND S/P amplitude ratios
4. **hash_driver4.f** - Driver with 5-character station name support
5. **hash_driver5.f** - Driver with simulated data for testing

## Key Subroutines

### FOCALMC (fmech_subs.f)
- Grid search using P-wave polarities only
- Input: azimuth, takeoff angle, polarity data
- Output: acceptable focal mechanisms

### FOCALAMP_MC (fmamp_subs.f)
- Grid search using P-wave polarities AND S/P amplitude ratios
- Input: azimuth, takeoff angle, polarity, S/P ratios
- Output: acceptable focal mechanisms
- Key formula: `sp_ratio = log10(4.9 * s_amp / p_amp)`

### GET_MISF_AMP (fmamp_subs.f)
- Calculate misfit for a given mechanism
- Returns: polarity misfit fraction, average S/P misfit, station distribution ratio

## Building

```bash
# Edit Makefile to set your Fortran compiler (gfortran, ifort, etc.)
make
```

## Running Examples

```bash
# Example 1: Polarity only
./hash_driver1 < examples/input/example1.inp

# Example 3: With S/P amplitude ratios
./hash_driver3 < examples/input/example3.inp
```

## S/P Amplitude Ratio Data Format

The amplitude file (e.g., north3.amp) contains:
```
event_id  num_stations
station_name  component  network  qns1  qns2  qpamp  qsamp
...
```

Where:
- qns1, qns2: S-wave amplitude quality
- qpamp: P-wave amplitude
- qsamp: S-wave amplitude

Station correction file (e.g., north3.statcor):
```
station_name  component  network  correction(log10)
```

## References

- Hardebeck, J. L., and P. M. Shearer (2002), A new method for determining first-motion focal mechanisms, Bull. Seism. Soc. Am., 92, 2264-2276.
- Hardebeck, J. L., and P. M. Shearer (2003), Using S/P amplitude ratios to constrain the focal mechanisms of small earthquakes, Bull. Seism. Soc. Am., 93, 2434-2444.

## License

See the HASH_manual_v1.2.pdf for license information.
