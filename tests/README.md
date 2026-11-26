# Validation against NTIA's ITM

Two suites check the port against NTIA's C++ Irregular Terrain Model.

## `ntia_reference.rs`: against NTIA's current C++

`reference/ntia-itm` is NTIA's C++ (https://github.com/NTIA/itm, under NTIA's terms), compiled by `cc` as a dev-dependency
only. The tests call it beside the port and require, per case, both to succeed with losses within 0.1 dB and equal
warning masks, or both to refuse the inputs:

- `published_point_to_point_cases_match_ntia`, `published_area_cases_match_ntia`: NTIA's published examples
  (`fixtures/ntia/p2p.csv`, `pfls.csv`, `area.csv`). The C++ reproduces their outputs within 0.05 dB (NTIA rounds
  them), and the port matches the C++.
- `generated_point_to_point_cases_match_ntia`, `generated_area_cases_match_ntia`: 5,000 deterministic cases per mode
  across ITM's accepted ranges, log-uniform where a range spans decades so every decade is as well represented
  (frequencies from 20 MHz to 20 GHz, antenna heights from 0.5 to 3,000 m, profile spacing from 1 m to 1 km, conductivity),
  with profiles of 2 to 1,000 intervals, every climate, both polarizations, refractivity from 250 to 400 N-units, every
  variability mode with and without its +10 and +20 flags, and time, location and situation from 1 to 99 %. Besides the
  loss and the warnings, each case compares the mode of propagation, and two refusals must carry the same error code.
  A sweep of 50,000 cases per mode agreed when this was written.

These tests caught the port's earlier errors of up to 36 dB in point-to-point and 27 dB in area mode (unfinished
profile analysis, an off-by-one in the receiver's ground height and in the terrain irregularity's sampling, a
fractional centre truncated in the least-squares fit, the area mode run as point-to-point in diffraction, the
crossover distance in troposcatter, NTIA's `DIM` ported as an absolute difference, a NaN reference attenuation passed on
where NTIA's `MAX` turns it to 0, the variability modes numbered in reverse, and warning bits and error codes that
differed from NTIA's).

## `ntia_validation.rs`: NTIA's command-line examples

The cases in NTIA's `cmd_examples` (`fixtures/pfl.txt` for the terrain). Area mode, TLS and CR, is checked against the
outputs NTIA printed, at ±0.2 dB. Point-to-point CR is checked against NTIA's current C++, loss within 0.01 dB and every
intermediate value: `o_p2pcr.txt` and `o_p2pcr_tbl.txt` were printed by NTIA's v1.4 DLL (January 2024), whose
intermediate values and some losses (up to 0.22 dB: 129.2 printed for 129.418 dB) differ from the current source.

## Running

```bash
cargo test --release --test ntia_reference
cargo test --release --test ntia_validation
```

The reference crate needs a C++ compiler (MSVC, GCC or Clang, as `cc` finds it); its build compiles copies of NTIA's
sources with the Windows-only include paths and `__declspec` made portable, leaving the vendored files as NTIA ships them.
