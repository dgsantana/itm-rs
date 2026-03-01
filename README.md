# ITM-RS: Irregular Terrain Model for Rust

A high-performance, scientifically-validated Rust implementation of the **ITS Irregular Terrain Model (ITM)**, also known as the **Longley-Rice** model.

## Overview

`itm-rs` provides a modern, memory-safe port of the standard atmospheric propagation model used for frequencies between 20 MHz and 20 GHz. This implementation is based directly on the **original NTIA C++ source code (v1.3)**, ensuring functional parity with the established scientific reference while providing the safety and performance benefits of the Rust ecosystem.

This library is designed for integration into large-scale simulations, including **Unreal Engine 5** and other C++ environments, via a robust [C-API](C_API_USAGE.md).

## Scientific Backing

The core algorithm predicts terrestrial radiowave propagation based on electromagnetic theory and extensive empirical data collected by Anita Longley and Phil Rice.

### Supported Propagation Mechanisms
- **Free Space Loss**: Standard LOS propagation.
- **Diffraction**: Multi-faceted implementation using **Vogler's smooth-earth approximation** and Bullington knife-edge methods.
- **Troposcatter**: Selection and calculation of tropospheric scatter for trans-horizon paths.
- **Variability Analysis**: Statistical modeling of signal reliability across time, location, and situation.

### Technical References
This implementation precisely follows the algorithms defined in the following primary references:
- **NTIA Technical Report TR-82-100**: *A Guide to the Use of the ITS Irregular Terrain Model in the Area Prediction Mode*.
- **NTIA Technical Report ERL 79-ITS 67**: *Prediction of Tropospheric Radio Transmission Loss Over Irregular Terrain: A Computer Method*.
- **Hufford, G. A. (1985)**: *Memorandum on the ITS Irregular Terrain Model, version 1.2.2*.

## Key Core Features

- **Point-to-Point (P2P) Mode**: Calculate path loss using actual terrain elevation profiles sampled from GIS data or game engine geometry.
- **Area Prediction Mode**: Generate statistical coverage estimates based on terrain roughness (`delta_h`) when specific profiles are unavailable.
- **Parallelized Signal Radius Calculation**: Efficiently determine the maximum distance a signal can travel before dropping below a specified sensitivity threshold. This is implemented as a multi-point parallel search using **Rayon** for high performance.
- **C-API Integration**: Built-in support for static linking into Unreal Engine, allowing real-time electromagnetic warfare (EW) and signal simulation.

## Heritage

This project is a pragmatic, high-fidelity port of the **NTIA/ITS C++ implementation**. It maintains the internal numeric heuristics and mathematical structures found in the canonical FORTRAN and C++ versions, ensuring that results match validated government benchmarks.

## Getting Started

Refer to the following documentation for integration details:
- [C-API Usage Guide](C_API_USAGE.md): Detailed instructions for Unreal Engine and C/C++ integration.
- [Error and Warning Reference](itm/ERRORS_AND_WARNINGS.md): Detailed breakdown of algorithm warning flags.

---
*For technical inquiries regarding the underlying ITM model, refer to the [official NTIA/ITS documentation](https://github.com/NTIA/itm).*
