# ITM Longley-Rice: Irregular Terrain Model for Rust

[![CI](https://github.com/dgsantana/itm-rs/workflows/CI/badge.svg)](https://github.com/dgsantana/itm-rs/actions/workflows/ci.yml)
[![NTIA Validation](https://github.com/dgsantana/itm-rs/workflows/NTIA%20Validation/badge.svg)](https://github.com/dgsantana/itm-rs/actions/workflows/validation.yml)
[![Security Audit](https://github.com/dgsantana/itm-rs/workflows/Security%20Audit/badge.svg)](https://github.com/dgsantana/itm-rs/actions/workflows/security.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Rust Version](https://img.shields.io/badge/rust-1.85%2B-orange.svg)](https://www.rust-lang.org/)
[![Docs](https://img.shields.io/badge/docs-rustdoc-blue.svg)](https://docs.rs/itm-longley-rice)
[![Crates.io](https://img.shields.io/crates/v/itm-longley-rice.svg)](https://crates.io/crates/itm-longley-rice)

A high-performance, scientifically-validated Rust implementation of the **ITS Irregular Terrain Model (ITM)**, also known as the **Longley-Rice** model.

## Overview

`itm-longley-rice` provides a modern, memory-safe port of the standard atmospheric propagation model used for frequencies between 20 MHz and 20 GHz. This implementation is based directly on the **original NTIA C++ source code (v1.3)**, ensuring functional parity with the established scientific reference while providing the safety and performance benefits of the Rust ecosystem.

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

**FFI feature note:** The C-API is enabled by default. Rust-only consumers can disable it with `--no-default-features` or by setting `default-features = false` in their dependency entry.

---
*For technical inquiries regarding the underlying ITM model, refer to the [official NTIA/ITS documentation](https://github.com/NTIA/itm).*

## Continuous Integration

This project uses GitHub Actions for continuous integration and validation:

- **CI Workflow**: Tests on Windows, Linux, and macOS with stable and beta Rust toolchains
- **NTIA Validation**: Automated testing against NTIA reference implementation outputs  
- **Security Audit**: Daily dependency vulnerability scanning with cargo-audit and cargo-deny
- **Code Coverage**: Automated coverage reports via Codecov
- **Release Automation**: Automatic cross-platform builds and crates.io publishing on version tags

All tests must pass before merging, ensuring scientific accuracy and code quality.

## Contributing

Contributions are welcome! Please read [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

**Before submitting a PR:**
1. Run `cargo test --all` to ensure all tests pass
2. Run `cargo test --test ntia_validation` to validate against NTIA reference
3. Run `cargo fmt` and `cargo clippy` to maintain code quality
4. Add tests for new functionality

## License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

### Original ITM License

The original ITS Irregular Terrain Model developed by NTIA is in the **public domain** as a work of the U.S. Federal Government (Title 15 USC § 105). This Rust implementation is a derivative work that maintains compatibility with the original while adding modern safety guarantees and performance optimizations.

### Attribution

When using this library, please acknowledge:
- **Original ITM**: National Telecommunications and Information Administration (NTIA)
- **Original Authors**: Anita Longley and Phil Rice
- **This Implementation**: Daniel Santana

### References

- Original NTIA ITM Repository: https://github.com/NTIA/itm
- NTIA ITM License: https://github.com/NTIA/itm/blob/master/LICENSE.md
