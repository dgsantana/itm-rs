# Contributing to ITM-RS

Thank you for your interest in contributing to ITM-RS! This document provides guidelines for contributing to this Rust implementation of the ITS Irregular Terrain Model.

## Project Goals

1. **Scientific Fidelity**: Maintain compatibility with the NTIA reference implementation
2. **Safety**: Leverage Rust's memory safety guarantees
3. **Performance**: Optimize for real-time applications (e.g., Unreal Engine integration)
4. **Usability**: Provide clean APIs for both Rust and C/C++ consumers

## Areas for Contribution

### High Priority
- **Validation Testing**: Add more test cases comparing against NTIA reference outputs
- **Documentation**: Improve code comments, API docs, and usage examples
- **Performance Optimization**: Profile and optimize hot paths
- **Bug Fixes**: Address issues in the GitHub issue tracker

### Medium Priority
- **Platform Support**: Test and improve compatibility across Windows/Linux/macOS
- **Terrain Profile Parsing**: Add utilities for common GIS formats (GeoTIFF, SRTM, etc.)
- **Examples**: Real-world usage examples and tutorials
- **Benchmarks**: Performance comparison with NTIA C++ implementation

### Research/Future Work
- **Algorithm Improvements**: Investigate discrepancies with older NTIA test files (see [tests/README.md](tests/README.md))
- **Extended Frequency Range**: Research applicability beyond 20 MHz - 20 GHz
- **GPU Acceleration**: Explore CUDA/OpenCL for massive parallel area calculations

## Development Setup

### Prerequisites
- Rust 1.85+ (with 2024 edition support)
- C compiler (MSVC on Windows, GCC/Clang on Linux/macOS) for FFI testing

### Building
```bash
# Clone the repository
git clone https://github.com/dgsantana/itm-rs
cd itm-rs

# Run tests
cargo test

# Build C static library
cargo build --release

# Run validation tests
cargo test --test ntia_validation
```

## Code Style

- Follow standard Rust conventions (`cargo fmt`)
- Run `cargo clippy` before submitting
- Maintain 100% safety: avoid `unsafe` unless absolutely necessary (currently only in FFI layer)
- Document all public APIs with rustdoc comments

## Testing Guidelines

### Unit Tests
- Add unit tests for new mathematical functions
- Maintain >80% code coverage for core modules

### Integration Tests
- Validate against NTIA reference outputs when adding new modes
- Document expected vs actual results with tolerances

### Validation Tests  
- See [tests/ntia_validation.rs](tests/ntia_validation.rs) for examples
- Use ±0.2 dB tolerance for path loss comparisons
- Note known issues with older NTIA reference files

## Submitting Changes

1. **Fork** the repository
2. **Create a branch**: `git checkout -b feature/your-feature-name`
3. **Make changes** with clear, atomic commits
4. **Add tests** for new functionality
5. **Update documentation** (README, code comments, etc.)
6. **Run tests**: `cargo test --all`
7. **Format code**: `cargo fmt`
8. **Check lints**: `cargo clippy`
9. **Submit PR** with clear description of changes

### Commit Messages
- Use present tense: "Add feature" not "Added feature"
- Reference issues: "Fix #123: Description"
- Keep first line under 50 characters
- Provide detailed context in commit body if needed

## Scientific Validation

When modifying core algorithms:
1. Document the mathematical basis with equation references
2. Compare outputs against NTIA C++ reference implementation
3. Note any intentional deviations with rationale
4. Update validation tests accordingly

## FFI/C-API Changes

When modifying the FFI layer ([src/ffi.rs](src/ffi.rs)):
1. Update [itm.h](itm.h) header file
2. Update [C_API_USAGE.md](C_API_USAGE.md) documentation
3. Test with actual C/C++ consumer (e.g., minimal test program)
4. Ensure Windows MSVC and Linux GCC compatibility

## Questions?

- **Algorithm Questions**: Refer to [official NTIA ITM](https://github.com/NTIA/itm) or NTIA technical reports
- **Implementation Questions**: Open a GitHub issue with the `question` label
- **Bug Reports**: Open a GitHub issue with detailed reproduction steps

## Code of Conduct

- Be respectful and constructive in all interactions
- Focus on technical merit and scientific accuracy
- Welcome newcomers and help them get started
- Acknowledge contributions and give credit appropriately

## License

By contributing to ITM-RS, you agree that your contributions will be licensed under the MIT License. All contributions must be your original work or properly attributed to their source.

---

Thank you for helping improve ITM-RS! 🦀📡
