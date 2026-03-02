# ITM Rust C-API Usage Guide

This document outlines how to integrate the `itm` Rust propagation model into any C or C++ project, with practical examples for Unreal Engine 5 game development. The C-API is framework-agnostic and works with any C/C++ build system.

## Building the Library

Because this library is set up to expose a C-API, compiling it via cargo will automatically generate both dynamic and static libraries.
For Unreal Engine, it is recommended to use the **Static Library** (`libitm.a` on Linux/Mac, or `itm.lib` on Windows) to avoid complex DLL copying and loading scripts.

Run the release build process:
```bash
cargo build --release
```
The resulting libraries will be located in the `target/release/` directory.

---

## Linking the Static Library

### Windows (MSVC)
When linking the static library (`itm.lib`) on Windows with MSVC, you must also link against the following system libraries required by the Rust standard library and threading runtime (used by `rayon`):

- `ws2_32.lib` - Windows Sockets
- `advapi32.lib` - Advanced Windows APIs
- `userenv.lib` - User environment
- `ntdll.lib` - NT Layer DLL
- `bcrypt.lib` - Cryptographic primitives

If you include `itm.h`, these are automatically linked via `#pragma comment(lib, ...)` directives. If manually configuring your build system, ensure these libraries are specified.

### Unreal Engine 5
In your `.Build.cs` file:
```csharp
PublicAdditionalLibraries.Add(Path.Combine(ModulePath, "ThirdParty", "itm", "itm.lib"));
```

The system libraries are automatically included by UE5 on Windows; no additional configuration needed.

### Linux / macOS
On Unix-like systems, the Rust standard library typically links against:
- `pthread` (POSIX threads)
- `dl` (dynamic linking)
- `m` (math library)

Example linker flags:
```bash
gcc my_app.c -L./target/release -litm -lpthread -ldl -lm -o my_app
```

---

## Exposed Functions

The C-API exposes three distinct propagation calculations depending on your fidelity requirements.

### 1. Area TLS (`itm_area_tls_c`)
Calculates signal attenuation over a generic area defined by a statistical terrain roughness (`delta_h`). 
Use this when you don't have exact terrain geography between the transmitter and receiver but want a fast, statistically accurate guess.

```c
extern "C" {
    int32_t itm_area_tls_c(
        double h_tx_m,       // Transmitter Height (meters)
        double h_rx_m,       // Receiver Height (meters)
        int32_t tx_site_idx, // Tx Siting Criteria (0 = Random, 1 = Careful, 2 = Very Careful)
        int32_t rx_site_idx, // Rx Siting Criteria (0 = Random, 1 = Careful, 2 = Very Careful)
        double d_km,         // Distance between points (kilometers)
        double delta_h_m,    // Terrain Roughness estimate (meters)
        int32_t climate_idx, // Radio Climate (1 = Equatorial, 5 = Continental Temperate, etc.)
        double n_0,          // Surface Refractivity (e.g., 301.0)
        double f_mhz,        // Frequency (MHz)
        int32_t pol_idx,     // Polarization (0 = Horizontal, 1 = Vertical)
        double epsilon,      // Earth Dielectric Constant (e.g., 15.0)
        double sigma,        // Earth Conductivity (e.g., 0.005)
        int32_t mdvar,       // Variability Mode (0 = Broadcast, 1 = Point-to-Point, 2 = Mobile)
        double time_pct,     // Time Reliability % (e.g., 50.0)
        double loc_pct,      // Location Reliability % (e.g., 50.0)
        double sit_pct,      // Situation Confidence % (e.g., 50.0)
        double* out_a_db,    // OUTPUT: Attenuation Loss in dB
        uint32_t* out_warn   // OUTPUT: Warning flags bitmask
    );
}
```

### 2. Point-to-Point TLS (`itm_p2p_tls_c`)
Calculates highly-accurate attenuation mapped against *actual* terrain geography. 
In Unreal Engine, use your terrain geometry (or Cesium tiles) to sample heights via Line Traces between the source and listener, pack them into a continuous array, and pass the pointer.

```c
extern "C" {
    int32_t itm_p2p_tls_c(
        double h_tx_m,          // Transmitter Height (meters)
        double h_rx_m,          // Receiver Height (meters)
        const double* pfl_data, // Pointer to terrain array
        size_t pfl_len,         // Total elements in the array
        int32_t climate_idx,    // Radio Climate (e.g., 5)
        double n_0,             // Surface Refractivity (e.g., 301.0)
        double f_mhz,           // Frequency (MHz)
        int32_t pol_idx,        // Polarization (0 = Horizontal, 1 = Vertical)
        double epsilon,         // Earth Dielectric Constant (e.g., 15.0)
        double sigma,           // Earth Conductivity (e.g., 0.005)
        int32_t mdvar,          // Variability Mode 
        double time_pct,        // Time Reliability % 
        double loc_pct,         // Location Reliability % 
        double sit_pct,         // Situation Confidence % 
        double* out_a_db,       // OUTPUT: Attenuation Loss in dB
        uint32_t* out_warn      // OUTPUT: Warning flags bitmask
    );
}
```

#### Structuring the `pfl_data` Array in UE:
Instead of sending distance, `pfl_data` expects a structured double array (`TArray<double>`):
- **Index `0`**: The number of line segments/intervals (`Total Points - 1`).
- **Index `1`**: The distance between each point in meters (`Spacing`).
- **Index `2` through `N+1`**: Expected terrain Elevation at that exact segment.

### 3. Signal Radius Bubble (`itm_calculate_signal_radius_c`)
Designed specifically for Electromagnetic Warfare (EW) simulations.
Performs a binary search internally over distances to find the exact radius at which the Signal Power drops below the receiver's Sensitivity Threshold.

```c
extern "C" {
    int32_t itm_calculate_signal_radius_c(
        double f_mhz,         // Frequency (MHz)
        double power_w,       // Transmission Power (Watts)
        double h_tx_m,        // Transmitter Height (meters)
        double h_rx_m,        // Receiver Height (meters)
        double rx_sens_dbm,   // Receiver Sensitivity limit (e.g., -90.0 dBm)
        double delta_h_m,     // Terrain Roughness (meters)
        int32_t climate_idx,  // Radio Climate (e.g., 5)
        double n_0,           // Surface Refractivity (e.g., 301.0)
        int32_t pol_idx,      // Polarization (0/1)
        double epsilon,       // Earth Dielectric (e.g., 15.0)
        double sigma,         // Earth Conductivity (e.g., 0.005)
        double* out_radius_m  // OUTPUT: Max Signal Radius in Meters
    );
}
```
*Note: The output is in meters. Multiply by 100 before using it as a sphere radius component in Unreal Engine (which operates in centimeters).*

## Error Handling
All C-API functions return an `int32_t`. 
A return of `0` means **Success**.
A return of `1` means **Success with Warnings** (Read the `out_warn` value).
Any return `> 1` means **Failure**. Map these numeric values back to the `ItmError` enum in `src/math/itm.rs` to debug input validation issues.
