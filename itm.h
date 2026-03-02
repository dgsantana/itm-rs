// Copyright (c) 2026 Daniel Santana
// SPDX-License-Identifier: MIT
//
// ITM-RS C API Header
// Rust implementation of the ITS Irregular Terrain Model (ITM)
//
// This is a derivative work of the original ITM developed by NTIA,
// which is in the public domain as a work of the U.S. Federal Government.
//
// See LICENSE file for details.

#ifndef ITM_H
#define ITM_H

#include <stdint.h>

// Windows MSVC: Link against required system libraries for Rust std + threading (rayon)
// Note: When using the static library, these are additional dependencies beyond
// what a typical C++ application already links (kernel32, msvcrt, etc.)
#if defined(_MSC_VER)
    #pragma comment(lib, "ws2_32.lib")      // Windows Sockets (required by std::net)
    #pragma comment(lib, "advapi32.lib")    // Advanced Windows APIs (required by std::thread)
    #pragma comment(lib, "userenv.lib")     // User environment (required by std::env)
    #pragma comment(lib, "bcrypt.lib")      // Cryptographic primitives (required by std::collections::HashMap)
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**
 * C-API wrapper for the ITM Area TLS function.
 * All pointers must be valid. Warnings are written to the `warnings` pointer.
 *
 * Returns 0 on success, or an error code > 0 on failure.
 */
int32_t itm_area_tls_c(
    double h_tx_m,
    double h_rx_m,
    int32_t tx_site_criteria_idx,
    int32_t rx_site_criteria_idx,
    double d_km,
    double delta_h_m,
    int32_t climate_idx,
    double n_0,
    double f_mhz,
    int32_t pol_idx,
    double epsilon,
    double sigma,
    int32_t mdvar,
    double time_pct,
    double location_pct,
    double situation_pct,
    double *out_a_db,
    uint32_t *out_warnings
);

/**
 * Calculates the maximum signal radius in meters for a given transmission.
 *
 * This performs a parallel multi-point search over distance using the ITM model
 * to find the point where the attenuation causes the signal to drop below the
 * receiver's sensitivity threshold.
 *
 * @param f_mhz Frequency in MHz.
 * @param power_w Transmission power in Watts.
 * @param h_tx_m Transmitter antenna height in meters.
 * @param h_rx_m Receiver antenna height in meters.
 * @param rx_sens_dbm Receiver sensitivity in dBm (e.g., -90.0). Defaults to -90.0 if 0.0 is passed.
 * @param delta_h_m Terrain irregularity (meters). Sensible default: 90.0.
 * @param climate_idx Climate type. Sensible default: 5 (Continental Temperate).
 * @param n_0 Surface refractivity. Sensible default: 301.0.
 * @param pol_idx Polarization (0:Horiz, 1:Vert). Sensible default: 1 (Vertical).
 * @param epsilon Earth dielectric constant. Sensible default: 15.0.
 * @param sigma Earth conductivity. Sensible default: 0.005.
 * @param out_radius_m Output pointer for radius in meters.
 *
 * Returns 0 on success, or an error code.
 */
int32_t itm_calculate_signal_radius_c(
    double f_mhz,
    double power_w,
    double h_tx_m,
    double h_rx_m,
    double rx_sens_dbm,
    double delta_h_m,
    int32_t climate_idx,
    double n_0,
    int32_t pol_idx,
    double epsilon,
    double sigma,
    double *out_radius_m
);

/**
 * C-API wrapper for the ITM Point-to-Point TLS function.
 * `pfl_data` must be a valid pointer to an array of `pfl_len` doubles.
 * The array format must match ITM expectations: [N-1, spacing_m, elev_1, elev_2, ... elev_N]
 *
 * Returns 0 on success, or an error code > 0 on failure.
 */
int32_t itm_p2p_tls_c(
    double h_tx_m,
    double h_rx_m,
    const double *pfl_data,
    uintptr_t pfl_len,
    int32_t climate_idx,
    double n_0,
    double f_mhz,
    int32_t pol_idx,
    double epsilon,
    double sigma,
    int32_t mdvar,
    double time_pct,
    double location_pct,
    double situation_pct,
    double *out_a_db,
    uint32_t *out_warnings
);

#ifdef __cplusplus
}
#endif

#endif // ITM_H
