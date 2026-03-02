// Copyright (c) 2026 Daniel Santana
// SPDX-License-Identifier: MIT
//
// FFI (Foreign Function Interface) layer for C/C++ integration.
// This module provides C-compatible exports for use in Unreal Engine and other
// C/C++ applications.

use crate::math::itm::{itm_area_tls, itm_p2p_tls};
use rayon::prelude::*;
use std::slice;

/// C-API wrapper for the ITM Area TLS function.
/// All pointers must be valid. Warnings are written to the `warnings` pointer.
///
/// Returns 0 on success, or an error code > 0 on failure.
///
/// # Safety
/// - `out_a_db` must be a non-null, properly aligned pointer to writable `f64` storage.
/// - `out_warnings` must be a non-null, properly aligned pointer to writable `u32` storage.
/// - Both output pointers must remain valid for the duration of the call.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn itm_area_tls_c(
    h_tx_m: f64,
    h_rx_m: f64,
    tx_site_criteria_idx: i32,
    rx_site_criteria_idx: i32,
    d_km: f64,
    delta_h_m: f64,
    climate_idx: i32,
    n_0: f64,
    f_mhz: f64,
    pol_idx: i32,
    epsilon: f64,
    sigma: f64,
    mdvar: i32,
    time_pct: f64,
    location_pct: f64,
    situation_pct: f64,
    out_a_db: *mut f64,
    out_warnings: *mut u32,
) -> i32 {
    if out_a_db.is_null() || out_warnings.is_null() {
        return 99; // ErrorInvalidPointer
    }

    match itm_area_tls(
        h_tx_m,
        h_rx_m,
        tx_site_criteria_idx,
        rx_site_criteria_idx,
        d_km,
        delta_h_m,
        climate_idx,
        n_0,
        f_mhz,
        pol_idx,
        epsilon,
        sigma,
        mdvar,
        time_pct,
        location_pct,
        situation_pct,
    ) {
        Ok((a_db, warnings, _inter)) => {
            unsafe {
                *out_a_db = a_db;
                *out_warnings = warnings;
            }
            0 // Success
        }
        Err(err) => err.code(),
    }
}

/// Calculates the maximum signal radius in meters for a given transmission.
///
/// This performs a binary search over distance (from 0.001 km to 2000 km) using the ITM model
/// to find the point where the `a_db` (attenuation loss) causes the signal to drop below the
/// receiver's sensitivity threshold.
///
/// # Arguments
/// * `f_mhz` - Frequency in MHz.
/// * `power_w` - Transmission power in Watts.
/// * `h_tx_m` - Transmitter antenna height in meters.
/// * `h_rx_m` - Receiver antenna height in meters.
/// * `rx_sens_dbm` - The sensitivity threshold of the receiver in dBm (e.g., -90.0). Defaults to -90.0 if you pass 0.0.
/// * `delta_h_m` - Terrain irregularity parameter in meters. **Default: 90.0** (Average hilly terrain).
/// * `climate_idx` - Climate type. **Default: 5** (Continental Temperate).
/// * `n_0` - Surface refractivity. **Default: 301.0** (Average atmospheric conditions).
/// * `pol_idx` - Polarization (0 for Horizontal, 1 for Vertical). **Default: 1** (Vertical).
/// * `epsilon` - Relative earth dielectric constant. **Default: 15.0** (Average/Standard ground).
/// * `sigma` - Earth surface conductivity. **Default: 0.005** (Average/Standard ground).
/// * `out_radius_m` - Output pointer where the resulting max radius in meters will be written.
///
/// Returns 0 on success, or an error code.
///
/// # Safety
/// - `out_radius_m` must be a non-null, properly aligned pointer to writable `f64` storage.
/// - `out_radius_m` must remain valid for the duration of the call.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn itm_calculate_signal_radius_c(
    f_mhz: f64,
    power_w: f64,
    h_tx_m: f64,
    h_rx_m: f64,
    mut rx_sens_dbm: f64,
    delta_h_m: f64,
    climate_idx: i32,
    n_0: f64,
    pol_idx: i32,
    epsilon: f64,
    sigma: f64,
    out_radius_m: *mut f64,
) -> i32 {
    if out_radius_m.is_null() {
        return 99; // ErrorInvalidPointer
    }

    if power_w <= 0.0 {
        unsafe {
            *out_radius_m = 0.0;
        }
        return 0;
    }

    // Default receiver sensitivity if unprovided
    if rx_sens_dbm == 0.0 {
        rx_sens_dbm = -90.0;
    }

    // Convert Watts to dBm
    // P(dBm) = 10 * log10(P(W) * 1000)
    let power_dbm = 10.0 * (power_w * 1000.0).log10();

    // Maximum allowable path loss (assuming 0 dBi antenna gains)
    // Loss = TxPower(dBm) - RxSens(dBm)
    let max_loss_db = power_dbm - rx_sens_dbm;

    let mut min_d_km = 0.001; // 1 meter
    let mut max_d_km = 2000.0; // ITM limit

    // Parallel multi-point search
    // We'll perform 10 iterations, each evaluating 8 points in parallel.
    // This gives roughly same or better precision than a 50rd-binary search but with fewer serial steps.
    for _ in 0..10 {
        let n_probes = 8;
        let step_size = (max_d_km - min_d_km) / (n_probes + 1) as f64;

        if step_size < 0.0001 {
            break; // Better than 1 meter precision reached
        }

        let probes: Vec<f64> = (1..=n_probes)
            .map(|i| min_d_km + step_size * i as f64)
            .collect();

        // Calculate losses for all probes in parallel
        let results: Vec<Option<f64>> = probes
            .par_iter()
            .map(|&d| {
                match itm_area_tls(
                    h_tx_m,
                    h_rx_m,
                    0,
                    0,
                    d,
                    delta_h_m,
                    climate_idx,
                    n_0,
                    f_mhz,
                    pol_idx,
                    epsilon,
                    sigma,
                    0,
                    50.0,
                    50.0,
                    50.0,
                ) {
                    Ok((loss, _, _)) => Some(loss),
                    Err(_) => None,
                }
            })
            .collect();

        // Find the furthest successful probe that's under the loss threshold
        let mut last_valid_idx: i32 = -1;
        for (i, res) in results.iter().enumerate() {
            if let Some(loss) = res {
                if *loss <= max_loss_db {
                    last_valid_idx = i as i32;
                } else {
                    break; // Losses usually increase with distance
                }
            } else {
                break;
            }
        }

        // Adjust search window
        if last_valid_idx == -1 {
            // All probes over threshold, must be in first segment
            max_d_km = probes[0];
        } else if last_valid_idx == (n_probes - 1) {
            // All probes under threshold, must be in last segment
            min_d_km = probes[last_valid_idx as usize];
        } else {
            // Threshold is between last_valid_idx and the next probe
            min_d_km = probes[last_valid_idx as usize];
            max_d_km = probes[last_valid_idx as usize + 1];
        }
    }

    unsafe {
        *out_radius_m = min_d_km * 1000.0;
    }
    0
}

/// C-API wrapper for the ITM Point-to-Point TLS function.
/// `pfl_data` must be a valid pointer to an array of `pfl_len` doubles.
/// The array format must match ITM expectations: [N-1, spacing_m, elev_1, elev_2, ... elev_N]
/// Returns 0 on success, or an error code > 0 on failure.
///
/// # Safety
/// - `pfl_data` must be non-null and point to a readable array of at least `pfl_len` `f64` values.
/// - `pfl_data` must be properly aligned for `f64` and remain valid for the duration of the call.
/// - `out_a_db` must be a non-null, properly aligned pointer to writable `f64` storage.
/// - `out_warnings` must be a non-null, properly aligned pointer to writable `u32` storage.
/// - Output pointers must remain valid for the duration of the call.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn itm_p2p_tls_c(
    h_tx_m: f64,
    h_rx_m: f64,
    pfl_data: *const f64,
    pfl_len: usize,
    climate_idx: i32,
    n_0: f64,
    f_mhz: f64,
    pol_idx: i32,
    epsilon: f64,
    sigma: f64,
    mdvar: i32,
    time_pct: f64,
    location_pct: f64,
    situation_pct: f64,
    out_a_db: *mut f64,
    out_warnings: *mut u32,
) -> i32 {
    if pfl_data.is_null() || out_a_db.is_null() || out_warnings.is_null() {
        return 99; // ErrorInvalidPointer
    }

    if pfl_len < 3 {
        return 98; // ErrorInvalidProfileArray
    }

    // Safely reconstruct the slice from the raw C pointer
    let pfl_slice = unsafe { slice::from_raw_parts(pfl_data, pfl_len) };

    match itm_p2p_tls(
        h_tx_m,
        h_rx_m,
        pfl_slice,
        climate_idx,
        n_0,
        f_mhz,
        pol_idx,
        epsilon,
        sigma,
        mdvar,
        time_pct,
        location_pct,
        situation_pct,
    ) {
        Ok((a_db, warnings, _inter)) => {
            unsafe {
                *out_a_db = a_db;
                *out_warnings = warnings;
            }
            0 // Success
        }
        Err(err) => err.code(),
    }
}

#[cfg(test)]
mod ffi_tests {
    use super::*;

    #[test]
    fn test_ffi_area_tls() {
        let mut a_db = 0.0;
        let mut warnings = 0u32;

        unsafe {
            let res = itm_area_tls_c(
                30.0,
                20.0,
                0,
                0,
                100.0,
                90.0,
                5,
                301.0,
                900.0,
                0,
                15.0,
                0.005,
                0,
                50.0,
                50.0,
                50.0,
                &mut a_db,
                &mut warnings,
            );
            assert_eq!(res, 0);
            assert!(a_db > 0.0);
        }
    }

    #[test]
    fn test_ffi_signal_radius() {
        let mut radius_m = 0.0;

        unsafe {
            let res = itm_calculate_signal_radius_c(
                900.0, // WiFi/cellular bandish 900 MHz
                1.0,   // 1 Watt transmission
                30.0,  // 30m Tx
                2.0,   // 2m Rx
                -90.0, // Sensitivity
                90.0,  // roughness
                5,     // Climate
                301.0,
                0,
                15.0,
                0.005,
                &mut radius_m,
            );

            assert_eq!(res, 0);
            assert!(radius_m > 0.0);
            // ~1 Watt at 900mhz to -90dBm is a fairly long distance (several km)
            assert!(radius_m > 1000.0);
        }
    }

    #[test]
    fn test_ffi_p2p_tls() {
        let mut a_db = 0.0;
        let mut warnings = 0u32;
        // Mock a basic array: 2 points, 100.0 spacing, 0.0 elev, 0.0 elev
        let pfl = [2.0, 100.0, 0.0, 0.0, 0.0];

        unsafe {
            let res = itm_p2p_tls_c(
                30.0,
                20.0,
                pfl.as_ptr(),
                pfl.len(),
                5,
                301.0,
                900.0,
                0,
                15.0,
                0.005,
                0,
                50.0,
                50.0,
                50.0,
                &mut a_db,
                &mut warnings,
            );
            assert_eq!(res, 0);
            assert!(a_db > 0.0);
        }
    }
}
