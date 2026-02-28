use std::f64;

use std::slice;

use crate::math::itm::{ItmError, itm_area_tls, itm_p2p_tls};

/// C-API wrapper for the ITM Area TLS function.
/// All pointers must be valid. Warnings are written to the `warnings` pointer.
///
/// Returns 0 on success, or an error code > 0 on failure.
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
        Err(ItmError::Success) => 0,
        Err(ItmError::SuccessWithWarnings) => 1,
        Err(ItmError::ErrorTxTerminalHeight) => 2,
        Err(ItmError::ErrorRxTerminalHeight) => 3,
        Err(ItmError::ErrorInvalidRadioClimate) => 4,
        Err(ItmError::ErrorRefractivity) => 5,
        Err(ItmError::ErrorFrequency) => 6,
        Err(ItmError::ErrorPolarization) => 7,
        Err(ItmError::ErrorEpsilon) => 8,
        Err(ItmError::ErrorSigma) => 9,
        Err(ItmError::ErrorMdvar) => 10,
        Err(ItmError::ErrorInvalidSituation) => 11,
        Err(ItmError::ErrorInvalidTime) => 12,
        Err(ItmError::ErrorInvalidLocation) => 13,
        Err(ItmError::ErrorSurfaceRefractivitySmall) => 14,
        Err(ItmError::ErrorSurfaceRefractivityLarge) => 15,
        Err(ItmError::ErrorEffectiveEarth) => 16,
        Err(ItmError::ErrorGroundImpedance) => 17,
        Err(ItmError::ErrorPathDistance) => 18,
        Err(ItmError::ErrorDeltaH) => 19,
        Err(ItmError::ErrorTxSitingCriteria) => 20,
        Err(ItmError::ErrorRxSitingCriteria) => 21,
        Err(ItmError::ErrorInvalidReliability) => 22,
        Err(ItmError::ErrorInvalidConfidence) => 23,
        Err(ItmError::Other(code)) => code,
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
/// * `rx_sens_dbm` - The sensitivity threshold of the receiver in dBm (e.g., -90.0). Defaults to -90.0 if you pass 0.0 for this value simply as a fallback.
/// * `delta_h_m` - Terrain irregularity parameter in meters (e.g. 90.0).
/// * `climate_idx` - Climate type (e.g., 5 for Continental Temperate).
/// * `n_0` - Surface refractivity (e.g. 301.0).
/// * `pol_idx` - Polarization (0 for Horizontal, 1 for Vertical).
/// * `epsilon` - Relative earth dielectric constant (e.g. 15.0).
/// * `sigma` - Earth surface conductivity (e.g. 0.005).
/// * `out_radius_m` - Output pointer where the resulting max radius in meters will be written.
///
/// Returns 0 on success, or an error code.
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

    let mut best_d_km = 0.0;

    // Binary search over distance (km) to match max_loss_db
    for _ in 0..50 {
        let mid_d_km = (min_d_km + max_d_km) / 2.0;

        match itm_area_tls(
            h_tx_m,
            h_rx_m,
            0, // Random tx siting
            0, // Random rx siting
            mid_d_km,
            delta_h_m,
            climate_idx,
            n_0,
            f_mhz,
            pol_idx,
            epsilon,
            sigma,
            0, // Broadcast mode
            50.0,
            50.0,
            50.0,
        ) {
            Ok((a_db, _, _)) => {
                // If attenuation is less than max allowable loss, signal still reaches
                if a_db <= max_loss_db {
                    best_d_km = mid_d_km;
                    min_d_km = mid_d_km;
                } else {
                    max_d_km = mid_d_km;
                }
            }
            Err(ItmError::Success) => return 0,
            Err(ItmError::SuccessWithWarnings) => return 1,
            Err(ItmError::ErrorTxTerminalHeight) => return 2,
            Err(ItmError::ErrorRxTerminalHeight) => return 3,
            Err(ItmError::ErrorInvalidRadioClimate) => return 4,
            Err(ItmError::ErrorRefractivity) => return 5,
            Err(ItmError::ErrorFrequency) => return 6,
            Err(ItmError::ErrorPolarization) => return 7,
            Err(ItmError::ErrorEpsilon) => return 8,
            Err(ItmError::ErrorSigma) => return 9,
            Err(ItmError::ErrorMdvar) => return 10,
            Err(ItmError::ErrorInvalidSituation) => return 11,
            Err(ItmError::ErrorInvalidTime) => return 12,
            Err(ItmError::ErrorInvalidLocation) => return 13,
            Err(ItmError::ErrorSurfaceRefractivitySmall) => return 14,
            Err(ItmError::ErrorSurfaceRefractivityLarge) => return 15,
            Err(ItmError::ErrorEffectiveEarth) => return 16,
            Err(ItmError::ErrorGroundImpedance) => return 17,
            Err(ItmError::ErrorPathDistance) => return 18,
            Err(ItmError::ErrorDeltaH) => return 19,
            Err(ItmError::ErrorTxSitingCriteria) => return 20,
            Err(ItmError::ErrorRxSitingCriteria) => return 21,
            Err(ItmError::ErrorInvalidReliability) => return 22,
            Err(ItmError::ErrorInvalidConfidence) => return 23,
            Err(ItmError::Other(code)) => return code,
        }

        // Close enough threshold (within 1 meter accuracy)
        if (max_d_km - min_d_km).abs() < 0.001 {
            break;
        }
    }

    unsafe {
        *out_radius_m = best_d_km * 1000.0;
    }
    0
}

/// C-API wrapper for the ITM Point-to-Point TLS function.
/// `pfl_data` must be a valid pointer to an array of `pfl_len` doubles.
/// The array format must match ITM expectations: [N-1, spacing_m, elev_1, elev_2, ... elev_N]
/// Returns 0 on success, or an error code > 0 on failure.
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
        Err(ItmError::Success) => 0,
        Err(ItmError::SuccessWithWarnings) => 1,
        Err(ItmError::ErrorTxTerminalHeight) => 2,
        Err(ItmError::ErrorRxTerminalHeight) => 3,
        Err(ItmError::ErrorInvalidRadioClimate) => 4,
        Err(ItmError::ErrorRefractivity) => 5,
        Err(ItmError::ErrorFrequency) => 6,
        Err(ItmError::ErrorPolarization) => 7,
        Err(ItmError::ErrorEpsilon) => 8,
        Err(ItmError::ErrorSigma) => 9,
        Err(ItmError::ErrorMdvar) => 10,
        Err(ItmError::ErrorInvalidSituation) => 11,
        Err(ItmError::ErrorInvalidTime) => 12,
        Err(ItmError::ErrorInvalidLocation) => 13,
        Err(ItmError::ErrorSurfaceRefractivitySmall) => 14,
        Err(ItmError::ErrorSurfaceRefractivityLarge) => 15,
        Err(ItmError::ErrorEffectiveEarth) => 16,
        Err(ItmError::ErrorGroundImpedance) => 17,
        Err(ItmError::ErrorPathDistance) => 18,
        Err(ItmError::ErrorDeltaH) => 19,
        Err(ItmError::ErrorTxSitingCriteria) => 20,
        Err(ItmError::ErrorRxSitingCriteria) => 21,
        Err(ItmError::ErrorInvalidReliability) => 22,
        Err(ItmError::ErrorInvalidConfidence) => 23,
        Err(ItmError::Other(code)) => code,
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
        let pfl = vec![2.0, 100.0, 0.0, 0.0, 0.0];

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
