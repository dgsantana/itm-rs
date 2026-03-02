#![allow(clippy::too_many_arguments)]
//! High-level ITM orchestration and input validation (Rust-first API).
//!
//! This module provides Rust-native implementations of the ITM wrapper
//! functions (area and point-to-point modes) and input validation. It
//! reuses lower-level math utilities defined elsewhere in the `math` module.

// use direct module functions from `crate::math::longley_rice` where needed
use crate::math::propagation::{free_space_loss, initialize_point_to_point};
use crate::math::terrain::{SitingCriteria, initialize_area};
use crate::math::variability::{Climate, VariabilityMode, variability_loss};

use std::fmt;

/// Error codes corresponding to ITM failures.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ItmError {
    Success = 0,
    SuccessWithWarnings = 1,
    ErrorTxTerminalHeight = 2,
    ErrorRxTerminalHeight = 3,
    ErrorInvalidRadioClimate = 4,
    ErrorRefractivity = 5,
    ErrorFrequency = 6,
    ErrorPolarization = 7,
    ErrorEpsilon = 8,
    ErrorSigma = 9,
    ErrorMdvar = 10,
    ErrorInvalidSituation = 11,
    ErrorInvalidTime = 12,
    ErrorInvalidLocation = 13,
    ErrorSurfaceRefractivitySmall = 14,
    ErrorSurfaceRefractivityLarge = 15,
    ErrorEffectiveEarth = 16,
    ErrorGroundImpedance = 17,
    ErrorPathDistance = 18,
    ErrorDeltaH = 19,
    ErrorTxSitingCriteria = 20,
    ErrorRxSitingCriteria = 21,
    ErrorInvalidReliability = 22,
    ErrorInvalidConfidence = 23,
    Other(i32) = 99,
}

impl ItmError {
    /// Convert an `ItmError` to its integer error code.
    pub const fn code(self) -> i32 {
        match self {
            ItmError::Success => 0,
            ItmError::SuccessWithWarnings => 1,
            ItmError::ErrorTxTerminalHeight => 2,
            ItmError::ErrorRxTerminalHeight => 3,
            ItmError::ErrorInvalidRadioClimate => 4,
            ItmError::ErrorRefractivity => 5,
            ItmError::ErrorFrequency => 6,
            ItmError::ErrorPolarization => 7,
            ItmError::ErrorEpsilon => 8,
            ItmError::ErrorSigma => 9,
            ItmError::ErrorMdvar => 10,
            ItmError::ErrorInvalidSituation => 11,
            ItmError::ErrorInvalidTime => 12,
            ItmError::ErrorInvalidLocation => 13,
            ItmError::ErrorSurfaceRefractivitySmall => 14,
            ItmError::ErrorSurfaceRefractivityLarge => 15,
            ItmError::ErrorEffectiveEarth => 16,
            ItmError::ErrorGroundImpedance => 17,
            ItmError::ErrorPathDistance => 18,
            ItmError::ErrorDeltaH => 19,
            ItmError::ErrorTxSitingCriteria => 20,
            ItmError::ErrorRxSitingCriteria => 21,
            ItmError::ErrorInvalidReliability => 22,
            ItmError::ErrorInvalidConfidence => 23,
            ItmError::Other(code) => code,
        }
    }
}

impl fmt::Display for ItmError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

/// Warning bitmask constants (u32).
pub mod warnings {
    pub const NO_WARNINGS: u32 = 0;
    pub const WARN_TX_TERMINAL_HEIGHT: u32 = 1 << 0;
    pub const WARN_RX_TERMINAL_HEIGHT: u32 = 1 << 1;
    pub const WARN_FREQUENCY: u32 = 1 << 2;
    pub const WARN_TX_HORIZON_ANGLE: u32 = 1 << 3;
    pub const WARN_RX_HORIZON_ANGLE: u32 = 1 << 4;
    pub const WARN_TX_HORIZON_DISTANCE_1: u32 = 1 << 5;
    pub const WARN_RX_HORIZON_DISTANCE_1: u32 = 1 << 6;
    pub const WARN_TX_HORIZON_DISTANCE_2: u32 = 1 << 7;
    pub const WARN_RX_HORIZON_DISTANCE_2: u32 = 1 << 8;
    pub const WARN_SURFACE_REFRACTIVITY: u32 = 1 << 9;
    pub const WARN_PATH_DISTANCE_TOO_SMALL_1: u32 = 1 << 10;
    pub const WARN_PATH_DISTANCE_TOO_SMALL_2: u32 = 1 << 11;
    pub const WARN_PATH_DISTANCE_TOO_BIG_1: u32 = 1 << 12;
    pub const WARN_PATH_DISTANCE_TOO_BIG_2: u32 = 1 << 13;
}

/// Intermediate values structure captured by the Ex variants.
///
/// Field names use idiomatic snake_case. Compatibility accessor methods for the
/// original field-style names are provided below to aid incremental migration.
#[derive(Debug, Clone)]
pub struct IntermediateValues {
    pub d_km: f64,
    pub a_ref_db: f64,
    pub a_fs_db: f64,
    pub delta_h_meter: f64,
    pub d_hzn_meter: [f64; 2],
    pub h_e_meter: [f64; 2],
    pub n_s: f64,
    pub theta_hzn: [f64; 2],
    pub mode: i32,
}

impl Default for IntermediateValues {
    fn default() -> Self {
        Self {
            d_km: 0.0,
            a_ref_db: 0.0,
            a_fs_db: 0.0,
            delta_h_meter: 0.0,
            d_hzn_meter: [0.0, 0.0],
            h_e_meter: [0.0, 0.0],
            n_s: 0.0,
            theta_hzn: [0.0, 0.0],
            mode: 0,
        }
    }
}

/// Validate inputs for both AREA and P2P modes.
///
/// Returns `Ok(warnings_mask)` if inputs pass (warnings may be set), or
/// `Err(ItmError)` for fatal validation failures.
pub fn validate_inputs(
    h_tx_meter: f64,
    h_rx_meter: f64,
    climate_idx_one_based: i32,
    time_pct: f64,
    location_pct: f64,
    situation_pct: f64,
    n_0: f64,
    f_mhz: f64,
    pol_idx: i32,
    epsilon: f64,
    sigma: f64,
    mdvar: i32,
) -> Result<u32, ItmError> {
    let mut warns = warnings::NO_WARNINGS;

    if !(1.0..=1000.0).contains(&h_tx_meter) {
        warns |= warnings::WARN_TX_TERMINAL_HEIGHT;
    }
    if !(0.5..=3000.0).contains(&h_tx_meter) {
        return Err(ItmError::ErrorTxTerminalHeight);
    }

    if !(1.0..=1000.0).contains(&h_rx_meter) {
        warns |= warnings::WARN_RX_TERMINAL_HEIGHT;
    }
    if !(0.5..=3000.0).contains(&h_rx_meter) {
        return Err(ItmError::ErrorRxTerminalHeight);
    }

    // climate indices are 1..=7 per NTIA/ITM
    if !(1..=7).contains(&climate_idx_one_based) {
        return Err(ItmError::ErrorInvalidRadioClimate);
    }

    if !(250.0..=400.0).contains(&n_0) {
        return Err(ItmError::ErrorRefractivity);
    }

    if !(40.0..=10000.0).contains(&f_mhz) {
        warns |= warnings::WARN_FREQUENCY;
    }
    if !(20.0..=20000.0).contains(&f_mhz) {
        return Err(ItmError::ErrorFrequency);
    }

    if !(pol_idx == 0 || pol_idx == 1) {
        return Err(ItmError::ErrorPolarization);
    }

    if epsilon < 1.0 {
        return Err(ItmError::ErrorEpsilon);
    }
    if sigma <= 0.0 {
        return Err(ItmError::ErrorSigma);
    }

    if (mdvar < 0)
        || (mdvar > 3 && mdvar < 10)
        || (mdvar > 13 && mdvar < 20)
        || (mdvar > 23 && mdvar < 30)
        || (mdvar > 33)
    {
        return Err(ItmError::ErrorMdvar);
    }

    if situation_pct <= 0.0 || situation_pct >= 100.0 {
        return Err(ItmError::ErrorInvalidSituation);
    }
    if time_pct <= 0.0 || time_pct >= 100.0 {
        return Err(ItmError::ErrorInvalidTime);
    }
    if location_pct <= 0.0 || location_pct >= 100.0 {
        return Err(ItmError::ErrorInvalidLocation);
    }

    Ok(warns)
}

/// Public Rust-first AREA TLS wrapper.
///
/// Returns (A_db, warnings_mask, IntermediateValues) or `ItmError`.
pub fn itm_area_tls(
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
) -> Result<(f64, u32, IntermediateValues), ItmError> {
    let mut warnings_mask = warnings::NO_WARNINGS;
    // validate
    match validate_inputs(
        h_tx_m,
        h_rx_m,
        climate_idx,
        time_pct,
        location_pct,
        situation_pct,
        n_0,
        f_mhz,
        pol_idx,
        epsilon,
        sigma,
        mdvar,
    ) {
        Ok(w) => warnings_mask |= w,
        Err(e) => return Err(e),
    }

    // area-specific validation
    if d_km <= 0.0 {
        return Err(ItmError::ErrorPathDistance);
    }
    if delta_h_m < 0.0 {
        return Err(ItmError::ErrorDeltaH);
    }
    if !(0..=2).contains(&tx_site_criteria_idx) {
        return Err(ItmError::ErrorTxSitingCriteria);
    }
    if !(0..=2).contains(&rx_site_criteria_idx) {
        return Err(ItmError::ErrorRxSitingCriteria);
    }

    let site_criteria = [tx_site_criteria_idx, rx_site_criteria_idx];
    let h_meter = [h_tx_m, h_rx_m];
    let mut inter = IntermediateValues {
        d_km,
        ..Default::default()
    };

    // Point-to-point initialization (ground impedance, gamma_e, N_s)
    let (z_g, gamma_e, n_s) = initialize_point_to_point(
        f_mhz,
        0.0,
        n_0,
        if pol_idx == 1 {
            crate::math::propagation::Polarization::Vertical
        } else {
            crate::math::propagation::Polarization::Horizontal
        },
        epsilon,
        sigma,
    );

    // initialize_area to compute effective heights and horizons
    let (h_e_m, d_hzn_m, theta_hzn) = initialize_area(
        [
            if site_criteria[0] == 0 {
                SitingCriteria::Random
            } else if site_criteria[0] == 1 {
                SitingCriteria::Careful
            } else {
                SitingCriteria::VeryCareful
            },
            if site_criteria[1] == 0 {
                SitingCriteria::Random
            } else if site_criteria[1] == 1 {
                SitingCriteria::Careful
            } else {
                SitingCriteria::VeryCareful
            },
        ],
        gamma_e,
        delta_h_m,
        h_meter,
    );

    // convert to meters
    let d_m = d_km * 1000.0;

    // Call Longley-Rice
    let mut warnings_u32 = warnings_mask;
    let lr_result = crate::math::longley_rice::longley_rice(
        theta_hzn, f_mhz, z_g, d_hzn_m, h_e_m, gamma_e, n_s, delta_h_m, h_meter, d_m, 0,
    )?;
    let a_ref_db = lr_result.a_ref;
    warnings_u32 |= lr_result.warnings;
    let propmode = lr_result.propmode;

    let a_fs_db = free_space_loss(d_m, f_mhz);

    // Extract mdvar flags for disabling variability components
    // mdvar >= 20: disable situation variability
    // mdvar >= 10 (after subtracting 20): disable location variability
    let disable_situation = mdvar >= 20;
    let disable_location = (mdvar % 20) >= 10;

    let result_variability = variability_loss(
        time_pct,
        location_pct,
        situation_pct,
        h_e_m,
        delta_h_m,
        f_mhz,
        d_m,
        a_ref_db,
        climate_from_idx(climate_idx)?,
        mdvar_to_mode(mdvar)?,
        disable_location,
        disable_situation,
        &mut warnings_u32,
    );

    let a_db = a_fs_db + result_variability;

    // fill inter values (snake_case fields)
    inter.a_ref_db = a_ref_db;
    inter.a_fs_db = a_fs_db;
    inter.delta_h_meter = delta_h_m;
    inter.d_hzn_meter = d_hzn_m;
    inter.h_e_meter = h_e_m;
    inter.n_s = n_s;
    inter.theta_hzn = theta_hzn;
    inter.mode = propmode;

    Ok((a_db, warnings_u32, inter))
}

/// AREA CR wrapper (confidence/reliability). Maps CR to TLS internals.
pub fn itm_area_cr(
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
    confidence_pct: f64,
    reliability_pct: f64,
) -> Result<(f64, u32, IntermediateValues), ItmError> {
    // CR: reliability maps to time, confidence maps to situation per original API
    let res = itm_area_tls(
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
        reliability_pct, // time
        50.0,            // location hard-coded per original wrapper
        confidence_pct,  // situation
    );

    match res {
        Err(ItmError::ErrorInvalidTime) => Err(ItmError::ErrorInvalidReliability),
        Err(ItmError::ErrorInvalidSituation) => Err(ItmError::ErrorInvalidConfidence),
        other => other,
    }
}

/// P2P TLS wrapper.
pub fn itm_p2p_tls(
    h_tx_m: f64,
    h_rx_m: f64,
    pfl: &[f64],
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
) -> Result<(f64, u32, IntermediateValues), ItmError> {
    // Perform validate inputs
    match validate_inputs(
        h_tx_m,
        h_rx_m,
        climate_idx,
        time_pct,
        location_pct,
        situation_pct,
        n_0,
        f_mhz,
        pol_idx,
        epsilon,
        sigma,
        mdvar,
    ) {
        Ok(w) => {
            // Continue
            let mut inter = IntermediateValues::default();
            // Use QuickPfl-like helper (delegate to extracted module)
            let q_pfl = crate::math::terrain::quick_pfl(pfl, 1.0 / 157e-9, [h_tx_m, h_rx_m]);
            let (theta_hzn, d_hzn_m, h_e_m, delta_h_m, d_m) = (
                q_pfl.theta_hzn,
                q_pfl.d_hzn,
                q_pfl.h_e,
                q_pfl.delta_h,
                q_pfl.d,
            );
            // initialize point-to-point
            let (z_g, gamma_e, n_s) = initialize_point_to_point(
                f_mhz,
                0.0,
                n_0,
                if pol_idx == 1 {
                    crate::math::propagation::Polarization::Vertical
                } else {
                    crate::math::propagation::Polarization::Horizontal
                },
                epsilon,
                sigma,
            );
            let mut warnings = w;
            let lr_result = crate::math::longley_rice::longley_rice(
                theta_hzn,
                f_mhz,
                z_g,
                d_hzn_m,
                h_e_m,
                gamma_e,
                n_s,
                delta_h_m,
                [h_tx_m, h_rx_m],
                d_m,
                0,
            )?;
            let a_ref_db = lr_result.a_ref;
            warnings |= lr_result.warnings;
            let propmode = lr_result.propmode;
            let a_fs_db = free_space_loss(d_m, f_mhz);

            // Extract mdvar flags for disabling variability components
            // mdvar >= 20: disable situation variability
            // mdvar >= 10 (after subtracting 20): disable location variability
            let disable_situation = mdvar >= 20;
            let disable_location = (mdvar % 20) >= 10;

            let variability = variability_loss(
                time_pct,
                location_pct,
                situation_pct,
                h_e_m,
                delta_h_m,
                f_mhz,
                d_m,
                a_ref_db,
                climate_from_idx(climate_idx)?,
                mdvar_to_mode(mdvar)?,
                disable_location,
                disable_situation,
                &mut warnings,
            );
            let a_db = a_fs_db + variability;
            inter.a_ref_db = a_ref_db;
            inter.a_fs_db = a_fs_db;
            inter.delta_h_meter = delta_h_m;
            inter.d_hzn_meter = d_hzn_m;
            inter.h_e_meter = h_e_m;
            inter.n_s = n_s;
            inter.theta_hzn = theta_hzn;
            inter.mode = propmode;
            Ok((a_db, warnings, inter))
        }
        Err(e) => Err(e),
    }
}

/// P2P CR wrapper (maps CR to TLS)
pub fn itm_p2p_cr(
    h_tx_m: f64,
    h_rx_m: f64,
    pfl: &[f64],
    climate_idx: i32,
    n_0: f64,
    f_mhz: f64,
    pol_idx: i32,
    epsilon: f64,
    sigma: f64,
    mdvar: i32,
    confidence_pct: f64,
    reliability_pct: f64,
) -> Result<(f64, u32, IntermediateValues), ItmError> {
    let res = itm_p2p_tls(
        h_tx_m,
        h_rx_m,
        pfl,
        climate_idx,
        n_0,
        f_mhz,
        pol_idx,
        epsilon,
        sigma,
        mdvar,
        reliability_pct,
        50.0,
        confidence_pct,
    );

    match res {
        Err(ItmError::ErrorInvalidTime) => Err(ItmError::ErrorInvalidReliability),
        Err(ItmError::ErrorInvalidSituation) => Err(ItmError::ErrorInvalidConfidence),
        other => other,
    }
}

/// Convert 1-based climate index to `Climate` enum.
fn climate_from_idx(idx_one_based: i32) -> Result<Climate, ItmError> {
    match idx_one_based {
        1 => Ok(Climate::Equatorial),
        2 => Ok(Climate::ContinentalSubtropical),
        3 => Ok(Climate::MaritimeSubtropical),
        4 => Ok(Climate::Desert),
        5 => Ok(Climate::ContinentalTemperate),
        6 => Ok(Climate::MaritimeTemperateOverLand),
        7 => Ok(Climate::MaritimeTemperateOverSea),
        _ => Err(ItmError::ErrorInvalidRadioClimate),
    }
}

/// Convert mdvar integer into VariabilityMode, preserving common encodings.
fn mdvar_to_mode(mdvar: i32) -> Result<VariabilityMode, ItmError> {
    // Simplified mapping: 0 -> Broadcast, 1 -> Mobile, 2 -> Accidental, 3 -> SingleMessage
    // Extended encodings (+10/+20) handled elsewhere in variability module; here we only map base mode.
    match mdvar % 10 {
        0 => Ok(VariabilityMode::Broadcast),
        1 => Ok(VariabilityMode::Mobile),
        2 => Ok(VariabilityMode::Accidental),
        3 => Ok(VariabilityMode::SingleMessage),
        _ => Ok(VariabilityMode::Broadcast),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_validate_inputs_ok() {
        let res = validate_inputs(
            10.0, 10.0, 5, 50.0, 50.0, 50.0, 301.0, 900.0, 0, 15.0, 0.005, 0,
        );
        assert!(res.is_ok());
    }

    #[test]
    fn test_itm_area_tls_runs() {
        let res = itm_area_tls(
            30.0, 20.0, 0, 0, 100.0, 90.0, 5, 301.0, 900.0, 0, 15.0, 0.005, 0, 50.0, 50.0, 50.0,
        );
        assert!(res.is_ok());
        let (a_db, _warnings, inter) = res.unwrap();
        assert!(a_db.is_finite());
        assert!(inter.a_ref_db.is_finite());
    }
}
