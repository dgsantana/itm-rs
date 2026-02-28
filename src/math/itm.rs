//! High-level ITM orchestration and input validation (Rust-first API).
//!
//! This module provides Rust-native implementations of the ITM wrapper
//! functions (area and point-to-point modes) and input validation. It
//! reuses lower-level math utilities defined elsewhere in the `math` module.

use crate::math::propagation::free_space_loss;
use crate::math::propagation::initialize_point_to_point;
use crate::math::terrain::{SitingCriteria, compute_delta_h, find_horizons, initialize_area};
use crate::math::variability::{Climate, VariabilityMode, variability_loss}; // variability_warnings removed

use num_complex::Complex;
use std::fmt;

/// Error codes corresponding to ITM failures.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ItmError {
    Success,
    SuccessWithWarnings,
    ErrorTxTerminalHeight,
    ErrorRxTerminalHeight,
    ErrorInvalidRadioClimate,
    ErrorRefractivity,
    ErrorFrequency,
    ErrorPolarization,
    ErrorEpsilon,
    ErrorSigma,
    ErrorMdvar,
    ErrorInvalidSituation,
    ErrorInvalidTime,
    ErrorInvalidLocation,
    ErrorSurfaceRefractivitySmall,
    ErrorSurfaceRefractivityLarge,
    ErrorEffectiveEarth,
    ErrorGroundImpedance,
    ErrorPathDistance,
    ErrorDeltaH,
    ErrorTxSitingCriteria,
    ErrorRxSitingCriteria,
    ErrorInvalidReliability,
    ErrorInvalidConfidence,
    Other(i32),
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
#[derive(Debug, Clone)]
pub struct IntermediateValues {
    pub d__km: f64,
    pub A_ref__db: f64,
    pub A_fs__db: f64,
    pub delta_h__meter: f64,
    pub d_hzn__meter: [f64; 2],
    pub h_e__meter: [f64; 2],
    pub N_s: f64,
    pub theta_hzn: [f64; 2],
    pub mode: i32,
}

impl Default for IntermediateValues {
    fn default() -> Self {
        Self {
            d__km: 0.0,
            A_ref__db: 0.0,
            A_fs__db: 0.0,
            delta_h__meter: 0.0,
            d_hzn__meter: [0.0, 0.0],
            h_e__meter: [0.0, 0.0],
            N_s: 0.0,
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

    if h_tx_meter < 1.0 || h_tx_meter > 1000.0 {
        warns |= warnings::WARN_TX_TERMINAL_HEIGHT;
    }
    if h_tx_meter < 0.5 || h_tx_meter > 3000.0 {
        return Err(ItmError::ErrorTxTerminalHeight);
    }

    if h_rx_meter < 1.0 || h_rx_meter > 1000.0 {
        warns |= warnings::WARN_RX_TERMINAL_HEIGHT;
    }
    if h_rx_meter < 0.5 || h_rx_meter > 3000.0 {
        return Err(ItmError::ErrorRxTerminalHeight);
    }

    // climate indices are 1..=7 per NTIA/ITM
    if !(1..=7).contains(&climate_idx_one_based) {
        return Err(ItmError::ErrorInvalidRadioClimate);
    }

    if n_0 < 250.0 || n_0 > 400.0 {
        return Err(ItmError::ErrorRefractivity);
    }

    if f_mhz < 40.0 || f_mhz > 10000.0 {
        warns |= warnings::WARN_FREQUENCY;
    }
    if f_mhz < 20.0 || f_mhz > 20000.0 {
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

// -- Helper stubs for components not yet ported --

/// QuickPfl: quick profile preprocess. For now, a lightweight implementation
/// that derives horizons and delta_h using already-ported functions.
fn quick_pfl(
    pfl: &[f64],
    gamma_e: f64,
    h_meter: [f64; 2],
) -> ([f64; 2], [f64; 2], [f64; 2], f64, f64) {
    // returns (theta_hzn, d_hzn, h_e, delta_h, d_meter)
    let d_meter = pfl[0] * pfl[1];
    let (theta_hzn, d_hzn) = find_horizons(pfl, gamma_e, h_meter);
    // use initialize_area with random siting to derive h_e approximation
    let site_criteria = [SitingCriteria::Random, SitingCriteria::Random];
    let (h_e, _d_hzn2, _theta_hzn2) = initialize_area(site_criteria, gamma_e, 0.0, h_meter);
    let delta_h = compute_delta_h(pfl, 0.0, d_meter);
    (theta_hzn, d_hzn, h_e, delta_h, d_meter)
}

/// Placeholder for Longley-Rice core algorithm. This must be replaced by the
/// full Rust port. For now it returns a plausible A_ref__db.
/// Wrapper for DiffractionLoss - currently maps to Vogler smooth-earth diffraction.
fn diffraction_loss(
    d_m: f64,
    d_hzn_m: [f64; 2],
    h_e_m: [f64; 2],
    z_g: Complex<f64>,
    a_e_m: f64,
    delta_h_m: f64,
    _h_m: [f64; 2],
    _mode: i32,
    theta_los: f64,
    d_sml_m: f64,
    f_mhz: f64,
) -> f64 {
    // Use Vogler 3-radii smooth-earth diffraction as an approximation for DiffractionLoss.
    crate::math::diffraction::smooth_earth_diffraction(
        d_m, f_mhz, a_e_m, theta_los, d_hzn_m, h_e_m, z_g,
    )
}

/// Wrapper for LineOfSightLoss. For now we approximate with diffraction wrapper at small distances.
fn line_of_sight_loss(
    d_m: f64,
    h_e_m: [f64; 2],
    z_g: Complex<f64>,
    delta_h_m: f64,
    m_d: f64,
    a_d0_db: f64,
    d_sml_m: f64,
    f_mhz: f64,
) -> f64 {
    // Simple approximation: use diffraction_loss as a proxy.
    diffraction_loss(
        d_m,
        [d_sml_m / 2.0, d_sml_m / 2.0],
        h_e_m,
        z_g,
        d_sml_m / 2.0,
        delta_h_m,
        [0.0, 0.0],
        0,
        0.0,
        d_sml_m,
        f_mhz,
    )
}

fn longley_rice_(
    theta_hzn: [f64; 2],
    f_mhz: f64,
    z_g: Complex<f64>,
    d_hzn_m: [f64; 2],
    h_e_m: [f64; 2],
    gamma_e: f64,
    n_s: f64,
    delta_h_m: f64,
    h_m: [f64; 2],
    d_m: f64,
    mode: i32,
    warnings: &mut u32,
    propmode: &mut i32,
) -> Result<f64, ItmError> {
    // Bring in troposcatter and diffraction functions
    use crate::math::scatter::troposcatter_loss;

    // effective earth radius
    let a_e_m = 1.0 / gamma_e;

    // Terrestrial smooth earth horizon distances
    let mut d_hzn_s_m = [0.0f64; 2];
    for i in 0..2 {
        d_hzn_s_m[i] = (2.0 * h_e_m[i] * a_e_m).sqrt();
    }

    let d_sml_m = d_hzn_s_m[0] + d_hzn_s_m[1];
    let d_ml_m = d_hzn_m[0] + d_hzn_m[1];

    let theta_los = -((theta_hzn[0] + theta_hzn[1]).max(-d_ml_m / a_e_m));

    if theta_hzn[0].abs() > 200e-3 {
        *warnings |= warnings::WARN_TX_HORIZON_ANGLE;
    }
    if theta_hzn[1].abs() > 200e-3 {
        *warnings |= warnings::WARN_RX_HORIZON_ANGLE;
    }

    if d_hzn_m[0] < 0.1 * d_hzn_s_m[0] {
        *warnings |= warnings::WARN_TX_HORIZON_DISTANCE_1;
    }
    if d_hzn_m[1] < 0.1 * d_hzn_s_m[1] {
        *warnings |= warnings::WARN_RX_HORIZON_DISTANCE_1;
    }

    if d_hzn_m[0] > 3.0 * d_hzn_s_m[0] {
        *warnings |= warnings::WARN_TX_HORIZON_DISTANCE_2;
    }
    if d_hzn_m[1] > 3.0 * d_hzn_s_m[1] {
        *warnings |= warnings::WARN_RX_HORIZON_DISTANCE_2;
    }

    if n_s < 150.0 {
        return Err(ItmError::ErrorSurfaceRefractivitySmall);
    }
    if n_s > 400.0 {
        return Err(ItmError::ErrorSurfaceRefractivityLarge);
    }
    if n_s < 250.0 {
        *warnings |= warnings::WARN_SURFACE_REFRACTIVITY;
    }

    if a_e_m < 4_000_000.0 || a_e_m > 13_333_333.0 {
        return Err(ItmError::ErrorEffectiveEarth);
    }

    if z_g.re <= z_g.im.abs() {
        return Err(ItmError::ErrorGroundImpedance);
    }

    let d_3_m = d_sml_m.max(d_ml_m + 5.0 * ((a_e_m.powi(2) / f_mhz).powf(1.0 / 3.0)));
    let d_4_m = d_3_m + 10.0 * ((a_e_m.powi(2) / f_mhz).powf(1.0 / 3.0));

    let a_3_db = diffraction_loss(
        d_3_m, d_hzn_m, h_e_m, z_g, a_e_m, delta_h_m, h_m, mode, theta_los, d_sml_m, f_mhz,
    );
    let a_4_db = diffraction_loss(
        d_4_m, d_hzn_m, h_e_m, z_g, a_e_m, delta_h_m, h_m, mode, theta_los, d_sml_m, f_mhz,
    );

    let m_d = (a_4_db - a_3_db) / (d_4_m - d_3_m);
    let a_d0_db = a_3_db - m_d * d_3_m;

    let d_min_m = (h_e_m[0] - h_e_m[1]).abs() / 200e-3;
    if d_m < d_min_m {
        *warnings |= warnings::WARN_PATH_DISTANCE_TOO_SMALL_1;
    }
    if d_m < 1e3 {
        *warnings |= warnings::WARN_PATH_DISTANCE_TOO_SMALL_2;
    }
    if d_m > 1000e3 {
        *warnings |= warnings::WARN_PATH_DISTANCE_TOO_BIG_1;
    }
    if d_m > 2000e3 {
        *warnings |= warnings::WARN_PATH_DISTANCE_TOO_BIG_2;
    }

    let mut a_ref_db: f64;
    let mut propm = 0i32;

    if d_m < d_sml_m {
        let a_sml_db = d_sml_m * m_d + a_d0_db;

        let mut d_0_m = 0.04 * f_mhz * h_e_m[0] * h_e_m[1];
        let d_1_m: f64;
        if a_d0_db >= 0.0 {
            d_0_m = d_0_m.min(0.5 * d_ml_m);
            d_1_m = d_0_m + 0.25 * (d_ml_m - d_0_m);
        } else {
            d_1_m = (-a_d0_db / m_d).max(0.25 * d_ml_m);
        }

        let a_1_db = line_of_sight_loss(d_1_m, h_e_m, z_g, delta_h_m, m_d, a_d0_db, d_sml_m, f_mhz);

        let mut flag = false;
        let mut khat1 = 0.0f64;
        let mut khat2 = 0.0f64;

        if d_0_m < d_1_m {
            let a_0_db =
                line_of_sight_loss(d_0_m, h_e_m, z_g, delta_h_m, m_d, a_d0_db, d_sml_m, f_mhz);
            let q = (d_sml_m / d_0_m).ln();

            let denom = (d_sml_m - d_0_m) * (d_1_m / d_0_m).ln() - (d_1_m - d_0_m) * q;
            if denom.abs() > 0.0 {
                khat2 = ((d_sml_m - d_0_m) * (a_1_db - a_0_db)
                    - (d_1_m - d_0_m) * (a_sml_db - a_0_db))
                    / denom;
            }
            if khat2.is_nan() || khat2 < 0.0 {
                khat2 = 0.0;
            }

            flag = a_d0_db > 0.0 || khat2 > 0.0;

            if flag {
                khat1 = (a_sml_db - a_0_db - khat2 * q) / (d_sml_m - d_0_m);
                if khat1 < 0.0 {
                    khat1 = 0.0;
                    khat2 = (a_sml_db - a_0_db).abs() / q;
                    if khat2 == 0.0 {
                        khat1 = m_d;
                    }
                }
            }
        }

        if !flag {
            khat1 = (a_sml_db - a_1_db).abs() / (d_sml_m - d_1_m);
            khat2 = 0.0;
            if khat1 == 0.0 {
                khat1 = m_d;
            }
        }

        let a_o_db = a_sml_db - khat1 * d_sml_m - khat2 * d_sml_m.ln();
        a_ref_db = a_o_db + khat1 * d_m + khat2 * d_m.ln();
        propm = 1; // LINE_OF_SIGHT
    } else {
        // trans-horizon path
        let d_5_m = d_ml_m + 200e3;
        let d_6_m = d_ml_m + 400e3;

        let mut h0 = -1.0f64;
        let a6_db = troposcatter_loss(
            d_6_m, &theta_hzn, &d_hzn_m, &h_e_m, a_e_m, n_s, f_mhz, theta_los, &mut h0,
        );
        let a5_db = troposcatter_loss(
            d_5_m, &theta_hzn, &d_hzn_m, &h_e_m, a_e_m, n_s, f_mhz, theta_los, &mut h0,
        );

        let (m_s, a_s0_db, d_x_m) = if a5_db < 1000.0 {
            let m_s = (a6_db - a5_db) / 200e3;
            let mut d_x = d_ml_m
                .max(d_ml_m + 1.088 * (a_e_m.powi(2) / f_mhz).powf(1.0 / 3.0) * f_mhz.ln().abs());
            d_x = d_x.max((a5_db - a_d0_db - m_s * d_5_m) / (m_d - m_s));
            let a_s0 = (m_d - m_s) * d_x + a_d0_db;
            (m_s, a_s0, d_x)
        } else {
            (m_d, a_d0_db, 10e6)
        };

        if d_m > d_x_m {
            a_ref_db = m_s * d_m + a_s0_db;
            propm = 3; // TROPOSCATTER
        } else {
            a_ref_db = m_d * d_m + a_d0_db;
            propm = 2; // DIFFRACTION
        }
    }

    if a_ref_db < 0.0 {
        a_ref_db = 0.0;
    }
    *propmode = propm;
    Ok(a_ref_db)
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

    let site_criteria = [tx_site_criteria_idx as i32, rx_site_criteria_idx as i32];
    let h_meter = [h_tx_m, h_rx_m];
    let mut inter = IntermediateValues::default();
    inter.d__km = d_km;

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
    let mut propmode = 0i32;
    let a_ref_db = longley_rice_(
        theta_hzn,
        f_mhz,
        z_g,
        d_hzn_m,
        h_e_m,
        gamma_e,
        n_s,
        delta_h_m,
        h_meter,
        d_m,
        0,
        &mut warnings_u32,
        &mut propmode,
    )?;

    let a_fs_db = free_space_loss(d_m, f_mhz);
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
        false,
        false,
        &mut warnings_u32,
    );

    let a_db = a_fs_db + result_variability;

    // fill inter values
    inter.A_ref__db = a_ref_db;
    inter.A_fs__db = a_fs_db;
    inter.delta_h__meter = delta_h_m;
    inter.d_hzn__meter = d_hzn_m;
    inter.h_e__meter = h_e_m;
    inter.N_s = n_s;
    inter.theta_hzn = theta_hzn;
    inter.mode = propmode;

    let final_warnings = warnings_u32;
    if final_warnings != warnings::NO_WARNINGS {
        Ok((a_db, final_warnings, inter))
    } else {
        Ok((a_db, final_warnings, inter))
    }
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
            // Use QuickPfl-like helper
            let (theta_hzn, d_hzn_m, h_e_m, delta_h_m, d_m) =
                quick_pfl(pfl, 1.0 / 157e-9, [h_tx_m, h_rx_m]);
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
            let mut propmode = 0i32;
            let a_ref_db = longley_rice_(
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
                &mut warnings,
                &mut propmode,
            )?;
            let a_fs_db = free_space_loss(d_m, f_mhz);
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
                false,
                false,
                &mut warnings,
            );
            let a_db = a_fs_db + variability;
            inter.A_ref__db = a_ref_db;
            inter.A_fs__db = a_fs_db;
            inter.delta_h__meter = delta_h_m;
            inter.d_hzn__meter = d_hzn_m;
            inter.h_e__meter = h_e_m;
            inter.N_s = n_s;
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
        let (a_db, warnings, inter) = res.unwrap();
        assert!(a_db.is_finite());
        assert!(inter.A_ref__db.is_finite());
    }
}
