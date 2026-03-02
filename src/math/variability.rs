#![allow(clippy::too_many_arguments)]
//!! Variability loss calculations for ITM.
//!!
//!! This module contains functions for time, location, and situation variability
//!! computations used to derive statistical propagation losses.
use super::statistics::inverse_ccdf;
use super::terrain::terrain_roughness;

const THIRD: f64 = 1.0 / 3.0;
const A_9000_M: f64 = 9_000_000.0;

/// Radio climate classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Climate {
    Equatorial,
    ContinentalSubtropical,
    MaritimeSubtropical,
    Desert,
    ContinentalTemperate,
    MaritimeTemperateOverLand,
    MaritimeTemperateOverSea,
}

impl Climate {
    fn index(self) -> usize {
        match self {
            Climate::Equatorial => 0,
            Climate::ContinentalSubtropical => 1,
            Climate::MaritimeSubtropical => 2,
            Climate::Desert => 3,
            Climate::ContinentalTemperate => 4,
            Climate::MaritimeTemperateOverLand => 5,
            Climate::MaritimeTemperateOverSea => 6,
        }
    }
}

/// Variability mode for statistical calculations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariabilityMode {
    Broadcast,
    Mobile,
    Accidental,
    SingleMessage,
}

/// Warning flags returned by variability calculations.
pub mod variability_warnings {
    /// Extreme variability quantiles detected (|z| > 3.10).
    pub const EXTREME_VARIABILITIES: u32 = 1 << 0;
}

/// Curve helper function for TN101v2 Eqn III.69 & III.70.
///
/// # Arguments
///
/// * `c1` - Curve fit parameter.
/// * `c2` - Curve fit parameter.
/// * `x1` - Curve fit parameter.
/// * `x2` - Curve fit parameter.
/// * `x3` - Curve fit parameter.
/// * `d_e_m` - Effective distance in meters.
///
/// # Returns
///
/// Curve value in dB.
///
/// # References
///
/// - ITM Technical Note 101 (TN101v2), Eqn III.69 and III.70.
pub fn curve(c1: f64, c2: f64, x1: f64, x2: f64, x3: f64, d_e_m: f64) -> f64 {
    let ratio = d_e_m / x1;
    let ratio_sq = ratio * ratio;
    let denom = 1.0 + ((d_e_m - x2) / x3).powi(2);
    (c1 + c2 / denom) * ratio_sq / (1.0 + ratio_sq)
}

/// Computes the variability loss.
///
/// # Arguments
///
/// * `time_percent` - Time percentage (0 < time < 100).
/// * `location_percent` - Location percentage (0 < location < 100).
/// * `situation_percent` - Situation percentage (0 < situation < 100).
/// * `h_e_m` - Effective antenna heights in meters `[tx, rx]`.
/// * `delta_h_m` - Terrain irregularity parameter in meters.
/// * `f_mhz` - Frequency in MHz.
/// * `d_m` - Path distance in meters.
/// * `a_ref_db` - Reference attenuation in dB.
/// * `climate` - Radio climate classification.
/// * `mode` - Variability mode.
/// * `disable_location_variability` - If true, location variability is suppressed.
/// * `disable_situation_variability` - If true, situation variability is suppressed.
/// * `warnings` - Bitmask of warning flags updated by this function.
///
/// # Returns
///
/// Variability loss in dB.
///
/// # References
///
/// - Longley-Rice ITM Algorithm, Eqn 5.3, 5.9, 5.10, 5.11.
/// - ITM Technical Note 101 (TN101), Fig 10.13, Eqn III.69 & III.70.
/// - Hufford (1982): Situation variability guidance.
///
/// # Examples
///
/// ```
/// use itm::math::variability::{variability_loss, Climate, VariabilityMode};
/// use itm::math::variability::variability_warnings;
///
/// let h_e_m = [30.0, 30.0];
/// let mut warnings = 0u32;
/// let loss = variability_loss(
///     50.0,
///     50.0,
///     50.0,
///     h_e_m,
///     90.0,
///     900.0,
///     100_000.0,
///     120.0,
///     Climate::ContinentalTemperate,
///     VariabilityMode::Broadcast,
///     false,
///     false,
///     &mut warnings,
/// );
/// assert!(loss.is_finite());
/// assert_eq!(warnings & variability_warnings::EXTREME_VARIABILITIES, 0);
/// ```
pub fn variability_loss(
    time_percent: f64,
    location_percent: f64,
    situation_percent: f64,
    h_e_m: [f64; 2],
    delta_h_m: f64,
    f_mhz: f64,
    d_m: f64,
    a_ref_db: f64,
    climate: Climate,
    mode: VariabilityMode,
    disable_location_variability: bool,
    disable_situation_variability: bool,
    warnings: &mut u32,
) -> f64 {
    // Asymptotic values from TN101, Fig 10.13 (TN101v2 Eqn III.69 & III.70).
    const ALL_YEAR: [[f64; 7]; 5] = [
        [-9.67, -0.62, 1.26, -9.21, -0.62, -0.39, 3.15],
        [12.7, 9.19, 15.5, 9.05, 9.19, 2.86, 857.9],
        [
            144.9e3, 228.9e3, 262.6e3, 84.1e3, 228.9e3, 141.7e3, 2222.0e3,
        ],
        [
            190.3e3, 205.2e3, 185.2e3, 101.1e3, 205.2e3, 315.9e3, 164.8e3,
        ],
        [133.8e3, 143.6e3, 99.8e3, 98.6e3, 143.6e3, 167.4e3, 116.3e3],
    ];

    const BSM1: [f64; 7] = [2.13, 2.66, 6.11, 1.98, 2.68, 6.86, 8.51];
    const BSM2: [f64; 7] = [159.5, 7.67, 6.65, 13.11, 7.16, 10.38, 169.8];
    const XSM1: [f64; 7] = [762.2e3, 100.4e3, 138.2e3, 139.1e3, 93.7e3, 187.8e3, 609.8e3];
    const XSM2: [f64; 7] = [
        123.6e3, 172.5e3, 242.2e3, 132.7e3, 186.8e3, 169.6e3, 119.9e3,
    ];
    const XSM3: [f64; 7] = [94.5e3, 136.4e3, 178.6e3, 193.5e3, 133.5e3, 108.9e3, 106.6e3];

    const BSP1: [f64; 7] = [2.11, 6.87, 10.08, 3.68, 4.75, 8.58, 8.43];
    const BSP2: [f64; 7] = [102.3, 15.53, 9.60, 159.3, 8.12, 13.97, 8.19];
    const XSP1: [f64; 7] = [636.9e3, 138.7e3, 165.3e3, 464.4e3, 93.2e3, 216.0e3, 136.2e3];
    const XSP2: [f64; 7] = [134.8e3, 143.7e3, 225.7e3, 93.1e3, 135.9e3, 152.0e3, 188.5e3];
    const XSP3: [f64; 7] = [95.6e3, 98.6e3, 129.7e3, 94.2e3, 113.4e3, 122.7e3, 122.9e3];

    const C_D: [f64; 7] = [1.224, 0.801, 1.380, 1.000, 1.224, 1.518, 1.518];
    const Z_D: [f64; 7] = [1.282, 2.161, 1.282, 20.0, 1.282, 1.282, 1.282];

    const BFM1: [f64; 7] = [1.0, 1.0, 1.0, 1.0, 0.92, 1.0, 1.0];
    const BFM2: [f64; 7] = [0.0, 0.0, 0.0, 0.0, 0.25, 0.0, 0.0];
    const BFM3: [f64; 7] = [0.0, 0.0, 0.0, 0.0, 1.77, 0.0, 0.0];

    const BFP1: [f64; 7] = [1.0, 0.93, 1.0, 0.93, 0.93, 1.0, 1.0];
    const BFP2: [f64; 7] = [0.0, 0.31, 0.0, 0.19, 0.31, 0.0, 0.0];
    const BFP3: [f64; 7] = [0.0, 2.00, 0.0, 1.79, 2.00, 0.0, 0.0];

    let mut z_t = inverse_ccdf(time_percent / 100.0);
    let mut z_l = inverse_ccdf(location_percent / 100.0);
    let z_s = inverse_ccdf(situation_percent / 100.0);

    let climate_idx = climate.index();
    let wn = f_mhz / 47.7;

    let d_ex_m = (2.0 * A_9000_M * h_e_m[0]).sqrt()
        + (2.0 * A_9000_M * h_e_m[1]).sqrt()
        + (575.7e12 / wn).powf(THIRD); // [Algorithm, Eqn 5.3]

    let d_e_m = if d_m < d_ex_m {
        130e3 * d_m / d_ex_m
    } else {
        130e3 + d_m - d_ex_m
    };

    // Situation variability calculations.
    let sigma_s = if disable_situation_variability {
        0.0
    } else {
        let d_scale_m = 100e3;
        5.0 + 3.0 * (-d_e_m / d_scale_m).exp() // [Algorithm, Eqn 5.10]
    };

    // Apply mode adjustments.
    match mode {
        VariabilityMode::SingleMessage => {
            z_t = z_s;
            z_l = z_s;
        }
        VariabilityMode::Accidental => {
            z_l = z_s;
        }
        VariabilityMode::Mobile => {
            z_l = z_t;
        }
        VariabilityMode::Broadcast => {}
    }

    if z_t.abs() > 3.10 || z_l.abs() > 3.10 || z_s.abs() > 3.10 {
        *warnings |= variability_warnings::EXTREME_VARIABILITIES;
    }

    // Location variability calculations.
    let sigma_l = if disable_location_variability {
        0.0
    } else {
        let delta_h_d_m = terrain_roughness(d_m, delta_h_m);
        10.0 * wn * delta_h_d_m / (wn * delta_h_d_m + 13.0) // [Algorithm, Eqn 5.9 context]
    };
    let y_l = sigma_l * z_l;

    // Time variability calculations.
    let q = (0.133 * wn).ln();
    let g_minus = BFM1[climate_idx] + BFM2[climate_idx] / ((BFM3[climate_idx] * q).powi(2) + 1.0);
    let g_plus = BFP1[climate_idx] + BFP2[climate_idx] / ((BFP3[climate_idx] * q).powi(2) + 1.0);

    let sigma_t_minus = curve(
        BSM1[climate_idx],
        BSM2[climate_idx],
        XSM1[climate_idx],
        XSM2[climate_idx],
        XSM3[climate_idx],
        d_e_m,
    ) * g_minus;

    let sigma_t_plus = curve(
        BSP1[climate_idx],
        BSP2[climate_idx],
        XSP1[climate_idx],
        XSP2[climate_idx],
        XSP3[climate_idx],
        d_e_m,
    ) * g_plus;

    let sigma_td = C_D[climate_idx] * sigma_t_plus;
    let tgtd = (sigma_t_plus - sigma_td) * Z_D[climate_idx];

    let sigma_t = if z_t < 0.0 {
        sigma_t_minus
    } else if z_t <= Z_D[climate_idx] {
        sigma_t_plus
    } else {
        sigma_td + tgtd / z_t
    };
    let y_t = sigma_t * z_t;

    let y_s_temp =
        sigma_s.powi(2) + y_t.powi(2) / (7.8 + z_s.powi(2)) + y_l.powi(2) / (24.0 + z_s.powi(2)); // [Algorithm, Eqn 5.11 part]

    let (y_r, y_s) = match mode {
        VariabilityMode::SingleMessage => (
            0.0,
            (sigma_t.powi(2) + sigma_l.powi(2) + y_s_temp).sqrt() * z_s,
        ),
        VariabilityMode::Accidental => (y_t, (sigma_l.powi(2) + y_s_temp).sqrt() * z_s),
        VariabilityMode::Mobile => (
            (sigma_t.powi(2) + sigma_l.powi(2)).sqrt() * z_t,
            y_s_temp.sqrt() * z_s,
        ),
        VariabilityMode::Broadcast => (y_t + y_l, y_s_temp.sqrt() * z_s),
    };

    let mut result = a_ref_db
        - curve(
            ALL_YEAR[0][climate_idx],
            ALL_YEAR[1][climate_idx],
            ALL_YEAR[2][climate_idx],
            ALL_YEAR[3][climate_idx],
            ALL_YEAR[4][climate_idx],
            d_e_m,
        )
        - y_r
        - y_s;

    // [Algorithm, Eqn 52]
    if result < 0.0 {
        result = result * (29.0 - result) / (29.0 - 10.0 * result);
    }

    result
}

#[cfg(test)]
mod variability_tests {
    use super::*;

    #[test]
    fn test_variability_loss_finite() {
        let h_e_m = [30.0, 30.0];
        let mut warnings = 0u32;

        let loss = variability_loss(
            50.0,
            50.0,
            50.0,
            h_e_m,
            90.0,
            900.0,
            100_000.0,
            120.0,
            Climate::ContinentalTemperate,
            VariabilityMode::Broadcast,
            false,
            false,
            &mut warnings,
        );

        assert!(loss.is_finite());
    }
}
