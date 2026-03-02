#![allow(clippy::too_many_arguments)]
//!! Tropospheric scatter propagation calculations.
//!!
//!! This module contains functions for modeling tropospheric scatter propagation,
//!! a mode where radio waves are scattered by atmospheric turbulence and
//!! irregularities in the troposphere.
use crate::math::constants::{
    HEIGHT_CONSTANT_H, SMOOTHING_DISTANCE_D0, TROPOSCATTER_UNDEFINED_SENTINEL, WAVENUMBER_DIVISOR,
};
use crate::math::height_gain::h0_function;
use std::f64::consts::SQRT_2;

/// Computes the attenuation function F(θ·d) for scatter propagation.
///
/// This function calculates the attenuation factor used in troposcatter loss predictions.
/// It uses different empirical curve fits depending on the path distance, accounting for
/// the fact that atmospheric scatter characteristics change with distance.
///
/// # Arguments
///
/// * `theta_distance` - Product of angular distance (theta) and path distance (d), in meter-radians.
///   This represents the angular path length through the scattering region.
///
/// # Returns
///
/// The attenuation function value F(θ·d) in dB. This value is used as a component
/// in the total troposcatter loss calculation.
///
/// # Algorithm Notes
///
/// The function uses three different sets of empirical parameters (a, b, c) based on
/// the distance range:
/// - Short paths (≤ 10 km): Parameters optimized for near-field scatter
/// - Medium paths (10-70 km): Transition region parameters
/// - Long paths (> 70 km): Parameters for far-field scatter
///
/// The formula used is: `F = a + b·(θ·d) + c·log₁₀(θ·d)`
///
/// # References
///
/// - Longley-Rice ITM Algorithm, Equation 6.9: Attenuation function parameters
/// - ITM Technical Note 101 (TN101): Troposcatter attenuation modeling
///
/// # Examples
///
/// ```
/// use itm::math::scatter::f_function;
///
/// // Calculate attenuation for a 5 km path
/// let f = f_function(5_000.0);
/// assert!(f > 0.0);
///
/// // Longer paths have different attenuation characteristics
/// let f_long = f_function(100_000.0);
/// ```
pub fn f_function(theta_distance: f64) -> f64 {
    // Constants from [Algorithm, 6.9]
    const A: [f64; 3] = [133.4, 104.6, 71.8];
    const B: [f64; 3] = [0.332e-3, 0.212e-3, 0.157e-3];
    const C: [f64; 3] = [-10.0, -2.5, 5.0];

    // Select the set of values based on distance
    let param_index = if theta_distance <= 10e3 {
        0 // Short paths: <= 10 km
    } else if theta_distance <= 70e3 {
        1 // Medium paths: 10 km to 70 km
    } else {
        2 // Long paths: > 70 km
    };

    // Calculate F_0 using the selected parameters [Algorithm, 6.9]
    A[param_index] + B[param_index] * theta_distance + C[param_index] * theta_distance.log10()
}

/// Computes the tropospheric scatter propagation loss.
///
/// This function calculates the signal loss due to tropospheric scatter, a propagation
/// mechanism where radio waves are scattered by atmospheric turbulence and irregularities.
/// This mode is particularly important for beyond-horizon communications and long-distance
/// radio links.
///
/// Troposcatter involves the scattering of radio waves by small-scale irregularities in
/// the refractive index of the troposphere, caused by turbulence and mixing of air masses
/// with different temperatures and humidity levels.
///
/// # Arguments
///
/// * `d` - Total path distance in meters from transmitter to receiver
/// * `theta_hzn` - Array of terminal horizon angles in radians [transmitter, receiver]
/// * `d_hzn` - Array of terminal horizon distances in meters [transmitter, receiver]
/// * `h_e` - Array of effective terminal heights in meters [transmitter, receiver],
///   accounting for atmospheric refraction effects
/// * `a_e` - Effective earth radius in meters, typically around 8.5e6 meters (accounts
///   for atmospheric refraction, approximately 4/3 of actual earth radius)
/// * `n_s` - Surface refractivity in N-Units (typically 250-400), which characterizes
///   the atmospheric refractive properties at the surface
/// * `f` - Frequency in MHz
/// * `theta_los` - Angular distance of the line-of-sight region in radians
/// * `h0` - Mutable reference to H_0 value (height gain factor), updated by this function
///
/// # Returns
///
/// The troposcatter loss in dB. Returns 1001.0 as a sentinel value if the calculation
/// is undefined (when both r_1 and r_2 are less than 0.2, indicating the scatter
/// function is not defined or infinite).
///
/// # Algorithm Notes
///
/// The calculation involves several key steps:
/// 1. Compute wavenumber from frequency
/// 2. Calculate scattering parameters (r_1, r_2) based on terminal heights
/// 3. Determine asymmetry parameter and cross-over height
/// 4. Calculate scattering efficiency factor (eta_s)
/// 5. Compute height gain factor H_0 with interpolation for low eta_s
/// 6. Combine components: F-function, frequency factor, refractivity correction, and H_0
///
/// Special handling:
/// - If H_0 > 15 dB, uses cached value to maintain continuity
/// - For eta_s < 1, interpolates between general case and special case (eta_s = 0)
/// - Constrains asymmetry parameters (s, q) to range [0.1, 10]
///
/// # References
///
/// - Longley-Rice ITM Algorithm, Equations 4.63, 4.66, 4.67, 6.8: Core troposcatter formulas
/// - ITM Technical Note 101 (TN101), Equations 9.3a, 9.3b, 9.4a, 9.5: Detailed derivations
/// - Rice, P.L., et al. "Transmission Loss Predictions for Tropospheric Communication Circuits"
///
/// # Examples
///
/// ```
/// use itm::math::scatter::troposcatter_loss;
///
/// let d = 100_000.0; // 100 km
/// let theta_hzn = [0.01, 0.01]; // horizon angles in radians
/// let d_hzn = [50_000.0, 50_000.0]; // horizon distances
/// let h_e = [100.0, 100.0]; // effective heights in meters
/// let a_e = 8.5e6; // effective earth radius
/// let n_s = 301.0; // surface refractivity
/// let f = 1000.0; // 1 GHz (1000 MHz)
/// let theta_los = 0.001;
/// let mut h0 = 0.0;
///
/// let loss = troposcatter_loss(d, &theta_hzn, &d_hzn, &h_e, a_e, n_s, f, theta_los, &mut h0);
/// assert!(loss > 0.0 && loss < 1000.0); // Valid loss value
/// ```
pub fn troposcatter_loss(
    d: f64,
    theta_hzn: &[f64; 2],
    d_hzn: &[f64; 2],
    h_e: &[f64; 2],
    a_e: f64,
    n_s: f64,
    f: f64,
    theta_los: f64,
    h0: &mut f64,
) -> f64 {
    // Wavenumber k = f / WAVENUMBER_DIVISOR
    let wavenumber = f / WAVENUMBER_DIVISOR;

    let height_gain = if *h0 > 15.0 {
        // Short-circuit calculations if already greater than 15 dB
        *h0
    } else {
        // Calculate asymmetry in horizon distances
        let mut distance_asymmetry = d_hzn[0] - d_hzn[1];
        let mut height_ratio = h_e[1] / h_e[0];

        // Ensure correct frame of reference
        if distance_asymmetry < 0.0 {
            distance_asymmetry = -distance_asymmetry;
            height_ratio = 1.0 / height_ratio;
        }

        // Total angular distance in radians
        let angular_distance = theta_hzn[0] + theta_hzn[1] + d / a_e;

        // Calculate scattering parameters [TN101, Eqn 9.4a]
        let r_1 = 2.0 * wavenumber * angular_distance * h_e[0];
        let r_2 = 2.0 * wavenumber * angular_distance * h_e[1];

        // Check if scatter function is defined
        if r_1 < 0.2 && r_2 < 0.2 {
            // "If both r_1 and r_2 are less than 0.2 the function A_scat is not defined
            // (or is infinite)" [Algorithm, page 11]
            return TROPOSCATTER_UNDEFINED_SENTINEL;
        }

        // Asymmetry parameter [TN101, Eqn 9.5]
        let mut asymmetry_s = (d - distance_asymmetry) / (d + distance_asymmetry);

        // "In all of this, we truncate the values of s and q at 0.1 and 10"
        // [Algorithm, page 16]
        let asymmetry_q = (height_ratio / asymmetry_s).clamp(0.1, 10.0);
        asymmetry_s = asymmetry_s.max(0.1);

        // Height of cross-over [Algorithm, 4.66] [TN101v1, 9.3b]
        let crossover_height =
            (d - distance_asymmetry) * (d + distance_asymmetry) * angular_distance * 0.25 / d;

        // Scale heights [Algorithm, 4.67]
        const SCALE_HEIGHT_Z0: f64 = 1.7556e3; // meters
        const SCALE_HEIGHT_Z1: f64 = 8.0e3; // meters

        // Scattering efficiency factor eta_s [TN101 Eqn 9.3a]
        let height_ratio_z1 = (crossover_height / SCALE_HEIGHT_Z1).min(1.7);
        let refractivity_factor = 0.031 - n_s * 2.32e-3 + n_s.powi(2) * 5.67e-6;
        let eta_s = (crossover_height / SCALE_HEIGHT_Z0)
            * (1.0 + refractivity_factor * (-height_ratio_z1.powi(6)).exp());

        // First term in TN101v1, Eqn 9.5
        let h_00 = (h0_function(r_1, eta_s) + h0_function(r_2, eta_s)) / 2.0;

        // Second term correction
        let delta_h0 =
            (6.0 * (0.6 - eta_s.max(1.0).log10()) * asymmetry_s.log10() * asymmetry_q.log10())
                .min(h_00);

        let mut h_0 = h_00 + delta_h0;

        // "If Delta_H_0 would make H_0 negative, use H_0 = 0" [TN101v1, p9.4]
        h_0 = h_0.max(0.0);

        // If eta_s <= 1, interpolate with the special case of eta_s = 0
        if eta_s < 1.0 {
            let sqrt2_term_1 = 1.0 + SQRT_2 / r_1;
            let sqrt2_term_2 = 1.0 + SQRT_2 / r_2;
            let special_case = 10.0
                * ((sqrt2_term_1 * sqrt2_term_2).powi(2) * (r_1 + r_2)
                    / (r_1 + r_2 + 2.0 * SQRT_2))
                    .log10();

            h_0 = eta_s * h_0 + (1.0 - eta_s) * special_case;
        }

        // "If, at d_5, calculations show that H_0 will exceed 15 dB, they are replaced
        // by the value it has at d_6" [Algorithm, page 12]
        if h_0 > 15.0 && *h0 >= 0.0 {
            h_0 = *h0;
        }

        h_0
    };

    // Update the H_0 value
    *h0 = height_gain;

    // Angular distance beyond line-of-sight
    let theta = d / a_e - theta_los;

    // Constants are provided by the centralized `crate::math::constants` module.
    // (SMOOTHING_DISTANCE_D0 and HEIGHT_CONSTANT_H are imported at top of file.)

    // Calculate total troposcatter loss [Algorithm, 4.63]
    f_function(theta * d) + 10.0 * (wavenumber * HEIGHT_CONSTANT_H * theta.powi(4)).log10()
        - 0.1 * (n_s - 301.0) * (-(theta * d) / SMOOTHING_DISTANCE_D0).exp()
        + height_gain
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_f_function_short_distance() {
        let f = f_function(5_000.0);
        assert!(f > 0.0, "F-function should be positive");
    }

    #[test]
    fn test_f_function_medium_distance() {
        let f = f_function(50_000.0);
        assert!(f > 0.0, "F-function should be positive");
    }

    #[test]
    fn test_f_function_long_distance() {
        let f = f_function(100_000.0);
        assert!(f > 0.0, "F-function should be positive");
    }

    #[test]
    fn test_f_function_transitions() {
        // Test that transitions between regions are reasonable
        let f1 = f_function(9_999.0);
        let f2 = f_function(10_001.0);
        let diff = (f1 - f2).abs();
        assert!(diff < 50.0, "Transition should be relatively smooth");

        let f3 = f_function(69_999.0);
        let f4 = f_function(70_001.0);
        let diff2 = (f3 - f4).abs();
        assert!(diff2 < 50.0, "Transition should be relatively smooth");
    }

    #[test]
    fn test_troposcatter_loss_basic() {
        let d = 100_000.0;
        let theta_hzn = [0.01, 0.01];
        let d_hzn = [50_000.0, 50_000.0];
        let h_e = [100.0, 100.0];
        let a_e = 8.5e6;
        let n_s = 301.0;
        let f = 1000.0;
        let theta_los = 0.001;
        let mut h0 = 0.0;

        let loss = troposcatter_loss(d, &theta_hzn, &d_hzn, &h_e, a_e, n_s, f, theta_los, &mut h0);

        assert!(loss > 0.0 && loss < 1000.0, "Loss should be in valid range");
        assert!(h0 >= 0.0, "H0 should be non-negative");
    }

    #[test]
    fn test_troposcatter_loss_undefined_case() {
        // Test case where r_1 and r_2 are too small
        let d = 100_000.0;
        let theta_hzn = [0.00001, 0.00001]; // Very small angles
        let d_hzn = [50_000.0, 50_000.0];
        let h_e = [0.1, 0.1]; // Very small heights
        let a_e = 8.5e6;
        let n_s = 301.0;
        let f = 100.0; // Low frequency
        let theta_los = 0.001;
        let mut h0 = 0.0;

        let loss = troposcatter_loss(d, &theta_hzn, &d_hzn, &h_e, a_e, n_s, f, theta_los, &mut h0);

        // Should return sentinel value
        assert_eq!(loss, 1001.0, "Should return undefined sentinel value");
    }

    #[test]
    fn test_troposcatter_loss_h0_caching() {
        let d = 100_000.0;
        let theta_hzn = [0.01, 0.01];
        let d_hzn = [50_000.0, 50_000.0];
        let h_e = [100.0, 100.0];
        let a_e = 8.5e6;
        let n_s = 301.0;
        let f = 1000.0;
        let theta_los = 0.001;
        let mut h0 = 20.0; // Already > 15

        let loss = troposcatter_loss(d, &theta_hzn, &d_hzn, &h_e, a_e, n_s, f, theta_los, &mut h0);

        // Should use cached h0 value
        assert_eq!(h0, 20.0, "H0 should remain unchanged when > 15");
        assert!(loss > 0.0 && loss < 1000.0, "Loss should be in valid range");
    }
}
