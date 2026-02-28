/// Diffraction modeling and calculations.
///
/// This module contains functions for modeling radio wave diffraction around obstacles,
/// including Fresnel integral approximations and knife-edge diffraction calculations.

/// Computes the approximate edge diffraction loss using the Fresnel integral.
///
/// This function approximates the diffraction loss in dB using the Fresnel-Kirchhoff diffraction parameter.
/// It provides a piecewise approximation of the Fresnel integral, which is commonly used in radio propagation
/// models to estimate signal attenuation due to obstacles.
///
/// # Arguments
///
/// * `v2` - The squared Fresnel-Kirchhoff diffraction parameter (v²). This parameter characterizes
///   the degree of diffraction and is typically derived from geometric and frequency considerations
///   in radio path analysis.
///
/// # Returns
///
/// The approximate diffraction loss in dB. Two different approximations are used depending on the value of v²:
/// - For v² < 5.76: A quadratic approximation
/// - For v² ≥ 5.76: A logarithmic approximation
///
/// # References
///
/// - Hufford, G.A. (1952). "An Integral Equation Approach to the Problem of Wave Propagation over an Irregular Surface"
/// - ITM Technical Note 101 (TN101), Section I: Fresnel integral approximations
/// - Deygout, J. (1966). "Multiple knife-edge diffraction of microwaves". IEEE Transactions on Antennas and Propagation, 14(4), 480-489.
///
/// # Examples
///
/// ```
/// use itm_rs::math::fresnel_integral;
///
/// // Small diffraction parameter
/// let loss1 = fresnel_integral(2.0);
/// assert!(loss1 > 0.0);
///
/// // Large diffraction parameter
/// let loss2 = fresnel_integral(10.0);
/// assert!(loss2 > loss1);
/// ```
pub fn fresnel_integral(v2: f64) -> f64 {
    if v2 < 5.76 {
        6.02 + 9.11 * v2.sqrt() - 1.27 * v2
    } else {
        12.953 + 10.0 * v2.log10()
    }
}

/// Computes the knife-edge diffraction loss over an obstruction.
///
/// This function calculates the diffraction loss when radio signals encounter a single
/// knife-edge obstruction (or an obstruction that can be approximated as one). It uses
/// the Fresnel diffraction theory to determine the signal attenuation caused by the
/// obstruction, which is particularly useful in terrain shadowing scenarios.
///
/// # Arguments
///
/// * `d` - Total path distance in meters. The distance from transmitter to receiver.
/// * `f` - Frequency of the signal in Hz. Must be positive.
/// * `a_e` - Earth's effective radius in meters, accounting for atmospheric refraction.
///   Typical value is around 8.5e6 meters (8500 km).
/// * `theta_los` - Line-of-sight grazing angle in radians. The angle between the
///   direct path and the Earth's surface.
/// * `d_hzn` - Array of two horizon distances in meters:
///   - `d_hzn[0]` = Distance from transmitter to horizon/obstruction point
///   - `d_hzn[1]` = Distance from obstruction point to receiver
///
/// # Returns
///
/// The knife-edge diffraction loss in dB. Represents the additional path loss due to
/// the obstruction compared to free space propagation.
///
/// # Algorithm Notes
///
/// This implementation uses the Fresnel-Kirchhoff diffraction model with the following approach:
/// - Calculates the Fresnel diffraction parameter (v) for both sides of the obstruction
/// - Accounts for asymmetric geometry where the obstruction may not be at the midpoint
/// - Uses pre-computed curve fits via `fresnel_integral()` for efficiency
/// - Follows the ITM (Irregular Terrain Model) implementation standards
///
/// # References
///
/// - Longley-Rice ITM Algorithm, Equation 4.12: Angular distance calculation
/// - ITM Technical Note 101 (TN101), Equation I.7: Fresnel parameter calculation (1/(4π) ≈ 0.0795775)
/// - ITM Technical Note 101 (TN101), Equation I.1: Total diffraction loss
/// - Bullington, K. (1947). "Radio Propagation at Frequencies above 30 Megacycles". Proceedings of the IRE, 35(10), 1122-1136.
pub(crate) fn knife_edge_diffraction(
    d: f64,
    f: f64,
    a_e: f64,
    theta_los: f64,
    d_hzn: [f64; 2],
) -> f64 {
    use std::f64;

    // Total line-of-sight distance coverage up to the obstruction
    let los_distance = d_hzn[0] + d_hzn[1];

    // Angular distance beyond the direct line-of-sight in the diffraction region
    let diffraction_angle = d / a_e - theta_los;

    // Distance beyond line-of-sight where diffraction occurs
    let beyond_los_distance = d - los_distance;

    // Fresnel diffraction parameter for receiver side
    // 1 / (4π) = 0.0795775 [TN101, Eqn I.7]
    const RECIP_FOUR_PI: f64 = 1.0 / (4.0 * std::f64::consts::PI);
    let diffraction_angle_squared = diffraction_angle.powi(2);
    let fresnel_param_receiver =
        RECIP_FOUR_PI * (f / 47.7) * diffraction_angle_squared * d_hzn[1] * beyond_los_distance
            / (beyond_los_distance + d_hzn[1]);

    // Fresnel diffraction parameter for transmitter side
    let fresnel_param_transmitter =
        RECIP_FOUR_PI * (f / 47.7) * diffraction_angle_squared * d_hzn[0] * beyond_los_distance
            / (beyond_los_distance + d_hzn[0]);

    // Total diffraction loss is sum of contributions from both sides [TN101, Eqn I.1]
    fresnel_integral(fresnel_param_transmitter) + fresnel_integral(fresnel_param_receiver)
}

use crate::math::constants::{EARTH_RADIUS_M, THIRD};
use num_complex::Complex;

/// Computes the smooth-earth diffraction loss using the Vogler 3-radii method.
///
/// This function implements the Vogler (1964) smooth-earth diffraction model with
/// three effective radii, combining a distance function with two height gain terms.
///
/// # Arguments
///
/// * `d_m` - Path distance in meters.
/// * `f_mhz` - Frequency in MHz.
/// * `a_e_m` - Effective earth radius in meters.
/// * `theta_los` - Angular distance of the line-of-sight region in radians.
/// * `d_hzn_m` - Terminal horizon distances in meters `[tx, rx]`.
/// * `h_e_m` - Effective terminal heights in meters `[tx, rx]`.
/// * `z_g` - Complex ground impedance.
///
/// # Returns
///
/// Smooth-earth diffraction loss in dB.
///
/// # References
///
/// - Longley-Rice ITM Algorithm, Eqn 4.12 and 4.20.
/// - ITM Technical Note 101 (TN101), Eqn 8.4.
/// - Vogler, L. (1964). "An Attenuation Function for Propagation over Irregular Terrain".
///
/// # Examples
///
/// ```
/// use itm_rs::math::smooth_earth_diffraction;
/// use num_complex::Complex;
///
/// let d_hzn_m = [50_000.0, 50_000.0];
/// let h_e_m = [30.0, 30.0];
/// let z_g = Complex::new(15.0, 0.1);
/// let loss = smooth_earth_diffraction(100_000.0, 900.0, 8.5e6, 0.001, d_hzn_m, h_e_m, z_g);
/// assert!(loss.is_finite());
/// ```
pub fn smooth_earth_diffraction(
    d_m: f64,
    f_mhz: f64,
    a_e_m: f64,
    theta_los: f64,
    d_hzn_m: [f64; 2],
    h_e_m: [f64; 2],
    z_g: Complex<f64>,
) -> f64 {
    let mut a_m = [0.0; 3];
    let mut d_km = [0.0; 3];
    let mut f_x_db = [0.0; 2];
    let mut k = [0.0; 3];
    let mut b_0 = [0.0; 3];
    let mut x_km = [0.0; 3];
    let mut c_0 = [0.0; 3];

    let theta_nlos = d_m / a_e_m - theta_los; // [Algorithm, Eqn 4.12]
    let d_ml_m = d_hzn_m[0] + d_hzn_m[1];

    // Compute 3 radii.
    a_m[0] = (d_m - d_ml_m) / (d_m / a_e_m - theta_los);
    a_m[1] = 0.5 * d_hzn_m[0].powi(2) / h_e_m[0];
    a_m[2] = 0.5 * d_hzn_m[1].powi(2) / h_e_m[1];

    d_km[0] = (a_m[0] * theta_nlos) / 1000.0;
    d_km[1] = d_hzn_m[0] / 1000.0;
    d_km[2] = d_hzn_m[1] / 1000.0;

    for i in 0..3 {
        // C_0 = (4/(3k))^(1/3) with k = a_i / a_0. [Vogler 1964, Eqn 2]
        c_0[i] = ((4.0 / 3.0) * EARTH_RADIUS_M / a_m[i]).powf(THIRD);

        // [Vogler 1964, Eqn 6a/7a]
        k[i] = 0.017778 * c_0[i] * f_mhz.powf(-THIRD) / z_g.norm();

        // [Vogler 1964, Fig 4]
        b_0[i] = 1.607 - k[i];
    }

    // Compute x_km for each radius. [Vogler 1964, Eqn 2]
    x_km[1] = b_0[1] * c_0[1].powi(2) * f_mhz.powf(THIRD) * d_km[1];
    x_km[2] = b_0[2] * c_0[2].powi(2) * f_mhz.powf(THIRD) * d_km[2];
    x_km[0] = b_0[0] * c_0[0].powi(2) * f_mhz.powf(THIRD) * d_km[0] + x_km[1] + x_km[2];

    // Height gain functions.
    f_x_db[0] = height_function(x_km[1], k[1]);
    f_x_db[1] = height_function(x_km[2], k[2]);

    // Distance function. [TN101, Eqn 8.4] & [Vogler 1964, Eqn 13]
    let g_x_db = 0.05751 * x_km[0] - 10.0 * x_km[0].log10();

    g_x_db - f_x_db[0] - f_x_db[1] - 20.0 // [Algorithm, Eqn 4.20]
}

/// Height function F(x, K) for smooth-earth diffraction.
///
/// # Arguments
///
/// * `x_km` - Normalized distance in kilometers.
/// * `k` - Vogler K value.
///
/// # Returns
///
/// Height function value in dB.
///
/// # References
///
/// - Vogler, L. (1964). Eqn 6a/7a and Fig 4.
/// - ITM Technical Note 101 (TN101), Eqn 8.4.
fn height_function(x_km: f64, k: f64) -> f64 {
    let mut result;

    if x_km < 200.0 {
        let w = -k.ln();

        if k < 1e-5 || x_km * w.powi(3) > 5495.0 {
            result = -117.0;

            if x_km > 1.0 {
                result = 17.372 * x_km.ln() + result;
            }
        } else {
            result = 2.5e-5 * x_km.powi(2) / k - 8.686 * w - 15.0;
        }
    } else {
        result = 0.05751 * x_km - 4.343 * x_km.ln();

        if x_km < 2000.0 {
            let w = 0.0134 * x_km * (-0.005 * x_km).exp();
            result = (1.0 - w) * result + w * (17.372 * x_km.ln() - 117.0);
        }
    }

    result
}

#[cfg(test)]
mod smooth_earth_diffraction_tests {
    use super::smooth_earth_diffraction;
    use num_complex::Complex;

    #[test]
    fn test_smooth_earth_diffraction_finite() {
        let d_hzn_m = [50_000.0, 50_000.0];
        let h_e_m = [30.0, 30.0];
        let z_g = Complex::new(15.0, 0.1);

        let loss = smooth_earth_diffraction(100_000.0, 900.0, 8.5e6, 0.001, d_hzn_m, h_e_m, z_g);
        assert!(loss.is_finite());
    }
}
