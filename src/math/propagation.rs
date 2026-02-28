use crate::math::constants::GAMMA_A;
/// Radio propagation models and calculations.
///
/// This module contains functions for calculating basic radio propagation parameters,
/// including free space path loss and related transmission characteristics.
use num_complex::Complex;

/// Computes the free space loss (path loss) between two points in free space.
///
/// This function calculates the power attenuation that occurs when electromagnetic waves
/// propagate in a straight line through free space with no obstructions, reflections, or
/// scattering. It is based on the Friis transmission equation, a fundamental formula in
/// radio propagation modeling.
///
/// # Arguments
///
/// * `d` - Distance between the two antennas in meters. Must be positive.
/// * `f` - Frequency of the signal in Hz. Must be positive.
///
/// # Returns
///
/// The free space loss in dB (decibels). Higher values indicate greater signal attenuation.
///
/// # Formula
///
/// The Friis transmission equation used is:
///
/// ```text
/// L = 20 * log10(d) + 20 * log10(f) + 32.45
/// ```
///
/// where:
/// - L = free space loss in dB
/// - d = distance between antennas in meters (converted to km internally)
/// - f = frequency of the signal in Hz
/// - 32.45 = constant derived from the Friis equation constants
///
/// # References
///
/// - Friis, H.T. (1946). "A Note on a Simple Transmission Formula". Proceedings of the IRE, 34(5), 254-256.
/// - ITU-R P.525: "Calculation of free-space attenuation"
///
/// # Examples
///
/// ```
/// use itm::math::free_space_loss;
///
/// // 1 km distance at 1 GHz
/// let loss1 = free_space_loss(1000.0, 1e9);
/// assert!(loss1 > 90.0); // Approximately 92 dB
///
/// // Shorter distance results in less loss
/// let loss2 = free_space_loss(100.0, 1e9);
/// assert!(loss2 < loss1);
/// ```
pub fn free_space_loss(d: f64, f: f64) -> f64 {
    20.0 * (d / 1000.0).log10() + 20.0 * (f).log10() + 32.45
}

/// Polarization of the radio wave.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Polarization {
    /// Horizontal polarization.
    Horizontal,
    /// Vertical polarization.
    Vertical,
}

/// Initializes point-to-point propagation parameters.
///
/// This function computes ground impedance, effective earth curvature, and surface
/// refractivity for point-to-point mode, based on frequency, path height, and ground
/// electrical properties. The implementation follows the ITM (Longley-Rice) model.
///
/// # Arguments
///
/// * `f_mhz` - Frequency in MHz.
/// * `h_sys_m` - Average height of the path above mean sea level in meters.
/// * `n_0` - Refractivity at sea level in N-Units.
/// * `polarization` - Wave polarization (horizontal or vertical).
/// * `epsilon` - Relative permittivity of the ground.
/// * `sigma` - Ground conductivity in S/m.
///
/// # Returns
///
/// A tuple `(z_g, gamma_e, n_s)` where:
/// - `z_g` is the complex ground impedance.
/// - `gamma_e` is the curvature of the effective earth (1/m).
/// - `n_s` is the surface refractivity in N-Units.
///
/// # References
///
/// - ITM Technical Note 101 (TN101), Eq 4.3: Refractivity scaling with height.
/// - ITM Technical Note 101 (TN101), Eq 4.4: Effective earth curvature.
/// - Longley-Rice ITM Algorithm: Ground impedance formulation.
///
/// # Examples
///
/// ```
/// use itm::math::{initialize_point_to_point, Polarization};
///
/// let (z_g, gamma_e, n_s) = initialize_point_to_point(
///     900.0,
///     100.0,
///     301.0,
///     Polarization::Horizontal,
///     15.0,
///     0.005,
/// );
/// assert!(z_g.re.is_finite());
/// assert!(gamma_e > 0.0);
/// assert!(n_s > 0.0);
/// ```
pub fn initialize_point_to_point(
    f_mhz: f64,
    h_sys_m: f64,
    n_0: f64,
    polarization: Polarization,
    epsilon: f64,
    sigma: f64,
) -> (Complex<f64>, f64, f64) {
    // Curvature of the actual earth (imported from `crate::math::constants::GAMMA_A`).

    // Scale refractivity based on elevation above mean sea level.
    let n_s = if h_sys_m == 0.0 {
        n_0
    } else {
        n_0 * (-h_sys_m / 9460.0).exp() // [TN101, Eq 4.3]
    };

    // Curvature of the effective earth. [TN101, Eq 4.4]
    let gamma_e = GAMMA_A * (1.0 - 0.04665 * (n_s / 179.3).exp());

    // Complex relative permittivity.
    let ep_r = Complex::new(epsilon, 18000.0 * sigma / f_mhz);

    // Ground impedance (horizontal polarization).
    let mut z_g = (ep_r - Complex::new(1.0, 0.0)).sqrt();

    // Adjust for vertical polarization.
    if polarization == Polarization::Vertical {
        z_g /= ep_r;
    }

    (z_g, gamma_e, n_s)
}

/// Line-of-sight loss wrapper.
///
/// For now this function uses the diffraction approximation as a proxy for LOS
/// loss computation at short distances.
///
/// # Arguments
///
/// * `d` - Total path distance.
/// * `h_e` - Effective heights [tx, rx].
/// * `z_g` - Complex ground impedance.
/// * `delta_h` - Terrain roughness parameter.
/// * `m_d` - Slope of diffraction loss.
/// * `a_d0` - Intercept of diffraction loss.
/// * `d_sml` - Smooth earth horizon distance sum.
/// * `f` - Frequency in MHz.
///
/// # Returns
///
/// The line-of-sight loss in dB.
pub fn line_of_sight_loss(
    d: f64,
    h_e: [f64; 2],
    z_g: Complex<f64>,
    delta_h: f64,
    m_d: f64,
    a_d0: f64,
    d_sml: f64,
    f: f64,
) -> f64 {
    // Implemented to match the C++ LineOfSightLoss algorithm.
    // delta_h_d = TerrainRoughness(d, delta_h)
    let delta_h_d = crate::math::terrain::terrain_roughness(d, delta_h);

    // sigma_h_d = SigmaHFunction(delta_h_d)
    let sigma_h_d = crate::math::terrain::sigma_h_function(delta_h_d);

    // wavenumber (wn = f / 47.7)
    let wn = f / crate::math::constants::WAVENUMBER_DIVISOR;

    // sin_psi per Algorithm Eqn 4.46
    let sum_h = h_e[0] + h_e[1];
    let sin_psi = sum_h / d.hypot(sum_h);

    // R_e = (sin_psi - Z_g) / (sin_psi + Z_g) * exp(-MIN(10.0, wn * sigma_h_d * sin_psi))
    let s_c = Complex::new(sin_psi, 0.0);
    let mut r_e = (s_c - z_g) / (s_c + z_g);
    // MIN(10, wn * sigma_h_d * sin_psi)
    let min_val = (wn * sigma_h_d * sin_psi).min(10.0);
    r_e *= (-min_val).exp();

    // q = |R_e|^2
    let mut q = r_e.re.powi(2) + r_e.im.powi(2);
    if q <= 0.0 {
        q = std::f64::EPSILON;
    }
    if q < 0.25 || q < sin_psi {
        let scale = (sin_psi / q).sqrt();
        r_e *= scale;
        q = r_e.re.powi(2) + r_e.im.powi(2);
    }

    // phase difference delta_phi = wn * 2 * h_e[0] * h_e[1] / d
    let mut delta_phi = wn * 2.0 * h_e[0] * h_e[1] / d;
    if delta_phi > std::f64::consts::PI / 2.0 {
        delta_phi = std::f64::consts::PI - (std::f64::consts::PI / 2.0).powi(2) / delta_phi;
    }

    // two-ray attenuation
    let rr = Complex::new(delta_phi.cos(), -delta_phi.sin()) + r_e;
    let a_t_db = -10.0 * (rr.re.powi(2) + rr.im.powi(2)).log10();

    // extended diffraction attenuation
    let a_d_db = m_d * d + a_d0;

    // weighting factor w = 1 / (1 + f * delta_h / max(10e3, d_sml_m))
    let denom = d_sml.max(10e3);
    let w = 1.0 / (1.0 + f * delta_h / denom);

    // final LOS loss
    w * a_t_db + (1.0 - w) * a_d_db
}
