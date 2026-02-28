/// Height gain factor calculations.
///
/// This module contains functions for computing height gain factors used in
/// radio propagation models, including curve fitting and interpolation for
/// different terrain roughness categories.

/// Curve fit helper function for h0 calculations.
///
/// This function computes a component of the h0 value using a curve fit model for
/// different terrain roughness indices. It is used internally by `h0_function()` to
/// calculate height gain factors in radio propagation analysis.
///
/// # Arguments
///
/// * `j` - Index into the curve fitting parameters (0-4). Typically corresponds to
///   a roughness or terrain category index.
/// * `r` - A normalized parameter (often related to ground reflection characteristics).
///   Must be positive.
///
/// # Returns
///
/// The computed h0 curve value in dB, calculated using pre-defined coefficients A and B
/// that are indexed by parameter `j`.
///
/// # Panics
///
/// Panics if `j` is out of bounds (j > 4).
///
/// # References
///
/// - Longley, A.G. & Rice, P.L. (1968). "Prediction of Tropospheric Radio Transmission Loss Over Irregular Terrain: A Computer Method-1968". ESSA Technical Report ERL 79-ITS 67.
/// - ITM Technical Note 101 (TN101): Height gain curve fit parameters
///
/// # Examples
///
/// ```
/// use itm_rs::math::h0_curve;
///
/// let value = h0_curve(0, 5.0);
/// assert!(value > 0.0);
/// ```
pub fn h0_curve(j: usize, r: f64) -> f64 {
    const A: [f64; 5] = [25.0, 80.0, 177.0, 395.0, 705.0];
    const B: [f64; 5] = [24.0, 45.0, 68.0, 80.0, 105.0];

    10.0 * (1.0 + A[j] * (r.recip()).powi(4) + B[j] * (r.recip()).powi(2)).log10()
}

/// Computes the h0 height gain factor with smooth interpolation.
///
/// This function calculates the h0 value, which represents a height gain factor used in
/// radio propagation models. It uses interpolation between pre-computed curve fit values
/// for different terrain roughness indices, allowing for smooth transitions between
/// discrete roughness categories.
///
/// # Arguments
///
/// * `r` - A normalized parameter related to ground reflection characteristics. Must be positive.
/// * `eta_s` - The terrain roughness index (typically between 1.0 and 5.0).
///   Values outside this range are clamped to [1.0, 5.0].
///
/// # Returns
///
/// The interpolated h0 value in dB. The function uses linear interpolation between
/// adjacent curve values when `eta_s` is not an integer.
///
/// # Behavior
///
/// - If `eta_s` is an integer between 1 and 5, the value from the corresponding
///   `h0_curve()` is returned directly.
/// - If `eta_s` is fractional, linear interpolation is performed between the two
///   nearest integer curve indices.
/// - Values of `eta_s` outside [1.0, 5.0] are clamped before processing.
///
/// # References
///
/// - Longley, A.G. & Rice, P.L. (1968). "Prediction of Tropospheric Radio Transmission Loss Over Irregular Terrain: A Computer Method-1968". ESSA Technical Report ERL 79-ITS 67.
/// - ITM Technical Note 101 (TN101): Height gain interpolation method for terrain roughness
///
/// # Examples
///
/// ```
/// use itm_rs::math::h0_function;
///
/// let value = h0_function(10.0, 2.5);
/// assert!(value > 0.0);
///
/// // Clamping example
/// let clamped_high = h0_function(10.0, 10.0);
/// let clamped_normal = h0_function(10.0, 5.0);
/// assert_eq!(clamped_high, clamped_normal);
/// ```
pub fn h0_function(r: f64, eta_s: f64) -> f64 {
    let eta_s = eta_s.clamp(1.0, 5.0);
    let i = eta_s as usize;
    let q = eta_s - i as f64;

    let mut result = h0_curve(i - 1, r);

    if q != 0.0 {
        result = (1.0 - q) * result + q * h0_curve(i, r);
    }
    result
}
