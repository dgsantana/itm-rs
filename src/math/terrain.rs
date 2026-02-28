/// Terrain profile analysis and calculations.
///
/// This module contains functions for analyzing terrain profiles, including
/// statistical fitting methods and terrain characterization for radio propagation models.
use std::f64::consts::PI;

/// Computes the terrain roughness parameter (interdecile range) for radio propagation.
///
/// This function calculates the terrain roughness, which quantifies terrain irregularity
/// over a given distance. The roughness parameter is essential for determining how terrain
/// variations affect radio wave propagation, particularly for scatter and diffraction modes.
///
/// The calculation uses an exponential smoothing function that accounts for the fact that
/// terrain features become less significant for very long paths, effectively modeling the
/// averaging effect of distance on perceived terrain roughness.
///
/// # Arguments
///
/// * `d` - Distance over which to calculate roughness in meters. Must be positive.
/// * `delta_h` - Interdecile range of terrain elevations in meters (difference between
///   the 90th and 10th percentile elevations along the path). This represents the
///   statistical measure of terrain height variation.
///
/// # Returns
///
/// The effective terrain roughness parameter in meters, accounting for distance-dependent
/// smoothing effects. The value is always non-negative.
///
/// # Formula
///
/// ```text
/// roughness = delta_h * (1 - 0.8 * exp(-d / 50000))
/// ```
///
/// where:
/// - `delta_h` = interdecile range (meters)
/// - `d` = distance (meters)
/// - 50000 = smoothing distance constant (50 km)
/// - 0.8 = empirical scaling factor
///
/// # Behavior
///
/// - For short distances (d << 50 km): Roughness approaches 0, as local variations dominate
/// - For long distances (d >> 50 km): Roughness approaches delta_h, as the full statistical
///   variation becomes relevant
/// - At d = 50 km: Roughness is approximately 0.43 * delta_h
///
/// # References
///
/// - Longley, A.G. & Rice, P.L. (1968). "Prediction of Tropospheric Radio Transmission Loss
///   Over Irregular Terrain: A Computer Method-1968". ESSA Technical Report ERL 79-ITS 67, Equation 3.
/// - ITM Technical Note 101 (TN101): Terrain irregularity parameter calculation
///
/// # Examples
///
/// ```
/// use itm_rs::math::terrain_roughness;
///
/// // Calculate roughness for a 10 km path with 100m elevation variation
/// let roughness = terrain_roughness(10_000.0, 100.0);
/// assert!(roughness > 0.0 && roughness < 100.0);
///
/// // Longer distances result in higher roughness for same delta_h
/// let roughness_short = terrain_roughness(10_000.0, 100.0);
/// let roughness_long = terrain_roughness(100_000.0, 100.0);
/// assert!(roughness_long > roughness_short);
/// ```
pub fn terrain_roughness(d: f64, delta_h: f64) -> f64 {
    // [ERL 79-ITS 67, Eqn 3], with distance in meters instead of kilometers
    delta_h * (1.0 - 0.8 * (-d / 50e3).exp())
}

/// Performs a linear least-squares fit on terrain profile data.
///
/// This function computes a best-fit line for terrain elevation data using the
/// least-squares method. It selects a subset of the profile data within the specified
/// distance range and calculates the fitted line's parameters (intercept and slope).
/// The fit minimizes the sum of squared vertical deviations from the line.
///
/// # Arguments
///
/// * `pfl` - Profile array with metadata:
///   - `pfl[0]` = Number of points in the profile
///   - `pfl[1]` = Distance increment between profile points (meters)
/// * `d_start` - Starting distance for the fit region in meters
/// * `d_end` - Ending distance for the fit region in meters
///
/// # Returns
///
/// A tuple `(intercept, slope)` representing the fitted line:
/// - `intercept` = Y-intercept of the fitted line
/// - `slope` = Slope of the fitted line (change in elevation per unit distance)
///
/// # Algorithm Notes
///
/// - Uses centered coordinate system to improve numerical stability
/// - Adjusts fit region if the resulting range is invalid (empty or inverted)
/// - Applies a scaling factor to convert the raw slope to the final slope value
/// - Extrapolates the line to both endpoints of the original profile
///
/// # References
///
/// - Longley, A.G. & Rice, P.L. (1968). "Prediction of Tropospheric Radio Transmission Loss Over Irregular Terrain: A Computer Method-1968". ESSA Technical Report ERL 79-ITS 67.
/// - Standard least-squares regression method with centered coordinates for numerical stability
/// - ITM Technical Note 101 (TN101): Terrain profile analysis methods
pub(crate) fn linear_least_square_fit(pfl: &[f64], d_start: f64, d_end: f64) -> (f64, f64) {
    let num_points = pfl[0] as usize;
    let distance_step = pfl[1];

    // Calculate starting and ending indices within the profile data
    let mut index_start = ((d_start / distance_step).max(0.0)) as usize;
    let mut index_end = num_points - (pfl[0].max(d_end / distance_step)) as usize;

    // Ensure valid range; if empty or invalid, adjust boundaries
    if index_end <= index_start {
        index_start = index_start.max(1);
        index_end = num_points - num_points.max(index_end + 1);
    }

    let sample_count = (index_end - index_start) as f64;

    // Centered coordinate system: shift the midpoint to origin for numerical stability
    let center_offset_i32 = (-0.5 * sample_count) as i32;
    let center_position = (index_end as f64) + (center_offset_i32 as f64);

    // Initialize sums for least-squares calculation
    // Start with average of endpoint values
    let mut sum_elevations = 0.5 * (pfl[index_start + 2] + pfl[index_end + 2]);
    let mut weighted_sum =
        0.5 * (pfl[index_start + 2] - pfl[index_end + 2]) * (center_offset_i32 as f64);

    // Accumulate elevation values and weighted values across the sample range
    let mut current_index = index_start;
    for iteration in 2..=(sample_count as usize) {
        current_index += 1;
        let current_offset = (center_offset_i32 + (iteration as i32)) as f64;

        sum_elevations += pfl[current_index + 2];
        weighted_sum += pfl[current_index + 2] * current_offset;
    }

    // Normalize the sums
    let mean_elevation = sum_elevations / sample_count;

    // Scale the weighted sum to get the slope
    // Factor: 12.0 / ((sample_count^2 + 2) * sample_count)
    let scaled_slope = weighted_sum * 12.0 / ((sample_count * sample_count + 2.0) * sample_count);

    // Calculate intercept and slope for the fitted line
    let intercept = mean_elevation - scaled_slope * center_position;
    let slope = mean_elevation + scaled_slope * ((num_points as f64) - center_position);

    (intercept, slope)
}

/// Siting criteria for terminal placement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SitingCriteria {
    /// Random siting.
    Random,
    /// Careful siting.
    Careful,
    /// Very careful siting.
    VeryCareful,
}

/// Initializes area mode calculations.
///
/// This function computes effective terminal heights, horizon distances, and
/// horizon angles for area mode propagation based on siting criteria, effective
/// earth curvature, terrain irregularity, and structural heights.
///
/// # Arguments
///
/// * `site_criteria` - Siting criteria for each terminal.
/// * `gamma_e` - Curvature of the effective earth (1/m).
/// * `delta_h_m` - Terrain irregularity parameter in meters.
/// * `h_m` - Terminal structural heights in meters.
///
/// # Returns
///
/// A tuple `(h_e, d_hzn, theta_hzn)` where:
/// - `h_e` is the effective terminal heights in meters.
/// - `d_hzn` is the terminal horizon distances in meters.
/// - `theta_hzn` is the terminal horizon angles in radians.
///
/// # References
///
/// - Longley-Rice ITM Algorithm, Eqn 3.2: Effective height adjustment.
/// - Longley-Rice ITM Algorithm, Eqn 3.3: Horizon distance calculation.
/// - Longley-Rice ITM Algorithm, Eqn 3.4: Horizon angle calculation.
///
/// # Examples
///
/// ```
/// use itm_rs::math::{initialize_area, SitingCriteria};
///
/// let site_criteria = [SitingCriteria::Random, SitingCriteria::Careful];
/// let h_m = [30.0, 20.0];
/// let (h_e, d_hzn, theta_hzn) = initialize_area(site_criteria, 1.57e-7, 90.0, h_m);
///
/// assert!(h_e[0].is_finite());
/// assert!(d_hzn[0].is_finite());
/// assert!(theta_hzn[0].is_finite());
/// ```
pub fn initialize_area(
    site_criteria: [SitingCriteria; 2],
    gamma_e: f64,
    delta_h_m: f64,
    h_m: [f64; 2],
) -> ([f64; 2], [f64; 2], [f64; 2]) {
    let mut h_e = [0.0; 2];
    let mut d_hzn = [0.0; 2];
    let mut theta_hzn = [0.0; 2];

    for i in 0..2 {
        match site_criteria[i] {
            SitingCriteria::Random => {
                h_e[i] = h_m[i];
            }
            SitingCriteria::Careful | SitingCriteria::VeryCareful => {
                let mut b = match site_criteria[i] {
                    SitingCriteria::Careful => 4.0,
                    SitingCriteria::VeryCareful => 9.0,
                    SitingCriteria::Random => unreachable!(),
                };

                if h_m[i] < 5.0 {
                    b *= (0.1 * PI * h_m[i]).sin();
                }

                // [Algorithm, Eqn 3.2]
                let denom = (2.0 * h_m[i] / delta_h_m.max(1e-3)).min(20.0);
                h_e[i] = h_m[i] + (1.0 + b) * (-denom).exp();
            }
        }

        let d_ls_m = (2.0 * h_e[i] / gamma_e).sqrt();

        // [Algorithm, Eqn 3.3]
        const H_3_M: f64 = 5.0;
        d_hzn[i] = d_ls_m * (-0.07 * (delta_h_m / h_e[i].max(H_3_M)).sqrt()).exp();

        // [Algorithm, Eqn 3.4]
        theta_hzn[i] = (0.65 * delta_h_m * (d_ls_m / d_hzn[i] - 1.0) - 2.0 * h_e[i]) / d_ls_m;
    }

    (h_e, d_hzn, theta_hzn)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_initialize_area_outputs_finite() {
        let site_criteria = [SitingCriteria::Random, SitingCriteria::VeryCareful];
        let h_m = [30.0, 20.0];

        let (h_e, d_hzn, theta_hzn) = initialize_area(site_criteria, 1.57e-7, 90.0, h_m);

        for i in 0..2 {
            assert!(h_e[i].is_finite());
            assert!(d_hzn[i].is_finite());
            assert!(theta_hzn[i].is_finite());
        }
    }
}

/// Computes the radio horizons for the terminals.
///
/// This function scans the terrain profile to determine the radio horizon angles
/// and distances for both terminals. It starts with the line-of-sight assumption
/// and updates the horizon angles when intervening terrain points obstruct the
/// line of sight.
///
/// # Arguments
///
/// * `pfl` - Terrain profile array with metadata:
///   - `pfl[0]` = Number of profile points
///   - `pfl[1]` = Distance increment between points (meters)
///   - `pfl[2..]` = Elevation samples along the path (meters)
/// * `a_e` - Effective earth radius in meters.
/// * `h_m` - Terminal structural heights in meters `[tx, rx]`.
///
/// # Returns
///
/// A tuple `(theta_hzn, d_hzn)` where:
/// - `theta_hzn` = Terminal radio horizon angles in radians `[tx, rx]`
/// - `d_hzn` = Terminal radio horizon distances in meters `[tx, rx]`
///
/// # References
///
/// - ITM Technical Note 101 (TN101), Eq 6.15: Initial horizon angles.
/// - Longley-Rice ITM Algorithm: Horizon search procedure.
///
/// # Examples
///
/// ```
/// use itm_rs::math::find_horizons;
///
/// // Minimal terrain profile: 3 points, 1 km spacing, flat terrain at 0 m.
/// let pfl = [3.0, 1000.0, 0.0, 0.0, 0.0];
/// let h_m = [10.0, 10.0];
/// let (theta_hzn, d_hzn) = find_horizons(&pfl, 8.5e6, h_m);
///
/// assert!(theta_hzn[0].is_finite());
/// assert!(d_hzn[0] > 0.0);
/// ```
///
pub fn find_horizons(pfl: &[f64], a_e: f64, h_m: [f64; 2]) -> ([f64; 2], [f64; 2]) {
    let np = pfl[0] as usize;
    let xi = pfl[1];

    let d_m = pfl[0] * pfl[1];

    // Compute terminal elevations including structural heights.
    let z_tx_m = pfl[2] + h_m[0];
    let z_rx_m = pfl[np + 2] + h_m[1];

    // Initial horizon angles assuming line-of-sight. [TN101, Eq 6.15]
    let mut theta_hzn = [
        (z_rx_m - z_tx_m) / d_m - d_m / (2.0 * a_e),
        -(z_rx_m - z_tx_m) / d_m - d_m / (2.0 * a_e),
    ];

    let mut d_hzn = [d_m, d_m];

    let mut d_tx_m = 0.0;
    let mut d_rx_m = d_m;

    for i in 1..np {
        d_tx_m += xi;
        d_rx_m -= xi;

        let theta_tx = (pfl[i + 2] - z_tx_m) / d_tx_m - d_tx_m / (2.0 * a_e);
        let theta_rx = -(z_rx_m - pfl[i + 2]) / d_rx_m - d_rx_m / (2.0 * a_e);

        if theta_tx > theta_hzn[0] {
            theta_hzn[0] = theta_tx;
            d_hzn[0] = d_tx_m;
        }

        if theta_rx > theta_hzn[1] {
            theta_hzn[1] = theta_rx;
            d_hzn[1] = d_rx_m;
        }
    }

    (theta_hzn, d_hzn)
}

#[cfg(test)]
mod find_horizons_tests {
    use super::find_horizons;

    #[test]
    fn test_find_horizons_flat_profile() {
        let pfl = [3.0, 1000.0, 0.0, 0.0, 0.0];
        let h_m = [10.0, 10.0];
        let (theta_hzn, d_hzn) = find_horizons(&pfl, 8.5e6, h_m);

        assert!(theta_hzn[0].is_finite());
        assert!(theta_hzn[1].is_finite());
        assert!(d_hzn[0] > 0.0);
        assert!(d_hzn[1] > 0.0);
    }
}

/// Computes the terrain irregularity parameter (delta_h).
///
/// This function calculates delta_h by sampling the terrain profile over a range,
/// removing the best-fit linear trend, and computing the interdecile range between
/// the 10th and 90th percentiles of the residuals. It then inverts the ITM smoothing
/// factor to recover the terrain irregularity parameter for the specified range.
///
/// # Arguments
///
/// * `pfl` - Terrain profile array with metadata:
///   - `pfl[0]` = Number of profile points
///   - `pfl[1]` = Distance increment between points (meters)
///   - `pfl[2..]` = Elevation samples along the path (meters)
/// * `d_start_m` - Distance into the profile to start considering data (meters).
/// * `d_end_m` - Distance into the profile to stop considering data (meters).
///
/// # Returns
///
/// Terrain irregularity parameter delta_h in meters.
///
/// # References
///
/// - Longley-Rice ITM Algorithm: Interdecile-based irregularity estimation.
/// - ERL 79-ITS 67, Eqn 3 (inverted): Delta_h smoothing correction.
///
/// # Examples
///
/// ```
/// use itm_rs::math::compute_delta_h;
///
/// let pfl = [5.0, 1000.0, 0.0, 10.0, 5.0, 12.0, 8.0];
/// let delta_h = compute_delta_h(&pfl, 0.0, 4000.0);
/// assert!(delta_h >= 0.0);
/// ```
pub fn compute_delta_h(pfl: &[f64], d_start_m: f64, d_end_m: f64) -> f64 {
    let np = pfl[0] as usize;
    let mut x_start = d_start_m / pfl[1];
    let mut x_end = d_end_m / pfl[1];

    // If there are fewer than 2 terrain points, delta_h = 0.
    if x_end - x_start < 2.0 {
        return 0.0;
    }

    let mut p10 = (0.1 * (x_end - x_start + 8.0)) as usize;
    p10 = p10.clamp(4, 25);

    let n = 10 * p10 - 5;
    let p90 = n - p10;

    let np_s = (n - 1) as f64;
    let mut s = vec![0.0; n + 2];
    s[0] = np_s;
    s[1] = 1.0;

    x_end = (x_end - x_start) / np_s;
    let mut i = x_start as usize;
    x_start -= i as f64 + 1.0;

    for j in 0..n {
        while x_start > 0.0 && (i + 1) < np {
            x_start -= 1.0;
            i += 1;
        }

        s[j + 2] = pfl[i + 3] + (pfl[i + 3] - pfl[i + 2]) * x_start;
        x_start += x_end;
    }

    let (mut fit_y1, mut fit_y2) = linear_least_square_fit(&s, 0.0, np_s);
    fit_y2 = (fit_y2 - fit_y1) / np_s;

    let mut diffs = Vec::with_capacity(n);
    for _ in 0..n {
        diffs.push(s[diffs.len() + 2] - fit_y1);
        fit_y1 += fit_y2;
    }

    diffs.select_nth_unstable_by(p10 - 1, |a, b| {
        b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal)
    });
    let q10 = diffs[p10 - 1];

    diffs.select_nth_unstable_by(p90, |a, b| {
        b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal)
    });
    let q90 = diffs[p90];

    let delta_h_d_m = q10 - q90;

    // [ERL 79-ITS 67, Eqn 3], inverted.
    let smoothing = 1.0 - 0.8 * (-(d_end_m - d_start_m) / 50e3).exp();
    delta_h_d_m / smoothing
}

#[cfg(test)]
mod compute_delta_h_tests {
    use super::compute_delta_h;

    #[test]
    fn test_compute_delta_h_non_negative() {
        let pfl = [5.0, 1000.0, 0.0, 10.0, 5.0, 12.0, 8.0];
        let delta_h = compute_delta_h(&pfl, 0.0, 4000.0);
        assert!(delta_h >= 0.0);
    }
}
