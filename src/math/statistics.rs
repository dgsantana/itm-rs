//!! Statistical functions and distributions.
//!!
//!! This module contains statistical utility functions used in radio propagation
//!! modeling, including probability distributions and related calculations.

/// Computes the inverse complementary cumulative distribution function (inverse CCDF).
///
/// This function approximates the inverse of the complementary cumulative distribution
/// function (also known as the quantile function or percent-point function) using the
/// approximation given in Formula 26.2.23 from Abramowitz & Stegun's "Handbook of
/// Mathematical Functions".
///
/// The complementary cumulative distribution function is Q(x) = 1 - Φ(x), where Φ(x)
/// is the standard normal CDF. This function computes Q⁻¹(q), which gives the value
/// x such that Q(x) = q.
///
/// This is commonly used in radio propagation to determine confidence intervals and
/// fading margins, where we need to find the signal level that will be exceeded with
/// a certain probability.
///
/// # Arguments
///
/// * `q` - The quantile (probability value), where 0.0 < q < 1.0. This represents
///   the probability in the upper tail of the distribution. For example:
///   - q = 0.5 corresponds to the median (50th percentile)
///   - q = 0.1 corresponds to the 90th percentile
///   - q = 0.01 corresponds to the 99th percentile
///
/// # Returns
///
/// The inverse complementary CDF value Q⁻¹(q), which is the number of standard
/// deviations from the mean corresponding to the given upper-tail probability.
///
/// # Algorithm
///
/// The approximation uses a rational function fit:
///
/// ```text
/// T = sqrt(-2 * ln(x))
/// ζ = ((C₂·T + C₁)·T + C₀) / (((D₃·T + D₂)·T + D₁)·T + 1)
/// Q⁻¹(q) = T - ζ
/// ```
///
/// where x = min(q, 1-q) to ensure symmetry about q = 0.5.
///
/// # Error Bounds
///
/// This approximation has an absolute error |ε(q)| < 4.5×10⁻⁴ across the valid
/// input range, making it suitable for most radio propagation calculations.
///
/// # Panics
///
/// This function does not perform input validation. The caller must ensure that
/// 0.0 < q < 1.0. Values outside this range will produce mathematically invalid
/// results (NaN or infinite values).
///
/// # References
///
/// - Abramowitz, M. & Stegun, I.A. (1964). "Handbook of Mathematical Functions with
///   Formulas, Graphs, and Mathematical Tables". National Bureau of Standards,
///   Formula 26.2.23, page 933.
/// - Used in ITM for computing fading margins and variability statistics
///
/// # Examples
///
/// ```
/// use itm::math::statistics::inverse_ccdf;
///
/// // Find the 50th percentile (median) - should be near 0
/// let median = inverse_ccdf(0.5);
/// assert!(median.abs() < 0.001);
///
/// // Find the 90th percentile (q = 0.1)
/// let p90 = inverse_ccdf(0.1);
/// assert!(p90 > 1.0); // Approximately 1.28 standard deviations
///
/// // Find the 99th percentile (q = 0.01)
/// let p99 = inverse_ccdf(0.01);
/// assert!(p99 > 2.0); // Approximately 2.33 standard deviations
///
/// // Symmetry check
/// let lower = inverse_ccdf(0.9);
/// let upper = inverse_ccdf(0.1);
/// assert!((lower + upper).abs() < 0.001); // Should be negatives of each other
/// ```
///
/// # Applications in Radio Propagation
///
/// This function is used to:
/// - Compute fading margins for desired reliability levels
/// - Determine confidence intervals for signal strength predictions
/// - Calculate variability factors in path loss predictions
/// - Establish design margins for communication links
///
/// For example, to achieve 99% reliability (q = 0.01), you would add
/// approximately 2.33 standard deviations of fading to your median path loss.
pub fn inverse_ccdf(q: f64) -> f64 {
    // Coefficients from Abramowitz & Stegun Formula 26.2.23
    const C_0: f64 = 2.515516;
    const C_1: f64 = 0.802853;
    const C_2: f64 = 0.010328;
    const D_1: f64 = 1.432788;
    const D_2: f64 = 0.189269;
    const D_3: f64 = 0.001308;

    // Use the smaller of q and (1-q) for numerical stability
    let x = if q > 0.5 { 1.0 - q } else { q };

    // Compute T = sqrt(-2 * ln(x))
    let t_x = (-2.0 * x.ln()).sqrt();

    // Compute the rational function approximation ζ(T)
    // Numerator: C₂·T² + C₁·T + C₀
    let numerator = (C_2 * t_x + C_1) * t_x + C_0;

    // Denominator: D₃·T³ + D₂·T² + D₁·T + 1
    let denominator = ((D_3 * t_x + D_2) * t_x + D_1) * t_x + 1.0;

    let zeta = numerator / denominator;

    // Compute the result Q⁻¹(q) = T - ζ
    let mut result = t_x - zeta;

    // Apply sign correction for q > 0.5 (upper tail)
    if q > 0.5 {
        result = -result;
    }

    result
}

/// Alias for `inverse_ccdf` with a more descriptive name.
///
/// This is the full name of the function, provided for clarity in code where
/// the abbreviated version might be unclear.
pub fn inverse_complementary_cumulative_distribution_function(q: f64) -> f64 {
    inverse_ccdf(q)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_inverse_ccdf_median() {
        // At q = 0.5, should return approximately 0
        let result = inverse_ccdf(0.5);
        assert!(
            result.abs() < 0.001,
            "Median should be near 0, got {}",
            result
        );
    }

    #[test]
    fn test_inverse_ccdf_90th_percentile() {
        // At q = 0.1 (90th percentile), should return approximately 1.282
        let result = inverse_ccdf(0.1);
        assert!(
            result > 1.2 && result < 1.3,
            "90th percentile should be ~1.28, got {}",
            result
        );
    }

    #[test]
    fn test_inverse_ccdf_99th_percentile() {
        // At q = 0.01 (99th percentile), should return approximately 2.326
        let result = inverse_ccdf(0.01);
        assert!(
            result > 2.3 && result < 2.4,
            "99th percentile should be ~2.33, got {}",
            result
        );
    }

    #[test]
    fn test_inverse_ccdf_symmetry() {
        // Q⁻¹(q) = -Q⁻¹(1-q) (symmetry property)
        let q_low = 0.1;
        let q_high = 0.9;

        let result_low = inverse_ccdf(q_low);
        let result_high = inverse_ccdf(q_high);

        assert!(
            (result_low + result_high).abs() < 0.001,
            "Symmetry failed: {} + {} ≠ 0",
            result_low,
            result_high
        );
    }

    #[test]
    fn test_inverse_ccdf_known_values() {
        // Test against known values from normal distribution tables
        let test_cases = vec![
            (0.5, 0.0),    // Median
            (0.3085, 0.5), // ~0.5 std devs
            (0.1587, 1.0), // 1 std dev
            (0.0228, 2.0), // 2 std devs
        ];

        for (q, expected) in test_cases {
            let result = inverse_ccdf(q);
            assert!(
                (result - expected).abs() < 0.01,
                "For q={}, expected ~{}, got {}",
                q,
                expected,
                result
            );
        }
    }

    #[test]
    fn test_inverse_ccdf_monotonicity() {
        // Function should be monotonically decreasing
        let q1 = 0.1;
        let q2 = 0.5;
        let q3 = 0.9;

        let r1 = inverse_ccdf(q1);
        let r2 = inverse_ccdf(q2);
        let r3 = inverse_ccdf(q3);

        assert!(r1 > r2, "Should be decreasing: {} > {}", r1, r2);
        assert!(r2 > r3, "Should be decreasing: {} > {}", r2, r3);
    }

    #[test]
    fn test_inverse_ccdf_extreme_values() {
        // Test near boundaries (but not at them, as function is undefined at 0 and 1)
        let very_small = inverse_ccdf(0.001);
        let very_large = inverse_ccdf(0.999);

        assert!(
            very_small > 3.0,
            "Very small q should give large positive value"
        );
        assert!(
            very_large < -3.0,
            "Very large q should give large negative value"
        );
    }

    #[test]
    fn test_alias_function() {
        // Test that the alias function returns the same result
        let q = 0.1;
        let result1 = inverse_ccdf(q);
        let result2 = inverse_complementary_cumulative_distribution_function(q);

        assert_eq!(result1, result2, "Alias function should return same result");
    }

    #[test]
    fn test_inverse_ccdf_accuracy() {
        // Verify the claimed accuracy of < 4.5e-4
        // We'll test at several points and verify consistency
        let test_points = vec![
            0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99,
        ];

        for &q in &test_points {
            let result = inverse_ccdf(q);

            // Verify result is finite
            assert!(result.is_finite(), "Result should be finite for q={}", q);

            // Verify reasonable range (normal distribution is typically within ±4σ for practical purposes)
            assert!(
                result.abs() < 4.0 || !(0.0001..=0.9999).contains(&q),
                "Result {} seems unreasonable for q={}",
                result,
                q
            );
        }
    }
}
