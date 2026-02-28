/*!
Centralized mathematical and physical constants for the `math` submodules.

Place commonly reused numerical constants here to avoid duplication, improve
readability and make future tuning easier. Prefer `f64` typed constants for
numeric stability across computation-heavy modules.

This file intentionally contains only pure constants and very small helper
`const` expressions. No functions or side effects: including this module is
cheap and safe in all math submodules.
*/

/// Mean earth radius in meters (approx). Used throughout the model.
///
/// Typical value: ~6.37e6 m (6370 km).
pub const EARTH_RADIUS_M: f64 = 6_370_000.0;

/// One-third constant used in Vogler formulas and other 1/3 exponents.
pub const THIRD: f64 = 1.0 / 3.0;

/// Divisor used to convert MHz frequency to an approximate wavenumber in
/// several legacy formulas (many places historically use `f / 47.7`).
pub const WAVENUMBER_DIVISOR: f64 = 47.7;

/// Curvature of the actual (geometric) earth: ~1 / (6370 km).
pub const GAMMA_A: f64 = 157e-9;

/// Reciprocal of 4π, used in Fresnel parameter calculations and related formulas.
pub const RECIP_FOUR_PI: f64 = 1.0 / (4.0 * std::f64::consts::PI);

/// Smoothing distance used in troposcatter smoothing terms (meters).
pub const SMOOTHING_DISTANCE_D0: f64 = 40_000.0; // 40 km

/// Height constant used in troposcatter formula (meters) — historical value used
/// in ITM implementations.
pub const HEIGHT_CONSTANT_H: f64 = 47.7;

/// Typical limits for refractivity `n0` used in input validation (N-units).
pub const REFRACTIVITY_MIN: f64 = 250.0;
pub const REFRACTIVITY_MAX: f64 = 400.0;

/// Typical acceptable frequency ranges (MHz) used for generating warnings/errors.
pub const FREQ_WARN_MIN_MHZ: f64 = 40.0;
pub const FREQ_WARN_MAX_MHZ: f64 = 10_000.0;
pub const FREQ_ERR_MIN_MHZ: f64 = 20.0;
pub const FREQ_ERR_MAX_MHZ: f64 = 20_000.0;

/// Default sentinel value returned by troposcatter when the scatter function is
/// not defined (kept consistent with current implementation).
pub const TROPOSCATTER_UNDEFINED_SENTINEL: f64 = 1001.0;

/// Commonly-used numeric tolerances.
pub const EPS_F64: f64 = 1e-12;
pub const LARGE_F64: f64 = 1.0e300;

/// A small helper const for converting kilometers to meters (1 km = 1000 m).
pub const KM_TO_M: f64 = 1000.0;

/// A small helper const for converting meters to kilometers.
pub const M_TO_KM: f64 = 0.001;
