//! Mathematical functions for the ITM (Irregular Terrain Model).
//!
//! This module provides mathematical utilities for radio propagation modeling,
//! organized into the following submodules:
//!
//! - [`propagation`] - Free space loss and basic propagation calculations
//! - [`diffraction`] - Diffraction modeling (Fresnel integral, knife-edge)
//! - [`scatter`] - Tropospheric scatter propagation
//! - [`height_gain`] - Height gain factor calculations
//! - [`terrain`] - Terrain profile analysis and fitting
//! - [`statistics`] - Statistical functions and distributions
//! - [`variability`] - Variability loss calculations
//!
//! # Overview
//!
//! These functions implement the core mathematical algorithms from the Longley-Rice
//! Irregular Terrain Model (ITM), a comprehensive radio propagation prediction method
//! developed by A.G. Longley and P.L. Rice in 1968.
//!
//! # References
//!
//! - Longley, A.G. & Rice, P.L. (1968). "Prediction of Tropospheric Radio Transmission
//!   Loss Over Irregular Terrain: A Computer Method-1968". ESSA Technical Report ERL 79-ITS 67.
//! - ITM Technical Note 101 (TN101): Detailed implementation notes and equations

// Submodules
pub mod constants;
pub mod diffraction;
pub mod height_gain;
pub mod itm;
pub mod longley_rice;
pub mod propagation;
pub mod scatter;
pub mod statistics;
pub mod terrain;
pub mod variability;

// Re-export public API for backward compatibility and convenience
pub use constants::{
    EARTH_RADIUS_M, EPS_F64, FREQ_ERR_MAX_MHZ, FREQ_ERR_MIN_MHZ, FREQ_WARN_MAX_MHZ,
    FREQ_WARN_MIN_MHZ, GAMMA_A, HEIGHT_CONSTANT_H, KM_TO_M, LARGE_F64, M_TO_KM, RECIP_FOUR_PI,
    REFRACTIVITY_MAX, REFRACTIVITY_MIN, SMOOTHING_DISTANCE_D0, THIRD,
    TROPOSCATTER_UNDEFINED_SENTINEL, WAVENUMBER_DIVISOR,
};
pub use diffraction::{fresnel_integral, smooth_earth_diffraction};
pub use height_gain::{h0_curve, h0_function};
pub use itm::{IntermediateValues, ItmError, itm_area_tls, validate_inputs};
pub use propagation::{Polarization, free_space_loss, initialize_point_to_point};
pub use scatter::{f_function, troposcatter_loss};
pub use statistics::{inverse_ccdf, inverse_complementary_cumulative_distribution_function};
pub use terrain::{
    SitingCriteria, compute_delta_h, find_horizons, initialize_area, terrain_roughness,
};
pub use variability::{Climate, VariabilityMode, curve, variability_loss};

// Re-export crate-internal functions
