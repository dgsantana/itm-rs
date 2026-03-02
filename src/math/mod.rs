//! Mathematical functions for the ITM (Irregular Terrain Model).
//!
//! This module provides mathematical utilities for radio propagation modeling,
//! organized into submodules by function category.
//!
//! # Overview
//!
//! These functions implement the core mathematical algorithms from the Longley-Rice
//! Irregular Terrain Model (ITM), a comprehensive radio propagation prediction method
//! developed by A.G. Longley and P.L. Rice in 1968.
//!
//! # Submodules
//!
//! - [`constants`] - Physical and mathematical constants
//! - [`diffraction`] - Diffraction loss calculations
//! - [`height_gain`] - Height gain (terrain) functions
//! - [`itm`] - Main ITM prediction functions (P2P and area modes)
//! - [`propagation`] - Basic propagation models (free space loss, etc.)
//! - [`scatter`] - Troposcatter calculations
//! - [`statistics`] - Statistical functions
//! - [`terrain`] - Terrain analysis and initialization
//! - [`variability`] - Variability and reliability calculations
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
pub(crate) mod longley_rice;
pub mod propagation;
pub mod scatter;
pub mod statistics;
pub mod terrain;
pub mod variability;

/// Prelude module for convenient imports of math module submodules.
pub mod prelude {
    pub use super::constants::*;
    pub use super::diffraction::*;
    pub use super::height_gain::*;
    pub use super::itm::*;
    pub use super::propagation::*;
    pub use super::scatter::*;
    pub use super::statistics::*;
    pub use super::terrain::*;
    pub use super::variability::*;
}
