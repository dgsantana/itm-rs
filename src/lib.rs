// Copyright (c) 2026 Daniel Santana
// SPDX-License-Identifier: MIT

//! # ITM-RS: Irregular Terrain Model for Rust
//!
//! A high-performance, scientifically-validated Rust implementation of the ITS Irregular
//! Terrain Model (ITM), also known as the Longley-Rice model.
//!
//! ## Overview
//!
//! This library provides radio propagation prediction for frequencies between 20 MHz and
//! 20 GHz, supporting both point-to-point (P2P) and area prediction modes.
//!
//! ## Original Work Attribution
//!
//! This is a derivative work of the ITS Irregular Terrain Model originally developed by
//! the National Telecommunications and Information Administration (NTIA). The original
//! ITM software is in the public domain as a work of the U.S. Federal Government.
//!
//! Original authors: Anita Longley and Phil Rice
//!
//! ## License
//!
//! This Rust implementation is licensed under the MIT License.
//! See LICENSE file for details.
//!
//! ## References
//!
//! - Original NTIA repository: <https://github.com/NTIA/itm>
//! - Hufford, G.A., Longley, A.G., Kissick, W.A. (1982). "A Guide to the Use of the
//!   ITS Irregular Terrain Model in the Area Prediction Mode." NTIA Technical Report
//!   TR-82-100.
//! - Longley, A.G., Rice, P.L. (1968). "Prediction of Tropospheric Radio Transmission
//!   Loss Over Irregular Terrain: A Computer Method - 1968." NTIA Technical Report
//!   ERL 79-ITS 67.

pub mod ffi;
pub mod math;
