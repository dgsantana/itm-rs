/*!
Core Longley-Rice helpers extracted from the ITM orchestration.

This module contains the lower-level Longley-Rice pieces used by the public
wrappers in `itm.rs`. It is intended to be an internal implementation module
that provides:

- `quick_pfl` — a lightweight profile preprocessor used by the P2P path
  interfaces (temporary helper until a full profile port exists)
- `diffraction_loss` — a wrapper over the diffraction implementation (currently
  using Vogler/smooth-earth approximation)
- `line_of_sight_loss` — a small wrapper that proxies to diffraction logic for
  short distances
- `longley_rice` — the (partial) Longley-Rice core algorithm used to compute
  the reference attenuation `A_ref` and select propagation mode.

Note: this file intentionally focuses on the algorithmic core and keeps public
exposure crate-local (`pub(crate)`) so the high-level API in `itm.rs` can stay
stable and small.
*/

use crate::math::diffraction::diffraction_loss;
use crate::math::itm::{ItmError, warnings as itm_warnings};
use crate::math::propagation::line_of_sight_loss;
use crate::math::scatter::troposcatter_loss;
use num_complex::Complex;

/// Result of the Longley-Rice core algorithm.
pub(crate) struct LongleyRiceResult {
    /// The reference attenuation `A_ref` in dB.
    pub a_ref: f64,
    /// Bitmask of warnings accumulated during computation.
    pub warnings: u32,
    /// Propagation mode selected (1 for LOS, 2 for Diffraction, 3 for Troposcatter).
    pub propmode: i32,
}

/// Partially-ported Longley-Rice core.
///
/// Computes the reference attenuation `A_ref` for a path and selects a propagation
/// mode. This is a pragmatic port intended to match the original algorithm structure:
/// diffraction (including smooth-earth), transition handling, and troposcatter selection.
///
/// Many internal numeric heuristics are directly transcribed from canonical ITM
/// references and from the initial project port.
///
/// # Arguments
///
/// * `theta_hzn` - Horizon angles [tx, rx].
/// * `f` - Frequency in MHz.
/// * `z_g` - Complex ground impedance.
/// * `d_hzn` - Horizon distances [tx, rx].
/// * `h_e` - Effective terminal heights [tx, rx].
/// * `gamma_e` - Curvature of the effective earth (1/a_e).
/// * `n_s` - Surface refractivity in N-Units.
/// * `delta_h` - Terrain roughness parameter.
/// * `h` - Terminal heights [tx, rx].
/// * `d` - Total path distance in meters.
/// * `mode` - Propagation mode flag (0 for P2P).
///
/// # Returns
///
/// A `LongleyRiceResult` struct indicating the loss, warnings, and mode, or an `ItmError`
/// representing fatal input/algorithm conditions.
pub(crate) fn longley_rice(
    theta_hzn: [f64; 2],
    f: f64,
    z_g: Complex<f64>,
    d_hzn: [f64; 2],
    h_e: [f64; 2],
    gamma_e: f64,
    n_s: f64,
    delta_h: f64,
    h: [f64; 2],
    d: f64,
    _mode: i32,
) -> Result<LongleyRiceResult, ItmError> {
    let mut warnings = 0u32;
    // effective earth radius
    let a_e = 1.0 / gamma_e;

    // Terrestrial smooth earth horizon distances
    let mut d_hzn_s = [0.0f64; 2];
    for i in 0..2 {
        d_hzn_s[i] = (2.0 * h_e[i] * a_e).sqrt();
    }

    let d_sml = d_hzn_s[0] + d_hzn_s[1];
    let d_ml = d_hzn[0] + d_hzn[1];

    let theta_los = -((theta_hzn[0] + theta_hzn[1]).max(-d_ml / a_e));

    if theta_hzn[0].abs() > 200e-3 {
        warnings |= itm_warnings::WARN_TX_HORIZON_ANGLE;
    }
    if theta_hzn[1].abs() > 200e-3 {
        warnings |= itm_warnings::WARN_RX_HORIZON_ANGLE;
    }

    if d_hzn[0] < 0.1 * d_hzn_s[0] {
        warnings |= itm_warnings::WARN_TX_HORIZON_DISTANCE_1;
    }
    if d_hzn[1] < 0.1 * d_hzn_s[1] {
        warnings |= itm_warnings::WARN_RX_HORIZON_DISTANCE_1;
    }

    if d_hzn[0] > 3.0 * d_hzn_s[0] {
        warnings |= itm_warnings::WARN_TX_HORIZON_DISTANCE_2;
    }
    if d_hzn[1] > 3.0 * d_hzn_s[1] {
        warnings |= itm_warnings::WARN_RX_HORIZON_DISTANCE_2;
    }

    if n_s < 150.0 {
        return Err(ItmError::ErrorSurfaceRefractivitySmall);
    }
    if n_s > 400.0 {
        return Err(ItmError::ErrorSurfaceRefractivityLarge);
    }
    if n_s < 250.0 {
        warnings |= itm_warnings::WARN_SURFACE_REFRACTIVITY;
    }

    if a_e < 4_000_000.0 || a_e > 13_333_333.0 {
        return Err(ItmError::ErrorEffectiveEarth);
    }

    if z_g.re <= z_g.im.abs() {
        return Err(ItmError::ErrorGroundImpedance);
    }

    let d_3 = d_sml.max(d_ml + 5.0 * ((a_e.powi(2) / f).powf(1.0 / 3.0)));
    let d_4 = d_3 + 10.0 * ((a_e.powi(2) / f).powf(1.0 / 3.0));

    let a_3_db = diffraction_loss(
        d_3, d_hzn, h_e, z_g, a_e, delta_h, h, 0, theta_los, d_sml, f,
    );
    let a_4_db = diffraction_loss(
        d_4, d_hzn, h_e, z_g, a_e, delta_h, h, 0, theta_los, d_sml, f,
    );

    let m_d = (a_4_db - a_3_db) / (d_4 - d_3);
    let a_d0 = a_3_db - m_d * d_3;

    let d_min = (h_e[0] - h_e[1]).abs() / 200e-3;
    if d < d_min {
        warnings |= itm_warnings::WARN_PATH_DISTANCE_TOO_SMALL_1;
    }
    if d < 1e3 {
        warnings |= itm_warnings::WARN_PATH_DISTANCE_TOO_SMALL_2;
    }
    if d > 1000e3 {
        warnings |= itm_warnings::WARN_PATH_DISTANCE_TOO_BIG_1;
    }
    if d > 2000e3 {
        warnings |= itm_warnings::WARN_PATH_DISTANCE_TOO_BIG_2;
    }

    let mut a_ref: f64;
    let mut propmode = 0i32;

    if d < d_sml {
        // small path / line-of-sight region
        let a_sml_db = d_sml * m_d + a_d0;

        let mut d_0 = 0.04 * f * h_e[0] * h_e[1];
        let d_1: f64;
        if a_d0 >= 0.0 {
            d_0 = d_0.min(0.5 * d_ml);
            d_1 = d_0 + 0.25 * (d_ml - d_0);
        } else {
            d_1 = (-a_d0 / m_d).max(0.25 * d_ml);
        }

        let a_1_db = line_of_sight_loss(d_1, h_e, z_g, delta_h, m_d, a_d0, d_sml, f);

        let mut flag = false;
        let mut khat1 = 0.0f64;
        let mut khat2 = 0.0f64;

        if d_0 < d_1 {
            let a_0_db = line_of_sight_loss(d_0, h_e, z_g, delta_h, m_d, a_d0, d_sml, f);
            let q = (d_sml / d_0).ln();

            let denom = (d_sml - d_0) * (d_1 / d_0).ln() - (d_1 - d_0) * q;
            if denom.abs() > 0.0 {
                khat2 =
                    ((d_sml - d_0) * (a_1_db - a_0_db) - (d_1 - d_0) * (a_sml_db - a_0_db)) / denom;
            }
            if khat2.is_nan() || khat2 < 0.0 {
                khat2 = 0.0;
            }

            flag = a_d0 > 0.0 || khat2 > 0.0;

            if flag {
                khat1 = (a_sml_db - a_0_db - khat2 * q) / (d_sml - d_0);
                if khat1 < 0.0 {
                    khat1 = 0.0;
                    khat2 = (a_sml_db - a_0_db).abs() / q;
                    if khat2 == 0.0 {
                        khat1 = m_d;
                    }
                }
            }
        }

        if !flag {
            khat1 = (a_sml_db - a_1_db).abs() / (d_sml - d_1);
            khat2 = 0.0;
            if khat1 == 0.0 {
                khat1 = m_d;
            }
        }

        let a_o_db = a_sml_db - khat1 * d_sml - khat2 * d_sml.ln();
        a_ref = a_o_db + khat1 * d + khat2 * d.ln();
        propmode = 1; // LINE_OF_SIGHT
    } else {
        // trans-horizon path
        let d_5 = d_ml + 200e3;
        let d_6 = d_ml + 400e3;

        let mut h0 = -1.0f64;
        let a6_db = troposcatter_loss(
            d_6, &theta_hzn, &d_hzn, &h_e, a_e, n_s, f, theta_los, &mut h0,
        );
        let a5_db = troposcatter_loss(
            d_5, &theta_hzn, &d_hzn, &h_e, a_e, n_s, f, theta_los, &mut h0,
        );

        let (m_s, a_s0_db, d_x_m) = if a5_db < 1000.0 {
            let m_s = (a6_db - a5_db) / 200e3;
            let mut d_x = d_ml.max(d_ml + 1.088 * (a_e.powi(2) / f).powf(1.0 / 3.0) * f.ln().abs());
            d_x = d_x.max((a5_db - a_d0 - m_s * d_5) / (m_d - m_s));
            let a_s0 = (m_d - m_s) * d_x + a_d0;
            (m_s, a_s0, d_x)
        } else {
            (m_d, a_d0, 10e6)
        };

        if d > d_x_m {
            a_ref = m_s * d + a_s0_db;
            propmode = 3; // TROPOSCATTER
        } else {
            a_ref = m_d * d + a_d0;
            propmode = 2; // DIFFRACTION
        }
    }

    if a_ref < 0.0 {
        a_ref = 0.0;
    }

    Ok(LongleyRiceResult {
        a_ref,
        warnings,
        propmode,
    })
}
