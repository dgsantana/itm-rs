//! Bindings to NTIA's C++ ITM, the reference the library's tests compare with. Test-only.

use std::ffi::{c_int, c_long};

/// NTIA's intermediate values, laid out as `IntermediateValues` in `include/itm.h`.
#[repr(C)]
#[derive(Debug, Default, Clone, Copy)]
pub struct Intermediate {
    /// Terminal horizon angles, radians.
    pub theta_hzn: [f64; 2],
    /// Terminal horizon distances, metres.
    pub d_hzn_meter: [f64; 2],
    /// Terminal effective heights, metres.
    pub h_e_meter: [f64; 2],
    /// Surface refractivity, N-units.
    pub n_s: f64,
    /// Terrain irregularity, metres.
    pub delta_h_meter: f64,
    /// Reference attenuation, dB.
    pub a_ref_db: f64,
    /// Free-space loss, dB.
    pub a_fs_db: f64,
    /// Path distance, km.
    pub d_km: f64,
    /// Mode of propagation (1 line of sight, 2 diffraction, 3 troposcatter).
    pub mode: c_int,
}

/// One reference result: NTIA's return code (0 success, 1 success with warnings, else an error), the loss and the
/// warning bits.
#[derive(Debug, Clone, Copy)]
pub struct Reference {
    /// NTIA's return code.
    pub code: i32,
    /// Basic transmission loss, dB.
    pub loss_db: f64,
    /// Warning bits.
    pub warnings: i64,
    /// Intermediate values.
    pub intermediate: Intermediate,
}

unsafe extern "C" {
    #[link_name = "ITM_P2P_TLS_Ex"]
    fn itm_p2p_tls_ex(
        h_tx: f64,
        h_rx: f64,
        pfl: *const f64,
        climate: c_int,
        n_0: f64,
        f_mhz: f64,
        pol: c_int,
        epsilon: f64,
        sigma: f64,
        mdvar: c_int,
        time: f64,
        location: f64,
        situation: f64,
        a_db: *mut f64,
        warnings: *mut c_long,
        intermediate: *mut Intermediate,
    ) -> c_int;

    #[link_name = "ITM_AREA_TLS_Ex"]
    fn itm_area_tls_ex(
        h_tx: f64,
        h_rx: f64,
        tx_site: c_int,
        rx_site: c_int,
        d_km: f64,
        delta_h: f64,
        climate: c_int,
        n_0: f64,
        f_mhz: f64,
        pol: c_int,
        epsilon: f64,
        sigma: f64,
        mdvar: c_int,
        time: f64,
        location: f64,
        situation: f64,
        a_db: *mut f64,
        warnings: *mut c_long,
        intermediate: *mut Intermediate,
    ) -> c_int;
}

/// NTIA's point-to-point time, location and situation mode. `pfl` is ITM's profile: point count less one, spacing in
/// metres, then the heights.
#[allow(clippy::too_many_arguments)]
pub fn p2p_tls(
    h_tx: f64,
    h_rx: f64,
    pfl: &[f64],
    climate: i32,
    n_0: f64,
    f_mhz: f64,
    pol: i32,
    epsilon: f64,
    sigma: f64,
    mdvar: i32,
    time: f64,
    location: f64,
    situation: f64,
) -> Reference {
    // NTIA reads pfl[0] + 3 values; a shorter slice would be read past its end.
    assert!(
        pfl.len() >= 2 && pfl[0] >= 0.0 && pfl.len() >= pfl[0] as usize + 3,
        "a profile of {} values cannot hold the {} points its first value declares",
        pfl.len(),
        pfl.first().copied().unwrap_or(0.0) + 1.0
    );
    let (mut loss_db, mut warnings, mut intermediate) = (0.0, 0 as c_long, Intermediate::default());
    // SAFETY: the assertion above makes `pfl` hold every value NTIA reads; the out pointers are valid locals.
    let code = unsafe {
        itm_p2p_tls_ex(
            h_tx, h_rx, pfl.as_ptr(), climate, n_0, f_mhz, pol, epsilon, sigma, mdvar, time, location, situation, &mut loss_db,
            &mut warnings, &mut intermediate,
        )
    };
    Reference { code, loss_db, warnings: warnings as i64, intermediate }
}

/// NTIA's area time, location and situation mode.
#[allow(clippy::too_many_arguments)]
pub fn area_tls(
    h_tx: f64,
    h_rx: f64,
    tx_site: i32,
    rx_site: i32,
    d_km: f64,
    delta_h: f64,
    climate: i32,
    n_0: f64,
    f_mhz: f64,
    pol: i32,
    epsilon: f64,
    sigma: f64,
    mdvar: i32,
    time: f64,
    location: f64,
    situation: f64,
) -> Reference {
    let (mut loss_db, mut warnings, mut intermediate) = (0.0, 0 as c_long, Intermediate::default());
    // SAFETY: plain values in, valid locals out.
    let code = unsafe {
        itm_area_tls_ex(
            h_tx, h_rx, tx_site, rx_site, d_km, delta_h, climate, n_0, f_mhz, pol, epsilon, sigma, mdvar, time, location, situation,
            &mut loss_db, &mut warnings, &mut intermediate,
        )
    };
    Reference { code, loss_db, warnings: warnings as i64, intermediate }
}
