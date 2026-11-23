//! The port against NTIA's current C++ ITM (`reference/ntia-itm`), built as this test's reference: NTIA's published
//! cases (`tests/fixtures/ntia`), then generated cases across ITM's accepted ranges, in both modes. A case agrees when
//! both succeed with losses within `TOLERANCE_DB` and equal warnings, or both refuse it.

use itm::math::itm::{itm_area_tls, itm_p2p_tls};
use ntia_itm_reference::{Reference, area_tls, p2p_tls};

/// Agreement asked of the port, decibels; NTIA rounds its published outputs to 0.01 dB.
const TOLERANCE_DB: f64 = 0.1;
/// Generated cases per mode.
const GENERATED: usize = 1000;

fn rows(csv: &str) -> Vec<Vec<f64>> {
    csv.lines()
        .skip(1)
        .filter(|line| !line.trim().is_empty())
        .map(|line| line.split(',').map(|value| value.trim().parse().expect("a number")).collect())
        .collect()
}

fn profiles(csv: &str) -> Vec<Vec<f64>> {
    csv.lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| line.split(',').filter(|value| !value.trim().is_empty()).map(|value| value.trim().parse().expect("a number")).collect())
        .collect()
}

/// A small deterministic generator (xorshift64*), so every run checks the same cases without a dependency.
struct Cases(u64);

impl Cases {
    fn unit(&mut self) -> f64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        (self.0.wrapping_mul(0x2545_F491_4F6C_DD1D) >> 11) as f64 / (1u64 << 53) as f64
    }
    fn range(&mut self, low: f64, high: f64) -> f64 {
        low + (high - low) * self.unit()
    }
    fn pick(&mut self, count: i32) -> i32 {
        (self.unit() * count as f64) as i32
    }
}

/// Whether the port's result agrees with NTIA's for one case; the reason when it does not.
fn compare(label: &str, port: Result<(f64, u32, itm::math::itm::IntermediateValues), itm::math::itm::ItmError>, reference: Reference) -> Option<String> {
    let reference_ok = reference.code == 0 || reference.code == 1;
    match port {
        Ok((loss, warnings, _)) if reference_ok => {
            let close = (loss - reference.loss_db).abs() <= TOLERANCE_DB;
            let same_warnings = i64::from(warnings) == reference.warnings;
            (!close || !same_warnings).then(|| {
                format!(
                    "{label}: port {loss:.3} dB w{warnings:#x}, NTIA {:.3} dB w{:#x}",
                    reference.loss_db, reference.warnings
                )
            })
        }
        Err(_) if !reference_ok => None,
        Ok((loss, _, _)) => Some(format!("{label}: port gave {loss:.3} dB, NTIA refused it (code {})", reference.code)),
        Err(error) => Some(format!("{label}: port refused it ({error:?}), NTIA gave {:.3} dB", reference.loss_db)),
    }
}

fn assert_agrees(failures: Vec<String>, what: &str) {
    assert!(failures.is_empty(), "{what}: {} case(s) disagree with NTIA:\n{}", failures.len(), failures.join("\n"));
}

#[test]
fn published_point_to_point_cases_match_ntia() {
    let cases = rows(include_str!("fixtures/ntia/p2p.csv"));
    let pfls = profiles(include_str!("fixtures/ntia/pfls.csv"));
    assert_eq!(cases.len(), pfls.len(), "one profile per case");
    let mut failures = Vec::new();
    for (index, (c, pfl)) in cases.iter().zip(&pfls).enumerate() {
        // h_tx, h_rx, epsilon, sigma, N_0, f, pol, climate, time, location, situation, mdvar, A
        let port = itm_p2p_tls(c[0], c[1], pfl, c[7] as i32, c[4], c[5], c[6] as i32, c[2], c[3], c[11] as i32, c[8], c[9], c[10]);
        let reference = p2p_tls(c[0], c[1], pfl, c[7] as i32, c[4], c[5], c[6] as i32, c[2], c[3], c[11] as i32, c[8], c[9], c[10]);
        assert!((reference.loss_db - c[12]).abs() <= 0.05, "case {index}: NTIA's C++ gives {} for a published {}", reference.loss_db, c[12]);
        failures.extend(compare(&format!("published p2p {index}"), port, reference));
    }
    assert_agrees(failures, "published point-to-point");
}

#[test]
fn published_area_cases_match_ntia() {
    let mut failures = Vec::new();
    for (index, c) in rows(include_str!("fixtures/ntia/area.csv")).iter().enumerate() {
        // h_tx, h_rx, delta_h, mdvar, d, tx_siting, rx_siting, epsilon, sigma, N_0, f, pol, climate, time, location,
        // situation, A
        let port = itm_area_tls(
            c[0], c[1], c[5] as i32, c[6] as i32, c[4], c[2], c[12] as i32, c[9], c[10], c[11] as i32, c[7], c[8], c[3] as i32, c[13], c[14], c[15],
        );
        let reference = area_tls(
            c[0], c[1], c[5] as i32, c[6] as i32, c[4], c[2], c[12] as i32, c[9], c[10], c[11] as i32, c[7], c[8], c[3] as i32, c[13], c[14], c[15],
        );
        assert!((reference.loss_db - c[16]).abs() <= 0.05, "case {index}: NTIA's C++ gives {} for a published {}", reference.loss_db, c[16]);
        failures.extend(compare(&format!("published area {index}"), port, reference));
    }
    assert_agrees(failures, "published area");
}

/// A terrain profile in ITM's form: hills and valleys from a few sine waves, 2 to 1,000 intervals, 10 m to 1 km apart.
fn generated_profile(cases: &mut Cases) -> Vec<f64> {
    let intervals = 2 + cases.pick(999) as usize;
    let spacing = cases.range(10.0, 1000.0);
    let base = cases.range(0.0, 1500.0);
    let waves: Vec<(f64, f64, f64)> = (0..3).map(|_| (cases.range(0.0, 300.0), cases.range(0.2, 6.0), cases.range(0.0, 6.3))).collect();
    let mut pfl = vec![intervals as f64, spacing];
    for i in 0..=intervals {
        let x = i as f64 / intervals as f64;
        pfl.push(base + waves.iter().map(|(height, cycles, phase)| height * (x * cycles * std::f64::consts::TAU + phase).sin()).sum::<f64>());
    }
    pfl
}

#[test]
fn generated_point_to_point_cases_match_ntia() {
    let mut cases = Cases(0x1717_2023_0F0E_5EED);
    let mut failures = Vec::new();
    for index in 0..GENERATED {
        let pfl = generated_profile(&mut cases);
        let (h_tx, h_rx) = (cases.range(0.5, 3000.0_f64.min(200.0)), cases.range(0.5, 200.0));
        let f_mhz = cases.range(20.0, 20_000.0);
        let climate = 1 + cases.pick(7);
        let pol = cases.pick(2);
        let (epsilon, sigma) = (cases.range(4.0, 80.0), cases.range(0.001, 5.0));
        let n_0 = cases.range(250.0, 400.0);
        let mdvar = cases.pick(4) + [0, 10, 20, 30][cases.pick(4) as usize];
        let (time, location, situation) = (cases.range(1.0, 99.0), cases.range(1.0, 99.0), cases.range(1.0, 99.0));
        let port = itm_p2p_tls(h_tx, h_rx, &pfl, climate, n_0, f_mhz, pol, epsilon, sigma, mdvar, time, location, situation);
        let reference = p2p_tls(h_tx, h_rx, &pfl, climate, n_0, f_mhz, pol, epsilon, sigma, mdvar, time, location, situation);
        failures.extend(compare(&format!("generated p2p {index}"), port, reference));
    }
    assert_agrees(failures, "generated point-to-point");
}

#[test]
fn generated_area_cases_match_ntia() {
    let mut cases = Cases(0x0A2E_A5EE_D000_1717);
    let mut failures = Vec::new();
    for index in 0..GENERATED {
        let (h_tx, h_rx) = (cases.range(0.5, 200.0), cases.range(0.5, 200.0));
        let (tx_site, rx_site) = (cases.pick(3), cases.pick(3));
        let d_km = cases.range(1.0, 1000.0);
        let delta_h = cases.range(0.0, 500.0);
        let climate = 1 + cases.pick(7);
        let n_0 = cases.range(250.0, 400.0);
        let f_mhz = cases.range(20.0, 20_000.0);
        let pol = cases.pick(2);
        let (epsilon, sigma) = (cases.range(4.0, 80.0), cases.range(0.001, 5.0));
        let mdvar = cases.pick(4) + [0, 10, 20, 30][cases.pick(4) as usize];
        let (time, location, situation) = (cases.range(1.0, 99.0), cases.range(1.0, 99.0), cases.range(1.0, 99.0));
        let port = itm_area_tls(
            h_tx, h_rx, tx_site, rx_site, d_km, delta_h, climate, n_0, f_mhz, pol, epsilon, sigma, mdvar, time, location, situation,
        );
        let reference = area_tls(
            h_tx, h_rx, tx_site, rx_site, d_km, delta_h, climate, n_0, f_mhz, pol, epsilon, sigma, mdvar, time, location, situation,
        );
        failures.extend(compare(&format!("generated area {index}"), port, reference));
    }
    assert_agrees(failures, "generated area");
}
