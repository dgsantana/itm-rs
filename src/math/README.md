# Math Module

Mathematical functions for the ITM (Irregular Terrain Model) radio propagation calculations.

## Module Structure

```
math/
├── mod.rs              # Module definitions and re-exports
├── propagation.rs      # Free space loss calculations
├── diffraction.rs      # Diffraction modeling
├── scatter.rs          # Tropospheric scatter propagation
├── height_gain.rs      # Height gain calculations
├── terrain.rs          # Terrain profile analysis
├── statistics.rs       # Statistical functions and distributions
└── variability.rs      # Variability loss calculations
```

## Quick Start

### Using Re-exported Functions (Recommended)

```rust
use itm_rs::math::{free_space_loss, fresnel_integral, h0_function, terrain_roughness, 
                     f_function, troposcatter_loss, inverse_ccdf};

// Calculate free space loss
let loss = free_space_loss(1000.0, 1e9); // 1km at 1 GHz

// Calculate diffraction loss
let diffraction = fresnel_integral(5.0);

// Calculate height gain
let gain = h0_function(10.0, 2.5);

// Calculate terrain roughness
let roughness = terrain_roughness(10_000.0, 100.0); // 10km path, 100m elevation variation

// Calculate troposcatter attenuation
let f = f_function(50_000.0); // F-function for 50km

// Calculate troposcatter loss
let theta_hzn = [0.01, 0.01];
let d_hzn = [50_000.0, 50_000.0];
let h_e = [100.0, 100.0];
let mut h0 = 0.0;
let scatter_loss = troposcatter_loss(100_000.0, &theta_hzn, &d_hzn, &h_e, 
                                      8.5e6, 301.0, 1000.0, 0.001, &mut h0);

// Calculate inverse CCDF for 90% reliability (q = 0.1)
let margin = inverse_ccdf(0.1); // Returns ~1.28 std deviations
```

### Using Submodules Directly

```rust
use itm_rs::math::propagation;
use itm_rs::math::diffraction;
use itm_rs::math::height_gain;

let loss = propagation::free_space_loss(1000.0, 1e9);
let diffraction = diffraction::fresnel_integral(5.0);
let gain = height_gain::h0_function(10.0, 2.5);
```

## Modules

### `propagation`

Basic radio propagation calculations.

**Functions:**
- `free_space_loss(d: f64, f: f64) -> f64` - Compute Friis free space path loss
- `initialize_point_to_point(...) -> (Complex<f64>, f64, f64)` - Initialize ground impedance, effective earth curvature, and surface refractivity

**Types:**
- `Polarization` - Wave polarization (horizontal or vertical)

### `diffraction`

Diffraction modeling for obstacles and terrain.

**Functions:**
- `fresnel_integral(v2: f64) -> f64` - Fresnel integral approximation for diffraction loss
- `smooth_earth_diffraction(...) -> f64` - Smooth-earth diffraction loss using Vogler 3-radii

### `scatter`

Tropospheric scatter propagation modeling.

**Functions:**
- `f_function(theta_distance: f64) -> f64` - Attenuation function F(θ·d) for scatter propagation
- `troposcatter_loss(...)` - Complete troposcatter loss calculation including height gain effects

### `height_gain`

Height gain factor calculations for different terrain roughness categories.

**Functions:**
- `h0_curve(j: usize, r: f64) -> f64` - Curve fit for height gain parameters
- `h0_function(r: f64, eta_s: f64) -> f64` - Interpolated height gain calculation

### `terrain`

Terrain profile analysis and statistical methods.

**Functions:**
- `terrain_roughness(d: f64, delta_h: f64) -> f64` - Calculate terrain roughness parameter from distance and elevation variation
- `initialize_area(site_criteria: [SitingCriteria; 2], gamma_e: f64, delta_h_m: f64, h_m: [f64; 2]) -> ([f64; 2], [f64; 2], [f64; 2])` - Initialize area mode parameters
- `find_horizons(pfl: &[f64], a_e: f64, h_m: [f64; 2]) -> ([f64; 2], [f64; 2])` - Compute terminal radio horizon angles and distances
- `compute_delta_h(pfl: &[f64], d_start_m: f64, d_end_m: f64) -> f64` - Compute terrain irregularity parameter
- Internal functions for terrain analysis (crate-private)

**Types:**
- `SitingCriteria` - Terminal siting criteria (random, careful, very careful)

### `statistics`

Statistical functions and probability distributions.

**Functions:**
- `inverse_ccdf(q: f64) -> f64` - Inverse complementary cumulative distribution function (quantile function)
- `inverse_complementary_cumulative_distribution_function(q: f64) -> f64` - Alias with full name

### `variability`

Variability loss calculations.

**Functions:**
- `curve(c1: f64, c2: f64, x1: f64, x2: f64, x3: f64, d_e_m: f64) -> f64` - Curve helper (TN101v2)
- `variability_loss(...) -> f64` - Compute statistical variability loss

**Types:**
- `Climate` - Radio climate classification
- `VariabilityMode` - Variability mode (broadcast, mobile, accidental, single message)

## Documentation

- See [`REFERENCES.md`](../../REFERENCES.md) for academic citations and equation sources
- See [`MODULE_ORGANIZATION.md`](../../MODULE_ORGANIZATION.md) for detailed organization rationale

## Testing

Run tests for the math module:

```bash
cargo test --lib math
```

Run tests for a specific submodule:

```bash
cargo test --lib math::propagation
```

## Contributing

When adding new functions:

1. **Choose the appropriate submodule** based on functionality
2. **Add comprehensive documentation** with:
   - Description
   - Arguments
   - Returns
   - Examples
   - References (academic papers, equations)
3. **Include unit tests**
4. **Update REFERENCES.md** with citations

### Example Function Template

```rust
/// Brief one-line description.
///
/// Detailed description of what the function does and when to use it.
///
/// # Arguments
///
/// * `param1` - Description of parameter 1
/// * `param2` - Description of parameter 2
///
/// # Returns
///
/// Description of return value
///
/// # References
///
/// - Author (Year). "Paper Title". Journal/Conference.
/// - ITM TN101, Equation X.Y: Description
///
/// # Examples
///
/// ```
/// use itm_rs::math::module::function;
///
/// let result = function(arg1, arg2);
/// assert!(result > 0.0);
/// ```
pub fn function(param1: f64, param2: f64) -> f64 {
    // Implementation
}
```

## License

See LICENSE file in the project root.
