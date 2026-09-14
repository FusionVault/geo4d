//! Probability of collision for a close approach — the short-term-encounter model with Chan's
//! analytic series.
//!
//! A miss distance alone does not say how dangerous a conjunction is; that depends on how well each
//! object's position is known. Given a position [`Covariance3`] for each object and a combined
//! hard-body radius, [`collision_probability`] reduces the 3D encounter to a 2D Gaussian integral
//! over a disk in the plane perpendicular to the relative velocity (the "B-plane") and evaluates it
//! with F. K. Chan's convergent series — the method operational conjunction assessment uses on
//! CCSDS Conjunction Data Messages.
//!
//! Assumptions (the short-term-encounter model): relative motion is rectilinear through the
//! encounter, position errors are Gaussian, velocity uncertainty is negligible, and the objects are
//! spheres. Covariances are taken in the same ECEF frame as the state vectors and treated as
//! constant across the brief encounter.

use crate::kinematics::closest_approach;
use crate::types::StateVector;
use geo3d::Vec3;

/// A 3×3 position covariance (metres²), symmetric. Add two together for a combined covariance.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Covariance3 {
    /// The symmetric covariance matrix, row-major, metres².
    pub m: [[f64; 3]; 3],
}

impl Covariance3 {
    /// From a full 3×3 matrix (assumed symmetric), metres².
    #[inline]
    pub const fn from_matrix(m: [[f64; 3]; 3]) -> Self {
        Covariance3 { m }
    }

    /// An isotropic covariance: the same 1σ position uncertainty (metres) in every direction.
    #[inline]
    pub fn isotropic(sigma_m: f64) -> Self {
        let v = sigma_m * sigma_m;
        Covariance3::from_matrix([[v, 0.0, 0.0], [0.0, v, 0.0], [0.0, 0.0, v]])
    }

    /// A diagonal covariance from three 1σ position uncertainties (metres) along the ECEF axes.
    #[inline]
    pub fn diagonal(sigma_x_m: f64, sigma_y_m: f64, sigma_z_m: f64) -> Self {
        Covariance3::from_matrix([
            [sigma_x_m * sigma_x_m, 0.0, 0.0],
            [0.0, sigma_y_m * sigma_y_m, 0.0],
            [0.0, 0.0, sigma_z_m * sigma_z_m],
        ])
    }

    /// The quadratic form `xᵀ · C · y`.
    fn quad(&self, x: Vec3, y: Vec3) -> f64 {
        let cy = Vec3::new(
            self.m[0][0] * y.x + self.m[0][1] * y.y + self.m[0][2] * y.z,
            self.m[1][0] * y.x + self.m[1][1] * y.y + self.m[1][2] * y.z,
            self.m[2][0] * y.x + self.m[2][1] * y.y + self.m[2][2] * y.z,
        );
        x.dot(cy)
    }
}

impl core::ops::Add for Covariance3 {
    type Output = Covariance3;
    /// The combined uncertainty of two independent objects — the element-wise sum.
    fn add(self, o: Covariance3) -> Covariance3 {
        let mut m = self.m;
        for (row, orow) in m.iter_mut().zip(o.m.iter()) {
            for (a, b) in row.iter_mut().zip(orow.iter()) {
                *a += *b;
            }
        }
        Covariance3 { m }
    }
}

/// The probability that two objects collide at their close approach, in `[0, 1]`. `cov_a` and
/// `cov_b` are each object's position covariance (metres², ECEF); `hard_body_radius_m` is the sum of
/// the two object radii. Uses Chan's series on the 2D encounter-plane integral.
///
/// Zero relative velocity has no encounter plane; there the result is a deterministic hit test
/// (`1.0` if the objects overlap at closest approach, else `0.0`).
///
/// ```
/// use geo4d::{collision_probability, Covariance3, Ecef, Epoch, StateVector, Vec3};
/// let t = Epoch::from_unix_seconds(1_700_000_000.0);
/// // Two objects on an exact head-on collision course (zero miss), 10 m 1σ each way, 5 m combined size.
/// let a = StateVector::new(t, Ecef::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
/// let b = StateVector::new(t, Ecef::new(100.0, 0.0, 0.0), Vec3::new(-1.0, 0.0, 0.0));
/// let pc = collision_probability(a, Covariance3::isotropic(10.0), b, Covariance3::isotropic(10.0), 5.0);
/// // Zero-miss, isotropic combined σ² = 200: Pc = 1 − exp(−R²/(2σ²)) = 1 − exp(−25/400) ≈ 0.0606.
/// assert!((pc - 0.0606).abs() < 1e-3);
/// ```
pub fn collision_probability(
    a: StateVector,
    cov_a: Covariance3,
    b: StateVector,
    cov_b: Covariance3,
    hard_body_radius_m: f64,
) -> f64 {
    let tca = closest_approach(a, b).tca;
    let miss = a.at(tca).position.vec() - b.at(tca).position.vec();
    let v_rel = a.velocity - b.velocity;

    let zhat = match v_rel.normalized() {
        Some(z) => z,
        // No relative motion: no encounter plane. Deterministic overlap test.
        None => return if miss.norm() < hard_body_radius_m { 1.0 } else { 0.0 },
    };
    // In-plane axes: x̂ along the projected miss (or an arbitrary perpendicular if the miss is along ẑ).
    let xhat = (miss - zhat * miss.dot(zhat))
        .normalized()
        .unwrap_or_else(|| any_perpendicular(zhat));
    let yhat = zhat.cross(xhat);

    // Project the combined covariance and the miss onto the plane.
    let cov = cov_a + cov_b;
    let (a2, b2, d2) = (cov.quad(xhat, xhat), cov.quad(xhat, yhat), cov.quad(yhat, yhat));
    let (mx, my) = (miss.dot(xhat), miss.dot(yhat));

    // Diagonalise the 2×2 [[a2,b2],[b2,d2]] and rotate the miss into its eigenframe.
    let mean = (a2 + d2) / 2.0;
    let radius = (((a2 - d2) / 2.0).powi(2) + b2 * b2).sqrt();
    let (lam1, lam2) = (mean + radius, mean - radius);
    if lam1 <= 0.0 || lam2 <= 0.0 {
        // Degenerate covariance: fall back to a deterministic test.
        return if miss.norm() < hard_body_radius_m { 1.0 } else { 0.0 };
    }
    let (s, c) = (0.5 * (2.0 * b2).atan2(a2 - d2)).sin_cos();
    let (xm, ym) = (mx * c + my * s, -mx * s + my * c);

    let u = hard_body_radius_m * hard_body_radius_m / (lam1 * lam2).sqrt();
    let v = xm * xm / lam1 + ym * ym / lam2;
    chan_series(u, v)
}

/// Chan's series for the offset-circle Gaussian integral:
/// `Pc = e^{-v/2} Σ_m (vᵐ/2ᵐm!)·[1 − e^{-u/2} Σ_{k≤m} uᵏ/2ᵏk!]`.
fn chan_series(u: f64, v: f64) -> f64 {
    let eu = (-u / 2.0).exp();
    let mut outer = 1.0; // vᵐ / (2ᵐ m!), m = 0
    let mut inner = 1.0; // Σ_{k=0}^{m} uᵏ/(2ᵏ k!), m = 0
    let mut inner_term = 1.0; // uᵏ/(2ᵏ k!), k = 0
    let mut sum = 0.0;
    for m in 0..1000usize {
        let contrib = outer * (1.0 - eu * inner);
        sum += contrib;
        if m > 2 && contrib < 1e-15 * sum.max(1e-300) {
            break;
        }
        let mp1 = 2.0 * (m as f64 + 1.0);
        outer *= v / mp1;
        inner_term *= u / mp1;
        inner += inner_term;
    }
    ((-v / 2.0).exp() * sum).clamp(0.0, 1.0)
}

/// Any unit vector perpendicular to `z` (which must be unit length).
fn any_perpendicular(z: Vec3) -> Vec3 {
    let axis = if z.x.abs() <= z.y.abs() && z.x.abs() <= z.z.abs() {
        Vec3::new(1.0, 0.0, 0.0)
    } else if z.y.abs() <= z.z.abs() {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    };
    z.cross(axis).normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Ecef, Epoch};

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn zero_miss_matches_the_rayleigh_closed_form() {
        // Head-on, zero miss: Pc = 1 − exp(−R²/(2σ_combined²)).
        let t = Epoch::from_unix_seconds(1_700_000_000.0);
        let a = StateVector::new(t, Ecef::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        let b = StateVector::new(t, Ecef::new(100.0, 0.0, 0.0), Vec3::new(-1.0, 0.0, 0.0));
        for (sig, r) in [(10.0_f64, 5.0_f64), (25.0, 10.0), (50.0, 20.0)] {
            let comb = 2.0 * sig * sig; // σ_a² + σ_b², both isotropic = sig
            let expected = 1.0 - (-(r * r) / (2.0 * comb)).exp();
            let pc = collision_probability(a, Covariance3::isotropic(sig), b, Covariance3::isotropic(sig), r);
            assert!(close(pc, expected, 1e-6), "σ={sig} R={r}: {pc} vs {expected}");
        }
    }

    #[test]
    fn probability_falls_with_miss_distance_and_stays_bounded() {
        let t = Epoch::from_unix_seconds(1_700_000_000.0);
        let a = StateVector::new(t, Ecef::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        let mut last = 1.0;
        for miss in [0.0, 5.0, 10.0, 30.0, 100.0] {
            // b passes offset in y (perpendicular to the x-axis closing motion) by `miss`.
            let b = StateVector::new(t, Ecef::new(100.0, miss, 0.0), Vec3::new(-1.0, 0.0, 0.0));
            let pc = collision_probability(a, Covariance3::isotropic(10.0), b, Covariance3::isotropic(10.0), 5.0);
            assert!((0.0..=1.0).contains(&pc));
            assert!(pc <= last + 1e-12, "Pc should not increase with miss: {pc} > {last}");
            last = pc;
        }
        assert!(last < 1e-3); // far miss → negligible
    }

    #[test]
    fn zero_relative_velocity_is_a_deterministic_hit_test() {
        let t = Epoch::from_unix_seconds(1_700_000_000.0);
        let a = StateVector::new(t, Ecef::new(0.0, 0.0, 0.0), Vec3::new(5.0, 0.0, 0.0));
        let touching = StateVector::new(t, Ecef::new(3.0, 0.0, 0.0), Vec3::new(5.0, 0.0, 0.0));
        assert_eq!(
            collision_probability(
                a,
                Covariance3::isotropic(1.0),
                touching,
                Covariance3::isotropic(1.0),
                5.0
            ),
            1.0
        );
        let apart = StateVector::new(t, Ecef::new(30.0, 0.0, 0.0), Vec3::new(5.0, 0.0, 0.0));
        assert_eq!(
            collision_probability(a, Covariance3::isotropic(1.0), apart, Covariance3::isotropic(1.0), 5.0),
            0.0
        );
    }
}
