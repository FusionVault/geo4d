//! 4D kinematics: closest approach, conflict windows and intercept — all in absolute time.
//!
//! `geo3d`'s [`closest_approach`](geo3d::closest_approach) works in *relative* time ("seconds from
//! now"), which is enough when both objects are sampled at the same instant. Here every object is a
//! [`StateVector`] carrying its own [`Epoch`], so two tracks fixed at *different* times can still be
//! compared: the answers come back as absolute epochs and [`TimeWindow`]s.
//!
//! All motion is constant-velocity (first order) — the standard assumption behind every closed form
//! below.

use crate::time::{Duration, Epoch, TimeWindow};
use crate::types::{Ecef4, StateVector};
use geo3d::{Ecef, LocalFrame, Vec3, WGS84};

/// The closest approach of two moving objects: when it happens and how near they pass.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Conjunction {
    /// The epoch of closest approach (may be in the past for diverging objects).
    pub tca: Epoch,
    /// The separation at that epoch, metres.
    pub distance_m: f64,
    /// The relative speed of the two objects, m/s.
    pub relative_speed_mps: f64,
}

/// The closest point of approach of two state vectors, in absolute time. Each is propagated to a
/// common epoch first, so they need not be sampled simultaneously. A first-order estimate: compare
/// [`distance_m`](Conjunction::distance_m) with a separation threshold.
///
/// ```
/// use geo4d::{closest_approach, Ecef, Epoch, StateVector, Vec3};
/// let t = Epoch::from_unix_seconds(1_700_000_000.0);
/// // Two objects closing along x, offset 5 m in y: closest 5 s from t, 5 m apart.
/// let a = StateVector::new(t, Ecef::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
/// let b = StateVector::new(t, Ecef::new(10.0, 5.0, 0.0), Vec3::new(-1.0, 0.0, 0.0));
/// let c = closest_approach(a, b);
/// assert!((c.distance_m - 5.0).abs() < 1e-9);
/// assert!((c.tca - t).as_seconds().abs() < 5.0 + 1e-9);
/// ```
pub fn closest_approach(a: StateVector, b: StateVector) -> Conjunction {
    let t_ref = a.epoch;
    let b_ref = b.at(t_ref);
    let prel = a.position.vec() - b_ref.position.vec();
    let vrel = a.velocity - b.velocity;
    let vv = vrel.dot(vrel);
    let tstar = if vv == 0.0 { 0.0 } else { -prel.dot(vrel) / vv };
    Conjunction {
        tca: t_ref + Duration::from_seconds(tstar),
        distance_m: (prel + vrel * tstar).norm(),
        relative_speed_mps: vv.sqrt(),
    }
}

/// The closest approach restricted to a time window — the nearest they pass **while the window is
/// open**. Because separation is convex in time, this is the unconstrained TCA clamped into the
/// window (no second solve).
pub fn closest_approach_within(a: StateVector, b: StateVector, window: TimeWindow) -> Conjunction {
    let free = closest_approach(a, b);
    let tca = window.clamp(free.tca);
    Conjunction {
        tca,
        distance_m: a.at(tca).position.distance_to(b.at(tca).position),
        relative_speed_mps: free.relative_speed_mps,
    }
}

/// A cylindrical protected zone: a horizontal separation minimum and a vertical separation minimum,
/// the shape aviation and drone-traffic (UTM) separation standards use.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ProtectedZone {
    /// Minimum horizontal separation, metres.
    pub horizontal_m: f64,
    /// Minimum vertical separation, metres.
    pub vertical_m: f64,
}

impl ProtectedZone {
    /// A protected zone from a horizontal and a vertical minimum, metres.
    #[inline]
    pub const fn new(horizontal_m: f64, vertical_m: f64) -> Self {
        ProtectedZone {
            horizontal_m,
            vertical_m,
        }
    }
}

/// A predicted loss of separation: the interval during which two objects are simultaneously within
/// the horizontal *and* vertical minima of a [`ProtectedZone`].
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Conflict {
    /// When separation is first lost.
    pub enters: Epoch,
    /// When separation is regained.
    pub exits: Epoch,
    /// The horizontal closest approach within the conflict, its epoch.
    pub tca: Epoch,
}

impl Conflict {
    /// The duration of the loss of separation.
    #[inline]
    pub fn duration(&self) -> Duration {
        self.exits - self.enters
    }
}

/// Whether and when two objects violate a cylindrical [`ProtectedZone`], searched over `window`
/// (pass [`TimeWindow::unbounded`] for all time). Horizontal and vertical are separated in a local
/// East-North-Up frame at the first object's position — an approximation good for the minutes-scale
/// look-ahead of conflict detection, where the frame's rotation is negligible.
///
/// Returns `None` if there is no conflict in the window (or if the conflict is unbounded because
/// there is no relative motion in some dimension and the window is open — pass a bounded window).
///
/// ```
/// use geo4d::{conflict, Ecef, Epoch, ProtectedZone, StateVector, TimeWindow, Vec3};
/// let t = Epoch::from_unix_seconds(1_700_000_000.0);
/// let frame = geo4d::Geodetic4::from_lat_lon(0.0, 0.0, 3_000.0, t).to_ecef4().position;
/// // Two aircraft converging head-on along the ECEF y-axis at 100 m/s each, 20 km apart.
/// let a = StateVector::new(t, frame, Vec3::new(0.0, 100.0, 0.0));
/// let b = StateVector::new(t, Ecef::new(frame.x, frame.y + 20_000.0, frame.z), Vec3::new(0.0, -100.0, 0.0));
/// // A 5 NM / 1000 ft protected zone is breached around the midpoint.
/// let c = conflict(a, b, ProtectedZone::new(9_260.0, 300.0), TimeWindow::unbounded()).unwrap();
/// assert!(c.enters < c.tca && c.tca < c.exits);
/// ```
pub fn conflict(a: StateVector, b: StateVector, zone: ProtectedZone, window: TimeWindow) -> Option<Conflict> {
    let t_ref = a.epoch;
    let b_ref = b.at(t_ref);
    let frame = LocalFrame::new(WGS84.to_geodetic(a.position));
    let p = frame.dir_to_enu(a.position.vec() - b_ref.position.vec());
    let v = frame.dir_to_enu(a.velocity - b.velocity);

    // Horizontal: ‖p_h + v_h·τ‖² < H² (a quadratic-in-τ interval); vertical: |p_z + v_z·τ| < V.
    let (ph, vh) = ((p.east, p.north), (v.east, v.north));
    let horiz = interval_below(
        vh.0 * vh.0 + vh.1 * vh.1,
        2.0 * (ph.0 * vh.0 + ph.1 * vh.1),
        ph.0 * ph.0 + ph.1 * ph.1 - zone.horizontal_m * zone.horizontal_m,
    )?;
    let vert = interval_abs_below(p.up, v.up, zone.vertical_m)?;

    let start_rel = window.start.map_or(f64::NEG_INFINITY, |s| (s - t_ref).as_seconds());
    let end_rel = window.end.map_or(f64::INFINITY, |e| (e - t_ref).as_seconds());
    let t_in = horiz.0.max(vert.0).max(start_rel);
    let t_out = horiz.1.min(vert.1).min(end_rel);
    if !t_in.is_finite() || !t_out.is_finite() || t_in >= t_out {
        return None;
    }
    // Worst horizontal moment inside the conflict.
    let denom = vh.0 * vh.0 + vh.1 * vh.1;
    let h_tca = if denom == 0.0 {
        t_in
    } else {
        (-(ph.0 * vh.0 + ph.1 * vh.1) / denom).clamp(t_in, t_out)
    };
    Some(Conflict {
        enters: t_ref + Duration::from_seconds(t_in),
        exits: t_ref + Duration::from_seconds(t_out),
        tca: t_ref + Duration::from_seconds(h_tca),
    })
}

/// The τ-interval where `a·τ² + b·τ + c < 0`. `None` if it is never negative; unbounded ends come
/// back as ±∞.
fn interval_below(a: f64, b: f64, c: f64) -> Option<(f64, f64)> {
    if a.abs() < 1e-12 {
        // Linear/constant: b·τ + c < 0.
        return if b.abs() < 1e-12 {
            (c < 0.0).then_some((f64::NEG_INFINITY, f64::INFINITY))
        } else {
            let r = -c / b;
            Some(if b > 0.0 {
                (f64::NEG_INFINITY, r)
            } else {
                (r, f64::INFINITY)
            })
        };
    }
    let disc = b * b - 4.0 * a * c;
    if disc <= 0.0 {
        return None; // parabola never dips below zero (a>0), or opens down (rare here)
    }
    let sq = disc.sqrt();
    let (r1, r2) = ((-b - sq) / (2.0 * a), (-b + sq) / (2.0 * a));
    Some((r1.min(r2), r1.max(r2)))
}

/// The τ-interval where `|p + v·τ| < half`. `None` if never; unbounded ends come back as ±∞.
fn interval_abs_below(p: f64, v: f64, half: f64) -> Option<(f64, f64)> {
    if v.abs() < 1e-12 {
        return (p.abs() < half).then_some((f64::NEG_INFINITY, f64::INFINITY));
    }
    let (ta, tb) = ((half - p) / v, (-half - p) / v);
    Some((ta.min(tb), ta.max(tb)))
}

/// A solved intercept: when and where a moving target is reached, and the velocity to do it with.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Intercept {
    /// The epoch of intercept.
    pub epoch: Epoch,
    /// Where the intercept happens, ECEF metres.
    pub point: Ecef,
    /// The constant velocity the interceptor must fly (its speed equals the requested speed), ECEF m/s.
    pub velocity: Vec3,
    /// Time from launch to intercept.
    pub time_to_go: Duration,
}

/// The earliest fixed-speed intercept of a moving `target` by an interceptor launching from `origin`
/// (at `origin`'s epoch) with the given constant `speed_mps`, heading free — the "lead collision"
/// solution. `None` if the target cannot be caught at that speed (it is faster and opening, or the
/// geometry is unreachable).
///
/// ```
/// use geo4d::{intercept, Ecef, Ecef4, Epoch, StateVector, Vec3};
/// let t = Epoch::from_unix_seconds(1_700_000_000.0);
/// let origin = Ecef4::new(Ecef::new(0.0, 0.0, 0.0), t);
/// // Target crossing at 100 m/s, 1 km ahead; interceptor at 300 m/s catches it.
/// let target = StateVector::new(t, Ecef::new(1_000.0, 0.0, 0.0), Vec3::new(0.0, 100.0, 0.0));
/// let ix = intercept(origin, target, 300.0).unwrap();
/// assert!(ix.epoch > t);
/// assert!((ix.velocity.norm() - 300.0).abs() < 1e-6);       // flies at the requested speed
/// // ...and actually meets the target there.
/// assert!(ix.point.distance_to(target.at(ix.epoch).position) < 1e-6);
/// ```
pub fn intercept(origin: Ecef4, target: StateVector, speed_mps: f64) -> Option<Intercept> {
    let t0 = origin.epoch;
    let tgt0 = target.at(t0).position;
    let d = tgt0.vec() - origin.position.vec();
    let vt = target.velocity;
    let tau = smallest_positive_root(vt.dot(vt) - speed_mps * speed_mps, 2.0 * d.dot(vt), d.dot(d))?;
    let point = Ecef::from(tgt0.vec() + vt * tau);
    Some(Intercept {
        epoch: t0 + Duration::from_seconds(tau),
        point,
        velocity: (point.vec() - origin.position.vec()) / tau,
        time_to_go: Duration::from_seconds(tau),
    })
}

/// The smallest strictly-positive real root of `a·τ² + b·τ + c`, or `None`.
fn smallest_positive_root(a: f64, b: f64, c: f64) -> Option<f64> {
    const EPS: f64 = 1e-12;
    if a.abs() < EPS {
        if b.abs() < EPS {
            return None;
        }
        let r = -c / b;
        return (r > EPS).then_some(r);
    }
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return None;
    }
    let sq = disc.sqrt();
    let (r1, r2) = ((-b - sq) / (2.0 * a), (-b + sq) / (2.0 * a));
    let (lo, hi) = (r1.min(r2), r1.max(r2));
    if lo > EPS {
        Some(lo)
    } else if hi > EPS {
        Some(hi)
    } else {
        None
    }
}

/// Whether two objects are on a collision course: their closest approach is within `miss_tol_m` and
/// still ahead (constant bearing, decreasing range). A cheap early-out before a full solve.
pub fn on_collision_course(a: StateVector, b: StateVector, miss_tol_m: f64) -> bool {
    let c = closest_approach(a, b);
    c.distance_m <= miss_tol_m && c.tca >= a.epoch
}

impl StateVector {
    /// The closest approach with another state vector (see [`closest_approach`]).
    #[inline]
    pub fn closest_approach(self, other: StateVector) -> Conjunction {
        closest_approach(self, other)
    }

    /// Whether the two objects breach a [`ProtectedZone`] within `window` (see [`conflict`]).
    #[inline]
    pub fn conflict(self, other: StateVector, zone: ProtectedZone, window: TimeWindow) -> Option<Conflict> {
        conflict(self, other, zone, window)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Geodetic4;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn closest_approach_in_absolute_time() {
        let t = Epoch::from_unix_seconds(1_700_000_000.0);
        let a = StateVector::new(t, Ecef::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        // b is sampled 10 s later, but still on the same physical line — the absolute answer agrees.
        let b_now = StateVector::new(t, Ecef::new(10.0, 5.0, 0.0), Vec3::new(-1.0, 0.0, 0.0));
        let b_later = b_now.at(t + Duration::from_seconds(10.0));
        let c1 = closest_approach(a, b_now);
        let c2 = closest_approach(a, b_later);
        assert!(close((c1.tca - t).as_seconds(), 5.0, 1e-9));
        assert!(close((c2.tca - t).as_seconds(), 5.0, 1e-6)); // same TCA despite different sample epoch
        assert!(close(c1.distance_m, 5.0, 1e-9) && close(c2.distance_m, 5.0, 1e-6));
        assert!(close(c1.relative_speed_mps, 2.0, 1e-9));
        // Windowed: force the answer before the true TCA → distance grows.
        let w = TimeWindow::between(t, t + Duration::from_seconds(2.0));
        let cw = closest_approach_within(a, b_now, w);
        assert_eq!(cw.tca, t + Duration::from_seconds(2.0));
        assert!(cw.distance_m > 5.0);
    }

    #[test]
    fn head_on_conflict_has_a_symmetric_window() {
        let t = Epoch::from_unix_seconds(1_700_000_000.0);
        let mid = Geodetic4::from_lat_lon(0.0, 0.0, 3_000.0, t).to_ecef4().position;
        let a = StateVector::new(t, mid, Vec3::new(0.0, 100.0, 0.0));
        let b = StateVector::new(
            t,
            Ecef::new(mid.x, mid.y + 20_000.0, mid.z),
            Vec3::new(0.0, -100.0, 0.0),
        );
        let c = a
            .conflict(b, ProtectedZone::new(9_260.0, 300.0), TimeWindow::unbounded())
            .unwrap();
        // They meet at the midpoint at t+100 s (20 km closing at 200 m/s); entered before, left after.
        assert!(close((c.tca - t).as_seconds(), 100.0, 0.5));
        assert!(c.enters < c.tca && c.tca < c.exits);
        assert!(c.duration().as_seconds() > 0.0);
        // Vertically separated by more than the zone → no conflict. At (lat 0, lon 0) the local "up"
        // is the ECEF x-axis, so a 1 km vertical gap is a +x offset.
        let high = StateVector::new(
            t,
            Ecef::new(mid.x + 1_000.0, mid.y + 20_000.0, mid.z),
            Vec3::new(0.0, -100.0, 0.0),
        );
        assert!(a
            .conflict(high, ProtectedZone::new(9_260.0, 300.0), TimeWindow::unbounded())
            .is_none());
    }

    #[test]
    fn intercept_solves_and_meets_the_target() {
        let t = Epoch::from_unix_seconds(1_700_000_000.0);
        let origin = crate::Ecef4::new(Ecef::new(0.0, 0.0, 0.0), t);
        let target = StateVector::new(t, Ecef::new(1_000.0, 0.0, 0.0), Vec3::new(0.0, 100.0, 0.0));
        let ix = intercept(origin, target, 300.0).unwrap();
        assert!(ix.epoch > t && close(ix.velocity.norm(), 300.0, 1e-6));
        assert!(ix.point.distance_to(target.at(ix.epoch).position) < 1e-6);
        // Too slow to catch a target opening at 200 m/s straight away.
        let fleeing = StateVector::new(t, Ecef::new(1_000.0, 0.0, 0.0), Vec3::new(200.0, 0.0, 0.0));
        assert!(intercept(origin, fleeing, 100.0).is_none());
    }

    #[test]
    fn collision_course_matches_zero_miss() {
        let t = Epoch::from_unix_seconds(1_700_000_000.0);
        let a = StateVector::new(t, Ecef::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        let head_on = StateVector::new(t, Ecef::new(100.0, 0.0, 0.0), Vec3::new(-1.0, 0.0, 0.0));
        assert!(on_collision_course(a, head_on, 1.0));
        let passing = StateVector::new(t, Ecef::new(100.0, 50.0, 0.0), Vec3::new(-1.0, 0.0, 0.0));
        assert!(!on_collision_course(a, passing, 1.0));
    }
}
