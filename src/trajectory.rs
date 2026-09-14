//! Trajectories: a time-ordered path sampled at [`Waypoint`]s, interpolated to any [`Epoch`].
//!
//! Where a [`StateVector`] is one moving object frozen at an instant, a [`Trajectory`] is the whole
//! flight plan, orbit arc or track history — a sequence of timestamped positions (optionally with
//! velocities) and a rule for reading a position *between* them. The interpolation choices mirror the
//! CCSDS Orbit Ephemeris Message: piecewise [`Linear`](Interpolation::Linear) in ECEF, cubic
//! [`Hermite`](Interpolation::Hermite) when velocities are known, or [`CatmullRom`](Interpolation::CatmullRom)
//! from positions alone. What happens off the ends is an explicit [`Extrapolation`] policy, never a
//! silent guess.
//!
//! This is the one part of `geo4d` that allocates (it owns a `Vec` of waypoints); everything else is
//! `Copy` and allocation-free.

use crate::time::{Duration, Epoch, TimeWindow};
use crate::types::{Ecef4, Geodetic4, StateVector};
use geo3d::{Ecef, Vec3};

/// One sample on a [`Trajectory`]: a position at an epoch, optionally with a known velocity (which
/// [`Hermite`](Interpolation::Hermite) interpolation uses as the tangent).
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Waypoint {
    /// When the sample is valid.
    pub epoch: Epoch,
    /// The position, ECEF metres.
    pub position: Ecef,
    /// The velocity if known, ECEF m/s.
    pub velocity: Option<Vec3>,
}

impl Waypoint {
    /// A waypoint with a known velocity.
    #[inline]
    pub const fn new(epoch: Epoch, position: Ecef, velocity: Vec3) -> Self {
        Waypoint {
            epoch,
            position,
            velocity: Some(velocity),
        }
    }

    /// A waypoint that is position-only (no velocity).
    #[inline]
    pub const fn position_only(epoch: Epoch, position: Ecef) -> Self {
        Waypoint {
            epoch,
            position,
            velocity: None,
        }
    }
}

impl From<StateVector> for Waypoint {
    #[inline]
    fn from(s: StateVector) -> Waypoint {
        Waypoint::new(s.epoch, s.position, s.velocity)
    }
}
impl From<Ecef4> for Waypoint {
    #[inline]
    fn from(p: Ecef4) -> Waypoint {
        Waypoint::position_only(p.epoch, p.position)
    }
}

/// How a position is read between samples.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Interpolation {
    /// Straight line in ECEF (the chord). Exact for uniform motion; below a curved path it sags by
    /// about `d²/8R` (≈ 2 cm over 1 km, ≈ 2 m over 10 km). The default.
    #[default]
    Linear,
    /// Cubic Hermite using each waypoint's velocity as the tangent (a finite-difference tangent
    /// where a velocity is missing). C¹-continuous; the right choice for ephemeris state vectors.
    Hermite,
    /// Cubic Catmull-Rom from positions alone, tangents estimated from the neighbouring samples in
    /// the sample-time parameterization. Smooth pass-through of the points without needing velocities.
    CatmullRom,
}

/// What a query outside the sampled span returns.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Extrapolation {
    /// Nothing — a query before the first or after the last sample returns `None`. The safe default.
    #[default]
    None,
    /// Hold the nearest endpoint position.
    Clamp,
    /// Continue in a straight line at the endpoint's velocity (or the end-segment secant).
    Linear,
}

/// Why a trajectory could not be built or extended.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TrajectoryError {
    /// A waypoint's epoch was not strictly later than the previous one (times must be increasing and
    /// distinct).
    NonMonotonic,
}

impl core::fmt::Display for TrajectoryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            TrajectoryError::NonMonotonic => f.write_str("waypoint epochs must be strictly increasing"),
        }
    }
}
impl std::error::Error for TrajectoryError {}

/// A time-ordered path with an interpolation and an extrapolation policy. Waypoints are kept sorted
/// and strictly increasing in time.
///
/// ```
/// use geo4d::{Duration, Ecef, Epoch, Interpolation, StateVector, Trajectory, Vec3};
///
/// let t0 = Epoch::from_unix_seconds(1_700_000_000.0);
/// // A climbing turn sampled every 60 s, with velocities → cubic Hermite.
/// let traj = Trajectory::from_states([
///     StateVector::new(t0, Ecef::new(6_378_137.0, 0.0, 0.0), Vec3::new(0.0, 250.0, 5.0)),
///     StateVector::new(t0 + Duration::from_minutes(1.0), Ecef::new(6_378_137.0, 15_000.0, 300.0), Vec3::new(-5.0, 250.0, 5.0)),
///     StateVector::new(t0 + Duration::from_minutes(2.0), Ecef::new(6_377_900.0, 30_000.0, 600.0), Vec3::new(-10.0, 245.0, 5.0)),
/// ]).unwrap().with_interpolation(Interpolation::Hermite);
///
/// // Read the state at any instant inside the span.
/// let mid = traj.state_at(t0 + Duration::from_seconds(90.0)).unwrap();
/// assert!(mid.position.y > 15_000.0 && mid.position.y < 30_000.0);
/// assert_eq!(traj.window().start, Some(t0));
/// ```
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Trajectory {
    points: Vec<Waypoint>,
    interpolation: Interpolation,
    extrapolation: Extrapolation,
}

impl Trajectory {
    /// An empty trajectory (linear interpolation, no extrapolation). Add points with
    /// [`push`](Self::push).
    #[inline]
    pub fn new() -> Self {
        Trajectory {
            points: Vec::new(),
            interpolation: Interpolation::Linear,
            extrapolation: Extrapolation::None,
        }
    }

    /// A trajectory from state vectors (position + velocity → Hermite-ready). The samples are sorted
    /// by epoch; a duplicate epoch is a [`TrajectoryError::NonMonotonic`].
    pub fn from_states<I: IntoIterator<Item = StateVector>>(states: I) -> Result<Self, TrajectoryError> {
        Self::from_waypoints(states.into_iter().map(Waypoint::from))
    }

    /// A trajectory from timestamped positions (no velocities → linear or Catmull-Rom).
    pub fn from_samples<I: IntoIterator<Item = Ecef4>>(samples: I) -> Result<Self, TrajectoryError> {
        Self::from_waypoints(samples.into_iter().map(Waypoint::from))
    }

    /// A trajectory from waypoints (sorted by epoch; duplicate epochs are rejected).
    pub fn from_waypoints<I: IntoIterator<Item = Waypoint>>(waypoints: I) -> Result<Self, TrajectoryError> {
        let mut points: Vec<Waypoint> = waypoints.into_iter().collect();
        points.sort_by_key(|w| w.epoch);
        if points.windows(2).any(|w| w[0].epoch == w[1].epoch) {
            return Err(TrajectoryError::NonMonotonic);
        }
        Ok(Trajectory {
            points,
            interpolation: Interpolation::Linear,
            extrapolation: Extrapolation::None,
        })
    }

    /// Set the interpolation policy.
    #[inline]
    pub fn with_interpolation(mut self, interpolation: Interpolation) -> Self {
        self.interpolation = interpolation;
        self
    }

    /// Set the extrapolation policy.
    #[inline]
    pub fn with_extrapolation(mut self, extrapolation: Extrapolation) -> Self {
        self.extrapolation = extrapolation;
        self
    }

    /// Append a waypoint. It must be strictly later than the current last, else
    /// [`TrajectoryError::NonMonotonic`].
    pub fn push(&mut self, waypoint: Waypoint) -> Result<(), TrajectoryError> {
        if self.points.last().is_some_and(|last| waypoint.epoch <= last.epoch) {
            return Err(TrajectoryError::NonMonotonic);
        }
        self.points.push(waypoint);
        Ok(())
    }

    /// The waypoints, in time order.
    #[inline]
    pub fn waypoints(&self) -> &[Waypoint] {
        &self.points
    }

    /// The number of waypoints.
    #[inline]
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Whether the trajectory has no waypoints.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// The epoch of the first waypoint.
    #[inline]
    pub fn start(&self) -> Option<Epoch> {
        self.points.first().map(|w| w.epoch)
    }

    /// The epoch of the last waypoint.
    #[inline]
    pub fn end(&self) -> Option<Epoch> {
        self.points.last().map(|w| w.epoch)
    }

    /// The time span covered by the samples.
    #[inline]
    pub fn window(&self) -> TimeWindow {
        TimeWindow::new(self.start(), self.end())
    }

    /// The total sampled duration, or `None` if there are fewer than two waypoints.
    #[inline]
    pub fn duration(&self) -> Option<Duration> {
        match (self.start(), self.end()) {
            (Some(s), Some(e)) if s != e => Some(e - s),
            _ => None,
        }
    }

    /// The interpolated ECEF position at an epoch, honouring the extrapolation policy off the ends.
    pub fn position_at(&self, epoch: Epoch) -> Option<Ecef> {
        self.evaluate(epoch).map(|(p, _)| p)
    }

    /// The interpolated position as a timestamped ECEF point.
    #[inline]
    pub fn sample(&self, epoch: Epoch) -> Option<Ecef4> {
        self.position_at(epoch).map(|p| Ecef4::new(p, epoch))
    }

    /// The interpolated state (position **and** velocity) at an epoch. The velocity is the analytic
    /// derivative of the interpolant.
    pub fn state_at(&self, epoch: Epoch) -> Option<StateVector> {
        self.evaluate(epoch).map(|(p, v)| StateVector::new(epoch, p, v))
    }

    /// The interpolated position as a geodetic fix.
    #[inline]
    pub fn geodetic_at(&self, epoch: Epoch) -> Option<Geodetic4> {
        self.sample(epoch).map(Ecef4::to_geodetic4)
    }

    /// Positions resampled at a fixed `step` from the first waypoint to the last (inclusive).
    pub fn resample(&self, step: Duration) -> Vec<Ecef4> {
        let (Some(start), Some(end)) = (self.start(), self.end()) else {
            return Vec::new();
        };
        let mut out = Vec::new();
        let mut t = start;
        while t <= end {
            if let Some(p) = self.sample(t) {
                out.push(p);
            }
            t = t + step;
        }
        out
    }

    // ── interpolation internals ──────────────────────────────────────────────────────────────────────

    /// Position and velocity at `epoch`, or `None` per the extrapolation policy.
    fn evaluate(&self, epoch: Epoch) -> Option<(Ecef, Vec3)> {
        let n = self.points.len();
        if n == 0 {
            return None;
        }
        let (start, end) = (self.points[0].epoch, self.points[n - 1].epoch);
        if epoch < start {
            return self.off_end(0, epoch);
        }
        if epoch > end {
            return self.off_end(n - 1, epoch);
        }
        if n == 1 {
            // A single-point trajectory, queried exactly at its epoch.
            return Some((self.points[0].position, self.tangent(0)));
        }
        // Inside: find the segment [i, i+1] with points[i].epoch <= epoch <= points[i+1].epoch.
        let hi = self.points.partition_point(|w| w.epoch <= epoch);
        let i = hi.saturating_sub(1).min(n - 2);
        Some(self.eval_segment(i, epoch))
    }

    /// Handle a query beyond an endpoint index per the extrapolation policy.
    fn off_end(&self, idx: usize, epoch: Epoch) -> Option<(Ecef, Vec3)> {
        let w = self.points[idx];
        match self.extrapolation {
            Extrapolation::None => None,
            Extrapolation::Clamp => Some((w.position, Vec3::ZERO)),
            Extrapolation::Linear => {
                let v = self.tangent(idx);
                let dt = (epoch - w.epoch).as_seconds();
                Some((w.position.offset(v * dt), v))
            }
        }
    }

    /// Evaluate segment `i` (`0 ≤ i < n-1`) at `epoch`.
    fn eval_segment(&self, i: usize, epoch: Epoch) -> (Ecef, Vec3) {
        let (w0, w1) = (self.points[i], self.points[i + 1]);
        let dt = (w1.epoch - w0.epoch).as_seconds();
        let u = (epoch - w0.epoch).as_seconds() / dt;
        match self.interpolation {
            Interpolation::Linear => {
                let secant = (w1.position.vec() - w0.position.vec()) / dt;
                (Ecef::from(w0.position.vec() + secant * (u * dt)), secant)
            }
            Interpolation::Hermite | Interpolation::CatmullRom => {
                let (p0, p1) = (w0.position.vec(), w1.position.vec());
                let (m0, m1) = (self.tangent(i), self.tangent(i + 1));
                hermite(p0, m0, p1, m1, dt, u)
            }
        }
    }

    /// The tangent (dp/dt, m/s) at waypoint `i` for the active interpolation.
    fn tangent(&self, i: usize) -> Vec3 {
        match self.interpolation {
            Interpolation::Linear => self.secant_tangent(i),
            Interpolation::Hermite => self.points[i].velocity.unwrap_or_else(|| self.secant_tangent(i)),
            Interpolation::CatmullRom => self.catmull_tangent(i),
        }
    }

    /// A one- or two-sided secant estimate of the tangent at `i` (m/s).
    fn secant_tangent(&self, i: usize) -> Vec3 {
        let n = self.points.len();
        if n < 2 {
            return Vec3::ZERO;
        }
        let (a, b) = (i.saturating_sub(1), (i + 1).min(n - 1));
        let (pa, pb) = (self.points[a].position.vec(), self.points[b].position.vec());
        let span = (self.points[b].epoch - self.points[a].epoch).as_seconds();
        if span == 0.0 {
            Vec3::ZERO
        } else {
            (pb - pa) / span
        }
    }

    /// The non-uniform Catmull-Rom tangent at `i` (m/s), using the sample times as knots.
    fn catmull_tangent(&self, i: usize) -> Vec3 {
        let n = self.points.len();
        if n < 2 {
            return Vec3::ZERO;
        }
        if i == 0 || i == n - 1 {
            return self.secant_tangent(i);
        }
        let (pm, p0, pp) = (
            self.points[i - 1].position.vec(),
            self.points[i].position.vec(),
            self.points[i + 1].position.vec(),
        );
        let (tm, t0, tp) = (self.points[i - 1].epoch, self.points[i].epoch, self.points[i + 1].epoch);
        let d_prev = (t0 - tm).as_seconds();
        let d_next = (tp - t0).as_seconds();
        let d_span = (tp - tm).as_seconds();
        (p0 - pm) / d_prev - (pp - pm) / d_span + (pp - p0) / d_next
    }
}

impl Default for Trajectory {
    #[inline]
    fn default() -> Self {
        Trajectory::new()
    }
}

/// Cubic Hermite on a segment: positions `p0,p1`, tangents `m0,m1` (dp/dt), segment length `dt`
/// seconds, normalized parameter `u ∈ [0,1]`. Returns position and its time-derivative (velocity).
fn hermite(p0: Vec3, m0: Vec3, p1: Vec3, m1: Vec3, dt: f64, u: f64) -> (Ecef, Vec3) {
    let u2 = u * u;
    let u3 = u2 * u;
    let (h00, h10, h01, h11) = (
        2.0 * u3 - 3.0 * u2 + 1.0,
        u3 - 2.0 * u2 + u,
        -2.0 * u3 + 3.0 * u2,
        u3 - u2,
    );
    let pos = p0 * h00 + m0 * (h10 * dt) + p1 * h01 + m1 * (h11 * dt);
    let (d00, d10, d01, d11) = (
        6.0 * u2 - 6.0 * u,
        3.0 * u2 - 4.0 * u + 1.0,
        -6.0 * u2 + 6.0 * u,
        3.0 * u2 - 2.0 * u,
    );
    let vel = (p0 * d00 + m0 * (d10 * dt) + p1 * d01 + m1 * (d11 * dt)) / dt;
    (Ecef::from(pos), vel)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    fn t0() -> Epoch {
        Epoch::from_unix_seconds(1_700_000_000.0)
    }

    #[test]
    fn linear_interpolation_is_the_chord() {
        let a = Ecef4::new(Ecef::new(0.0, 0.0, 0.0), t0());
        let b = Ecef4::new(Ecef::new(100.0, 0.0, 0.0), t0() + Duration::from_seconds(10.0));
        let traj = Trajectory::from_samples([a, b]).unwrap();
        let mid = traj.position_at(t0() + Duration::from_seconds(5.0)).unwrap();
        assert_eq!(mid, Ecef::new(50.0, 0.0, 0.0));
        // Endpoints exact; velocity is the constant secant.
        assert_eq!(traj.position_at(t0()).unwrap(), a.position);
        assert_eq!(traj.position_at(b.epoch).unwrap(), b.position);
        assert!(close(
            traj.state_at(t0() + Duration::from_seconds(2.0)).unwrap().velocity.x,
            10.0,
            1e-9
        ));
    }

    #[test]
    fn hermite_passes_through_nodes_with_given_velocity() {
        let states = [
            StateVector::new(t0(), Ecef::new(0.0, 0.0, 0.0), Vec3::new(10.0, 0.0, 0.0)),
            StateVector::new(
                t0() + Duration::from_seconds(10.0),
                Ecef::new(100.0, 50.0, 0.0),
                Vec3::new(10.0, 0.0, 0.0),
            ),
        ];
        let traj = Trajectory::from_states(states)
            .unwrap()
            .with_interpolation(Interpolation::Hermite);
        // Nodes are hit exactly, with their stated velocity.
        for s in states {
            let got = traj.state_at(s.epoch).unwrap();
            assert!(got.position.distance_to(s.position) < 1e-6);
            assert!((got.velocity - s.velocity).norm() < 1e-6);
        }
        // Between nodes it curves (y rises then the tangent pulls it back), staying inside the box.
        let mid = traj.position_at(t0() + Duration::from_seconds(5.0)).unwrap();
        assert!(mid.x > 0.0 && mid.x < 100.0 && mid.y > 0.0 && mid.y <= 50.0);
    }

    #[test]
    fn extrapolation_policies() {
        let a = Ecef4::new(Ecef::new(0.0, 0.0, 0.0), t0());
        let b = Ecef4::new(Ecef::new(100.0, 0.0, 0.0), t0() + Duration::from_seconds(10.0));
        let before = t0() - Duration::from_seconds(5.0);
        let after = b.epoch + Duration::from_seconds(5.0);
        // None → out of range is None.
        let none = Trajectory::from_samples([a, b]).unwrap();
        assert!(none.position_at(before).is_none() && none.position_at(after).is_none());
        // Clamp → holds endpoints.
        let clamp = none.clone().with_extrapolation(Extrapolation::Clamp);
        assert_eq!(clamp.position_at(before).unwrap(), a.position);
        assert_eq!(clamp.position_at(after).unwrap(), b.position);
        // Linear → continues at the end secant (10 m/s): 5 s past the end is x = 150.
        let lin = none.with_extrapolation(Extrapolation::Linear);
        assert!(close(lin.position_at(after).unwrap().x, 150.0, 1e-6));
        assert!(close(lin.position_at(before).unwrap().x, -50.0, 1e-6));
    }

    #[test]
    fn catmull_rom_is_smooth_through_interior_points() {
        let pts: Vec<Ecef4> = (0..5)
            .map(|k| {
                Ecef4::new(
                    Ecef::new(k as f64 * 100.0, (k as f64 * 100.0).sin() * 50.0, 0.0),
                    t0() + Duration::from_seconds(k as f64 * 10.0),
                )
            })
            .collect();
        let traj = Trajectory::from_samples(pts.clone())
            .unwrap()
            .with_interpolation(Interpolation::CatmullRom);
        // Interpolant passes through every sample.
        for p in &pts {
            let got = traj.position_at(p.epoch).unwrap();
            assert!(got.distance_to(p.position) < 1e-6, "{got:?} vs {:?}", p.position);
        }
    }

    #[test]
    fn build_validation_and_resample() {
        // Duplicate epochs are rejected.
        let a = Ecef4::new(Ecef::new(0.0, 0.0, 0.0), t0());
        assert_eq!(
            Trajectory::from_samples([a, a]).err(),
            Some(TrajectoryError::NonMonotonic)
        );
        // Out-of-order push is rejected; in-order is accepted; unsorted input is sorted.
        let mut tr = Trajectory::new();
        assert!(tr.push(Waypoint::from(a)).is_ok());
        assert_eq!(
            tr.push(Waypoint::position_only(
                t0() - Duration::SECOND,
                Ecef::new(1.0, 0.0, 0.0)
            ))
            .err(),
            Some(TrajectoryError::NonMonotonic)
        );
        let b = Ecef4::new(Ecef::new(100.0, 0.0, 0.0), t0() + Duration::from_seconds(10.0));
        let traj = Trajectory::from_samples([b, a]).unwrap(); // reversed input
        assert_eq!(traj.start(), Some(t0()));
        assert_eq!(traj.duration(), Some(Duration::from_seconds(10.0)));
        // Resample every 2 s → 6 points (t0 .. t0+10 inclusive).
        let s = traj.resample(Duration::from_seconds(2.0));
        assert_eq!(s.len(), 6);
        assert_eq!(s[0].position, a.position);
        assert!(close(s[3].position.x, 60.0, 1e-9));
    }
}
