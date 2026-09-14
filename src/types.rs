//! Timestamped coordinates and state vectors — a 3D position bound to an [`Epoch`].
//!
//! A point in space is not a point in 4D space-time until it carries a time. [`Geodetic4`] and
//! [`Ecef4`] are a `geo3d` position plus an [`Epoch`]; [`StateVector`] adds a velocity, which is the
//! minimum needed to say where something *will* be. The 3D machinery — ellipsoids, ECEF ↔ geodetic,
//! local frames — is `geo3d`'s and is re-exported from the crate root.

use crate::time::{Duration, Epoch};
use geo3d::{eci, Ecef, Eci, Geodetic, LocalFrame, Track, Vec3, WGS84};

/// A geodetic position (latitude, longitude, height) at an instant. The "where and when" of a fix,
/// a waypoint or a report.
///
/// ```
/// use geo4d::{Epoch, Geodetic, Geodetic4};
/// let fix = Geodetic4::new(Geodetic::new(-33.87, 151.21, 20.0), Epoch::from_unix_seconds(1_700_000_000.0));
/// assert_eq!(fix.to_ecef4().to_geodetic4().position.lat_deg, fix.position.lat_deg);
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Geodetic4 {
    /// The geodetic position.
    pub position: Geodetic,
    /// When the position was valid.
    pub epoch: Epoch,
}

impl Geodetic4 {
    /// A geodetic position at an epoch.
    #[inline]
    pub const fn new(position: Geodetic, epoch: Epoch) -> Self {
        Geodetic4 { position, epoch }
    }

    /// From latitude and longitude in degrees, height in metres, and an epoch.
    #[inline]
    pub const fn from_lat_lon(lat_deg: f64, lon_deg: f64, height_m: f64, epoch: Epoch) -> Self {
        Geodetic4::new(Geodetic::new(lat_deg, lon_deg, height_m), epoch)
    }

    /// The same instant in ECEF on WGS84.
    #[inline]
    pub fn to_ecef4(self) -> Ecef4 {
        Ecef4::new(WGS84.to_ecef(self.position), self.epoch)
    }

    /// The same position at a different epoch (the position is unchanged; only the timestamp moves).
    #[inline]
    pub const fn at_epoch(self, epoch: Epoch) -> Self {
        Geodetic4 { epoch, ..self }
    }
}

/// An Earth-Centred Earth-Fixed position at an instant.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Ecef4 {
    /// The ECEF position, metres.
    pub position: Ecef,
    /// When the position was valid.
    pub epoch: Epoch,
}

impl Ecef4 {
    /// An ECEF position at an epoch.
    #[inline]
    pub const fn new(position: Ecef, epoch: Epoch) -> Self {
        Ecef4 { position, epoch }
    }

    /// The same instant as a geodetic position on WGS84.
    #[inline]
    pub fn to_geodetic4(self) -> Geodetic4 {
        Geodetic4::new(WGS84.to_geodetic(self.position), self.epoch)
    }

    /// The Earth-Centred Inertial position at this instant, rotating the fixed frame by the epoch's
    /// [Earth Rotation Angle](Epoch::earth_rotation_angle_rad). The time-aware counterpart of
    /// `geo3d`'s manual `ecef_to_eci`: here the sidereal angle is taken straight from the timestamp.
    #[inline]
    pub fn to_eci(self) -> Eci {
        eci::ecef_to_eci(self.position, self.epoch.earth_rotation_angle_rad())
    }

    /// The straight-line separation from another ECEF position, metres (the epoch is ignored).
    #[inline]
    pub fn distance_to(self, o: Ecef4) -> f64 {
        self.position.distance_to(o.position)
    }
}

/// A moving object at an instant: an ECEF position (metres) and velocity (metres per second),
/// stamped with the [`Epoch`] they are valid at. Constant velocity is assumed between updates, so a
/// state vector dead-reckons to any other epoch.
///
/// This is `geo3d`'s [`Track`] with the one thing a track lacks — an absolute time — so two of them
/// sampled at *different* epochs can still be compared (see [`closest_approach`](crate::closest_approach)).
///
/// ```
/// use geo4d::{Duration, Ecef, Epoch, StateVector, Vec3};
/// let t0 = Epoch::from_unix_seconds(1_700_000_000.0);
/// let s = StateVector::new(t0, Ecef::new(7_000_000.0, 0.0, 0.0), Vec3::new(0.0, 7_500.0, 0.0));
/// let later = s.at(t0 + Duration::from_seconds(10.0));           // dead reckon 10 s
/// assert!((later.position.y - 75_000.0).abs() < 1e-6);
/// assert_eq!(later.velocity, s.velocity);
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct StateVector {
    /// The epoch the position and velocity are valid at.
    pub epoch: Epoch,
    /// Position now, ECEF metres.
    pub position: Ecef,
    /// Velocity, ECEF m/s (constant between updates).
    pub velocity: Vec3,
}

impl StateVector {
    /// A state vector from an epoch, ECEF position and ECEF velocity.
    #[inline]
    pub const fn new(epoch: Epoch, position: Ecef, velocity: Vec3) -> Self {
        StateVector {
            epoch,
            position,
            velocity,
        }
    }

    /// A state vector inferred from two timestamped fixes (the later one becomes "now"). Coincident
    /// or out-of-order epochs give zero velocity.
    pub fn from_fixes(earlier: Ecef4, later: Ecef4) -> Self {
        let dt = (later.epoch - earlier.epoch).as_seconds();
        let velocity = if dt > 0.0 {
            (later.position.vec() - earlier.position.vec()) / dt
        } else {
            Vec3::ZERO
        };
        StateVector::new(later.epoch, later.position, velocity)
    }

    /// The state dead-reckoned to another epoch (past or future), assuming constant velocity. The
    /// velocity is carried over unchanged; the epoch and position advance.
    #[inline]
    pub fn at(self, epoch: Epoch) -> Self {
        let dt = (epoch - self.epoch).as_seconds();
        StateVector::new(epoch, self.position.offset(self.velocity * dt), self.velocity)
    }

    /// The position dead-reckoned to another epoch, as a timestamped ECEF point.
    #[inline]
    pub fn position_at(self, epoch: Epoch) -> Ecef4 {
        let s = self.at(epoch);
        Ecef4::new(s.position, s.epoch)
    }

    /// Speed through space, m/s.
    #[inline]
    pub fn speed_mps(self) -> f64 {
        self.velocity.norm()
    }

    /// The timestamped position (velocity dropped).
    #[inline]
    pub fn ecef4(self) -> Ecef4 {
        Ecef4::new(self.position, self.epoch)
    }

    /// As a `geo3d` [`Track`] (position + velocity), dropping the absolute epoch — for the
    /// relative-time kinematics in `geo3d`.
    #[inline]
    pub fn to_track(self) -> Track {
        Track::new(self.position, self.velocity)
    }

    /// The span from this state's epoch to another (`other − self`).
    #[inline]
    pub fn dt_to(self, other: StateVector) -> Duration {
        other.epoch - self.epoch
    }
}

impl From<StateVector> for Track {
    #[inline]
    fn from(s: StateVector) -> Track {
        s.to_track()
    }
}

impl Geodetic4 {
    /// A state vector at this position and time from a local ENU velocity given as a `geo3d`
    /// [`Course`](geo3d::Course) (course over ground, ground speed, climb) — the shape ADS-B/AIS/GNSS
    /// report. The course is resolved in the local frame at this position and turned into an ECEF
    /// velocity.
    pub fn state_from_course(self, course: geo3d::Course) -> StateVector {
        let frame = LocalFrame::new(self.position);
        StateVector::new(
            self.epoch,
            WGS84.to_ecef(self.position),
            frame.velocity_of_course(course),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::TimeScale;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn geodetic4_ecef4_round_trip() {
        let t = Epoch::from_gregorian_utc(2026, 1, 2, 3, 4, 5.0);
        let g = Geodetic4::from_lat_lon(-33.87, 151.21, 20.0, t);
        let back = g.to_ecef4().to_geodetic4();
        assert!(close(back.position.lat_deg, g.position.lat_deg, 1e-9));
        assert!(close(back.position.height_m, g.position.height_m, 1e-4));
        assert_eq!(back.epoch, t);
    }

    #[test]
    fn state_vector_dead_reckons_and_infers() {
        let t0 = Epoch::from_unix_seconds(1_700_000_000.0);
        let s = StateVector::new(t0, Ecef::new(0.0, 0.0, 0.0), Vec3::new(10.0, 0.0, 0.0));
        let s10 = s.at(t0 + Duration::from_seconds(10.0));
        assert_eq!(s10.position, Ecef::new(100.0, 0.0, 0.0));
        assert!(close(s10.speed_mps(), 10.0, 1e-12));
        // Inferring velocity from two fixes 5 s apart.
        let a = Ecef4::new(Ecef::new(0.0, 0.0, 0.0), t0);
        let b = Ecef4::new(Ecef::new(50.0, 0.0, 0.0), t0 + Duration::from_seconds(5.0));
        let inferred = StateVector::from_fixes(a, b);
        assert_eq!(inferred.velocity, Vec3::new(10.0, 0.0, 0.0));
        assert_eq!(inferred.epoch, b.epoch);
        // Out-of-order fixes → zero velocity.
        assert_eq!(StateVector::from_fixes(b, a).velocity, Vec3::ZERO);
    }

    #[test]
    fn eci_conversion_uses_the_epoch() {
        // At an epoch whose ERA ≈ 0 the ECI and ECEF frames nearly coincide; a quarter turn later
        // they don't. Here we just confirm the round trip and that the epoch drives the rotation.
        let t = Epoch::from_julian_date(TimeScale::Utc, 2_451_545.0);
        let p = Ecef4::new(Ecef::new(7_000_000.0, 0.0, 0.0), t);
        let eci = p.to_eci();
        let back = eci::eci_to_ecef(eci, t.earth_rotation_angle_rad());
        assert!(close(back.x, p.position.x, 1e-6) && close(back.y, p.position.y, 1e-6));
    }

    #[test]
    fn state_from_course_matches_speed() {
        let t = Epoch::now();
        let g = Geodetic4::from_lat_lon(-33.9, 151.2, 100.0, t);
        let s = g.state_from_course(geo3d::Course::new(90.0, 200.0, 0.0)); // due east, 200 m/s
        assert!(close(s.speed_mps(), 200.0, 1e-6));
        assert_eq!(s.epoch, t);
    }
}
