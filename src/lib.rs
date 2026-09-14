//! 4D (spatio-temporal) geodesy in pure `f64`, built on [`geo3d`].
//!
//! The fourth dimension is **time**. `geo3d` answers *where is this?*; `geo4d` answers *where is
//! this, when is it there, and how does the Earth itself move underneath it in the meantime?* It adds
//! the time axis to the 3D geometry and nothing else — the ellipsoids, ECEF/ENU frames, geodesics and
//! static volumes are `geo3d`'s, re-exported here so a program needs only one dependency.
//!
//! - [`Epoch`] / [`Duration`] — a high-precision atomic clock: nanosecond-exact, `Ord`, stored on
//!   [`TAI`](TimeScale::Tai) so arithmetic never trips over leap seconds. Converts to and from
//!   [`Gps`](TimeScale::Gps), [`Tt`](TimeScale::Tt), [`Utc`](TimeScale::Utc), Julian dates, GPS
//!   week/time-of-week, calendars and ITRF decimal years.
//! - [`Geodetic4`] / [`Ecef4`] / [`StateVector`] — a position, or a position and velocity, bound to
//!   an [`Epoch`]. A point is not a point in space-time until it carries a time.
//! - [`Trajectory`] — a time-ordered path with linear, cubic-[`Hermite`](Interpolation::Hermite) or
//!   [`CatmullRom`](Interpolation::CatmullRom) interpolation and an explicit [`Extrapolation`] policy.
//! - [`PlateMotion`] / [`Helmert14`] — the Earth moving: propagate a coordinate across epochs along
//!   its tectonic plate (ITRF2020 model), and transform between time-dependent reference frames
//!   (ITRF, GDA2020) with a 14-parameter Helmert.
//! - [`closest_approach`], [`conflict`], [`intercept`], [`collision_probability`] — 4D kinematics in
//!   absolute time: time of closest approach, cylindrical conflict windows, fixed-speed intercept,
//!   and the probability of collision (Chan's series on the encounter plane).
//! - [`Timed`] / [`Volume4`] — box any `geo3d` volume into a [`TimeWindow`] for 4D geofencing, and
//!   test whether a [`Trajectory`] enters it.
//!
//! Angles are **degrees**, lengths **metres**, times **seconds** (or a [`Duration`]) at the API.
//! Everything except [`Trajectory`] is `Copy` and allocation-free, and the whole crate is
//! `#![forbid(unsafe_code)]`.
//!
//! ```
//! use geo4d::{conflict, Duration, Epoch, Geodetic4, ProtectedZone, TimeWindow};
//!
//! // Two aircraft, each a position + course-over-ground at a shared time.
//! let now = Epoch::from_gregorian_utc(2026, 3, 14, 9, 30, 0.0);
//! let a = Geodetic4::from_lat_lon(-33.80, 151.30, 2_500.0, now)
//!     .state_from_course(geo4d::Course::new(225.0, 240.0, 0.0));   // SW at 240 m/s
//! let b = Geodetic4::from_lat_lon(-33.95, 151.10, 2_450.0, now)
//!     .state_from_course(geo4d::Course::new(45.0, 250.0, 0.0));    // NE at 250 m/s
//!
//! // When do they pass closest, and is 5 NM / 1000 ft ever breached in the next 10 minutes?
//! let cpa = a.closest_approach(b);
//! let zone = ProtectedZone::new(9_260.0, 300.0);
//! let window = TimeWindow::for_duration(now, Duration::from_minutes(10.0));
//! if let Some(c) = a.conflict(b, zone, window) {
//!     assert!(c.enters >= now && c.duration().as_seconds() > 0.0);
//! }
//! assert!(cpa.tca >= now);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

/// Every Rust code block in the README is compiled and run as a doctest.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
pub struct ReadmeDoctests;

mod conjunction;
mod datum;
mod kinematics;
mod time;
mod trajectory;
mod types;
mod volume;

pub use conjunction::{collision_probability, Covariance3};
pub use datum::{Helmert14, PlateMotion};
pub use kinematics::{
    closest_approach, closest_approach_within, conflict, intercept, on_collision_course, Conflict, Conjunction,
    Intercept, ProtectedZone,
};
pub use time::{Duration, Epoch, TimeScale, TimeWindow};
pub use trajectory::{Extrapolation, Interpolation, Trajectory, TrajectoryError, Waypoint};
pub use types::{Ecef4, Geodetic4, StateVector};
pub use volume::{Timed, Volume4};

// ── Re-exports of the geo3d 3D geometry geo4d builds on ──────────────────────────────────────────────
//
// So a program can `use geo4d::{Geodetic, WGS84, Sphere, ...}` and never name geo3d directly. The
// whole crate is also available as `geo4d::geo3d` for anything not surfaced here.
pub use geo3d;
pub use geo3d::Aer;
pub use geo3d::{
    eci, AltitudeBand, Cone, Course, Cpa, Cylinder, Ecef, Eci, Ellipsoid, Enu, Geodesic, Geodetic, Helmert7,
    LocalFrame, Ned, Sphere, Track, Vec3, Volume, AIRY_1830, BESSEL_1841, CLARKE_1866, GRS80, INTERNATIONAL_1924, PZ90,
    WGS72, WGS84,
};
