# geo4d

[![crates.io](https://img.shields.io/crates/v/geo4d.svg)](https://crates.io/crates/geo4d)
[![docs.rs](https://docs.rs/geo4d/badge.svg)](https://docs.rs/geo4d)
[![ci](https://github.com/FusionVault/geo4d/actions/workflows/ci.yml/badge.svg)](https://github.com/FusionVault/geo4d/actions/workflows/ci.yml)

4D (spatio-temporal) geodesy in pure `f64`, built on [`geo3d`](https://crates.io/crates/geo3d).

The fourth dimension is **time**. Where `geo3d` answers *where is this?*, `geo4d` answers *where is
this, when is it there, and how does the Earth itself move underneath it in the meantime?* It takes
the 3D geometry — ellipsoids, ECEF/ENU frames, geodesics, static volumes — and binds it to a
high-precision atomic clock, so you can dead-reckon a track to an exact instant, propagate a survey
mark across tectonic time, predict when two aircraft lose separation, or ask whether a drone enters a
temporary no-fly zone *while it is active*.

`geo3d` is a dependency and its whole API is re-exported, so `cargo add geo4d` gives you the 3D stack
(`Geodetic`, `Ecef`, `LocalFrame`, `WGS84`, `Sphere`, …) **and** the time axis, from one crate.
Angles are degrees, lengths metres, times seconds (or a `Duration`); everything except `Trajectory`
is `Copy`, allocation-free and `#![forbid(unsafe_code)]`.

## Install

```bash
cargo add geo4d
```

Optional: `cargo add geo4d --features serde` derives `Serialize`/`Deserialize` on every type (and
turns on `geo3d`'s serde too).

## Quick start

```rust
use geo4d::{Duration, Epoch, Geodetic4, ProtectedZone, TimeWindow};

let now = Epoch::from_gregorian_utc(2026, 3, 14, 9, 30, 0.0);

// Two aircraft, each a position + course-over-ground (course°, ground speed, climb) at that time.
let a = Geodetic4::from_lat_lon(-33.80, 151.30, 2_500.0, now)
    .state_from_course(geo4d::Course::new(225.0, 240.0, 0.0));
let b = Geodetic4::from_lat_lon(-33.95, 151.10, 2_450.0, now)
    .state_from_course(geo4d::Course::new(45.0, 250.0, 0.0));

// When is their closest approach, and how near?
let cpa = a.closest_approach(b);
assert!(cpa.tca >= now);

// Is a 5 NM / 1000 ft protected zone breached in the next 10 minutes?
let zone = ProtectedZone::new(9_260.0, 300.0);
let window = TimeWindow::for_duration(now, Duration::from_minutes(10.0));
let alert = a.conflict(b, zone, window);       // Some(Conflict { enters, exits, tca }) or None
assert!(alert.is_none() || alert.unwrap().enters >= now);
```

Every Rust example in this README is compiled and run as a doctest by CI, so they are safe to copy.

## Guide

### Time: epochs, durations and scales

An `Epoch` is one instant, stored to the nanosecond on the **TAI** (atomic) scale, so it is `Ord` and
epoch arithmetic never trips over leap seconds. A `Duration` is a signed span of whole nanoseconds.
`Epoch − Epoch` is a `Duration`; `Epoch ± Duration` is an `Epoch`.

Convert to and from the scales that matter — `Tai`, `Gps`, `Tt`, `Utc` — plus Julian dates, GPS
week/time-of-week, calendars and ITRF decimal years.

```rust
use geo4d::{Duration, Epoch, TimeScale};

// The J2000.0 epoch: 2000-01-01T12:00:00 TT = Julian Date 2451545.0.
let j2000 = Epoch::from_gregorian(TimeScale::Tt, 2000, 1, 1, 12, 0, 0.0);
assert_eq!(j2000, Epoch::J2000);
assert!((j2000.julian_date(TimeScale::Tt) - 2_451_545.0).abs() < 1e-9);
assert!((j2000.leap_seconds() - 32.0).abs() < 1e-9);   // TAI − UTC was 32 s in 2000

// Build from Unix time, read GPS week / time-of-week, do exact arithmetic.
let t = Epoch::from_unix_seconds(1_700_000_000.0);
let (week, tow) = t.gps_week_seconds();
assert_eq!(Epoch::from_gps_week_seconds(week, tow), t);
assert_eq!((t + Duration::from_minutes(5.0)) - t, Duration::from_minutes(5.0));

// Decimal years use the ITRF/GDA calendar convention: 2020.0 is 2020-01-01T00:00:00 UTC.
assert_eq!(Epoch::from_decimal_year(2020.0), Epoch::from_gregorian_utc(2020, 1, 1, 0, 0, 0.0));
```

`TimeWindow` is a span between two optional epochs — the temporal half of a geofence and the search
bound for conflict prediction.

```rust
use geo4d::{Duration, Epoch, TimeWindow};

let open = Epoch::from_gregorian_utc(2026, 7, 1, 19, 0, 0.0);
let match_time = TimeWindow::for_duration(open, Duration::from_hours(3.0));
assert!(match_time.contains(open + Duration::from_hours(1.0)));
assert!(!match_time.contains(open + Duration::from_hours(4.0)));
assert_eq!(match_time.duration(), Some(Duration::from_hours(3.0)));
```

### Timestamped coordinates and state vectors

`Geodetic4` and `Ecef4` are a `geo3d` position plus an `Epoch`. `StateVector` adds a velocity — the
minimum needed to say where something *will* be — and dead-reckons to any epoch.

```rust
use geo4d::{Duration, Ecef, Epoch, Geodetic4, StateVector, Vec3};

let t0 = Epoch::from_unix_seconds(1_700_000_000.0);

// A satellite state, propagated 10 s along its velocity.
let s = StateVector::new(t0, Ecef::new(7_000_000.0, 0.0, 0.0), Vec3::new(0.0, 7_500.0, 0.0));
let later = s.at(t0 + Duration::from_seconds(10.0));
assert!((later.position.y - 75_000.0).abs() < 1e-6);

// A geodetic fix and its ECEF / ECI counterparts at that instant.
let fix = Geodetic4::from_lat_lon(-33.87, 151.21, 20.0, t0);
let ecef = fix.to_ecef4();
let eci = ecef.to_eci();                       // rotated by the epoch's Earth-rotation angle
assert_eq!(ecef.to_geodetic4().position.lat_deg, fix.position.lat_deg);
assert!(eci.vec().norm() > 6_000_000.0);
```

### Trajectories

A `Trajectory` is a time-ordered path — a flight plan, an orbit arc, a track history — with a rule
for reading a position *between* the samples. The interpolation choices mirror the CCSDS Orbit
Ephemeris Message: piecewise `Linear` in ECEF, cubic `Hermite` when velocities are known, or
`CatmullRom` from positions alone. What happens off the ends is an explicit `Extrapolation` policy,
never a silent guess. This is the one part of the crate that allocates.

```rust
use geo4d::{Duration, Ecef, Epoch, Extrapolation, Interpolation, StateVector, Trajectory, Vec3};

let t0 = Epoch::from_unix_seconds(1_700_000_000.0);

// Sampled every minute, with velocities → cubic Hermite (C¹, passes through each state exactly).
let traj = Trajectory::from_states([
    StateVector::new(t0, Ecef::new(6_378_137.0, 0.0, 0.0), Vec3::new(0.0, 250.0, 5.0)),
    StateVector::new(t0 + Duration::from_minutes(1.0), Ecef::new(6_378_137.0, 15_000.0, 300.0), Vec3::new(-5.0, 250.0, 5.0)),
    StateVector::new(t0 + Duration::from_minutes(2.0), Ecef::new(6_377_900.0, 30_000.0, 600.0), Vec3::new(-10.0, 245.0, 5.0)),
]).unwrap().with_interpolation(Interpolation::Hermite);

let state = traj.state_at(t0 + Duration::from_seconds(90.0)).unwrap();  // position and velocity
assert!(state.position.y > 15_000.0 && state.position.y < 30_000.0);

// Resample onto a fixed grid, or read a geodetic fix at any instant.
let every_10s = traj.resample(Duration::from_seconds(10.0));
assert_eq!(every_10s.len(), 13);                         // 0 s .. 120 s inclusive
assert!(traj.geodetic_at(t0 + Duration::from_seconds(30.0)).is_some());

// Off the ends: None by default, or hold / extend if you opt in.
let held = traj.clone().with_extrapolation(Extrapolation::Clamp);
assert!(traj.position_at(t0 - Duration::from_seconds(1.0)).is_none());
assert!(held.position_at(t0 - Duration::from_seconds(1.0)).is_some());
```

### Datum shifts and plate motion

The ground moves. A mark in Sydney drifts about 5.7 cm/yr north-east, so its ITRF coordinate in 2026
is nearly two metres from its 1994 one. Two operations capture this, and they must not be confused:

- **Same frame, different epoch** — `PlateMotion`: propagate a coordinate along its tectonic plate,
  `X(t₂) = X(t₁) + V·(t₂ − t₁)`, using the ITRF2020 plate-motion model (13 named plates).
- **Same epoch, different frame** — `Helmert14`: a 14-parameter time-dependent Helmert, evaluated at
  the coordinate's epoch before it is applied. Named transforms for ITRF↔ITRF and Australia's
  GDA2020 are built in; build any other from the published parameters.

```rust
use geo4d::{Ecef, Ecef4, Epoch, Geodetic4, Helmert14, PlateMotion};

// Propagate a Sydney coordinate from 2020.0 to 2026.0 along the Australian plate.
let mark = Geodetic4::from_lat_lon(-33.87, 151.21, 20.0, Epoch::from_decimal_year(2020.0)).to_ecef4();
let in_2026 = PlateMotion::AUSTRALIA.propagate(mark, Epoch::from_decimal_year(2026.0));
assert!((in_2026.position.distance_to(mark.position) - 0.34).abs() < 0.05);   // ~34 cm in 6 yr

// Transform a GNSS position from ITRF2020 to ITRF2014 at its own epoch (a few mm apart).
let p = Ecef4::new(Ecef::new(-4_052_051.0, 4_212_836.0, -2_545_106.0), Epoch::from_decimal_year(2024.6));
let itrf2014 = Helmert14::itrf2020_to_itrf2014().apply(p);
assert!(itrf2014.position.distance_to(p.position) < 0.01);
assert_eq!(itrf2014.epoch, p.epoch);

// GDA94 → GDA2020 shifts a coordinate by 1.4–1.8 m across Australia (26 years of plate motion
// plus a frame change) — about 1.5 m at Sydney.
let g20 = Helmert14::gda94_to_gda2020().apply(mark);
assert!((g20.position.distance_to(mark.position) - 1.5).abs() < 0.2);
```

### Closest approach, conflict and intercept

4D kinematics in absolute time. Because each `StateVector` carries its own epoch, two tracks fixed at
*different* times can still be compared — the answers come back as epochs and windows.

```rust
use geo4d::{closest_approach, conflict, intercept, Ecef, Ecef4, Epoch, ProtectedZone, StateVector, TimeWindow, Vec3};

let t = Epoch::from_unix_seconds(1_700_000_000.0);

// Time and distance of closest approach.
let a = StateVector::new(t, Ecef::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
let b = StateVector::new(t, Ecef::new(10.0, 5.0, 0.0), Vec3::new(-1.0, 0.0, 0.0));
let cpa = closest_approach(a, b);
assert!((cpa.distance_m - 5.0).abs() < 1e-9);        // 5 m apart, 5 s from now

// Cylindrical loss-of-separation window (horizontal + vertical minima).
let conf = conflict(a, b, ProtectedZone::new(6.0, 100.0), TimeWindow::unbounded());
assert!(conf.is_some());

// Fixed-speed intercept of a moving target ("lead collision").
let origin = Ecef4::new(Ecef::new(0.0, 0.0, 0.0), t);
let target = StateVector::new(t, Ecef::new(1_000.0, 0.0, 0.0), Vec3::new(0.0, 100.0, 0.0));
let ix = intercept(origin, target, 300.0).unwrap();
assert!((ix.velocity.norm() - 300.0).abs() < 1e-6);  // flies at the requested speed
assert!(ix.point.distance_to(target.at(ix.epoch).position) < 1e-6);   // and meets it
```

### Probability of collision

A miss distance alone does not say how dangerous a conjunction is — that depends on how well each
object's position is known. Given a position covariance for each object and a combined hard-body
radius, `collision_probability` evaluates the short-term-encounter integral with Chan's series, the
method operational conjunction assessment uses.

```rust
use geo4d::{collision_probability, Covariance3, Ecef, Epoch, StateVector, Vec3};

let t = Epoch::from_unix_seconds(1_700_000_000.0);
let a = StateVector::new(t, Ecef::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
let b = StateVector::new(t, Ecef::new(100.0, 0.0, 0.0), Vec3::new(-1.0, 0.0, 0.0));

// Head-on, zero miss, 10 m 1σ each way, 5 m combined size → Pc = 1 − exp(−R²/2σ²) ≈ 6 %.
let pc = collision_probability(a, Covariance3::isotropic(10.0), b, Covariance3::isotropic(10.0), 5.0);
assert!((pc - 0.0606).abs() < 1e-3);
```

### 4D geofencing

A `Timed<V>` boxes any `geo3d` volume — `Sphere`, `Cone`, `Cylinder`, `AltitudeBand` — into a
`TimeWindow`, so a point is inside only when it is inside the shape **and** the window is open. A
`Trajectory` can be tested against one to find whether, and when, it first enters.

```rust
use geo4d::{AltitudeBand, Cylinder, Duration, Ecef, Epoch, Geodetic, Geodetic4, StateVector, Timed, Trajectory, Volume4};

// A temporary flight restriction over a stadium: 2 km radius, surface to 3 km, for a 3-hour match.
let kickoff = Epoch::from_gregorian_utc(2026, 7, 1, 19, 0, 0.0);
let tfr = Timed::for_duration(
    Cylinder::new(Geodetic::new(-33.85, 151.06, 0.0), 2_000.0, AltitudeBand::new(Some(0.0), Some(3_000.0))),
    kickoff,
    Duration::from_hours(3.0),
);

// Inside the cylinder, but only while the match is on.
let overhead = Geodetic::new(-33.85, 151.06, 1_500.0);
assert!(tfr.contains(Geodetic4::new(overhead, kickoff + Duration::from_hours(1.0))));
assert!(!tfr.contains(Geodetic4::new(overhead, kickoff + Duration::from_hours(5.0))));

// Does a drone's flight path enter it while it is active?
let start = Geodetic4::from_lat_lon(-33.85, 150.90, 1_500.0, kickoff).to_ecef4();
let vel = geo4d::LocalFrame::new(start.to_geodetic4().position).velocity_of_course(geo4d::Course::new(90.0, 60.0, 0.0));
let drone = Trajectory::from_states([StateVector::new(kickoff, start.position, vel).at(kickoff),
                                     StateVector::new(kickoff, start.position, vel).at(kickoff + Duration::from_hours(1.0))]).unwrap();
let _entry = tfr.first_entry(&drone, Duration::from_seconds(5.0));   // Some(epoch) or None
```

### Serde

With the `serde` feature every type — `Epoch`, `Duration`, `StateVector`, `Trajectory`, the volumes,
and the re-exported `geo3d` coordinate types — serialises to plain fields.

```rust,ignore
// cargo add geo4d --features serde
let json = serde_json::to_string(&geo4d::Epoch::from_unix_seconds(1_700_000_000.0))?;
```

### Interoperability

The re-exported coordinate types (`Ecef`, `Vec3`, `Enu`, …) bridge to the common math and GIS crates
through `geo3d`'s optional conversions, forwarded here as features: `mint`, `glam`, `nalgebra` and
`geo-types`. All default-off, so the default build pulls nothing extra.

```bash
cargo add geo4d --features glam
# then: let v: glam::DVec3 = ecef4.position.into();
```

## Precision and accuracy classes

Each function states its accuracy class, and the tests hold it:

- **Exact / integer** — `Epoch` and `Duration` are nanosecond-exact integers on the TAI timeline; the
  scale offsets (TT − TAI = 32.184 s, GPS = TAI − 19 s) and the IERS leap-second table (10 s in 1972
  to 37 s since 2017) are exact. Calendar and GPS-week round trips hold to the nanosecond.
- **First order (constant velocity)** — `StateVector` dead reckoning, `closest_approach`, `conflict`,
  `intercept` and `PlateMotion` propagation are linear-motion models, documented as such.
- **Model-bounded** — `collision_probability` is the short-term-encounter (rectilinear, Gaussian)
  model with Chan's series; `Helmert14` is the linearised 14-parameter Helmert; `Trajectory` linear
  interpolation sags below a curved path by about `d²/8R` (≈ 2 cm over 1 km).
- **Reference data** — plate poles are the ITRF2020 plate-motion model (Altamimi et al. 2023); the
  ITRF and GDA2020 transformation parameters are the published IERS / ICSM values.

The sidereal angle used for ECI is the IAU 2000 Earth Rotation Angle (with UT1 ≈ UTC); sub-arcsecond
precession, nutation and polar motion are out of scope, as in `geo3d`.

## Design

- Builds on `geo3d`; re-exports its whole API so downstream code needs one dependency. Cargo's
  dead-code elimination means you pay only for what you use.
- Time lives on a single continuous TAI counter; leap seconds, scale offsets and calendars are
  applied only at the boundary, so arithmetic can never land "inside" a leap second.
- Everything except `Trajectory` is `Copy` and allocation-free; `#![forbid(unsafe_code)]`.
- Degrees, metres and seconds at the API; radians and internal nanoseconds never leak out.

Not in scope, by design: an orbit propagator (use [`sgp4`](https://crates.io/crates/sgp4) and feed
its states to a `Trajectory`), full IERS Earth-orientation (UT1/polar motion) modelling, map
projections, and the 2D/planar and terrain geometry that lives in the wider
[`geo`](https://crates.io/crates/geo) ecosystem.

## Develop

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-features       # unit and doc tests, including this README
cargo run --example tracking    # a 4D conflict / geofence walkthrough
cargo publish --dry-run
```

CI runs all of the above plus a build on the declared minimum Rust version (1.85).

## Versioning and license

Plain semver; tags `vX.Y.Z` at the published commit; see [CHANGELOG.md](CHANGELOG.md).
Apache-2.0, see [LICENSE](LICENSE).
