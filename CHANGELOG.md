# Changelog

All notable changes to `geo4d`. Plain semver; tags `vX.Y.Z` at the published commit.

## 0.1.0 — 2026-09-14

Initial release. 4D (spatio-temporal) geodesy on top of `geo3d`.

- Time: `Epoch` and `Duration` (nanosecond-exact, `Ord`, TAI-internal) with `TimeScale` (TAI/GPS/TT/UTC), the IERS leap-second table, Julian date / MJD, GPS week and time-of-week, Gregorian calendars, ITRF decimal years and the IAU 2000 Earth Rotation Angle. `TimeWindow` with containment, clamping and intersection.
- Timestamped coordinates: `Geodetic4`, `Ecef4` (with time-aware ECI conversion) and `StateVector` (position + velocity + epoch), including `from_fixes`, dead reckoning and a course-over-ground bridge.
- `Trajectory`: time-ordered path with `Linear`, `Hermite` and `CatmullRom` interpolation, an `Extrapolation` policy, resampling and state/geodetic queries.
- Datum: `PlateMotion` (ITRF2020 plate-motion model, 13 named plates) for cross-epoch propagation, and `Helmert14` (14-parameter time-dependent Helmert) in both rotation conventions, with named ITRF2020↔ITRF2014/2008 and GDA2020 transforms.
- Kinematics in absolute time: `closest_approach` (+ windowed), `conflict` against a cylindrical `ProtectedZone`, `intercept` (fixed-speed lead-collision) and `on_collision_course`.
- `collision_probability`: short-term-encounter probability of collision (Chan's series) with a `Covariance3` type.
- `Volume4` / `Timed<V>`: box any `geo3d` volume into a `TimeWindow` for 4D geofencing, with sampled trajectory intersection and a bisected entry epoch.
- Re-exports the whole `geo3d` API; optional `serde` feature (forwarded to `geo3d`).
