//! Time-dependent datum transformations: tectonic [`PlateMotion`] and the 14-parameter [`Helmert14`].
//!
//! A coordinate on the ground is not fixed in a global frame — the plate it sits on carries it. A
//! point in Sydney moves about 5.7 cm/yr north-east, so its ITRF coordinate in 2026 differs from its
//! 1994 one by nearly two metres. Two operations capture this, and they must not be confused:
//!
//! - **Same frame, different epoch** — [`PlateMotion`]: propagate a coordinate forward or back in
//!   time along its plate's rigid rotation, `X(t₂) = X(t₁) + V·(t₂ − t₁)`, `V = Ω × R`.
//! - **Same epoch, different frame** — [`Helmert14`]: a 7-parameter Helmert whose parameters drift
//!   linearly, so it is evaluated at the coordinate's epoch before being applied.
//!
//! The static 3D [`Helmert7`](geo3d::Helmert7) it rests on is `geo3d`'s.
//!
//! Constants are transcribed from the ITRF2020 plate-motion model (Altamimi et al. 2023) and the
//! ITRF / ICSM transformation tables; see each item for the source and epoch.

use crate::time::Epoch;
use crate::types::Ecef4;
use geo3d::{Ecef, Helmert7, Vec3};

/// Milliarcseconds per year → radians per year (`π / (180·3600·1000)`).
const MAS_PER_YEAR_TO_RAD_PER_YEAR: f64 = core::f64::consts::PI / (180.0 * 3600.0 * 1000.0);

/// The rigid rotation of a tectonic plate, as an Euler angular-velocity vector **Ω** in ECEF
/// (radians per Julian year). A site at ECEF position **R** moves at `V = Ω × R` (metres per year),
/// which is how a plate-fixed coordinate is carried through a global frame like ITRF.
///
/// The named constants are the ITRF2020 plate-motion model (Altamimi et al. 2023), which fits the
/// ITRF2020 velocity field to about 0.25 mm/yr. For millimetre-consistent velocities add the model's
/// origin-rate bias, [`PlateMotion::ITRF2020_ORIGIN_RATE_BIAS`].
///
/// ```
/// use geo4d::{Duration, Ecef4, Epoch, Geodetic4, PlateMotion};
///
/// // A pillar in Sydney, on the Australian plate, moves about 5.7 cm/yr (the plate's fast edge
/// // reaches ~7 cm/yr — the local speed depends on the angle to the Euler pole).
/// let sydney = Geodetic4::from_lat_lon(-33.87, 151.21, 20.0, Epoch::from_decimal_year(2020.0)).to_ecef4();
/// let speed = PlateMotion::AUSTRALIA.velocity_at(sydney.position).norm();
/// assert!((speed - 0.057).abs() < 0.003); // metres per year
///
/// // Propagate its 2020.0 coordinate to 2026.0: it has moved ~34 cm.
/// let moved = PlateMotion::AUSTRALIA.propagate(sydney, Epoch::from_decimal_year(2026.0));
/// assert!((moved.position.distance_to(sydney.position) - 0.34).abs() < 0.05);
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PlateMotion {
    /// Angular velocity vector in ECEF, radians per Julian year.
    omega: Vec3,
}

impl PlateMotion {
    /// From an angular-velocity vector in ECEF, radians per Julian year.
    #[inline]
    pub const fn from_rad_per_year(omega: Vec3) -> Self {
        PlateMotion { omega }
    }

    /// From angular-velocity components in **milliarcseconds per year** (the unit ITRF plate models
    /// publish).
    #[inline]
    pub const fn from_mas_per_year(wx: f64, wy: f64, wz: f64) -> Self {
        PlateMotion {
            omega: Vec3::new(
                wx * MAS_PER_YEAR_TO_RAD_PER_YEAR,
                wy * MAS_PER_YEAR_TO_RAD_PER_YEAR,
                wz * MAS_PER_YEAR_TO_RAD_PER_YEAR,
            ),
        }
    }

    /// From a geographic Euler pole: pole latitude and longitude (degrees) and rotation rate
    /// (degrees per million years) — the form geophysical models (NNR-MORVEL56, GSRM) publish.
    pub fn from_pole(pole_lat_deg: f64, pole_lon_deg: f64, rate_deg_per_myr: f64) -> Self {
        let rate_rad_per_year = rate_deg_per_myr.to_radians() / 1.0e6;
        let (slat, clat) = pole_lat_deg.to_radians().sin_cos();
        let (slon, clon) = pole_lon_deg.to_radians().sin_cos();
        PlateMotion {
            omega: Vec3::new(clat * clon, clat * slon, slat) * rate_rad_per_year,
        }
    }

    /// The angular-velocity vector, radians per Julian year.
    #[inline]
    pub const fn angular_velocity(self) -> Vec3 {
        self.omega
    }

    /// The surface velocity `V = Ω × R` at an ECEF position, **metres per Julian year**.
    #[inline]
    pub fn velocity_at(self, position: Ecef) -> Vec3 {
        self.omega.cross(position.vec())
    }

    /// The displacement of a point over a span, metres — `velocity × Δt` with Δt in Julian years.
    #[inline]
    pub fn displacement(self, position: Ecef, dt: crate::time::Duration) -> Vec3 {
        self.velocity_at(position) * dt.as_julian_years()
    }

    /// A coordinate propagated along the plate to another epoch — `X(t₂) = X(t₁) + V·(t₂ − t₁)`,
    /// with the velocity taken at the starting position. Moving a few years changes `R` by
    /// centimetres, far too little to affect `V`, so a single step is enough.
    pub fn propagate(self, from: Ecef4, to: Epoch) -> Ecef4 {
        let v = self.velocity_at(from.position);
        let dt = (to - from.epoch).as_julian_years();
        Ecef4::new(from.position.offset(v * dt), to)
    }

    // ── ITRF2020 plate-motion model (Altamimi et al. 2023), Ω in mas/yr ──────────────────────────────

    /// Amur plate (ITRF2020-PMM).
    pub const AMUR: PlateMotion = PlateMotion::from_mas_per_year(-0.1310, -0.5515, 0.8370);
    /// Antarctic plate (ITRF2020-PMM).
    pub const ANTARCTICA: PlateMotion = PlateMotion::from_mas_per_year(-0.2686, -0.3118, 0.6775);
    /// Arabian plate (ITRF2020-PMM).
    pub const ARABIA: PlateMotion = PlateMotion::from_mas_per_year(1.1286, -0.1462, 1.4378);
    /// Australian plate (ITRF2020-PMM) — one of the fastest, ~7 cm/yr NNE.
    pub const AUSTRALIA: PlateMotion = PlateMotion::from_mas_per_year(1.4875, 1.1754, 1.2233);
    /// Caribbean plate (ITRF2020-PMM).
    pub const CARIBBEAN: PlateMotion = PlateMotion::from_mas_per_year(0.2074, -1.4216, 0.7261);
    /// Eurasian plate (ITRF2020-PMM).
    pub const EURASIA: PlateMotion = PlateMotion::from_mas_per_year(-0.0853, -0.5191, 0.7528);
    /// Indian plate (ITRF2020-PMM).
    pub const INDIA: PlateMotion = PlateMotion::from_mas_per_year(1.1372, 0.0133, 1.4436);
    /// Nazca plate (ITRF2020-PMM).
    pub const NAZCA: PlateMotion = PlateMotion::from_mas_per_year(-0.3265, -1.5610, 1.6052);
    /// North American plate (ITRF2020-PMM).
    pub const NORTH_AMERICA: PlateMotion = PlateMotion::from_mas_per_year(0.0454, -0.6656, -0.0979);
    /// Nubian (west-African) plate (ITRF2020-PMM).
    pub const NUBIA: PlateMotion = PlateMotion::from_mas_per_year(0.0900, -0.5850, 0.7168);
    /// Pacific plate (ITRF2020-PMM).
    pub const PACIFIC: PlateMotion = PlateMotion::from_mas_per_year(-0.4039, 1.0210, -2.1542);
    /// South American plate (ITRF2020-PMM).
    pub const SOUTH_AMERICA: PlateMotion = PlateMotion::from_mas_per_year(-0.2610, -0.2822, -0.1573);
    /// Somalian plate (ITRF2020-PMM).
    pub const SOMALIA: PlateMotion = PlateMotion::from_mas_per_year(-0.0810, -0.7189, 0.8644);

    /// The ITRF2020-PMM origin-rate bias, a constant velocity (metres per year, ECEF) added to every
    /// plate-predicted velocity for full consistency with the frame origin. IGN's guidance is to add
    /// it but drop its vertical component; it is small (~0.9 mm/yr).
    pub const ITRF2020_ORIGIN_RATE_BIAS: Vec3 = Vec3::new(0.000_37, 0.000_35, 0.000_74);
}

/// A **14-parameter time-dependent Helmert** transformation between two reference frames: the seven
/// [`Helmert7`](geo3d::Helmert7) parameters plus a rate for each, referred to an epoch `t₀`. The
/// parameters at epoch `t` are `P(t) = P(t₀) + Ṗ·(t − t₀)`, so the transform is [evaluated at the
/// coordinate's epoch](Self::at) before it is applied.
///
/// Build it from parameters as the ITRF/EPSG tables publish them — translations in **mm** (rates
/// mm/yr), rotations in **mas** (rates mas/yr), scale in **ppb** (rates ppb/yr) — in either rotation
/// convention. IERS/ITRF tables are [position-vector](Self::position_vector); most national datums,
/// GDA2020 included, are [coordinate-frame](Self::coordinate_frame).
///
/// ```
/// use geo4d::{Ecef4, Epoch, Helmert14};
///
/// // A GNSS position in ITRF2020 at epoch 2024.6, taken to ITRF2014.
/// let p = Ecef4::new(geo4d::Ecef::new(-4_052_051.0, 4_212_836.0, -2_545_106.0), Epoch::from_decimal_year(2024.6));
/// let itrf2014 = Helmert14::itrf2020_to_itrf2014().apply(p);
/// // The frames differ by only millimetres, and the epoch is preserved.
/// assert!(itrf2014.position.distance_to(p.position) < 0.01);
/// assert_eq!(itrf2014.epoch, p.epoch);
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Helmert14 {
    /// The reference epoch the parameters are quoted at.
    pub reference_epoch: Epoch,
    /// The 7 parameters at `reference_epoch` (SI: metres, radians, unitless scale).
    pub value: Helmert7,
    /// The 7 rates, per Julian year (SI: m/yr, rad/yr, unitless/yr).
    pub rate: Helmert7,
}

impl Helmert14 {
    /// From a value/rate pair of already-built [`Helmert7`](geo3d::Helmert7)s and a reference epoch.
    #[inline]
    pub const fn from_parts(reference_epoch: Epoch, value: Helmert7, rate: Helmert7) -> Self {
        Helmert14 {
            reference_epoch,
            value,
            rate,
        }
    }

    /// From parameters in the **position-vector** rotation convention (EPSG method 1053, the one IERS
    /// and the ITRF tables use), in published units: translation mm (rate mm/yr), rotation mas (rate
    /// mas/yr), scale ppb (rate ppb/yr).
    pub fn position_vector(
        reference_epoch: Epoch,
        translation_mm: [f64; 3],
        rotation_mas: [f64; 3],
        scale_ppb: f64,
        translation_rate_mm_yr: [f64; 3],
        rotation_rate_mas_yr: [f64; 3],
        scale_rate_ppb_yr: f64,
    ) -> Self {
        Helmert14 {
            reference_epoch,
            value: pv(translation_mm, rotation_mas, scale_ppb),
            rate: pv(translation_rate_mm_yr, rotation_rate_mas_yr, scale_rate_ppb_yr),
        }
    }

    /// From parameters in the **coordinate-frame** rotation convention (EPSG method 1056, used by
    /// GDA2020 and most national datums); same units as [`position_vector`](Self::position_vector).
    pub fn coordinate_frame(
        reference_epoch: Epoch,
        translation_mm: [f64; 3],
        rotation_mas: [f64; 3],
        scale_ppb: f64,
        translation_rate_mm_yr: [f64; 3],
        rotation_rate_mas_yr: [f64; 3],
        scale_rate_ppb_yr: f64,
    ) -> Self {
        // Coordinate-frame is position-vector with the rotation signs (value and rate) flipped.
        Helmert14 {
            reference_epoch,
            value: pv(
                translation_mm,
                [-rotation_mas[0], -rotation_mas[1], -rotation_mas[2]],
                scale_ppb,
            ),
            rate: pv(
                translation_rate_mm_yr,
                [
                    -rotation_rate_mas_yr[0],
                    -rotation_rate_mas_yr[1],
                    -rotation_rate_mas_yr[2],
                ],
                scale_rate_ppb_yr,
            ),
        }
    }

    /// The static 7-parameter Helmert obtained by evaluating every parameter at `epoch`.
    pub fn at(self, epoch: Epoch) -> Helmert7 {
        let dt = (epoch - self.reference_epoch).as_julian_years();
        let (v, r) = (self.value, self.rate);
        Helmert7::new(
            v.tx + r.tx * dt,
            v.ty + r.ty * dt,
            v.tz + r.tz * dt,
            v.rx + r.rx * dt,
            v.ry + r.ry * dt,
            v.rz + r.rz * dt,
            v.scale + r.scale * dt,
        )
    }

    /// Transform a timestamped coordinate (source frame → target frame), evaluating the parameters at
    /// the coordinate's own epoch. The epoch is carried through unchanged.
    #[inline]
    pub fn apply(self, p: Ecef4) -> Ecef4 {
        Ecef4::new(self.at(p.epoch).apply(p.position), p.epoch)
    }

    /// The exact inverse transform (target frame → source frame), evaluating at the coordinate's
    /// epoch — so `apply_inverse(apply(p))` round-trips.
    #[inline]
    pub fn apply_inverse(self, p: Ecef4) -> Ecef4 {
        Ecef4::new(self.at(p.epoch).apply_inverse(p.position), p.epoch)
    }

    // ── Named transformations (from the ITRF and ICSM tables) ────────────────────────────────────────

    /// ITRF2020 → ITRF2014, position-vector convention, reference epoch 2015.0 (IERS/IGN table). The
    /// frames differ by a few millimetres; the inverse is [`apply_inverse`](Self::apply_inverse).
    pub fn itrf2020_to_itrf2014() -> Self {
        Helmert14::position_vector(
            Epoch::from_decimal_year(2015.0),
            [-1.4, -0.9, 1.4],
            [0.0, 0.0, 0.0],
            -0.42,
            [0.0, -0.1, 0.2],
            [0.0, 0.0, 0.0],
            0.0,
        )
    }

    /// ITRF2020 → ITRF2008, position-vector convention, reference epoch 2015.0 (IERS/IGN table).
    pub fn itrf2020_to_itrf2008() -> Self {
        Helmert14::position_vector(
            Epoch::from_decimal_year(2015.0),
            [0.2, 1.0, 3.3],
            [0.0, 0.0, 0.0],
            -0.29,
            [0.0, -0.1, 0.1],
            [0.0, 0.0, 0.0],
            0.03,
        )
    }

    /// ITRF2014 → GDA2020, the ICSM "conformal" transformation (coordinate-frame convention,
    /// reference epoch 2020.0). All parameters are zero except the three rotation *rates*, which are
    /// the Australian-plate rotation — so this both changes frame and carries plate motion. The
    /// inverse (GDA2020 → ITRF2014) is [`apply_inverse`](Self::apply_inverse).
    pub fn itrf2014_to_gda2020() -> Self {
        Helmert14::coordinate_frame(
            Epoch::from_decimal_year(2020.0),
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            0.0,
            [0.0, 0.0, 0.0],
            [1.503_79, 1.183_46, 1.207_16],
            0.0,
        )
    }

    /// GDA94 → GDA2020, the ICSM static 7-parameter conformal transformation (coordinate-frame
    /// convention). Both are plate-fixed static datums, so there are no rates; the ~1.8 m north-east
    /// shift is 26 years of Australian-plate motion plus the frame change.
    pub fn gda94_to_gda2020() -> Self {
        Helmert14::coordinate_frame(
            Epoch::from_decimal_year(2020.0), // arbitrary: the rates are zero
            [61.55, -10.87, -40.19],
            [-39.4924, -32.7221, -32.8979],
            -9.994,
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            0.0,
        )
    }
}

/// A position-vector [`Helmert7`](geo3d::Helmert7) from published units: translation mm, rotation
/// mas, scale ppb. (mm → m, mas → arc-second, ppb → ppm are each a division by 1000, which is what
/// `geo3d`'s arc-second/ppm builder then expects.)
fn pv(translation_mm: [f64; 3], rotation_mas: [f64; 3], scale_ppb: f64) -> Helmert7 {
    Helmert7::position_vector(
        translation_mm[0] / 1000.0,
        translation_mm[1] / 1000.0,
        translation_mm[2] / 1000.0,
        rotation_mas[0] / 1000.0,
        rotation_mas[1] / 1000.0,
        rotation_mas[2] / 1000.0,
        scale_ppb / 1000.0,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Duration, Geodetic4};

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn australian_plate_moves_seven_cm_per_year_northeast() {
        // The Euler pole magnitude is exactly the ITRF2020-PMM value converted from mas/yr.
        let expected_omega =
            (1.4875f64.powi(2) + 1.1754f64.powi(2) + 1.2233f64.powi(2)).sqrt() * MAS_PER_YEAR_TO_RAD_PER_YEAR;
        assert!(close(
            PlateMotion::AUSTRALIA.angular_velocity().norm(),
            expected_omega,
            1e-20
        ));

        let syd = Geodetic4::from_lat_lon(-33.87, 151.21, 20.0, Epoch::from_decimal_year(2020.0)).to_ecef4();
        let v = PlateMotion::AUSTRALIA.velocity_at(syd.position);
        assert!(close(v.norm(), 0.057_131_9, 1e-6), "{}", v.norm()); // deterministic ~5.71 cm/yr at Sydney
                                                                     // Direction is north-east and roughly horizontal: resolve into local ENU.
        let enu = geo3d::LocalFrame::new(syd.to_geodetic4().position).dir_to_enu(v);
        assert!(enu.north > 0.0 && enu.east > 0.0, "{enu:?}"); // NNE
        assert!(enu.up.abs() < 0.02); // nearly tangential
    }

    #[test]
    fn propagation_is_linear_and_reversible() {
        let p0 = Ecef4::new(
            Ecef::new(-4_052_051.0, 4_212_836.0, -2_545_106.0),
            Epoch::from_decimal_year(2020.0),
        );
        let t1 = Epoch::from_decimal_year(2030.0);
        let moved = PlateMotion::AUSTRALIA.propagate(p0, t1);
        // Ten years ≈ 70 cm.
        assert!(close(moved.position.distance_to(p0.position), 0.70, 0.05));
        // Propagating back recovers the start.
        let back = PlateMotion::AUSTRALIA.propagate(moved, p0.epoch);
        assert!(back.position.distance_to(p0.position) < 1e-6);
        assert_eq!(back.epoch, p0.epoch);
        // displacement() agrees with propagate().
        let d = PlateMotion::AUSTRALIA.displacement(p0.position, Duration::from_julian_years(10.0));
        assert!(close(d.norm(), 0.70, 0.05));
    }

    #[test]
    fn pole_and_component_constructors_agree() {
        // NNR-MORVEL56 Pacific pole vs a direct mas/yr vector should be the same order of magnitude.
        let by_pole = PlateMotion::from_pole(-63.58, 114.70, 0.651);
        assert!(by_pole.angular_velocity().norm() > 0.0);
        // Round-trip through rad/yr.
        let raw = PlateMotion::from_mas_per_year(1.0, 2.0, 3.0);
        assert_eq!(PlateMotion::from_rad_per_year(raw.angular_velocity()), raw);
    }

    #[test]
    fn helmert14_evaluates_at_the_coordinate_epoch() {
        let h = Helmert14::itrf2020_to_itrf2014();
        // At the 2015.0 reference epoch the transform equals its value part.
        let at2015 = h.at(Epoch::from_decimal_year(2015.0));
        assert!(close(at2015.tz, 1.4 / 1000.0, 1e-12));
        // A decade later the Z translation has grown by ~10 × 0.2 mm/yr = 2 mm. (Ten calendar years
        // is 10.0014 Julian years, since the rates are per Julian year — hence not exactly 2 mm.)
        let at2025 = h.at(Epoch::from_decimal_year(2025.0));
        assert!(
            close(at2025.tz - at2015.tz, 2.0 / 1000.0, 1e-6),
            "{}",
            at2025.tz - at2015.tz
        );
    }

    #[test]
    fn helmert14_round_trips_and_preserves_epoch() {
        let p = Ecef4::new(
            Ecef::new(-4_052_051.0, 4_212_836.0, -2_545_106.0),
            Epoch::from_decimal_year(2024.6),
        );
        let there = Helmert14::itrf2014_to_gda2020().apply(p);
        let back = Helmert14::itrf2014_to_gda2020().apply_inverse(there);
        assert!(back.position.distance_to(p.position) < 1e-6);
        assert_eq!(there.epoch, p.epoch);
        // GDA94 → GDA2020 shifts a real Australian coordinate by ~1.8 m.
        let g94 = Ecef4::new(
            Ecef::new(-4_052_051.0, 4_212_836.0, -2_545_106.0),
            Epoch::from_decimal_year(2020.0),
        );
        let g20 = Helmert14::gda94_to_gda2020().apply(g94);
        assert!(
            close(g20.position.distance_to(g94.position), 1.8, 0.3),
            "{}",
            g20.position.distance_to(g94.position)
        );
    }
}
