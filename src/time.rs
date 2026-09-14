//! High-precision time: [`Duration`], [`Epoch`] and the [`TimeScale`]s.
//!
//! The fourth dimension. A [`Duration`] is a signed span stored as whole nanoseconds (`i128`), so it
//! is exact, `Ord` and never drifts. An [`Epoch`] is one instant, stored internally on the **TAI**
//! (International Atomic Time) scale as nanoseconds from 1970-01-01T00:00:00 TAI — a uniform count
//! with no leap seconds, so epoch arithmetic is plain integer arithmetic. The messiness (leap
//! seconds, scale offsets, calendars) lives only at the boundary, in the constructors and accessors.
//!
//! Four time scales convert at the API — [`TimeScale::Tai`], [`Gps`](TimeScale::Gps),
//! [`Tt`](TimeScale::Tt) and [`Utc`](TimeScale::Utc):
//!
//! - `TT = TAI + 32.184 s` (Terrestrial Time — the scale ephemerides and J2000 are defined on).
//! - `GPS = TAI − 19 s` (GPS Time — no leap seconds; `GPS − UTC` grows as leaps are added).
//! - `UTC = TAI − (leap seconds)`, a step function tabulated from the IERS bulletins (10 s at
//!   1972-01-01 up to 37 s at 2017-01-01). Unix time is UTC seconds since 1970.
//!
//! ```
//! use geo4d::{Epoch, TimeScale};
//!
//! // The J2000.0 epoch: 2000-01-01T12:00:00 TT = Julian Date 2451545.0 (TT).
//! let j2000 = Epoch::from_gregorian(TimeScale::Tt, 2000, 1, 1, 12, 0, 0.0);
//! assert!((j2000.julian_date(TimeScale::Tt) - 2_451_545.0).abs() < 1e-9);
//! assert_eq!(j2000, Epoch::J2000);
//!
//! // Same instant read on other scales, and as a calendar date.
//! assert!((j2000.leap_seconds() - 32.0).abs() < 1e-9);            // TAI−UTC in 2000 was 32 s
//! let (y, mo, d, h, mi, s) = j2000.to_gregorian_utc();           // ...so 11:58:55.816 UTC
//! assert_eq!((y, mo, d, h, mi), (2000, 1, 1, 11, 58));
//! assert!((s - 55.816).abs() < 1e-3);
//! ```

use core::ops::{Add, Sub};

/// A time scale that [`Epoch`]s convert to and from at the API. The internal storage is always TAI.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum TimeScale {
    /// International Atomic Time — the uniform reference; no leap seconds.
    Tai,
    /// GPS Time — `TAI − 19 s`; no leap seconds. Epoch 1980-01-06T00:00:00 UTC.
    Gps,
    /// Terrestrial Time — `TAI + 32.184 s`; the scale J2000 and planetary ephemerides use.
    Tt,
    /// Coordinated Universal Time — `TAI − (leap seconds)`. Unix time is UTC seconds since 1970.
    Utc,
}

const NS_PER_SEC: i128 = 1_000_000_000;
const SEC_PER_DAY: i64 = 86_400;
/// `TT − TAI`, nanoseconds (32.184 s exactly).
const TT_MINUS_TAI_NS: i128 = 32_184_000_000;
/// `TAI − GPS`, nanoseconds (19 s exactly).
const TAI_MINUS_GPS_NS: i128 = 19 * NS_PER_SEC;
/// GPS epoch (1980-01-06T00:00:00 UTC) as Unix seconds.
const GPS_EPOCH_UNIX: i64 = 315_964_800;
/// Julian Date of the Unix epoch 1970-01-01T00:00:00.
const UNIX_JD: f64 = 2_440_587.5;

// ── The IERS leap-second table ─────────────────────────────────────────────────────────────────────

/// `(year, month, day, TAI−UTC after this UTC instant)` — every leap step from 1972 to 2017. Before
/// 1972-01-01 the value 10 is used (pre-1972 rubber-band UTC is not modelled). As of the 2035 CGPM
/// resolution no leap second has been introduced since 2017-01-01; update this table if one is.
const LEAP_TABLE: [(i64, i64, i64, i64); 28] = [
    (1972, 1, 1, 10),
    (1972, 7, 1, 11),
    (1973, 1, 1, 12),
    (1974, 1, 1, 13),
    (1975, 1, 1, 14),
    (1976, 1, 1, 15),
    (1977, 1, 1, 16),
    (1978, 1, 1, 17),
    (1979, 1, 1, 18),
    (1980, 1, 1, 19),
    (1981, 7, 1, 20),
    (1982, 7, 1, 21),
    (1983, 7, 1, 22),
    (1985, 7, 1, 23),
    (1988, 1, 1, 24),
    (1990, 1, 1, 25),
    (1991, 1, 1, 26),
    (1992, 7, 1, 27),
    (1993, 7, 1, 28),
    (1994, 7, 1, 29),
    (1996, 1, 1, 30),
    (1997, 7, 1, 31),
    (1999, 1, 1, 32),
    (2006, 1, 1, 33),
    (2009, 1, 1, 34),
    (2012, 7, 1, 35),
    (2015, 7, 1, 36),
    (2017, 1, 1, 37),
];

/// Days from 1970-01-01 to the civil date `(y, m, d)` in the proleptic Gregorian calendar. Exact for
/// all dates (Howard Hinnant's algorithm).
const fn days_from_civil(y: i64, m: i64, d: i64) -> i64 {
    let y = if m <= 2 { y - 1 } else { y };
    let era = (if y >= 0 { y } else { y - 399 }) / 400;
    let yoe = y - era * 400; // [0, 399]
    let doy = (153 * (if m > 2 { m - 3 } else { m + 9 }) + 2) / 5 + d - 1; // [0, 365]
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy; // [0, 146096]
    era * 146097 + doe - 719_468
}

/// The civil date `(y, m, d)` a given number of days from 1970-01-01 (inverse of [`days_from_civil`]).
const fn civil_from_days(z: i64) -> (i64, i64, i64) {
    let z = z + 719_468;
    let era = (if z >= 0 { z } else { z - 146_096 }) / 146_097;
    let doe = z - era * 146_097; // [0, 146096]
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365; // [0, 399]
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100); // [0, 365]
    let mp = (5 * doy + 2) / 153; // [0, 11]
    let d = doy - (153 * mp + 2) / 5 + 1; // [1, 31]
    let m = if mp < 10 { mp + 3 } else { mp - 9 }; // [1, 12]
    (if m <= 2 { y + 1 } else { y }, m, d)
}

/// `TAI − UTC` (leap seconds) applicable to a UTC instant given as Unix seconds.
fn leap_secs_at_unix(unix_s: i64) -> i64 {
    let mut leap = LEAP_TABLE[0].3;
    for (y, m, d, l) in LEAP_TABLE {
        if unix_s >= days_from_civil(y, m, d) * SEC_PER_DAY {
            leap = l;
        }
    }
    leap
}

/// `TAI − UTC` (leap seconds) applicable to a TAI instant given as nanoseconds from the TAI origin.
fn leap_secs_at_tai(tai_ns: i128) -> i64 {
    let mut leap = LEAP_TABLE[0].3;
    for (y, m, d, l) in LEAP_TABLE {
        // The transition expressed on the TAI scale is the UTC instant plus the new offset.
        let tai_transition = (days_from_civil(y, m, d) * SEC_PER_DAY + l) as i128 * NS_PER_SEC;
        if tai_ns >= tai_transition {
            leap = l;
        }
    }
    leap
}

// ── Duration ───────────────────────────────────────────────────────────────────────────────────────

/// A signed span of time, stored as whole nanoseconds. Exact, `Ord`, and about ±5.4×10²¹ years of
/// range — so it never overflows or drifts the way `f64` seconds would over long baselines.
///
/// ```
/// use geo4d::Duration;
/// let a = Duration::from_hours(1.5);
/// assert_eq!(a, Duration::from_minutes(90.0));
/// assert_eq!(a.as_seconds(), 5400.0);
/// assert_eq!((a + Duration::from_seconds(600.0)).as_hours(), 1.6666666666666667);
/// ```
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Duration {
    nanos: i128,
}

impl Duration {
    /// The zero span.
    pub const ZERO: Duration = Duration { nanos: 0 };
    /// One second.
    pub const SECOND: Duration = Duration { nanos: NS_PER_SEC };
    /// One minute (60 s).
    pub const MINUTE: Duration = Duration { nanos: 60 * NS_PER_SEC };
    /// One hour (3600 s).
    pub const HOUR: Duration = Duration {
        nanos: 3_600 * NS_PER_SEC,
    };
    /// One day (86 400 s).
    pub const DAY: Duration = Duration {
        nanos: SEC_PER_DAY as i128 * NS_PER_SEC,
    };
    /// One week (7 days).
    pub const WEEK: Duration = Duration {
        nanos: 7 * SEC_PER_DAY as i128 * NS_PER_SEC,
    };
    /// One **Julian year** — 365.25 days exactly. The year that geodetic velocities (mm/yr, mas/yr)
    /// and ITRF epochs are quoted in, so it is the one to divide by when propagating coordinates.
    pub const JULIAN_YEAR: Duration = Duration {
        nanos: 31_557_600 * NS_PER_SEC,
    };

    /// From an exact number of nanoseconds.
    #[inline]
    pub const fn from_nanos(nanos: i128) -> Self {
        Duration { nanos }
    }

    /// From microseconds.
    #[inline]
    pub const fn from_micros(micros: i128) -> Self {
        Duration { nanos: micros * 1_000 }
    }

    /// From milliseconds.
    #[inline]
    pub const fn from_millis(millis: i128) -> Self {
        Duration {
            nanos: millis * 1_000_000,
        }
    }

    /// From a floating-point number of seconds (rounded to the nearest nanosecond).
    #[inline]
    pub fn from_seconds(seconds: f64) -> Self {
        Duration {
            nanos: (seconds * 1e9).round() as i128,
        }
    }

    /// From a floating-point number of minutes.
    #[inline]
    pub fn from_minutes(minutes: f64) -> Self {
        Duration::from_seconds(minutes * 60.0)
    }

    /// From a floating-point number of hours.
    #[inline]
    pub fn from_hours(hours: f64) -> Self {
        Duration::from_seconds(hours * 3_600.0)
    }

    /// From a floating-point number of days (86 400 s each).
    #[inline]
    pub fn from_days(days: f64) -> Self {
        Duration::from_seconds(days * SEC_PER_DAY as f64)
    }

    /// From a floating-point number of Julian years (365.25 days each).
    #[inline]
    pub fn from_julian_years(years: f64) -> Self {
        Duration::from_seconds(years * 31_557_600.0)
    }

    /// The exact whole-nanosecond count.
    #[inline]
    pub const fn total_nanoseconds(self) -> i128 {
        self.nanos
    }

    /// As seconds (may lose precision for very large spans).
    #[inline]
    pub fn as_seconds(self) -> f64 {
        self.nanos as f64 / 1e9
    }

    /// As minutes.
    #[inline]
    pub fn as_minutes(self) -> f64 {
        self.as_seconds() / 60.0
    }

    /// As hours.
    #[inline]
    pub fn as_hours(self) -> f64 {
        self.as_seconds() / 3_600.0
    }

    /// As days (86 400 s each).
    #[inline]
    pub fn as_days(self) -> f64 {
        self.as_seconds() / SEC_PER_DAY as f64
    }

    /// As Julian years (365.25 days each) — the unit geodetic rates are per, so this is what to
    /// multiply a velocity (m/yr) by when propagating across epochs.
    #[inline]
    pub fn as_julian_years(self) -> f64 {
        self.as_seconds() / 31_557_600.0
    }

    /// The absolute value.
    #[inline]
    pub const fn abs(self) -> Self {
        Duration {
            nanos: self.nanos.abs(),
        }
    }

    /// −1, 0 or 1 as the span is negative, zero or positive.
    #[inline]
    pub const fn signum(self) -> i32 {
        self.nanos.signum() as i32
    }

    /// Whether the span is zero.
    #[inline]
    pub const fn is_zero(self) -> bool {
        self.nanos == 0
    }
}

impl Add for Duration {
    type Output = Duration;
    #[inline]
    fn add(self, o: Duration) -> Duration {
        Duration {
            nanos: self.nanos + o.nanos,
        }
    }
}
impl Sub for Duration {
    type Output = Duration;
    #[inline]
    fn sub(self, o: Duration) -> Duration {
        Duration {
            nanos: self.nanos - o.nanos,
        }
    }
}
impl core::ops::Neg for Duration {
    type Output = Duration;
    #[inline]
    fn neg(self) -> Duration {
        Duration { nanos: -self.nanos }
    }
}
impl core::ops::Mul<f64> for Duration {
    type Output = Duration;
    #[inline]
    fn mul(self, k: f64) -> Duration {
        Duration {
            nanos: (self.nanos as f64 * k).round() as i128,
        }
    }
}
impl core::ops::Div<f64> for Duration {
    type Output = Duration;
    #[inline]
    fn div(self, k: f64) -> Duration {
        Duration {
            nanos: (self.nanos as f64 / k).round() as i128,
        }
    }
}
impl core::ops::Div<Duration> for Duration {
    type Output = f64;
    /// The dimensionless ratio of two spans.
    #[inline]
    fn div(self, o: Duration) -> f64 {
        self.nanos as f64 / o.nanos as f64
    }
}

// ── Epoch ──────────────────────────────────────────────────────────────────────────────────────────

/// One instant in time, to the nanosecond. Stored internally on the uniform **TAI** scale, so
/// comparisons and differences are exact integer operations regardless of leap seconds; scales,
/// calendars and Julian dates are produced only at the boundary.
///
/// `Epoch − Epoch` is a [`Duration`]; `Epoch ± Duration` is an [`Epoch`]. Epochs are `Ord`, so they
/// sort and compare directly.
///
/// ```
/// use geo4d::{Duration, Epoch, TimeScale};
///
/// let launch = Epoch::from_gregorian_utc(2026, 3, 14, 9, 30, 0.0);
/// let later = launch + Duration::from_minutes(10.0);
/// assert_eq!(later - launch, Duration::from_minutes(10.0));
/// assert!(later > launch);
///
/// // GPS week / time-of-week round-trips.
/// let (week, tow) = launch.gps_week_seconds();
/// assert_eq!(Epoch::from_gps_week_seconds(week, tow), launch);
/// ```
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Epoch {
    /// Nanoseconds from 1970-01-01T00:00:00 TAI.
    tai_ns: i128,
}

impl Epoch {
    /// The J2000.0 epoch: 2000-01-01T12:00:00 TT (Julian Date 2451545.0 TT), the origin most
    /// astronomy and geodesy reckon from.
    pub const J2000: Epoch = Epoch {
        // 2000-01-01T12:00:00 TT = 2000-01-01T11:59:27.816 TAI.
        tai_ns: (days_from_civil(2000, 1, 1) * SEC_PER_DAY + 43_200) as i128 * NS_PER_SEC - TT_MINUS_TAI_NS,
    };
    /// The Unix epoch, 1970-01-01T00:00:00 UTC.
    pub const UNIX_EPOCH: Epoch = Epoch {
        tai_ns: LEAP_TABLE[0].3 as i128 * NS_PER_SEC, // TAI−UTC = 10 s here (pre-1972 convention)
    };
    /// The GPS epoch, 1980-01-06T00:00:00 UTC (GPS week 0, time-of-week 0).
    pub const GPS_EPOCH: Epoch = Epoch {
        tai_ns: (GPS_EPOCH_UNIX as i128 + 19) * NS_PER_SEC,
    };

    /// From nanoseconds on the given scale, measured from that scale's 1970-01-01T00:00:00. The most
    /// precise constructor; the others reduce to this.
    pub fn from_scale_nanos(scale: TimeScale, nanos: i128) -> Self {
        let tai_ns = match scale {
            TimeScale::Tai => nanos,
            TimeScale::Tt => nanos - TT_MINUS_TAI_NS,
            TimeScale::Gps => nanos + TAI_MINUS_GPS_NS,
            // `nanos` here is a UTC/Unix-nanosecond count; add the leap seconds for that instant.
            TimeScale::Utc => {
                let unix_s = (nanos.div_euclid(NS_PER_SEC)) as i64;
                nanos + leap_secs_at_unix(unix_s) as i128 * NS_PER_SEC
            }
        };
        Epoch { tai_ns }
    }

    /// From a floating-point number of seconds on the given scale, from that scale's
    /// 1970-01-01T00:00:00. Convenience; large values carry only `f64` precision.
    #[inline]
    pub fn from_scale_seconds(scale: TimeScale, seconds: f64) -> Self {
        Epoch::from_scale_nanos(scale, (seconds * 1e9).round() as i128)
    }

    /// From Unix time — UTC seconds since 1970-01-01T00:00:00, the usual `time_t`.
    #[inline]
    pub fn from_unix_seconds(seconds: f64) -> Self {
        Epoch::from_scale_seconds(TimeScale::Utc, seconds)
    }

    /// From Unix time in exact nanoseconds.
    #[inline]
    pub fn from_unix_nanos(nanos: i128) -> Self {
        Epoch::from_scale_nanos(TimeScale::Utc, nanos)
    }

    /// From GPS seconds since the GPS epoch (1980-01-06T00:00:00 UTC).
    #[inline]
    pub fn from_gps_seconds(seconds: f64) -> Self {
        Epoch::from_scale_seconds(TimeScale::Gps, seconds + GPS_EPOCH_UNIX as f64)
    }

    /// From a GPS week number (continuous, **not** the 1024-week-rollover value) and time-of-week in
    /// seconds.
    #[inline]
    pub fn from_gps_week_seconds(week: u32, tow_seconds: f64) -> Self {
        Epoch::from_gps_seconds(week as f64 * 604_800.0 + tow_seconds)
    }

    /// From a calendar date and time on the given scale. Month and day are 1-based; `seconds` may be
    /// fractional. Exact to the nanosecond (integer path).
    pub fn from_gregorian(
        scale: TimeScale,
        year: i64,
        month: i64,
        day: i64,
        hour: i64,
        minute: i64,
        seconds: f64,
    ) -> Self {
        let whole =
            days_from_civil(year, month, day) * SEC_PER_DAY + hour * 3_600 + minute * 60 + seconds.trunc() as i64;
        let frac_ns = ((seconds - seconds.trunc()) * 1e9).round() as i128;
        Epoch::from_scale_nanos(scale, whole as i128 * NS_PER_SEC + frac_ns)
    }

    /// From a UTC calendar date and time (the common case).
    #[inline]
    pub fn from_gregorian_utc(year: i64, month: i64, day: i64, hour: i64, minute: i64, seconds: f64) -> Self {
        Epoch::from_gregorian(TimeScale::Utc, year, month, day, hour, minute, seconds)
    }

    /// From a Julian Date on the given scale.
    #[inline]
    pub fn from_julian_date(scale: TimeScale, jd: f64) -> Self {
        Epoch::from_scale_seconds(scale, (jd - UNIX_JD) * SEC_PER_DAY as f64)
    }

    /// From a Modified Julian Date (`MJD = JD − 2400000.5`) on the given scale.
    #[inline]
    pub fn from_mjd(scale: TimeScale, mjd: f64) -> Self {
        Epoch::from_julian_date(scale, mjd + 2_400_000.5)
    }

    /// From a decimal year such as ITRF's `2020.0`, in the **calendar-fraction** convention geodetic
    /// reference frames label their epochs with: `year = Y + (day-of-year − 1 + fraction of day) /
    /// (days in Y)`, on UTC. So `2020.0` is exactly 2020-01-01T00:00:00 UTC (GDA2020's reference
    /// epoch) and `2020.5` is mid-year. To *propagate* a coordinate across epochs, do not subtract
    /// two decimal years (the 365/366 divisor jumps at New Year); use
    /// [`Duration::as_julian_years`] of the epoch difference instead.
    pub fn from_decimal_year(year: f64) -> Self {
        let y = year.floor() as i64;
        let start = Epoch::from_gregorian_utc(y, 1, 1, 0, 0, 0.0);
        let next = Epoch::from_gregorian_utc(y + 1, 1, 1, 0, 0, 0.0);
        start + (next - start) * (year - y as f64)
    }

    /// The current instant from the system clock (UTC). Requires the standard library.
    pub fn now() -> Self {
        let d = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("system clock before 1970");
        Epoch::from_unix_nanos(d.as_nanos() as i128)
    }

    /// The span since another epoch (`self − other`).
    #[inline]
    pub const fn duration_since(self, other: Epoch) -> Duration {
        Duration {
            nanos: self.tai_ns - other.tai_ns,
        }
    }

    /// `TAI − UTC`, the leap-second count in effect at this instant, in seconds.
    #[inline]
    pub fn leap_seconds(self) -> f64 {
        leap_secs_at_tai(self.tai_ns) as f64
    }

    /// This instant as nanoseconds on the given scale, from that scale's 1970-01-01T00:00:00.
    #[inline]
    pub fn to_scale_nanos(self, scale: TimeScale) -> i128 {
        match scale {
            TimeScale::Tai => self.tai_ns,
            TimeScale::Tt => self.tai_ns + TT_MINUS_TAI_NS,
            TimeScale::Gps => self.tai_ns - TAI_MINUS_GPS_NS,
            TimeScale::Utc => self.tai_ns - leap_secs_at_tai(self.tai_ns) as i128 * NS_PER_SEC,
        }
    }

    /// This instant as seconds on the given scale, from that scale's 1970-01-01T00:00:00.
    #[inline]
    pub fn to_scale_seconds(self, scale: TimeScale) -> f64 {
        self.to_scale_nanos(scale) as f64 / 1e9
    }

    /// Unix time — UTC seconds since 1970-01-01T00:00:00.
    #[inline]
    pub fn to_unix_seconds(self) -> f64 {
        self.to_scale_seconds(TimeScale::Utc)
    }

    /// GPS seconds since the GPS epoch.
    #[inline]
    pub fn gps_seconds(self) -> f64 {
        self.to_scale_seconds(TimeScale::Gps) - GPS_EPOCH_UNIX as f64
    }

    /// GPS week number (continuous) and time-of-week in seconds.
    pub fn gps_week_seconds(self) -> (u32, f64) {
        let g = self.gps_seconds();
        let week = (g / 604_800.0).floor();
        ((week as u32), g - week * 604_800.0)
    }

    /// The Julian Date on the given scale.
    #[inline]
    pub fn julian_date(self, scale: TimeScale) -> f64 {
        UNIX_JD + self.to_scale_nanos(scale) as f64 / (1e9 * SEC_PER_DAY as f64)
    }

    /// The Modified Julian Date (`JD − 2400000.5`) on the given scale.
    #[inline]
    pub fn mjd(self, scale: TimeScale) -> f64 {
        self.julian_date(scale) - 2_400_000.5
    }

    /// The decimal year in the calendar-fraction convention (see
    /// [`from_decimal_year`](Self::from_decimal_year)), e.g. `2020.0` for 2020-01-01T00:00:00 UTC.
    /// Inverse of [`from_decimal_year`](Self::from_decimal_year). For epoch differences used in
    /// coordinate propagation prefer [`Duration::as_julian_years`].
    pub fn decimal_year(self) -> f64 {
        let (y, ..) = self.to_gregorian_utc();
        let start = Epoch::from_gregorian_utc(y, 1, 1, 0, 0, 0.0);
        let next = Epoch::from_gregorian_utc(y + 1, 1, 1, 0, 0, 0.0);
        y as f64 + (self - start) / (next - start)
    }

    /// The calendar date and time on the given scale as `(year, month, day, hour, minute, seconds)`,
    /// month and day 1-based, `seconds` fractional.
    pub fn to_gregorian(self, scale: TimeScale) -> (i64, u8, u8, u8, u8, f64) {
        let ns = self.to_scale_nanos(scale);
        let total_secs = ns.div_euclid(NS_PER_SEC) as i64;
        let frac = (ns.rem_euclid(NS_PER_SEC)) as f64 / 1e9;
        let days = total_secs.div_euclid(SEC_PER_DAY);
        let sod = total_secs.rem_euclid(SEC_PER_DAY);
        let (y, m, d) = civil_from_days(days);
        let hour = sod / 3_600;
        let minute = (sod % 3_600) / 60;
        let second = sod % 60;
        (y, m as u8, d as u8, hour as u8, minute as u8, second as f64 + frac)
    }

    /// The UTC calendar date and time (the common case).
    #[inline]
    pub fn to_gregorian_utc(self) -> (i64, u8, u8, u8, u8, f64) {
        self.to_gregorian(TimeScale::Utc)
    }

    /// Earth Rotation Angle (radians in `[0, 2π)`), the IAU 2000 linear function of UT1 that
    /// superseded the older sidereal-time series; UT1 is approximated by UTC here. This is the angle
    /// to rotate ECEF into a true-equator inertial frame at this instant.
    pub fn earth_rotation_angle_rad(self) -> f64 {
        let tu = self.julian_date(TimeScale::Utc) - 2_451_545.0;
        // IERS Conventions 2010, eq. (5.15): 0.7790572732640 turns at J2000 + 1.00273781191135448
        // turns per UT1 day (both shown here as the nearest f64).
        let turns = 0.779_057_273_264_0 + 1.002_737_811_911_354_6 * tu;
        turns.rem_euclid(1.0) * core::f64::consts::TAU
    }
}

impl Add<Duration> for Epoch {
    type Output = Epoch;
    #[inline]
    fn add(self, d: Duration) -> Epoch {
        Epoch {
            tai_ns: self.tai_ns + d.total_nanoseconds(),
        }
    }
}
impl Sub<Duration> for Epoch {
    type Output = Epoch;
    #[inline]
    fn sub(self, d: Duration) -> Epoch {
        Epoch {
            tai_ns: self.tai_ns - d.total_nanoseconds(),
        }
    }
}
impl Sub<Epoch> for Epoch {
    type Output = Duration;
    #[inline]
    fn sub(self, o: Epoch) -> Duration {
        self.duration_since(o)
    }
}

// ── TimeWindow ───────────────────────────────────────────────────────────────────────────────────

/// A span of time between two [`Epoch`]s, either end optionally open. The temporal half of 4D
/// geofencing and the search bound for conflict prediction: a restriction that only holds
/// `[start, end)`.
///
/// ```
/// use geo4d::{Duration, Epoch, TimeWindow};
/// let open = Epoch::from_gregorian_utc(2026, 7, 1, 19, 0, 0.0);
/// let game = TimeWindow::for_duration(open, Duration::from_hours(3.0));
/// assert!(game.contains(open + Duration::from_hours(1.0)));
/// assert!(!game.contains(open + Duration::from_hours(4.0)));
/// assert_eq!(game.duration(), Some(Duration::from_hours(3.0)));
/// ```
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TimeWindow {
    /// The start, inclusive; `None` = open (unbounded in the past).
    pub start: Option<Epoch>,
    /// The end, exclusive; `None` = open (unbounded in the future).
    pub end: Option<Epoch>,
}

impl TimeWindow {
    /// A window between two optional bounds.
    #[inline]
    pub const fn new(start: Option<Epoch>, end: Option<Epoch>) -> Self {
        TimeWindow { start, end }
    }

    /// A closed window `[start, end)`.
    #[inline]
    pub const fn between(start: Epoch, end: Epoch) -> Self {
        TimeWindow::new(Some(start), Some(end))
    }

    /// A window from `start` running for `duration`.
    #[inline]
    pub fn for_duration(start: Epoch, duration: Duration) -> Self {
        TimeWindow::between(start, start + duration)
    }

    /// Everything at or after `start`.
    #[inline]
    pub const fn from(start: Epoch) -> Self {
        TimeWindow::new(Some(start), None)
    }

    /// Everything before `end`.
    #[inline]
    pub const fn until(end: Epoch) -> Self {
        TimeWindow::new(None, Some(end))
    }

    /// The all-of-time window (both ends open).
    #[inline]
    pub const fn unbounded() -> Self {
        TimeWindow::new(None, None)
    }

    /// Whether an epoch is inside `[start, end)`.
    #[inline]
    pub fn contains(&self, epoch: Epoch) -> bool {
        self.start.is_none_or(|s| epoch >= s) && self.end.is_none_or(|e| epoch < e)
    }

    /// The epoch clamped into the window (unchanged if already inside; the nearer bound otherwise).
    #[inline]
    pub fn clamp(&self, epoch: Epoch) -> Epoch {
        let mut e = epoch;
        if let Some(s) = self.start {
            if e < s {
                e = s;
            }
        }
        if let Some(end) = self.end {
            if e > end {
                e = end;
            }
        }
        e
    }

    /// The length of the window, or `None` if either end is open.
    #[inline]
    pub fn duration(&self) -> Option<Duration> {
        match (self.start, self.end) {
            (Some(s), Some(e)) => Some(e - s),
            _ => None,
        }
    }

    /// The overlap of two windows (open ends propagate), or `None` if they are disjoint.
    pub fn intersect(&self, other: &TimeWindow) -> Option<TimeWindow> {
        let start = max_opt(self.start, other.start, true);
        let end = max_opt(self.end, other.end, false);
        if let (Some(s), Some(e)) = (start, end) {
            if s >= e {
                return None;
            }
        }
        Some(TimeWindow::new(start, end))
    }

    /// Whether two windows overlap at all.
    #[inline]
    pub fn overlaps(&self, other: &TimeWindow) -> bool {
        self.intersect(other).is_some()
    }
}

/// The later of two optional starts (`take_later = true`) or the earlier of two optional ends
/// (`take_later = false`); `None` means unbounded on that side and yields to the bounded one.
fn max_opt(a: Option<Epoch>, b: Option<Epoch>, take_later: bool) -> Option<Epoch> {
    match (a, b) {
        (Some(a), Some(b)) => Some(if (a > b) == take_later { a } else { b }),
        (x, None) | (None, x) => x,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn time_window_contains_clamp_and_intersect() {
        let t0 = Epoch::from_unix_seconds(1_000_000.0);
        let w = TimeWindow::for_duration(t0, Duration::from_hours(2.0));
        assert!(w.contains(t0) && w.contains(t0 + Duration::from_hours(1.0)));
        assert!(!w.contains(t0 + Duration::from_hours(2.0))); // end exclusive
        assert!(!w.contains(t0 - Duration::SECOND));
        assert_eq!(w.clamp(t0 - Duration::HOUR), t0);
        assert_eq!(w.clamp(t0 + Duration::from_hours(5.0)), t0 + Duration::from_hours(2.0));
        assert_eq!(w.clamp(t0 + Duration::HOUR), t0 + Duration::HOUR);
        // Intersection with an overlapping and a disjoint window.
        let later = TimeWindow::from(t0 + Duration::HOUR);
        assert_eq!(w.intersect(&later).unwrap().duration(), Some(Duration::HOUR));
        assert!(w.intersect(&TimeWindow::from(t0 + Duration::from_hours(3.0))).is_none());
        assert!(TimeWindow::unbounded().contains(t0));
        assert_eq!(TimeWindow::unbounded().intersect(&w), Some(w));
    }

    #[test]
    fn civil_day_algorithm_round_trips() {
        for (y, m, d) in [
            (1970, 1, 1),
            (2000, 1, 1),
            (2026, 9, 14),
            (1972, 6, 30),
            (1600, 12, 31),
            (2400, 2, 29),
        ] {
            let z = days_from_civil(y, m, d);
            assert_eq!(civil_from_days(z), (y, m, d), "{y}-{m}-{d}");
        }
        assert_eq!(days_from_civil(1970, 1, 1), 0);
        assert_eq!(days_from_civil(2000, 1, 1), 10_957);
        assert_eq!(days_from_civil(1980, 1, 6) * SEC_PER_DAY, GPS_EPOCH_UNIX);
    }

    #[test]
    fn leap_seconds_step_correctly() {
        // Integer-valued and therefore exact.
        assert_eq!(Epoch::from_gregorian_utc(1971, 1, 1, 0, 0, 0.0).leap_seconds(), 10.0);
        assert_eq!(Epoch::from_gregorian_utc(1972, 1, 1, 0, 0, 0.0).leap_seconds(), 10.0);
        assert_eq!(Epoch::from_gregorian_utc(1972, 7, 1, 0, 0, 0.0).leap_seconds(), 11.0);
        assert_eq!(Epoch::from_gregorian_utc(2000, 1, 1, 0, 0, 0.0).leap_seconds(), 32.0);
        assert_eq!(Epoch::from_gregorian_utc(2017, 1, 1, 0, 0, 0.0).leap_seconds(), 37.0);
        assert_eq!(Epoch::from_gregorian_utc(2026, 1, 1, 0, 0, 0.0).leap_seconds(), 37.0);
        // One second before a step still reads the old value.
        assert_eq!(
            (Epoch::from_gregorian_utc(1972, 7, 1, 0, 0, 0.0) - Duration::SECOND).leap_seconds(),
            10.0
        );
    }

    #[test]
    fn scale_offsets_are_exact_to_the_nanosecond() {
        let e = Epoch::from_gregorian_utc(2020, 6, 15, 12, 0, 0.0);
        // Compared in integer nanoseconds, the scale offsets are exact — no float drift.
        // TT − TAI = 32.184 s, TAI − GPS = 19 s, TAI − UTC = 37 s (2020), GPS − UTC = 18 s.
        assert_eq!(
            e.to_scale_nanos(TimeScale::Tt) - e.to_scale_nanos(TimeScale::Tai),
            32_184_000_000
        );
        assert_eq!(
            e.to_scale_nanos(TimeScale::Tai) - e.to_scale_nanos(TimeScale::Gps),
            19_000_000_000
        );
        assert_eq!(
            e.to_scale_nanos(TimeScale::Tai) - e.to_scale_nanos(TimeScale::Utc),
            37_000_000_000
        );
        assert_eq!(
            e.to_scale_nanos(TimeScale::Gps) - e.to_scale_nanos(TimeScale::Utc),
            18_000_000_000
        );
    }

    #[test]
    fn j2000_is_consistent_everywhere() {
        let j = Epoch::J2000;
        assert!(close(j.julian_date(TimeScale::Tt), 2_451_545.0, 1e-9));
        assert!(close(j.mjd(TimeScale::Tt), 51_544.5, 1e-9));
        // J2000 in UTC is 11:58:55.816 (TT is ahead of UTC by 32.184 + 32 s = 64.184 s).
        let (y, mo, d, h, mi, s) = j.to_gregorian_utc();
        assert_eq!((y, mo, d, h, mi), (2000, 1, 1, 11, 58));
        assert!(close(s, 55.816, 1e-6), "{s}");
        assert_eq!(j, Epoch::from_gregorian(TimeScale::Tt, 2000, 1, 1, 12, 0, 0.0));
    }

    #[test]
    fn unix_and_gps_reference_values() {
        // 2000-01-01T12:00:00 UTC is Unix 946 728 000 and JD(UTC) 2451545.0.
        let noon2000 = Epoch::from_gregorian_utc(2000, 1, 1, 12, 0, 0.0);
        assert!(close(noon2000.to_unix_seconds(), 946_728_000.0, 1e-6));
        assert!(close(noon2000.julian_date(TimeScale::Utc), 2_451_545.0, 1e-9));
        assert_eq!(Epoch::from_unix_seconds(0.0), Epoch::UNIX_EPOCH);
        // The GPS epoch: week 0, tow 0, gps_seconds 0.
        assert_eq!(Epoch::GPS_EPOCH, Epoch::from_gregorian_utc(1980, 1, 6, 0, 0, 0.0));
        assert!(close(Epoch::GPS_EPOCH.gps_seconds(), 0.0, 1e-9));
        let (w, tow) = Epoch::GPS_EPOCH.gps_week_seconds();
        assert_eq!((w, tow), (0, 0.0));
        // A known GPS week: 2020-06-15 is week 2110.
        let (w, _) = Epoch::from_gregorian_utc(2020, 6, 15, 0, 0, 0.0).gps_week_seconds();
        assert_eq!(w, 2110);
    }

    #[test]
    fn round_trips_hold_to_the_nanosecond() {
        let e = Epoch::from_gregorian_utc(2026, 9, 14, 17, 33, 12.123_456_789);
        let (y, mo, d, h, mi, s) = e.to_gregorian_utc();
        assert_eq!((y, mo, d, h, mi), (2026, 9, 14, 17, 33));
        assert!(close(s, 12.123_456_789, 1e-9), "{s}");
        // Nanosecond arithmetic is exact and Ord works.
        let later = e + Duration::from_nanos(1);
        assert!(later > e && (later - e) == Duration::from_nanos(1));
        // Scale and calendar round-trips.
        for scale in [TimeScale::Tai, TimeScale::Gps, TimeScale::Tt, TimeScale::Utc] {
            let (y, mo, d, h, mi, s) = e.to_gregorian(scale);
            let back = Epoch::from_gregorian(scale, y, mo as i64, d as i64, h as i64, mi as i64, s);
            assert!((back - e).abs() <= Duration::from_nanos(1), "{scale:?}");
        }
    }

    #[test]
    fn decimal_year_uses_the_calendar_convention() {
        // ITRF/GDA labels: 2020.0 is exactly 2020-01-01T00:00:00 UTC.
        let e = Epoch::from_decimal_year(2020.0);
        assert_eq!(e, Epoch::from_gregorian_utc(2020, 1, 1, 0, 0, 0.0));
        assert!(close(e.decimal_year(), 2020.0, 1e-9));
        // Mid-year round-trips, and 2000 is a leap year (366 days).
        assert!(close(Epoch::from_decimal_year(2020.5).decimal_year(), 2020.5, 1e-9));
        let mid = Epoch::from_decimal_year(2000.5);
        let (_, mo, d, ..) = mid.to_gregorian_utc();
        assert_eq!((mo, d), (7, 2)); // 183 days after 1 Jan of a 366-day year
    }

    #[test]
    fn era_matches_reference() {
        // At J2000.0 UT1 the ERA is exactly 0.7790572732640 turns (IERS Conventions 2010, eq. 5.15).
        let era = Epoch::from_julian_date(TimeScale::Utc, 2_451_545.0).earth_rotation_angle_rad();
        assert!(close(era, 0.779_057_273_264_0 * core::f64::consts::TAU, 1e-12), "{era}");
        assert!((0.0..core::f64::consts::TAU).contains(&era));
    }

    #[test]
    fn durations_convert_and_compose() {
        assert_eq!(Duration::from_days(1.0), Duration::HOUR * 24.0);
        assert!(close(Duration::JULIAN_YEAR.as_days(), 365.25, 1e-9));
        assert!(close(Duration::from_hours(3.0) / Duration::from_hours(2.0), 1.5, 1e-12));
        assert_eq!(Duration::from_seconds(1.0).total_nanoseconds(), 1_000_000_000);
        assert_eq!((-Duration::SECOND).signum(), -1);
        assert_eq!(Duration::from_millis(1_500).as_seconds(), 1.5);
    }
}
