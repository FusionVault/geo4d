//! 4D geofencing: a [`Volume4`] is a 3D volume that only exists for a stretch of time.
//!
//! [`Timed`] boxes any `geo3d` [`Volume`](geo3d::Volume) — a [`Sphere`](geo3d::Sphere),
//! [`Cone`](geo3d::Cone), [`Cylinder`](geo3d::Cylinder) or [`AltitudeBand`](geo3d::AltitudeBand) —
//! into a [`TimeWindow`], so a query is inside only when it is inside the shape **and** the window is
//! open. That is exactly a temporary flight restriction, a stadium no-fly bubble during a match, or a
//! drone corridor active for one delivery: `Timed<Cylinder>` is a temporary cylinder, `Timed<Sphere>`
//! a temporary bubble. A [`Trajectory`] can be tested against one to find whether — and when — it
//! first enters.

use crate::time::{Duration, Epoch, TimeWindow};
use crate::trajectory::Trajectory;
use crate::types::Geodetic4;
use geo3d::Volume;

/// Something that can decide whether a **timestamped** point is inside it — the 4D counterpart of
/// `geo3d`'s [`Volume`](geo3d::Volume).
pub trait Volume4 {
    /// Whether the point is inside in both space and time.
    fn contains(&self, point: Geodetic4) -> bool;
}

/// A 3D [`Volume`](geo3d::Volume) that is only active during a [`TimeWindow`] — the building block of
/// 4D geofencing.
///
/// ```
/// use geo4d::{AltitudeBand, Cylinder, Duration, Epoch, Geodetic, Geodetic4, TimeWindow, Timed, Trajectory, Volume4};
///
/// // A temporary cylinder over a stadium: 2 km radius, surface to 3 km, active for a 3-hour match.
/// let kickoff = Epoch::from_gregorian_utc(2026, 7, 1, 19, 0, 0.0);
/// let tfr = Timed::for_duration(
///     Cylinder::new(Geodetic::new(-33.85, 151.06, 0.0), 2_000.0, AltitudeBand::new(Some(0.0), Some(3_000.0))),
///     kickoff,
///     Duration::from_hours(3.0),
/// );
///
/// // A point inside the cylinder is only "inside" while the match is on.
/// let overhead = Geodetic::new(-33.85, 151.06, 1_500.0);
/// assert!(tfr.contains(Geodetic4::new(overhead, kickoff + Duration::from_hours(1.0))));
/// assert!(!tfr.contains(Geodetic4::new(overhead, kickoff + Duration::from_hours(5.0)))); // match over
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Timed<V> {
    /// The 3D volume.
    pub volume: V,
    /// The window during which it is active.
    pub window: TimeWindow,
}

impl<V> Timed<V> {
    /// A volume active during `window`.
    #[inline]
    pub const fn new(volume: V, window: TimeWindow) -> Self {
        Timed { volume, window }
    }

    /// A volume active at all times (an always-on geofence).
    #[inline]
    pub const fn always(volume: V) -> Self {
        Timed::new(volume, TimeWindow::unbounded())
    }

    /// A volume active from `start` for `duration`.
    #[inline]
    pub fn for_duration(volume: V, start: Epoch, duration: Duration) -> Self {
        Timed::new(volume, TimeWindow::for_duration(start, duration))
    }

    /// A volume active over an explicit `[start, end)`.
    #[inline]
    pub const fn between(volume: V, start: Epoch, end: Epoch) -> Self {
        Timed::new(volume, TimeWindow::between(start, end))
    }
}

impl<V: Volume> Timed<V> {
    /// Whether a timestamped point is inside in both space and time.
    #[inline]
    pub fn contains4(&self, point: Geodetic4) -> bool {
        self.window.contains(point.epoch) && self.volume.contains(point.position)
    }

    /// Whether a trajectory ever enters the active volume, sampling at `step`. A sampled test:
    /// features shorter than `step` can be missed, so choose `step` well below the time the path
    /// spends crossing the volume.
    #[inline]
    pub fn intersects(&self, trajectory: &Trajectory, step: Duration) -> bool {
        self.first_entry(trajectory, step).is_some()
    }

    /// The epoch at which a trajectory first enters the active volume, or `None` if it never does
    /// within the overlap of the window and the trajectory's own span. Sampled at `step`, then the
    /// crossing is refined by bisection to effectively the nanosecond.
    pub fn first_entry(&self, trajectory: &Trajectory, step: Duration) -> Option<Epoch> {
        let overlap = self.window.intersect(&trajectory.window())?;
        let start = overlap.start?;
        let end = overlap.end?;
        // The overlap is non-empty, so `end > start`; fall back to a single span if `step` is unusable.
        let step = if step.total_nanoseconds() > 0 && step <= end - start {
            step
        } else {
            end - start
        };

        let inside = |t: Epoch| -> bool {
            self.window.contains(t)
                && trajectory
                    .geodetic_at(t)
                    .is_some_and(|g| self.volume.contains(g.position))
        };

        let mut prev_out: Option<Epoch> = None;
        let mut t = start;
        loop {
            let tc = if t > end { end } else { t };
            if inside(tc) {
                return Some(match prev_out {
                    Some(out) => bisect_entry(out, tc, inside),
                    None => tc, // already inside at the first sample
                });
            }
            if tc == end {
                return None;
            }
            prev_out = Some(tc);
            t = t + step;
        }
    }
}

impl<V: Volume> Volume4 for Timed<V> {
    #[inline]
    fn contains(&self, point: Geodetic4) -> bool {
        self.contains4(point)
    }
}

/// Bisect a bracket whose start is outside and end is inside the volume, returning the crossing
/// epoch to ~nanosecond precision.
fn bisect_entry(mut out: Epoch, mut inside_epoch: Epoch, inside: impl Fn(Epoch) -> bool) -> Epoch {
    for _ in 0..40 {
        let mid = out + (inside_epoch - out) / 2.0;
        if mid == out || mid == inside_epoch {
            break;
        }
        if inside(mid) {
            inside_epoch = mid;
        } else {
            out = mid;
        }
    }
    inside_epoch
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Ecef4, StateVector, Waypoint};
    use geo3d::{AltitudeBand, Cylinder, Ecef, Geodetic, LocalFrame, Sphere};

    fn t0() -> Epoch {
        Epoch::from_gregorian_utc(2026, 7, 1, 19, 0, 0.0)
    }

    #[test]
    fn timed_volume_is_space_and_time() {
        let tfr = Timed::for_duration(
            Sphere::new(Geodetic::new(0.0, 0.0, 1_000.0), 5_000.0),
            t0(),
            Duration::from_hours(2.0),
        );
        let center = Geodetic::new(0.0, 0.0, 1_000.0);
        assert!(tfr.contains(Geodetic4::new(center, t0() + Duration::from_hours(1.0))));
        assert!(!tfr.contains(Geodetic4::new(center, t0() - Duration::SECOND))); // before it opens
        assert!(!tfr.contains(Geodetic4::new(center, t0() + Duration::from_hours(3.0)))); // after it closes
                                                                                          // Outside spatially, during the window: still out.
        assert!(!tfr.contains(Geodetic4::new(
            Geodetic::new(0.1, 0.0, 1_000.0),
            t0() + Duration::from_hours(1.0)
        )));
        // always() ignores time.
        let perma = Timed::always(Sphere::new(center, 5_000.0));
        assert!(perma.contains(Geodetic4::new(center, t0() + Duration::from_hours(100.0))));
    }

    #[test]
    fn trajectory_enters_only_while_active() {
        // A drone flying due east at 1 km altitude, passing over a geofenced cylinder at the origin.
        let start = Geodetic4::from_lat_lon(0.0, -0.2, 1_000.0, t0()).to_ecef4();
        let vel =
            LocalFrame::new(start.to_geodetic4().position).velocity_of_course(geo3d::Course::new(90.0, 200.0, 0.0));
        let s0 = StateVector::new(t0(), start.position, vel);
        // 400 s of flight → travels ~80 km east, crossing the origin around the midpoint.
        let traj = crate::Trajectory::from_states([s0, s0.at(t0() + Duration::from_seconds(400.0))]).unwrap();

        let cyl = Cylinder::new(
            Geodetic::new(0.0, 0.0, 0.0),
            3_000.0,
            AltitudeBand::new(Some(0.0), Some(3_000.0)),
        );
        // Active window covers the crossing → it enters ~96 s in (3 km before it reaches the origin,
        // which it passes at ~111 s: 0.2° ≈ 22.3 km east at 200 m/s).
        let active = Timed::for_duration(cyl, t0(), Duration::from_hours(1.0));
        let entry = active.first_entry(&traj, Duration::from_seconds(1.0)).unwrap();
        assert!(
            (entry - t0()).as_seconds() > 90.0 && (entry - t0()).as_seconds() < 100.0,
            "{}",
            (entry - t0()).as_seconds()
        );
        assert!(active.intersects(&traj, Duration::from_seconds(1.0)));

        // Same geofence, but only active AFTER the drone has gone → no intersection.
        let later = Timed::for_duration(cyl, t0() + Duration::from_hours(1.0), Duration::from_hours(1.0));
        assert!(!later.intersects(&traj, Duration::from_seconds(1.0)));
        assert!(later.first_entry(&traj, Duration::from_seconds(1.0)).is_none());
    }

    #[test]
    fn entry_epoch_is_refined_by_bisection() {
        // A straight run in ECEF crossing a sphere; the coarse-sampled entry is still sharp.
        let a = Ecef4::new(Ecef::new(6_378_137.0 - 10_000.0, 0.0, 0.0), t0());
        let b = Ecef4::new(
            Ecef::new(6_378_137.0 + 10_000.0, 0.0, 0.0),
            t0() + Duration::from_seconds(100.0),
        );
        let traj = crate::Trajectory::from_waypoints([Waypoint::from(a), Waypoint::from(b)]).unwrap();
        let sphere = Sphere::new(Geodetic::new(0.0, 0.0, 0.0), 1_000.0); // ±1 km around the surface point
        let geo = Timed::always(sphere);
        let entry = geo.first_entry(&traj, Duration::from_seconds(10.0)).unwrap();
        // The sphere spans x ∈ [a-1000, a+1000]; entry when x reaches 6378137-1000, i.e. at t ≈ 45 s.
        let g = traj.geodetic_at(entry).unwrap();
        assert!(sphere.contains(g.position));
        // Just before the refined entry we are outside.
        assert!(!sphere.contains(traj.geodetic_at(entry - Duration::from_millis(50)).unwrap().position));
    }
}
