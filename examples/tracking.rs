//! A 4D walkthrough: two aircraft in conflict, a time-boxed geofence, an intercept, tectonic drift
//! and a reference-frame change — everything bound to an absolute clock.
//!
//! `cargo run --example tracking`

use geo4d::{
    closest_approach, collision_probability, conflict, intercept, Course, Covariance3, Duration, Ecef, Ecef4, Epoch,
    Geodetic4, Helmert14, PlateMotion, ProtectedZone, StateVector, TimeScale, TimeWindow, Vec3,
};

fn utc(e: Epoch) -> String {
    let (y, mo, d, h, mi, s) = e.to_gregorian_utc();
    format!("{y:04}-{mo:02}-{d:02}T{h:02}:{mi:02}:{s:06.3}Z")
}

fn main() {
    let now = Epoch::from_gregorian_utc(2026, 3, 14, 9, 30, 0.0);
    println!(
        "scenario epoch: {}  (GPS week {:?})",
        utc(now),
        now.gps_week_seconds().0
    );

    // Two aircraft near Sydney, each a position + course-over-ground turned into an ECEF state.
    let a = Geodetic4::from_lat_lon(-33.80, 151.30, 2_500.0, now).state_from_course(Course::new(225.0, 240.0, 0.0));
    let b = Geodetic4::from_lat_lon(-33.95, 151.10, 2_450.0, now).state_from_course(Course::new(45.0, 250.0, 0.0));

    let cpa = closest_approach(a, b);
    println!(
        "\nclosest approach: {} — {:.0} m apart, closing at {:.0} m/s",
        utc(cpa.tca),
        cpa.distance_m,
        cpa.relative_speed_mps
    );

    let zone = ProtectedZone::new(9_260.0, 300.0); // 5 NM / 1000 ft
    match conflict(a, b, zone, TimeWindow::for_duration(now, Duration::from_minutes(10.0))) {
        Some(c) => println!(
            "CONFLICT: separation lost {} → {} ({:.0} s), worst at {}",
            utc(c.enters),
            utc(c.exits),
            c.duration().as_seconds(),
            utc(c.tca)
        ),
        None => println!("no loss of separation in the next 10 minutes"),
    }

    // A fixed-speed intercept of aircraft B from a ground point.
    let base = Ecef4::new(
        geo4d::WGS84.to_ecef(Geodetic4::from_lat_lon(-33.9, 151.0, 0.0, now).position),
        now,
    );
    if let Some(ix) = intercept(base, b, 340.0) {
        println!(
            "\nintercept of B at 340 m/s: {} (t+{:.0} s), {:.2} km downrange",
            utc(ix.epoch),
            ix.time_to_go.as_seconds(),
            ix.point.distance_to(base.position) / 1000.0
        );
    }

    // Tectonic drift and a reference-frame change on a ground station.
    let station_2020 = Geodetic4::from_lat_lon(-33.87, 151.21, 20.0, Epoch::from_decimal_year(2020.0)).to_ecef4();
    let station_2026 = PlateMotion::AUSTRALIA.propagate(station_2020, Epoch::from_decimal_year(2026.0));
    println!(
        "\nAustralian plate: the station moved {:.1} cm from 2020.0 to 2026.0",
        station_2026.position.distance_to(station_2020.position) * 100.0
    );
    let itrf2014 = Helmert14::itrf2020_to_itrf2014().apply(station_2026);
    println!(
        "ITRF2020 → ITRF2014 at epoch 2026.0: {:.1} mm",
        itrf2014.position.distance_to(station_2026.position) * 1000.0
    );

    // A satellite conjunction: probability of collision from each object's covariance.
    let epoch = Epoch::from_julian_date(TimeScale::Utc, 2_461_000.5);
    let sat = StateVector::new(epoch, Ecef::new(7_000_000.0, 0.0, 0.0), Vec3::new(0.0, 7_500.0, 0.0));
    let debris = StateVector::new(epoch, Ecef::new(7_000_120.0, 0.0, 0.0), Vec3::new(0.0, -7_500.0, 30.0));
    let pc = collision_probability(
        sat,
        Covariance3::isotropic(40.0),
        debris,
        Covariance3::isotropic(60.0),
        8.0,
    );
    println!("\nsatellite vs debris: probability of collision ≈ {pc:.2e}");
}
