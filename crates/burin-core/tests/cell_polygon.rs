//! Cell polygons for drawing: closed, counterclockwise, around their own nucleus, across the
//! antimeridian and around the poles, for every cell at levels 0-3 under three profiles.

use burin_core::geo::Shape;
use burin_core::hierarchy::{cid_to_suid, SPACE};
use burin_core::profile::Profile;
use geo::{Contains, LineString, Point, Polygon};

fn profiles() -> [Profile; 3] {
    [Profile::ogc(), Profile::burin_1(), Profile { lon_0_udeg: 0, north_square: 1, south_square: 2, ..Profile::default() }]
}

fn signed_area(ring: &[(f64, f64)]) -> f64 {
    ring.windows(2).map(|w| w[0].0 * w[1].1 - w[1].0 * w[0].1).sum::<f64>() / 2.0
}

#[test]
fn every_polygon_is_closed_counterclockwise_and_holds_its_nucleus() {
    for profile in profiles() {
        let grid = profile.grid();
        let (mut crossing, mut caps) = (0, 0);
        for level in 0..=3u32 {
            for c in (0..6u64).flat_map(|b| SPACE.descendants(9 + b, level)) {
                let path = SPACE.path(c).unwrap();
                let ring = grid.cell_polygon(&path, 5).unwrap();
                let name = cid_to_suid(c).unwrap();
                assert!(ring.len() >= 5 && ring.first() == ring.last(), "{name}: not a closed ring");
                assert!(signed_area(&ring) > 0.0, "{name}: not counterclockwise ({} under {})", signed_area(&ring), profile.describe());
                let poly = Polygon::new(LineString::from(ring.clone()), vec![]);
                let (mut lon, lat) = grid.nucleus(&path).unwrap();
                if grid.shape(&path) == Shape::Cap {
                    caps += 1;
                    let seam = ring[0].1;
                    let pole = if path[0] == 0 { 90.0 } else { -90.0 };
                    assert!(poly.contains(&Point::new(10.0, (seam + pole) / 2.0)), "{name}: cap misses its pole");
                    continue;
                }
                if ring.iter().any(|p| p.0 > 180.0) {
                    crossing += 1;
                    if lon < 0.0 {
                        lon += 360.0;
                    }
                }
                assert!(poly.contains(&Point::new(lon, lat)), "{name}: nucleus ({lon}, {lat}) outside its polygon");
            }
        }
        assert_eq!(caps, 8, "one cap per polar base cell per level");
        assert!(crossing > 0, "some cells cross the antimeridian under {}", profile.describe());
    }
}

#[test]
fn equatorial_cells_are_their_four_corners() {
    let grid = Profile::ogc().grid();
    let ring = grid.cell_polygon(&SPACE.path(SPACE.cid(&[3, 4, 5]).unwrap()).unwrap(), 9).unwrap();
    assert_eq!(ring.len(), 5);
    let polar = grid.cell_polygon(&SPACE.path(SPACE.cid(&[0, 1, 5]).unwrap()).unwrap(), 9).unwrap();
    assert_eq!(polar.len(), 4 * 8 + 1, "n points per edge on a polar cell");
}
