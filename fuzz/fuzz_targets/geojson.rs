//! Any GeoJSON geometry either is refused or covers without panicking, and its cover is one
//! canonical tree: the cells it lists rebuild the same root.
#![no_main]
use burin_core::profile::Profile;
use burin_core::{cover, geometry_from_geojson, polyfill, Tree};
use libfuzzer_sys::fuzz_target;
use std::sync::Arc;

fuzz_target!(|data: &[u8]| {
    let Some((&r, json)) = data.split_first() else { return };
    let Ok(v) = serde_json::from_slice::<serde_json::Value>(json) else { return };
    let Ok(g) = geometry_from_geojson(&v) else { return };
    let (p, resolution) = (Profile::ogc(), (r % 5) as u32);
    let ctx = Arc::new(p.ctx(resolution).unwrap());
    let Ok(t) = cover(&g, resolution, &p, ctx.clone()) else { return };
    let cells = polyfill(&g, resolution, &p).expect("polyfill agrees with cover on what it accepts");
    let rebuilt = Tree::from_cells(p.hierarchy(), resolution, cells, ctx).unwrap();
    assert_eq!(rebuilt.root(), t.root());
});
