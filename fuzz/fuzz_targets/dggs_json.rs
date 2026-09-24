//! Zone data reads without panicking, and a coverage it accepts writes back to zone data that
//! reads as the same coverage.
#![no_main]
use burin_core::profile::Profile;
use burin_core::zone_data::{from_dggs_json, to_dggs_json, Presence};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let Ok(v) = serde_json::from_slice::<serde_json::Value>(data) else { return };
    let p = Profile::ogc();
    for presence in [Presence::NonNull, Presence::NonZero] {
        let Ok(t) = from_dggs_json(&v, &p, None, None, presence) else { continue };
        let zone = burin_core::suid_to_cid(v["zoneId"].as_str().unwrap()).unwrap();
        let doc = to_dggs_json(&t, &p, zone, "f").expect("an accepted coverage writes back");
        let back = from_dggs_json(&doc, &p, None, None, Presence::NonNull).expect("and reads again");
        assert_eq!(back.root(), t.root());
    }
});
