//! Cell identifiers and points from arbitrary bits: parsing, topology and the point rule refuse
//! or answer, never panic, and what they answer is consistent.
#![no_main]
use burin_core::hierarchy::{cid_to_suid, suid_to_cid, SPACE};
use burin_core::profile::Profile;
use burin_core::zone::{cell_from_point, neighbours, position, subzone_at, subzone_index};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    if let Ok(s) = std::str::from_utf8(data) {
        if let Ok(c) = suid_to_cid(s) {
            assert_eq!(cid_to_suid(c).unwrap(), s, "a suid has one spelling");
        }
    }
    let word = |i: usize| data.get(i..i + 8).map(|b| u64::from_le_bytes(b.try_into().unwrap()));
    let p = Profile::ogc();
    if let Some(cid) = word(0) {
        if let Ok(pos) = position(&SPACE, cid) {
            let n = neighbours(&p, cid).expect("a cell has neighbours");
            assert!(n.iter().all(|&m| position(&SPACE, m).is_ok_and(|q| q.level == pos.level)));
            if let (Some(k), Some(parent)) = (data.get(8), SPACE.parent(cid)) {
                let i = subzone_index(&SPACE, parent, cid).unwrap();
                assert_eq!(subzone_at(&SPACE, parent, 1, i).unwrap(), cid);
                let _ = subzone_at(&SPACE, cid, (*k % 4) as u32, u64::from(*k) * 7);
            }
        }
    }
    if let (Some(lon), Some(lat), Some(&level)) = (word(0), word(8), data.get(16)) {
        let (lon, lat) = (f64::from_bits(lon), f64::from_bits(lat));
        if let Ok(c) = cell_from_point(&p, lon, lat, level as u32) {
            assert_eq!(SPACE.level(c), level as u32);
        }
    }
});
