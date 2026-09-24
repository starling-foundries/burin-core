//! A parsed opening record verifies without panicking, and it has one spelling: what it writes,
//! it reads back unchanged.
#![no_main]
use burin_core::opening::OpeningRecord;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let Ok(v) = serde_json::from_slice::<serde_json::Value>(data) else { return };
    let Ok(r) = OpeningRecord::from_json(&v) else { return };
    let ok = r.verify();
    let _ = r.cid();
    let again = OpeningRecord::from_json(&r.to_json()).expect("a record reads back what it writes");
    assert_eq!(again, r);
    assert_eq!(again.verify(), ok);
});
