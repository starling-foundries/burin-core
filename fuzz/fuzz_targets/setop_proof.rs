//! A parsed set-operation proof verifies without panicking and reads back what it writes.
#![no_main]
use burin_core::setops::SetOpProof;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let Ok(v) = serde_json::from_slice::<serde_json::Value>(data) else { return };
    let Ok(p) = SetOpProof::from_json(&v) else { return };
    let ok = p.verify();
    let _ = p.size();
    let again = SetOpProof::from_json(&p.to_json()).expect("a proof reads back what it writes");
    assert_eq!(again, p);
    assert_eq!(again.verify(), ok);
});
