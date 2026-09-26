//! A fingerprint read from arbitrary text re-spells to that text, and a check on an arbitrary
//! record answers or refuses, never panics.
#![no_main]
use burin_core::fingerprint::Fingerprint;
use burin_core::opening::OpeningRecord;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let Ok(text) = std::str::from_utf8(data) else { return };
    let (line, rest) = text.split_once('\n').unwrap_or((text, ""));
    let Ok(fp) = Fingerprint::parse(line, &[]) else { return };
    assert_eq!(fp.to_text(), line, "one spelling");
    let Ok(v) = serde_json::from_str::<serde_json::Value>(rest) else { return };
    let Ok(r) = OpeningRecord::from_json(&v) else { return };
    let cid = r.cid().unwrap_or(9);
    let _ = fp.check_cell(&r, cid);
    let _ = fp.check_point(&r, (cid % 360) as f64 - 180.0, (cid % 180) as f64 - 90.0);
    let _ = fp.check_instant(&r, cid as i64);
});
