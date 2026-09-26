//! Fingerprints and their checks: one spelling, and every way a proof can fail to answer the
//! question is refused with a reason, never answered.

use burin_core::fingerprint::{Axis, Fingerprint};
use burin_core::hierarchy::{suid_to_cid, SPACE};
use burin_core::opening::{open_path, OpeningRecord};
use burin_core::profile::Profile;
use burin_core::setops::{holds, relation_transcript, SetRelation};
use burin_core::time::{self, EPOCH_TICK};
use burin_core::tree::Tree;
use burin_core::Error;
use std::sync::Arc;

fn space(d: u32, suids: &[&str], p: &Profile) -> Tree {
    Tree::from_cells(SPACE, d, suids.iter().map(|s| suid_to_cid(s).unwrap()), Arc::new(p.ctx(d).unwrap())).unwrap()
}

fn hours(spans: &[(u64, u64)]) -> Tree {
    let shifted: Vec<(u64, u64)> = spans.iter().map(|&(a, b)| (EPOCH_TICK + a, EPOCH_TICK + b)).collect();
    time::intervals_tree(&shifted, Arc::new(time::ctx(&Profile::ogc()).unwrap())).unwrap()
}

fn record(t: &Tree, cell: u64) -> OpeningRecord {
    OpeningRecord::new(t, open_path(t, cell).unwrap().expect("a constant node"))
}

fn is_refused<T: std::fmt::Debug>(r: burin_core::Result<T>, why: &str) {
    match r {
        Err(Error::InvalidProof(m)) => assert!(m.contains(why), "refused for {m:?}, expected {why:?}"),
        other => panic!("expected a refusal ({why}), got {other:?}"),
    }
}

#[test]
fn a_fingerprint_has_one_spelling() {
    let p = Profile::ogc();
    let t = space(4, &["Q4"], &p);
    let fp = Fingerprint::of(&t, &p).unwrap();
    let text = fp.to_text();
    assert!(text.starts_with("burin:ogc-rhealpix:space:4:"));
    assert_eq!(Fingerprint::parse(&text, &[]).unwrap(), fp);
    let root = &text[text.len() - 64..];
    for bad in [
        text.to_uppercase(),
        text.replace(":4:", ":04:"),
        text.replace(":4:", ":+4:"),
        text.replace("ogc-rhealpix", "OGC-rhealpix"),
        format!("burin:ogc-rhealpix:space:19:{root}"),
        format!("burin:ogc-rhealpix:time:61:{root}"),
        format!("burin:ogc-rhealpix:space:4:{}", &root[..62]),
        format!("burin:ogc-rhealpix:space:4:{root}:"),
        format!("brin:ogc-rhealpix:space:4:{root}"),
        format!("burin:somewhere:space:4:{root}"),
    ] {
        assert!(Fingerprint::parse(&bad, &[]).is_err(), "{bad:?} parsed");
    }
    let tf = Fingerprint::of(&hours(&[(0, 10)]), &p).unwrap();
    assert_eq!(tf.axis, Axis::Time);
    assert!(tf.to_text().starts_with("burin:ogc-rhealpix:time:"));
    assert_eq!(Fingerprint::parse(&tf.to_text(), &[]).unwrap(), tf);
    // a profile that is not registered is named by its id, and must be supplied to be read
    let odd = Profile { lon_0_udeg: 10_000_000, ..Profile::ogc() };
    let of = Fingerprint::of(&space(4, &["Q4"], &odd), &odd).unwrap();
    assert!(of.to_text().contains(&odd.id_hex()));
    assert!(Fingerprint::parse(&of.to_text(), &[]).is_err());
    assert_eq!(Fingerprint::parse(&of.to_text(), &[odd]).unwrap(), of);
}

#[test]
fn a_cell_is_answered_only_by_a_proof_about_it() {
    let p = Profile::ogc();
    let t = space(4, &["Q4", "N0"], &p);
    let fp = Fingerprint::of(&t, &p).unwrap();
    let (inside, outside) = (suid_to_cid("Q412").unwrap(), suid_to_cid("P000").unwrap());
    assert!(fp.check_cell(&record(&t, inside), inside).unwrap());
    assert!(!fp.check_cell(&record(&t, outside), outside).unwrap());
    assert!(fp.check_cell(&record(&t, suid_to_cid("Q4").unwrap()), suid_to_cid("Q4").unwrap()).unwrap(), "a whole coarse cell");

    is_refused(fp.check_cell(&record(&t, inside), outside), "is about cell");
    let other = space(4, &["Q4"], &p);
    is_refused(fp.check_cell(&record(&other, inside), inside), "another fingerprint");
    let deeper = space(5, &["Q4", "N0"], &p);
    is_refused(fp.check_cell(&record(&deeper, inside), inside), "depth 5");
    let tick = time::tick_cell(EPOCH_TICK + 3).unwrap();
    is_refused(fp.check_cell(&record(&hours(&[(0, 10)]), tick), inside), "H(2,1)");
    let mut forged = record(&t, inside);
    forged.root[31] ^= 1;
    is_refused(fp.check_cell(&forged, inside), "another fingerprint");
    let mut forged = record(&t, inside);
    if let burin_core::Opening::Entries(es) = &mut forged.opening {
        if let Some(burin_core::Entry::Hash(h)) = es.iter_mut().find(|e| matches!(e, burin_core::Entry::Hash(_))) {
            h[0] ^= 1;
        }
    }
    is_refused(fp.check_cell(&forged, inside), "does not verify");
    assert!(matches!(fp.check_cell(&record(&t, inside), suid_to_cid("Q41234").unwrap()), Err(Error::Invalid(_))), "a question deeper than the fingerprint");
}

#[test]
fn points_and_instants_are_read_on_the_fingerprint_s_own_grid_and_line() {
    let (ogc, b1) = (Profile::ogc(), Profile::burin_1());
    let (lon, lat) = (-73.9712, 40.7740);
    for p in [&ogc, &b1] {
        let cell = burin_core::zone::cell_from_point(p, lon, lat, 11).unwrap();
        let t = Tree::from_cells(SPACE, 11, [cell], Arc::new(p.ctx(11).unwrap())).unwrap();
        let fp = Fingerprint::of(&t, p).unwrap();
        assert!(fp.check_point(&record(&t, cell), lon, lat).unwrap());
        assert!(fp.check_instant(&record(&t, cell), 0).is_err(), "a space fingerprint has no instants");
    }
    let ogc_cell = burin_core::zone::cell_from_point(&ogc, lon, lat, 11).unwrap();
    let b1_cell = burin_core::zone::cell_from_point(&b1, lon, lat, 11).unwrap();
    assert_ne!(ogc_cell, b1_cell, "the grids differ, so the same place has different cells");
    let t = Tree::from_cells(SPACE, 11, [ogc_cell], Arc::new(ogc.ctx(11).unwrap())).unwrap();
    let as_b1 = Fingerprint { profile: b1.clone(), ..Fingerprint::of(&t, &ogc).unwrap() };
    is_refused(as_b1.check_point(&record(&t, ogc_cell), lon, lat), "is about cell");

    let open = hours(&[(1_000, 5_000)]);
    let fp = Fingerprint::of(&open, &ogc).unwrap();
    for (us, inside) in [(1_000, true), (4_999, true), (5_000, false), (-1, false)] {
        let tick = time::tick_cell(time::tick(&ogc, us).unwrap()).unwrap();
        assert_eq!(fp.check_instant(&record(&open, tick), us).unwrap(), inside, "{us}");
    }
    let tick = time::tick_cell(time::tick(&ogc, 1_000).unwrap()).unwrap();
    is_refused(fp.check_instant(&record(&open, tick), 1_001), "is about cell");
    assert!(fp.check_point(&record(&open, tick), lon, lat).is_err(), "a time fingerprint has no points");
}

#[test]
fn relations_are_decided_by_their_own_transcript_only() {
    let p = Profile::ogc();
    let sets = [space(3, &["Q4"], &p), space(3, &["Q41", "Q418"], &p), space(3, &["S"], &p), space(3, &["Q4", "N"], &p)];
    let fps: Vec<Fingerprint> = sets.iter().map(|t| Fingerprint::of(t, &p).unwrap()).collect();
    for (i, x) in sets.iter().enumerate() {
        for (j, y) in sets.iter().enumerate() {
            for rel in SetRelation::ALL {
                let tr = relation_transcript(x, y, rel).unwrap();
                assert_eq!(fps[i].check_relation(rel, &fps[j], tr.as_ref()).unwrap(), holds(x, y, rel).unwrap(), "{i} {j} {}", rel.name());
                if let Some(tr) = &tr {
                    is_refused(fps[i].check_relation(SetRelation::Equals, &fps[j], Some(tr)), "equals takes no transcript");
                    if fps[i] != fps[j] {
                        is_refused(fps[j].check_relation(rel, &fps[i], Some(tr)), "is not the");
                    }
                    let k = (i + 1) % sets.len();
                    if fps[k] != fps[i] {
                        is_refused(fps[k].check_relation(rel, &fps[j], Some(tr)), "is not the");
                    }
                } else {
                    is_refused(fps[i].check_relation(SetRelation::Within, &fps[j], None), "needs a transcript");
                }
            }
        }
    }
    let tf = Fingerprint::of(&hours(&[(0, 10)]), &p).unwrap();
    assert!(fps[0].check_relation(SetRelation::Equals, &tf, None).is_err(), "space and time are not comparable");
}
