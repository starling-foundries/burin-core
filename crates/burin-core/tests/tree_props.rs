//! Structural properties of the tree, openings and set algebra (no fixtures needed).

use burin_core::hash::Ctx;
use burin_core::hierarchy::{suid_to_cid, Cid, Hierarchy, SPACE};
use burin_core::opening::{open_path, verify_opening, Claim, Entry, Opening, OpeningRecord};
use burin_core::setops::{difference, divergence, intersect, prove, union, Op, SetOpProof};
use burin_core::tree::Tree;
use std::collections::BTreeSet;
use std::sync::Arc;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

fn ctx(d: u32) -> Arc<Ctx> {
    Arc::new(Ctx::sha256(9, 6, d))
}

fn random_cells(rng: &mut Rng, h: &Hierarchy, d: u32, n: usize) -> Vec<Cid> {
    (0..n)
        .map(|_| {
            let level = rng.below(d as u64 + 1) as u32;
            let mut path = vec![rng.below(h.b as u64) as u32];
            for _ in 0..level {
                path.push(rng.below(h.a as u64) as u32);
            }
            h.cid(&path).unwrap()
        })
        .collect()
}

fn leaf_set(h: &Hierarchy, d: u32, cells: &[Cid]) -> BTreeSet<Cid> {
    let mut s = BTreeSet::new();
    for &c in cells {
        s.extend(h.descendants(c, d - h.level(c)));
    }
    s
}

#[test]
fn roll_up_is_an_identity() {
    let c = ctx(4);
    let q45 = suid_to_cid("Q45").unwrap();
    let whole = Tree::from_cells(SPACE, 4, [q45], c.clone()).unwrap();
    let spelled: Vec<Cid> = SPACE.descendants(q45, 2).collect();
    let fine = Tree::from_cells(SPACE, 4, spelled, c.clone()).unwrap();
    assert_eq!(whole.root(), fine.root());
    assert_eq!(fine.cells(), vec![q45]);
    // heterogeneous coverage does not collapse
    let partial = Tree::from_cells(SPACE, 4, SPACE.descendants(q45, 2).skip(1), c).unwrap();
    assert_ne!(partial.root(), whole.root());
}

#[test]
fn order_and_duplicates_do_not_matter() {
    let mut rng = Rng(0x9E3779B97F4A7C15);
    let c = ctx(4);
    for _ in 0..20 {
        let cells = random_cells(&mut rng, &SPACE, 4, 40);
        let a = Tree::from_cells(SPACE, 4, cells.iter().copied(), c.clone()).unwrap();
        let mut shuffled = cells.clone();
        for i in (1..shuffled.len()).rev() {
            let j = rng.below(i as u64 + 1) as usize;
            shuffled.swap(i, j);
        }
        shuffled.extend(cells.iter().take(7).copied());
        let b = Tree::from_cells(SPACE, 4, shuffled, c.clone()).unwrap();
        assert_eq!(a.root(), b.root());
        assert_eq!(a.leaves().into_iter().collect::<BTreeSet<_>>(), leaf_set(&SPACE, 4, &cells));
    }
}

#[test]
fn set_algebra_matches_a_set_oracle() {
    let mut rng = Rng(42);
    let d = 3;
    let c = ctx(d);
    for _ in 0..30 {
        let ca = random_cells(&mut rng, &SPACE, d, 12);
        let cb = random_cells(&mut rng, &SPACE, d, 12);
        let (sa, sb) = (leaf_set(&SPACE, d, &ca), leaf_set(&SPACE, d, &cb));
        let a = Tree::from_cells(SPACE, d, ca, c.clone()).unwrap();
        let b = Tree::from_cells(SPACE, d, cb, c.clone()).unwrap();
        let expect = |s: BTreeSet<Cid>| Tree::from_cells(SPACE, d, s, c.clone()).unwrap().root();
        assert_eq!(union(&a, &b).unwrap().root(), expect(&sa | &sb));
        assert_eq!(intersect(&a, &b).unwrap().root(), expect(&sa & &sb));
        assert_eq!(difference(&a, &b).unwrap().root(), expect(&sa - &sb));
        for op in [Op::Union, Op::Intersect, Op::Difference] {
            let p = prove(&a, &b, op).unwrap();
            assert!(p.verify(), "{op:?} transcript");
            let back = SetOpProof::from_json(&p.to_json()).unwrap();
            assert_eq!(back, p);
            assert!(back.verify());
            // tamper: flip the claimed result root
            let mut t = p.clone();
            t.root_c[0] ^= 1;
            assert!(!t.verify());
            // tamper: expand a decided node (a decided node must be a leaf step)
            let mut t = p.clone();
            if let Some(s) = t.steps.iter_mut().find(|s| s.children.is_none()) {
                s.children = Some(vec![s.clone(); 9]);
                assert!(!t.verify());
            }
        }
        let div = divergence(&a, &b).unwrap();
        assert_eq!(div.is_empty(), a.root() == b.root());
    }
}

#[test]
fn every_cell_opens_and_verifies() {
    let mut rng = Rng(7);
    let d = 3;
    let c = ctx(d);
    let cells = random_cells(&mut rng, &SPACE, d, 15);
    let t = Tree::from_cells(SPACE, d, cells, c.clone()).unwrap();
    let root = t.root();
    let mut n_full = 0;
    let mut n_empty = 0;
    let mut n_partial = 0;
    for b in 0..6u32 {
        for level in 0..=d {
            let mut path = vec![b];
            for _ in 0..level {
                path.push(rng.below(9) as u32);
            }
            let cid = SPACE.cid(&path).unwrap();
            match open_path(&t, cid).unwrap() {
                None => {
                    n_partial += 1;
                    assert!(t.node_at(cid).unwrap().0.is_branch());
                }
                Some(op) => {
                    assert!(verify_opening(&op, &root, d, &c));
                    assert_eq!(op.path().unwrap(), path);
                    let claim = op.terminal().unwrap();
                    assert_eq!(claim == Claim::Full, t.covers(cid).unwrap());
                    if claim == Claim::Full {
                        n_full += 1
                    } else {
                        n_empty += 1
                    }
                    let rec = OpeningRecord::new(&t, op.clone());
                    assert_eq!(rec.cid(), Some(cid));
                    let back = OpeningRecord::from_json(&rec.to_json()).unwrap();
                    assert_eq!(back, rec);
                    assert!(back.verify());
                    // tamper: flip the terminal claim
                    let flipped = flip(&op);
                    assert!(!verify_opening(&flipped, &root, d, &c), "flipped claim must not verify");
                    // tamper: a sibling hash
                    if let Opening::Entries(es) = &op {
                        let mut es2 = es.clone();
                        if let Some(Entry::Hash(h)) = es2.iter_mut().find(|e| matches!(e, Entry::Hash(_))) {
                            h[0] ^= 1;
                        }
                        assert!(!verify_opening(&Opening::Entries(es2), &root, d, &c));
                    }
                    // wrong depth
                    assert!(!verify_opening(&op, &root, d + 1, &Ctx::sha256(9, 6, d + 1)));
                }
            }
        }
    }
    assert!(n_full > 0 && n_empty > 0 && n_partial > 0, "{n_full} {n_empty} {n_partial}");
}

fn flip(op: &Opening) -> Opening {
    match op {
        Opening::Claim(Claim::Full) => Opening::Claim(Claim::Empty),
        Opening::Claim(Claim::Empty) => Opening::Claim(Claim::Full),
        Opening::Entries(es) => Opening::Entries(
            es.iter()
                .map(|e| match e {
                    Entry::Open(o) => Entry::Open(Box::new(flip(o))),
                    h => h.clone(),
                })
                .collect(),
        ),
    }
}

#[test]
fn malformed_openings_verify_false() {
    let c = ctx(2);
    let t = Tree::from_cells(SPACE, 2, [suid_to_cid("Q4").unwrap()], c.clone()).unwrap();
    let root = t.root();
    assert!(!verify_opening(&Opening::Claim(Claim::Full), &root, 2, &c));
    assert!(!verify_opening(&Opening::Entries(vec![]), &root, 2, &c));
    let too_few = Opening::Entries(vec![Entry::Hash([0; 32]); 5]);
    assert!(!verify_opening(&too_few, &root, 2, &c));
    let bad_arity = Opening::Entries(vec![Entry::Open(Box::new(Opening::Entries(vec![Entry::Hash([0; 32]); 8]))); 6]);
    assert!(!verify_opening(&bad_arity, &root, 2, &c));
    assert!(OpeningRecord::from_json(&serde_json::json!({"v": 2})).is_err());
}

#[test]
fn hierarchy_round_trips() {
    for s in ["N", "S8", "Q453", "O00000000000000", "R888888888888888"] {
        let cid = suid_to_cid(s).unwrap();
        assert_eq!(burin_core::hierarchy::cid_to_suid(cid).unwrap(), s);
        assert_eq!(SPACE.level(cid) as usize, s.len() - 1);
    }
    assert_eq!(suid_to_cid("Q453").unwrap(), (9 + 3) * 729 + 4 * 81 + 5 * 9 + 3);
    assert!(suid_to_cid("X1").is_err());
    assert!(SPACE.check(5, None).is_err());
    assert!(SPACE.check(15 * 9, None).is_err(), "base 15-9=6 is out of range");
}

#[test]
fn bulk_build_equals_cell_by_cell_insertion() {
    let mut rng = Rng(0x2545_f491_4f6c_dd1d);
    for d in [0u32, 1, 3, 5] {
        for n in [0usize, 1, 7, 60, 400] {
            let mut cells = random_cells(&mut rng, &SPACE, d, n);
            if n > 1 {
                cells.push(cells[0]); // duplicates
                cells.push(SPACE.parent(cells[1]).unwrap_or(cells[1])); // an ancestor of a listed cell
            }
            let bulk = Tree::from_cells(SPACE, d, cells.clone(), ctx(d)).unwrap();
            let mut one = Tree::empty(SPACE, d, ctx(d)).unwrap();
            for &c in &cells {
                one = one.set_full(c).unwrap();
            }
            assert_eq!(bulk.root(), one.root(), "d={d} n={n}");
            assert_eq!(bulk.cells(), one.cells(), "d={d} n={n}: canonical form");
        }
    }
    let planet: Vec<Cid> = (0..6).map(|b| SPACE.cid(&[b]).unwrap()).collect();
    assert!(Tree::from_cells(SPACE, 4, planet, ctx(4)).unwrap().is_full());
    assert!(Tree::from_cells(SPACE, 2, [SPACE.cid(&[0, 0, 0, 0]).unwrap()], ctx(2)).is_err(), "a cell below the depth");
}

#[test]
fn large_lists_build_the_same_tree_on_any_schedule() {
    // enough leaf ranges to take the parallel path when the `parallel` feature is on
    let mut rng = Rng(0x0dd_ba11);
    let d = 7;
    let cells = random_cells(&mut rng, &SPACE, d, 30_000);
    let bulk = Tree::from_cells(SPACE, d, cells.clone(), ctx(d)).unwrap();
    let one = cells.iter().fold(Tree::empty(SPACE, d, ctx(d)).unwrap(), |t, &c| t.set_full(c).unwrap());
    assert_eq!(bulk.root(), one.root());
    assert_eq!(bulk.cells(), one.cells());
    for _ in 0..3 {
        assert_eq!(Tree::from_cells(SPACE, d, cells.clone(), ctx(d)).unwrap().root(), bulk.root(), "a rebuild differs");
    }
}
