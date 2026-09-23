use burin_core::hash::Ctx;
use burin_core::hierarchy::{suid_to_cid, SPACE};
use burin_core::index::{Index, Relation};
use burin_core::tree::Tree;
use std::sync::Arc;

fn t(cells: &[&str], ctx: &Arc<Ctx>) -> Tree {
    Tree::from_cells(SPACE, 4, cells.iter().map(|s| suid_to_cid(s).unwrap()), ctx.clone()).unwrap()
}

#[test]
fn query_fractions_and_relations() {
    let ctx = Arc::new(Ctx::sha256(9, 6, 4));
    let mut ix = Index::new();
    ix.add("a-whole-Q4", &t(&["Q4"], &ctx)).unwrap();
    ix.add("b-Q45-and-Q47", &t(&["Q45", "Q47"], &ctx)).unwrap();
    ix.add("c-same-as-b", &t(&["Q47", "Q45"], &ctx)).unwrap();
    ix.add("d-elsewhere", &t(&["R1"], &ctx)).unwrap();
    assert_eq!(ix.len(), 4);
    assert_eq!(ix.distinct_roots(), 3, "b and c share a root and a tree");
    let q = t(&["Q45"], &ctx);
    let hits = ix.query(&q).unwrap();
    let keys: Vec<&str> = hits.iter().map(|h| h.key.as_str()).collect();
    assert_eq!(keys, vec!["a-whole-Q4", "b-Q45-and-Q47", "c-same-as-b"]);
    let a = &hits[0];
    assert_eq!(a.relation, Relation::Contains);
    assert_eq!(a.fraction_of_query, 1.0);
    assert!((a.fraction_of_item - 1.0 / 9.0).abs() < 1e-12);
    let b = &hits[1];
    assert_eq!(b.relation, Relation::Contains);
    assert_eq!(b.fraction_of_item, 0.5);
    let q2 = t(&["Q4", "Q5"], &ctx);
    let hits = ix.query(&q2).unwrap();
    assert_eq!(hits[0].key, "a-whole-Q4");
    assert_eq!(hits[0].relation, Relation::Within);
    assert_eq!(hits[0].fraction_of_query, 0.5);
    assert_eq!(ix.query(&t(&["S0"], &ctx)).unwrap().len(), 0);
    assert_eq!(ix.query(&t(&["Q4"], &ctx)).unwrap()[0].relation, Relation::Equal);
}

#[test]
fn attention_counts_touched_cells() {
    let ctx = Arc::new(Ctx::sha256(9, 6, 4));
    let mut ix = Index::new();
    ix.add("x", &t(&["Q45"], &ctx)).unwrap();
    ix.add("y", &t(&["Q451", "Q47"], &ctx)).unwrap();
    ix.add("z", &t(&["R"], &ctx)).unwrap();
    let at = ix.attention(2).unwrap();
    assert_eq!(at[&suid_to_cid("Q45").unwrap()], 2, "x fully, y partially");
    assert_eq!(at[&suid_to_cid("Q47").unwrap()], 1);
    assert_eq!(at[&suid_to_cid("R00").unwrap()], 1, "a full base expands to every level-2 cell");
    assert_eq!(at.values().filter(|&&n| n == 1).count(), 81 + 1);
    assert!(ix.attention(5).is_err());
}
