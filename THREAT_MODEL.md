# Threat model

What a burin-core check establishes, what it assumes, and what it does not establish.

## What a check shows

`Fingerprint.check_cell`, `check_point` and `check_instant` return `True` or `False`: the set
whose fingerprint is F does, or does not, contain that cell, or the cell of that point or
instant, on F's grid and at F's depth. `check_relation` returns whether the relation holds
between the sets of two fingerprints. Any proof that does not settle exactly that question
against exactly that fingerprint raises `InvalidProof`.

The low-level verifiers (`verify_opening`, `verify_setop`, `verify_relation`) check a proof
against its own root only. A caller who uses them must also compare the root with a fingerprint
it trusts and the cell with the one it asked about. The checks above do both.

## Assumptions

- **SHA-256 is collision resistant.** Every node hash commits to its children in order, under a
  domain tag (§4), and a fixed number of children means no length ambiguity. A proof that
  verifies for a false statement therefore yields two different inputs with the same node hash.
- **The verifier obtained the fingerprint authentically.** A fingerprint says what set; it does
  not say who published it. Until records are signed (SPEC-time.md §15), authenticity is carried
  by how the fingerprint was obtained.
- **The verifier computes, or trusts, the point or instant it asks about.** The check reads it on
  the fingerprint's own grid and time line, so a verifier cannot mismatch them.

## What an attacker holding a proof can and cannot do

| attempt | result |
|---|---|
| change any byte of a proof | the proof does not verify |
| present a valid proof for another cell, point or instant | refused: the proof is about another cell |
| present a proof against another set, depth, grid or axis | refused, with which of those differs |
| present a record with a second spelling of a field | refused on reading (SPEC §5) |
| send an oversized or malformed record | refused with bounded work: apertures, depths and listings are capped, and parsers are fuzzed |
| move a proof between identical siblings | verifies, and states something equally true (SPEC §5) |
| offer a transcript for another relation or another pair of sets | refused |

## What is not shown

- **Truth.** A point or instant is whatever was supplied; a set is whatever the publisher
  computed. A check shows consistency with a published fingerprint, not that anything happened.
- **Who, and when.** Nothing yet identifies a publisher or binds a time of issue.
- **Completeness.** A publisher may decline to prove a cell. Withholding cannot be detected.
- **Privacy.** Fingerprints are deterministic, so a small or guessable set can be found by trying
  candidates, and an opening shows which of its path's siblings are empty, full or identical.
  Hiding a set is outside this crate.

## The implementation

There is no `unsafe` code. Every reader and verifier is fuzzed (`fuzz/`), and every single edit of
real records is tested to be refused or to state something true (`wire_mutations.rs`). Point
lookup is defined to the bit and checked on five platforms (SPEC §6).

The soundness of the set-algebra transcripts (SPEC §5) rests on the argument above: each step's
three hashes chain to the three roots, decided steps follow from the constant ladders, and equal
hashes of two subtrees mean equal sets unless SHA-256 collides. That argument has not been
reviewed outside the project.
