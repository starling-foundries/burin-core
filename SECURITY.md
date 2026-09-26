# Security

## Reporting

Report a vulnerability privately through GitHub: **Security → Report a vulnerability** on this
repository. Please do not open a public issue for it.

## In scope

- a verifier or check accepting a proof of a false statement;
- a check answering a question other than the one asked;
- a reader or verifier that panics, hangs, or allocates without bound on any input;
- a result that differs between platforms where SPEC.md defines it to the bit;
- two spellings of one record field that both read.

## Not in scope

The limits stated in [THREAT_MODEL.md](THREAT_MODEL.md): the truth of supplied points and
instants, the authenticity of a fingerprint obtained out of band, withheld proofs, and privacy.

## Supported versions

The latest minor release receives fixes.
