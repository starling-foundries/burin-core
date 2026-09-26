//! One error type for the crate. The low-level verifiers (`verify_opening`, `SetOpProof::verify`)
//! return `false` on any defect; the checks of `fingerprint` return `InvalidProof` with the reason.

#[derive(Debug, thiserror::Error, PartialEq)]
pub enum Error {
    #[error("{0}")]
    Invalid(String),
    #[error("({0:.20},{1:.20}) is not a point in the image of the projection")]
    OutOfImage(f64, f64),
    /// The proof proves nothing about the question asked.
    #[error("invalid proof: {0}")]
    InvalidProof(String),
}

pub type Result<T> = core::result::Result<T, Error>;

pub(crate) fn invalid<T>(msg: impl Into<String>) -> Result<T> {
    Err(Error::Invalid(msg.into()))
}
