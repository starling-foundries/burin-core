//! One error type for the crate. Verifiers never return it: a malformed proof verifies `false`.

#[derive(Debug, thiserror::Error, PartialEq)]
pub enum Error {
    #[error("{0}")]
    Invalid(String),
    #[error("({0:.20},{1:.20}) is not a point in the image of the projection")]
    OutOfImage(f64, f64),
}

pub type Result<T> = core::result::Result<T, Error>;

pub(crate) fn invalid<T>(msg: impl Into<String>) -> Result<T> {
    Err(Error::Invalid(msg.into()))
}
