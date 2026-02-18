//! Crate-wide error types.

use thiserror::Error;

/// Unified error type for all `bio_utils_rs` operations.
#[derive(Debug, Error)]
pub enum BioError {
    /// A function argument was out of range or otherwise invalid.
    #[error("Invalid parameter: {0}")]
    InvalidParameterError(String),

    /// Wrapper around [`std::io::Error`].
    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),

    /// JSON serialization or deserialization failed.
    #[cfg(feature = "io")]
    #[error("Serialization error: {0}")]
    SerializationError(#[from] serde_json::Error),

    /// File path does not end with a recognized sequence file extension.
    #[error("File has invalid extension: {0}")]
    InvalidFileExtensionError(String),

    /// The specified file path does not exist on disk.
    #[error("File does not exist: {0}")]
    FileDoesNotExistError(String),

    /// Needletail failed to open or parse a sequence file.
    #[cfg(feature = "io")]
    #[error("Needletail failed to parse file: {0}")]
    NeedletailParseError(#[from] needletail::errors::ParseError),
}
