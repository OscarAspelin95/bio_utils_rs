use thiserror::Error;

#[derive(Debug, Error)]
pub enum BioError {
    #[error("Invalid parameter: {0}")]
    InvalidParameterError(String),

    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),

    #[error("Serialization error: {0}")]
    SerializationError(#[from] serde_json::Error),

    #[error("File has invalid extension: {0}")]
    InvalidFileExtensionError(String),

    #[error("File does not exist: {0}")]
    FileDoesNotExistError(String),

    #[error("Needletail failed to parse file: {0}")]
    NeedletailParseError(#[from] needletail::errors::ParseError),
}
