//! File type classification for sequence files.

use crate::errors::BioError;

/// Compression type of a sequence file, inferred from its extension.
///
/// Recognized extensions:
/// - **Gzip**: `.fastq.gz`, `.fq.gz`, `.fasta.gz`, `.fa.gz`
/// - **Plain**: `.fastq`, `.fq`, `.fasta`, `.fa`
#[derive(Debug, PartialEq)]
pub enum SeqFileType {
    /// Gzip-compressed file.
    Gzip,
    /// Uncompressed plain-text file.
    Plain,
}

impl TryFrom<String> for SeqFileType {
    type Error = BioError;

    /// Determines the file type from the file path string.
    ///
    /// # Errors
    ///
    /// Returns [`BioError::InvalidFileExtensionError`] if the path does not end
    /// with a recognized FASTQ or FASTA extension.
    fn try_from(value: String) -> Result<Self, Self::Error> {
        if value.ends_with(".fastq.gz")
            || value.ends_with(".fq.gz")
            || value.ends_with(".fasta.gz")
            || value.ends_with(".fa.gz")
        {
            return Ok(Self::Gzip);
        }

        if value.ends_with(".fastq")
            || value.ends_with(".fq")
            || value.ends_with(".fasta")
            || value.ends_with(".fa")
        {
            return Ok(Self::Plain);
        }

        Err(BioError::InvalidFileExtensionError(value))
    }
}
