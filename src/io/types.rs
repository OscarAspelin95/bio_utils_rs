use crate::errors::BioError;

#[derive(Debug, PartialEq)]
pub enum FastqType {
    Gzip,
    Plain,
}

impl TryFrom<String> for FastqType {
    type Error = BioError;

    fn try_from(value: String) -> Result<Self, Self::Error> {
        if value.ends_with(".fastq.gz") || value.ends_with(".fq.gz") {
            return Ok(Self::Gzip);
        }

        if value.ends_with(".fastq") || value.ends_with(".fq") {
            return Ok(Self::Plain);
        }

        Err(BioError::InvalidFileExtensionError(value))
    }
}
