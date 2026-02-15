use super::types::SeqFileType;
use crate::errors::BioError;
use bio::io::fastq::Reader;
use flate2::read::MultiGzDecoder;
use needletail::{FastxReader, parse_fastx_file, parse_fastx_stdin};
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

/// Validates that `path` exists and has a recognized sequence file extension.
fn validate_seq_file(path: &Path) -> Result<(&Path, SeqFileType), BioError> {
    if !path.exists() {
        return Err(BioError::FileDoesNotExistError(path.display().to_string()));
    }

    let file_type = SeqFileType::try_from(path.display().to_string())?;
    Ok((path, file_type))
}

/// Creates a [`bio::io::fastq::Reader`] for a FASTQ file.
///
/// Pass `Some(path)` for a file (plain or gzip), or `None` to read from stdin.
///
/// # Errors
///
/// Returns [`BioError`] if the file does not exist, has an unrecognized
/// extension, or cannot be opened.
pub fn bio_fastq_reader(
    fastq: Option<PathBuf>,
) -> Result<Reader<BufReader<Box<dyn Read + Send>>>, BioError> {
    let reader = match fastq {
        Some(fastq) => {
            let (fastq_file, file_type) = validate_seq_file(&fastq)?;

            let f = File::open(fastq_file)?;

            let reader: Box<dyn Read + Send> = match file_type {
                SeqFileType::Gzip => Box::new(MultiGzDecoder::new(f)),
                SeqFileType::Plain => Box::new(f),
            };

            reader
        }

        None => Box::new(std::io::stdin()),
    };

    Ok(Reader::new(reader))
}

/// Creates a [`bio::io::fasta::Reader`] for a FASTA file.
///
/// Pass `Some(path)` for a file (plain or gzip), or `None` to read from stdin.
///
/// # Errors
///
/// Returns [`BioError`] if the file does not exist, has an unrecognized
/// extension, or cannot be opened.
pub fn bio_fasta_reader(
    fasta: Option<PathBuf>,
) -> Result<bio::io::fasta::Reader<BufReader<Box<dyn Read + Send>>>, BioError> {
    let reader = match fasta {
        Some(fasta) => {
            let (fasta_file, file_type) = validate_seq_file(&fasta)?;

            let f = File::open(fasta_file)?;

            let reader: Box<dyn Read + Send> = match file_type {
                SeqFileType::Gzip => Box::new(MultiGzDecoder::new(f)),
                SeqFileType::Plain => Box::new(f),
            };

            reader
        }

        None => Box::new(std::io::stdin()),
    };

    Ok(bio::io::fasta::Reader::new(reader))
}

/// Creates a needletail [`FastxReader`] that handles both FASTQ and FASTA formats.
///
/// Pass `Some(path)` for a file, or `None` to read from stdin. Format is
/// auto-detected by needletail.
///
/// # Errors
///
/// Returns [`BioError`] if the file does not exist, has an unrecognized
/// extension, or needletail fails to parse it.
pub fn needletail_reader(path: Option<PathBuf>) -> Result<Box<dyn FastxReader>, BioError> {
    let reader = match path {
        Some(path) => {
            let (seq_file, _) = validate_seq_file(&path)?;
            parse_fastx_file(seq_file)?
        }
        None => parse_fastx_stdin()?,
    };

    Ok(reader)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::temp_seq_file;
    use rstest::*;

    #[rstest]
    #[case("valid.fastq.gz", SeqFileType::Gzip)]
    #[case("valid.fastq", SeqFileType::Plain)]
    #[case("valid.fq.gz", SeqFileType::Gzip)]
    #[case("valid.fq", SeqFileType::Plain)]
    #[case("valid.fasta.gz", SeqFileType::Gzip)]
    #[case("valid.fasta", SeqFileType::Plain)]
    #[case("valid.fa.gz", SeqFileType::Gzip)]
    #[case("valid.fa", SeqFileType::Plain)]
    fn test_valid_seq_file(#[case] file_name: &str, #[case] expected_type: SeqFileType) {
        let (_tmp_dir, tmp_file) = temp_seq_file(file_name);
        let (_, file_type) = validate_seq_file(&tmp_file).expect("file should pass validation");
        assert_eq!(file_type, expected_type);
    }

    #[rstest]
    #[case("invalid.txt")]
    #[case("invalid.csv")]
    fn test_invalid_seq_file(#[case] file_name: &str) {
        let (_tmp_dir, tmp_file) = temp_seq_file(file_name);
        let result = validate_seq_file(&tmp_file);
        assert!(result.is_err());
    }

    #[test]
    fn test_nonexistent_file() {
        let result = validate_seq_file(Path::new("/nonexistent/file.fastq"));
        assert!(result.is_err());
    }
}
