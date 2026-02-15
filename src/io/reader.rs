use super::types::FastqType;
use crate::errors::BioError;
use bio::io::fastq::Reader;
use flate2::read::MultiGzDecoder;
use needletail::{FastxReader, parse_fastx_file, parse_fastx_stdin};
use rstest::*;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};
use tempfile::TempDir;

fn validate_fastq(fastq: &Path) -> Result<(&Path, FastqType), BioError> {
    if !fastq.exists() {
        return Err(BioError::FileDoesNotExistError(fastq.display().to_string()));
    }

    let fastq_type = FastqType::try_from(fastq.display().to_string())?;
    Ok((fastq, fastq_type))
}

pub fn bio_fastq_reader(
    fastq: Option<PathBuf>,
) -> Result<Reader<BufReader<Box<dyn Read + Send>>>, BioError> {
    let reader = match fastq {
        Some(fastq) => {
            let (fastq_file, fastq_type) = validate_fastq(&fastq)?;

            let f = File::open(fastq_file)?;

            let reader: Box<dyn Read + Send> = match fastq_type {
                FastqType::Gzip => Box::new(MultiGzDecoder::new(f)),
                FastqType::Plain => Box::new(f),
            };

            reader
        }

        None => Box::new(std::io::stdin()),
    };

    Ok(Reader::new(reader))
}

pub fn needletail_fastq_reader(fastq: Option<PathBuf>) -> Result<Box<dyn FastxReader>, BioError> {
    let reader = match fastq {
        Some(fastq) => {
            let (fastq_file, _) = validate_fastq(&fastq)?;
            parse_fastx_file(fastq_file)?
        }
        None => parse_fastx_stdin()?,
    };

    Ok(reader)
}

#[rstest]
#[case("valid.fastq.gz".into(), FastqType::Gzip)]
#[case("valid.fastq".into(), FastqType::Plain)]
#[case("valid.fq.gz".into(), FastqType::Gzip)]
#[case("valid.fq".into(), FastqType::Plain)]
fn test_valid_fastq(#[case] file_name: String, #[case] expected_type: FastqType) {
    // Create tmp file for validation.
    let tmp_dir = TempDir::new().expect("failed to generate tmp dir.");
    let tmp_file = tmp_dir.path().join(file_name);
    let _ = File::create(&tmp_file).expect("failed to create file");

    let (_, fastq_type) = validate_fastq(&tmp_file).expect("file should pass validation");
    assert_eq!(fastq_type, expected_type);
}

#[rstest]
#[case("invalid.fasta.gz".into())]
fn test_invalid_fastq(#[case] file_name: String) {
    // Create tmp file for validation.
    // NOTE - we should create a separate function for generating temp dirs/files to avoid duplication.
    let tmp_dir = TempDir::new().expect("failed to generate tmp dir.");
    let tmp_file = tmp_dir.path().join(file_name);
    let _ = File::create(&tmp_file).expect("failed to create file");

    let result = validate_fastq(&tmp_file);
    assert!(result.is_err())
}
