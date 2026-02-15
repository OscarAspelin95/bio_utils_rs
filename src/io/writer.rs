use bio::io::fastq::Writer;
use flate2::Compression;
use flate2::write::GzEncoder;
use serde::Serialize;
use std::io::Write;
use std::path::PathBuf;
use std::{fs::File, io::BufWriter};

use crate::errors::BioError;

/// Serializes `s` as JSON to a file or stdout.
///
/// Pass `Some(path)` to write to a file, or `None` for stdout.
/// Files ending in `.gz` are gzip-compressed automatically.
///
/// # Errors
///
/// Returns [`BioError`] on I/O or serialization failure.
pub fn write_json<T: Serialize>(outfile: Option<PathBuf>, s: T) -> Result<(), BioError> {
    let writer = get_bufwriter(outfile)?;
    serde_json::to_writer(writer, &s)?;

    Ok(())
}

/// Creates a buffered writer for a file or stdout.
///
/// Files ending in `.gz` are wrapped in a gzip encoder with fast compression.
/// Pass `None` to write to stdout.
///
/// # Errors
///
/// Returns [`BioError`] if the file cannot be created or has no extension.
pub fn get_bufwriter(outfile: Option<PathBuf>) -> Result<Box<dyn Write + Send>, BioError> {
    match outfile {
        Some(outfile) => {
            let f = File::create(&outfile)?;

            let extension = outfile.extension().map(|e| e.display().to_string()).ok_or(
                BioError::InvalidFileExtensionError(outfile.display().to_string()),
            )?;

            let writer = match extension.as_str() {
                "gz" => Box::new(BufWriter::new(GzEncoder::new(f, Compression::fast())))
                    as Box<dyn Write + Send>,
                _ => Box::new(BufWriter::new(f)) as Box<dyn Write + Send>,
            };

            Ok(writer)
        }
        None => {
            let writer = BufWriter::new(std::io::stdout());
            Ok(Box::new(writer))
        }
    }
}

/// Creates a [`bio::io::fastq::Writer`] for writing FASTQ records.
///
/// Output is always gzip-compressed when `Some(path)` is provided.
/// Pass `None` to write plain FASTQ to stdout.
///
/// # Errors
///
/// Returns [`BioError`] if the output file cannot be created.
pub fn bio_fastq_writer(outfile: Option<PathBuf>) -> Result<Writer<Box<dyn Write>>, BioError> {
    let writer: Box<dyn Write> = match outfile {
        Some(path) => {
            let f = File::create(path)?;
            Box::new(BufWriter::new(GzEncoder::new(f, Compression::fast())))
        }
        None => Box::new(BufWriter::new(std::io::stdout())),
    };

    Ok(Writer::new(writer))
}

/// Creates a [`bio::io::fasta::Writer`] for writing FASTA records.
///
/// Output is always gzip-compressed when `Some(path)` is provided.
/// Pass `None` to write plain FASTA to stdout.
///
/// # Errors
///
/// Returns [`BioError`] if the output file cannot be created.
pub fn bio_fasta_writer(
    outfile: Option<PathBuf>,
) -> Result<bio::io::fasta::Writer<Box<dyn Write>>, BioError> {
    let writer: Box<dyn Write> = match outfile {
        Some(path) => {
            let f = File::create(path)?;
            Box::new(BufWriter::new(GzEncoder::new(f, Compression::fast())))
        }
        None => Box::new(BufWriter::new(std::io::stdout())),
    };

    Ok(bio::io::fasta::Writer::new(writer))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use tempfile::TempDir;

    #[test]
    fn test_write_json_to_file() {
        let tmp_dir = TempDir::new().unwrap();
        let outfile = tmp_dir.path().join("test.json");

        let data = serde_json::json!({"key": "value", "num": 42});
        write_json(Some(outfile.clone()), &data).unwrap();

        let mut content = String::new();
        File::open(&outfile)
            .unwrap()
            .read_to_string(&mut content)
            .unwrap();
        let parsed: serde_json::Value = serde_json::from_str(&content).unwrap();
        assert_eq!(parsed["key"], "value");
        assert_eq!(parsed["num"], 42);
    }

    #[test]
    fn test_get_bufwriter_plain_file() {
        let tmp_dir = TempDir::new().unwrap();
        let outfile = tmp_dir.path().join("test.txt");

        let mut writer = get_bufwriter(Some(outfile.clone())).unwrap();
        writer.write_all(b"hello").unwrap();
        writer.flush().unwrap();
        drop(writer);

        let mut content = String::new();
        File::open(&outfile)
            .unwrap()
            .read_to_string(&mut content)
            .unwrap();
        assert_eq!(content, "hello");
    }

    #[test]
    fn test_get_bufwriter_gz_file() {
        let tmp_dir = TempDir::new().unwrap();
        let outfile = tmp_dir.path().join("test.gz");

        let mut writer = get_bufwriter(Some(outfile.clone())).unwrap();
        writer.write_all(b"compressed").unwrap();
        drop(writer);

        let metadata = std::fs::metadata(&outfile).unwrap();
        assert!(metadata.len() > 0);
    }

    #[test]
    fn test_bio_fastq_writer_to_file() {
        let tmp_dir = TempDir::new().unwrap();
        let outfile = tmp_dir.path().join("test.fastq.gz");

        let mut writer = bio_fastq_writer(Some(outfile.clone())).unwrap();
        writer
            .write("read1", None, b"ACGT", b"IIII")
            .unwrap();
        drop(writer);

        let metadata = std::fs::metadata(&outfile).unwrap();
        assert!(metadata.len() > 0);
    }

    #[test]
    fn test_bio_fasta_writer_to_file() {
        let tmp_dir = TempDir::new().unwrap();
        let outfile = tmp_dir.path().join("test.fasta.gz");

        let mut writer = bio_fasta_writer(Some(outfile.clone())).unwrap();
        writer.write("seq1", None, b"ACGT").unwrap();
        drop(writer);

        let metadata = std::fs::metadata(&outfile).unwrap();
        assert!(metadata.len() > 0);
    }
}
