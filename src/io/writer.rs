use bio::io::fastq::Writer;
use flate2::Compression;
use flate2::write::GzEncoder;
use serde::Serialize;
use serde_json;
use std::io::Write;
use std::path::PathBuf;
use std::{fs::File, io::BufWriter};

use crate::errors::BioError;

pub fn write_json<T: Serialize>(outfile: Option<PathBuf>, s: T) -> Result<(), BioError> {
    let writer = get_bufwriter(outfile)?;
    serde_json::to_writer(writer, &s)?;

    Ok(())
}

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

/// Writer specifically for bio::io::fastq::Records:
/// * Some(outfile) -> will write gzip fastq.
/// * None -> will write plain fastq to stdout.
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
