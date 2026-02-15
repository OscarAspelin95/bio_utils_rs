//! Readers and writers for FASTQ and FASTA sequence files.
//!
//! Supports both plain-text and gzip-compressed files. When `None` is passed
//! as the file path, functions default to stdin (readers) or stdout (writers).

mod reader;
pub use reader::*;

mod writer;
pub use writer::*;

pub mod types;
