//! Shared bioinformatics utilities for sequence I/O, nucleotide operations, and k-mer hashing.
//!
//! # Modules
//!
//! - [`io`] — Readers and writers for FASTQ/FASTA files (plain and gzip-compressed).
//! - [`nucleotide`] — Sequence operations, quality metrics, entropy, homopolymer detection, and pattern search.
//! - [`aminoacid`] - Nucleotide to aminoacid translations.
//! - [`kmers`] — K-mer encoding and FracMinHash sketching.
//! - [`errors`] — Shared error types used across the crate.

pub mod aminoacid;
pub mod errors;
pub mod io;
pub mod kmers;
pub mod nucleotide;

#[cfg(test)]
pub mod test_utils;
