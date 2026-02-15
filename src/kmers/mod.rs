//! K-mer encoding and sketching.
//!
//! Provides a FracMinHash implementation for generating compact sequence
//! sketches from canonical (strand-aware) k-mers.

mod hash;
mod kmerize;
pub use kmerize::frac_min_hash;
