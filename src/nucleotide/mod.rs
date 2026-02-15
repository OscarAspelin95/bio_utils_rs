//! Nucleotide sequence operations and metrics.
//!
//! Provides functions for:
//! - Reverse complement and base conversions ([`reverse_complement`], [`error_to_phred`])
//! - Quality and composition metrics ([`gc_content`], [`nucleotide_counts`], [`mean_error_and_phred`])
//! - Shannon entropy ([`shannon_entropy`], [`nucleotide_probabilities`])
//! - Homopolymer detection ([`find_homopolymers`])
//! - Exact and fuzzy pattern search ([`search_exact`], [`search_fuzzy`])
//! - Static lookup tables ([`NT_LOOKUP`], [`PHRED_TO_ERROR`])

mod seq;
pub use seq::*;

mod metrics;
pub use metrics::*;

mod statics;
pub use statics::*;

mod entropy;
pub use entropy::*;

mod homopolymer;
pub use homopolymer::*;

mod search;
pub use search::*;
