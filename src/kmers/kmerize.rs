use super::hash::mm_hash64;
use crate::errors::BioError;
use crate::nucleotide::NT_LOOKUP;
use std::collections::HashSet;

/// Computes a FracMinHash sketch of canonical k-mers from a DNA sequence.
///
/// Encodes each k-mer as a 2-bit packed `u64`, selects the canonical
/// (lexicographically smaller) orientation of forward/reverse complement,
/// and retains hashes that fall below `u64::MAX / ds_factor`.
///
/// Ambiguous bases (anything not `A`/`C`/`G`/`T`/`a`/`c`/`g`/`t`/`u`/`U`)
/// reset the current k-mer window.
///
/// # Errors
///
/// Returns [`BioError::InvalidParameterError`] if:
/// - `kmer_size` exceeds `seq.len()`
/// - `ds_factor` is `0` or greater than `200`
pub fn frac_min_hash(
    kmer_size: usize,
    ds_factor: u64,
    seq: &[u8],
) -> Result<HashSet<u64>, BioError> {
    if kmer_size > seq.len() {
        return Err(BioError::InvalidParameterError(format!(
            "kmer size {} cannot be longer than sequence len {}.",
            kmer_size,
            seq.len()
        )));
    }

    if ds_factor == 0 || ds_factor > 200 {
        return Err(BioError::InvalidParameterError(format!(
            "downsampling factor {} must be in range 1-200.",
            ds_factor
        )));
    }

    // -- fwd strand
    let mut kmer_forward: u64 = 0;
    let nbits = kmer_size << 1;
    let mask: u64 = (1 << nbits) - 1;

    // -- rev strand
    let mut kmer_reverse: u64 = 0;
    let shift = ((kmer_size - 1) * 2) as u64;

    // -- hashes.
    let mut canonical_hashes: HashSet<u64> = HashSet::with_capacity(seq.len() - kmer_size + 1);
    let mut valid_kmer_index: usize = 0;

    seq.iter().for_each(|nt_char| {
        let nt = NT_LOOKUP[*nt_char as usize] as u64;

        if nt >= 4 {
            valid_kmer_index = 0;
            kmer_forward = 0;
            kmer_reverse = 0;
            return;
        }

        // -- fwd strand
        kmer_forward = (kmer_forward << 2 | nt) & mask;

        // -- rev strand
        let nt_rev = 3 - nt;
        kmer_reverse = kmer_reverse >> 2 | nt_rev << shift;

        if valid_kmer_index >= kmer_size - 1 {
            let canonical = match kmer_forward < kmer_reverse {
                true => kmer_forward,
                false => kmer_reverse,
            };
            if canonical <= u64::MAX / ds_factor {
                canonical_hashes.insert(mm_hash64(canonical));
            }
        }

        valid_kmer_index += 1;
    });

    Ok(canonical_hashes)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;

    #[rstest]
    #[case(b"AAAAAAAA", 3, 1)]
    #[case(b"AAAAAAAC", 3, 2)]
    #[case(b"ATCGATCGATCG", 4, 3)]
    #[case(b"ATCNATCNATCN", 4, 0)]
    fn test_kmerize(
        #[case] seq: &[u8],
        #[case] kmer_size: usize,
        #[case] expected_num_hashes: usize,
    ) {
        let result = frac_min_hash(kmer_size, 1, seq).expect("Failed to run frac_min_hash");
        assert_eq!(result.len(), expected_num_hashes);
    }

    #[rstest]
    #[case(b"AAAAAAAA", b"TTTTTTTT", 3)]
    #[case(b"AAAAAAAAAT", b"ATTTTTTTTT", 3)]
    fn test_kmerize_reverse(#[case] seq1: &[u8], #[case] seq2: &[u8], #[case] kmer_size: usize) {
        let result1 = frac_min_hash(kmer_size, 1, seq1).expect("Failed to run frac_min_hash");
        let result2 = frac_min_hash(kmer_size, 1, seq2).expect("Failed to run frac_min_hash");

        assert_eq!(result1, result2);
    }
}
