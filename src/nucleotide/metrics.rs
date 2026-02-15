use super::seq::error_to_phred;
use super::statics::PHRED_TO_ERROR;
use memchr::memchr_iter;

/// Computes the mean error probability and corresponding Phred score for a quality string.
///
/// Each byte in `qual` is treated as a raw Phred+33 encoded quality score and
/// looked up in [`PHRED_TO_ERROR`](super::statics::PHRED_TO_ERROR).
/// Returns `(0.0, 0)` for empty input.
#[inline]
pub fn mean_error_and_phred(qual: &[u8]) -> (f64, u8) {
    if qual.is_empty() {
        return (0.0, 0);
    }

    let error_sum: f64 = qual
        .iter()
        .map(|phred| PHRED_TO_ERROR[*phred as usize])
        .sum::<f64>();

    let error_mean = error_sum / qual.len() as f64;
    (error_mean, error_to_phred(error_mean))
}

/// Returns the truncated mean of a slice of lengths.
///
/// Returns `0` for empty input.
#[inline]
pub fn mean_len(lengths: &[usize]) -> usize {
    if lengths.is_empty() {
        return 0;
    }

    lengths.iter().sum::<usize>() / lengths.len()
}

/// Counts canonical nucleotides in a DNA sequence.
///
/// Returns a tuple of:
/// - `[usize; 4]` — counts for `[A, C, G, T]` (uppercase only).
/// - `usize` — number of softmasked bases (`a`, `c`, `g`, `t`).
/// - `usize` — number of ambiguous/unknown bases (everything else).
///
/// # Examples
///
/// ```
/// use bio_utils_rs::nucleotide::nucleotide_counts;
///
/// let (counts, soft, ambig) = nucleotide_counts(b"AACGttNN");
/// assert_eq!(counts, [2, 1, 1, 0]);
/// assert_eq!(soft, 2);
/// assert_eq!(ambig, 2);
/// ```
#[inline]
pub fn nucleotide_counts(seq: &[u8]) -> ([usize; 4], usize, usize) {
    if seq.is_empty() {
        return ([0; 4], 0, 0);
    }

    let mut canonical = [0usize; 4]; // A=0, C=1, G=2, T=3

    let mut softmasked_count: usize = 0;
    let mut ambiguous_count: usize = 0;

    for &nt in seq {
        match nt {
            b'A' => canonical[0] += 1,
            b'C' => canonical[1] += 1,
            b'G' => canonical[2] += 1,
            b'T' => canonical[3] += 1,
            b'a' | b'c' | b'g' | b't' => softmasked_count += 1,
            _ => ambiguous_count += 1,
        }
    }

    (canonical, softmasked_count, ambiguous_count)
}

/// Computes the GC content of a DNA sequence as a fraction in `[0.0, 1.0]`.
///
/// Counts both uppercase (`G`, `C`) and lowercase (`g`, `c`) bases.
/// Uses SIMD-accelerated byte search via [`memchr`].
/// Returns `0.0` for empty input.
#[inline]
pub fn gc_content(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }

    let gc_count = memchr_iter(b'G', seq).count()
        + memchr_iter(b'C', seq).count()
        + memchr_iter(b'g', seq).count()
        + memchr_iter(b'c', seq).count();

    match gc_count {
        0 => 0.0,
        _ => gc_count as f64 / seq.len() as f64,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(b"", 0.0_f64)]
    #[case(b"A", 0.0_f64)]
    #[case(b"G", 1.0_f64)]
    #[case(b"ATCG", 0.5_f64)]
    #[case(b"AATTC", 1.0_f64 / 5.0_f64)]
    #[case(b"AAAAAAG", 1.0_f64 / 7.0_f64)]
    fn test_gc_content(#[case] seq: &[u8], #[case] expected: f64) {
        assert_eq!(gc_content(seq), expected);
    }

    #[rstest]
    #[case(vec![], 0)]
    #[case(vec![10, 20, 30], 20)]
    fn test_mean_len(#[case] lengths: Vec<usize>, #[case] expected: usize) {
        assert_eq!(mean_len(&lengths), expected);
    }

    #[rstest]
    #[case(b"", [0, 0, 0, 0], 0, 0)]
    #[case(b"aaAA", [2, 0, 0, 0], 2, 0)]
    #[case(b"aaAAttTTccCCggGGNN", [2, 2, 2, 2], 8, 2)]
    fn test_nucleotide_counts(
        #[case] seq: &[u8],
        #[case] expected_counts: [usize; 4],
        #[case] expected_softmasked: usize,
        #[case] expected_ambiguous: usize,
    ) {
        let (counts, softmasked, ambiguous) = nucleotide_counts(seq);

        assert_eq!(counts, expected_counts);
        assert_eq!(softmasked, expected_softmasked);
        assert_eq!(ambiguous, expected_ambiguous);
    }
}
