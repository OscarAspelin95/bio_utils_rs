use crate::errors::BioError;

/// Checks whether the run at `[i, j)` qualifies as a homopolymer.
#[inline]
fn valid_homopolymer(
    i: usize,
    j: usize,
    nt: &u8,
    min_hp_len: usize,
    include_softmask: bool,
) -> bool {
    let long_enough = j - i >= min_hp_len;

    if !long_enough {
        return false;
    };

    match include_softmask {
        true => true,
        false => nt.is_ascii_uppercase(),
    }
}

/// Finds all homopolymer runs in a DNA sequence.
///
/// Returns a list of `(start, end, nucleotide, length)` tuples for every
/// run of identical bases at least `min_len` long. Coordinates are
/// zero-based half-open intervals `[start, end)`.
///
/// When `include_softmask` is `false`, lowercase runs are ignored.
///
/// # Errors
///
/// Currently infallible but returns `Result` for forward compatibility.
#[inline]
pub fn find_homopolymers(
    seq: &[u8],
    min_len: usize,
    include_softmask: bool,
) -> Result<Vec<(usize, usize, u8, usize)>, BioError> {
    let mut hps: Vec<(usize, usize, u8, usize)> = Vec::new();

    let seq_len = seq.len();
    if seq_len < min_len {
        return Ok(hps);
    }

    let mut i = 0;
    let mut j = 1;

    while i <= seq_len - min_len {
        while j < seq_len && seq[j] == seq[i] {
            j += 1;
        }

        if valid_homopolymer(i, j, &seq[i], min_len, include_softmask) {
            hps.push((i, j, seq[i], j - i));
        }

        i = j;
        j += 1;
    }

    Ok(hps)
}

/// Inspired by https://github.com/bluenote-1577/myloasm/blob/main/src/kmer_comp.rs#L22
///
/// Applies hard homopolymer compression for a given input nt sequence.
/// e.g., b"AAAATTCGCG" -> vec![b'A', b'T', b'C', b'G', b'C', b'G']
///
/// TODO - try find an empirical value for the capacity for c_seq.
/// Query some genomes from NCBI and calculate before/after to find a suitable empirical value.
#[inline]
pub fn homopolymer_compression(seq: &[u8]) -> Vec<u8> {
    if seq.is_empty() {
        return vec![];
    }

    // capacity can be at most seq.len(), but is in most cases lower.
    let mut c_seq = Vec::with_capacity(seq.len());
    let mut last = seq[0];
    c_seq.push(last);

    for &nt in &seq[1..] {
        if nt != last {
            c_seq.push(nt);
            last = nt;
        }
    }

    c_seq
}

/// Applies soft homopolymer compression for a given input nt sequence and max allowed hp length.
/// e.g., b"AAAATTCGCG" (max_len=2) -> vec![b'A', b'A', b'T', b'T', b'C', b'G', b'C', b'G']
///
/// TODO - try find an empirical value for the capacity for c_seq.
/// Query some genomes from NCBI and calculate before/after to find a suitable empirical value.
fn homopolymer_compression_soft(seq: &[u8], max_len: usize) -> Vec<u8> {
    let mut hp_comp: Vec<u8> = Vec::new();

    if seq.is_empty() {
        return hp_comp;
    }

    let mut current_nt = b'_';
    let mut current_hp_len: usize = 0;

    for nt in seq {
        // Start of new hp
        if nt != &current_nt {
            current_nt = *nt;
            hp_comp.push(*nt);
            current_hp_len = 1;
            continue;
        }

        // Otherwise same nt
        if current_hp_len < max_len {
            hp_comp.push(*nt);
            current_hp_len += 1;
            continue;
        }
    }

    hp_comp
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;

    #[rstest]
    #[case(b"ATCG", 2, false, vec![])]
    #[case(b"AAAA", 4, false, vec![(0, 4, b'A', 4)])]
    #[case(b"aaaa", 4, true, vec![(0, 4, b'a', 4)])]
    #[case(b"aaaa", 4, false, vec![])]
    #[case(b"ATGGGGGCGccccAGT", 4, true, vec![(2, 7, b'G', 5), (9, 13, b'c', 4)])]
    fn test_find_homopolymer(
        #[case] seq: &[u8],
        #[case] min_len: usize,
        #[case] include_softmask: bool,
        #[case] expected: Vec<(usize, usize, u8, usize)>,
    ) {
        let hps = find_homopolymers(seq, min_len, include_softmask).expect("expected valid result");
        assert_eq!(hps, expected);
    }

    #[rstest]
    #[case(b"", b"")]
    #[case(b"T", b"T")]
    #[case(b"AAAA", b"A")]
    #[case(b"CCCGTTT", b"CGT")]
    #[case(b"AAANNGGT", b"ANGT")]
    fn test_homopolymer_compression(#[case] seq: &[u8], #[case] expected: &[u8]) {
        let c_seq = homopolymer_compression(&seq[..]);
        assert_eq!(&c_seq[..], expected);
    }

    #[rstest]
    #[case(b"", 10, b"")]
    #[case(b"ATCG", 1, b"ATCG")]
    #[case(b"AAAAAAAA", 1, b"A")]
    #[case(b"CCCCCCCCCCC", 2, b"CC")]
    #[case(b"GGG", 5, b"GGG")]
    #[case(b"TTT", 0, b"T")]
    #[case(b"ATCAAAGTCCCCCCCCGT", 2, b"ATCAAGTCCGT")]
    #[case(b"AAGGCCTT", 1, b"AGCT")]
    #[case(b"AGCTTTT", 2, b"AGCTT")]
    fn test_hp_comp_soft(#[case] seq: &[u8], #[case] max_len: usize, #[case] expected: &[u8]) {
        let hp_comp = homopolymer_compression_soft(seq, max_len);
        assert_eq!(&hp_comp[..], expected);
    }
}
