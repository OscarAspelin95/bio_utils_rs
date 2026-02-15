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
}
