/// Converts an error probability to a Phred quality score.
///
/// Applies the standard formula: `Q = -10 * log10(error)`.
#[inline]
pub fn error_to_phred(error: f64) -> u8 {
    (-10_f64 * error.log10()) as u8
}

/// Returns the reverse complement of a DNA sequence.
///
/// Handles all IUPAC ambiguity codes. Unrecognized bytes are mapped to `N`.
///
/// # Examples
///
/// ```
/// use bio_utils_rs::nucleotide::reverse_complement;
///
/// assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
/// assert_eq!(reverse_complement(b"AACG"), b"CGTT");
/// ```
#[inline]
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    let mut rc = Vec::with_capacity(seq.len());

    for &nt in seq.iter().rev() {
        rc.push(match nt {
            // Canonical
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            // Ambiguous
            b'R' => b'Y', // AG <-> CT
            b'Y' => b'R', // CT <-> AG
            b'S' => b'S', // GC
            b'W' => b'W', // AT
            b'K' => b'M', // GT <-> AC
            b'M' => b'K', // AC <-> GT
            b'B' => b'V', // CGT <-> ACG
            b'D' => b'H', // AGT <-> ACT
            b'H' => b'D', // ACT <-> AGT
            b'V' => b'B', // ACG <-> CGT
            b'N' => b'N',
            // Unknown nucleotides map to N
            _ => b'N',
        });
    }

    rc
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(b"A", vec![b'T'])]
    #[case(b"ATA", vec![b'T', b'A', b'T'])]
    #[case(b"AACGTT", vec![b'A', b'A', b'C', b'G', b'T', b'T'])]
    fn test_reverse_complement(#[case] seq: &[u8], #[case] expected: Vec<u8>) {
        assert_eq!(reverse_complement(seq), expected);
    }

    #[rstest]
    #[case(b"N", vec![b'N'])]
    #[case(b"TNT", vec![b'A', b'N', b'A'])]
    fn test_reverse_complement_ambig(#[case] seq: &[u8], #[case] expected: Vec<u8>) {
        assert_eq!(reverse_complement(seq), expected);
    }

    #[rstest]
    #[case(b"X", vec![b'N'])]
    #[case(b"AXA", vec![b'T', b'N', b'T'])]
    fn test_reverse_complement_unknown(#[case] seq: &[u8], #[case] expected: Vec<u8>) {
        assert_eq!(reverse_complement(seq), expected);
    }

    #[rstest]
    #[case(0.1, 10)]
    #[case(0.01, 20)]
    #[case(0.001, 30)]
    fn test_error_to_phred(#[case] error: f64, #[case] expected: u8) {
        assert_eq!(error_to_phred(error), expected);
    }
}
