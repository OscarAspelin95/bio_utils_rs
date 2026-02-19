use std::sync::LazyLock;

/// ASCII offset for Phred+33 quality encoding (Sanger/Illumina 1.8+).
pub const PHRED_OFFSET: usize = 33;

/// Maximum Phred index stored in [`PHRED_TO_ERROR`] (Phred 60 = index 93).
const MAX_PHRED_INDEX: usize = 93;

/// 2-bit nucleotide encoding table indexed by ASCII byte value.
///
/// Encodes `A`/`a` = 0, `C`/`c` = 1, `G`/`g` = 2, `T`/`t`/`U`/`u` = 3.
/// All other bytes (ambiguous/invalid) map to 4.
///
/// Inspired by [minimap2](https://github.com/lh3/minimap2/blob/master/sketch.c).
pub const NT_LOOKUP: [u8; 256] = [
    0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
];

/// Phred score to error probability lookup table.
///
/// Indexed by raw quality byte (Phred+33 encoded). Indices below
/// [`PHRED_OFFSET`] default to `1.0`. The table is capped at Phred 60
/// (index 93) since higher scores represent negligible error rates.
pub static PHRED_TO_ERROR: LazyLock<[f64; MAX_PHRED_INDEX + 1]> = LazyLock::new(|| {
    let mut error_lookup = [1.0; MAX_PHRED_INDEX + 1];

    for (i, entry) in error_lookup.iter_mut().enumerate().skip(PHRED_OFFSET) {
        *entry = 10_f64.powf(-((i - PHRED_OFFSET) as f64) / 10.0);
    }

    error_lookup
});

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    // A
    #[case(0, 0)]
    #[case(b'a', 0)]
    #[case(b'A', 0)]
    // C
    #[case(1, 1)]
    #[case(b'c', 1)]
    #[case(b'C', 1)]
    // G
    #[case(2, 2)]
    #[case(b'g', 2)]
    #[case(b'G', 2)]
    // T
    #[case(3, 3)]
    #[case(b't', 3)]
    #[case(b'T', 3)]
    // U
    #[case(b'u', 3)]
    #[case(b'U', 3)]
    // Ambig
    #[case(4, 4)]
    #[case(b'N', 4)]
    #[case(b'Y', 4)]
    #[case(b'R', 4)]
    fn test_lookup_table(#[case] nt: u8, #[case] expected: u8) {
        assert_eq!(NT_LOOKUP[nt as usize], expected);
    }

    #[rstest]
    #[case(30 + PHRED_OFFSET, 0.001)]
    #[case(40 + PHRED_OFFSET, 0.0001)]
    #[case(50 + PHRED_OFFSET, 0.00001)]
    fn test_phred_to_error(#[case] phred: usize, #[case] expected: f64) {
        assert_eq!(PHRED_TO_ERROR[phred], expected);
    }

    #[test]
    fn test_phred_table_size() {
        assert_eq!(PHRED_TO_ERROR.len(), MAX_PHRED_INDEX + 1);
    }
}
