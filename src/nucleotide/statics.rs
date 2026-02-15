use lazy_static::lazy_static;
use rstest::rstest;

pub const PHRED_OFFSET: usize = 33;

/// Inspired by https://github.com/lh3/minimap2/blob/master/sketch.c
///
/// let (A, C, G, T) = (0, 1, 2, 3)
///
/// Lookup table maps:
/// 	(a|A)		=>	A
/// 	(t|T|u|U)	=>	T
/// 	(g|G)		=>	G
/// 	(c|C)		=>	C
/// 	Ambig		=>	Ambig
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

// We should use a better approach here, which is to cap the phred score at e.g., 60.
// This makes a smaller table and avoids calculating non-sensically low error rates.
lazy_static! {
    pub static ref PHRED_TO_ERROR: [f64; 128] = {
        let mut error_lookup: [f64; 128] = [1.0; 128];

        for i in 0..128 {
            if i >= 33 {
                error_lookup[i] = 10_f64.powf(-((i - PHRED_OFFSET) as f64) / 10.0);
            };
        }

        return error_lookup;
    };
}

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
