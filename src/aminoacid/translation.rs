use super::codon_table::{CodonTable, NT_CODON_MAP};
use super::utils::Frame;

pub fn translate(codon_table_type: CodonTable, frame: &Frame, seq: &[u8]) -> Vec<u8> {
    let start_pos = frame.start_pos();

    if seq.len() < 3 {
        return vec![];
    }

    let codon_table = codon_table_type.table();

    let mut translated: Vec<u8> = Vec::with_capacity(seq.len() / 3);

    for codon in seq[start_pos..].chunks_exact(3) {
        let b1 = NT_CODON_MAP[codon[0] as usize] as usize;
        let b2 = NT_CODON_MAP[codon[1] as usize] as usize;
        let b3 = NT_CODON_MAP[codon[2] as usize] as usize;

        let index = (b1 << 4) | (b2 << 2) | b3;

        let aa = codon_table[index];

        translated.push(aa);

        if aa == b'*' {
            break;
        }
    }

    translated
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use crate::aminoacid::utils::Frame;
    use rstest::*;

    #[rstest]
    // frame 1
    #[case(b"", b"", Frame::First)]
    #[case(b"A", b"", Frame::First)]
    #[case(b"GT", b"", Frame::First)]
    #[case(b"ATG", b"M", Frame::First)]
    #[case(b"ATGTGAAAA", b"M*", Frame::First)] // terminates after `GTG` due to stop codon.
    // frame 2
    #[case(b"A", b"", Frame::Second)]
    #[case(b"ATG", b"", Frame::Second)]
    #[case(b"AATG", b"M", Frame::Second)] // first codon is `ATG`.
    #[case(b"TTGA", b"*", Frame::Second)] // first codon is `TGA`.

    fn test_translate(#[case] seq: &[u8], #[case] expected: &[u8], #[case] frame: Frame) {
        let translated = translate(CodonTable::Standard, &frame, seq);
        assert_eq!(&translated[..], expected);
    }
}
