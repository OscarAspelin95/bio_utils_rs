/// https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=cgencodes#SG1
pub const CODON_STANDARD: &[u8; 64] =
    b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

/// The NCBI codon tables are designed so that the mapping
/// t/T/U	->	0
/// c/C		->	1
/// a/A		->	2
/// g/G		->	3
/// enables extremely efficient codon lookups by converting the codon to a 6-bit representation. Essentially, it is a kmer of length 3.
///
/// With this 2-bit encoding (00, 01, 10, 11) we can create the index for our codon table array as index = (base_1 << 4) | (base_2 << 2) | base_3.
///
/// Example
/// ```
/// use bio_utils_rs::aminoacid::codon_table::{CODON_STANDARD, NT_CODON_MAP};
///
/// let codon = b"ATG";
///
/// let b1 = NT_CODON_MAP[codon[0] as usize] as usize; // =b10
/// let b2 = NT_CODON_MAP[codon[1] as usize] as usize; // =b00
/// let b3 = NT_CODON_MAP[codon[2] as usize] as usize; // =b11
///
/// let index = (b1 << 4) | (b2 << 2) | b3; // b100011
///
/// let aa = CODON_STANDARD[index]; // aminoacid at index 35 -> b'M'
/// assert_eq!(aa, b'M');
/// ```
pub const NT_CODON_MAP: [u8; 256] = {
    let mut map = [0u8; 256];
    map[b'T' as usize] = 0;
    map[b'C' as usize] = 1;
    map[b'A' as usize] = 2;
    map[b'G' as usize] = 3;
    // softmask
    map[b't' as usize] = 0;
    map[b'c' as usize] = 1;
    map[b'a' as usize] = 2;
    map[b'g' as usize] = 3;
    // not sure softmasked `u` is common
    map[b'U' as usize] = 0;

    map
};

/// Types of codon tables.
pub enum CodonTable {
    Standard,
}

impl CodonTable {
    /// Map codon table type to actual array.
    pub fn table(&self) -> &[u8; 64] {
        match self {
            CodonTable::Standard => CODON_STANDARD,
        }
    }
}
