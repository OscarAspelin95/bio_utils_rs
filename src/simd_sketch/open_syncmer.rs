use super::traits::Sketcher;
use packed_seq::{PackedSeqVec, SeqVec};
use simd_minimizers::*;
use std::collections::HashSet;

pub struct OpenSyncmerSketch {
    pub kmer_size: usize,
    pub window_size: usize,
}

impl Sketcher for OpenSyncmerSketch {
    fn sketch(&self, seq: &[u8]) -> HashSet<u64> {
        let packed_seq = PackedSeqVec::from_ascii(seq);
        let mut syncmer_positions = Vec::new();

        canonical_open_syncmers(self.kmer_size, self.window_size)
            .run(packed_seq.as_slice(), &mut syncmer_positions)
            .values_u64()
            .collect()
    }
}
