use super::traits::Sketcher;
use packed_seq::{PackedSeqVec, SeqVec};
use simd_minimizers::{canonical_minimizers, seq_hash};
use std::collections::HashSet;

pub struct MinimizerSketch {
    pub kmer_size: usize,
    pub window_size: usize,
}

impl Sketcher for MinimizerSketch {
    fn sketch(&self, seq: &[u8]) -> HashSet<u64> {
        let packed_seq = PackedSeqVec::from_ascii(seq);
        let hasher = <seq_hash::NtHasher>::new(self.kmer_size);

        let capacity = seq.len() * 2 / (self.window_size + 1);

        let mut minimizer_positions = Vec::with_capacity(capacity);
        let mut super_kmers = Vec::with_capacity(capacity);

        canonical_minimizers(self.kmer_size, self.window_size)
            .hasher(&hasher)
            .super_kmers(&mut super_kmers)
            .run(packed_seq.as_slice(), &mut minimizer_positions)
            .values_u64()
            .collect()
    }
}
