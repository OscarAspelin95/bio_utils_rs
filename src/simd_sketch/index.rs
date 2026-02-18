use super::traits::Sketcher;
use dashmap::DashMap;
use fixedbitset::FixedBitSet;
use rayon::prelude::*;
use rustc_hash::FxBuildHasher;

pub fn build_reverse_index(
    seqs: &[&[u8]],
    sketcher: &dyn Sketcher,
) -> DashMap<u64, FixedBitSet, FxBuildHasher> {
    let num_seqs = seqs.len();
    let map = DashMap::with_capacity_and_hasher(num_seqs, FxBuildHasher);

    seqs.par_iter().enumerate().for_each(|(i, seq)| {
        let hashes = sketcher.sketch(seq);

        for h in &hashes {
            map.entry(*h)
                .and_modify(|bitset: &mut FixedBitSet| bitset.set(i, true))
                .or_insert_with(|| {
                    let mut bitset = FixedBitSet::with_capacity(num_seqs);
                    bitset.set(i, true);
                    bitset
                });
        }
    });

    map
}
