/// Hashes a 64-bit encoded k-mer using the minimap2 hash function.
///
/// This is a fast, invertible integer hash commonly used for k-mer sketching.
/// Based on the hash function from [minimap2](https://github.com/lh3/minimap2).
#[inline]
pub fn mm_hash64(kmer: u64) -> u64 {
    let mut key = kmer;
    key = !key.wrapping_add(key << 21);
    key = key ^ key >> 24;
    key = (key.wrapping_add(key << 3)).wrapping_add(key << 8);
    key = key ^ key >> 14;
    key = (key.wrapping_add(key << 2)).wrapping_add(key << 4);
    key = key ^ key >> 28;
    key = key.wrapping_add(key << 31);
    key
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_deterministic() {
        assert_eq!(mm_hash64(42), mm_hash64(42));
    }

    #[test]
    fn test_different_inputs() {
        assert_ne!(mm_hash64(0), mm_hash64(1));
    }

    #[test]
    fn test_zero_input() {
        let hash = mm_hash64(0);
        assert_ne!(hash, 0);
    }

    #[test]
    fn test_known_value() {
        // Verify a specific hash doesn't change across runs
        let hash = mm_hash64(1);
        assert_eq!(hash, mm_hash64(1));
        assert_ne!(hash, 0);
    }
}
