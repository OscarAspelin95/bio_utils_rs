use bio::pattern_matching::myers::MyersBuilder;
use memchr::memmem;

/// Builds a Myers matcher with all IUPAC ambiguity codes pre-configured.
#[inline]
fn myers_builder(pattern: &[u8]) -> bio::pattern_matching::myers::Myers {
    MyersBuilder::new()
        .ambig(b'N', b"ACGT")
        .ambig(b'R', b"AG")
        .ambig(b'Y', b"CT")
        .ambig(b'S', b"GC")
        .ambig(b'W', b"AT")
        .ambig(b'K', b"GT")
        .ambig(b'M', b"AC")
        .ambig(b'B', b"CGT")
        .ambig(b'D', b"AGT")
        .ambig(b'H', b"ACT")
        .ambig(b'V', b"ACG")
        .build_64(pattern)
}

/// Searches `seq` for approximate matches of `pattern` using the Myers bit-parallel algorithm.
///
/// IUPAC ambiguity codes in the pattern are expanded automatically.
/// Returns a vec of `(end_position, edit_distance)` for every match with
/// at most `max_mismatches` edits. The `end_position` is the inclusive end
/// index of each match in `seq`.
pub fn search_fuzzy(seq: &[u8], pattern: &[u8], max_mismatches: u8) -> Vec<(usize, u8)> {
    let mut myers = myers_builder(pattern);
    myers.find_all_lazy(seq, max_mismatches).collect()
}

/// Searches `seq` for all exact occurrences of `pattern`.
///
/// Uses SIMD-accelerated substring search via [`memchr::memmem`].
/// Returns the starting byte offset of each non-overlapping match.
pub fn search_exact(seq: &[u8], pattern: &[u8]) -> Vec<usize> {
    memmem::find_iter(seq, pattern).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_search_exact_found() {
        let hits = search_exact(b"ACGTACGTACGT", b"ACGT");
        assert_eq!(hits, vec![0, 4, 8]);
    }

    #[test]
    fn test_search_exact_not_found() {
        let hits = search_exact(b"AAAA", b"CG");
        assert!(hits.is_empty());
    }

    #[test]
    fn test_search_exact_single() {
        let hits = search_exact(b"AACGTAA", b"CGT");
        assert_eq!(hits, vec![2]);
    }

    #[test]
    fn test_search_fuzzy_exact_match() {
        let hits = search_fuzzy(b"AACGTAA", b"CGT", 0);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].1, 0); // edit distance 0
    }

    #[test]
    fn test_search_fuzzy_with_mismatches() {
        let hits = search_fuzzy(b"AACCTAA", b"CGT", 2);
        assert!(!hits.is_empty());
        for (_end, dist) in &hits {
            assert!(*dist <= 2);
        }
    }

    #[test]
    fn test_search_fuzzy_no_match() {
        let hits = search_fuzzy(b"AAAAAAA", b"CGT", 0);
        assert!(hits.is_empty());
    }
}
