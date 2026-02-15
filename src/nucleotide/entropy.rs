/// Converts canonical nucleotide counts to relative frequencies.
///
/// Takes a `[A, C, G, T]` count array (as returned by [`nucleotide_counts`](super::nucleotide_counts))
/// and returns a probability vector of the same order.
/// Returns an empty vec if all counts are zero.
#[inline]
pub fn nucleotide_probabilities(canonical: &[usize; 4]) -> Vec<f32> {
    let canonical_count: usize = canonical.iter().sum();

    if canonical_count == 0 {
        return vec![];
    }

    canonical
        .iter()
        .map(|&count| count as f32 / canonical_count as f32)
        .collect()
}

/// Computes Shannon entropy (in bits) from a probability distribution.
///
/// Returns `0.0` for an empty slice. Maximum entropy for four equiprobable
/// nucleotides is `2.0` bits.
#[inline]
pub fn shannon_entropy(probs: &[f32]) -> f32 {
    let shannon: f32 = probs
        .iter()
        .map(|prob| match prob {
            0_f32 => 0.0_f32,
            _ => prob * prob.log2(),
        })
        .sum();

    -shannon
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;

    #[rstest]
    #[case(vec![], 0.0)]
    #[case(vec![0.0_f32, 0.0_f32, 0.0_f32, 0.0_f32], 0.0)]
    #[case(vec![1.0_f32], 0.0)]
    #[case(vec![1.0_f32, 0.0_f32, 0.0_f32, 0.0_f32], 0.0)]
    #[case(vec![0.5_f32, 0.5_f32], 1.0)]
    #[case(vec![0.5_f32, 0.5_f32, 0.0_f32, 0.0_f32], 1.0)]
    fn test_shannon_entropy(#[case] probs: Vec<f32>, #[case] expected: f32) {
        assert_eq!(shannon_entropy(&probs), expected);
    }

    #[test]
    fn test_nucleotide_probabilities_equal() {
        let counts = [10, 10, 10, 10];
        let probs = nucleotide_probabilities(&counts);
        assert_eq!(probs, vec![0.25, 0.25, 0.25, 0.25]);
    }

    #[test]
    fn test_nucleotide_probabilities_single() {
        let counts = [100, 0, 0, 0];
        let probs = nucleotide_probabilities(&counts);
        assert_eq!(probs, vec![1.0, 0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_nucleotide_probabilities_empty() {
        let counts = [0, 0, 0, 0];
        let probs = nucleotide_probabilities(&counts);
        assert!(probs.is_empty());
    }

    #[test]
    fn test_entropy_from_counts() {
        let counts = [25, 25, 25, 25];
        let probs = nucleotide_probabilities(&counts);
        assert_eq!(shannon_entropy(&probs), 2.0);
    }
}
