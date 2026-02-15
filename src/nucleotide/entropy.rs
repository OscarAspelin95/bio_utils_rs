use rstest::*;
use std::collections::HashMap;

#[inline]
pub fn nucleotide_probabilities(canonical: &HashMap<&u8, usize>) -> Vec<f32> {
    let canonical_count = canonical.values().sum::<usize>();

    let probs: Vec<f32> = canonical
        .values()
        .map(|count| *count as f32 / canonical_count as f32)
        .collect();

    probs
}

#[inline]
pub fn shannon_entropy(probs: &[f32]) -> f32 {
    // Probabilities of each nucleotide.

    let shannon: f32 = probs
        .iter()
        .map(|prob| match prob {
            0_f32 => 0 as f32,
            // This is safe because prob is never negative since (both count and sum_count are usize)
            _ => prob * prob.log2(),
        })
        .sum();

    -shannon
}

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
