//! 审计指标增强测试
//!
//! 测试 Task 5 的所有子任务：MSE、精确匹配率、ARI、NMI

use crosscell::validation::accuracy::{
    AccuracyMetrics, calculate_mse, calculate_exact_match_ratio,
    calculate_ari, calculate_nmi,
};

// ============================================================================
// 5.9 单元测试
// ============================================================================

#[test]
fn test_mse_identical() {
    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let mse = calculate_mse(&x, &x);
    assert!((mse - 0.0).abs() < 1e-15, "MSE of identical vectors should be 0");
}

#[test]
fn test_mse_known_value() {
    // x = [1, 2, 3], y = [2, 3, 4] => diffs = [1, 1, 1] => MSE = 1.0
    let x = vec![1.0, 2.0, 3.0];
    let y = vec![2.0, 3.0, 4.0];
    let mse = calculate_mse(&x, &y);
    assert!((mse - 1.0).abs() < 1e-15, "MSE should be 1.0, got {}", mse);
}

#[test]
fn test_mse_known_value_2() {
    // x = [0, 0], y = [3, 4] => diffs^2 = [9, 16] => MSE = 12.5
    let x = vec![0.0, 0.0];
    let y = vec![3.0, 4.0];
    let mse = calculate_mse(&x, &y);
    assert!((mse - 12.5).abs() < 1e-15, "MSE should be 12.5, got {}", mse);
}

#[test]
fn test_exact_match_ratio_identical() {
    let x = vec![1.0, 2.0, 3.0];
    let ratio = calculate_exact_match_ratio(&x, &x);
    assert!((ratio - 1.0).abs() < 1e-15, "Exact match ratio of identical vectors should be 1.0");
}

#[test]
fn test_exact_match_ratio_none_match() {
    let x = vec![1.0, 2.0, 3.0];
    let y = vec![4.0, 5.0, 6.0];
    let ratio = calculate_exact_match_ratio(&x, &y);
    assert!((ratio - 0.0).abs() < 1e-15, "No matches should give 0.0");
}

#[test]
fn test_exact_match_ratio_partial() {
    let x = vec![1.0, 2.0, 3.0, 4.0];
    let y = vec![1.0, 2.0, 9.0, 9.0];
    let ratio = calculate_exact_match_ratio(&x, &y);
    assert!((ratio - 0.5).abs() < 1e-15, "2/4 matches should give 0.5, got {}", ratio);
}

#[test]
fn test_ari_identical_labels() {
    let labels = vec!["A".to_string(), "B".to_string(), "A".to_string(), "B".to_string(), "C".to_string()];
    let ari = calculate_ari(&labels, &labels);
    assert!((ari - 1.0).abs() < 1e-10, "ARI of identical labels should be 1.0, got {}", ari);
}

#[test]
fn test_ari_completely_different() {
    // Two completely different clusterings
    let a = vec!["A".to_string(), "A".to_string(), "B".to_string(), "B".to_string()];
    let b = vec!["X".to_string(), "Y".to_string(), "X".to_string(), "Y".to_string()];
    let ari = calculate_ari(&a, &b);
    // ARI should be near -0.5 for this anti-correlated case (or at least < 0.5)
    assert!(ari < 0.5, "ARI of anti-correlated labels should be low, got {}", ari);
}

#[test]
fn test_nmi_identical_labels() {
    let labels = vec!["A".to_string(), "B".to_string(), "A".to_string(), "B".to_string(), "C".to_string()];
    let nmi = calculate_nmi(&labels, &labels);
    assert!((nmi - 1.0).abs() < 1e-10, "NMI of identical labels should be 1.0, got {}", nmi);
}

#[test]
fn test_nmi_single_cluster() {
    // All same label => H(A) = 0, H(B) = 0 => NMI = 1.0 (degenerate)
    let a = vec!["A".to_string(), "A".to_string(), "A".to_string()];
    let b = vec!["A".to_string(), "A".to_string(), "A".to_string()];
    let nmi = calculate_nmi(&a, &b);
    assert!((nmi - 1.0).abs() < 1e-10, "NMI of single cluster should be 1.0, got {}", nmi);
}

#[test]
fn test_ari_none_when_no_cluster_column() {
    // AccuracyMetrics::perfect() should have ari=None, nmi=None
    let metrics = AccuracyMetrics::perfect();
    assert!(metrics.ari.is_none(), "ARI should be None by default");
    assert!(metrics.nmi.is_none(), "NMI should be None by default");
}

#[test]
fn test_mse_empty_vectors() {
    let mse = calculate_mse(&[], &[]);
    assert!((mse - 0.0).abs() < 1e-15, "MSE of empty vectors should be 0");
}

#[test]
fn test_summary_contains_new_metrics() {
    let metrics = AccuracyMetrics::perfect();
    let summary = metrics.summary();
    assert!(summary.contains("MSE"), "Summary should contain MSE, got: {}", summary);
    assert!(summary.contains("Exact Match"), "Summary should contain Exact Match, got: {}", summary);
}

// ============================================================================
// 5.6 属性测试 - Property 5: 指标数学性质
// ============================================================================

#[cfg(test)]
mod prop_tests {
    use super::*;
    use proptest::prelude::*;

    // Feature: reviewer-enhancements, Property 5: 指标数学性质
    proptest! {
        #[test]
        fn prop_mse_self_is_zero(values in proptest::collection::vec(-1000.0f64..1000.0, 1..200)) {
            let mse = calculate_mse(&values, &values);
            prop_assert!((mse - 0.0).abs() < 1e-10, "MSE(x, x) should be 0, got {}", mse);
        }

        #[test]
        fn prop_exact_match_self_is_one(values in proptest::collection::vec(-1000.0f64..1000.0, 1..200)) {
            let ratio = calculate_exact_match_ratio(&values, &values);
            prop_assert!((ratio - 1.0).abs() < 1e-10, "ExactMatch(x, x) should be 1.0, got {}", ratio);
        }

        #[test]
        fn prop_mse_non_negative(
            x in proptest::collection::vec(-100.0f64..100.0, 1..100),
            y in proptest::collection::vec(-100.0f64..100.0, 1..100),
        ) {
            let len = x.len().min(y.len());
            let mse = calculate_mse(&x[..len], &y[..len]);
            prop_assert!(mse >= 0.0, "MSE should be non-negative, got {}", mse);
        }

        #[test]
        fn prop_exact_match_ratio_in_range(
            x in proptest::collection::vec(-100.0f64..100.0, 1..100),
            y in proptest::collection::vec(-100.0f64..100.0, 1..100),
        ) {
            let len = x.len().min(y.len());
            let ratio = calculate_exact_match_ratio(&x[..len], &y[..len]);
            prop_assert!(ratio >= 0.0 && ratio <= 1.0, "Ratio should be in [0, 1], got {}", ratio);
        }
    }

    // Feature: reviewer-enhancements, Property 6: ARI 数学性质
    proptest! {
        #[test]
        fn prop_ari_self_is_one(
            labels in proptest::collection::vec("[A-E]", 3..50)
        ) {
            let ari = calculate_ari(&labels, &labels);
            prop_assert!((ari - 1.0).abs() < 1e-10, "ARI(x, x) should be 1.0, got {}", ari);
        }

        #[test]
        fn prop_ari_symmetric(
            a in proptest::collection::vec("[A-D]", 5..30),
            b in proptest::collection::vec("[X-Z]", 5..30),
        ) {
            let len = a.len().min(b.len());
            let a = &a[..len];
            let b = &b[..len];
            let ari_ab = calculate_ari(a, b);
            let ari_ba = calculate_ari(b, a);
            prop_assert!((ari_ab - ari_ba).abs() < 1e-10, "ARI should be symmetric: {} vs {}", ari_ab, ari_ba);
        }
    }

    // Feature: reviewer-enhancements, Property 7: Summary 包含所有指标名称
    proptest! {
        #[test]
        fn prop_summary_contains_metric_names(
            corr in 0.0f64..1.0,
            mse_val in 0.0f64..100.0,
        ) {
            let metrics = AccuracyMetrics {
                correlation: corr,
                mean_relative_error: 0.0,
                mean_absolute_error: 0.0,
                max_absolute_error: 0.0,
                sparsity_preserved: 1.0,
                metadata_match: true,
                categorical_match: true,
                embedding_error: 0.0,
                nnz_match: true,
                original_nnz: 0,
                converted_nnz: 0,
                mse: mse_val,
                exact_match_ratio: 0.5,
                ari: None,
                nmi: None,
            };
            let summary = metrics.summary();
            prop_assert!(summary.contains("MSE"), "Summary must contain MSE");
            prop_assert!(summary.contains("Exact Match"), "Summary must contain Exact Match");
        }
    }
}
