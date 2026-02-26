//! CELLxGENE Schema 5.0 验证测试
//!
//! 测试 Task 7 的所有子任务：Schema 验证正确性

use arrow::array::{ArrayRef, BooleanArray, StringArray};
use crosscell::ir::*;
use crosscell::validation::schema::{
    validate_cellxgene_schema, SchemaValidationResult, CELLXGENE_V5_OBS_REQUIRED,
    CELLXGENE_V5_VAR_REQUIRED,
};
use std::sync::Arc;

// ============================================================================
// Helpers
// ============================================================================

fn make_obs(n_rows: usize, col_names: &[&str]) -> DataFrame {
    let columns: Vec<String> = col_names.iter().map(|s| s.to_string()).collect();
    let data: Vec<ArrayRef> = col_names
        .iter()
        .map(|_| Arc::new(StringArray::from(vec!["val"; n_rows])) as ArrayRef)
        .collect();
    DataFrame::new(columns, data, n_rows).unwrap()
}

fn make_var(n_rows: usize, col_names: &[&str]) -> DataFrame {
    let columns: Vec<String> = col_names.iter().map(|s| s.to_string()).collect();
    let data: Vec<ArrayRef> = col_names
        .iter()
        .map(|&name| {
            if name == "feature_is_filtered" {
                Arc::new(BooleanArray::from(vec![false; n_rows])) as ArrayRef
            } else {
                Arc::new(StringArray::from(vec!["val"; n_rows])) as ArrayRef
            }
        })
        .collect();
    DataFrame::new(columns, data, n_rows).unwrap()
}

fn make_data(obs_cols: &[&str], var_cols: &[&str]) -> SingleCellData {
    let n_cells = 5;
    let n_genes = 10;
    let obs = make_obs(n_cells, obs_cols);
    let var = make_var(n_genes, var_cols);
    let expr = ExpressionMatrix::Dense(DenseMatrix::zeros(n_cells, n_genes));
    let metadata = DatasetMetadata::new(n_cells, n_genes, "test".to_string());
    SingleCellData::new(expr, obs, var, metadata).unwrap()
}

// ============================================================================
// 7.4 单元测试
// ============================================================================

#[test]
fn test_all_required_fields_pass() {
    let data = make_data(CELLXGENE_V5_OBS_REQUIRED, CELLXGENE_V5_VAR_REQUIRED);
    let result = validate_cellxgene_schema(&data, "cellxgene-v5");
    assert!(
        result.passed,
        "Should pass with all required fields: {:?}",
        result
    );
    assert!(result.missing_obs_fields.is_empty());
    assert!(result.missing_var_fields.is_empty());
}

#[test]
fn test_missing_single_obs_field() {
    // Remove last obs field
    let mut obs_cols: Vec<&str> = CELLXGENE_V5_OBS_REQUIRED.to_vec();
    let removed = obs_cols.pop().unwrap(); // "tissue_type"
    let data = make_data(&obs_cols, CELLXGENE_V5_VAR_REQUIRED);
    let result = validate_cellxgene_schema(&data, "cellxgene-v5");
    assert!(!result.passed);
    assert_eq!(result.missing_obs_fields.len(), 1);
    assert_eq!(result.missing_obs_fields[0], removed);
}

#[test]
fn test_missing_all_obs_fields() {
    let data = make_data(&[], CELLXGENE_V5_VAR_REQUIRED);
    let result = validate_cellxgene_schema(&data, "cellxgene-v5");
    assert!(!result.passed);
    assert_eq!(
        result.missing_obs_fields.len(),
        CELLXGENE_V5_OBS_REQUIRED.len()
    );
}

#[test]
fn test_missing_var_feature_is_filtered() {
    let data = make_data(CELLXGENE_V5_OBS_REQUIRED, &[]);
    let result = validate_cellxgene_schema(&data, "cellxgene-v5");
    assert!(!result.passed);
    assert!(result.missing_obs_fields.is_empty());
    assert_eq!(result.missing_var_fields, vec!["feature_is_filtered"]);
}

#[test]
fn test_both_obs_and_var_missing() {
    let data = make_data(&[], &[]);
    let result = validate_cellxgene_schema(&data, "cellxgene-v5");
    assert!(!result.passed);
    assert_eq!(result.missing_obs_fields.len(), 12);
    assert_eq!(result.missing_var_fields.len(), 1);
}

#[test]
fn test_extra_columns_still_pass() {
    let mut obs_cols: Vec<&str> = CELLXGENE_V5_OBS_REQUIRED.to_vec();
    obs_cols.push("my_custom_column");
    obs_cols.push("another_extra");
    let mut var_cols: Vec<&str> = CELLXGENE_V5_VAR_REQUIRED.to_vec();
    var_cols.push("feature_name");
    let data = make_data(&obs_cols, &var_cols);
    let result = validate_cellxgene_schema(&data, "cellxgene-v5");
    assert!(result.passed, "Extra columns should not cause failure");
}

#[test]
fn test_unsupported_schema_version_fails() {
    let data = make_data(CELLXGENE_V5_OBS_REQUIRED, CELLXGENE_V5_VAR_REQUIRED);
    let result = validate_cellxgene_schema(&data, "cellxgene-v99");
    assert!(!result.passed);
    assert!(result.missing_obs_fields[0].contains("Unsupported"));
}

#[test]
fn test_summary_passed() {
    let result = SchemaValidationResult {
        passed: true,
        missing_obs_fields: vec![],
        missing_var_fields: vec![],
    };
    let summary = result.summary();
    assert!(
        summary.contains("passed"),
        "Summary should say passed: {}",
        summary
    );
}

#[test]
fn test_summary_failed_lists_fields() {
    let result = SchemaValidationResult {
        passed: false,
        missing_obs_fields: vec!["donor_id".to_string(), "tissue_type".to_string()],
        missing_var_fields: vec!["feature_is_filtered".to_string()],
    };
    let summary = result.summary();
    assert!(summary.contains("donor_id"));
    assert!(summary.contains("tissue_type"));
    assert!(summary.contains("feature_is_filtered"));
    assert!(summary.contains("failed"));
}

#[test]
fn test_obs_required_has_12_fields() {
    assert_eq!(
        CELLXGENE_V5_OBS_REQUIRED.len(),
        12,
        "CELLxGENE v5 requires 12 obs fields"
    );
}

#[test]
fn test_var_required_has_1_field() {
    assert_eq!(
        CELLXGENE_V5_VAR_REQUIRED.len(),
        1,
        "CELLxGENE v5 requires 1 var field"
    );
}

// ============================================================================
// 7.3 属性测试 - Property 8: CELLxGENE Schema 验证正确性
// ============================================================================

#[cfg(test)]
mod prop_tests {
    use super::*;
    use proptest::prelude::*;

    // Feature: reviewer-enhancements, Property 8: CELLxGENE Schema 验证正确性
    proptest! {
        /// For any subset of required obs fields, the missing fields should be exactly
        /// the complement of the provided subset.
        #[test]
        fn prop_missing_obs_is_exact_complement(
            // Generate a bitmask for which obs fields to include
            mask in proptest::collection::vec(proptest::bool::ANY, 12..=12)
        ) {
            let included: Vec<&str> = CELLXGENE_V5_OBS_REQUIRED.iter()
                .zip(mask.iter())
                .filter(|(_, &include)| include)
                .map(|(&f, _)| f)
                .collect();
            let expected_missing: Vec<String> = CELLXGENE_V5_OBS_REQUIRED.iter()
                .zip(mask.iter())
                .filter(|(_, &include)| !include)
                .map(|(&f, _)| f.to_string())
                .collect();

            let data = make_data(&included, CELLXGENE_V5_VAR_REQUIRED);
            let result = validate_cellxgene_schema(&data, "cellxgene-v5");

            prop_assert_eq!(&result.missing_obs_fields, &expected_missing,
                "Missing obs fields should be exact complement of provided fields");
            prop_assert_eq!(result.passed, expected_missing.is_empty() && result.missing_var_fields.is_empty());
        }

        /// For any subset of required var fields, the missing fields should be exactly
        /// the complement of the provided subset.
        #[test]
        fn prop_missing_var_is_exact_complement(
            include_feature_is_filtered in proptest::bool::ANY,
        ) {
            let var_cols: Vec<&str> = if include_feature_is_filtered {
                vec!["feature_is_filtered"]
            } else {
                vec![]
            };

            let data = make_data(CELLXGENE_V5_OBS_REQUIRED, &var_cols);
            let result = validate_cellxgene_schema(&data, "cellxgene-v5");

            if include_feature_is_filtered {
                prop_assert!(result.missing_var_fields.is_empty());
                prop_assert!(result.passed);
            } else {
                prop_assert_eq!(result.missing_var_fields, vec!["feature_is_filtered".to_string()]);
                prop_assert!(!result.passed);
            }
        }

        /// Validation with all fields should always pass
        #[test]
        fn prop_all_fields_always_pass(
            n_extra_obs in 0usize..5,
            n_extra_var in 0usize..5,
        ) {
            let mut obs_cols: Vec<&str> = CELLXGENE_V5_OBS_REQUIRED.to_vec();
            for i in 0..n_extra_obs {
                // We can't dynamically create &str with different lifetimes in proptest easily,
                // so we use a fixed set of extra names
                let extras = ["extra_a", "extra_b", "extra_c", "extra_d", "extra_e"];
                if i < extras.len() {
                    obs_cols.push(extras[i]);
                }
            }
            let mut var_cols: Vec<&str> = CELLXGENE_V5_VAR_REQUIRED.to_vec();
            let var_extras = ["feature_name", "feature_biotype", "feature_reference", "feature_length", "extra_var"];
            for i in 0..n_extra_var {
                if i < var_extras.len() {
                    var_cols.push(var_extras[i]);
                }
            }

            let data = make_data(&obs_cols, &var_cols);
            let result = validate_cellxgene_schema(&data, "cellxgene-v5");
            prop_assert!(result.passed, "Should always pass with all required fields + extras");
        }
    }
}
