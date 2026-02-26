//! CELLxGENE Schema 验证
//!
//! 检查 SingleCellData 是否符合 CELLxGENE Schema 5.0 的 Curator 必需字段要求。
//! 参考: refs/single-cell-curation/schema/5.0.0/schema.md

use crate::ir::SingleCellData;

/// CELLxGENE Schema 5.0 Curator 必须标注的 obs 字段
pub const CELLXGENE_V5_OBS_REQUIRED: &[&str] = &[
    "assay_ontology_term_id",
    "cell_type_ontology_term_id",
    "development_stage_ontology_term_id",
    "disease_ontology_term_id",
    "donor_id",
    "is_primary_data",
    "organism_ontology_term_id",
    "self_reported_ethnicity_ontology_term_id",
    "sex_ontology_term_id",
    "suspension_type",
    "tissue_ontology_term_id",
    "tissue_type",
];

/// CELLxGENE Schema 5.0 Curator 必须标注的 var 字段
/// feature_is_filtered 是唯一的 curator 必须标注的 var 列
pub const CELLXGENE_V5_VAR_REQUIRED: &[&str] = &["feature_is_filtered"];

/// Schema 验证结果
#[derive(Debug, Clone)]
pub struct SchemaValidationResult {
    /// 是否通过验证
    pub passed: bool,
    /// obs 中缺失的必需字段
    pub missing_obs_fields: Vec<String>,
    /// var 中缺失的必需字段
    pub missing_var_fields: Vec<String>,
}

impl SchemaValidationResult {
    /// 生成摘要报告
    pub fn summary(&self) -> String {
        if self.passed {
            return "✅ CELLxGENE Schema 5.0 validation passed".to_string();
        }
        let mut s = "❌ CELLxGENE Schema 5.0 validation failed\n".to_string();
        if !self.missing_obs_fields.is_empty() {
            s.push_str(&format!(
                "  Missing obs fields ({}):\n",
                self.missing_obs_fields.len()
            ));
            for f in &self.missing_obs_fields {
                s.push_str(&format!("    - {}\n", f));
            }
        }
        if !self.missing_var_fields.is_empty() {
            s.push_str(&format!(
                "  Missing var fields ({}):\n",
                self.missing_var_fields.len()
            ));
            for f in &self.missing_var_fields {
                s.push_str(&format!("    - {}\n", f));
            }
        }
        s
    }
}

/// 验证 SingleCellData 是否符合 CELLxGENE Schema
///
/// 目前支持 schema_version = "cellxgene-v5"
pub fn validate_cellxgene_schema(
    data: &SingleCellData,
    schema_version: &str,
) -> SchemaValidationResult {
    let (obs_required, var_required) = match schema_version {
        "cellxgene-v5" => (CELLXGENE_V5_OBS_REQUIRED, CELLXGENE_V5_VAR_REQUIRED),
        _ => {
            return SchemaValidationResult {
                passed: false,
                missing_obs_fields: vec![format!("Unsupported schema version: {}", schema_version)],
                missing_var_fields: vec![],
            };
        }
    };

    let obs_columns: std::collections::HashSet<&str> = data
        .cell_metadata
        .columns
        .iter()
        .map(|s| s.as_str())
        .collect();

    let var_columns: std::collections::HashSet<&str> = data
        .gene_metadata
        .columns
        .iter()
        .map(|s| s.as_str())
        .collect();

    let missing_obs: Vec<String> = obs_required
        .iter()
        .filter(|&&f| !obs_columns.contains(f))
        .map(|&f| f.to_string())
        .collect();

    let missing_var: Vec<String> = var_required
        .iter()
        .filter(|&&f| !var_columns.contains(f))
        .map(|&f| f.to_string())
        .collect();

    let passed = missing_obs.is_empty() && missing_var.is_empty();

    SchemaValidationResult {
        passed,
        missing_obs_fields: missing_obs,
        missing_var_fields: missing_var,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::*;
    use arrow::array::{BooleanArray, StringArray};
    use std::sync::Arc;

    /// Helper: create a DataFrame with given column names (all string columns)
    fn make_df(n_rows: usize, col_names: &[&str]) -> DataFrame {
        let columns: Vec<String> = col_names.iter().map(|s| s.to_string()).collect();
        let data: Vec<arrow::array::ArrayRef> = col_names
            .iter()
            .map(|_| Arc::new(StringArray::from(vec!["x"; n_rows])) as arrow::array::ArrayRef)
            .collect();
        DataFrame::new(columns, data, n_rows).unwrap()
    }

    /// Helper: create a DataFrame with given column names, using bool for feature_is_filtered
    fn make_var_df(n_rows: usize, col_names: &[&str]) -> DataFrame {
        let columns: Vec<String> = col_names.iter().map(|s| s.to_string()).collect();
        let data: Vec<arrow::array::ArrayRef> = col_names
            .iter()
            .map(|&name| {
                if name == "feature_is_filtered" {
                    Arc::new(BooleanArray::from(vec![false; n_rows])) as arrow::array::ArrayRef
                } else {
                    Arc::new(StringArray::from(vec!["x"; n_rows])) as arrow::array::ArrayRef
                }
            })
            .collect();
        DataFrame::new(columns, data, n_rows).unwrap()
    }

    #[test]
    fn test_all_fields_present_passes() {
        let obs = make_df(5, CELLXGENE_V5_OBS_REQUIRED);
        let var = make_var_df(10, CELLXGENE_V5_VAR_REQUIRED);
        let expr = ExpressionMatrix::Dense(DenseMatrix::zeros(5, 10));
        let metadata = DatasetMetadata::new(5, 10, "test".to_string());
        let data = SingleCellData::new(expr, obs, var, metadata).unwrap();

        let result = validate_cellxgene_schema(&data, "cellxgene-v5");
        assert!(
            result.passed,
            "Should pass when all fields present: {:?}",
            result
        );
        assert!(result.missing_obs_fields.is_empty());
        assert!(result.missing_var_fields.is_empty());
    }

    #[test]
    fn test_missing_obs_fields() {
        // Only provide first 3 obs fields
        let partial_obs = &CELLXGENE_V5_OBS_REQUIRED[..3];
        let obs = make_df(5, partial_obs);
        let var = make_var_df(10, CELLXGENE_V5_VAR_REQUIRED);
        let expr = ExpressionMatrix::Dense(DenseMatrix::zeros(5, 10));
        let metadata = DatasetMetadata::new(5, 10, "test".to_string());
        let data = SingleCellData::new(expr, obs, var, metadata).unwrap();

        let result = validate_cellxgene_schema(&data, "cellxgene-v5");
        assert!(!result.passed);
        assert_eq!(
            result.missing_obs_fields.len(),
            CELLXGENE_V5_OBS_REQUIRED.len() - 3
        );
        assert!(result.missing_var_fields.is_empty());
    }

    #[test]
    fn test_missing_var_fields() {
        let obs = make_df(5, CELLXGENE_V5_OBS_REQUIRED);
        let var = make_var_df(10, &[]); // no var columns
        let expr = ExpressionMatrix::Dense(DenseMatrix::zeros(5, 10));
        let metadata = DatasetMetadata::new(5, 10, "test".to_string());
        let data = SingleCellData::new(expr, obs, var, metadata).unwrap();

        let result = validate_cellxgene_schema(&data, "cellxgene-v5");
        assert!(!result.passed);
        assert!(result.missing_obs_fields.is_empty());
        assert_eq!(result.missing_var_fields, vec!["feature_is_filtered"]);
    }

    #[test]
    fn test_empty_dataframes() {
        let obs = DataFrame::empty(0);
        let var = DataFrame::empty(0);
        let expr = ExpressionMatrix::Dense(DenseMatrix::zeros(0, 0));
        let metadata = DatasetMetadata::new(0, 0, "test".to_string());
        let data = SingleCellData::new(expr, obs, var, metadata).unwrap();

        let result = validate_cellxgene_schema(&data, "cellxgene-v5");
        assert!(!result.passed);
        assert_eq!(result.missing_obs_fields.len(), 12);
        assert_eq!(result.missing_var_fields.len(), 1);
    }

    #[test]
    fn test_unsupported_schema_version() {
        let obs = DataFrame::empty(0);
        let var = DataFrame::empty(0);
        let expr = ExpressionMatrix::Dense(DenseMatrix::zeros(0, 0));
        let metadata = DatasetMetadata::new(0, 0, "test".to_string());
        let data = SingleCellData::new(expr, obs, var, metadata).unwrap();

        let result = validate_cellxgene_schema(&data, "cellxgene-v99");
        assert!(!result.passed);
        assert!(result.missing_obs_fields[0].contains("Unsupported"));
    }

    #[test]
    fn test_extra_columns_dont_affect_validation() {
        // All required + extra columns
        let mut cols: Vec<&str> = CELLXGENE_V5_OBS_REQUIRED.to_vec();
        cols.push("extra_column_1");
        cols.push("extra_column_2");
        let obs = make_df(5, &cols);
        let var = make_var_df(10, &["feature_is_filtered", "feature_name"]);
        let expr = ExpressionMatrix::Dense(DenseMatrix::zeros(5, 10));
        let metadata = DatasetMetadata::new(5, 10, "test".to_string());
        let data = SingleCellData::new(expr, obs, var, metadata).unwrap();

        let result = validate_cellxgene_schema(&data, "cellxgene-v5");
        assert!(result.passed, "Extra columns should not cause failure");
    }

    #[test]
    fn test_summary_output() {
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
}
