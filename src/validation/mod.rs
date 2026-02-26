//! 验证引擎

pub mod accuracy;
pub mod compare;
pub mod cross_format;
pub mod dimensions;
pub mod required_fields;
pub mod roundtrip;
pub mod schema;
pub mod type_compatibility;

pub use accuracy::{
    AccuracyMetrics, AccuracyReport, MetadataAccuracy, EmbeddingAccuracy,
    calculate_matrix_accuracy, calculate_metadata_accuracy, calculate_embedding_accuracy,
    calculate_full_accuracy, pearson_correlation, calculate_sparsity,
    calculate_mse, calculate_exact_match_ratio, calculate_ari, calculate_nmi,
};
pub use cross_format::{
    CrossFormatReport, compare_h5ad_rds, compare_ir_cross_format,
    ConversionDirection, CrossFormatValidationResult, RoundtripValidationResult,
    validate_h5ad_to_rds, validate_rds_to_h5ad, validate_roundtrip,
    compare_h5ad_rds_default, validate_h5ad_to_rds_default,
    validate_rds_to_h5ad_default, validate_roundtrip_default,
    DEFAULT_TOLERANCE,
};
pub use dimensions::{dimension_report, validate_dimensions, validate_expression_matrix_dimensions};
pub use schema::{SchemaValidationResult, validate_cellxgene_schema, CELLXGENE_V5_OBS_REQUIRED, CELLXGENE_V5_VAR_REQUIRED};

use std::collections::HashMap;

/// 验证报告
#[derive(Debug)]
pub struct ValidationReport {
    pub results: HashMap<String, ComparisonResult>,
}

impl ValidationReport {
    /// 创建新的验证报告
    pub fn new() -> Self {
        Self {
            results: HashMap::new(),
        }
    }
    
    /// 添加比较结果
    pub fn add_result(&mut self, component: String, result: ComparisonResult) {
        self.results.insert(component, result);
    }
    
    /// 检查是否所有比较都通过
    pub fn passed(&self) -> bool {
        self.results.values().all(|r| r.passed)
    }
    
    /// 获取失败的组件列表
    pub fn failed_components(&self) -> Vec<&String> {
        self.results
            .iter()
            .filter(|(_, r)| !r.passed)
            .map(|(k, _)| k)
            .collect()
    }
    
    /// 获取最大数值差异
    pub fn max_difference(&self) -> Option<f64> {
        self.results
            .values()
            .filter_map(|r| r.max_difference)
            .fold(None, |acc, diff| match acc {
                None => Some(diff),
                Some(max) => Some(max.max(diff)),
            })
    }
    
    /// 生成摘要报告
    pub fn summary(&self) -> String {
        let total = self.results.len();
        let passed = self.results.values().filter(|r| r.passed).count();
        let failed = total - passed;
        
        let mut summary = format!(
            "Validation Summary: {}/{} components passed",
            passed, total
        );
        
        if failed > 0 {
            summary.push_str(&format!(" ({} failed)", failed));
        }
        
        if let Some(max_diff) = self.max_difference() {
            summary.push_str(&format!("\nMax numerical difference: {:.2e}", max_diff));
        }
        
        summary
    }
    
    /// 生成详细报告
    pub fn detailed_report(&self) -> String {
        let mut report = self.summary();
        report.push_str("\n\nDetailed Results:\n");
        
        // 按组件名排序
        let mut components: Vec<_> = self.results.iter().collect();
        components.sort_by_key(|(k, _)| *k);
        
        for (component, result) in components {
            let status = if result.passed { "✓" } else { "✗" };
            report.push_str(&format!("\n{} {}: {}", status, component, result.message));
            
            if let Some(diff) = result.max_difference {
                report.push_str(&format!(" (max diff: {:.2e})", diff));
            }
        }
        
        report
    }
}

impl Default for ValidationReport {
    fn default() -> Self {
        Self::new()
    }
}

/// 比较结果
#[derive(Debug)]
pub struct ComparisonResult {
    pub passed: bool,
    pub message: String,
    pub max_difference: Option<f64>,
}
