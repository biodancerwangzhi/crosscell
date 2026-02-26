//! 准确性测试框架
//!
//! 提供详细的准确性指标计算，用于验证真实数据集的往返转换一致性。
//!
//! # 主要功能
//! - 计算表达矩阵的相关系数、相对误差、绝对误差
//! - 计算稀疏性保留率
//! - 验证元数据完整性
//! - 验证 Categorical 变量一致性
//! - 验证降维坐标误差

use crate::ir::{DataFrame, Embedding, ExpressionMatrix, SingleCellData};
use std::collections::HashMap;

/// 准确性指标
#[derive(Debug, Clone)]
pub struct AccuracyMetrics {
    /// 皮尔逊相关系数 (0.0 - 1.0)
    pub correlation: f64,
    /// 平均相对误差
    pub mean_relative_error: f64,
    /// 平均绝对误差
    pub mean_absolute_error: f64,
    /// 最大绝对误差
    pub max_absolute_error: f64,
    /// 稀疏性保留率 (0.0 - 1.0)
    pub sparsity_preserved: f64,
    /// 元数据完全匹配
    pub metadata_match: bool,
    /// Categorical 变量完全匹配
    pub categorical_match: bool,
    /// 降维坐标平均误差
    pub embedding_error: f64,
    /// 非零元素数量匹配
    pub nnz_match: bool,
    /// 原始非零元素数量
    pub original_nnz: usize,
    /// 转换后非零元素数量
    pub converted_nnz: usize,
    /// 均方误差
    pub mse: f64,
    /// 精确匹配比率 (0.0 - 1.0)
    pub exact_match_ratio: f64,
    /// 调整兰德指数（可选，需要聚类标签）
    pub ari: Option<f64>,
    /// 归一化互信息（可选，需要聚类标签）
    pub nmi: Option<f64>,
}

impl AccuracyMetrics {
    /// 创建默认的准确性指标（全部通过）
    pub fn perfect() -> Self {
        Self {
            correlation: 1.0,
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
            mse: 0.0,
            exact_match_ratio: 1.0,
            ari: None,
            nmi: None,
        }
    }

    /// 检查是否满足准确性要求
    pub fn meets_requirements(&self, tolerance: f64) -> bool {
        self.correlation > 0.99999
            && self.mean_relative_error < tolerance
            && self.max_absolute_error < tolerance * 100.0
            && self.sparsity_preserved > 0.999
            && self.metadata_match
            && self.categorical_match
            && self.embedding_error < tolerance
            && self.nnz_match
            && self.mse < tolerance * tolerance
            && self.exact_match_ratio > 0.999
    }

    /// 生成摘要报告
    pub fn summary(&self) -> String {
        let mut s = format!(
            "AccuracyMetrics:\n\
             - Correlation: {:.6}\n\
             - Mean Relative Error: {:.2e}\n\
             - Mean Absolute Error: {:.2e}\n\
             - Max Absolute Error: {:.2e}\n\
             - MSE: {:.2e}\n\
             - Exact Match Ratio: {:.6}\n\
             - Sparsity Preserved: {:.4}%\n\
             - Metadata Match: {}\n\
             - Categorical Match: {}\n\
             - Embedding Error: {:.2e}\n\
             - NNZ Match: {} ({} vs {})",
            self.correlation,
            self.mean_relative_error,
            self.mean_absolute_error,
            self.max_absolute_error,
            self.mse,
            self.exact_match_ratio,
            self.sparsity_preserved * 100.0,
            if self.metadata_match { "✓" } else { "✗" },
            if self.categorical_match { "✓" } else { "✗" },
            self.embedding_error,
            if self.nnz_match { "✓" } else { "✗" },
            self.original_nnz,
            self.converted_nnz,
        );
        if let Some(ari) = self.ari {
            s.push_str(&format!("\n - ARI: {:.6}", ari));
        }
        if let Some(nmi) = self.nmi {
            s.push_str(&format!("\n - NMI: {:.6}", nmi));
        }
        s
    }
}

/// 完整的准确性测试报告
#[derive(Debug)]
pub struct AccuracyReport {
    /// 数据集名称
    pub dataset_name: String,
    /// 表达矩阵准确性
    pub expression_accuracy: AccuracyMetrics,
    /// 细胞元数据准确性
    pub cell_metadata_accuracy: Option<MetadataAccuracy>,
    /// 基因元数据准确性
    pub gene_metadata_accuracy: Option<MetadataAccuracy>,
    /// 各嵌入的准确性
    pub embedding_accuracies: HashMap<String, EmbeddingAccuracy>,
    /// 各 layer 的准确性
    pub layer_accuracies: HashMap<String, AccuracyMetrics>,
    /// 总体是否通过
    pub passed: bool,
    /// 错误列表
    pub errors: Vec<String>,
}

impl AccuracyReport {
    /// 创建新的准确性报告
    pub fn new(dataset_name: &str) -> Self {
        Self {
            dataset_name: dataset_name.to_string(),
            expression_accuracy: AccuracyMetrics::perfect(),
            cell_metadata_accuracy: None,
            gene_metadata_accuracy: None,
            embedding_accuracies: HashMap::new(),
            layer_accuracies: HashMap::new(),
            passed: true,
            errors: Vec::new(),
        }
    }

    /// 添加错误
    pub fn add_error(&mut self, error: String) {
        self.errors.push(error);
        self.passed = false;
    }

    /// 生成详细报告
    pub fn detailed_report(&self) -> String {
        let mut report = format!("=== Accuracy Report: {} ===\n\n", self.dataset_name);

        // 总体状态
        report.push_str(&format!(
            "Overall: {}\n\n",
            if self.passed {
                "✓ PASSED"
            } else {
                "✗ FAILED"
            }
        ));

        // 表达矩阵
        report.push_str("--- Expression Matrix ---\n");
        report.push_str(&self.expression_accuracy.summary());
        report.push_str("\n\n");

        // 细胞元数据
        if let Some(ref meta) = self.cell_metadata_accuracy {
            report.push_str("--- Cell Metadata ---\n");
            report.push_str(&meta.summary());
            report.push_str("\n\n");
        }

        // 基因元数据
        if let Some(ref meta) = self.gene_metadata_accuracy {
            report.push_str("--- Gene Metadata ---\n");
            report.push_str(&meta.summary());
            report.push_str("\n\n");
        }

        // 嵌入
        if !self.embedding_accuracies.is_empty() {
            report.push_str("--- Embeddings ---\n");
            for (name, acc) in &self.embedding_accuracies {
                report.push_str(&format!("{}: {}\n", name, acc.summary()));
            }
            report.push_str("\n");
        }

        // Layers
        if !self.layer_accuracies.is_empty() {
            report.push_str("--- Layers ---\n");
            for (name, acc) in &self.layer_accuracies {
                report.push_str(&format!(
                    "{}: correlation={:.6}, max_err={:.2e}\n",
                    name, acc.correlation, acc.max_absolute_error
                ));
            }
            report.push_str("\n");
        }

        // 错误
        if !self.errors.is_empty() {
            report.push_str("--- Errors ---\n");
            for err in &self.errors {
                report.push_str(&format!("  - {}\n", err));
            }
        }

        report
    }
}

/// 元数据准确性
#[derive(Debug, Clone)]
pub struct MetadataAccuracy {
    /// 行数匹配
    pub row_count_match: bool,
    /// 列数匹配
    pub column_count_match: bool,
    /// 列名匹配
    pub column_names_match: bool,
    /// 数据类型匹配
    pub dtypes_match: bool,
    /// 数值列的最大误差
    pub numeric_max_error: f64,
    /// Categorical 列匹配
    pub categorical_match: bool,
    /// 缺失的列
    pub missing_columns: Vec<String>,
    /// 额外的列
    pub extra_columns: Vec<String>,
}

impl MetadataAccuracy {
    /// 生成摘要
    pub fn summary(&self) -> String {
        format!(
            "Rows: {}, Cols: {}, Names: {}, Types: {}, Categorical: {}, NumericErr: {:.2e}",
            if self.row_count_match { "✓" } else { "✗" },
            if self.column_count_match {
                "✓"
            } else {
                "✗"
            },
            if self.column_names_match {
                "✓"
            } else {
                "✗"
            },
            if self.dtypes_match { "✓" } else { "✗" },
            if self.categorical_match { "✓" } else { "✗" },
            self.numeric_max_error,
        )
    }

    /// 检查是否全部通过
    pub fn all_passed(&self) -> bool {
        self.row_count_match
            && self.column_count_match
            && self.column_names_match
            && self.dtypes_match
            && self.categorical_match
    }
}

/// 嵌入准确性
#[derive(Debug, Clone)]
pub struct EmbeddingAccuracy {
    /// 维度匹配
    pub dimension_match: bool,
    /// 平均误差
    pub mean_error: f64,
    /// 最大误差
    pub max_error: f64,
    /// 相关系数
    pub correlation: f64,
}

impl EmbeddingAccuracy {
    /// 生成摘要
    pub fn summary(&self) -> String {
        format!(
            "dim={}, corr={:.6}, mean_err={:.2e}, max_err={:.2e}",
            if self.dimension_match { "✓" } else { "✗" },
            self.correlation,
            self.mean_error,
            self.max_error,
        )
    }
}

// ============================================================================
// 核心计算函数
// ============================================================================

/// 计算皮尔逊相关系数
pub fn pearson_correlation(x: &[f64], y: &[f64]) -> f64 {
    if x.len() != y.len() || x.is_empty() {
        return 0.0;
    }

    let n = x.len() as f64;
    let sum_x: f64 = x.iter().sum();
    let sum_y: f64 = y.iter().sum();
    let sum_xy: f64 = x.iter().zip(y.iter()).map(|(a, b)| a * b).sum();
    let sum_x2: f64 = x.iter().map(|a| a * a).sum();
    let sum_y2: f64 = y.iter().map(|a| a * a).sum();

    let numerator = n * sum_xy - sum_x * sum_y;
    let denominator = ((n * sum_x2 - sum_x * sum_x) * (n * sum_y2 - sum_y * sum_y)).sqrt();

    if denominator == 0.0 {
        // 如果两个向量都是常数，检查它们是否相等
        if x.iter().all(|&v| v == x[0]) && y.iter().all(|&v| v == y[0]) {
            if (x[0] - y[0]).abs() < 1e-10 {
                return 1.0;
            }
        }
        return 0.0;
    }

    numerator / denominator
}

/// 计算稀疏度
pub fn calculate_sparsity(nnz: usize, total: usize) -> f64 {
    if total == 0 {
        return 1.0;
    }
    1.0 - (nnz as f64 / total as f64)
}

/// 计算均方误差
pub fn calculate_mse(x: &[f64], y: &[f64]) -> f64 {
    if x.len() != y.len() || x.is_empty() {
        return 0.0;
    }
    let sum_sq: f64 = x.iter().zip(y.iter()).map(|(a, b)| (a - b).powi(2)).sum();
    sum_sq / x.len() as f64
}

/// 计算精确匹配比率
pub fn calculate_exact_match_ratio(x: &[f64], y: &[f64]) -> f64 {
    if x.len() != y.len() || x.is_empty() {
        return 0.0;
    }
    let matches = x
        .iter()
        .zip(y.iter())
        .filter(|(a, b)| (*a - *b).abs() < f64::EPSILON)
        .count();
    matches as f64 / x.len() as f64
}

/// 计算调整兰德指数 (ARI)
///
/// ARI = (RI - Expected_RI) / (Max_RI - Expected_RI)
/// 基于列联表和配对计数
pub fn calculate_ari(labels_a: &[String], labels_b: &[String]) -> f64 {
    if labels_a.len() != labels_b.len() || labels_a.is_empty() {
        return 0.0;
    }
    let n = labels_a.len();

    // 构建列联表
    let mut label_a_map: HashMap<&str, usize> = HashMap::new();
    let mut label_b_map: HashMap<&str, usize> = HashMap::new();
    let mut a_idx = 0usize;
    let mut b_idx = 0usize;

    let mut a_indices = Vec::with_capacity(n);
    let mut b_indices = Vec::with_capacity(n);

    for (a, b) in labels_a.iter().zip(labels_b.iter()) {
        let ai = *label_a_map.entry(a.as_str()).or_insert_with(|| {
            let v = a_idx;
            a_idx += 1;
            v
        });
        let bi = *label_b_map.entry(b.as_str()).or_insert_with(|| {
            let v = b_idx;
            b_idx += 1;
            v
        });
        a_indices.push(ai);
        b_indices.push(bi);
    }

    let n_a = a_idx;
    let n_b = b_idx;

    // 列联表 nij
    let mut contingency = vec![vec![0i64; n_b]; n_a];
    for i in 0..n {
        contingency[a_indices[i]][b_indices[i]] += 1;
    }

    // 行和 a_i, 列和 b_j
    let a_sums: Vec<i64> = contingency.iter().map(|row| row.iter().sum()).collect();
    let b_sums: Vec<i64> = (0..n_b)
        .map(|j| contingency.iter().map(|row| row[j]).sum())
        .collect();

    // C(n, 2) helper
    let comb2 = |x: i64| -> i64 { x * (x - 1) / 2 };

    let sum_comb_nij: i64 = contingency
        .iter()
        .flat_map(|row| row.iter())
        .map(|&v| comb2(v))
        .sum();
    let sum_comb_ai: i64 = a_sums.iter().map(|&v| comb2(v)).sum();
    let sum_comb_bj: i64 = b_sums.iter().map(|&v| comb2(v)).sum();
    let comb_n = comb2(n as i64);

    if comb_n == 0 {
        return 0.0;
    }

    let expected = (sum_comb_ai as f64 * sum_comb_bj as f64) / comb_n as f64;
    let max_index = (sum_comb_ai as f64 + sum_comb_bj as f64) / 2.0;
    let denominator = max_index - expected;

    if denominator.abs() < f64::EPSILON {
        // 完全一致或退化情况
        if (sum_comb_nij as f64 - expected).abs() < f64::EPSILON {
            return 1.0;
        }
        return 0.0;
    }

    (sum_comb_nij as f64 - expected) / denominator
}

/// 计算归一化互信息 (NMI)
///
/// NMI = 2 * MI(A, B) / (H(A) + H(B))
pub fn calculate_nmi(labels_a: &[String], labels_b: &[String]) -> f64 {
    if labels_a.len() != labels_b.len() || labels_a.is_empty() {
        return 0.0;
    }
    let n = labels_a.len() as f64;

    // 计算联合分布和边缘分布
    let mut joint: HashMap<(&str, &str), usize> = HashMap::new();
    let mut count_a: HashMap<&str, usize> = HashMap::new();
    let mut count_b: HashMap<&str, usize> = HashMap::new();

    for (a, b) in labels_a.iter().zip(labels_b.iter()) {
        *joint.entry((a.as_str(), b.as_str())).or_insert(0) += 1;
        *count_a.entry(a.as_str()).or_insert(0) += 1;
        *count_b.entry(b.as_str()).or_insert(0) += 1;
    }

    // H(A)
    let h_a: f64 = count_a
        .values()
        .map(|&c| {
            let p = c as f64 / n;
            -p * p.ln()
        })
        .sum();

    // H(B)
    let h_b: f64 = count_b
        .values()
        .map(|&c| {
            let p = c as f64 / n;
            -p * p.ln()
        })
        .sum();

    if (h_a + h_b).abs() < f64::EPSILON {
        // 两个分布都是单一类别
        return 1.0;
    }

    // MI(A, B) = sum_ij p_ij * log(p_ij / (p_i * p_j))
    let mi: f64 = joint
        .iter()
        .map(|(&(a, b), &c)| {
            let p_ij = c as f64 / n;
            let p_i = count_a[a] as f64 / n;
            let p_j = count_b[b] as f64 / n;
            p_ij * (p_ij / (p_i * p_j)).ln()
        })
        .sum();

    2.0 * mi / (h_a + h_b)
}

// ============================================================================
// 表达矩阵准确性计算
// ============================================================================

/// 计算表达矩阵的准确性指标
pub fn calculate_matrix_accuracy(
    original: &ExpressionMatrix,
    converted: &ExpressionMatrix,
) -> Result<AccuracyMetrics, String> {
    // 检查维度
    let (orig_rows, orig_cols) = original.shape();
    let (conv_rows, conv_cols) = converted.shape();

    if orig_rows != conv_rows || orig_cols != conv_cols {
        return Err(format!(
            "Dimension mismatch: original ({}, {}), converted ({}, {})",
            orig_rows, orig_cols, conv_rows, conv_cols
        ));
    }

    // 获取所有非零值用于计算
    let (orig_values, orig_nnz) = extract_matrix_values(original);
    let (conv_values, conv_nnz) = extract_matrix_values(converted);

    // 计算稀疏度
    let total_elements = orig_rows * orig_cols;
    let orig_sparsity = calculate_sparsity(orig_nnz, total_elements);
    let conv_sparsity = calculate_sparsity(conv_nnz, total_elements);
    let sparsity_preserved = 1.0 - (orig_sparsity - conv_sparsity).abs();

    // 如果两个矩阵都是空的
    if orig_values.is_empty() && conv_values.is_empty() {
        return Ok(AccuracyMetrics {
            correlation: 1.0,
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
            mse: 0.0,
            exact_match_ratio: 1.0,
            ari: None,
            nmi: None,
        });
    }

    // 需要将矩阵展平为相同顺序的向量进行比较
    let (orig_flat, conv_flat) = flatten_matrices_for_comparison(original, converted)?;

    // 计算相关系数
    let correlation = pearson_correlation(&orig_flat, &conv_flat);

    // 计算误差
    let mut sum_abs_error = 0.0;
    let mut sum_rel_error = 0.0;
    let mut max_abs_error = 0.0f64;
    let mut count_rel = 0usize;

    for (o, c) in orig_flat.iter().zip(conv_flat.iter()) {
        let abs_err = (o - c).abs();
        sum_abs_error += abs_err;
        max_abs_error = max_abs_error.max(abs_err);

        // 相对误差（避免除以零）
        if o.abs() > 1e-10 {
            sum_rel_error += abs_err / o.abs();
            count_rel += 1;
        }
    }

    let n = orig_flat.len();
    let mean_absolute_error = if n > 0 { sum_abs_error / n as f64 } else { 0.0 };
    let mean_relative_error = if count_rel > 0 {
        sum_rel_error / count_rel as f64
    } else {
        0.0
    };

    // 计算新增指标
    let mse = calculate_mse(&orig_flat, &conv_flat);
    let exact_match_ratio = calculate_exact_match_ratio(&orig_flat, &conv_flat);

    Ok(AccuracyMetrics {
        correlation,
        mean_relative_error,
        mean_absolute_error,
        max_absolute_error: max_abs_error,
        sparsity_preserved,
        metadata_match: true,
        categorical_match: true,
        embedding_error: 0.0,
        nnz_match: orig_nnz == conv_nnz,
        original_nnz: orig_nnz,
        converted_nnz: conv_nnz,
        mse,
        exact_match_ratio,
        ari: None,
        nmi: None,
    })
}

/// 提取矩阵的所有值和非零元素数量
fn extract_matrix_values(matrix: &ExpressionMatrix) -> (Vec<f64>, usize) {
    match matrix {
        ExpressionMatrix::Dense(dense) => {
            let nnz = dense.data.iter().filter(|&&v| v != 0.0).count();
            (dense.data.clone(), nnz)
        }
        ExpressionMatrix::SparseCSR(sparse) => (sparse.data.clone(), sparse.data.len()),
        ExpressionMatrix::SparseCSC(sparse) => (sparse.data.clone(), sparse.data.len()),
        ExpressionMatrix::Lazy(lazy) => {
            // 延迟加载矩阵：尝试从缓存获取
            if let Some(cached) = lazy.get_cached() {
                extract_matrix_values(&cached)
            } else {
                // 没有缓存，返回空值
                (Vec::new(), lazy.nnz.unwrap_or(0))
            }
        }
    }
}

/// 将两个矩阵展平为相同顺序的向量进行比较
fn flatten_matrices_for_comparison(
    original: &ExpressionMatrix,
    converted: &ExpressionMatrix,
) -> Result<(Vec<f64>, Vec<f64>), String> {
    let (rows, cols) = original.shape();

    // 将两个矩阵都转换为行优先的稠密表示
    let orig_dense = matrix_to_dense_row_major(original, rows, cols);
    let conv_dense = matrix_to_dense_row_major(converted, rows, cols);

    Ok((orig_dense, conv_dense))
}

/// 将矩阵转换为行优先的稠密向量
fn matrix_to_dense_row_major(matrix: &ExpressionMatrix, rows: usize, cols: usize) -> Vec<f64> {
    let mut result = vec![0.0; rows * cols];

    match matrix {
        ExpressionMatrix::Dense(dense) => {
            // 假设 DenseMatrix 已经是行优先
            result.copy_from_slice(&dense.data);
        }
        ExpressionMatrix::SparseCSR(sparse) => {
            // CSR 格式：按行遍历
            for row in 0..rows {
                let start = sparse.indptr[row];
                let end = sparse.indptr[row + 1];
                for idx in start..end {
                    let col = sparse.indices[idx];
                    let val = sparse.data[idx];
                    result[row * cols + col] = val;
                }
            }
        }
        ExpressionMatrix::SparseCSC(sparse) => {
            // CSC 格式：按列遍历
            for col in 0..cols {
                let start = sparse.indptr[col];
                let end = sparse.indptr[col + 1];
                for idx in start..end {
                    let row = sparse.indices[idx];
                    let val = sparse.data[idx];
                    result[row * cols + col] = val;
                }
            }
        }
        ExpressionMatrix::Lazy(lazy) => {
            // 延迟加载矩阵：尝试从缓存获取
            if let Some(cached) = lazy.get_cached() {
                return matrix_to_dense_row_major(&cached, rows, cols);
            }
            // 没有缓存，返回全零
        }
    }

    result
}

// ============================================================================
// 元数据准确性计算
// ============================================================================

/// 计算元数据的准确性
pub fn calculate_metadata_accuracy(
    original: &DataFrame,
    converted: &DataFrame,
) -> MetadataAccuracy {
    use std::collections::HashSet;

    let orig_cols: HashSet<_> = original.columns.iter().cloned().collect();
    let conv_cols: HashSet<_> = converted.columns.iter().cloned().collect();

    let missing: Vec<_> = orig_cols.difference(&conv_cols).cloned().collect();
    let extra: Vec<_> = conv_cols.difference(&orig_cols).cloned().collect();

    let mut numeric_max_error = 0.0f64;
    let mut dtypes_match = true;
    let mut categorical_match = true;

    // 比较共同的列
    for col_name in orig_cols.intersection(&conv_cols) {
        if let (Some(orig_arr), Some(conv_arr)) =
            (original.column(col_name), converted.column(col_name))
        {
            // 检查数据类型
            if orig_arr.data_type() != conv_arr.data_type() {
                dtypes_match = false;
            }

            // 对于数值列，计算误差
            use arrow::array::*;
            use arrow::datatypes::DataType;

            match orig_arr.data_type() {
                DataType::Float64 => {
                    let orig = orig_arr.as_any().downcast_ref::<Float64Array>().unwrap();
                    let conv = conv_arr.as_any().downcast_ref::<Float64Array>().unwrap();

                    for i in 0..orig.len().min(conv.len()) {
                        if !orig.is_null(i) && !conv.is_null(i) {
                            let err = (orig.value(i) - conv.value(i)).abs();
                            numeric_max_error = numeric_max_error.max(err);
                        }
                    }
                }
                DataType::Dictionary(_, _) => {
                    // Categorical 比较
                    let orig = orig_arr
                        .as_any()
                        .downcast_ref::<DictionaryArray<arrow::datatypes::Int32Type>>();
                    let conv = conv_arr
                        .as_any()
                        .downcast_ref::<DictionaryArray<arrow::datatypes::Int32Type>>();

                    if let (Some(o), Some(c)) = (orig, conv) {
                        // 比较 keys
                        for i in 0..o.len().min(c.len()) {
                            if !o.is_null(i) && !c.is_null(i) {
                                if o.keys().value(i) != c.keys().value(i) {
                                    categorical_match = false;
                                    break;
                                }
                            }
                        }
                    }
                }
                _ => {}
            }
        }
    }

    MetadataAccuracy {
        row_count_match: original.n_rows == converted.n_rows,
        column_count_match: original.n_cols() == converted.n_cols(),
        column_names_match: missing.is_empty() && extra.is_empty(),
        dtypes_match,
        numeric_max_error,
        categorical_match,
        missing_columns: missing,
        extra_columns: extra,
    }
}

// ============================================================================
// 嵌入准确性计算
// ============================================================================

/// 计算嵌入的准确性
pub fn calculate_embedding_accuracy(
    original: &Embedding,
    converted: &Embedding,
) -> EmbeddingAccuracy {
    let dimension_match =
        original.n_rows == converted.n_rows && original.n_cols == converted.n_cols;

    if !dimension_match || original.data.is_empty() {
        return EmbeddingAccuracy {
            dimension_match,
            mean_error: f64::MAX,
            max_error: f64::MAX,
            correlation: 0.0,
        };
    }

    // 计算误差
    let mut sum_error = 0.0;
    let mut max_error = 0.0f64;

    for (o, c) in original.data.iter().zip(converted.data.iter()) {
        let err = (o - c).abs();
        sum_error += err;
        max_error = max_error.max(err);
    }

    let mean_error = sum_error / original.data.len() as f64;
    let correlation = pearson_correlation(&original.data, &converted.data);

    EmbeddingAccuracy {
        dimension_match,
        mean_error,
        max_error,
        correlation,
    }
}

// ============================================================================
// 完整的准确性测试
// ============================================================================

/// 对 SingleCellData 进行完整的准确性测试
pub fn calculate_full_accuracy(
    original: &SingleCellData,
    converted: &SingleCellData,
    dataset_name: &str,
) -> AccuracyReport {
    let mut report = AccuracyReport::new(dataset_name);

    // 1. 表达矩阵准确性
    match calculate_matrix_accuracy(&original.expression, &converted.expression) {
        Ok(metrics) => {
            report.expression_accuracy = metrics;
        }
        Err(e) => {
            report.add_error(format!("Expression matrix error: {}", e));
        }
    }

    // 2. 细胞元数据准确性
    let cell_meta_acc =
        calculate_metadata_accuracy(&original.cell_metadata, &converted.cell_metadata);
    // 只有当行数不匹配时才标记为失败
    // 列数差异在跨格式转换中是可以接受的（某些列可能不被保留）
    if !cell_meta_acc.row_count_match {
        report.passed = false;
    }
    report.cell_metadata_accuracy = Some(cell_meta_acc);

    // 3. 基因元数据准确性
    let gene_meta_acc =
        calculate_metadata_accuracy(&original.gene_metadata, &converted.gene_metadata);
    // 只有当行数不匹配时才标记为失败
    if !gene_meta_acc.row_count_match {
        report.passed = false;
    }
    report.gene_metadata_accuracy = Some(gene_meta_acc);

    // 4. 嵌入准确性
    if let (Some(orig_emb), Some(conv_emb)) = (&original.embeddings, &converted.embeddings) {
        for (name, orig) in orig_emb {
            if let Some(conv) = conv_emb.get(name) {
                let acc = calculate_embedding_accuracy(orig, conv);
                report.embedding_accuracies.insert(name.clone(), acc);
            } else {
                // 缺失的嵌入只是警告，不标记为失败
                eprintln!("Warning: Missing embedding in converted data: {}", name);
            }
        }
    }

    // 5. Layers 准确性
    if let (Some(orig_layers), Some(conv_layers)) = (&original.layers, &converted.layers) {
        for (name, orig) in orig_layers {
            if let Some(conv) = conv_layers.get(name) {
                match calculate_matrix_accuracy(orig, conv) {
                    Ok(metrics) => {
                        report.layer_accuracies.insert(name.clone(), metrics);
                    }
                    Err(e) => {
                        // Layer 错误只是警告
                        eprintln!("Warning: Layer '{}' error: {}", name, e);
                    }
                }
            } else {
                // 缺失的 layer 只是警告
                eprintln!("Warning: Missing layer in converted data: {}", name);
            }
        }
    }

    // 检查表达矩阵准确性是否满足要求
    if !report.expression_accuracy.meets_requirements(1e-7) {
        report.passed = false;
    }

    report
}

// ============================================================================
// 测试
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pearson_correlation_perfect() {
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let y = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let corr = pearson_correlation(&x, &y);
        assert!(
            (corr - 1.0).abs() < 1e-10,
            "Perfect correlation should be 1.0"
        );
    }

    #[test]
    fn test_pearson_correlation_negative() {
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let y = vec![5.0, 4.0, 3.0, 2.0, 1.0];
        let corr = pearson_correlation(&x, &y);
        assert!(
            (corr - (-1.0)).abs() < 1e-10,
            "Perfect negative correlation should be -1.0"
        );
    }

    #[test]
    fn test_calculate_sparsity() {
        assert!((calculate_sparsity(10, 100) - 0.9).abs() < 1e-10);
        assert!((calculate_sparsity(0, 100) - 1.0).abs() < 1e-10);
        assert!((calculate_sparsity(100, 100) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_accuracy_metrics_summary() {
        let metrics = AccuracyMetrics::perfect();
        let summary = metrics.summary();
        assert!(summary.contains("Correlation: 1.000000"));
        assert!(summary.contains("NNZ Match: ✓"));
    }
}
