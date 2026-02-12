//! 跨格式验证工具
//!
//! 提供 .h5ad 和 .rds 格式之间的比较和验证功能。
//! 用于验证 CrossCell 的格式转换正确性。

use crate::anndata::reader::read_h5ad;
use crate::error::CrossCellError;
use crate::ir::{ExpressionMatrix, SingleCellData};
use crate::seurat::seurat_to_ir::seurat_rds_to_ir;
use std::path::Path;

use super::ComparisonResult;
use super::compare::{
    compare_dataframe, compare_embeddings, compare_layers,
};
use super::accuracy::pearson_correlation;

/// 跨格式比较报告
#[derive(Debug)]
pub struct CrossFormatReport {
    /// 表达矩阵比较结果
    pub expression_match: ComparisonResult,
    /// 元数据比较结果
    pub metadata_match: ComparisonResult,
    /// 嵌入比较结果
    pub embeddings_match: Option<ComparisonResult>,
    /// Layers 比较结果
    pub layers_match: Option<ComparisonResult>,
    /// 表达矩阵相关系数
    pub expression_correlation: f64,
    /// H5AD 稀疏度
    pub h5ad_sparsity: f64,
    /// RDS 稀疏度
    pub rds_sparsity: f64,
    /// 总体是否通过
    pub passed: bool,
}

impl CrossFormatReport {
    /// 生成摘要
    pub fn summary(&self) -> String {
        let status = if self.passed { "✅ 通过" } else { "❌ 失败" };
        
        let mut summary = format!(
            "跨格式验证报告 {}\n\n\
             📊 表达矩阵\n\
             - 相关系数: {:.6}\n\
             - H5AD 稀疏度: {:.2}%\n\
             - RDS 稀疏度: {:.2}%\n\
             - 状态: {}\n",
            status,
            self.expression_correlation,
            self.h5ad_sparsity * 100.0,
            self.rds_sparsity * 100.0,
            if self.expression_match.passed { "✓" } else { "✗" }
        );
        
        summary.push_str(&format!(
            "\n📋 元数据\n- 状态: {}\n",
            if self.metadata_match.passed { "✓" } else { "✗" }
        ));
        
        if let Some(ref emb) = self.embeddings_match {
            summary.push_str(&format!(
                "\n🎯 嵌入\n- 状态: {}\n",
                if emb.passed { "✓" } else { "✗" }
            ));
        }
        
        if let Some(ref layers) = self.layers_match {
            summary.push_str(&format!(
                "\n📚 Layers\n- 状态: {}\n",
                if layers.passed { "✓" } else { "✗" }
            ));
        }
        
        summary
    }
}

/// 比较 .h5ad 和 .rds 文件
///
/// # 参数
/// - `h5ad_path`: H5AD 文件路径
/// - `rds_path`: RDS 文件路径
/// - `tolerance`: 数值容差（默认 1e-7）
///
/// # 返回
/// - `CrossFormatReport`: 跨格式比较报告
pub fn compare_h5ad_rds<P: AsRef<Path>>(
    h5ad_path: P,
    rds_path: P,
    tolerance: f64,
) -> Result<CrossFormatReport, CrossCellError> {
    // 读取 H5AD 文件
    let h5ad_ir = read_h5ad(h5ad_path.as_ref())?;
    
    // 读取 RDS 文件
    let rds_path_str = rds_path.as_ref().to_str()
        .ok_or_else(|| CrossCellError::InvalidFormat("Invalid RDS path".to_string()))?;
    let rds_ir = seurat_rds_to_ir(rds_path_str)?;
    
    // 比较两个 IR
    compare_ir_cross_format(&h5ad_ir, &rds_ir, tolerance)
}

/// 比较两个 IR（跨格式）
///
/// 考虑格式差异（如 CSR vs CSC）进行比较
pub fn compare_ir_cross_format(
    h5ad_ir: &SingleCellData,
    rds_ir: &SingleCellData,
    tolerance: f64,
) -> Result<CrossFormatReport, CrossCellError> {
    // 计算稀疏度
    let h5ad_sparsity = calculate_matrix_sparsity(&h5ad_ir.expression);
    let rds_sparsity = calculate_matrix_sparsity(&rds_ir.expression);
    
    // 计算相关系数
    let expression_correlation = calculate_expression_correlation(
        &h5ad_ir.expression,
        &rds_ir.expression,
    )?;
    
    // 比较表达矩阵（考虑格式差异）
    let expression_match = compare_expression_cross_format(
        &h5ad_ir.expression,
        &rds_ir.expression,
        tolerance,
    );
    
    // 比较元数据
    let metadata_match = compare_dataframe(
        &h5ad_ir.cell_metadata,
        &rds_ir.cell_metadata,
        tolerance,
    );
    
    // 比较嵌入
    let embeddings_match = match (&h5ad_ir.embeddings, &rds_ir.embeddings) {
        (Some(h5ad_emb), Some(rds_emb)) => {
            Some(compare_embeddings(h5ad_emb, rds_emb, tolerance))
        }
        (None, None) => None,
        _ => Some(ComparisonResult {
            passed: false,
            message: "嵌入存在性不匹配".to_string(),
            max_difference: None,
        }),
    };
    
    // 比较 layers
    let layers_match = match (&h5ad_ir.layers, &rds_ir.layers) {
        (Some(h5ad_layers), Some(rds_layers)) => {
            Some(compare_layers(h5ad_layers, rds_layers, tolerance))
        }
        (None, None) => None,
        _ => Some(ComparisonResult {
            passed: false,
            message: "Layers 存在性不匹配".to_string(),
            max_difference: None,
        }),
    };
    
    // 判断总体是否通过
    let passed = expression_match.passed
        && metadata_match.passed
        && embeddings_match.as_ref().map_or(true, |e| e.passed)
        && layers_match.as_ref().map_or(true, |l| l.passed)
        && expression_correlation > 0.9999;
    
    Ok(CrossFormatReport {
        expression_match,
        metadata_match,
        embeddings_match,
        layers_match,
        expression_correlation,
        h5ad_sparsity,
        rds_sparsity,
        passed,
    })
}

/// 比较表达矩阵（跨格式）
///
/// 考虑 CSR 和 CSC 格式差异，将两者转换为相同格式后比较
fn compare_expression_cross_format(
    h5ad_expr: &ExpressionMatrix,
    rds_expr: &ExpressionMatrix,
    tolerance: f64,
) -> ComparisonResult {
    // 检查维度
    let (h5ad_rows, h5ad_cols) = h5ad_expr.shape();
    let (rds_rows, rds_cols) = rds_expr.shape();
    
    if h5ad_rows != rds_rows || h5ad_cols != rds_cols {
        return ComparisonResult {
            passed: false,
            message: format!(
                "维度不匹配: H5AD ({}, {}), RDS ({}, {})",
                h5ad_rows, h5ad_cols, rds_rows, rds_cols
            ),
            max_difference: None,
        };
    }
    
    // 如果格式相同，直接比较
    match (h5ad_expr, rds_expr) {
        (ExpressionMatrix::SparseCSR(h5ad), ExpressionMatrix::SparseCSR(rds)) => {
            compare_sparse_csr_values(h5ad, rds, tolerance)
        }
        (ExpressionMatrix::SparseCSC(h5ad), ExpressionMatrix::SparseCSC(rds)) => {
            compare_sparse_csc_values(h5ad, rds, tolerance)
        }
        (ExpressionMatrix::Dense(h5ad), ExpressionMatrix::Dense(rds)) => {
            compare_dense_values(h5ad, rds, tolerance)
        }
        // 不同格式：转换后比较
        (ExpressionMatrix::SparseCSR(h5ad), ExpressionMatrix::SparseCSC(rds)) => {
            // 将 CSR 转换为 CSC 后比较
            let h5ad_csc = crate::sparse::convert::csr_to_csc(h5ad);
            compare_sparse_csc_values(&h5ad_csc, rds, tolerance)
        }
        (ExpressionMatrix::SparseCSC(h5ad), ExpressionMatrix::SparseCSR(rds)) => {
            // 将 CSC 转换为 CSR 后比较
            let h5ad_csr = crate::sparse::convert::csc_to_csr(h5ad);
            compare_sparse_csr_values(&h5ad_csr, rds, tolerance)
        }
        _ => ComparisonResult {
            passed: false,
            message: "不支持的矩阵格式组合".to_string(),
            max_difference: None,
        },
    }
}

/// 比较 CSR 矩阵值
fn compare_sparse_csr_values(
    a: &crate::ir::SparseMatrixCSR,
    b: &crate::ir::SparseMatrixCSR,
    tolerance: f64,
) -> ComparisonResult {
    // 比较非零元素数量
    if a.data.len() != b.data.len() {
        return ComparisonResult {
            passed: false,
            message: format!(
                "非零元素数量不匹配: {} vs {}",
                a.data.len(),
                b.data.len()
            ),
            max_difference: None,
        };
    }
    
    // 比较索引
    if a.indices != b.indices {
        return ComparisonResult {
            passed: false,
            message: "列索引不匹配".to_string(),
            max_difference: None,
        };
    }
    
    if a.indptr != b.indptr {
        return ComparisonResult {
            passed: false,
            message: "行指针不匹配".to_string(),
            max_difference: None,
        };
    }
    
    // 比较数值
    let max_diff = a
        .data
        .iter()
        .zip(b.data.iter())
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max);
    
    let passed = max_diff < tolerance;
    
    ComparisonResult {
        passed,
        message: if passed {
            format!("CSR 矩阵匹配 (最大差异: {:.2e})", max_diff)
        } else {
            format!(
                "CSR 矩阵不匹配 (最大差异: {:.2e} > 容差 {:.2e})",
                max_diff, tolerance
            )
        },
        max_difference: Some(max_diff),
    }
}

/// 比较 CSC 矩阵值
fn compare_sparse_csc_values(
    a: &crate::ir::SparseMatrixCSC,
    b: &crate::ir::SparseMatrixCSC,
    tolerance: f64,
) -> ComparisonResult {
    // 比较非零元素数量
    if a.data.len() != b.data.len() {
        return ComparisonResult {
            passed: false,
            message: format!(
                "非零元素数量不匹配: {} vs {}",
                a.data.len(),
                b.data.len()
            ),
            max_difference: None,
        };
    }
    
    // 比较索引
    if a.indices != b.indices {
        return ComparisonResult {
            passed: false,
            message: "行索引不匹配".to_string(),
            max_difference: None,
        };
    }
    
    if a.indptr != b.indptr {
        return ComparisonResult {
            passed: false,
            message: "列指针不匹配".to_string(),
            max_difference: None,
        };
    }
    
    // 比较数值
    let max_diff = a
        .data
        .iter()
        .zip(b.data.iter())
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max);
    
    let passed = max_diff < tolerance;
    
    ComparisonResult {
        passed,
        message: if passed {
            format!("CSC 矩阵匹配 (最大差异: {:.2e})", max_diff)
        } else {
            format!(
                "CSC 矩阵不匹配 (最大差异: {:.2e} > 容差 {:.2e})",
                max_diff, tolerance
            )
        },
        max_difference: Some(max_diff),
    }
}

/// 比较稠密矩阵值
fn compare_dense_values(
    a: &crate::ir::DenseMatrix,
    b: &crate::ir::DenseMatrix,
    tolerance: f64,
) -> ComparisonResult {
    if a.data.len() != b.data.len() {
        return ComparisonResult {
            passed: false,
            message: format!(
                "数据长度不匹配: {} vs {}",
                a.data.len(),
                b.data.len()
            ),
            max_difference: None,
        };
    }
    
    let max_diff = a
        .data
        .iter()
        .zip(b.data.iter())
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max);
    
    let passed = max_diff < tolerance;
    
    ComparisonResult {
        passed,
        message: if passed {
            format!("稠密矩阵匹配 (最大差异: {:.2e})", max_diff)
        } else {
            format!(
                "稠密矩阵不匹配 (最大差异: {:.2e} > 容差 {:.2e})",
                max_diff, tolerance
            )
        },
        max_difference: Some(max_diff),
    }
}

/// 计算表达矩阵相关系数
///
/// 将两个矩阵转换为相同格式后，按位置比较数据
fn calculate_expression_correlation(
    a: &ExpressionMatrix,
    b: &ExpressionMatrix,
) -> Result<f64, CrossCellError> {
    let (a_rows, a_cols) = a.shape();
    let (b_rows, b_cols) = b.shape();
    
    if a_rows != b_rows || a_cols != b_cols {
        return Err(CrossCellError::ValidationFailed(
            format!("矩阵维度不匹配: ({}, {}) vs ({}, {})", a_rows, a_cols, b_rows, b_cols),
        ));
    }
    
    if a_rows == 0 || a_cols == 0 {
        return Ok(1.0); // 空矩阵视为完全匹配
    }
    
    // 将两个矩阵都转换为 CSC 格式，然后按位置提取数据
    let a_csc = to_csc(a);
    let b_csc = to_csc(b);
    
    // 构建位置到值的映射
    let a_values = sparse_to_dense_vector(&a_csc);
    let b_values = sparse_to_dense_vector(&b_csc);
    
    if a_values.len() != b_values.len() {
        return Err(CrossCellError::ValidationFailed(
            "转换后数据长度不匹配".to_string(),
        ));
    }
    
    Ok(pearson_correlation(&a_values, &b_values))
}

/// 将表达矩阵转换为 CSC 格式
fn to_csc(matrix: &ExpressionMatrix) -> crate::ir::SparseMatrixCSC {
    match matrix {
        ExpressionMatrix::SparseCSC(m) => m.clone(),
        ExpressionMatrix::SparseCSR(m) => crate::sparse::convert::csr_to_csc(m),
        ExpressionMatrix::Dense(m) => dense_to_csc(m),
        ExpressionMatrix::Lazy(_) => {
            // Lazy 矩阵返回空的 CSC
            crate::ir::SparseMatrixCSC {
                n_rows: 0,
                n_cols: 0,
                data: Vec::new(),
                indices: Vec::new(),
                indptr: vec![0],
            }
        }
    }
}

/// 将稠密矩阵转换为 CSC 格式
fn dense_to_csc(dense: &crate::ir::DenseMatrix) -> crate::ir::SparseMatrixCSC {
    let mut data = Vec::new();
    let mut indices = Vec::new();
    let mut indptr = vec![0usize];
    
    // 按列遍历（CSC 是列优先）
    // DenseMatrix 是行优先存储
    for col in 0..dense.n_cols {
        for row in 0..dense.n_rows {
            // 行优先索引
            let idx = row * dense.n_cols + col;
            let val = dense.data[idx];
            if val != 0.0 {
                data.push(val);
                indices.push(row);
            }
        }
        indptr.push(data.len());
    }
    
    crate::ir::SparseMatrixCSC {
        n_rows: dense.n_rows,
        n_cols: dense.n_cols,
        data,
        indices,
        indptr,
    }
}

/// 将稀疏 CSC 矩阵转换为稠密向量（按列优先顺序）
///
/// 对于大矩阵，只采样部分数据以避免内存问题
fn sparse_to_dense_vector(csc: &crate::ir::SparseMatrixCSC) -> Vec<f64> {
    let total_elements = csc.n_rows * csc.n_cols;
    
    // 如果矩阵太大，使用采样方法
    if total_elements > 10_000_000 {
        // 对于大矩阵，只比较非零元素
        // 但需要确保两个矩阵的非零元素在相同位置
        return sample_sparse_values(csc, 1_000_000);
    }
    
    // 对于小矩阵，转换为完整的稠密向量
    let mut dense = vec![0.0; total_elements];
    
    for col in 0..csc.n_cols {
        let start = csc.indptr[col];
        let end = csc.indptr[col + 1];
        
        for i in start..end {
            let row = csc.indices[i];
            let val = csc.data[i];
            // 列优先索引
            let idx = col * csc.n_rows + row;
            dense[idx] = val;
        }
    }
    
    dense
}

/// 从稀疏矩阵采样值（用于大矩阵）
fn sample_sparse_values(csc: &crate::ir::SparseMatrixCSC, max_samples: usize) -> Vec<f64> {
    let nnz = csc.data.len();
    
    if nnz <= max_samples {
        // 如果非零元素不多，直接返回所有值
        // 但需要按位置排序
        let mut indexed_values: Vec<(usize, f64)> = Vec::with_capacity(nnz);
        
        for col in 0..csc.n_cols {
            let start = csc.indptr[col];
            let end = csc.indptr[col + 1];
            
            for i in start..end {
                let row = csc.indices[i];
                let val = csc.data[i];
                let idx = col * csc.n_rows + row;
                indexed_values.push((idx, val));
            }
        }
        
        // 按位置排序
        indexed_values.sort_by_key(|(idx, _)| *idx);
        indexed_values.into_iter().map(|(_, v)| v).collect()
    } else {
        // 均匀采样
        let step = nnz / max_samples;
        csc.data.iter().step_by(step).copied().collect()
    }
}

/// 计算矩阵稀疏度
fn calculate_matrix_sparsity(matrix: &ExpressionMatrix) -> f64 {
    let (rows, cols) = matrix.shape();
    let total = rows * cols;
    
    if total == 0 {
        return 1.0;
    }
    
    let nnz = match matrix {
        ExpressionMatrix::SparseCSR(m) => m.data.len(),
        ExpressionMatrix::SparseCSC(m) => m.data.len(),
        ExpressionMatrix::Dense(m) => m.data.iter().filter(|&&x| x != 0.0).count(),
        ExpressionMatrix::Lazy(_) => 0, // Lazy 矩阵暂不支持
    };
    
    1.0 - (nnz as f64 / total as f64)
}


// ============================================================================
// 跨格式验证函数
// ============================================================================

/// 验证 .h5ad → .rds 转换正确性
///
/// # 参数
/// - `original_h5ad`: 原始 H5AD 文件路径
/// - `converted_rds`: 转换后的 RDS 文件路径
/// - `tolerance`: 数值容差
///
/// # 返回
/// - `CrossFormatValidationResult`: 验证结果
pub fn validate_h5ad_to_rds<P: AsRef<Path>>(
    original_h5ad: P,
    converted_rds: P,
    tolerance: f64,
) -> Result<CrossFormatValidationResult, CrossCellError> {
    let report = compare_h5ad_rds(&original_h5ad, &converted_rds, tolerance)?;
    
    Ok(CrossFormatValidationResult {
        direction: ConversionDirection::H5adToRds,
        original_path: original_h5ad.as_ref().to_string_lossy().to_string(),
        converted_path: converted_rds.as_ref().to_string_lossy().to_string(),
        report,
    })
}

/// 验证 .rds → .h5ad 转换正确性
///
/// # 参数
/// - `original_rds`: 原始 RDS 文件路径
/// - `converted_h5ad`: 转换后的 H5AD 文件路径
/// - `tolerance`: 数值容差
///
/// # 返回
/// - `CrossFormatValidationResult`: 验证结果
pub fn validate_rds_to_h5ad<P: AsRef<Path>>(
    original_rds: P,
    converted_h5ad: P,
    tolerance: f64,
) -> Result<CrossFormatValidationResult, CrossCellError> {
    // 注意：这里交换参数顺序，因为 compare_h5ad_rds 期望 h5ad 在前
    let report = compare_h5ad_rds(&converted_h5ad, &original_rds, tolerance)?;
    
    Ok(CrossFormatValidationResult {
        direction: ConversionDirection::RdsToH5ad,
        original_path: original_rds.as_ref().to_string_lossy().to_string(),
        converted_path: converted_h5ad.as_ref().to_string_lossy().to_string(),
        report,
    })
}

/// 验证完整往返转换
///
/// .h5ad → .rds → .h5ad 或 .rds → .h5ad → .rds
///
/// # 参数
/// - `original`: 原始文件路径
/// - `roundtrip`: 往返后的文件路径
/// - `tolerance`: 数值容差
///
/// # 返回
/// - `RoundtripValidationResult`: 往返验证结果
pub fn validate_roundtrip<P: AsRef<Path>>(
    original: P,
    roundtrip: P,
    tolerance: f64,
) -> Result<RoundtripValidationResult, CrossCellError> {
    let original_path = original.as_ref();
    let roundtrip_path = roundtrip.as_ref();
    
    // 根据文件扩展名确定格式
    let original_ext = original_path.extension()
        .and_then(|e| e.to_str())
        .unwrap_or("");
    let roundtrip_ext = roundtrip_path.extension()
        .and_then(|e| e.to_str())
        .unwrap_or("");
    
    // 验证扩展名匹配
    if original_ext != roundtrip_ext {
        return Err(CrossCellError::InvalidFormat(format!(
            "往返验证要求相同格式: {} vs {}",
            original_ext, roundtrip_ext
        )));
    }
    
    // 读取两个文件
    let (original_ir, roundtrip_ir) = match original_ext {
        "h5ad" => {
            let orig = read_h5ad(original_path)?;
            let rt = read_h5ad(roundtrip_path)?;
            (orig, rt)
        }
        "rds" => {
            let orig_str = original_path.to_str()
                .ok_or_else(|| CrossCellError::InvalidFormat("Invalid path".to_string()))?;
            let rt_str = roundtrip_path.to_str()
                .ok_or_else(|| CrossCellError::InvalidFormat("Invalid path".to_string()))?;
            let orig = seurat_rds_to_ir(orig_str)?;
            let rt = seurat_rds_to_ir(rt_str)?;
            (orig, rt)
        }
        _ => {
            return Err(CrossCellError::UnsupportedFormat {
                expected: "h5ad or rds".to_string(),
                actual: original_ext.to_string(),
            });
        }
    };
    
    // 比较两个 IR
    let report = compare_ir_cross_format(&original_ir, &roundtrip_ir, tolerance)?;
    
    Ok(RoundtripValidationResult {
        format: original_ext.to_string(),
        original_path: original_path.to_string_lossy().to_string(),
        roundtrip_path: roundtrip_path.to_string_lossy().to_string(),
        report,
    })
}

/// 转换方向
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConversionDirection {
    /// H5AD → RDS
    H5adToRds,
    /// RDS → H5AD
    RdsToH5ad,
}

impl std::fmt::Display for ConversionDirection {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ConversionDirection::H5adToRds => write!(f, "H5AD → RDS"),
            ConversionDirection::RdsToH5ad => write!(f, "RDS → H5AD"),
        }
    }
}

/// 跨格式验证结果
#[derive(Debug)]
pub struct CrossFormatValidationResult {
    /// 转换方向
    pub direction: ConversionDirection,
    /// 原始文件路径
    pub original_path: String,
    /// 转换后文件路径
    pub converted_path: String,
    /// 比较报告
    pub report: CrossFormatReport,
}

impl CrossFormatValidationResult {
    /// 是否通过验证
    pub fn passed(&self) -> bool {
        self.report.passed
    }
    
    /// 生成摘要
    pub fn summary(&self) -> String {
        let status = if self.passed() { "✅ 通过" } else { "❌ 失败" };
        format!(
            "跨格式验证 {} {}\n\
             原始文件: {}\n\
             转换文件: {}\n\n\
             {}",
            self.direction,
            status,
            self.original_path,
            self.converted_path,
            self.report.summary()
        )
    }
}

/// 往返验证结果
#[derive(Debug)]
pub struct RoundtripValidationResult {
    /// 文件格式
    pub format: String,
    /// 原始文件路径
    pub original_path: String,
    /// 往返后文件路径
    pub roundtrip_path: String,
    /// 比较报告
    pub report: CrossFormatReport,
}

impl RoundtripValidationResult {
    /// 是否通过验证
    pub fn passed(&self) -> bool {
        self.report.passed
    }
    
    /// 生成摘要
    pub fn summary(&self) -> String {
        let status = if self.passed() { "✅ 通过" } else { "❌ 失败" };
        format!(
            "往返验证 (.{}) {}\n\
             原始文件: {}\n\
             往返文件: {}\n\n\
             {}",
            self.format,
            status,
            self.original_path,
            self.roundtrip_path,
            self.report.summary()
        )
    }
}

// ============================================================================
// 便捷函数
// ============================================================================

/// 默认容差值
pub const DEFAULT_TOLERANCE: f64 = 1e-7;

/// 使用默认容差比较 H5AD 和 RDS 文件
pub fn compare_h5ad_rds_default<P: AsRef<Path>>(
    h5ad_path: P,
    rds_path: P,
) -> Result<CrossFormatReport, CrossCellError> {
    compare_h5ad_rds(h5ad_path, rds_path, DEFAULT_TOLERANCE)
}

/// 使用默认容差验证 H5AD → RDS 转换
pub fn validate_h5ad_to_rds_default<P: AsRef<Path>>(
    original_h5ad: P,
    converted_rds: P,
) -> Result<CrossFormatValidationResult, CrossCellError> {
    validate_h5ad_to_rds(original_h5ad, converted_rds, DEFAULT_TOLERANCE)
}

/// 使用默认容差验证 RDS → H5AD 转换
pub fn validate_rds_to_h5ad_default<P: AsRef<Path>>(
    original_rds: P,
    converted_h5ad: P,
) -> Result<CrossFormatValidationResult, CrossCellError> {
    validate_rds_to_h5ad(original_rds, converted_h5ad, DEFAULT_TOLERANCE)
}

/// 使用默认容差验证往返转换
pub fn validate_roundtrip_default<P: AsRef<Path>>(
    original: P,
    roundtrip: P,
) -> Result<RoundtripValidationResult, CrossCellError> {
    validate_roundtrip(original, roundtrip, DEFAULT_TOLERANCE)
}
