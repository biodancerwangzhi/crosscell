//! AI-Ready 数据变换
//!
//! 提供库大小归一化、高变基因筛选、基因标识符替换等变换，
//! 用于将单细胞数据转换为 AI 基础模型可直接使用的格式。

use crate::ir::{DenseMatrix, ExpressionMatrix, SingleCellData, SparseMatrixCSC, SparseMatrixCSR};
use log::warn;

/// 库大小归一化 + log1p
///
/// 对每个细胞（行）：
/// 1. 计算总计数 total
/// 2. 缩放: x_norm = x * 10000 / total
/// 3. 变换: x_final = ln(1 + x_norm)
/// 全零行保持为零（避免除零）
pub fn normalize_library_size(matrix: &ExpressionMatrix) -> Result<ExpressionMatrix, String> {
    match matrix {
        ExpressionMatrix::SparseCSR(csr) => normalize_csr(csr).map(ExpressionMatrix::SparseCSR),
        ExpressionMatrix::SparseCSC(csc) => normalize_csc(csc).map(ExpressionMatrix::SparseCSC),
        ExpressionMatrix::Dense(dense) => normalize_dense(dense).map(ExpressionMatrix::Dense),
        ExpressionMatrix::Lazy(_) => Err("Cannot normalize lazy matrix; load it first".to_string()),
    }
}

fn normalize_csr(csr: &SparseMatrixCSR) -> Result<SparseMatrixCSR, String> {
    let mut new_data = csr.data.clone();
    for row in 0..csr.n_rows {
        let start = csr.indptr[row];
        let end = csr.indptr[row + 1];
        let total: f64 = csr.data[start..end].iter().sum();
        if total > 0.0 {
            for idx in start..end {
                new_data[idx] = (csr.data[idx] * 10000.0 / total + 1.0).ln();
            }
        }
    }
    SparseMatrixCSR::new(
        new_data,
        csr.indices.clone(),
        csr.indptr.clone(),
        csr.n_rows,
        csr.n_cols,
    )
}

fn normalize_csc(csc: &SparseMatrixCSC) -> Result<SparseMatrixCSC, String> {
    // For CSC we need row totals first, then apply per-element
    let mut row_totals = vec![0.0f64; csc.n_rows];
    for col in 0..csc.n_cols {
        let start = csc.indptr[col];
        let end = csc.indptr[col + 1];
        for idx in start..end {
            row_totals[csc.indices[idx]] += csc.data[idx];
        }
    }
    let mut new_data = csc.data.clone();
    for col in 0..csc.n_cols {
        let start = csc.indptr[col];
        let end = csc.indptr[col + 1];
        for idx in start..end {
            let row = csc.indices[idx];
            let total = row_totals[row];
            if total > 0.0 {
                new_data[idx] = (csc.data[idx] * 10000.0 / total + 1.0).ln();
            }
        }
    }
    SparseMatrixCSC::new(
        new_data,
        csc.indices.clone(),
        csc.indptr.clone(),
        csc.n_rows,
        csc.n_cols,
    )
}

fn normalize_dense(dense: &DenseMatrix) -> Result<DenseMatrix, String> {
    let mut new_data = dense.data.clone();
    for row in 0..dense.n_rows {
        let offset = row * dense.n_cols;
        let total: f64 = new_data[offset..offset + dense.n_cols].iter().sum();
        if total > 0.0 {
            for col in 0..dense.n_cols {
                let v = new_data[offset + col];
                new_data[offset + col] = (v * 10000.0 / total + 1.0).ln();
            }
        }
    }
    DenseMatrix::new(new_data, dense.n_rows, dense.n_cols)
}

/// 选择方差最高的前 N 个基因
///
/// 计算每列（基因）方差，按方差降序排序，保留前 N 列。
/// 同时裁剪表达矩阵和 var DataFrame。
/// N > 基因数时输出警告并保留所有基因。
pub fn select_top_variable_genes(data: &mut SingleCellData, n_genes: usize) -> Result<(), String> {
    let (n_rows, n_cols) = data.expression.shape();

    if n_genes >= n_cols {
        warn!(
            "Requested {} top genes but only {} genes available; keeping all",
            n_genes, n_cols
        );
        return Ok(());
    }

    // Compute per-column variance
    let variances = compute_column_variances(&data.expression, n_rows, n_cols)?;

    // Sort indices by variance descending
    let mut indices: Vec<usize> = (0..n_cols).collect();
    indices.sort_by(|&a, &b| {
        variances[b]
            .partial_cmp(&variances[a])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let keep: Vec<usize> = indices[..n_genes].to_vec();

    // Sort kept indices for consistent column ordering
    let mut keep_sorted = keep.clone();
    keep_sorted.sort();

    // Filter expression matrix columns
    data.expression = filter_matrix_columns(&data.expression, &keep_sorted, n_rows, n_cols)?;

    // Filter var DataFrame rows
    data.gene_metadata = filter_dataframe_rows(&data.gene_metadata, &keep_sorted)?;

    // Update metadata
    data.metadata.n_genes = n_genes;

    // Filter layers if present
    if let Some(ref mut layers) = data.layers {
        let mut new_layers = std::collections::HashMap::new();
        for (name, matrix) in layers.iter() {
            let (lr, lc) = matrix.shape();
            if lc == n_cols {
                new_layers.insert(
                    name.clone(),
                    filter_matrix_columns(matrix, &keep_sorted, lr, lc)?,
                );
            } else {
                new_layers.insert(name.clone(), matrix.clone());
            }
        }
        *layers = new_layers;
    }

    Ok(())
}

fn compute_column_variances(
    matrix: &ExpressionMatrix,
    n_rows: usize,
    n_cols: usize,
) -> Result<Vec<f64>, String> {
    let mut means = vec![0.0f64; n_cols];
    let mut m2 = vec![0.0f64; n_cols]; // sum of squared deviations

    if n_rows == 0 {
        return Ok(vec![0.0; n_cols]);
    }

    match matrix {
        ExpressionMatrix::Dense(dense) => {
            for row in 0..n_rows {
                let offset = row * n_cols;
                for col in 0..n_cols {
                    means[col] += dense.data[offset + col];
                }
            }
            for col in 0..n_cols {
                means[col] /= n_rows as f64;
            }
            for row in 0..n_rows {
                let offset = row * n_cols;
                for col in 0..n_cols {
                    let diff = dense.data[offset + col] - means[col];
                    m2[col] += diff * diff;
                }
            }
        }
        ExpressionMatrix::SparseCSR(csr) => {
            // First pass: compute sums
            for row in 0..n_rows {
                let start = csr.indptr[row];
                let end = csr.indptr[row + 1];
                for idx in start..end {
                    means[csr.indices[idx]] += csr.data[idx];
                }
            }
            for col in 0..n_cols {
                means[col] /= n_rows as f64;
            }
            // Second pass: compute sum of squared deviations
            // For sparse: non-zero elements contribute (val - mean)^2,
            // zero elements contribute mean^2
            let mut nnz_per_col = vec![0usize; n_cols];
            for row in 0..n_rows {
                let start = csr.indptr[row];
                let end = csr.indptr[row + 1];
                for idx in start..end {
                    let col = csr.indices[idx];
                    let diff = csr.data[idx] - means[col];
                    m2[col] += diff * diff;
                    nnz_per_col[col] += 1;
                }
            }
            // Add contribution from zero elements
            for col in 0..n_cols {
                let n_zeros = n_rows - nnz_per_col[col];
                m2[col] += means[col] * means[col] * n_zeros as f64;
            }
        }
        ExpressionMatrix::SparseCSC(csc) => {
            for col in 0..n_cols {
                let start = csc.indptr[col];
                let end = csc.indptr[col + 1];
                let nnz = end - start;
                let sum: f64 = csc.data[start..end].iter().sum();
                let mean = sum / n_rows as f64;
                means[col] = mean;
                let mut sq_dev: f64 = csc.data[start..end]
                    .iter()
                    .map(|&v| (v - mean).powi(2))
                    .sum();
                sq_dev += mean * mean * (n_rows - nnz) as f64;
                m2[col] = sq_dev;
            }
        }
        ExpressionMatrix::Lazy(_) => {
            return Err("Cannot compute variance on lazy matrix".to_string());
        }
    }

    // variance = m2 / n_rows
    Ok(m2.iter().map(|&v| v / n_rows as f64).collect())
}

/// Filter matrix to keep only specified columns (sorted indices)
fn filter_matrix_columns(
    matrix: &ExpressionMatrix,
    keep_cols: &[usize],
    n_rows: usize,
    _n_cols: usize,
) -> Result<ExpressionMatrix, String> {
    let new_n_cols = keep_cols.len();
    // Build reverse mapping: old_col -> new_col
    let mut col_map = std::collections::HashMap::new();
    for (new_idx, &old_idx) in keep_cols.iter().enumerate() {
        col_map.insert(old_idx, new_idx);
    }

    match matrix {
        ExpressionMatrix::Dense(dense) => {
            let mut new_data = Vec::with_capacity(n_rows * new_n_cols);
            for row in 0..n_rows {
                let offset = row * dense.n_cols;
                for &col in keep_cols {
                    new_data.push(dense.data[offset + col]);
                }
            }
            Ok(ExpressionMatrix::Dense(DenseMatrix::new(
                new_data, n_rows, new_n_cols,
            )?))
        }
        ExpressionMatrix::SparseCSR(csr) => {
            let mut new_data = Vec::new();
            let mut new_indices = Vec::new();
            let mut new_indptr = vec![0usize];
            for row in 0..n_rows {
                let start = csr.indptr[row];
                let end = csr.indptr[row + 1];
                for idx in start..end {
                    if let Some(&new_col) = col_map.get(&csr.indices[idx]) {
                        new_data.push(csr.data[idx]);
                        new_indices.push(new_col);
                    }
                }
                new_indptr.push(new_data.len());
            }
            Ok(ExpressionMatrix::SparseCSR(SparseMatrixCSR::new(
                new_data,
                new_indices,
                new_indptr,
                n_rows,
                new_n_cols,
            )?))
        }
        ExpressionMatrix::SparseCSC(csc) => {
            let mut new_data = Vec::new();
            let mut new_indices = Vec::new();
            let mut new_indptr = vec![0usize];
            for &old_col in keep_cols {
                let start = csc.indptr[old_col];
                let end = csc.indptr[old_col + 1];
                for idx in start..end {
                    new_data.push(csc.data[idx]);
                    new_indices.push(csc.indices[idx]);
                }
                new_indptr.push(new_data.len());
            }
            Ok(ExpressionMatrix::SparseCSC(SparseMatrixCSC::new(
                new_data,
                new_indices,
                new_indptr,
                n_rows,
                new_n_cols,
            )?))
        }
        ExpressionMatrix::Lazy(_) => Err("Cannot filter columns on lazy matrix".to_string()),
    }
}

/// Filter DataFrame rows by indices (keep only specified rows)
fn filter_dataframe_rows(
    df: &crate::ir::DataFrame,
    keep_indices: &[usize],
) -> Result<crate::ir::DataFrame, String> {
    use arrow::array::*;
    use arrow::compute::take;
    use std::sync::Arc;

    let indices_array =
        UInt32Array::from(keep_indices.iter().map(|&i| i as u32).collect::<Vec<u32>>());

    let new_data: Vec<arrow::array::ArrayRef> = df
        .data
        .iter()
        .map(|arr| {
            take(arr.as_ref(), &indices_array, None)
                .map(|a| Arc::from(a) as ArrayRef)
                .map_err(|e| format!("Failed to filter column: {}", e))
        })
        .collect::<Result<Vec<_>, _>>()?;

    crate::ir::DataFrame::new(df.columns.clone(), new_data, keep_indices.len())
}

/// 替换基因标识符
///
/// 从 var DataFrame 中读取指定列的值，将其设为基因标识符。
/// 具体做法：将指定列重命名为 "_index"（如果存在旧的 _index 列则替换），
/// 使该列成为基因名称的来源。
///
/// 如果指定列不存在，返回错误。
pub fn apply_gene_id_column(data: &mut SingleCellData, column_name: &str) -> Result<(), String> {
    use arrow::array::*;

    // Check if column exists
    let col_idx = data
        .gene_metadata
        .column_index(column_name)
        .ok_or_else(|| {
            format!(
                "Column '{}' not found in gene metadata. Available columns: {:?}",
                column_name, data.gene_metadata.columns
            )
        })?;

    // Get the column values as strings
    let col_data = &data.gene_metadata.data[col_idx];
    let string_values: Vec<String> =
        if let Some(arr) = col_data.as_any().downcast_ref::<StringArray>() {
            (0..arr.len())
                .map(|i| {
                    if arr.is_null(i) {
                        String::new()
                    } else {
                        arr.value(i).to_string()
                    }
                })
                .collect()
        } else if let Some(arr) = col_data.as_any().downcast_ref::<LargeStringArray>() {
            (0..arr.len())
                .map(|i| {
                    if arr.is_null(i) {
                        String::new()
                    } else {
                        arr.value(i).to_string()
                    }
                })
                .collect()
        } else {
            // Try to convert other types to string representation
            let arr = col_data.as_ref();
            (0..arr.len())
                .map(|i| {
                    if arr.is_null(i) {
                        String::new()
                    } else {
                        // Use arrow's display formatting
                        arrow::util::display::array_value_to_string(arr, i).unwrap_or_default()
                    }
                })
                .collect()
        };

    // Replace or add "_index" column with the new gene identifiers
    let new_index_array: ArrayRef = std::sync::Arc::new(StringArray::from(string_values));

    if let Some(idx) = data.gene_metadata.column_index("_index") {
        // Replace existing _index column
        data.gene_metadata.data[idx] = new_index_array;
    } else {
        // Add _index column at the beginning
        data.gene_metadata.columns.insert(0, "_index".to_string());
        data.gene_metadata.data.insert(0, new_index_array);
    }

    Ok(())
}
