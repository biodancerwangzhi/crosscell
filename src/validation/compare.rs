//! 数据比较工具
//!
//! 提供 IR 数据结构的比较功能，用于验证往返转换的一致性。

use crate::ir::{SingleCellData, ExpressionMatrix, DataFrame, Embedding, PairwiseMatrix};
use std::collections::HashMap;

/// 比较两个 SingleCellData 对象
///
/// # 参数
/// - `original`: 原始数据
/// - `converted`: 转换后的数据
/// - `tolerance`: 数值容差（默认 1e-7）
///
/// # 返回
/// - `ValidationReport`: 验证报告
pub fn compare_single_cell_data(
    original: &SingleCellData,
    converted: &SingleCellData,
    tolerance: f64,
) -> super::ValidationReport {
    let mut results = HashMap::new();
    
    // 比较表达矩阵
    let expr_result = compare_expression_matrix(&original.expression, &converted.expression, tolerance);
    results.insert("expression_matrix".to_string(), expr_result);
    
    // 比较细胞元数据
    let cell_meta_result = compare_dataframe(&original.cell_metadata, &converted.cell_metadata, tolerance);
    results.insert("cell_metadata".to_string(), cell_meta_result);
    
    // 比较基因元数据
    let gene_meta_result = compare_dataframe(&original.gene_metadata, &converted.gene_metadata, tolerance);
    results.insert("gene_metadata".to_string(), gene_meta_result);
    
    // 比较嵌入（如果存在）
    if let (Some(orig_emb), Some(conv_emb)) = (&original.embeddings, &converted.embeddings) {
        let emb_result = compare_embeddings(orig_emb, conv_emb, tolerance);
        results.insert("embeddings".to_string(), emb_result);
    } else if original.embeddings.is_some() != converted.embeddings.is_some() {
        results.insert("embeddings".to_string(), super::ComparisonResult {
            passed: false,
            message: "Embeddings presence mismatch".to_string(),
            max_difference: None,
        });
    }
    
    // 比较 layers（如果存在）
    if let (Some(orig_layers), Some(conv_layers)) = (&original.layers, &converted.layers) {
        let layers_result = compare_layers(orig_layers, conv_layers, tolerance);
        results.insert("layers".to_string(), layers_result);
    } else if original.layers.is_some() != converted.layers.is_some() {
        results.insert("layers".to_string(), super::ComparisonResult {
            passed: false,
            message: "Layers presence mismatch".to_string(),
            max_difference: None,
        });
    }
    
    // 比较 cell_pairwise（如果存在）
    if let (Some(orig_pw), Some(conv_pw)) = (&original.cell_pairwise, &converted.cell_pairwise) {
        let pw_result = compare_pairwise_matrices(orig_pw, conv_pw, tolerance);
        results.insert("cell_pairwise".to_string(), pw_result);
    } else if original.cell_pairwise.is_some() != converted.cell_pairwise.is_some() {
        results.insert("cell_pairwise".to_string(), super::ComparisonResult {
            passed: false,
            message: "Cell pairwise matrices presence mismatch".to_string(),
            max_difference: None,
        });
    }
    
    // 比较 gene_pairwise（如果存在）
    if let (Some(orig_pw), Some(conv_pw)) = (&original.gene_pairwise, &converted.gene_pairwise) {
        let pw_result = compare_pairwise_matrices(orig_pw, conv_pw, tolerance);
        results.insert("gene_pairwise".to_string(), pw_result);
    } else if original.gene_pairwise.is_some() != converted.gene_pairwise.is_some() {
        results.insert("gene_pairwise".to_string(), super::ComparisonResult {
            passed: false,
            message: "Gene pairwise matrices presence mismatch".to_string(),
            max_difference: None,
        });
    }
    
    super::ValidationReport { results }
}

/// 比较表达矩阵
///
/// # 参数
/// - `original`: 原始矩阵
/// - `converted`: 转换后的矩阵
/// - `tolerance`: 数值容差
///
/// # 返回
/// - `ComparisonResult`: 比较结果
pub fn compare_expression_matrix(
    original: &ExpressionMatrix,
    converted: &ExpressionMatrix,
    tolerance: f64,
) -> super::ComparisonResult {
    // 比较维度
    let (orig_rows, orig_cols) = original.shape();
    let (conv_rows, conv_cols) = converted.shape();
    
    if orig_rows != conv_rows || orig_cols != conv_cols {
        return super::ComparisonResult {
            passed: false,
            message: format!(
                "Dimension mismatch: original ({}, {}), converted ({}, {})",
                orig_rows, orig_cols, conv_rows, conv_cols
            ),
            max_difference: None,
        };
    }
    
    // 根据矩阵类型进行比较
    match (original, converted) {
        (ExpressionMatrix::SparseCSR(orig), ExpressionMatrix::SparseCSR(conv)) => {
            compare_sparse_csr(orig, conv, tolerance)
        }
        (ExpressionMatrix::SparseCSC(orig), ExpressionMatrix::SparseCSC(conv)) => {
            compare_sparse_csc(orig, conv, tolerance)
        }
        (ExpressionMatrix::Dense(orig), ExpressionMatrix::Dense(conv)) => {
            compare_dense(orig, conv, tolerance)
        }
        _ => super::ComparisonResult {
            passed: false,
            message: "Matrix type mismatch".to_string(),
            max_difference: None,
        },
    }
}

/// 比较稀疏 CSR 矩阵
fn compare_sparse_csr(
    original: &crate::ir::SparseMatrixCSR,
    converted: &crate::ir::SparseMatrixCSR,
    tolerance: f64,
) -> super::ComparisonResult {
    // 比较非零元素数量
    if original.data.len() != converted.data.len() {
        return super::ComparisonResult {
            passed: false,
            message: format!(
                "Non-zero element count mismatch: {} vs {}",
                original.data.len(),
                converted.data.len()
            ),
            max_difference: None,
        };
    }
    
    // 比较索引
    if original.indices != converted.indices {
        return super::ComparisonResult {
            passed: false,
            message: "Indices mismatch".to_string(),
            max_difference: None,
        };
    }
    
    if original.indptr != converted.indptr {
        return super::ComparisonResult {
            passed: false,
            message: "Indptr mismatch".to_string(),
            max_difference: None,
        };
    }
    
    // 比较数值
    let max_diff = original
        .data
        .iter()
        .zip(converted.data.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    
    let passed = max_diff < tolerance;
    
    super::ComparisonResult {
        passed,
        message: if passed {
            format!("Sparse CSR matrices match (max diff: {:.2e})", max_diff)
        } else {
            format!("Sparse CSR matrices differ (max diff: {:.2e} > tolerance {:.2e})", max_diff, tolerance)
        },
        max_difference: Some(max_diff),
    }
}

/// 比较稀疏 CSC 矩阵
fn compare_sparse_csc(
    original: &crate::ir::SparseMatrixCSC,
    converted: &crate::ir::SparseMatrixCSC,
    tolerance: f64,
) -> super::ComparisonResult {
    // 比较非零元素数量
    if original.data.len() != converted.data.len() {
        return super::ComparisonResult {
            passed: false,
            message: format!(
                "Non-zero element count mismatch: {} vs {}",
                original.data.len(),
                converted.data.len()
            ),
            max_difference: None,
        };
    }
    
    // 比较索引
    if original.indices != converted.indices {
        return super::ComparisonResult {
            passed: false,
            message: "Indices mismatch".to_string(),
            max_difference: None,
        };
    }
    
    if original.indptr != converted.indptr {
        return super::ComparisonResult {
            passed: false,
            message: "Indptr mismatch".to_string(),
            max_difference: None,
        };
    }
    
    // 比较数值
    let max_diff = original
        .data
        .iter()
        .zip(converted.data.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    
    let passed = max_diff < tolerance;
    
    super::ComparisonResult {
        passed,
        message: if passed {
            format!("Sparse CSC matrices match (max diff: {:.2e})", max_diff)
        } else {
            format!("Sparse CSC matrices differ (max diff: {:.2e} > tolerance {:.2e})", max_diff, tolerance)
        },
        max_difference: Some(max_diff),
    }
}

/// 比较稠密矩阵
fn compare_dense(
    original: &crate::ir::DenseMatrix,
    converted: &crate::ir::DenseMatrix,
    tolerance: f64,
) -> super::ComparisonResult {
    // 比较数据长度
    if original.data.len() != converted.data.len() {
        return super::ComparisonResult {
            passed: false,
            message: format!(
                "Data length mismatch: {} vs {}",
                original.data.len(),
                converted.data.len()
            ),
            max_difference: None,
        };
    }
    
    // 比较数值
    let max_diff = original
        .data
        .iter()
        .zip(converted.data.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    
    let passed = max_diff < tolerance;
    
    super::ComparisonResult {
        passed,
        message: if passed {
            format!("Dense matrices match (max diff: {:.2e})", max_diff)
        } else {
            format!("Dense matrices differ (max diff: {:.2e} > tolerance {:.2e})", max_diff, tolerance)
        },
        max_difference: Some(max_diff),
    }
}

/// 比较 DataFrame
pub fn compare_dataframe(
    original: &DataFrame,
    converted: &DataFrame,
    tolerance: f64,
) -> super::ComparisonResult {
    // 比较列数
    if original.n_cols() != converted.n_cols() {
        return super::ComparisonResult {
            passed: false,
            message: format!(
                "Column count mismatch: {} vs {}",
                original.n_cols(),
                converted.n_cols()
            ),
            max_difference: None,
        };
    }
    
    // 比较行数
    if original.n_rows != converted.n_rows {
        return super::ComparisonResult {
            passed: false,
            message: format!(
                "Row count mismatch: {} vs {}",
                original.n_rows,
                converted.n_rows
            ),
            max_difference: None,
        };
    }
    
    // 比较列名（排序后比较，因为顺序可能不同）
    let mut orig_cols = original.columns.clone();
    let mut conv_cols = converted.columns.clone();
    orig_cols.sort();
    conv_cols.sort();
    
    if orig_cols != conv_cols {
        return super::ComparisonResult {
            passed: false,
            message: format!(
                "Column names mismatch: {:?} vs {:?}",
                orig_cols, conv_cols
            ),
            max_difference: None,
        };
    }
    
    // 比较每列的数据类型和值
    use arrow::datatypes::DataType;
    use arrow::array::*;
    
    let mut max_diff_overall: f64 = 0.0;
    
    for col_name in &original.columns {
        let orig_array = original.column(col_name).unwrap();
        let conv_array = converted.column(col_name).unwrap();
        
        // 比较数据类型
        if orig_array.data_type() != conv_array.data_type() {
            return super::ComparisonResult {
                passed: false,
                message: format!(
                    "Column '{}' type mismatch: {:?} vs {:?}",
                    col_name,
                    orig_array.data_type(),
                    conv_array.data_type()
                ),
                max_difference: None,
            };
        }
        
        // 根据数据类型比较值
        match orig_array.data_type() {
            DataType::Int64 => {
                let orig = orig_array.as_any().downcast_ref::<Int64Array>().unwrap();
                let conv = conv_array.as_any().downcast_ref::<Int64Array>().unwrap();
                
                for i in 0..orig.len() {
                    match (orig.is_null(i), conv.is_null(i)) {
                        (true, true) => continue,
                        (true, false) | (false, true) => {
                            return super::ComparisonResult {
                                passed: false,
                                message: format!(
                                    "Column '{}' null mismatch at row {}",
                                    col_name, i
                                ),
                                max_difference: None,
                            };
                        }
                        (false, false) => {
                            if orig.value(i) != conv.value(i) {
                                return super::ComparisonResult {
                                    passed: false,
                                    message: format!(
                                        "Column '{}' value mismatch at row {}: {} vs {}",
                                        col_name, i, orig.value(i), conv.value(i)
                                    ),
                                    max_difference: None,
                                };
                            }
                        }
                    }
                }
            }
            DataType::Float64 => {
                let orig = orig_array.as_any().downcast_ref::<Float64Array>().unwrap();
                let conv = conv_array.as_any().downcast_ref::<Float64Array>().unwrap();
                
                for i in 0..orig.len() {
                    match (orig.is_null(i), conv.is_null(i)) {
                        (true, true) => continue,
                        (true, false) | (false, true) => {
                            return super::ComparisonResult {
                                passed: false,
                                message: format!(
                                    "Column '{}' null mismatch at row {}",
                                    col_name, i
                                ),
                                max_difference: None,
                            };
                        }
                        (false, false) => {
                            let diff = (orig.value(i) - conv.value(i)).abs();
                            max_diff_overall = max_diff_overall.max(diff);
                            
                            if diff >= tolerance {
                                return super::ComparisonResult {
                                    passed: false,
                                    message: format!(
                                        "Column '{}' value mismatch at row {}: {} vs {} (diff: {:.2e} > tolerance {:.2e})",
                                        col_name, i, orig.value(i), conv.value(i), diff, tolerance
                                    ),
                                    max_difference: Some(diff),
                                };
                            }
                        }
                    }
                }
            }
            DataType::Utf8 => {
                let orig = orig_array.as_any().downcast_ref::<StringArray>().unwrap();
                let conv = conv_array.as_any().downcast_ref::<StringArray>().unwrap();
                
                for i in 0..orig.len() {
                    match (orig.is_null(i), conv.is_null(i)) {
                        (true, true) => continue,
                        (true, false) | (false, true) => {
                            return super::ComparisonResult {
                                passed: false,
                                message: format!(
                                    "Column '{}' null mismatch at row {}",
                                    col_name, i
                                ),
                                max_difference: None,
                            };
                        }
                        (false, false) => {
                            if orig.value(i) != conv.value(i) {
                                return super::ComparisonResult {
                                    passed: false,
                                    message: format!(
                                        "Column '{}' value mismatch at row {}: '{}' vs '{}'",
                                        col_name, i, orig.value(i), conv.value(i)
                                    ),
                                    max_difference: None,
                                };
                            }
                        }
                    }
                }
            }
            DataType::Boolean => {
                let orig = orig_array.as_any().downcast_ref::<BooleanArray>().unwrap();
                let conv = conv_array.as_any().downcast_ref::<BooleanArray>().unwrap();
                
                for i in 0..orig.len() {
                    match (orig.is_null(i), conv.is_null(i)) {
                        (true, true) => continue,
                        (true, false) | (false, true) => {
                            return super::ComparisonResult {
                                passed: false,
                                message: format!(
                                    "Column '{}' null mismatch at row {}",
                                    col_name, i
                                ),
                                max_difference: None,
                            };
                        }
                        (false, false) => {
                            if orig.value(i) != conv.value(i) {
                                return super::ComparisonResult {
                                    passed: false,
                                    message: format!(
                                        "Column '{}' value mismatch at row {}: {} vs {}",
                                        col_name, i, orig.value(i), conv.value(i)
                                    ),
                                    max_difference: None,
                                };
                            }
                        }
                    }
                }
            }
            DataType::Dictionary(_, _) => {
                // Dictionary (Categorical) 类型比较
                // 比较编码和类别标签
                let orig = orig_array.as_any().downcast_ref::<DictionaryArray<arrow::datatypes::Int32Type>>().unwrap();
                let conv = conv_array.as_any().downcast_ref::<DictionaryArray<arrow::datatypes::Int32Type>>().unwrap();
                
                // 比较 keys (编码)
                for i in 0..orig.len() {
                    match (orig.is_null(i), conv.is_null(i)) {
                        (true, true) => continue,
                        (true, false) | (false, true) => {
                            return super::ComparisonResult {
                                passed: false,
                                message: format!(
                                    "Column '{}' null mismatch at row {}",
                                    col_name, i
                                ),
                                max_difference: None,
                            };
                        }
                        (false, false) => {
                            if orig.keys().value(i) != conv.keys().value(i) {
                                return super::ComparisonResult {
                                    passed: false,
                                    message: format!(
                                        "Column '{}' categorical code mismatch at row {}: {} vs {}",
                                        col_name, i, orig.keys().value(i), conv.keys().value(i)
                                    ),
                                    max_difference: None,
                                };
                            }
                        }
                    }
                }
                
                // 比较 values (类别标签)
                let orig_values = orig.values().as_any().downcast_ref::<StringArray>().unwrap();
                let conv_values = conv.values().as_any().downcast_ref::<StringArray>().unwrap();
                
                if orig_values.len() != conv_values.len() {
                    return super::ComparisonResult {
                        passed: false,
                        message: format!(
                            "Column '{}' categorical levels count mismatch: {} vs {}",
                            col_name, orig_values.len(), conv_values.len()
                        ),
                        max_difference: None,
                    };
                }
                
                for i in 0..orig_values.len() {
                    if orig_values.value(i) != conv_values.value(i) {
                        return super::ComparisonResult {
                            passed: false,
                            message: format!(
                                "Column '{}' categorical level mismatch at index {}: '{}' vs '{}'",
                                col_name, i, orig_values.value(i), conv_values.value(i)
                            ),
                            max_difference: None,
                        };
                    }
                }
            }
            _ => {
                // 其他类型暂不支持详细比较
                // 只检查数据类型匹配即可
            }
        }
    }
    
    super::ComparisonResult {
        passed: true,
        message: if max_diff_overall > 0.0 {
            format!("DataFrames match (max numeric diff: {:.2e})", max_diff_overall)
        } else {
            "DataFrames match (exact)".to_string()
        },
        max_difference: if max_diff_overall > 0.0 { Some(max_diff_overall) } else { None },
    }
}

/// 比较嵌入
pub fn compare_embeddings(
    original: &HashMap<String, Embedding>,
    converted: &HashMap<String, Embedding>,
    tolerance: f64,
) -> super::ComparisonResult {
    // 比较键集合
    let orig_keys: std::collections::HashSet<_> = original.keys().collect();
    let conv_keys: std::collections::HashSet<_> = converted.keys().collect();
    
    if orig_keys != conv_keys {
        return super::ComparisonResult {
            passed: false,
            message: format!(
                "Embedding keys mismatch: {:?} vs {:?}",
                orig_keys, conv_keys
            ),
            max_difference: None,
        };
    }
    
    // 比较每个嵌入
    let mut max_diff_overall: f64 = 0.0;
    for (key, orig_emb) in original {
        let conv_emb = &converted[key];
        
        // 比较维度
        if orig_emb.n_rows != conv_emb.n_rows || orig_emb.n_cols != conv_emb.n_cols {
            return super::ComparisonResult {
                passed: false,
                message: format!(
                    "Embedding '{}' dimension mismatch: ({}, {}) vs ({}, {})",
                    key, orig_emb.n_rows, orig_emb.n_cols, conv_emb.n_rows, conv_emb.n_cols
                ),
                max_difference: None,
            };
        }
        
        // 比较数值
        let max_diff = orig_emb
            .data
            .iter()
            .zip(conv_emb.data.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        
        max_diff_overall = max_diff_overall.max(max_diff);
        
        if max_diff >= tolerance {
            return super::ComparisonResult {
                passed: false,
                message: format!(
                    "Embedding '{}' differs (max diff: {:.2e} > tolerance {:.2e})",
                    key, max_diff, tolerance
                ),
                max_difference: Some(max_diff),
            };
        }
    }
    
    super::ComparisonResult {
        passed: true,
        message: format!("All embeddings match (max diff: {:.2e})", max_diff_overall),
        max_difference: Some(max_diff_overall),
    }
}

/// 比较 layers
pub fn compare_layers(
    original: &HashMap<String, ExpressionMatrix>,
    converted: &HashMap<String, ExpressionMatrix>,
    tolerance: f64,
) -> super::ComparisonResult {
    // 比较键集合
    let orig_keys: std::collections::HashSet<_> = original.keys().collect();
    let conv_keys: std::collections::HashSet<_> = converted.keys().collect();
    
    if orig_keys != conv_keys {
        return super::ComparisonResult {
            passed: false,
            message: format!(
                "Layer keys mismatch: {:?} vs {:?}",
                orig_keys, conv_keys
            ),
            max_difference: None,
        };
    }
    
    // 比较每个 layer
    let mut max_diff_overall: f64 = 0.0;
    for (key, orig_layer) in original {
        let conv_layer = &converted[key];
        
        let result = compare_expression_matrix(orig_layer, conv_layer, tolerance);
        
        if !result.passed {
            return super::ComparisonResult {
                passed: false,
                message: format!("Layer '{}': {}", key, result.message),
                max_difference: result.max_difference,
            };
        }
        
        if let Some(diff) = result.max_difference {
            max_diff_overall = max_diff_overall.max(diff);
        }
    }
    
    super::ComparisonResult {
        passed: true,
        message: format!("All layers match (max diff: {:.2e})", max_diff_overall),
        max_difference: Some(max_diff_overall),
    }
}

/// 比较成对矩阵
pub fn compare_pairwise_matrices(
    original: &HashMap<String, PairwiseMatrix>,
    converted: &HashMap<String, PairwiseMatrix>,
    tolerance: f64,
) -> super::ComparisonResult {
    // 比较键集合
    let orig_keys: std::collections::HashSet<_> = original.keys().collect();
    let conv_keys: std::collections::HashSet<_> = converted.keys().collect();
    
    if orig_keys != conv_keys {
        return super::ComparisonResult {
            passed: false,
            message: format!(
                "Pairwise matrix keys mismatch: {:?} vs {:?}",
                orig_keys, conv_keys
            ),
            max_difference: None,
        };
    }
    
    // 比较每个成对矩阵
    let mut max_diff_overall: f64 = 0.0;
    for (key, orig_pw) in original {
        let conv_pw = &converted[key];
        
        let result = compare_expression_matrix(&orig_pw.matrix, &conv_pw.matrix, tolerance);
        
        if !result.passed {
            return super::ComparisonResult {
                passed: false,
                message: format!("Pairwise matrix '{}': {}", key, result.message),
                max_difference: result.max_difference,
            };
        }
        
        if let Some(diff) = result.max_difference {
            max_diff_overall = max_diff_overall.max(diff);
        }
    }
    
    super::ComparisonResult {
        passed: true,
        message: format!("All pairwise matrices match (max diff: {:.2e})", max_diff_overall),
        max_difference: Some(max_diff_overall),
    }
}
