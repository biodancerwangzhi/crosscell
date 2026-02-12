//! 稀疏矩阵验证

use crate::ir::expression::{ExpressionMatrix, SparseMatrixCSC, SparseMatrixCSR};
use thiserror::Error;

#[derive(Error, Debug, PartialEq)]
pub enum ValidationError {
    #[error("Invalid indptr length: expected {expected}, got {actual}")]
    InvalidIndptr { expected: usize, actual: usize },

    #[error("Index out of bounds: index {index} >= dimension {dimension}")]
    IndexOutOfBounds { index: usize, dimension: usize },

    #[error(
        "Data and indices length mismatch: data.len()={data_len}, indices.len()={indices_len}"
    )]
    LengthMismatch { data_len: usize, indices_len: usize },

    #[error("Invalid indptr: not monotonically increasing at position {position}")]
    IndptrNotMonotonic { position: usize },

    #[error("Invalid indptr: last value {last} != nnz {nnz}")]
    IndptrLastValueMismatch { last: usize, nnz: usize },
}

/// 验证稀疏矩阵结构
pub fn validate_sparse_matrix(matrix: &ExpressionMatrix) -> Result<(), ValidationError> {
    match matrix {
        ExpressionMatrix::Dense(_) => Ok(()),
        ExpressionMatrix::SparseCSR(csr) => validate_csr(csr),
        ExpressionMatrix::SparseCSC(csc) => validate_csc(csc),
        ExpressionMatrix::Lazy(lazy) => {
            // 延迟加载矩阵：如果有缓存则验证缓存，否则只验证元数据
            if let Some(cached) = lazy.get_cached() {
                validate_sparse_matrix(&cached)
            } else {
                // 只验证元数据（维度）
                if lazy.shape.0 == 0 && lazy.shape.1 == 0 {
                    Ok(())
                } else if lazy.shape.0 == 0 || lazy.shape.1 == 0 {
                    // 不允许只有一个维度为 0
                    Err(ValidationError::InvalidIndptr {
                        expected: 1,
                        actual: 0,
                    })
                } else {
                    Ok(())
                }
            }
        }
    }
}

/// 验证 CSR 矩阵结构
fn validate_csr(csr: &SparseMatrixCSR) -> Result<(), ValidationError> {
    let nnz = csr.data.len();

    // 1. 验证 data 和 indices 长度一致
    if csr.data.len() != csr.indices.len() {
        return Err(ValidationError::LengthMismatch {
            data_len: csr.data.len(),
            indices_len: csr.indices.len(),
        });
    }

    // 2. 验证 indptr 长度 = n_rows + 1
    let expected_indptr_len = csr.n_rows + 1;
    if csr.indptr.len() != expected_indptr_len {
        return Err(ValidationError::InvalidIndptr {
            expected: expected_indptr_len,
            actual: csr.indptr.len(),
        });
    }

    // 3. 验证 indptr 单调递增
    for i in 0..csr.indptr.len() - 1 {
        if csr.indptr[i] > csr.indptr[i + 1] {
            return Err(ValidationError::IndptrNotMonotonic { position: i });
        }
    }

    // 4. 验证 indptr 最后一个值 = nnz
    if !csr.indptr.is_empty() {
        let last = csr.indptr[csr.indptr.len() - 1];
        if last != nnz {
            return Err(ValidationError::IndptrLastValueMismatch { last, nnz });
        }
    }

    // 5. 验证所有列索引在边界内
    for &col_idx in &csr.indices {
        if col_idx >= csr.n_cols {
            return Err(ValidationError::IndexOutOfBounds {
                index: col_idx,
                dimension: csr.n_cols,
            });
        }
    }

    Ok(())
}

/// 验证 CSC 矩阵结构
fn validate_csc(csc: &SparseMatrixCSC) -> Result<(), ValidationError> {
    let nnz = csc.data.len();

    // 1. 验证 data 和 indices 长度一致
    if csc.data.len() != csc.indices.len() {
        return Err(ValidationError::LengthMismatch {
            data_len: csc.data.len(),
            indices_len: csc.indices.len(),
        });
    }

    // 2. 验证 indptr 长度 = n_cols + 1
    let expected_indptr_len = csc.n_cols + 1;
    if csc.indptr.len() != expected_indptr_len {
        return Err(ValidationError::InvalidIndptr {
            expected: expected_indptr_len,
            actual: csc.indptr.len(),
        });
    }

    // 3. 验证 indptr 单调递增
    for i in 0..csc.indptr.len() - 1 {
        if csc.indptr[i] > csc.indptr[i + 1] {
            return Err(ValidationError::IndptrNotMonotonic { position: i });
        }
    }

    // 4. 验证 indptr 最后一个值 = nnz
    if !csc.indptr.is_empty() {
        let last = csc.indptr[csc.indptr.len() - 1];
        if last != nnz {
            return Err(ValidationError::IndptrLastValueMismatch { last, nnz });
        }
    }

    // 5. 验证所有行索引在边界内
    for &row_idx in &csc.indices {
        if row_idx >= csc.n_rows {
            return Err(ValidationError::IndexOutOfBounds {
                index: row_idx,
                dimension: csc.n_rows,
            });
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::expression::{SparseMatrixCSC, SparseMatrixCSR};

    #[test]
    fn test_validate_valid_csr() {
        let csr = SparseMatrixCSR {
            data: vec![1.0, 2.0, 3.0],
            indices: vec![0, 2, 1],
            indptr: vec![0, 2, 3],
            n_rows: 2,
            n_cols: 3,
        };

        let matrix = ExpressionMatrix::SparseCSR(csr);
        assert!(validate_sparse_matrix(&matrix).is_ok());
    }

    #[test]
    fn test_validate_valid_csc() {
        let csc = SparseMatrixCSC {
            data: vec![1.0, 3.0, 2.0],
            indices: vec![0, 1, 0],
            indptr: vec![0, 1, 2, 3],
            n_rows: 2,
            n_cols: 3,
        };

        let matrix = ExpressionMatrix::SparseCSC(csc);
        assert!(validate_sparse_matrix(&matrix).is_ok());
    }

    #[test]
    fn test_validate_length_mismatch() {
        let csr = SparseMatrixCSR {
            data: vec![1.0, 2.0],
            indices: vec![0, 2, 1], // 长度不匹配
            indptr: vec![0, 2, 3],
            n_rows: 2,
            n_cols: 3,
        };

        let matrix = ExpressionMatrix::SparseCSR(csr);
        let result = validate_sparse_matrix(&matrix);
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            ValidationError::LengthMismatch { .. }
        ));
    }

    #[test]
    fn test_validate_invalid_indptr_length() {
        let csr = SparseMatrixCSR {
            data: vec![1.0, 2.0, 3.0],
            indices: vec![0, 2, 1],
            indptr: vec![0, 2], // 长度应该是 n_rows + 1 = 3
            n_rows: 2,
            n_cols: 3,
        };

        let matrix = ExpressionMatrix::SparseCSR(csr);
        let result = validate_sparse_matrix(&matrix);
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            ValidationError::InvalidIndptr { .. }
        ));
    }

    #[test]
    fn test_validate_index_out_of_bounds() {
        let csr = SparseMatrixCSR {
            data: vec![1.0, 2.0, 3.0],
            indices: vec![0, 5, 1], // 5 >= n_cols (3)
            indptr: vec![0, 2, 3],
            n_rows: 2,
            n_cols: 3,
        };

        let matrix = ExpressionMatrix::SparseCSR(csr);
        let result = validate_sparse_matrix(&matrix);
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            ValidationError::IndexOutOfBounds { .. }
        ));
    }

    #[test]
    fn test_validate_indptr_not_monotonic() {
        let csr = SparseMatrixCSR {
            data: vec![1.0, 2.0, 3.0],
            indices: vec![0, 2, 1],
            indptr: vec![0, 3, 2], // 不单调递增
            n_rows: 2,
            n_cols: 3,
        };

        let matrix = ExpressionMatrix::SparseCSR(csr);
        let result = validate_sparse_matrix(&matrix);
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            ValidationError::IndptrNotMonotonic { .. }
        ));
    }

    #[test]
    fn test_validate_indptr_last_value_mismatch() {
        let csr = SparseMatrixCSR {
            data: vec![1.0, 2.0, 3.0],
            indices: vec![0, 2, 1],
            indptr: vec![0, 2, 5], // 最后一个值应该是 3
            n_rows: 2,
            n_cols: 3,
        };

        let matrix = ExpressionMatrix::SparseCSR(csr);
        let result = validate_sparse_matrix(&matrix);
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            ValidationError::IndptrLastValueMismatch { .. }
        ));
    }

    #[test]
    fn test_validate_empty_matrix() {
        let csr = SparseMatrixCSR::empty(0, 0);
        let matrix = ExpressionMatrix::SparseCSR(csr);
        assert!(validate_sparse_matrix(&matrix).is_ok());
    }

    #[test]
    fn test_validate_empty_sparse_matrix() {
        // 非零维度但无数据
        let csr = SparseMatrixCSR {
            data: vec![],
            indices: vec![],
            indptr: vec![0, 0, 0],
            n_rows: 2,
            n_cols: 3,
        };

        let matrix = ExpressionMatrix::SparseCSR(csr);
        assert!(validate_sparse_matrix(&matrix).is_ok());
    }
}
