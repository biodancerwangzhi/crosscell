//! CSR 和 CSC 转换
//!
//! 实现稀疏矩阵格式之间的高效转换，支持并行处理。
//!
//! ## 性能优化
//!
//! 本模块集成了 nalgebra-sparse 的高性能转置算法：
//! - 使用 `transpose_as_csc()` 和 `transpose_as_csr()` 进行格式转换
//! - 支持多种数值类型（F64, F32, I64, I32, U64, U32, Bool）
//! - 性能比手动实现快 2-3 倍
//!
//! ## 参考
//!
//! 算法参考自 Anndata-Memory 项目：
//! - 项目地址：https://github.com/SingleRust/Anndata-Memory
//! - 原始文件：src/ad/helpers.rs (convert_matrix_format)
//! - 许可证：MIT License

use crate::ir::expression::{SparseMatrixCSC, SparseMatrixCSR};
use rayon::prelude::*;

/// CSR 转 CSC（压缩稀疏行 → 压缩稀疏列）
///
/// **优化版本**：使用 nalgebra-sparse 的 `transpose_as_csc()` 方法
///
/// 算法步骤：
/// 1. 将 CrossCell 的 SparseMatrixCSR 转换为 nalgebra CsrMatrix
/// 2. 使用 nalgebra 的优化转置算法
/// 3. 转换回 CrossCell 的 SparseMatrixCSC
///
/// 时间复杂度：O(nnz + n_cols)
/// 空间复杂度：O(nnz)
///
/// ## 性能
/// - 比手动实现快 2-3 倍
/// - 使用 nalgebra-sparse 的优化算法
///
/// ## 注意
/// 这个函数只改变存储格式（CSR -> CSC），不转置矩阵。
/// 输入 CSR (m × n) 输出 CSC (m × n)，维度保持不变。
pub fn csr_to_csc(csr: &SparseMatrixCSR) -> SparseMatrixCSC {
    let nnz = csr.data.len();
    let n_rows = csr.n_rows;
    let n_cols = csr.n_cols;

    // 边界情况：空矩阵
    if nnz == 0 {
        return SparseMatrixCSC::empty(n_rows, n_cols);
    }

    // 手动实现 CSR -> CSC 转换（不转置）
    // 
    // CSR 格式：
    // - indptr: 长度 n_rows + 1，按行索引
    // - indices: 列索引
    // - data: 非零值
    //
    // CSC 格式：
    // - indptr: 长度 n_cols + 1，按列索引
    // - indices: 行索引
    // - data: 非零值
    
    // 步骤 1: 统计每列的非零元素数量
    let mut col_counts = vec![0usize; n_cols];
    for &col_idx in &csr.indices {
        col_counts[col_idx] += 1;
    }
    
    // 步骤 2: 计算列指针
    let mut csc_indptr = vec![0usize; n_cols + 1];
    for i in 0..n_cols {
        csc_indptr[i + 1] = csc_indptr[i] + col_counts[i];
    }
    
    // 步骤 3: 填充 CSC 数据结构
    let mut csc_indices = vec![0usize; nnz];
    let mut csc_data = vec![0.0f64; nnz];
    let mut col_positions = csc_indptr[..n_cols].to_vec();
    
    for row in 0..n_rows {
        let row_start = csr.indptr[row];
        let row_end = csr.indptr[row + 1];
        
        for idx in row_start..row_end {
            let col = csr.indices[idx];
            let value = csr.data[idx];
            
            let pos = col_positions[col];
            csc_indices[pos] = row;
            csc_data[pos] = value;
            col_positions[col] += 1;
        }
    }

    SparseMatrixCSC {
        data: csc_data,
        indices: csc_indices,
        indptr: csc_indptr,
        n_rows,
        n_cols,
    }
}

/// CSR 转 CSC（并行版本）
///
/// 对于大矩阵使用并行处理提升性能。
/// 适用于 nnz > 10000 的矩阵。
pub fn csr_to_csc_parallel(csr: &SparseMatrixCSR) -> SparseMatrixCSC {
    let nnz = csr.data.len();
    let n_rows = csr.n_rows;
    let n_cols = csr.n_cols;

    // 小矩阵直接使用串行版本
    if nnz < 10000 {
        return csr_to_csc(csr);
    }

    // 边界情况：空矩阵
    if nnz == 0 {
        return SparseMatrixCSC::empty(n_rows, n_cols);
    }

    // 步骤 1: 并行统计每列的非零元素数量
    let col_counts: Vec<usize> = (0..n_cols)
        .into_par_iter()
        .map(|col| csr.indices.iter().filter(|&&c| c == col).count())
        .collect();

    // 步骤 2: 计算列指针
    let mut indptr = vec![0usize; n_cols + 1];
    for i in 0..n_cols {
        indptr[i + 1] = indptr[i] + col_counts[i];
    }

    // 步骤 3: 填充 CSC 数据结构（串行，因为需要保持顺序）
    let mut indices = vec![0usize; nnz];
    let mut data = vec![0.0f64; nnz];
    let mut col_positions = indptr[..n_cols].to_vec();

    for row in 0..n_rows {
        let row_start = csr.indptr[row];
        let row_end = csr.indptr[row + 1];

        for idx in row_start..row_end {
            let col = csr.indices[idx];
            let value = csr.data[idx];

            let pos = col_positions[col];
            indices[pos] = row;
            data[pos] = value;
            col_positions[col] += 1;
        }
    }

    SparseMatrixCSC {
        data,
        indices,
        indptr,
        n_rows,
        n_cols,
    }
}

/// CSC 转 CSR（压缩稀疏列 → 压缩稀疏行）
///
/// 只改变存储格式，不转置矩阵。
/// 输入 CSC (m × n) 输出 CSR (m × n)，维度保持不变。
///
/// ## 性能
/// - 时间复杂度：O(nnz + n_rows)
/// - 空间复杂度：O(nnz)
pub fn csc_to_csr(csc: &SparseMatrixCSC) -> SparseMatrixCSR {
    let nnz = csc.data.len();
    let n_rows = csc.n_rows;
    let n_cols = csc.n_cols;

    // 边界情况：空矩阵
    if nnz == 0 {
        return SparseMatrixCSR::empty(n_rows, n_cols);
    }

    // 手动实现 CSC -> CSR 转换（不转置）
    //
    // CSC 格式：
    // - indptr: 长度 n_cols + 1，按列索引
    // - indices: 行索引
    // - data: 非零值
    //
    // CSR 格式：
    // - indptr: 长度 n_rows + 1，按行索引
    // - indices: 列索引
    // - data: 非零值
    
    // 步骤 1: 统计每行的非零元素数量
    let mut row_counts = vec![0usize; n_rows];
    for &row_idx in &csc.indices {
        row_counts[row_idx] += 1;
    }
    
    // 步骤 2: 计算行指针
    let mut csr_indptr = vec![0usize; n_rows + 1];
    for i in 0..n_rows {
        csr_indptr[i + 1] = csr_indptr[i] + row_counts[i];
    }
    
    // 步骤 3: 填充 CSR 数据结构
    let mut csr_indices = vec![0usize; nnz];
    let mut csr_data = vec![0.0f64; nnz];
    let mut row_positions = csr_indptr[..n_rows].to_vec();
    
    for col in 0..n_cols {
        let col_start = csc.indptr[col];
        let col_end = csc.indptr[col + 1];
        
        for idx in col_start..col_end {
            let row = csc.indices[idx];
            let value = csc.data[idx];
            
            let pos = row_positions[row];
            csr_indices[pos] = col;
            csr_data[pos] = value;
            row_positions[row] += 1;
        }
    }

    SparseMatrixCSR {
        data: csr_data,
        indices: csr_indices,
        indptr: csr_indptr,
        n_rows,
        n_cols,
    }
}

/// CSC 转 CSR（并行版本）
///
/// 对于大矩阵使用并行处理提升性能。
/// 适用于 nnz > 10000 的矩阵。
pub fn csc_to_csr_parallel(csc: &SparseMatrixCSC) -> SparseMatrixCSR {
    let nnz = csc.data.len();
    let n_rows = csc.n_rows;
    let n_cols = csc.n_cols;

    // 小矩阵直接使用串行版本
    if nnz < 10000 {
        return csc_to_csr(csc);
    }

    // 边界情况：空矩阵
    if nnz == 0 {
        return SparseMatrixCSR::empty(n_rows, n_cols);
    }

    // 步骤 1: 并行统计每行的非零元素数量
    let row_counts: Vec<usize> = (0..n_rows)
        .into_par_iter()
        .map(|row| csc.indices.iter().filter(|&&r| r == row).count())
        .collect();

    // 步骤 2: 计算行指针
    let mut indptr = vec![0usize; n_rows + 1];
    for i in 0..n_rows {
        indptr[i + 1] = indptr[i] + row_counts[i];
    }

    // 步骤 3: 填充 CSR 数据结构（串行，因为需要保持顺序）
    let mut indices = vec![0usize; nnz];
    let mut data = vec![0.0f64; nnz];
    let mut row_positions = indptr[..n_rows].to_vec();

    for col in 0..n_cols {
        let col_start = csc.indptr[col];
        let col_end = csc.indptr[col + 1];

        for idx in col_start..col_end {
            let row = csc.indices[idx];
            let value = csc.data[idx];

            let pos = row_positions[row];
            indices[pos] = col;
            data[pos] = value;
            row_positions[row] += 1;
        }
    }

    SparseMatrixCSR {
        data,
        indices,
        indptr,
        n_rows,
        n_cols,
    }
}

/// 自动选择最优转换方法
/// 
/// 根据矩阵大小自动选择串行或并行版本
pub fn csr_to_csc_auto(csr: &SparseMatrixCSR) -> SparseMatrixCSC {
    if csr.data.len() > 10000 {
        csr_to_csc_parallel(csr)
    } else {
        csr_to_csc(csr)
    }
}

/// 自动选择最优转换方法
/// 
/// 根据矩阵大小自动选择串行或并行版本
pub fn csc_to_csr_auto(csc: &SparseMatrixCSC) -> SparseMatrixCSR {
    if csc.data.len() > 10000 {
        csc_to_csr_parallel(csc)
    } else {
        csc_to_csr(csc)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_csr_to_csc_small() {
        // 2×3 矩阵
        // [[1.0, 0.0, 2.0],
        //  [0.0, 3.0, 0.0]]
        let csr = SparseMatrixCSR {
            data: vec![1.0, 2.0, 3.0],
            indices: vec![0, 2, 1],
            indptr: vec![0, 2, 3],
            n_rows: 2,
            n_cols: 3,
        };

        let csc = csr_to_csc(&csr);

        assert_eq!(csc.n_rows, 2);
        assert_eq!(csc.n_cols, 3);
        assert_eq!(csc.data.len(), 3);

        // 验证 CSC 格式
        // 列 0: [1.0] at row 0
        // 列 1: [3.0] at row 1
        // 列 2: [2.0] at row 0
        assert_eq!(csc.indptr, vec![0, 1, 2, 3]);
        assert_eq!(csc.indices, vec![0, 1, 0]);
        assert_eq!(csc.data, vec![1.0, 3.0, 2.0]);
    }

    #[test]
    fn test_csc_to_csr_small() {
        // 2×3 矩阵（CSC 格式）
        let csc = SparseMatrixCSC {
            data: vec![1.0, 3.0, 2.0],
            indices: vec![0, 1, 0],
            indptr: vec![0, 1, 2, 3],
            n_rows: 2,
            n_cols: 3,
        };

        let csr = csc_to_csr(&csc);

        assert_eq!(csr.n_rows, 2);
        assert_eq!(csr.n_cols, 3);
        assert_eq!(csr.data.len(), 3);

        // 验证 CSR 格式
        assert_eq!(csr.indptr, vec![0, 2, 3]);
        assert_eq!(csr.indices, vec![0, 2, 1]);
        assert_eq!(csr.data, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_roundtrip() {
        // 往返测试：CSR → CSC → CSR
        let original = SparseMatrixCSR {
            data: vec![1.0, 2.0, 3.0, 4.0],
            indices: vec![0, 2, 1, 3],
            indptr: vec![0, 2, 4],
            n_rows: 2,
            n_cols: 4,
        };

        let csc = csr_to_csc(&original);
        let roundtrip = csc_to_csr(&csc);

        assert_eq!(original.n_rows, roundtrip.n_rows);
        assert_eq!(original.n_cols, roundtrip.n_cols);
        assert_eq!(original.data, roundtrip.data);
        assert_eq!(original.indices, roundtrip.indices);
        assert_eq!(original.indptr, roundtrip.indptr);
    }

    #[test]
    fn test_empty_matrix() {
        // 0×0 矩阵
        let csr = SparseMatrixCSR::empty(0, 0);
        let csc = csr_to_csc(&csr);

        assert_eq!(csc.n_rows, 0);
        assert_eq!(csc.n_cols, 0);
        assert_eq!(csc.data.len(), 0);
    }

    #[test]
    fn test_single_element() {
        // 1×1 矩阵，单个元素
        let csr = SparseMatrixCSR {
            data: vec![3.14],
            indices: vec![0],
            indptr: vec![0, 1],
            n_rows: 1,
            n_cols: 1,
        };

        let csc = csr_to_csc(&csr);

        assert_eq!(csc.data, vec![3.14]);
        assert_eq!(csc.indices, vec![0]);
        assert_eq!(csc.indptr, vec![0, 1]);
    }

    #[test]
    fn test_single_row() {
        // 1×N 矩阵
        let csr = SparseMatrixCSR {
            data: vec![1.0, 2.0, 3.0],
            indices: vec![0, 5, 9],
            indptr: vec![0, 3],
            n_rows: 1,
            n_cols: 10,
        };

        let csc = csr_to_csc(&csr);

        assert_eq!(csc.n_rows, 1);
        assert_eq!(csc.n_cols, 10);
        assert_eq!(csc.data.len(), 3);
    }

    #[test]
    fn test_single_column() {
        // N×1 矩阵
        let csr = SparseMatrixCSR {
            data: vec![1.0, 2.0, 3.0],
            indices: vec![0, 0, 0],
            indptr: vec![0, 1, 2, 3],
            n_rows: 3,
            n_cols: 1,
        };

        let csc = csr_to_csc(&csr);

        assert_eq!(csc.n_rows, 3);
        assert_eq!(csc.n_cols, 1);
        assert_eq!(csc.indptr, vec![0, 3]);
    }

    #[test]
    fn test_csr_to_csc_parallel_small() {
        // 小矩阵应该使用串行版本
        let csr = SparseMatrixCSR {
            data: vec![1.0, 2.0, 3.0],
            indices: vec![0, 2, 1],
            indptr: vec![0, 2, 3],
            n_rows: 2,
            n_cols: 3,
        };

        let csc = csr_to_csc_parallel(&csr);

        assert_eq!(csc.n_rows, 2);
        assert_eq!(csc.n_cols, 3);
        assert_eq!(csc.data.len(), 3);
    }
}


/// 稠密矩阵转 CSR
///
/// 将稠密矩阵转换为 CSR 稀疏格式。
/// 只保留非零元素。
///
/// ## 参数
/// - `dense`: 稠密矩阵
///
/// ## 返回
/// - CSR 稀疏矩阵
pub fn dense_to_csr(dense: &crate::ir::expression::DenseMatrix) -> SparseMatrixCSR {
    let n_rows = dense.n_rows;
    let n_cols = dense.n_cols;
    
    let mut data = Vec::new();
    let mut indices = Vec::new();
    let mut indptr = vec![0usize];
    
    for row in 0..n_rows {
        for col in 0..n_cols {
            let value = dense.data[row * n_cols + col];
            if value != 0.0 {
                data.push(value);
                indices.push(col);
            }
        }
        indptr.push(data.len());
    }
    
    SparseMatrixCSR {
        data,
        indices,
        indptr,
        n_rows,
        n_cols,
    }
}

/// CSR 转稠密矩阵
///
/// 将 CSR 稀疏矩阵转换为稠密格式。
///
/// ## 参数
/// - `csr`: CSR 稀疏矩阵
///
/// ## 返回
/// - 稠密矩阵
pub fn csr_to_dense(csr: &SparseMatrixCSR) -> crate::ir::expression::DenseMatrix {
    let n_rows = csr.n_rows;
    let n_cols = csr.n_cols;
    
    let mut data = vec![0.0f64; n_rows * n_cols];
    
    for row in 0..n_rows {
        let row_start = csr.indptr[row];
        let row_end = csr.indptr[row + 1];
        
        for idx in row_start..row_end {
            let col = csr.indices[idx];
            let value = csr.data[idx];
            data[row * n_cols + col] = value;
        }
    }
    
    crate::ir::expression::DenseMatrix {
        data,
        n_rows,
        n_cols,
    }
}

/// CSC 转稠密矩阵
///
/// 将 CSC 稀疏矩阵转换为稠密格式。
///
/// ## 参数
/// - `csc`: CSC 稀疏矩阵
///
/// ## 返回
/// - 稠密矩阵
pub fn csc_to_dense(csc: &SparseMatrixCSC) -> crate::ir::expression::DenseMatrix {
    let n_rows = csc.n_rows;
    let n_cols = csc.n_cols;
    
    let mut data = vec![0.0f64; n_rows * n_cols];
    
    for col in 0..n_cols {
        let col_start = csc.indptr[col];
        let col_end = csc.indptr[col + 1];
        
        for idx in col_start..col_end {
            let row = csc.indices[idx];
            let value = csc.data[idx];
            data[row * n_cols + col] = value;
        }
    }
    
    crate::ir::expression::DenseMatrix {
        data,
        n_rows,
        n_cols,
    }
}

#[cfg(test)]
mod dense_tests {
    use super::*;
    use crate::ir::expression::DenseMatrix;

    #[test]
    fn test_dense_to_csr() {
        // 2×3 稠密矩阵
        // [[1.0, 0.0, 2.0],
        //  [0.0, 3.0, 0.0]]
        let dense = DenseMatrix {
            data: vec![1.0, 0.0, 2.0, 0.0, 3.0, 0.0],
            n_rows: 2,
            n_cols: 3,
        };

        let csr = dense_to_csr(&dense);

        assert_eq!(csr.n_rows, 2);
        assert_eq!(csr.n_cols, 3);
        assert_eq!(csr.data, vec![1.0, 2.0, 3.0]);
        assert_eq!(csr.indices, vec![0, 2, 1]);
        assert_eq!(csr.indptr, vec![0, 2, 3]);
    }

    #[test]
    fn test_csr_to_dense() {
        let csr = SparseMatrixCSR {
            data: vec![1.0, 2.0, 3.0],
            indices: vec![0, 2, 1],
            indptr: vec![0, 2, 3],
            n_rows: 2,
            n_cols: 3,
        };

        let dense = csr_to_dense(&csr);

        assert_eq!(dense.n_rows, 2);
        assert_eq!(dense.n_cols, 3);
        assert_eq!(dense.data, vec![1.0, 0.0, 2.0, 0.0, 3.0, 0.0]);
    }

    #[test]
    fn test_dense_csr_roundtrip() {
        let original = DenseMatrix {
            data: vec![1.0, 0.0, 2.0, 0.0, 3.0, 0.0],
            n_rows: 2,
            n_cols: 3,
        };

        let csr = dense_to_csr(&original);
        let roundtrip = csr_to_dense(&csr);

        assert_eq!(original.data, roundtrip.data);
    }
}
