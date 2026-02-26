//! Unit tests for R binding converters
//! 
//! Note: Tests that require R runtime (R! macro) cannot be run as regular Rust tests.
//! These tests verify the pure Rust logic only.

#[cfg(test)]
mod tests {
    use crosscell::ir::{DenseMatrix, SparseMatrixCSC, ExpressionMatrix, DataFrame, Embedding};
    use arrow::array::{Int64Array, Float64Array, StringArray, BooleanArray};
    use std::sync::Arc;

    #[test]
    fn test_dense_matrix_creation() {
        // Test that DenseMatrix can be created with row-major data
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let matrix = DenseMatrix::new(data, 2, 3).unwrap();
        assert_eq!(matrix.n_rows, 2);
        assert_eq!(matrix.n_cols, 3);
        assert_eq!(matrix.data.len(), 6);
    }

    #[test]
    fn test_sparse_csc_creation() {
        // Test CSC matrix creation (like dgCMatrix)
        // Matrix:
        // [1, 0, 2]
        // [0, 3, 0]
        let data = vec![1.0, 3.0, 2.0];
        let indices = vec![0, 1, 0]; // row indices
        let indptr = vec![0, 1, 2, 3]; // column pointers
        let matrix = SparseMatrixCSC::new(data, indices, indptr, 2, 3).unwrap();
        assert_eq!(matrix.n_rows, 2);
        assert_eq!(matrix.n_cols, 3);
        assert_eq!(matrix.data.len(), 3);
    }

    #[test]
    fn test_dataframe_creation() {
        // Test DataFrame creation with Arrow arrays
        let int_col = Arc::new(Int64Array::from(vec![1, 2, 3]));
        let float_col = Arc::new(Float64Array::from(vec![1.1, 2.2, 3.3]));
        let str_col = Arc::new(StringArray::from(vec!["a", "b", "c"]));
        
        let df = DataFrame::new(
            vec!["int".to_string(), "float".to_string(), "str".to_string()],
            vec![int_col, float_col, str_col],
            3
        ).unwrap();
        
        assert_eq!(df.n_rows, 3);
        assert_eq!(df.columns.len(), 3);
    }

    #[test]
    fn test_embedding_creation() {
        // Test Embedding creation
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let emb = Embedding::new("PCA".to_string(), data, 2, 3).unwrap();
        assert_eq!(emb.name, "PCA");
        assert_eq!(emb.n_rows, 2);
        assert_eq!(emb.n_cols, 3);
    }

    #[test]
    fn test_column_major_to_row_major_conversion() {
        // R matrices are column-major, we need row-major
        // Column-major [1,2,3,4,5,6] for 2x3 matrix means:
        // col1: [1,2], col2: [3,4], col3: [5,6]
        // So the matrix is:
        // [1, 3, 5]
        // [2, 4, 6]
        
        let col_major = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let n_rows = 2;
        let n_cols = 3;
        
        let mut row_major = vec![0.0; n_rows * n_cols];
        for row in 0..n_rows {
            for col in 0..n_cols {
                row_major[row * n_cols + col] = col_major[col * n_rows + row];
            }
        }
        
        // Row-major should be [1,3,5,2,4,6]
        assert_eq!(row_major, vec![1.0, 3.0, 5.0, 2.0, 4.0, 6.0]);
    }

    #[test]
    fn test_factor_code_conversion() {
        // R factors are 1-indexed, Arrow DictionaryArray is 0-indexed
        let r_codes = vec![1, 2, 1, 2, 3]; // R factor codes (1-indexed)
        let arrow_codes: Vec<i32> = r_codes.iter()
            .map(|&v| v - 1) // Convert to 0-indexed
            .collect();
        
        assert_eq!(arrow_codes, vec![0, 1, 0, 1, 2]);
    }

    #[test]
    fn test_na_handling() {
        // R uses special values for NA
        // Integer NA: i32::MIN
        // Double NA: NaN (but we treat NaN as NA)
        
        let r_int_na = i32::MIN;
        let is_na = r_int_na == i32::MIN;
        assert!(is_na);
        
        let r_double_na = f64::NAN;
        let is_na = r_double_na.is_nan();
        assert!(is_na);
    }
}
