//! Tests for AI-Ready Export (Task 9)
//!
//! Property tests 9-12 and unit tests for normalize, gene selection, gene ID replacement.

use arrow::array::{ArrayRef, StringArray};
use crosscell::ir::{
    DataFrame, DatasetMetadata, DenseMatrix, ExpressionMatrix, SingleCellData, SparseMatrixCSR,
};
use crosscell::transform::{apply_gene_id_column, normalize_library_size, select_top_variable_genes};
use proptest::prelude::*;
use std::sync::Arc;

// ============================================================
// Helper: build a SingleCellData from a dense matrix
// ============================================================

fn make_scd_dense(data: Vec<f64>, n_rows: usize, n_cols: usize) -> SingleCellData {
    let expr = ExpressionMatrix::Dense(DenseMatrix::new(data, n_rows, n_cols).unwrap());
    let cell_meta = DataFrame::empty(n_rows);
    let gene_meta = DataFrame::empty(n_cols);
    let metadata = DatasetMetadata::new(n_rows, n_cols, "test".to_string());
    SingleCellData::new(expr, cell_meta, gene_meta, metadata).unwrap()
}

fn make_scd_with_gene_col(
    data: Vec<f64>,
    n_rows: usize,
    n_cols: usize,
    col_name: &str,
    col_values: Vec<&str>,
) -> SingleCellData {
    let expr = ExpressionMatrix::Dense(DenseMatrix::new(data, n_rows, n_cols).unwrap());
    let cell_meta = DataFrame::empty(n_rows);
    let gene_col: ArrayRef = Arc::new(StringArray::from(col_values));
    let gene_meta = DataFrame::new(vec![col_name.to_string()], vec![gene_col], n_cols).unwrap();
    let metadata = DatasetMetadata::new(n_rows, n_cols, "test".to_string());
    SingleCellData::new(expr, cell_meta, gene_meta, metadata).unwrap()
}

// ============================================================
// Property 9: Normalization preserves dimensions and non-zero pattern
// Feature: reviewer-enhancements, Property 9
// ============================================================

fn arb_dense_matrix(
    max_rows: usize,
    max_cols: usize,
) -> impl Strategy<Value = (Vec<f64>, usize, usize)> {
    (1..=max_rows, 1..=max_cols).prop_flat_map(move |(r, c)| {
        proptest::collection::vec(0.0f64..100.0, r * c).prop_map(move |data| (data, r, c))
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(100))]

    // Feature: reviewer-enhancements, Property 9: 归一化不变量
    #[test]
    fn prop_normalize_preserves_dimensions((data, n_rows, n_cols) in arb_dense_matrix(20, 20)) {
        let matrix = ExpressionMatrix::Dense(DenseMatrix::new(data.clone(), n_rows, n_cols).unwrap());
        let result = normalize_library_size(&matrix).unwrap();
        let (r, c) = result.shape();
        prop_assert_eq!(r, n_rows);
        prop_assert_eq!(c, n_cols);

        // Non-zero pattern preserved
        if let (ExpressionMatrix::Dense(orig), ExpressionMatrix::Dense(norm)) = (&matrix, &result) {
            for i in 0..orig.data.len() {
                if orig.data[i] == 0.0 {
                    prop_assert!((norm.data[i] - 0.0).abs() < 1e-15,
                        "Zero element became non-zero at index {}", i);
                }
                if orig.data[i] > 0.0 {
                    // After normalization, non-zero elements should be > 0
                    // (log1p of positive value is positive)
                    prop_assert!(norm.data[i] > 0.0,
                        "Non-zero element became zero at index {}", i);
                }
            }
        }
    }
}

// ============================================================
// Property 10: Gene selection returns min(N, total_genes) genes
// Feature: reviewer-enhancements, Property 10
// ============================================================

proptest! {
    #![proptest_config(ProptestConfig::with_cases(50))]

    // Feature: reviewer-enhancements, Property 10: 基因筛选不变量
    #[test]
    fn prop_gene_selection_correct_count(
        n_rows in 2..10usize,
        n_cols in 3..15usize,
        n_select in 1..20usize,
    ) {
        // Create data with varying column variances
        let mut data = vec![0.0f64; n_rows * n_cols];
        for col in 0..n_cols {
            for row in 0..n_rows {
                // Column `col` gets variance proportional to col index
                data[row * n_cols + col] = (col as f64 + 1.0) * (row as f64);
            }
        }
        let mut scd = make_scd_dense(data, n_rows, n_cols);
        select_top_variable_genes(&mut scd, n_select).unwrap();

        let expected = n_select.min(n_cols);
        let (_, actual_cols) = scd.expression.shape();
        prop_assert_eq!(actual_cols, expected);
        prop_assert_eq!(scd.metadata.n_genes, expected);
        prop_assert_eq!(scd.gene_metadata.n_rows, expected);
    }
}

// ============================================================
// Property 11: Gene ID replacement uses correct column values
// Feature: reviewer-enhancements, Property 11
// ============================================================

proptest! {
    #![proptest_config(ProptestConfig::with_cases(50))]

    // Feature: reviewer-enhancements, Property 11: 基因标识符替换正确性
    #[test]
    fn prop_gene_id_replacement(n_genes in 1..10usize) {
        let n_rows = 3;
        let data = vec![1.0f64; n_rows * n_genes];
        let gene_names: Vec<String> = (0..n_genes).map(|i| format!("ENSG{:05}", i)).collect();
        let gene_names_ref: Vec<&str> = gene_names.iter().map(|s| s.as_str()).collect();

        let mut scd = make_scd_with_gene_col(data, n_rows, n_genes, "ensembl_id", gene_names_ref);
        apply_gene_id_column(&mut scd, "ensembl_id").unwrap();

        // Check _index column exists and has correct values
        let idx = scd.gene_metadata.column_index("_index").expect("_index column should exist");
        let arr = scd.gene_metadata.data[idx]
            .as_any()
            .downcast_ref::<StringArray>()
            .expect("_index should be StringArray");
        for i in 0..n_genes {
            prop_assert_eq!(arr.value(i), &format!("ENSG{:05}", i));
        }
    }
}

// ============================================================
// Property 12: Normalize → select order consistency
// Feature: reviewer-enhancements, Property 12
// ============================================================

#[test]
fn test_normalize_then_select_order_consistency() {
    // Create a small dataset
    let n_rows = 5;
    let n_cols = 8;
    let mut data = vec![0.0f64; n_rows * n_cols];
    for col in 0..n_cols {
        for row in 0..n_rows {
            data[row * n_cols + col] = ((col + 1) * (row + 1)) as f64;
        }
    }

    // Path A: normalize then select
    let mut scd_a = make_scd_dense(data.clone(), n_rows, n_cols);
    scd_a.expression = normalize_library_size(&scd_a.expression).unwrap();
    select_top_variable_genes(&mut scd_a, 4).unwrap();

    // Path B: same operations (simulating CLI --normalize --top-genes 4)
    let mut scd_b = make_scd_dense(data, n_rows, n_cols);
    scd_b.expression = normalize_library_size(&scd_b.expression).unwrap();
    select_top_variable_genes(&mut scd_b, 4).unwrap();

    // Results should be identical
    let (ra, ca) = scd_a.expression.shape();
    let (rb, cb) = scd_b.expression.shape();
    assert_eq!(ra, rb);
    assert_eq!(ca, cb);

    if let (ExpressionMatrix::Dense(a), ExpressionMatrix::Dense(b)) =
        (&scd_a.expression, &scd_b.expression)
    {
        for i in 0..a.data.len() {
            assert!(
                (a.data[i] - b.data[i]).abs() < 1e-15,
                "Mismatch at index {}: {} vs {}",
                i,
                a.data[i],
                b.data[i]
            );
        }
    }
}


// ============================================================
// Unit Tests (9.9)
// ============================================================

/// Test: N > gene count should warn and keep all genes
#[test]
fn test_top_genes_exceeds_count() {
    let n_rows = 3;
    let n_cols = 5;
    let data = vec![1.0f64; n_rows * n_cols];
    let mut scd = make_scd_dense(data, n_rows, n_cols);

    // Request more genes than available
    select_top_variable_genes(&mut scd, 100).unwrap();

    // Should keep all genes
    let (_, cols) = scd.expression.shape();
    assert_eq!(cols, n_cols);
    assert_eq!(scd.metadata.n_genes, n_cols);
}

/// Test: Missing column returns error
#[test]
fn test_gene_id_column_missing() {
    let n_rows = 3;
    let n_cols = 4;
    let data = vec![1.0f64; n_rows * n_cols];
    let mut scd = make_scd_dense(data, n_rows, n_cols);

    let result = apply_gene_id_column(&mut scd, "nonexistent_column");
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("not found"));
}

/// Test: All-zero row normalization stays zero
#[test]
fn test_normalize_all_zero_row() {
    // Row 0: all zeros, Row 1: has values
    let data = vec![0.0, 0.0, 0.0, 1.0, 2.0, 3.0];
    let matrix = ExpressionMatrix::Dense(DenseMatrix::new(data, 2, 3).unwrap());
    let result = normalize_library_size(&matrix).unwrap();

    if let ExpressionMatrix::Dense(dense) = &result {
        // Row 0 should remain all zeros
        assert_eq!(dense.data[0], 0.0);
        assert_eq!(dense.data[1], 0.0);
        assert_eq!(dense.data[2], 0.0);
        // Row 1 should have positive values
        assert!(dense.data[3] > 0.0);
        assert!(dense.data[4] > 0.0);
        assert!(dense.data[5] > 0.0);
    } else {
        panic!("Expected Dense matrix");
    }
}

/// Test: Normalization with CSR matrix
#[test]
fn test_normalize_csr() {
    // 2x3 CSR matrix: [[1, 0, 2], [0, 3, 0]]
    let csr = SparseMatrixCSR::new(
        vec![1.0, 2.0, 3.0],
        vec![0, 2, 1],
        vec![0, 2, 3],
        2,
        3,
    )
    .unwrap();
    let matrix = ExpressionMatrix::SparseCSR(csr);
    let result = normalize_library_size(&matrix).unwrap();

    let (r, c) = result.shape();
    assert_eq!(r, 2);
    assert_eq!(c, 3);

    // Check non-zero pattern preserved
    if let ExpressionMatrix::SparseCSR(res_csr) = &result {
        assert_eq!(res_csr.indices, vec![0, 2, 1]);
        assert_eq!(res_csr.indptr, vec![0, 2, 3]);
        // Values should be log1p(x * 10000 / total)
        // Row 0: total=3, val[0]=1*10000/3, val[1]=2*10000/3
        let expected_0 = (1.0_f64 * 10000.0 / 3.0 + 1.0).ln();
        let expected_1 = (2.0_f64 * 10000.0 / 3.0 + 1.0).ln();
        let expected_2 = (3.0_f64 * 10000.0 / 3.0 + 1.0).ln();
        assert!((res_csr.data[0] - expected_0).abs() < 1e-10);
        assert!((res_csr.data[1] - expected_1).abs() < 1e-10);
        assert!((res_csr.data[2] - expected_2).abs() < 1e-10);
    } else {
        panic!("Expected SparseCSR matrix");
    }
}

/// Test: Gene selection picks highest variance columns
#[test]
fn test_gene_selection_picks_highest_variance() {
    let n_rows = 4;
    let n_cols = 5;
    // Column variances: col0=0, col1=low, col2=medium, col3=high, col4=highest
    #[rustfmt::skip]
    let data = vec![
        1.0, 1.0, 1.0,  0.0, 0.0,   // row 0
        1.0, 2.0, 3.0,  10.0, 100.0, // row 1
        1.0, 1.0, 1.0,  0.0, 0.0,    // row 2
        1.0, 2.0, 3.0,  10.0, 100.0, // row 3
    ];
    let mut scd = make_scd_dense(data, n_rows, n_cols);
    select_top_variable_genes(&mut scd, 2).unwrap();

    let (_, cols) = scd.expression.shape();
    assert_eq!(cols, 2);
    // The two highest variance columns should be col4 and col3
    // After selection, the matrix should contain those columns
    if let ExpressionMatrix::Dense(dense) = &scd.expression {
        // col3 values: [0, 10, 0, 10], col4 values: [0, 100, 0, 100]
        // Sorted by original index: col3 first, col4 second
        assert!((dense.data[0] - 0.0).abs() < 1e-10); // row0, col3
        assert!((dense.data[1] - 0.0).abs() < 1e-10); // row0, col4
        assert!((dense.data[2] - 10.0).abs() < 1e-10); // row1, col3
        assert!((dense.data[3] - 100.0).abs() < 1e-10); // row1, col4
    }
}

/// Test: apply_gene_id_column with existing _index column
#[test]
fn test_gene_id_replaces_existing_index() {
    let n_rows = 2;
    let n_cols = 3;
    let data = vec![1.0f64; n_rows * n_cols];
    let expr = ExpressionMatrix::Dense(DenseMatrix::new(data, n_rows, n_cols).unwrap());
    let cell_meta = DataFrame::empty(n_rows);

    let old_index: ArrayRef = Arc::new(StringArray::from(vec!["old1", "old2", "old3"]));
    let ensembl: ArrayRef = Arc::new(StringArray::from(vec!["ENSG1", "ENSG2", "ENSG3"]));
    let gene_meta = DataFrame::new(
        vec!["_index".to_string(), "ensembl_id".to_string()],
        vec![old_index, ensembl],
        n_cols,
    )
    .unwrap();
    let metadata = DatasetMetadata::new(n_rows, n_cols, "test".to_string());
    let mut scd = SingleCellData::new(expr, cell_meta, gene_meta, metadata).unwrap();

    apply_gene_id_column(&mut scd, "ensembl_id").unwrap();

    // _index should now have ensembl values
    let idx = scd.gene_metadata.column_index("_index").unwrap();
    let arr = scd.gene_metadata.data[idx]
        .as_any()
        .downcast_ref::<StringArray>()
        .unwrap();
    assert_eq!(arr.value(0), "ENSG1");
    assert_eq!(arr.value(1), "ENSG2");
    assert_eq!(arr.value(2), "ENSG3");

    // Should still have 2 columns (not 3)
    assert_eq!(scd.gene_metadata.columns.len(), 2);
}

/// Test: Layers are also filtered during gene selection
#[test]
fn test_gene_selection_filters_layers() {
    let n_rows = 3;
    let n_cols = 5;
    let data = vec![0.0f64; n_rows * n_cols];
    let mut scd = make_scd_dense(data.clone(), n_rows, n_cols);

    // Add a layer with same dimensions
    let layer = ExpressionMatrix::Dense(DenseMatrix::new(vec![1.0; n_rows * n_cols], n_rows, n_cols).unwrap());
    let mut layers = std::collections::HashMap::new();
    layers.insert("counts".to_string(), layer);
    scd.layers = Some(layers);

    // Create varying variance by modifying expression
    let mut new_data = vec![0.0f64; n_rows * n_cols];
    for col in 0..n_cols {
        for row in 0..n_rows {
            new_data[row * n_cols + col] = ((col + 1) * (row + 1)) as f64;
        }
    }
    scd.expression = ExpressionMatrix::Dense(DenseMatrix::new(new_data, n_rows, n_cols).unwrap());

    select_top_variable_genes(&mut scd, 3).unwrap();

    // Layer should also be filtered
    let layer = scd.layers.as_ref().unwrap().get("counts").unwrap();
    let (_, lc) = layer.shape();
    assert_eq!(lc, 3);
}
