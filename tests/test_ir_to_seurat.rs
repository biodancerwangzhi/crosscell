//! Integration tests for IR to Seurat conversion

use crosscell::ir::{
    DataFrame, DatasetMetadata, Embedding, ExpressionMatrix, SingleCellData, SparseMatrixCSC,
};
use crosscell::seurat::{ir_to_seurat_rds, seurat_rds_to_ir, write_seurat_rds};
use std::collections::HashMap;

#[test]
fn test_ir_to_seurat_minimal() {
    // Create minimal IR with sparse matrix
    let sparse = SparseMatrixCSC {
        n_rows: 10, // 10 cells
        n_cols: 20, // 20 genes
        indptr: vec![
            0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42,
        ],
        indices: vec![
            0, 5, 1, 6, 2, 7, 3, 8, 4, 9, // First 10 genes
            0, 5, 1, 6, 2, 7, 3, 8, 4, 9, // Next 10 genes
            0, 5, 1, 6, 2, 7, 3, 8, 4, 9, // Next 10 genes
            0, 5, 1, 6, 2, 7, 3, 8, 4, 9, // Last 10 genes
            0, 5, // Last 2 genes
        ],
        data: vec![1.0; 42],
    };

    let expression = ExpressionMatrix::SparseCSC(sparse);
    let cell_metadata = DataFrame::empty(10);
    let gene_metadata = DataFrame::empty(20);
    let metadata = DatasetMetadata::new(10, 20, "test".to_string());

    let ir = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata)
        .expect("Failed to create IR");

    // Convert to Seurat RObject
    let result = ir_to_seurat_rds(&ir);
    assert!(
        result.is_ok(),
        "Failed to convert IR to Seurat: {:?}",
        result.err()
    );
}

#[test]
fn test_write_and_read_seurat_rds() {
    // Create minimal IR
    let sparse = SparseMatrixCSC {
        n_rows: 5,  // 5 cells
        n_cols: 10, // 10 genes
        indptr: vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
        indices: vec![0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0],
        data: vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0],
    };

    let expression = ExpressionMatrix::SparseCSC(sparse);
    let cell_metadata = DataFrame::empty(5);
    let gene_metadata = DataFrame::empty(10);
    let metadata = DatasetMetadata::new(5, 10, "test".to_string());

    let ir = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata)
        .expect("Failed to create IR");

    // Write to RDS file
    let output_path = "tests/data/rust_generated_seurat.rds";
    let result = write_seurat_rds(&ir, output_path);
    assert!(
        result.is_ok(),
        "Failed to write Seurat RDS: {:?}",
        result.err()
    );

    // Verify file exists
    use std::path::Path;
    assert!(Path::new(output_path).exists(), "Output file not created");

    // Try to read it back
    let read_result = seurat_rds_to_ir(output_path);
    assert!(
        read_result.is_ok(),
        "Failed to read back Seurat RDS: {:?}",
        read_result.err()
    );

    let ir_read = read_result.unwrap();
    assert_eq!(ir_read.metadata.n_cells, 5);
    assert_eq!(ir_read.metadata.n_genes, 10);
}

#[test]
fn test_ir_to_seurat_with_embeddings() {
    // Create IR with embeddings
    let sparse = SparseMatrixCSC {
        n_rows: 10,
        n_cols: 20,
        indptr: vec![
            0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42,
        ],
        indices: vec![
            0, 5, 1, 6, 2, 7, 3, 8, 4, 9, 0, 5, 1, 6, 2, 7, 3, 8, 4, 9, 0, 5, 1, 6, 2, 7, 3, 8, 4,
            9, 0, 5, 1, 6, 2, 7, 3, 8, 4, 9, 0, 5,
        ],
        data: vec![1.0; 42],
    };

    let expression = ExpressionMatrix::SparseCSC(sparse);
    let cell_metadata = DataFrame::empty(10);
    let gene_metadata = DataFrame::empty(20);
    let metadata = DatasetMetadata::new(10, 20, "test".to_string());

    // Create embeddings
    let mut embeddings = HashMap::new();
    let pca_data = vec![
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0,
        17.0, 18.0, 19.0, 20.0,
    ];
    let pca =
        Embedding::new("pca".to_string(), pca_data, 10, 2).expect("Failed to create PCA embedding");
    embeddings.insert("pca".to_string(), pca);

    let mut ir = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata)
        .expect("Failed to create IR");
    ir.embeddings = Some(embeddings);

    // Convert to Seurat
    let result = ir_to_seurat_rds(&ir);
    assert!(
        result.is_ok(),
        "Failed to convert IR with embeddings to Seurat: {:?}",
        result.err()
    );

    // Write and read back
    let output_path = "tests/data/rust_generated_seurat_with_embeddings.rds";
    let write_result = write_seurat_rds(&ir, output_path);
    assert!(
        write_result.is_ok(),
        "Failed to write Seurat with embeddings: {:?}",
        write_result.err()
    );

    let read_result = seurat_rds_to_ir(output_path);
    assert!(
        read_result.is_ok(),
        "Failed to read back Seurat with embeddings: {:?}",
        read_result.err()
    );

    let ir_read = read_result.unwrap();
    assert!(ir_read.embeddings.is_some(), "Embeddings not preserved");
    assert!(
        ir_read.embeddings.as_ref().unwrap().contains_key("pca"),
        "PCA embedding not found"
    );
}

#[test]
fn test_ir_to_seurat_with_layers() {
    // Create IR with multiple layers
    let sparse = SparseMatrixCSC {
        n_rows: 5,
        n_cols: 10,
        indptr: vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
        indices: vec![0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0],
        data: vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0],
    };

    let expression = ExpressionMatrix::SparseCSC(sparse.clone());
    let cell_metadata = DataFrame::empty(5);
    let gene_metadata = DataFrame::empty(10);
    let metadata = DatasetMetadata::new(5, 10, "test".to_string());

    // Create additional layer
    let mut layers = HashMap::new();
    let layer_sparse = SparseMatrixCSC {
        n_rows: 5,
        n_cols: 10,
        indptr: vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
        indices: vec![0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0],
        data: vec![2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0],
    };
    layers.insert("SCT".to_string(), ExpressionMatrix::SparseCSC(layer_sparse));

    let mut ir = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata)
        .expect("Failed to create IR");
    ir.layers = Some(layers);

    // Convert to Seurat
    let result = ir_to_seurat_rds(&ir);
    assert!(
        result.is_ok(),
        "Failed to convert IR with layers to Seurat: {:?}",
        result.err()
    );

    // Write and read back
    let output_path = "tests/data/rust_generated_seurat_with_layers.rds";
    let write_result = write_seurat_rds(&ir, output_path);
    assert!(
        write_result.is_ok(),
        "Failed to write Seurat with layers: {:?}",
        write_result.err()
    );

    let read_result = seurat_rds_to_ir(output_path);
    assert!(
        read_result.is_ok(),
        "Failed to read back Seurat with layers: {:?}",
        read_result.err()
    );

    let ir_read = read_result.unwrap();
    assert!(ir_read.layers.is_some(), "Layers not preserved");
    assert!(
        ir_read.layers.as_ref().unwrap().contains_key("SCT"),
        "SCT layer not found"
    );
}
