//! 集成测试：.h5ad 文件读取
//!
//! 测试从 HDF5 格式的 AnnData 文件读取表达矩阵。
//!
//! 运行测试前，需要先运行 Python 脚本生成测试数据：
//! ```bash
//! python3 tests/create_test_h5ad.py
//! ```

use crosscell::anndata::read_h5ad;
use crosscell::ir::ExpressionMatrix;
use std::path::Path;

/// 检查测试数据文件是否存在
fn test_file_exists(filename: &str) -> bool {
    Path::new("tests/data").join(filename).exists()
}

#[test]
fn test_read_small_sparse_h5ad() {
    if !test_file_exists("small_sparse.h5ad") {
        eprintln!("Skipping test: tests/data/small_sparse.h5ad not found");
        eprintln!("Run: python3 tests/create_test_h5ad.py");
        return;
    }

    let result = read_h5ad("tests/data/small_sparse.h5ad");
    assert!(
        result.is_ok(),
        "Failed to read small_sparse.h5ad: {:?}",
        result.err()
    );

    let data = result.unwrap();

    // 验证维度
    assert_eq!(data.metadata.n_cells, 100);
    assert_eq!(data.metadata.n_genes, 50);

    // 验证表达矩阵类型
    match data.expression {
        ExpressionMatrix::SparseCSR(ref csr) => {
            assert_eq!(csr.n_rows, 100);
            assert_eq!(csr.n_cols, 50);
            // 验证稀疏性（应该约 95% 稀疏，即约 250 个非零元素）
            let nnz = csr.data.len();
            assert!(
                nnz > 200 && nnz < 300,
                "Expected ~250 non-zero elements, got {}",
                nnz
            );
        }
        _ => panic!("Expected SparseCSR matrix"),
    }

    // 验证细胞元数据
    assert_eq!(data.cell_metadata.n_rows, 100);
    let n_cell_cols = data.cell_metadata.n_cols();
    println!("  Cell metadata columns: {}", n_cell_cols);

    // 检查特定列是否存在（如果有的话）
    if data.cell_metadata.column("n_genes").is_some() {
        println!("    ✓ n_genes column found");
    }
    if data.cell_metadata.column("cell_type").is_some() {
        println!("    ✓ cell_type column found");
    }

    // 验证基因元数据
    assert_eq!(data.gene_metadata.n_rows, 50);
    let n_gene_cols = data.gene_metadata.n_cols();
    println!("  Gene metadata columns: {}", n_gene_cols);

    if data.gene_metadata.column("highly_variable").is_some() {
        println!("    ✓ highly_variable column found");
    }

    // 验证降维嵌入
    if let Some(ref embeddings) = data.embeddings {
        println!("  Embeddings: {} found", embeddings.len());

        // 验证 PCA
        if let Some(pca) = embeddings.get("X_pca") {
            assert_eq!(pca.n_rows, 100, "PCA should have 100 cells");
            assert_eq!(pca.n_cols, 20, "PCA should have 20 components");
            println!("    ✓ X_pca: {} × {}", pca.n_rows, pca.n_cols);
        }

        // 验证 UMAP
        if let Some(umap) = embeddings.get("X_umap") {
            assert_eq!(umap.n_rows, 100, "UMAP should have 100 cells");
            assert_eq!(umap.n_cols, 2, "UMAP should have 2 components");
            println!("    ✓ X_umap: {} × {}", umap.n_rows, umap.n_cols);
        }
    } else {
        println!("  ⚠ No embeddings found");
    }

    // 验证 layers
    if let Some(ref layers) = data.layers {
        println!("  Layers: {} found", layers.len());

        // 验证 counts layer
        if let Some(counts) = layers.get("counts") {
            let (n_rows, n_cols) = counts.shape();
            assert_eq!(n_rows, 100, "counts layer should have 100 cells");
            assert_eq!(n_cols, 50, "counts layer should have 50 genes");
            println!("    ✓ counts: {} × {}", n_rows, n_cols);
        }

        // 验证 log1p layer
        if let Some(log1p) = layers.get("log1p") {
            let (n_rows, n_cols) = log1p.shape();
            assert_eq!(n_rows, 100, "log1p layer should have 100 cells");
            assert_eq!(n_cols, 50, "log1p layer should have 50 genes");
            println!("    ✓ log1p: {} × {}", n_rows, n_cols);
        }
    } else {
        println!("  ⚠ No layers found");
    }

    // 验证细胞-细胞 pairwise 矩阵
    if let Some(ref cell_pairwise) = data.cell_pairwise {
        println!("  Cell pairwise matrices: {} found", cell_pairwise.len());

        // 验证 connectivities
        if let Some(conn) = cell_pairwise.get("connectivities") {
            let (n_rows, n_cols) = conn.matrix.shape();
            assert_eq!(n_rows, 100, "connectivities should be 100×100");
            assert_eq!(n_cols, 100, "connectivities should be 100×100");
            println!("    ✓ connectivities: {} × {}", n_rows, n_cols);
        }

        // 验证 distances
        if let Some(dist) = cell_pairwise.get("distances") {
            let (n_rows, n_cols) = dist.matrix.shape();
            assert_eq!(n_rows, 100, "distances should be 100×100");
            assert_eq!(n_cols, 100, "distances should be 100×100");
            println!("    ✓ distances: {} × {}", n_rows, n_cols);
        }
    } else {
        println!("  ⚠ No cell pairwise matrices found");
    }

    // 验证基因-基因 pairwise 矩阵
    if let Some(ref gene_pairwise) = data.gene_pairwise {
        println!("  Gene pairwise matrices: {} found", gene_pairwise.len());

        // 验证 gene_correlation
        if let Some(corr) = gene_pairwise.get("gene_correlation") {
            let (n_rows, n_cols) = corr.matrix.shape();
            assert_eq!(n_rows, 50, "gene_correlation should be 50×50");
            assert_eq!(n_cols, 50, "gene_correlation should be 50×50");
            println!("    ✓ gene_correlation: {} × {}", n_rows, n_cols);
        }
    } else {
        println!("  ⚠ No gene pairwise matrices found");
    }

    println!(
        "✓ Successfully read small_sparse.h5ad: {} cells × {} genes",
        data.metadata.n_cells, data.metadata.n_genes
    );
    println!("  Cell metadata: {} columns", data.cell_metadata.n_cols());
    println!("  Gene metadata: {} columns", data.gene_metadata.n_cols());
}

#[test]
fn test_read_small_dense_h5ad() {
    if !test_file_exists("small_dense.h5ad") {
        eprintln!("Skipping test: tests/data/small_dense.h5ad not found");
        eprintln!("Run: python3 tests/create_test_h5ad.py");
        return;
    }

    let result = read_h5ad("tests/data/small_dense.h5ad");
    assert!(
        result.is_ok(),
        "Failed to read small_dense.h5ad: {:?}",
        result.err()
    );

    let data = result.unwrap();

    // 验证维度
    assert_eq!(data.metadata.n_cells, 10);
    assert_eq!(data.metadata.n_genes, 20);

    // 验证表达矩阵类型
    match data.expression {
        ExpressionMatrix::Dense(ref dense) => {
            assert_eq!(dense.n_rows, 10);
            assert_eq!(dense.n_cols, 20);
            assert_eq!(dense.data.len(), 200);
        }
        _ => panic!("Expected Dense matrix"),
    }

    // 验证元数据
    assert_eq!(data.cell_metadata.n_rows, 10);
    let n_cell_cols = data.cell_metadata.n_cols();
    println!("  Cell metadata columns: {}", n_cell_cols);

    assert_eq!(data.gene_metadata.n_rows, 20);
    let n_gene_cols = data.gene_metadata.n_cols();
    println!("  Gene metadata columns: {}", n_gene_cols);

    // 验证降维嵌入
    if let Some(ref embeddings) = data.embeddings {
        println!("  Embeddings: {} found", embeddings.len());

        // 验证 PCA
        if let Some(pca) = embeddings.get("X_pca") {
            assert_eq!(pca.n_rows, 10, "PCA should have 10 cells");
            assert_eq!(pca.n_cols, 5, "PCA should have 5 components");
            println!("    ✓ X_pca: {} × {}", pca.n_rows, pca.n_cols);
        }
    } else {
        println!("  ⚠ No embeddings found");
    }

    // 验证 layers
    if let Some(ref layers) = data.layers {
        println!("  Layers: {} found", layers.len());

        // 验证 raw layer
        if let Some(raw) = layers.get("raw") {
            let (n_rows, n_cols) = raw.shape();
            assert_eq!(n_rows, 10, "raw layer should have 10 cells");
            assert_eq!(n_cols, 20, "raw layer should have 20 genes");
            println!("    ✓ raw: {} × {}", n_rows, n_cols);
        }
    } else {
        println!("  ⚠ No layers found");
    }

    println!(
        "✓ Successfully read small_dense.h5ad: {} cells × {} genes",
        data.metadata.n_cells, data.metadata.n_genes
    );
    println!("  Cell metadata: {} columns", data.cell_metadata.n_cols());
    println!("  Gene metadata: {} columns", data.gene_metadata.n_cols());
}

#[test]
fn test_read_empty_h5ad() {
    if !test_file_exists("empty.h5ad") {
        eprintln!("Skipping test: tests/data/empty.h5ad not found");
        eprintln!("Run: python3 tests/create_test_h5ad.py");
        return;
    }

    let result = read_h5ad("tests/data/empty.h5ad");
    assert!(
        result.is_ok(),
        "Failed to read empty.h5ad: {:?}",
        result.err()
    );

    let data = result.unwrap();

    // 验证维度
    assert_eq!(data.metadata.n_cells, 0);
    assert_eq!(data.metadata.n_genes, 0);

    println!(
        "✓ Successfully read empty.h5ad: {} cells × {} genes",
        data.metadata.n_cells, data.metadata.n_genes
    );
}

#[test]
fn test_read_single_row_h5ad() {
    if !test_file_exists("single_row.h5ad") {
        eprintln!("Skipping test: tests/data/single_row.h5ad not found");
        eprintln!("Run: python3 tests/create_test_h5ad.py");
        return;
    }

    let result = read_h5ad("tests/data/single_row.h5ad");
    assert!(
        result.is_ok(),
        "Failed to read single_row.h5ad: {:?}",
        result.err()
    );

    let data = result.unwrap();

    // 验证维度
    assert_eq!(data.metadata.n_cells, 1);
    assert_eq!(data.metadata.n_genes, 1000);

    println!(
        "✓ Successfully read single_row.h5ad: {} cells × {} genes",
        data.metadata.n_cells, data.metadata.n_genes
    );
}

#[test]
fn test_read_single_column_h5ad() {
    if !test_file_exists("single_column.h5ad") {
        eprintln!("Skipping test: tests/data/single_column.h5ad not found");
        eprintln!("Run: python3 tests/create_test_h5ad.py");
        return;
    }

    let result = read_h5ad("tests/data/single_column.h5ad");
    assert!(
        result.is_ok(),
        "Failed to read single_column.h5ad: {:?}",
        result.err()
    );

    let data = result.unwrap();

    // 验证维度
    assert_eq!(data.metadata.n_cells, 1000);
    assert_eq!(data.metadata.n_genes, 1);

    println!(
        "✓ Successfully read single_column.h5ad: {} cells × {} genes",
        data.metadata.n_cells, data.metadata.n_genes
    );
}

#[test]
fn test_read_medium_sparse_h5ad() {
    if !test_file_exists("medium_sparse.h5ad") {
        eprintln!("Skipping test: tests/data/medium_sparse.h5ad not found");
        eprintln!("Run: python3 tests/create_test_h5ad.py");
        return;
    }

    let result = read_h5ad("tests/data/medium_sparse.h5ad");
    assert!(
        result.is_ok(),
        "Failed to read medium_sparse.h5ad: {:?}",
        result.err()
    );

    let data = result.unwrap();

    // 验证维度
    assert_eq!(data.metadata.n_cells, 1000);
    assert_eq!(data.metadata.n_genes, 500);

    // 验证稀疏性
    match data.expression {
        ExpressionMatrix::SparseCSR(ref csr) => {
            let nnz = csr.data.len();
            let total = 1000 * 500;
            let sparsity = 1.0 - (nnz as f64 / total as f64);
            assert!(
                sparsity > 0.94 && sparsity < 0.96,
                "Expected ~95% sparsity, got {:.2}%",
                sparsity * 100.0
            );
        }
        _ => panic!("Expected SparseCSR matrix"),
    }

    println!(
        "✓ Successfully read medium_sparse.h5ad: {} cells × {} genes",
        data.metadata.n_cells, data.metadata.n_genes
    );
}

#[test]
fn test_read_nonexistent_file() {
    let result = read_h5ad("tests/data/nonexistent.h5ad");
    assert!(result.is_err(), "Should fail to read nonexistent file");

    println!("✓ Correctly handled nonexistent file");
}
