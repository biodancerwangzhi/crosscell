//! 集成测试：.h5ad 文件部分加载功能
//!
//! 测试 inspect_h5ad 和 read_h5ad_partial 函数。
//! 借鉴 scDIOR 的设计，支持仅加载元数据而不加载表达矩阵。

use crosscell::anndata::{inspect_h5ad, read_h5ad_partial, PartialLoadOptions};
use crosscell::ir::ExpressionMatrix;
use std::path::Path;

/// 检查测试数据文件是否存在
fn test_file_exists(filename: &str) -> bool {
    Path::new("tests/data").join(filename).exists()
}

#[test]
fn test_inspect_h5ad_sparse() {
    if !test_file_exists("small_sparse.h5ad") {
        eprintln!("Skipping test: tests/data/small_sparse.h5ad not found");
        return;
    }

    let result = inspect_h5ad("tests/data/small_sparse.h5ad");
    assert!(
        result.is_ok(),
        "Failed to inspect small_sparse.h5ad: {:?}",
        result.err()
    );

    let info = result.unwrap();

    // 验证基本信息
    assert_eq!(info.n_cells, 100);
    assert_eq!(info.n_genes, 50);
    assert!(info.is_sparse, "Expected sparse matrix");

    // 验证稀疏格式
    assert!(info.sparse_format.is_some());

    // 验证 nnz
    assert!(info.nnz.is_some());
    let nnz = info.nnz.unwrap();
    assert!(
        nnz > 200 && nnz < 300,
        "Expected ~250 non-zero elements, got {}",
        nnz
    );

    // 验证元数据列数
    assert!(
        info.n_obs_columns > 0,
        "Expected some cell metadata columns"
    );
    assert!(
        info.n_var_columns > 0,
        "Expected some gene metadata columns"
    );

    // 验证嵌入名称
    assert!(!info.embedding_names.is_empty(), "Expected some embeddings");

    // 验证层名称
    assert!(!info.layer_names.is_empty(), "Expected some layers");

    // 打印信息
    println!("{}", info.display());

    println!("✓ Successfully inspected small_sparse.h5ad");
}

#[test]
fn test_inspect_h5ad_dense() {
    if !test_file_exists("small_dense.h5ad") {
        eprintln!("Skipping test: tests/data/small_dense.h5ad not found");
        return;
    }

    let result = inspect_h5ad("tests/data/small_dense.h5ad");
    assert!(
        result.is_ok(),
        "Failed to inspect small_dense.h5ad: {:?}",
        result.err()
    );

    let info = result.unwrap();

    // 验证基本信息
    assert_eq!(info.n_cells, 10);
    assert_eq!(info.n_genes, 20);
    assert!(!info.is_sparse, "Expected dense matrix");

    // 验证稀疏格式为 None
    assert!(info.sparse_format.is_none());
    assert!(info.nnz.is_none());

    println!("{}", info.display());

    println!("✓ Successfully inspected small_dense.h5ad");
}

#[test]
fn test_partial_load_metadata_only() {
    if !test_file_exists("small_sparse.h5ad") {
        eprintln!("Skipping test: tests/data/small_sparse.h5ad not found");
        return;
    }

    let options = PartialLoadOptions::metadata_only();
    let result = read_h5ad_partial("tests/data/small_sparse.h5ad", &options);
    assert!(
        result.is_ok(),
        "Failed to read with metadata_only: {:?}",
        result.err()
    );

    let data = result.unwrap();

    // 验证维度
    assert_eq!(data.metadata.n_cells, 100);
    assert_eq!(data.metadata.n_genes, 50);

    // 验证表达矩阵是 Lazy（未加载）
    match &data.expression {
        ExpressionMatrix::Lazy(lazy) => {
            assert_eq!(lazy.shape, (100, 50));
            assert!(lazy.is_sparse);
            println!("  ✓ Expression matrix is lazy (not loaded)");
        }
        _ => panic!("Expected Lazy matrix for metadata_only mode"),
    }

    // 验证元数据已加载
    assert_eq!(data.cell_metadata.n_rows, 100);
    assert!(
        data.cell_metadata.n_cols() > 0,
        "Expected cell metadata to be loaded"
    );

    assert_eq!(data.gene_metadata.n_rows, 50);
    assert!(
        data.gene_metadata.n_cols() > 0,
        "Expected gene metadata to be loaded"
    );

    // 验证嵌入未加载
    assert!(
        data.embeddings.is_none(),
        "Expected embeddings to not be loaded"
    );

    // 验证 layers 未加载
    assert!(data.layers.is_none(), "Expected layers to not be loaded");

    println!("✓ Successfully loaded metadata only from small_sparse.h5ad");
}

#[test]
fn test_partial_load_lazy() {
    if !test_file_exists("small_sparse.h5ad") {
        eprintln!("Skipping test: tests/data/small_sparse.h5ad not found");
        return;
    }

    let options = PartialLoadOptions::lazy();
    let result = read_h5ad_partial("tests/data/small_sparse.h5ad", &options);
    assert!(
        result.is_ok(),
        "Failed to read with lazy mode: {:?}",
        result.err()
    );

    let data = result.unwrap();

    // 验证表达矩阵是 Lazy
    match &data.expression {
        ExpressionMatrix::Lazy(lazy) => {
            assert_eq!(lazy.shape, (100, 50));
            assert!(lazy.is_sparse);
            assert!(lazy.nnz.is_some());
            println!(
                "  ✓ Expression matrix is lazy with nnz={}",
                lazy.nnz.unwrap()
            );
        }
        _ => panic!("Expected Lazy matrix for lazy mode"),
    }

    // 验证元数据已加载
    assert!(data.cell_metadata.n_cols() > 0);
    assert!(data.gene_metadata.n_cols() > 0);

    // 验证嵌入已加载
    assert!(
        data.embeddings.is_some(),
        "Expected embeddings to be loaded"
    );

    // 验证 layers 已加载
    assert!(data.layers.is_some(), "Expected layers to be loaded");

    println!("✓ Successfully loaded with lazy mode from small_sparse.h5ad");
}

#[test]
fn test_partial_load_full() {
    if !test_file_exists("small_sparse.h5ad") {
        eprintln!("Skipping test: tests/data/small_sparse.h5ad not found");
        return;
    }

    let options = PartialLoadOptions::full();
    let result = read_h5ad_partial("tests/data/small_sparse.h5ad", &options);
    assert!(
        result.is_ok(),
        "Failed to read with full mode: {:?}",
        result.err()
    );

    let data = result.unwrap();

    // 验证表达矩阵已完整加载
    match &data.expression {
        ExpressionMatrix::SparseCSR(csr) => {
            assert_eq!(csr.n_rows, 100);
            assert_eq!(csr.n_cols, 50);
            println!(
                "  ✓ Expression matrix fully loaded: {} × {}",
                csr.n_rows, csr.n_cols
            );
        }
        _ => panic!("Expected SparseCSR matrix for full mode"),
    }

    // 验证所有数据已加载
    assert!(data.cell_metadata.n_cols() > 0);
    assert!(data.gene_metadata.n_cols() > 0);
    assert!(data.embeddings.is_some());
    assert!(data.layers.is_some());

    println!("✓ Successfully loaded with full mode from small_sparse.h5ad");
}

#[test]
fn test_partial_load_custom_options() {
    if !test_file_exists("small_sparse.h5ad") {
        eprintln!("Skipping test: tests/data/small_sparse.h5ad not found");
        return;
    }

    // 自定义选项：只加载表达矩阵和元数据，不加载嵌入和 layers
    let options = PartialLoadOptions {
        load_expression: true,
        load_cell_metadata: true,
        load_gene_metadata: true,
        load_embeddings: false,
        load_layers: false,
        load_pairwise: false,
        load_spatial: false,
        lazy_load: false,
    };

    let result = read_h5ad_partial("tests/data/small_sparse.h5ad", &options);
    assert!(
        result.is_ok(),
        "Failed to read with custom options: {:?}",
        result.err()
    );

    let data = result.unwrap();

    // 验证表达矩阵已加载
    match &data.expression {
        ExpressionMatrix::SparseCSR(_) => {
            println!("  ✓ Expression matrix loaded");
        }
        _ => panic!("Expected SparseCSR matrix"),
    }

    // 验证元数据已加载
    assert!(data.cell_metadata.n_cols() > 0);
    assert!(data.gene_metadata.n_cols() > 0);

    // 验证嵌入未加载
    assert!(
        data.embeddings.is_none(),
        "Expected embeddings to not be loaded"
    );

    // 验证 layers 未加载
    assert!(data.layers.is_none(), "Expected layers to not be loaded");

    // 验证 pairwise 未加载
    assert!(
        data.cell_pairwise.is_none(),
        "Expected cell_pairwise to not be loaded"
    );
    assert!(
        data.gene_pairwise.is_none(),
        "Expected gene_pairwise to not be loaded"
    );

    println!("✓ Successfully loaded with custom options from small_sparse.h5ad");
}

#[test]
fn test_h5ad_info_estimate_memory() {
    if !test_file_exists("small_sparse.h5ad") {
        eprintln!("Skipping test: tests/data/small_sparse.h5ad not found");
        return;
    }

    let info = inspect_h5ad("tests/data/small_sparse.h5ad").unwrap();

    let estimated_memory = info.estimate_full_memory();

    // 验证估算值合理（应该大于 0）
    assert!(estimated_memory > 0, "Expected positive memory estimate");

    println!(
        "  Estimated full memory: {:.2} MB",
        estimated_memory as f64 / 1024.0 / 1024.0
    );
    println!("✓ Memory estimation works correctly");
}

#[test]
fn test_inspect_real_dataset() {
    // 测试真实数据集（如果存在）
    let pbmc_path = "data/pbmc3k_raw.h5ad";
    if !Path::new(pbmc_path).exists() {
        eprintln!("Skipping test: {} not found", pbmc_path);
        return;
    }

    let result = inspect_h5ad(pbmc_path);
    assert!(
        result.is_ok(),
        "Failed to inspect pbmc3k: {:?}",
        result.err()
    );

    let info = result.unwrap();

    println!("📊 PBMC3K Dataset Info:");
    println!("{}", info.display());

    // 验证基本信息
    assert!(info.n_cells > 0);
    assert!(info.n_genes > 0);

    println!(
        "✓ Successfully inspected real dataset: {} cells × {} genes",
        info.n_cells, info.n_genes
    );
}
