//! 集成测试：.h5ad 文件写入
//!
//! 测试将 IR 数据写入 HDF5 格式的 AnnData 文件。

use crosscell::anndata::{read_h5ad, write_h5ad};
use crosscell::ir::{
    DatasetMetadata, ExpressionMatrix, SparseMatrixCSR, DenseMatrix, DataFrame, SingleCellData,
};
use std::path::Path;
use tempfile::NamedTempFile;

#[test]
fn test_write_sparse_csr_matrix() {
    // 创建一个简单的稀疏 CSR 矩阵
    // 3 细胞 × 4 基因
    // [[1.0, 0.0, 2.0, 0.0],
    //  [0.0, 3.0, 0.0, 4.0],
    //  [5.0, 0.0, 0.0, 6.0]]
    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    let indices = vec![0, 2, 1, 3, 0, 3];
    let indptr = vec![0, 2, 4, 6];
    let csr = SparseMatrixCSR::new(data, indices, indptr, 3, 4).unwrap();
    
    // 创建 IR 数据
    let expression = ExpressionMatrix::SparseCSR(csr);
    let cell_metadata = DataFrame::empty(3);
    let gene_metadata = DataFrame::empty(4);
    let metadata = DatasetMetadata::new(3, 4, "test".to_string());
    
    let data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata).unwrap();
    
    // 写入临时文件
    let temp_file = NamedTempFile::new().unwrap();
    let temp_path = temp_file.path().with_extension("h5ad");
    
    let result = write_h5ad(&data, &temp_path);
    assert!(result.is_ok(), "Failed to write h5ad: {:?}", result.err());
    
    // 读取并验证
    let read_result = read_h5ad(&temp_path);
    assert!(read_result.is_ok(), "Failed to read written h5ad: {:?}", read_result.err());
    
    let read_data = read_result.unwrap();
    assert_eq!(read_data.metadata.n_cells, 3);
    assert_eq!(read_data.metadata.n_genes, 4);
    
    // 验证矩阵内容
    match read_data.expression {
        ExpressionMatrix::SparseCSR(ref read_csr) => {
            assert_eq!(read_csr.n_rows, 3);
            assert_eq!(read_csr.n_cols, 4);
            assert_eq!(read_csr.data.len(), 6);
            
            // 验证数值（允许小误差）
            for (i, &val) in read_csr.data.iter().enumerate() {
                let expected = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0][i];
                assert!((val - expected).abs() < 1e-10, 
                    "Value mismatch at index {}: expected {}, got {}", i, expected, val);
            }
            
            // 验证索引
            assert_eq!(read_csr.indices, vec![0, 2, 1, 3, 0, 3]);
            assert_eq!(read_csr.indptr, vec![0, 2, 4, 6]);
        }
        _ => panic!("Expected SparseCSR matrix"),
    }
    
    // 清理临时文件
    std::fs::remove_file(&temp_path).ok();
    
    println!("✓ Successfully wrote and verified sparse CSR matrix");
}

#[test]
fn test_write_dense_matrix() {
    // 创建一个简单的稠密矩阵
    // 2 细胞 × 3 基因
    // [[1.0, 2.0, 3.0],
    //  [4.0, 5.0, 6.0]]
    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    let dense = DenseMatrix::new(data, 2, 3).unwrap();
    
    // 创建 IR 数据
    let expression = ExpressionMatrix::Dense(dense);
    let cell_metadata = DataFrame::empty(2);
    let gene_metadata = DataFrame::empty(3);
    let metadata = DatasetMetadata::new(2, 3, "test".to_string());
    
    let data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata).unwrap();
    
    // 写入临时文件
    let temp_file = NamedTempFile::new().unwrap();
    let temp_path = temp_file.path().with_extension("h5ad");
    
    let result = write_h5ad(&data, &temp_path);
    assert!(result.is_ok(), "Failed to write h5ad: {:?}", result.err());
    
    // 读取并验证
    let read_result = read_h5ad(&temp_path);
    assert!(read_result.is_ok(), "Failed to read written h5ad: {:?}", read_result.err());
    
    let read_data = read_result.unwrap();
    assert_eq!(read_data.metadata.n_cells, 2);
    assert_eq!(read_data.metadata.n_genes, 3);
    
    // 验证矩阵内容
    match read_data.expression {
        ExpressionMatrix::Dense(ref read_dense) => {
            assert_eq!(read_dense.n_rows, 2);
            assert_eq!(read_dense.n_cols, 3);
            assert_eq!(read_dense.data.len(), 6);
            
            // 验证数值（允许小误差）
            for (i, &val) in read_dense.data.iter().enumerate() {
                let expected = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0][i];
                assert!((val - expected).abs() < 1e-10, 
                    "Value mismatch at index {}: expected {}, got {}", i, expected, val);
            }
        }
        _ => panic!("Expected Dense matrix"),
    }
    
    // 清理临时文件
    std::fs::remove_file(&temp_path).ok();
    
    println!("✓ Successfully wrote and verified dense matrix");
}

#[test]
fn test_write_empty_matrix() {
    // 创建空矩阵
    let data = vec![];
    let indices = vec![];
    let indptr = vec![0];
    let csr = SparseMatrixCSR::new(data, indices, indptr, 0, 0).unwrap();
    
    // 创建 IR 数据
    let expression = ExpressionMatrix::SparseCSR(csr);
    let cell_metadata = DataFrame::empty(0);
    let gene_metadata = DataFrame::empty(0);
    let metadata = DatasetMetadata::new(0, 0, "test".to_string());
    
    let data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata).unwrap();
    
    // 写入临时文件
    let temp_file = NamedTempFile::new().unwrap();
    let temp_path = temp_file.path().with_extension("h5ad");
    
    let result = write_h5ad(&data, &temp_path);
    assert!(result.is_ok(), "Failed to write empty h5ad: {:?}", result.err());
    
    // 读取并验证
    let read_result = read_h5ad(&temp_path);
    assert!(read_result.is_ok(), "Failed to read written empty h5ad: {:?}", read_result.err());
    
    let read_data = read_result.unwrap();
    assert_eq!(read_data.metadata.n_cells, 0);
    assert_eq!(read_data.metadata.n_genes, 0);
    
    // 清理临时文件
    std::fs::remove_file(&temp_path).ok();
    
    println!("✓ Successfully wrote and verified empty matrix");
}

#[test]
fn test_roundtrip_sparse_matrix() {
    // 读取测试文件
    if !Path::new("tests/data/small_sparse.h5ad").exists() {
        eprintln!("Skipping test: tests/data/small_sparse.h5ad not found");
        return;
    }
    
    let original_data = read_h5ad("tests/data/small_sparse.h5ad").unwrap();
    
    // 写入临时文件
    let temp_file = NamedTempFile::new().unwrap();
    let temp_path = temp_file.path().with_extension("h5ad");
    
    write_h5ad(&original_data, &temp_path).unwrap();
    
    // 读取并比较
    let roundtrip_data = read_h5ad(&temp_path).unwrap();
    
    // 验证维度
    assert_eq!(roundtrip_data.metadata.n_cells, original_data.metadata.n_cells);
    assert_eq!(roundtrip_data.metadata.n_genes, original_data.metadata.n_genes);
    
    // 验证矩阵类型和内容
    match (&original_data.expression, &roundtrip_data.expression) {
        (ExpressionMatrix::SparseCSR(orig_csr), ExpressionMatrix::SparseCSR(rt_csr)) => {
            assert_eq!(rt_csr.n_rows, orig_csr.n_rows);
            assert_eq!(rt_csr.n_cols, orig_csr.n_cols);
            assert_eq!(rt_csr.data.len(), orig_csr.data.len());
            
            // 验证数值相关性
            let mut sum_orig = 0.0;
            let mut sum_rt = 0.0;
            let mut sum_sq_orig = 0.0;
            let mut sum_sq_rt = 0.0;
            let mut sum_prod = 0.0;
            
            for i in 0..orig_csr.data.len() {
                let o = orig_csr.data[i];
                let r = rt_csr.data[i];
                sum_orig += o;
                sum_rt += r;
                sum_sq_orig += o * o;
                sum_sq_rt += r * r;
                sum_prod += o * r;
            }
            
            let n = orig_csr.data.len() as f64;
            let correlation = (n * sum_prod - sum_orig * sum_rt) / 
                ((n * sum_sq_orig - sum_orig * sum_orig).sqrt() * 
                 (n * sum_sq_rt - sum_rt * sum_rt).sqrt());
            
            assert!(correlation > 0.9999, "Correlation too low: {}", correlation);
            println!("  Correlation: {:.6}", correlation);
        }
        _ => panic!("Matrix type mismatch"),
    }
    
    // 清理临时文件
    std::fs::remove_file(&temp_path).ok();
    
    println!("✓ Successfully completed roundtrip test for sparse matrix");
}

#[test]
fn test_write_metadata() {
    use arrow::array::{Float64Array, StringArray, DictionaryArray, Int32Array};
    use arrow::datatypes::Int32Type;
    use std::sync::Arc;
    
    // 创建包含元数据的 IR 数据
    // 3 细胞 × 2 基因
    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    let csr = SparseMatrixCSR::new(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0], vec![0, 1, 0, 1, 0, 1], vec![0, 2, 4, 6], 3, 2).unwrap();
    let expression = ExpressionMatrix::SparseCSR(csr);
    
    // 创建细胞元数据（包含不同类型的列）
    let cell_col1 = Arc::new(Float64Array::from(vec![1.5, 2.5, 3.5])) as arrow::array::ArrayRef;
    let cell_col2 = Arc::new(StringArray::from(vec!["A", "B", "C"])) as arrow::array::ArrayRef;
    
    // 创建 Categorical 列
    let keys = Int32Array::from(vec![0, 1, 0]);
    let values = StringArray::from(vec!["type1", "type2"]);
    let dict_array = DictionaryArray::<Int32Type>::try_new(keys, Arc::new(values)).unwrap();
    let cell_col3 = Arc::new(dict_array) as arrow::array::ArrayRef;
    
    let cell_metadata = DataFrame::new(
        vec!["n_genes".to_string(), "cell_id".to_string(), "cell_type".to_string()],
        vec![cell_col1, cell_col2, cell_col3],
        3,
    ).unwrap();
    
    // 创建基因元数据
    let gene_col1 = Arc::new(StringArray::from(vec!["gene1", "gene2"])) as arrow::array::ArrayRef;
    let gene_metadata = DataFrame::new(
        vec!["gene_name".to_string()],
        vec![gene_col1],
        2,
    ).unwrap();
    
    let metadata = DatasetMetadata::new(3, 2, "test".to_string());
    let data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata).unwrap();
    
    // 写入临时文件
    let temp_file = NamedTempFile::new().unwrap();
    let temp_path = temp_file.path().with_extension("h5ad");
    
    let result = write_h5ad(&data, &temp_path);
    assert!(result.is_ok(), "Failed to write h5ad with metadata: {:?}", result.err());
    
    // 读取并验证
    let read_result = read_h5ad(&temp_path);
    assert!(read_result.is_ok(), "Failed to read written h5ad: {:?}", read_result.err());
    
    let read_data = read_result.unwrap();
    
    // 验证细胞元数据
    assert_eq!(read_data.cell_metadata.n_cols(), 3);
    assert_eq!(read_data.cell_metadata.columns, vec!["cell_id", "cell_type", "n_genes"]);
    
    // 验证基因元数据
    assert_eq!(read_data.gene_metadata.n_cols(), 1);
    assert_eq!(read_data.gene_metadata.columns, vec!["gene_name"]);
    
    // 清理临时文件
    std::fs::remove_file(&temp_path).ok();
    
    println!("✓ Successfully wrote and verified metadata");
}

#[test]
fn test_write_metadata_roundtrip_with_python() {
    use arrow::array::{Float64Array, StringArray, DictionaryArray, Int32Array};
    use arrow::datatypes::Int32Type;
    use std::sync::Arc;
    
    // 创建包含元数据的 IR 数据
    let csr = SparseMatrixCSR::new(
        vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        vec![0, 1, 0, 1, 0, 1],
        vec![0, 2, 4, 6],
        3,
        2,
    ).unwrap();
    let expression = ExpressionMatrix::SparseCSR(csr);
    
    // 创建细胞元数据
    let cell_col1 = Arc::new(Float64Array::from(vec![1.5, 2.5, 3.5])) as arrow::array::ArrayRef;
    let cell_col2 = Arc::new(StringArray::from(vec!["A", "B", "C"])) as arrow::array::ArrayRef;
    
    // 创建 Categorical 列
    let keys = Int32Array::from(vec![0, 1, 0]);
    let values = StringArray::from(vec!["type1", "type2"]);
    let dict_array = DictionaryArray::<Int32Type>::try_new(keys, Arc::new(values)).unwrap();
    let cell_col3 = Arc::new(dict_array) as arrow::array::ArrayRef;
    
    let cell_metadata = DataFrame::new(
        vec!["n_genes".to_string(), "cell_id".to_string(), "cell_type".to_string()],
        vec![cell_col1, cell_col2, cell_col3],
        3,
    ).unwrap();
    
    // 创建基因元数据
    let gene_col1 = Arc::new(StringArray::from(vec!["gene1", "gene2"])) as arrow::array::ArrayRef;
    let gene_metadata = DataFrame::new(
        vec!["gene_name".to_string()],
        vec![gene_col1],
        2,
    ).unwrap();
    
    let metadata = DatasetMetadata::new(3, 2, "test".to_string());
    let data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata).unwrap();
    
    // 写入到固定路径，方便 Python 验证
    let test_path = "tests/data/rust_written_metadata.h5ad";
    
    let result = write_h5ad(&data, test_path);
    assert!(result.is_ok(), "Failed to write h5ad: {:?}", result.err());
    
    println!("✓ Wrote test file to: {}", test_path);
    println!("  Run Python verification:");
    println!("  python3 -c \"import anndata; ad = anndata.read_h5ad('{}'); print(ad)\"", test_path);
    
    // 读取并验证往返
    let read_result = read_h5ad(test_path);
    assert!(read_result.is_ok(), "Failed to read written h5ad: {:?}", read_result.err());
    
    let read_data = read_result.unwrap();
    assert_eq!(read_data.cell_metadata.n_cols(), 3);
    assert_eq!(read_data.gene_metadata.n_cols(), 1);
    
    println!("✓ Roundtrip test passed");
}


#[test]
fn test_write_embeddings() {
    use crosscell::ir::Embedding;
    use std::collections::HashMap;
    
    // 创建基础数据
    let csr = SparseMatrixCSR::new(
        vec![1.0, 2.0, 3.0, 4.0],
        vec![0, 1, 0, 1],
        vec![0, 2, 4],
        2,
        2,
    ).unwrap();
    let expression = ExpressionMatrix::SparseCSR(csr);
    let cell_metadata = DataFrame::empty(2);
    let gene_metadata = DataFrame::empty(2);
    let metadata = DatasetMetadata::new(2, 2, "test".to_string());
    
    let mut data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata).unwrap();
    
    // 添加嵌入
    let mut embeddings = HashMap::new();
    
    // PCA: 2 细胞 × 3 成分
    let pca_data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    let pca = Embedding::new("X_pca".to_string(), pca_data, 2, 3).unwrap();
    embeddings.insert("X_pca".to_string(), pca);
    
    // UMAP: 2 细胞 × 2 成分
    let umap_data = vec![0.1, 0.2, 0.3, 0.4];
    let umap = Embedding::new("X_umap".to_string(), umap_data, 2, 2).unwrap();
    embeddings.insert("X_umap".to_string(), umap);
    
    data.embeddings = Some(embeddings);
    
    // 写入临时文件
    let temp_file = NamedTempFile::new().unwrap();
    let temp_path = temp_file.path().with_extension("h5ad");
    
    let result = write_h5ad(&data, &temp_path);
    assert!(result.is_ok(), "Failed to write h5ad with embeddings: {:?}", result.err());
    
    // 读取并验证
    let read_result = read_h5ad(&temp_path);
    assert!(read_result.is_ok(), "Failed to read written h5ad: {:?}", read_result.err());
    
    let read_data = read_result.unwrap();
    
    // 验证嵌入存在
    assert!(read_data.embeddings.is_some());
    let read_embeddings = read_data.embeddings.unwrap();
    assert_eq!(read_embeddings.len(), 2);
    
    // 验证 PCA
    let read_pca = read_embeddings.get("X_pca").unwrap();
    assert_eq!(read_pca.n_rows, 2);
    assert_eq!(read_pca.n_cols, 3);
    assert_eq!(read_pca.data.len(), 6);
    
    // 验证 UMAP
    let read_umap = read_embeddings.get("X_umap").unwrap();
    assert_eq!(read_umap.n_rows, 2);
    assert_eq!(read_umap.n_cols, 2);
    assert_eq!(read_umap.data.len(), 4);
    
    // 清理临时文件
    std::fs::remove_file(&temp_path).ok();
    
    println!("✓ Successfully wrote and verified embeddings");
}

#[test]
fn test_write_layers() {
    use std::collections::HashMap;
    
    // 创建基础数据
    let csr = SparseMatrixCSR::new(
        vec![1.0, 2.0, 3.0, 4.0],
        vec![0, 1, 0, 1],
        vec![0, 2, 4],
        2,
        3,
    ).unwrap();
    let expression = ExpressionMatrix::SparseCSR(csr);
    let cell_metadata = DataFrame::empty(2);
    let gene_metadata = DataFrame::empty(3);
    let metadata = DatasetMetadata::new(2, 3, "test".to_string());
    
    let mut data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata).unwrap();
    
    // 添加 layers
    let mut layers = HashMap::new();
    
    // counts layer (稀疏)
    let counts_csr = SparseMatrixCSR::new(
        vec![10.0, 20.0, 30.0, 40.0],
        vec![0, 1, 0, 1],
        vec![0, 2, 4],
        2,
        3,
    ).unwrap();
    layers.insert("counts".to_string(), ExpressionMatrix::SparseCSR(counts_csr));
    
    // log1p layer (稠密)
    let log1p_dense = DenseMatrix::new(
        vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
        2,
        3,
    ).unwrap();
    layers.insert("log1p".to_string(), ExpressionMatrix::Dense(log1p_dense));
    
    data.layers = Some(layers);
    
    // 写入临时文件
    let temp_file = NamedTempFile::new().unwrap();
    let temp_path = temp_file.path().with_extension("h5ad");
    
    let result = write_h5ad(&data, &temp_path);
    assert!(result.is_ok(), "Failed to write h5ad with layers: {:?}", result.err());
    
    // 读取并验证
    let read_result = read_h5ad(&temp_path);
    assert!(read_result.is_ok(), "Failed to read written h5ad: {:?}", read_result.err());
    
    let read_data = read_result.unwrap();
    
    // 验证 layers 存在
    assert!(read_data.layers.is_some());
    let read_layers = read_data.layers.unwrap();
    assert_eq!(read_layers.len(), 2);
    
    // 验证 counts layer
    let read_counts = read_layers.get("counts").unwrap();
    let (n_rows, n_cols) = read_counts.shape();
    assert_eq!(n_rows, 2);
    assert_eq!(n_cols, 3);
    
    // 验证 log1p layer
    let read_log1p = read_layers.get("log1p").unwrap();
    let (n_rows, n_cols) = read_log1p.shape();
    assert_eq!(n_rows, 2);
    assert_eq!(n_cols, 3);
    
    // 清理临时文件
    std::fs::remove_file(&temp_path).ok();
    
    println!("✓ Successfully wrote and verified layers");
}

#[test]
fn test_write_pairwise_matrices() {
    use crosscell::ir::PairwiseMatrix;
    use std::collections::HashMap;
    
    // 创建基础数据
    let csr = SparseMatrixCSR::new(
        vec![1.0, 2.0, 3.0, 4.0],
        vec![0, 1, 0, 1],
        vec![0, 2, 4],
        2,
        2,
    ).unwrap();
    let expression = ExpressionMatrix::SparseCSR(csr);
    let cell_metadata = DataFrame::empty(2);
    let gene_metadata = DataFrame::empty(2);
    let metadata = DatasetMetadata::new(2, 2, "test".to_string());
    
    let mut data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata).unwrap();
    
    // 添加细胞-细胞成对矩阵 (obsp)
    let mut cell_pairwise = HashMap::new();
    
    // connectivities (稀疏方阵)
    let conn_csr = SparseMatrixCSR::new(
        vec![1.0, 0.5, 0.5, 1.0],
        vec![0, 1, 0, 1],
        vec![0, 2, 4],
        2,
        2,
    ).unwrap();
    let conn_pw = PairwiseMatrix::new(
        "connectivities".to_string(),
        ExpressionMatrix::SparseCSR(conn_csr),
    ).unwrap();
    cell_pairwise.insert("connectivities".to_string(), conn_pw);
    
    // distances (稠密方阵)
    let dist_dense = DenseMatrix::new(
        vec![0.0, 1.5, 1.5, 0.0],
        2,
        2,
    ).unwrap();
    let dist_pw = PairwiseMatrix::new(
        "distances".to_string(),
        ExpressionMatrix::Dense(dist_dense),
    ).unwrap();
    cell_pairwise.insert("distances".to_string(), dist_pw);
    
    data.cell_pairwise = Some(cell_pairwise);
    
    // 添加基因-基因成对矩阵 (varp)
    let mut gene_pairwise = HashMap::new();
    
    let corr_dense = DenseMatrix::new(
        vec![1.0, 0.8, 0.8, 1.0],
        2,
        2,
    ).unwrap();
    let corr_pw = PairwiseMatrix::new(
        "correlation".to_string(),
        ExpressionMatrix::Dense(corr_dense),
    ).unwrap();
    gene_pairwise.insert("correlation".to_string(), corr_pw);
    
    data.gene_pairwise = Some(gene_pairwise);
    
    // 写入临时文件
    let temp_file = NamedTempFile::new().unwrap();
    let temp_path = temp_file.path().with_extension("h5ad");
    
    let result = write_h5ad(&data, &temp_path);
    assert!(result.is_ok(), "Failed to write h5ad with pairwise matrices: {:?}", result.err());
    
    // 读取并验证
    let read_result = read_h5ad(&temp_path);
    assert!(read_result.is_ok(), "Failed to read written h5ad: {:?}", read_result.err());
    
    let read_data = read_result.unwrap();
    
    // 验证 cell_pairwise 存在
    assert!(read_data.cell_pairwise.is_some());
    let read_cell_pw = read_data.cell_pairwise.unwrap();
    assert_eq!(read_cell_pw.len(), 2);
    
    // 验证 connectivities
    let read_conn = read_cell_pw.get("connectivities").unwrap();
    let (n_rows, n_cols) = read_conn.matrix.shape();
    assert_eq!(n_rows, 2);
    assert_eq!(n_cols, 2);
    
    // 验证 distances
    let read_dist = read_cell_pw.get("distances").unwrap();
    let (n_rows, n_cols) = read_dist.matrix.shape();
    assert_eq!(n_rows, 2);
    assert_eq!(n_cols, 2);
    
    // 验证 gene_pairwise 存在
    assert!(read_data.gene_pairwise.is_some());
    let read_gene_pw = read_data.gene_pairwise.unwrap();
    assert_eq!(read_gene_pw.len(), 1);
    
    // 验证 correlation
    let read_corr = read_gene_pw.get("correlation").unwrap();
    let (n_rows, n_cols) = read_corr.matrix.shape();
    assert_eq!(n_rows, 2);
    assert_eq!(n_cols, 2);
    
    // 清理临时文件
    std::fs::remove_file(&temp_path).ok();
    
    println!("✓ Successfully wrote and verified pairwise matrices");
}

#[test]
fn test_roundtrip_with_embeddings_and_layers() {
    use crosscell::ir::{Embedding, PairwiseMatrix};
    use std::collections::HashMap;
    
    // 创建完整的数据集
    let csr = SparseMatrixCSR::new(
        vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        vec![0, 1, 2, 0, 1, 2],
        vec![0, 3, 6],
        2,
        3,
    ).unwrap();
    let expression = ExpressionMatrix::SparseCSR(csr);
    let cell_metadata = DataFrame::empty(2);
    let gene_metadata = DataFrame::empty(3);
    let metadata = DatasetMetadata::new(2, 3, "test".to_string());
    
    let mut data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata).unwrap();
    
    // 添加嵌入
    let mut embeddings = HashMap::new();
    let pca_data = vec![1.0, 2.0, 3.0, 4.0];
    let pca = Embedding::new("X_pca".to_string(), pca_data, 2, 2).unwrap();
    embeddings.insert("X_pca".to_string(), pca);
    data.embeddings = Some(embeddings);
    
    // 添加 layers
    let mut layers = HashMap::new();
    let counts_csr = SparseMatrixCSR::new(
        vec![10.0, 20.0, 30.0, 40.0, 50.0, 60.0],
        vec![0, 1, 2, 0, 1, 2],
        vec![0, 3, 6],
        2,
        3,
    ).unwrap();
    layers.insert("counts".to_string(), ExpressionMatrix::SparseCSR(counts_csr));
    data.layers = Some(layers);
    
    // 添加成对矩阵
    let mut cell_pairwise = HashMap::new();
    let conn_csr = SparseMatrixCSR::new(
        vec![1.0, 0.5, 0.5, 1.0],
        vec![0, 1, 0, 1],
        vec![0, 2, 4],
        2,
        2,
    ).unwrap();
    let conn_pw = PairwiseMatrix::new(
        "connectivities".to_string(),
        ExpressionMatrix::SparseCSR(conn_csr),
    ).unwrap();
    cell_pairwise.insert("connectivities".to_string(), conn_pw);
    data.cell_pairwise = Some(cell_pairwise);
    
    // 写入临时文件
    let temp_file = NamedTempFile::new().unwrap();
    let temp_path = temp_file.path().with_extension("h5ad");
    
    let result = write_h5ad(&data, &temp_path);
    assert!(result.is_ok(), "Failed to write complete h5ad: {:?}", result.err());
    
    // 读取并验证
    let read_result = read_h5ad(&temp_path);
    assert!(read_result.is_ok(), "Failed to read written h5ad: {:?}", read_result.err());
    
    let read_data = read_result.unwrap();
    
    // 验证所有组件
    assert_eq!(read_data.metadata.n_cells, 2);
    assert_eq!(read_data.metadata.n_genes, 3);
    assert!(read_data.embeddings.is_some());
    assert!(read_data.layers.is_some());
    assert!(read_data.cell_pairwise.is_some());
    
    let read_embeddings = read_data.embeddings.unwrap();
    assert_eq!(read_embeddings.len(), 1);
    assert!(read_embeddings.contains_key("X_pca"));
    
    let read_layers = read_data.layers.unwrap();
    assert_eq!(read_layers.len(), 1);
    assert!(read_layers.contains_key("counts"));
    
    let read_cell_pw = read_data.cell_pairwise.unwrap();
    assert_eq!(read_cell_pw.len(), 1);
    assert!(read_cell_pw.contains_key("connectivities"));
    
    // 清理临时文件
    std::fs::remove_file(&temp_path).ok();
    
    println!("✓ Successfully completed roundtrip test with embeddings, layers, and pairwise matrices");
}

// ============================================================================
// 任务 7.7: 额外的单元测试 - 边界情况、错误处理、特殊值
// ============================================================================

#[test]
fn test_write_single_cell_single_gene() {
    // 边界情况：1 细胞 × 1 基因
    let data = vec![42.0];
    let dense = DenseMatrix::new(data, 1, 1).unwrap();
    
    let expression = ExpressionMatrix::Dense(dense);
    let cell_metadata = DataFrame::empty(1);
    let gene_metadata = DataFrame::empty(1);
    let metadata = DatasetMetadata::new(1, 1, "test".to_string());
    
    let data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata).unwrap();
    
    let temp_file = NamedTempFile::new().unwrap();
    let temp_path = temp_file.path().with_extension("h5ad");
    
    let result = write_h5ad(&data, &temp_path);
    assert!(result.is_ok(), "Failed to write 1×1 h5ad: {:?}", result.err());
    
    let read_result = read_h5ad(&temp_path);
    assert!(read_result.is_ok(), "Failed to read 1×1 h5ad: {:?}", read_result.err());
    
    let read_data = read_result.unwrap();
    assert_eq!(read_data.metadata.n_cells, 1);
    assert_eq!(read_data.metadata.n_genes, 1);
    
    match read_data.expression {
        ExpressionMatrix::Dense(ref d) => {
            assert_eq!(d.data[0], 42.0);
        }
        _ => panic!("Expected Dense matrix"),
    }
    
    std::fs::remove_file(&temp_path).ok();
    println!("✓ Successfully wrote and verified 1×1 matrix");
}

#[test]
fn test_write_single_row_matrix() {
    // 边界情况：1 细胞 × 多基因
    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let dense = DenseMatrix::new(data, 1, 5).unwrap();
    
    let expression = ExpressionMatrix::Dense(dense);
    let cell_metadata = DataFrame::empty(1);
    let gene_metadata = DataFrame::empty(5);
    let metadata = DatasetMetadata::new(1, 5, "test".to_string());
    
    let data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata).unwrap();
    
    let temp_file = NamedTempFile::new().unwrap();
    let temp_path = temp_file.path().with_extension("h5ad");
    
    write_h5ad(&data, &temp_path).unwrap();
    let read_data = read_h5ad(&temp_path).unwrap();
    
    assert_eq!(read_data.metadata.n_cells, 1);
    assert_eq!(read_data.metadata.n_genes, 5);
    
    std::fs::remove_file(&temp_path).ok();
    println!("✓ Successfully wrote and verified single row matrix");
}

#[test]
fn test_write_single_column_matrix() {
    // 边界情况：多细胞 × 1 基因
    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let dense = DenseMatrix::new(data, 5, 1).unwrap();
    
    let expression = ExpressionMatrix::Dense(dense);
    let cell_metadata = DataFrame::empty(5);
    let gene_metadata = DataFrame::empty(1);
    let metadata = DatasetMetadata::new(5, 1, "test".to_string());
    
    let data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata).unwrap();
    
    let temp_file = NamedTempFile::new().unwrap();
    let temp_path = temp_file.path().with_extension("h5ad");
    
    write_h5ad(&data, &temp_path).unwrap();
    let read_data = read_h5ad(&temp_path).unwrap();
    
    assert_eq!(read_data.metadata.n_cells, 5);
    assert_eq!(read_data.metadata.n_genes, 1);
    
    std::fs::remove_file(&temp_path).ok();
    println!("✓ Successfully wrote and verified single column matrix");
}

#[test]
fn test_write_sparse_with_zeros() {
    // 测试包含显式零值的稀疏矩阵
    let data = vec![1.0, 0.0, 2.0, 0.0, 3.0];
    let indices = vec![0, 1, 2, 3, 4];
    let indptr = vec![0, 2, 5];
    let csr = SparseMatrixCSR::new(data, indices, indptr, 2, 5).unwrap();
    
    let expression = ExpressionMatrix::SparseCSR(csr);
    let cell_metadata = DataFrame::empty(2);
    let gene_metadata = DataFrame::empty(5);
    let metadata = DatasetMetadata::new(2, 5, "test".to_string());
    
    let data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata).unwrap();
    
    let temp_file = NamedTempFile::new().unwrap();
    let temp_path = temp_file.path().with_extension("h5ad");
    
    write_h5ad(&data, &temp_path).unwrap();
    let read_data = read_h5ad(&temp_path).unwrap();
    
    match read_data.expression {
        ExpressionMatrix::SparseCSR(ref csr) => {
            assert_eq!(csr.data.len(), 5);
            // 验证零值被保留
            assert_eq!(csr.data[1], 0.0);
            assert_eq!(csr.data[3], 0.0);
        }
        _ => panic!("Expected SparseCSR matrix"),
    }
    
    std::fs::remove_file(&temp_path).ok();
    println!("✓ Successfully wrote and verified sparse matrix with explicit zeros");
}

#[test]
fn test_write_negative_values() {
    // 测试负值
    let data = vec![-1.0, -2.0, -3.0, 4.0, 5.0, 6.0];
    let dense = DenseMatrix::new(data, 2, 3).unwrap();
    
    let expression = ExpressionMatrix::Dense(dense);
    let cell_metadata = DataFrame::empty(2);
    let gene_metadata = DataFrame::empty(3);
    let metadata = DatasetMetadata::new(2, 3, "test".to_string());
    
    let data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata).unwrap();
    
    let temp_file = NamedTempFile::new().unwrap();
    let temp_path = temp_file.path().with_extension("h5ad");
    
    write_h5ad(&data, &temp_path).unwrap();
    let read_data = read_h5ad(&temp_path).unwrap();
    
    match read_data.expression {
        ExpressionMatrix::Dense(ref d) => {
            assert_eq!(d.data[0], -1.0);
            assert_eq!(d.data[1], -2.0);
            assert_eq!(d.data[2], -3.0);
        }
        _ => panic!("Expected Dense matrix"),
    }
    
    std::fs::remove_file(&temp_path).ok();
    println!("✓ Successfully wrote and verified negative values");
}

#[test]
fn test_write_very_sparse_matrix() {
    // 测试极高稀疏度矩阵（99.9% 稀疏）
    // 1000 细胞 × 500 基因，只有 500 个非零元素
    let data: Vec<f64> = (1..=500).map(|x| x as f64).collect();
    let indices: Vec<usize> = (0..500).collect();
    let mut indptr: Vec<usize> = vec![0];
    for i in 1..=1000 {
        if i <= 500 {
            indptr.push(i);
        } else {
            indptr.push(500);
        }
    }
    
    let csr = SparseMatrixCSR::new(data, indices, indptr, 1000, 500).unwrap();
    let expression = ExpressionMatrix::SparseCSR(csr);
    let cell_metadata = DataFrame::empty(1000);
    let gene_metadata = DataFrame::empty(500);
    let metadata = DatasetMetadata::new(1000, 500, "test".to_string());
    
    let data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata).unwrap();
    
    let temp_file = NamedTempFile::new().unwrap();
    let temp_path = temp_file.path().with_extension("h5ad");
    
    write_h5ad(&data, &temp_path).unwrap();
    let read_data = read_h5ad(&temp_path).unwrap();
    
    assert_eq!(read_data.metadata.n_cells, 1000);
    assert_eq!(read_data.metadata.n_genes, 500);
    
    match read_data.expression {
        ExpressionMatrix::SparseCSR(ref csr) => {
            assert_eq!(csr.data.len(), 500);
            // 验证稀疏度保留（修复：使用 >= 而不是 >）
            let sparsity = 1.0 - (csr.data.len() as f64) / (1000.0 * 500.0);
            assert!(sparsity >= 0.999, "Sparsity should be >= 99.9%, got {}", sparsity);
        }
        _ => panic!("Expected SparseCSR matrix"),
    }
    
    std::fs::remove_file(&temp_path).ok();
    println!("✓ Successfully wrote and verified very sparse matrix (99.9% sparse)");
}

#[test]
fn test_write_large_values() {
    // 测试大数值
    let data = vec![1e10, 1e15, 1e20, 1e-10, 1e-15, 1e-20];
    let dense = DenseMatrix::new(data, 2, 3).unwrap();
    
    let expression = ExpressionMatrix::Dense(dense);
    let cell_metadata = DataFrame::empty(2);
    let gene_metadata = DataFrame::empty(3);
    let metadata = DatasetMetadata::new(2, 3, "test".to_string());
    
    let data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata).unwrap();
    
    let temp_file = NamedTempFile::new().unwrap();
    let temp_path = temp_file.path().with_extension("h5ad");
    
    write_h5ad(&data, &temp_path).unwrap();
    let read_data = read_h5ad(&temp_path).unwrap();
    
    match read_data.expression {
        ExpressionMatrix::Dense(ref d) => {
            // 验证大数值保留精度
            assert!((d.data[0] - 1e10).abs() < 1e3);
            assert!((d.data[1] - 1e15).abs() < 1e8);
            assert!((d.data[3] - 1e-10).abs() < 1e-17);
        }
        _ => panic!("Expected Dense matrix"),
    }
    
    std::fs::remove_file(&temp_path).ok();
    println!("✓ Successfully wrote and verified large values");
}

#[test]
fn test_write_empty_embeddings() {
    // 测试空嵌入 HashMap
    use std::collections::HashMap;
    
    let csr = SparseMatrixCSR::new(vec![1.0], vec![0], vec![0, 1], 1, 1).unwrap();
    let expression = ExpressionMatrix::SparseCSR(csr);
    let cell_metadata = DataFrame::empty(1);
    let gene_metadata = DataFrame::empty(1);
    let metadata = DatasetMetadata::new(1, 1, "test".to_string());
    
    let mut data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata).unwrap();
    data.embeddings = Some(HashMap::new());
    
    let temp_file = NamedTempFile::new().unwrap();
    let temp_path = temp_file.path().with_extension("h5ad");
    
    write_h5ad(&data, &temp_path).unwrap();
    let read_data = read_h5ad(&temp_path).unwrap();
    
    // 空 embeddings 应该被读取为 Some(empty HashMap)
    assert!(read_data.embeddings.is_some());
    assert_eq!(read_data.embeddings.unwrap().len(), 0);
    
    std::fs::remove_file(&temp_path).ok();
    println!("✓ Successfully wrote and verified empty embeddings");
}

#[test]
fn test_write_empty_layers() {
    // 测试空 layers HashMap
    // 空 HashMap 写入后不会创建 /layers group，读回来是 None
    // 这是合理行为：Some(空HashMap) 和 None 语义等价
    use std::collections::HashMap;
    
    let csr = SparseMatrixCSR::new(vec![1.0], vec![0], vec![0, 1], 1, 1).unwrap();
    let expression = ExpressionMatrix::SparseCSR(csr);
    let cell_metadata = DataFrame::empty(1);
    let gene_metadata = DataFrame::empty(1);
    let metadata = DatasetMetadata::new(1, 1, "test".to_string());
    
    let mut data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata).unwrap();
    data.layers = Some(HashMap::new());
    
    let temp_file = NamedTempFile::new().unwrap();
    let temp_path = temp_file.path().with_extension("h5ad");
    
    write_h5ad(&data, &temp_path).unwrap();
    let read_data = read_h5ad(&temp_path).unwrap();
    
    // 空 layers 往返后变为 None（writer 不写空 /layers group）
    assert!(read_data.layers.is_none() || read_data.layers.as_ref().map_or(true, |l| l.is_empty()));
    
    std::fs::remove_file(&temp_path).ok();
    println!("✓ Successfully wrote and verified empty layers");
}

#[test]
fn test_write_mixed_sparse_dense_layers() {
    // 测试混合稀疏和稠密 layers
    use std::collections::HashMap;
    
    let csr = SparseMatrixCSR::new(
        vec![1.0, 2.0],
        vec![0, 1],
        vec![0, 1, 2],
        2,
        2,
    ).unwrap();
    let expression = ExpressionMatrix::SparseCSR(csr);
    let cell_metadata = DataFrame::empty(2);
    let gene_metadata = DataFrame::empty(2);
    let metadata = DatasetMetadata::new(2, 2, "test".to_string());
    
    let mut data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata).unwrap();
    
    let mut layers = HashMap::new();
    
    // 稀疏 layer
    let sparse_csr = SparseMatrixCSR::new(
        vec![10.0, 20.0],
        vec![0, 1],
        vec![0, 1, 2],
        2,
        2,
    ).unwrap();
    layers.insert("sparse_layer".to_string(), ExpressionMatrix::SparseCSR(sparse_csr));
    
    // 稠密 layer
    let dense = DenseMatrix::new(vec![1.0, 2.0, 3.0, 4.0], 2, 2).unwrap();
    layers.insert("dense_layer".to_string(), ExpressionMatrix::Dense(dense));
    
    data.layers = Some(layers);
    
    let temp_file = NamedTempFile::new().unwrap();
    let temp_path = temp_file.path().with_extension("h5ad");
    
    write_h5ad(&data, &temp_path).unwrap();
    let read_data = read_h5ad(&temp_path).unwrap();
    
    assert!(read_data.layers.is_some());
    let read_layers = read_data.layers.unwrap();
    assert_eq!(read_layers.len(), 2);
    assert!(read_layers.contains_key("sparse_layer"));
    assert!(read_layers.contains_key("dense_layer"));
    
    std::fs::remove_file(&temp_path).ok();
    println!("✓ Successfully wrote and verified mixed sparse/dense layers");
}

#[test]
fn test_roundtrip_consistency_numerical() {
    // 测试往返数值一致性
    let data = vec![1.23456789, 9.87654321, 3.14159265, 2.71828182];
    let dense = DenseMatrix::new(data.clone(), 2, 2).unwrap();
    
    let expression = ExpressionMatrix::Dense(dense);
    let cell_metadata = DataFrame::empty(2);
    let gene_metadata = DataFrame::empty(2);
    let metadata = DatasetMetadata::new(2, 2, "test".to_string());
    
    let original_data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata).unwrap();
    
    let temp_file = NamedTempFile::new().unwrap();
    let temp_path = temp_file.path().with_extension("h5ad");
    
    write_h5ad(&original_data, &temp_path).unwrap();
    let read_data = read_h5ad(&temp_path).unwrap();
    
    match read_data.expression {
        ExpressionMatrix::Dense(ref d) => {
            for (i, &original_val) in data.iter().enumerate() {
                let read_val = d.data[i];
                let rel_error = ((read_val - original_val) / original_val).abs();
                assert!(rel_error < 1e-10, 
                    "Value mismatch at index {}: original={}, read={}, rel_error={}",
                    i, original_val, read_val, rel_error);
            }
        }
        _ => panic!("Expected Dense matrix"),
    }
    
    std::fs::remove_file(&temp_path).ok();
    println!("✓ Successfully verified numerical consistency in roundtrip");
}

#[test]
fn test_write_multiple_embeddings() {
    // 测试多个不同维度的嵌入
    use crosscell::ir::Embedding;
    use std::collections::HashMap;
    
    // 修复：创建 3 细胞 × 1 基因的矩阵，与 metadata 一致
    let csr = SparseMatrixCSR::new(vec![1.0, 2.0, 3.0], vec![0, 0, 0], vec![0, 1, 2, 3], 3, 1).unwrap();
    let expression = ExpressionMatrix::SparseCSR(csr);
    let cell_metadata = DataFrame::empty(3);
    let gene_metadata = DataFrame::empty(1);
    let metadata = DatasetMetadata::new(3, 1, "test".to_string());
    
    let mut data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata).unwrap();
    
    let mut embeddings = HashMap::new();
    
    // PCA: 3 细胞 × 50 成分
    let pca_data: Vec<f64> = (0..150).map(|x| x as f64).collect();
    let pca = Embedding::new("X_pca".to_string(), pca_data, 3, 50).unwrap();
    embeddings.insert("X_pca".to_string(), pca);
    
    // UMAP: 3 细胞 × 2 成分
    let umap_data = vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
    let umap = Embedding::new("X_umap".to_string(), umap_data, 3, 2).unwrap();
    embeddings.insert("X_umap".to_string(), umap);
    
    // tSNE: 3 细胞 × 2 成分
    let tsne_data = vec![1.1, 1.2, 1.3, 1.4, 1.5, 1.6];
    let tsne = Embedding::new("X_tsne".to_string(), tsne_data, 3, 2).unwrap();
    embeddings.insert("X_tsne".to_string(), tsne);
    
    data.embeddings = Some(embeddings);
    
    let temp_file = NamedTempFile::new().unwrap();
    let temp_path = temp_file.path().with_extension("h5ad");
    
    write_h5ad(&data, &temp_path).unwrap();
    let read_data = read_h5ad(&temp_path).unwrap();
    
    assert!(read_data.embeddings.is_some());
    let read_embeddings = read_data.embeddings.unwrap();
    assert_eq!(read_embeddings.len(), 3);
    
    // 验证 PCA
    let pca = read_embeddings.get("X_pca").unwrap();
    assert_eq!(pca.n_rows, 3);
    assert_eq!(pca.n_cols, 50);
    
    // 验证 UMAP
    let umap = read_embeddings.get("X_umap").unwrap();
    assert_eq!(umap.n_rows, 3);
    assert_eq!(umap.n_cols, 2);
    
    // 验证 tSNE
    let tsne = read_embeddings.get("X_tsne").unwrap();
    assert_eq!(tsne.n_rows, 3);
    assert_eq!(tsne.n_cols, 2);
    
    std::fs::remove_file(&temp_path).ok();
    println!("✓ Successfully wrote and verified multiple embeddings with different dimensions");
}
