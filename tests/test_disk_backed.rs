//! 磁盘支持矩阵的集成测试

use crosscell::ir::expression::{ExpressionMatrix, SparseFormat, SparseMatrixCSC, SparseMatrixCSR};
use crosscell::storage::{ChunkStore, DiskBackedConfig, DiskBackedMatrix};
use tempfile::tempdir;

#[test]
fn test_chunk_store_create_and_open() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("test.chunks");

    // 创建存储
    let store = ChunkStore::create(&path, 1000, 500, 100, SparseFormat::CSR, true).unwrap();

    assert_eq!(store.total_rows, 1000);
    assert_eq!(store.total_cols, 500);
    assert_eq!(store.chunk_size, 100);
    assert_eq!(store.num_chunks(), 10);

    // 重新打开
    let store2 = ChunkStore::open(&path).unwrap();
    assert_eq!(store2.total_rows, 1000);
    assert_eq!(store2.total_cols, 500);
}

#[test]
fn test_chunk_store_csr_write_read() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("test_csr.chunks");

    // 创建存储
    let mut store = ChunkStore::create(&path, 100, 50, 25, SparseFormat::CSR, true).unwrap();

    // 创建测试矩阵 (5 行 × 50 列)
    let matrix = SparseMatrixCSR::new(
        vec![1.0, 2.0, 3.0, 4.0, 5.0],
        vec![0, 10, 20, 30, 40],
        vec![0, 1, 2, 3, 4, 5],
        5,
        50,
    )
    .unwrap();

    // 写入
    store.write_csr_chunk(0, &matrix).unwrap();

    // 读取
    let matrix2 = store.read_csr_chunk(0).unwrap();

    assert_eq!(matrix.data, matrix2.data);
    assert_eq!(matrix.indices, matrix2.indices);
    assert_eq!(matrix.indptr, matrix2.indptr);
    assert_eq!(matrix.n_rows, matrix2.n_rows);
    assert_eq!(matrix.n_cols, matrix2.n_cols);
}

#[test]
fn test_chunk_store_csc_write_read() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("test_csc.chunks");

    // 创建存储
    let mut store = ChunkStore::create(&path, 100, 50, 25, SparseFormat::CSC, true).unwrap();

    // 创建测试矩阵 (100 行 × 5 列)
    let matrix = SparseMatrixCSC::new(
        vec![1.0, 2.0, 3.0, 4.0, 5.0],
        vec![0, 10, 20, 30, 40],
        vec![0, 1, 2, 3, 4, 5],
        100,
        5,
    )
    .unwrap();

    // 写入
    store.write_csc_chunk(0, &matrix).unwrap();

    // 读取
    let matrix2 = store.read_csc_chunk(0).unwrap();

    assert_eq!(matrix.data, matrix2.data);
    assert_eq!(matrix.indices, matrix2.indices);
    assert_eq!(matrix.indptr, matrix2.indptr);
}

#[test]
fn test_chunk_store_multiple_chunks() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("test_multi.chunks");

    // 创建存储 (100 行，每块 25 行 = 4 块)
    let mut store = ChunkStore::create(&path, 100, 50, 25, SparseFormat::CSR, true).unwrap();

    // 写入 4 个块
    for i in 0..4 {
        let matrix =
            SparseMatrixCSR::new(vec![i as f64 + 1.0], vec![i], vec![0, 1], 1, 50).unwrap();
        store.write_csr_chunk(i, &matrix).unwrap();
    }

    // 验证
    assert_eq!(store.written_chunks(), 4);

    // 读取并验证
    for i in 0..4 {
        let matrix = store.read_csr_chunk(i).unwrap();
        assert_eq!(matrix.data[0], i as f64 + 1.0);
    }
}

#[test]
fn test_chunk_store_compression() {
    let dir = tempdir().unwrap();
    let path_compressed = dir.path().join("compressed.chunks");
    let path_uncompressed = dir.path().join("uncompressed.chunks");

    // 创建大矩阵
    let n_rows = 100;
    let n_cols = 1000;
    let data: Vec<f64> = (0..500).map(|i| i as f64).collect();
    let indices: Vec<usize> = (0..500).map(|i| i % n_cols).collect();
    let mut indptr = vec![0usize];
    for i in 1..=n_rows {
        indptr.push((i * 5).min(500));
    }

    let matrix = SparseMatrixCSR::new(data, indices, indptr, n_rows, n_cols).unwrap();

    // 写入压缩版本
    let mut store_c = ChunkStore::create(
        &path_compressed,
        n_rows,
        n_cols,
        n_rows,
        SparseFormat::CSR,
        true,
    )
    .unwrap();
    store_c.write_csr_chunk(0, &matrix).unwrap();

    // 写入未压缩版本
    let mut store_u = ChunkStore::create(
        &path_uncompressed,
        n_rows,
        n_cols,
        n_rows,
        SparseFormat::CSR,
        false,
    )
    .unwrap();
    store_u.write_csr_chunk(0, &matrix).unwrap();

    // 压缩版本应该更小
    let size_c = std::fs::metadata(&path_compressed).unwrap().len();
    let size_u = std::fs::metadata(&path_uncompressed).unwrap().len();

    println!("Compressed: {} bytes", size_c);
    println!("Uncompressed: {} bytes", size_u);

    // 压缩应该有效果（至少减少 10%）
    assert!(size_c < size_u, "Compression should reduce file size");

    // 验证数据完整性
    let matrix_c = store_c.read_csr_chunk(0).unwrap();
    let matrix_u = store_u.read_csr_chunk(0).unwrap();

    assert_eq!(matrix.data, matrix_c.data);
    assert_eq!(matrix.data, matrix_u.data);
}

#[test]
fn test_disk_backed_matrix_from_csr() {
    let dir = tempdir().unwrap();

    // 创建测试矩阵 (100 行 × 50 列)
    let n_rows = 100;
    let n_cols = 50;
    let data: Vec<f64> = (0..200).map(|i| i as f64).collect();
    let indices: Vec<usize> = (0..200).map(|i| i % n_cols).collect();
    let mut indptr = vec![0usize];
    for i in 1..=n_rows {
        indptr.push(i * 2);
    }

    let matrix = SparseMatrixCSR::new(
        data.clone(),
        indices.clone(),
        indptr.clone(),
        n_rows,
        n_cols,
    )
    .unwrap();
    let expr = ExpressionMatrix::SparseCSR(matrix);

    // 创建磁盘支持矩阵
    let config = DiskBackedConfig::new()
        .with_temp_dir(dir.path())
        .with_chunk_size(25)
        .with_compression(true)
        .with_cleanup(false); // 不自动清理，方便测试

    let disk_matrix = DiskBackedMatrix::from_matrix(&expr, config).unwrap();

    assert_eq!(disk_matrix.shape(), (100, 50));
    assert_eq!(disk_matrix.num_chunks(), 4); // 100 / 25 = 4

    // 读取每个块并验证
    for chunk_idx in 0..4 {
        let chunk = disk_matrix.read_csr_chunk(chunk_idx).unwrap();
        assert_eq!(chunk.n_rows, 25);
        assert_eq!(chunk.n_cols, 50);
    }

    // 转换回内存矩阵
    let restored = disk_matrix.to_matrix().unwrap();

    if let ExpressionMatrix::SparseCSR(m) = restored {
        assert_eq!(m.n_rows, n_rows);
        assert_eq!(m.n_cols, n_cols);
        assert_eq!(m.data, data);
        assert_eq!(m.indices, indices);
        assert_eq!(m.indptr, indptr);
    } else {
        panic!("Expected SparseCSR");
    }
}

#[test]
fn test_disk_backed_matrix_process_chunks() {
    let dir = tempdir().unwrap();

    // 创建测试矩阵
    let n_rows = 100;
    let n_cols = 50;
    let data: Vec<f64> = (0..200).map(|i| i as f64).collect();
    let indices: Vec<usize> = (0..200).map(|i| i % n_cols).collect();
    let mut indptr = vec![0usize];
    for i in 1..=n_rows {
        indptr.push(i * 2);
    }

    let matrix = SparseMatrixCSR::new(data, indices, indptr, n_rows, n_cols).unwrap();
    let expr = ExpressionMatrix::SparseCSR(matrix);

    let config = DiskBackedConfig::new()
        .with_temp_dir(dir.path())
        .with_chunk_size(25)
        .with_cleanup(false);

    let disk_matrix = DiskBackedMatrix::from_matrix(&expr, config).unwrap();

    // 流式处理每个块
    let results = disk_matrix
        .process_chunks(|idx, chunk| {
            let (rows, cols) = chunk.shape();
            Ok((idx, rows, cols))
        })
        .unwrap();

    assert_eq!(results.len(), 4);
    for (idx, rows, cols) in results {
        assert_eq!(rows, 25);
        assert_eq!(cols, 50);
        assert!(idx < 4);
    }
}

#[test]
fn test_disk_backed_matrix_memory_estimate() {
    let dir = tempdir().unwrap();

    let config = DiskBackedConfig::new()
        .with_temp_dir(dir.path())
        .with_chunk_size(5000);

    let disk_matrix = DiskBackedMatrix::create(
        1_000_000, // 1M 细胞
        30_000,    // 30k 基因
        SparseFormat::CSR,
        config,
    )
    .unwrap();

    let chunk_memory = disk_matrix.estimate_chunk_memory();

    // 单个块的内存应该远小于完整矩阵
    // 5000 × 30000 × 5% 稀疏度 × 16 bytes ≈ 120 MB
    println!("Estimated chunk memory: {} bytes", chunk_memory);
    assert!(chunk_memory < 200 * 1024 * 1024); // < 200 MB
}
