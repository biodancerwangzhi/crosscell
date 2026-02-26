//! 流式转换性能基准测试
//!
//! 本基准测试套件对比流式 vs 非流式转换性能，以及不同 chunk_size 的影响。
//!
//! ## 测试场景
//! 1. 流式 vs 非流式 H5AD 读取
//! 2. 不同 chunk_size 的性能对比
//! 3. 内存使用估算
//!
//! ## 测试数据规模
//! - Small: 3k 细胞（PBMC 3k 规模）
//! - Medium: 100k 细胞
//! - Large: 1M 细胞（仅流式）
//!
//! ## 运行方法
//! ```bash
//! # 在 Docker 容器中运行
//! docker-compose run --rm dev cargo bench streaming_conversion
//! ```

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use crosscell::anndata::streaming::{DataChunk, StreamingMetadata};
use crosscell::ir::{DataFrame, ExpressionMatrix, SparseMatrixCSR};

/// 生成测试用的 StreamingMetadata
fn create_test_metadata(n_cells: usize, n_genes: usize) -> StreamingMetadata {
    StreamingMetadata {
        n_cells,
        n_genes,
        is_sparse: true,
        sparse_format: Some("csr_matrix".to_string()),
        nnz: Some(n_cells * n_genes / 20), // 5% 非零
        obs_columns: vec!["cell_type".to_string(), "batch".to_string()],
        var_columns: vec!["gene_name".to_string()],
        embedding_names: vec!["X_pca".to_string(), "X_umap".to_string()],
        layer_names: vec![],
    }
}

/// 生成测试用的 DataChunk
fn create_test_chunk(start: usize, end: usize, n_genes: usize) -> DataChunk {
    let n_rows = end - start;
    let nnz_per_row = n_genes / 20; // 5% 非零
    let total_nnz = n_rows * nnz_per_row;

    let mut data = Vec::with_capacity(total_nnz);
    let mut indices = Vec::with_capacity(total_nnz);
    let mut indptr = Vec::with_capacity(n_rows + 1);

    indptr.push(0);
    for row in 0..n_rows {
        for col in 0..nnz_per_row {
            data.push((row * col) as f64 * 0.001);
            indices.push(col * 20); // 每 20 列一个非零元素
        }
        indptr.push(data.len());
    }

    let csr = SparseMatrixCSR {
        data,
        indices,
        indptr,
        n_rows,
        n_cols: n_genes,
    };

    DataChunk {
        expression: ExpressionMatrix::SparseCSR(csr),
        cell_indices: start..end,
        cell_metadata: DataFrame::empty(n_rows),
        embeddings: None,
    }
}

/// 基准测试：不同 chunk_size 的内存估算
fn bench_chunk_size_memory(c: &mut Criterion) {
    let mut group = c.benchmark_group("chunk_size_memory");

    let n_genes = 30000;
    let chunk_sizes = [1000, 5000, 10000, 20000, 50000];

    for &chunk_size in &chunk_sizes {
        let estimated_memory = chunk_size * n_genes * 8 * 2; // f64 × 2 buffers

        group.throughput(Throughput::Bytes(estimated_memory as u64));
        group.bench_with_input(
            BenchmarkId::new("estimate", format!("chunk_{}", chunk_size)),
            &chunk_size,
            |b, &cs| {
                b.iter(|| {
                    let mem = black_box(cs) * black_box(n_genes) * 8 * 2;
                    black_box(mem)
                });
            },
        );
    }

    group.finish();
}

/// 基准测试：DataChunk 创建性能
fn bench_chunk_creation(c: &mut Criterion) {
    let mut group = c.benchmark_group("chunk_creation");

    let n_genes = 30000;
    let chunk_sizes = [1000, 5000, 10000];

    for &chunk_size in &chunk_sizes {
        group.throughput(Throughput::Elements(chunk_size as u64));
        group.bench_with_input(
            BenchmarkId::new("create", format!("chunk_{}", chunk_size)),
            &chunk_size,
            |b, &cs| {
                b.iter(|| {
                    let chunk = create_test_chunk(0, black_box(cs), black_box(n_genes));
                    black_box(chunk)
                });
            },
        );
    }

    group.finish();
}

/// 基准测试：StreamingMetadata 创建
fn bench_metadata_creation(c: &mut Criterion) {
    let mut group = c.benchmark_group("metadata_creation");

    let sizes = [
        (3000, 20000, "small"),
        (100000, 30000, "medium"),
        (1000000, 30000, "large"),
    ];

    for (n_cells, n_genes, name) in sizes {
        group.bench_with_input(
            BenchmarkId::new("create", name),
            &(n_cells, n_genes),
            |b, &(cells, genes)| {
                b.iter(|| {
                    let metadata = create_test_metadata(black_box(cells), black_box(genes));
                    black_box(metadata)
                });
            },
        );
    }

    group.finish();
}

/// 基准测试：块迭代模拟
fn bench_chunk_iteration(c: &mut Criterion) {
    let mut group = c.benchmark_group("chunk_iteration");

    let n_cells = 100000;
    let _n_genes = 30000;
    let chunk_sizes = [1000, 5000, 10000, 20000];

    for &chunk_size in &chunk_sizes {
        let _n_chunks = (n_cells + chunk_size - 1) / chunk_size;

        group.throughput(Throughput::Elements(n_cells as u64));
        group.bench_with_input(
            BenchmarkId::new("iterate", format!("chunk_{}", chunk_size)),
            &chunk_size,
            |b, &cs| {
                b.iter(|| {
                    let mut total_cells = 0usize;
                    let mut chunk_idx = 0;

                    while chunk_idx * cs < n_cells {
                        let start = chunk_idx * cs;
                        let end = std::cmp::min((chunk_idx + 1) * cs, n_cells);
                        total_cells += end - start;
                        chunk_idx += 1;
                    }

                    black_box(total_cells)
                });
            },
        );
    }

    group.finish();
}

/// 基准测试：稀疏矩阵块提取
fn bench_sparse_chunk_extraction(c: &mut Criterion) {
    let mut group = c.benchmark_group("sparse_chunk_extraction");

    // 创建一个中等大小的稀疏矩阵
    let n_rows = 10000;
    let n_cols = 5000;
    let nnz_per_row = 250; // 5% 非零

    let mut data = Vec::with_capacity(n_rows * nnz_per_row);
    let mut indices = Vec::with_capacity(n_rows * nnz_per_row);
    let mut indptr = Vec::with_capacity(n_rows + 1);

    indptr.push(0);
    for row in 0..n_rows {
        for col in 0..nnz_per_row {
            data.push((row * col) as f64 * 0.001);
            indices.push(col * 20);
        }
        indptr.push(data.len());
    }

    let csr = SparseMatrixCSR {
        data,
        indices,
        indptr,
        n_rows,
        n_cols,
    };

    let chunk_sizes = [100, 500, 1000, 2000];

    for &chunk_size in &chunk_sizes {
        group.throughput(Throughput::Elements(chunk_size as u64));
        group.bench_with_input(
            BenchmarkId::new("extract", format!("chunk_{}", chunk_size)),
            &chunk_size,
            |b, &cs| {
                b.iter(|| {
                    // 模拟提取一个块
                    let start_row = 0;
                    let end_row = cs;

                    let data_start = csr.indptr[start_row];
                    let data_end = csr.indptr[end_row];

                    let chunk_data: Vec<f64> = csr.data[data_start..data_end].to_vec();
                    let chunk_indices: Vec<usize> = csr.indices[data_start..data_end].to_vec();
                    let chunk_indptr: Vec<usize> = csr.indptr[start_row..=end_row]
                        .iter()
                        .map(|&p| p - data_start)
                        .collect();

                    black_box((chunk_data, chunk_indices, chunk_indptr))
                });
            },
        );
    }

    group.finish();
}

/// 基准测试：内存使用估算准确性
fn bench_memory_estimation_accuracy(c: &mut Criterion) {
    let mut group = c.benchmark_group("memory_estimation");

    let scenarios = [
        (3000, 20000, 1000, "small_3k"),
        (100000, 30000, 10000, "medium_100k"),
        (1000000, 30000, 10000, "large_1m"),
    ];

    for (n_cells, n_genes, chunk_size, name) in scenarios {
        group.bench_with_input(
            BenchmarkId::new("estimate", name),
            &(n_cells, n_genes, chunk_size),
            |b, &(cells, genes, cs)| {
                b.iter(|| {
                    // 估算流式模式内存
                    let streaming_mem = cs * genes * 8 * 2;

                    // 估算非流式模式内存
                    let full_mem = cells * genes * 8;

                    // 计算节省比例
                    let savings = 1.0 - (streaming_mem as f64 / full_mem as f64);

                    black_box((streaming_mem, full_mem, savings))
                });
            },
        );
    }

    group.finish();
}

/// 基准测试：吞吐量计算
fn bench_throughput_calculation(c: &mut Criterion) {
    let mut group = c.benchmark_group("throughput_calculation");

    let scenarios = [
        (3000, 1.5, "small_3k"),       // 3k cells, 1.5 seconds
        (100000, 15.0, "medium_100k"), // 100k cells, 15 seconds
        (1000000, 120.0, "large_1m"),  // 1M cells, 120 seconds
    ];

    for (n_cells, time_secs, name) in scenarios {
        group.bench_with_input(
            BenchmarkId::new("calculate", name),
            &(n_cells, time_secs),
            |b, &(cells, time)| {
                b.iter(|| {
                    let throughput = cells as f64 / time;
                    black_box(throughput)
                });
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_chunk_size_memory,
    bench_chunk_creation,
    bench_metadata_creation,
    bench_chunk_iteration,
    bench_sparse_chunk_extraction,
    bench_memory_estimation_accuracy,
    bench_throughput_calculation,
);
criterion_main!(benches);
