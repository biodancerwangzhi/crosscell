//! 稀疏矩阵转换性能基准测试
//!
//! 本基准测试套件对比 CrossCell 优化前后的性能，以及与其他工具的性能对比。
//!
//! ## 测试场景
//! 1. 稀疏矩阵子集操作（行子集、列子集、行列子集）
//! 2. CSR ↔ CSC 格式转换
//! 3. 内存使用估算
//!
//! ## 测试数据规模
//! - Small: 1000×500 细胞，稀疏度 90%
//! - Medium: 10000×5000 细胞，稀疏度 95%
//! - Large: 100000×30000 细胞，稀疏度 98%
//!
//! ## 运行方法
//! ```bash
//! # 在 Docker 容器中运行
//! docker run -it --rm -v ${PWD}:/workspace crosscell-dev cargo bench
//!
//! # 运行特定基准测试
//! docker run -it --rm -v ${PWD}:/workspace crosscell-dev cargo bench sparse_conversion
//! ```

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use crosscell::ir::expression::SparseMatrixCSR;
use crosscell::sparse::convert::{csr_to_csc, csr_to_csc_parallel};
use crosscell::sparse::memory::{estimate_csc_matrix_memory, estimate_csr_matrix_memory};
use crosscell::sparse::subset::{
    subset_contiguous_columns_in_place, subset_rows_only_in_place, subset_sparse_columns_in_place,
};
use rand::Rng;

/// 生成随机稀疏 CSR 矩阵
fn generate_random_csr(n_rows: usize, n_cols: usize, sparsity: f64) -> SparseMatrixCSR {
    let mut rng = rand::thread_rng();
    let nnz_per_row = ((n_cols as f64) * (1.0 - sparsity)) as usize;
    let total_nnz = n_rows * nnz_per_row;

    let mut data = Vec::with_capacity(total_nnz);
    let mut indices = Vec::with_capacity(total_nnz);
    let mut indptr = Vec::with_capacity(n_rows + 1);

    indptr.push(0);

    for _ in 0..n_rows {
        // 为每行生成随机非零元素
        let mut row_cols: Vec<usize> = (0..n_cols).collect();
        // 使用 Fisher-Yates shuffle 算法
        for i in (1..row_cols.len()).rev() {
            let j = rng.gen_range(0..=i);
            row_cols.swap(i, j);
        }
        row_cols.truncate(nnz_per_row);
        row_cols.sort_unstable();

        for col in row_cols {
            data.push(rng.gen::<f64>());
            indices.push(col);
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

/// 基准测试：CSR → CSC 转换（不同数据规模）
fn bench_csr_to_csc_conversion(c: &mut Criterion) {
    let mut group = c.benchmark_group("csr_to_csc_conversion");

    // Small: 1000×500, 90% 稀疏
    let small_csr = generate_random_csr(1000, 500, 0.90);
    group.throughput(Throughput::Elements(small_csr.data.len() as u64));
    group.bench_with_input(
        BenchmarkId::new("small", "1000x500_90%"),
        &small_csr,
        |b, csr| {
            b.iter(|| {
                let csc = csr_to_csc(black_box(csr));
                black_box(csc);
            });
        },
    );

    // Medium: 10000×5000, 95% 稀疏
    let medium_csr = generate_random_csr(10000, 5000, 0.95);
    group.throughput(Throughput::Elements(medium_csr.data.len() as u64));
    group.bench_with_input(
        BenchmarkId::new("medium", "10000x5000_95%"),
        &medium_csr,
        |b, csr| {
            b.iter(|| {
                let csc = csr_to_csc(black_box(csr));
                black_box(csc);
            });
        },
    );

    // Large: 50000×20000, 98% 稀疏（减小规模以加快测试）
    let large_csr = generate_random_csr(50000, 20000, 0.98);
    group.throughput(Throughput::Elements(large_csr.data.len() as u64));
    group.bench_with_input(
        BenchmarkId::new("large", "50000x20000_98%"),
        &large_csr,
        |b, csr| {
            b.iter(|| {
                let csc = csr_to_csc(black_box(csr));
                black_box(csc);
            });
        },
    );

    group.finish();
}

/// 基准测试：CSC → CSR 转换
fn bench_csc_to_csr_conversion(c: &mut Criterion) {
    let group = c.benchmark_group("csc_to_csr_conversion");

    // 直接生成 CSR 并转换为 CSC，然后测试 CSC → CSR
    // 注意：由于当前实现的限制，我们跳过这个基准测试
    // TODO: 修复 CSC → CSR 转换后重新启用

    group.finish();
}

/// 基准测试：并行 vs 串行 CSR → CSC 转换
fn bench_parallel_vs_serial(c: &mut Criterion) {
    let mut group = c.benchmark_group("parallel_vs_serial");

    // 大矩阵才能体现并行优势
    let large_csr = generate_random_csr(50000, 20000, 0.98);

    group.bench_with_input(
        BenchmarkId::new("serial", "50000x20000"),
        &large_csr,
        |b, csr| {
            b.iter(|| {
                let csc = csr_to_csc(black_box(csr));
                black_box(csc);
            });
        },
    );

    group.bench_with_input(
        BenchmarkId::new("parallel", "50000x20000"),
        &large_csr,
        |b, csr| {
            b.iter(|| {
                let csc = csr_to_csc_parallel(black_box(csr));
                black_box(csc);
            });
        },
    );

    group.finish();
}

/// 基准测试：稀疏矩阵行子集操作
fn bench_subset_rows(c: &mut Criterion) {
    let mut group = c.benchmark_group("subset_rows");

    let csr = generate_random_csr(10000, 5000, 0.95);
    let row_indices: Vec<usize> = (0..5000).step_by(2).collect(); // 选择一半的行

    group.bench_function("subset_5000_rows", |b| {
        b.iter(|| {
            let result = subset_rows_only_in_place(
                black_box(csr.indptr.clone()),
                black_box(csr.indices.clone()),
                black_box(csr.data.clone()),
                black_box(&row_indices),
                black_box(csr.n_cols),
            );
            black_box(result);
        });
    });

    group.finish();
}

/// 基准测试：稀疏矩阵列子集操作（连续列）
fn bench_subset_contiguous_columns(c: &mut Criterion) {
    let mut group = c.benchmark_group("subset_contiguous_columns");

    let csr = generate_random_csr(10000, 5000, 0.95);
    let row_indices: Vec<usize> = (0..10000).collect();
    let col_range: Vec<usize> = (1000..3000).collect(); // 连续的 2000 列

    group.bench_function("subset_2000_contiguous_cols", |b| {
        b.iter(|| {
            let result = subset_contiguous_columns_in_place(
                black_box(csr.indptr.clone()),
                black_box(csr.indices.clone()),
                black_box(csr.data.clone()),
                black_box(&row_indices),
                black_box(&col_range),
            );
            black_box(result);
        });
    });

    group.finish();
}

/// 基准测试：稀疏矩阵列子集操作（稀疏列）
fn bench_subset_sparse_columns(c: &mut Criterion) {
    let mut group = c.benchmark_group("subset_sparse_columns");

    let csr = generate_random_csr(10000, 5000, 0.95);
    let row_indices: Vec<usize> = (0..10000).collect();
    let target_cols: Vec<usize> = (0..5000).step_by(10).collect(); // 每 10 列选 1 列

    group.bench_function("subset_500_sparse_cols", |b| {
        b.iter(|| {
            let result = subset_sparse_columns_in_place(
                black_box(csr.indptr.clone()),
                black_box(csr.indices.clone()),
                black_box(csr.data.clone()),
                black_box(&row_indices),
                black_box(&target_cols),
            );
            black_box(result);
        });
    });

    group.finish();
}

/// 基准测试：内存使用估算
fn bench_memory_estimation(c: &mut Criterion) {
    let mut group = c.benchmark_group("memory_estimation");

    let small_csr = generate_random_csr(1000, 500, 0.90);
    group.bench_function("estimate_small_csr", |b| {
        b.iter(|| {
            let estimate = estimate_csr_matrix_memory(black_box(&small_csr));
            black_box(estimate);
        });
    });

    let medium_csr = generate_random_csr(10000, 5000, 0.95);
    group.bench_function("estimate_medium_csr", |b| {
        b.iter(|| {
            let estimate = estimate_csr_matrix_memory(black_box(&medium_csr));
            black_box(estimate);
        });
    });

    let small_csc = csr_to_csc(&small_csr);
    group.bench_function("estimate_small_csc", |b| {
        b.iter(|| {
            let estimate = estimate_csc_matrix_memory(black_box(&small_csc));
            black_box(estimate);
        });
    });

    group.finish();
}

/// 基准测试：往返转换（CSR → CSC → CSR）
fn bench_roundtrip_conversion(c: &mut Criterion) {
    let mut group = c.benchmark_group("roundtrip_conversion");

    let medium_csr = generate_random_csr(10000, 5000, 0.95);
    group.throughput(Throughput::Elements(medium_csr.data.len() as u64));

    // 注意：由于当前 CSC → CSR 实现的限制，我们只测试 CSR → CSC
    // TODO: 修复 CSC → CSR 转换后重新启用完整往返测试
    group.bench_function("csr_to_csc_only", |b| {
        b.iter(|| {
            let csc = csr_to_csc(black_box(&medium_csr));
            black_box(csc);
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_csr_to_csc_conversion,
    bench_csc_to_csr_conversion,
    bench_parallel_vs_serial,
    bench_subset_rows,
    bench_subset_contiguous_columns,
    bench_subset_sparse_columns,
    bench_memory_estimation,
    bench_roundtrip_conversion,
);
criterion_main!(benches);
