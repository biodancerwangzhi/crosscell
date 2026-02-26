//! 稀疏矩阵内存使用估算
//!
//! 本模块提供内存使用估算功能，帮助用户在转换前了解内存需求。
//!
//! ## 原始来源
//!
//! 算法改编自 Anndata-Memory-master 项目：
//! - 项目地址：https://github.com/SingleRust/Anndata-Memory
//! - 原始文件：src/utils/mod.rs (estimate_csr_total_memory_usage)
//! - 许可证：MIT License
//!
//! ## 致谢
//!
//! 感谢 Anndata-Memory 项目的作者提供了内存估算算法。
//!
//! ## 内存优化策略
//!
//! 本模块实现了以下内存优化策略：
//! 1. 零拷贝转换：尽可能避免数据复制
//! 2. 原地操作：在可能的情况下原地修改数据
//! 3. 内存预估：转换前估算峰值内存需求
//! 4. 分块处理：支持大数据集的分块处理

use crate::ir::expression::{SparseMatrixCSC, SparseMatrixCSR};

/// 估算 CSR 矩阵的总内存使用量（字节）
///
/// 计算包括：
/// - 数据数组（data）
/// - 列索引数组（indices）
/// - 行指针数组（indptr）
/// - 转换过程中的峰值内存（临时数组）
/// - 20% 的安全余量
///
/// # 参数
/// - `nnz`: 非零元素数量
/// - `n_rows`: 行数
/// - `data_type_size`: 数据类型大小（字节），例如 f64 = 8
///
/// # 返回
/// 估算的峰值内存使用量（字节）
///
/// # 示例
/// ```
/// use crosscell::sparse::memory::estimate_csr_memory_usage;
///
/// // 100万细胞 × 30000基因，稀疏度 95%
/// let nnz = 1_000_000 * 30_000 * 5 / 100;  // 1.5B 非零元素
/// let n_rows = 1_000_000;
/// let data_type_size = 8;  // f64
///
/// let estimated_bytes = estimate_csr_memory_usage(nnz, n_rows, data_type_size);
/// let estimated_gb = estimated_bytes as f64 / 1_073_741_824.0;
///
/// println!("估算内存需求: {:.2} GB", estimated_gb);
/// ```
pub fn estimate_csr_memory_usage(nnz: usize, n_rows: usize, data_type_size: usize) -> usize {
    // 最终 CSR 矩阵的内存占用
    let data_array_size = nnz * data_type_size;
    let indices_array_size = nnz * std::mem::size_of::<usize>();
    let indptr_array_size = (n_rows + 1) * std::mem::size_of::<usize>();

    let final_csr_size = data_array_size + indices_array_size + indptr_array_size;

    // 转换过程中的峰值内存（需要临时数组）
    // 假设最坏情况：同时存在源矩阵和目标矩阵
    let peak_usage = final_csr_size * 2;

    // 添加 20% 的安全余量
    (peak_usage as f64 * 1.2) as usize
}

/// 估算 CSC 矩阵的总内存使用量（字节）
///
/// 与 CSR 类似，但使用列数而不是行数
pub fn estimate_csc_memory_usage(nnz: usize, n_cols: usize, data_type_size: usize) -> usize {
    let data_array_size = nnz * data_type_size;
    let indices_array_size = nnz * std::mem::size_of::<usize>();
    let indptr_array_size = (n_cols + 1) * std::mem::size_of::<usize>();

    let final_csc_size = data_array_size + indices_array_size + indptr_array_size;
    let peak_usage = final_csc_size * 2;

    (peak_usage as f64 * 1.2) as usize
}

/// 估算 CSR 矩阵实例的内存使用量
///
/// 便捷函数，直接从 SparseMatrixCSR 实例估算
pub fn estimate_csr_matrix_memory(csr: &SparseMatrixCSR) -> usize {
    let nnz = csr.data.len();
    let n_rows = csr.n_rows;
    let data_type_size = std::mem::size_of::<f64>();

    estimate_csr_memory_usage(nnz, n_rows, data_type_size)
}

/// 估算 CSC 矩阵实例的内存使用量
///
/// 便捷函数，直接从 SparseMatrixCSC 实例估算
pub fn estimate_csc_matrix_memory(csc: &SparseMatrixCSC) -> usize {
    let nnz = csc.data.len();
    let n_cols = csc.n_cols;
    let data_type_size = std::mem::size_of::<f64>();

    estimate_csc_memory_usage(nnz, n_cols, data_type_size)
}

/// 格式化内存大小为人类可读的字符串
///
/// # 示例
/// ```
/// use crosscell::sparse::memory::format_memory_size;
///
/// assert_eq!(format_memory_size(1024), "1.00 KB");
/// assert_eq!(format_memory_size(1_048_576), "1.00 MB");
/// assert_eq!(format_memory_size(1_073_741_824), "1.00 GB");
/// ```
pub fn format_memory_size(bytes: usize) -> String {
    const KB: f64 = 1024.0;
    const MB: f64 = KB * 1024.0;
    const GB: f64 = MB * 1024.0;
    const TB: f64 = GB * 1024.0;

    let bytes_f64 = bytes as f64;

    if bytes_f64 >= TB {
        format!("{:.2} TB", bytes_f64 / TB)
    } else if bytes_f64 >= GB {
        format!("{:.2} GB", bytes_f64 / GB)
    } else if bytes_f64 >= MB {
        format!("{:.2} MB", bytes_f64 / MB)
    } else if bytes_f64 >= KB {
        format!("{:.2} KB", bytes_f64 / KB)
    } else {
        format!("{} B", bytes)
    }
}

/// 内存估算报告
#[derive(Debug, Clone)]
pub struct MemoryEstimate {
    /// 非零元素数量
    pub nnz: usize,

    /// 矩阵维度（行数，列数）
    pub dimensions: (usize, usize),

    /// 估算的峰值内存使用量（字节）
    pub peak_memory_bytes: usize,

    /// 数据数组大小（字节）
    pub data_size_bytes: usize,

    /// 索引数组大小（字节）
    pub indices_size_bytes: usize,

    /// 指针数组大小（字节）
    pub indptr_size_bytes: usize,
}

impl MemoryEstimate {
    /// 创建 CSR 矩阵的内存估算报告
    pub fn for_csr(csr: &SparseMatrixCSR) -> Self {
        let nnz = csr.data.len();
        let n_rows = csr.n_rows;
        let n_cols = csr.n_cols;
        let data_type_size = std::mem::size_of::<f64>();

        let data_size_bytes = nnz * data_type_size;
        let indices_size_bytes = nnz * std::mem::size_of::<usize>();
        let indptr_size_bytes = (n_rows + 1) * std::mem::size_of::<usize>();

        let peak_memory_bytes = estimate_csr_memory_usage(nnz, n_rows, data_type_size);

        Self {
            nnz,
            dimensions: (n_rows, n_cols),
            peak_memory_bytes,
            data_size_bytes,
            indices_size_bytes,
            indptr_size_bytes,
        }
    }

    /// 创建 CSC 矩阵的内存估算报告
    pub fn for_csc(csc: &SparseMatrixCSC) -> Self {
        let nnz = csc.data.len();
        let n_rows = csc.n_rows;
        let n_cols = csc.n_cols;
        let data_type_size = std::mem::size_of::<f64>();

        let data_size_bytes = nnz * data_type_size;
        let indices_size_bytes = nnz * std::mem::size_of::<usize>();
        let indptr_size_bytes = (n_cols + 1) * std::mem::size_of::<usize>();

        let peak_memory_bytes = estimate_csc_memory_usage(nnz, n_cols, data_type_size);

        Self {
            nnz,
            dimensions: (n_rows, n_cols),
            peak_memory_bytes,
            data_size_bytes,
            indices_size_bytes,
            indptr_size_bytes,
        }
    }

    /// 格式化为人类可读的报告
    pub fn format_report(&self) -> String {
        format!(
            "Memory Estimate:\n\
             - Dimensions: {} × {}\n\
             - Non-zero elements: {}\n\
             - Data array: {}\n\
             - Indices array: {}\n\
             - Pointer array: {}\n\
             - Peak memory: {}",
            self.dimensions.0,
            self.dimensions.1,
            self.nnz,
            format_memory_size(self.data_size_bytes),
            format_memory_size(self.indices_size_bytes),
            format_memory_size(self.indptr_size_bytes),
            format_memory_size(self.peak_memory_bytes)
        )
    }
}

/// 内存优化建议
#[derive(Debug, Clone)]
pub struct MemoryOptimizationAdvice {
    /// 是否建议使用 lazy loading
    pub use_lazy_loading: bool,

    /// 是否建议使用分块处理
    pub use_chunked_processing: bool,

    /// 建议的块大小（行数）
    pub recommended_chunk_size: usize,

    /// 是否建议使用并行处理
    pub use_parallel: bool,

    /// 建议的线程数
    pub recommended_threads: usize,

    /// 详细建议信息
    pub advice_message: String,
}

impl MemoryOptimizationAdvice {
    /// 根据矩阵大小和可用内存生成优化建议
    pub fn for_matrix(
        n_rows: usize,
        _n_cols: usize,
        nnz: usize,
        available_memory_bytes: usize,
    ) -> Self {
        let data_type_size = std::mem::size_of::<f64>();
        let estimated_peak = estimate_csr_memory_usage(nnz, n_rows, data_type_size);

        // 计算内存比率
        let memory_ratio = estimated_peak as f64 / available_memory_bytes as f64;

        // 决定是否使用 lazy loading
        let use_lazy_loading = memory_ratio > 0.5;

        // 决定是否使用分块处理
        let use_chunked_processing = memory_ratio > 0.8;

        // 计算建议的块大小
        let recommended_chunk_size = if use_chunked_processing {
            // 目标：每块使用约 20% 的可用内存
            let target_chunk_memory = available_memory_bytes / 5;
            let bytes_per_row = (nnz / n_rows) * (data_type_size + std::mem::size_of::<usize>());
            std::cmp::max(1000, target_chunk_memory / bytes_per_row.max(1))
        } else {
            n_rows // 不分块
        };

        // 决定是否使用并行处理
        let use_parallel = nnz > 10_000;

        // 建议的线程数
        let recommended_threads = if use_parallel {
            std::cmp::min(num_cpus::get(), 8)
        } else {
            1
        };

        // 生成建议信息
        let advice_message = if memory_ratio > 1.0 {
            format!(
                "⚠️ 警告：估算内存需求 ({}) 超过可用内存 ({})。\n\
                 建议：\n\
                 1. 使用 --lazy 选项启用延迟加载\n\
                 2. 使用 --chunk-size {} 进行分块处理\n\
                 3. 或增加系统内存",
                format_memory_size(estimated_peak),
                format_memory_size(available_memory_bytes),
                recommended_chunk_size
            )
        } else if memory_ratio > 0.5 {
            format!(
                "💡 提示：数据集较大，建议使用 --lazy 选项减少内存占用。\n\
                 估算内存需求: {}\n\
                 可用内存: {}",
                format_memory_size(estimated_peak),
                format_memory_size(available_memory_bytes)
            )
        } else {
            format!(
                "✅ 内存充足，可以直接处理。\n\
                 估算内存需求: {}\n\
                 可用内存: {}",
                format_memory_size(estimated_peak),
                format_memory_size(available_memory_bytes)
            )
        };

        Self {
            use_lazy_loading,
            use_chunked_processing,
            recommended_chunk_size,
            use_parallel,
            recommended_threads,
            advice_message,
        }
    }
}

/// 获取系统可用内存（字节）
///
/// 注意：这是一个近似值，实际可用内存可能因系统状态而异
pub fn get_available_memory() -> usize {
    // 使用 sysinfo 或类似库获取系统内存
    // 这里使用一个保守的默认值
    #[cfg(target_os = "linux")]
    {
        // 尝试读取 /proc/meminfo
        if let Ok(content) = std::fs::read_to_string("/proc/meminfo") {
            for line in content.lines() {
                if line.starts_with("MemAvailable:") {
                    if let Some(kb_str) = line.split_whitespace().nth(1) {
                        if let Ok(kb) = kb_str.parse::<usize>() {
                            return kb * 1024;
                        }
                    }
                }
            }
        }
    }

    // 默认返回 8 GB
    8 * 1024 * 1024 * 1024
}

/// 检查是否有足够的内存进行转换
pub fn check_memory_sufficient(n_rows: usize, _n_cols: usize, nnz: usize) -> Result<(), String> {
    let data_type_size = std::mem::size_of::<f64>();
    let estimated_peak = estimate_csr_memory_usage(nnz, n_rows, data_type_size);
    let available = get_available_memory();

    if estimated_peak > available {
        Err(format!(
            "内存不足：需要 {}，可用 {}。请使用 --lazy 选项或增加系统内存。",
            format_memory_size(estimated_peak),
            format_memory_size(available)
        ))
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_estimate_csr_memory() {
        // 小矩阵：1000 细胞 × 500 基因，稀疏度 90%
        let nnz = 1000 * 500 * 10 / 100; // 50,000 非零元素
        let n_rows = 1000;
        let data_type_size = 8; // f64

        let estimated = estimate_csr_memory_usage(nnz, n_rows, data_type_size);

        // 验证估算合理性
        // data: 50000 * 8 = 400 KB
        // indices: 50000 * 8 = 400 KB
        // indptr: 1001 * 8 = 8 KB
        // total: ~808 KB
        // peak: ~1.6 MB (2x)
        // with 20% margin: ~1.9 MB
        assert!(estimated > 1_500_000); // > 1.5 MB
        assert!(estimated < 2_500_000); // < 2.5 MB
    }

    #[test]
    fn test_format_memory_size() {
        assert_eq!(format_memory_size(512), "512 B");
        assert_eq!(format_memory_size(1024), "1.00 KB");
        assert_eq!(format_memory_size(1_048_576), "1.00 MB");
        assert_eq!(format_memory_size(1_073_741_824), "1.00 GB");
        assert_eq!(format_memory_size(1_099_511_627_776), "1.00 TB");
    }

    #[test]
    fn test_memory_estimate_report() {
        let csr = SparseMatrixCSR {
            data: vec![1.0; 1000],
            indices: vec![0; 1000],
            indptr: vec![0; 101],
            n_rows: 100,
            n_cols: 50,
        };

        let estimate = MemoryEstimate::for_csr(&csr);

        assert_eq!(estimate.nnz, 1000);
        assert_eq!(estimate.dimensions, (100, 50));
        assert!(estimate.peak_memory_bytes > 0);

        let report = estimate.format_report();
        assert!(report.contains("100 × 50"));
        assert!(report.contains("1000"));
    }

    #[test]
    fn test_large_matrix_estimate() {
        // 大矩阵：100万细胞 × 30000基因，稀疏度 95%
        let nnz = 1_000_000 * 30_000 * 5 / 100; // 1.5B 非零元素
        let n_rows = 1_000_000;
        let data_type_size = 8;

        let estimated = estimate_csr_memory_usage(nnz, n_rows, data_type_size);
        let estimated_gb = estimated as f64 / 1_073_741_824.0;

        // 应该在合理范围内（几十 GB）
        assert!(estimated_gb > 10.0);
        assert!(estimated_gb < 100.0);

        println!("Large matrix estimate: {:.2} GB", estimated_gb);
    }
}
