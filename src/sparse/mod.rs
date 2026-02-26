//! 稀疏矩阵操作

pub mod convert;
pub mod memory;
pub mod subset;
pub mod validate;

pub use convert::{csc_to_csr, csr_to_csc};
pub use memory::{
    estimate_csc_matrix_memory, estimate_csc_memory_usage, estimate_csr_matrix_memory,
    estimate_csr_memory_usage, format_memory_size, MemoryEstimate,
};
pub use subset::is_contiguous;
pub use validate::validate_sparse_matrix;
