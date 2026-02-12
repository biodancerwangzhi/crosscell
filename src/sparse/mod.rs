//! 稀疏矩阵操作

pub mod convert;
pub mod validate;
pub mod subset;
pub mod memory;

pub use convert::{csc_to_csr, csr_to_csc};
pub use validate::validate_sparse_matrix;
pub use subset::{is_contiguous};
pub use memory::{
    estimate_csr_memory_usage, estimate_csc_memory_usage,
    estimate_csr_matrix_memory, estimate_csc_matrix_memory,
    format_memory_size, MemoryEstimate,
};
