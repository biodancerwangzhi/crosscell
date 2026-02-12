//! CrossCell 统一错误处理
//!
//! 定义所有 CrossCell 操作的统一错误类型，提供清晰的错误消息和上下文信息。

use thiserror::Error;

/// CrossCell 统一错误类型
///
/// 整合所有子模块的错误类型，提供统一的错误处理接口。
#[derive(Error, Debug)]
pub enum CrossCellError {
    // ========== I/O 错误 ==========
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("File not found: {path}")]
    FileNotFound { path: String },

    #[error("Failed to read file '{path}': {source}")]
    FileReadError {
        path: String,
        source: std::io::Error,
    },

    #[error("Failed to write file '{path}': {source}")]
    FileWriteError {
        path: String,
        source: std::io::Error,
    },

    // ========== 格式错误 ==========
    #[error("Invalid file format: {0}")]
    InvalidFormat(String),

    #[error("Unsupported file format: expected {expected}, got {actual}")]
    UnsupportedFormat { expected: String, actual: String },

    #[error("Invalid magic bytes: expected {expected:?}, got {actual:?}")]
    InvalidMagic { expected: Vec<u8>, actual: Vec<u8> },

    // ========== 数据验证错误 ==========
    #[error("Missing required field: {field}")]
    MissingField { field: String },

    #[error("Missing required fields: {fields:?}")]
    MissingFields { fields: Vec<String> },

    #[error("Dimension mismatch: {context}")]
    DimensionMismatch { context: String },

    #[error("Dimension mismatch: expected {expected}, got {actual} for {field}")]
    DimensionMismatchDetailed {
        field: String,
        expected: String,
        actual: String,
    },

    #[error("Invalid dimensions: {0}")]
    InvalidDimensions(String),

    // ========== 类型错误 ==========
    #[error("Unsupported data type: {type_name}")]
    UnsupportedType { type_name: String },

    #[error("Type mismatch: expected {expected}, got {actual} for field '{field}'")]
    TypeMismatch {
        field: String,
        expected: String,
        actual: String,
    },

    #[error("Incompatible types: cannot convert {from} to {to}")]
    IncompatibleTypes { from: String, to: String },

    // ========== 稀疏矩阵错误 ==========
    #[error("Sparse matrix validation error: {0}")]
    SparseValidation(#[from] crate::sparse::validate::ValidationError),

    #[error("Invalid sparse matrix structure: {0}")]
    InvalidSparseMatrix(String),

    // ========== AnnData 特定错误 ==========
    #[error("AnnData error: {0}")]
    AnnData(#[from] crate::anndata::AnnDataError),

    #[error("HDF5 error: {0}")]
    Hdf5(#[from] hdf5::Error),

    // ========== Seurat 特定错误 ==========
    #[error("Seurat error: {0}")]
    Seurat(#[from] crate::seurat::error::SeuratError),

    // ========== RDS 特定错误 ==========
    #[error("RDS error: {0}")]
    Rds(#[from] crate::rds::RdsError),

    // ========== 转换错误 ==========
    #[error("Conversion error: {0}")]
    Conversion(String),

    #[error("Failed to convert {from} to {to}: {reason}")]
    ConversionFailed {
        from: String,
        to: String,
        reason: String,
    },

    // ========== 验证错误 ==========
    #[error("Validation failed: {0}")]
    ValidationFailed(String),

    #[error("Roundtrip validation failed: {details}")]
    RoundtripValidationFailed { details: String },

    // ========== 解析错误 ==========
    #[error("Parse error: {0}")]
    Parse(String),

    #[error("Failed to parse {field}: {reason}")]
    ParseFailed { field: String, reason: String },

    // ========== 内存错误 ==========
    #[error("Out of memory: {context}")]
    OutOfMemory { context: String },

    #[error("Memory allocation failed: requested {requested} bytes")]
    AllocationFailed { requested: usize },

    // ========== 其他错误 ==========
    #[error("Internal error: {0}")]
    Internal(String),

    #[error("Not implemented: {0}")]
    NotImplemented(String),

    #[error("Operation not supported: {0}")]
    NotSupported(String),

    #[error("Not found: {0}")]
    NotFound(String),

    #[error("Unsupported operation: {0}")]
    UnsupportedOperation(String),
}

/// CrossCell Result 类型别名
pub type Result<T> = std::result::Result<T, CrossCellError>;

impl CrossCellError {
    /// 创建文件未找到错误
    pub fn file_not_found(path: impl Into<String>) -> Self {
        Self::FileNotFound { path: path.into() }
    }

    /// 创建文件读取错误
    pub fn file_read_error(path: impl Into<String>, source: std::io::Error) -> Self {
        Self::FileReadError {
            path: path.into(),
            source,
        }
    }

    /// 创建文件写入错误
    pub fn file_write_error(path: impl Into<String>, source: std::io::Error) -> Self {
        Self::FileWriteError {
            path: path.into(),
            source,
        }
    }

    /// 创建缺失字段错误
    pub fn missing_field(field: impl Into<String>) -> Self {
        Self::MissingField {
            field: field.into(),
        }
    }

    /// 创建缺失多个字段错误
    pub fn missing_fields(fields: Vec<String>) -> Self {
        Self::MissingFields { fields }
    }

    /// 创建维度不匹配错误
    pub fn dimension_mismatch(context: impl Into<String>) -> Self {
        Self::DimensionMismatch {
            context: context.into(),
        }
    }

    /// 创建详细的维度不匹配错误
    pub fn dimension_mismatch_detailed(
        field: impl Into<String>,
        expected: impl Into<String>,
        actual: impl Into<String>,
    ) -> Self {
        Self::DimensionMismatchDetailed {
            field: field.into(),
            expected: expected.into(),
            actual: actual.into(),
        }
    }

    /// 创建类型不匹配错误
    pub fn type_mismatch(
        field: impl Into<String>,
        expected: impl Into<String>,
        actual: impl Into<String>,
    ) -> Self {
        Self::TypeMismatch {
            field: field.into(),
            expected: expected.into(),
            actual: actual.into(),
        }
    }

    /// 创建不兼容类型错误
    pub fn incompatible_types(from: impl Into<String>, to: impl Into<String>) -> Self {
        Self::IncompatibleTypes {
            from: from.into(),
            to: to.into(),
        }
    }

    /// 创建转换失败错误
    pub fn conversion_failed(
        from: impl Into<String>,
        to: impl Into<String>,
        reason: impl Into<String>,
    ) -> Self {
        Self::ConversionFailed {
            from: from.into(),
            to: to.into(),
            reason: reason.into(),
        }
    }

    /// 创建解析失败错误
    pub fn parse_failed(field: impl Into<String>, reason: impl Into<String>) -> Self {
        Self::ParseFailed {
            field: field.into(),
            reason: reason.into(),
        }
    }

    /// 检查错误是否是 I/O 错误
    pub fn is_io_error(&self) -> bool {
        matches!(
            self,
            Self::Io(_) | Self::FileReadError { .. } | Self::FileWriteError { .. }
        )
    }

    /// 检查错误是否是验证错误
    pub fn is_validation_error(&self) -> bool {
        matches!(
            self,
            Self::ValidationFailed(_)
                | Self::RoundtripValidationFailed { .. }
                | Self::SparseValidation(_)
                | Self::DimensionMismatch { .. }
                | Self::DimensionMismatchDetailed { .. }
        )
    }

    /// 检查错误是否是类型错误
    pub fn is_type_error(&self) -> bool {
        matches!(
            self,
            Self::UnsupportedType { .. }
                | Self::TypeMismatch { .. }
                | Self::IncompatibleTypes { .. }
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_creation() {
        let err = CrossCellError::file_not_found("test.h5ad");
        assert!(matches!(err, CrossCellError::FileNotFound { .. }));

        let err = CrossCellError::missing_field("X");
        assert!(matches!(err, CrossCellError::MissingField { .. }));

        let err = CrossCellError::dimension_mismatch("matrix shape mismatch");
        assert!(matches!(err, CrossCellError::DimensionMismatch { .. }));
    }

    #[test]
    fn test_error_display() {
        let err = CrossCellError::file_not_found("test.h5ad");
        assert_eq!(err.to_string(), "File not found: test.h5ad");

        let err = CrossCellError::missing_field("X");
        assert_eq!(err.to_string(), "Missing required field: X");

        let err = CrossCellError::dimension_mismatch_detailed(
            "expression_matrix",
            "100x50",
            "100x60",
        );
        assert_eq!(
            err.to_string(),
            "Dimension mismatch: expected 100x50, got 100x60 for expression_matrix"
        );
    }

    #[test]
    fn test_error_classification() {
        let io_err = CrossCellError::Io(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            "file not found",
        ));
        assert!(io_err.is_io_error());
        assert!(!io_err.is_validation_error());

        let val_err = CrossCellError::dimension_mismatch("test");
        assert!(val_err.is_validation_error());
        assert!(!val_err.is_io_error());

        let type_err = CrossCellError::type_mismatch("field", "int", "string");
        assert!(type_err.is_type_error());
        assert!(!type_err.is_validation_error());
    }
}
