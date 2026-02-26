//! 类型兼容性验证
//!
//! 验证数据类型在格式转换时的兼容性

use crate::error::{CrossCellError, Result};
use arrow::datatypes::DataType;

/// 验证 Arrow DataType 可以转换为 HDF5
pub fn validate_arrow_to_hdf5_type(dtype: &DataType) -> Result<()> {
    match dtype {
        // 支持的数值类型
        DataType::Int8
        | DataType::Int16
        | DataType::Int32
        | DataType::Int64
        | DataType::UInt8
        | DataType::UInt16
        | DataType::UInt32
        | DataType::UInt64
        | DataType::Float32
        | DataType::Float64 => Ok(()),

        // 支持的字符串类型
        DataType::Utf8 | DataType::LargeUtf8 => Ok(()),

        // 支持的布尔类型
        DataType::Boolean => Ok(()),

        // 支持的 Categorical 类型
        DataType::Dictionary(key_type, value_type) => {
            // 验证 key 是整数类型
            match **key_type {
                DataType::Int8
                | DataType::Int16
                | DataType::Int32
                | DataType::Int64
                | DataType::UInt8
                | DataType::UInt16
                | DataType::UInt32
                | DataType::UInt64 => {}
                _ => {
                    return Err(CrossCellError::type_mismatch(
                        "Dictionary key",
                        "integer type",
                        format!("{:?}", key_type),
                    ))
                }
            }

            // 验证 value 是字符串类型
            match **value_type {
                DataType::Utf8 | DataType::LargeUtf8 => Ok(()),
                _ => Err(CrossCellError::type_mismatch(
                    "Dictionary value",
                    "string type",
                    format!("{:?}", value_type),
                )),
            }
        }

        // 支持的日期时间类型
        DataType::Date32 | DataType::Date64 | DataType::Timestamp(_, _) => Ok(()),

        // 不支持的类型
        _ => Err(CrossCellError::incompatible_types(
            format!("Arrow {:?}", dtype),
            "HDF5".to_string(),
        )),
    }
}

/// 验证 Arrow DataType 可以转换为 R 类型
pub fn validate_arrow_to_r_type(dtype: &DataType) -> Result<()> {
    match dtype {
        // 支持的数值类型
        DataType::Int8 | DataType::Int16 | DataType::Int32 | DataType::Int64 => Ok(()),

        DataType::Float32 | DataType::Float64 => Ok(()),

        // 支持的字符串类型
        DataType::Utf8 | DataType::LargeUtf8 => Ok(()),

        // 支持的布尔类型
        DataType::Boolean => Ok(()),

        // 支持的 Categorical 类型（转换为 R Factor）
        DataType::Dictionary(key_type, value_type) => {
            // 验证 key 是整数类型
            match **key_type {
                DataType::Int8
                | DataType::Int16
                | DataType::Int32
                | DataType::Int64
                | DataType::UInt8
                | DataType::UInt16
                | DataType::UInt32
                | DataType::UInt64 => {}
                _ => {
                    return Err(CrossCellError::type_mismatch(
                        "Dictionary key",
                        "integer type",
                        format!("{:?}", key_type),
                    ))
                }
            }

            // 验证 value 是字符串类型
            match **value_type {
                DataType::Utf8 | DataType::LargeUtf8 => Ok(()),
                _ => Err(CrossCellError::type_mismatch(
                    "Dictionary value",
                    "string type",
                    format!("{:?}", value_type),
                )),
            }
        }

        // 支持的日期时间类型
        DataType::Date32 | DataType::Date64 | DataType::Timestamp(_, _) => Ok(()),

        // 不支持的类型
        _ => Err(CrossCellError::incompatible_types(
            format!("Arrow {:?}", dtype),
            "R".to_string(),
        )),
    }
}

/// 验证 HDF5 类型可以转换为 Arrow
pub fn validate_hdf5_to_arrow_compatible(hdf5_type: &str) -> Result<()> {
    match hdf5_type {
        // 支持的数值类型
        "i8" | "i16" | "i32" | "i64" | "u8" | "u16" | "u32" | "u64" | "f32" | "f64" => Ok(()),

        // 支持的字符串类型
        "string" | "utf8" => Ok(()),

        // 支持的布尔类型
        "bool" => Ok(()),

        // 不支持的类型
        _ => Err(CrossCellError::incompatible_types(
            format!("HDF5 {}", hdf5_type),
            "Arrow".to_string(),
        )),
    }
}

/// 验证 R 类型可以转换为 Arrow
pub fn validate_r_to_arrow_compatible(r_type: &str) -> Result<()> {
    match r_type {
        // 支持的数值类型
        "integer" | "numeric" | "double" => Ok(()),

        // 支持的字符串类型
        "character" => Ok(()),

        // 支持的布尔类型
        "logical" => Ok(()),

        // 支持的 Factor 类型
        "factor" => Ok(()),

        // 支持的日期时间类型
        "Date" | "POSIXct" | "POSIXlt" => Ok(()),

        // 不支持的类型
        _ => Err(CrossCellError::incompatible_types(
            format!("R {}", r_type),
            "Arrow".to_string(),
        )),
    }
}

/// 生成类型兼容性报告
pub fn type_compatibility_report(dtype: &DataType, target_format: &str) -> String {
    let mut report = format!(
        "Type Compatibility Report for {:?} → {}:\n",
        dtype, target_format
    );

    let result = match target_format {
        "hdf5" | "anndata" => validate_arrow_to_hdf5_type(dtype),
        "rds" | "seurat" => validate_arrow_to_r_type(dtype),
        _ => Err(CrossCellError::NotSupported(format!(
            "Unknown target format: {}",
            target_format
        ))),
    };

    match result {
        Ok(_) => {
            report.push_str("✓ Type is compatible\n");
        }
        Err(e) => {
            report.push_str(&format!("✗ Type is incompatible: {}\n", e));
        }
    }

    report
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow::datatypes::TimeUnit;

    #[test]
    fn test_validate_arrow_to_hdf5_numeric_types() {
        assert!(validate_arrow_to_hdf5_type(&DataType::Int32).is_ok());
        assert!(validate_arrow_to_hdf5_type(&DataType::Int64).is_ok());
        assert!(validate_arrow_to_hdf5_type(&DataType::Float32).is_ok());
        assert!(validate_arrow_to_hdf5_type(&DataType::Float64).is_ok());
    }

    #[test]
    fn test_validate_arrow_to_hdf5_string_types() {
        assert!(validate_arrow_to_hdf5_type(&DataType::Utf8).is_ok());
        assert!(validate_arrow_to_hdf5_type(&DataType::LargeUtf8).is_ok());
    }

    #[test]
    fn test_validate_arrow_to_hdf5_boolean() {
        assert!(validate_arrow_to_hdf5_type(&DataType::Boolean).is_ok());
    }

    #[test]
    fn test_validate_arrow_to_hdf5_dictionary() {
        let dict_type = DataType::Dictionary(Box::new(DataType::Int32), Box::new(DataType::Utf8));
        assert!(validate_arrow_to_hdf5_type(&dict_type).is_ok());
    }

    #[test]
    fn test_validate_arrow_to_hdf5_invalid_dictionary_key() {
        let dict_type = DataType::Dictionary(
            Box::new(DataType::Float32), // Invalid key type
            Box::new(DataType::Utf8),
        );
        assert!(validate_arrow_to_hdf5_type(&dict_type).is_err());
    }

    #[test]
    fn test_validate_arrow_to_hdf5_datetime() {
        assert!(validate_arrow_to_hdf5_type(&DataType::Date32).is_ok());
        assert!(validate_arrow_to_hdf5_type(&DataType::Date64).is_ok());
        assert!(
            validate_arrow_to_hdf5_type(&DataType::Timestamp(TimeUnit::Millisecond, None)).is_ok()
        );
    }

    #[test]
    fn test_validate_arrow_to_hdf5_unsupported() {
        use arrow::datatypes::Fields;
        use std::sync::Arc;

        // List type is not supported
        assert!(validate_arrow_to_hdf5_type(&DataType::List(Arc::new(
            arrow::datatypes::Field::new("item", DataType::Int32, true)
        )))
        .is_err());

        // Struct type is not supported
        let fields: Vec<Arc<arrow::datatypes::Field>> = vec![];
        assert!(validate_arrow_to_hdf5_type(&DataType::Struct(Fields::from(fields))).is_err());
    }

    #[test]
    fn test_validate_arrow_to_r_numeric_types() {
        assert!(validate_arrow_to_r_type(&DataType::Int32).is_ok());
        assert!(validate_arrow_to_r_type(&DataType::Int64).is_ok());
        assert!(validate_arrow_to_r_type(&DataType::Float64).is_ok());
    }

    #[test]
    fn test_validate_arrow_to_r_string() {
        assert!(validate_arrow_to_r_type(&DataType::Utf8).is_ok());
    }

    #[test]
    fn test_validate_arrow_to_r_boolean() {
        assert!(validate_arrow_to_r_type(&DataType::Boolean).is_ok());
    }

    #[test]
    fn test_validate_arrow_to_r_dictionary() {
        let dict_type = DataType::Dictionary(Box::new(DataType::Int32), Box::new(DataType::Utf8));
        assert!(validate_arrow_to_r_type(&dict_type).is_ok());
    }

    #[test]
    fn test_validate_hdf5_to_arrow() {
        assert!(validate_hdf5_to_arrow_compatible("i32").is_ok());
        assert!(validate_hdf5_to_arrow_compatible("f64").is_ok());
        assert!(validate_hdf5_to_arrow_compatible("string").is_ok());
        assert!(validate_hdf5_to_arrow_compatible("bool").is_ok());
    }

    #[test]
    fn test_validate_hdf5_to_arrow_unsupported() {
        assert!(validate_hdf5_to_arrow_compatible("complex").is_err());
        assert!(validate_hdf5_to_arrow_compatible("unknown").is_err());
    }

    #[test]
    fn test_validate_r_to_arrow() {
        assert!(validate_r_to_arrow_compatible("integer").is_ok());
        assert!(validate_r_to_arrow_compatible("numeric").is_ok());
        assert!(validate_r_to_arrow_compatible("character").is_ok());
        assert!(validate_r_to_arrow_compatible("logical").is_ok());
        assert!(validate_r_to_arrow_compatible("factor").is_ok());
        assert!(validate_r_to_arrow_compatible("Date").is_ok());
    }

    #[test]
    fn test_validate_r_to_arrow_unsupported() {
        assert!(validate_r_to_arrow_compatible("closure").is_err());
        assert!(validate_r_to_arrow_compatible("environment").is_err());
        assert!(validate_r_to_arrow_compatible("unknown").is_err());
    }

    #[test]
    fn test_type_compatibility_report() {
        let report = type_compatibility_report(&DataType::Int32, "hdf5");
        assert!(report.contains("compatible"));

        let report = type_compatibility_report(&DataType::Utf8, "seurat");
        assert!(report.contains("compatible"));
    }
}
