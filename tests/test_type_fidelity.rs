//! Task 3: 数据类型保真度测试
//!
//! 验证 crosscell_original_dtype 属性的写入和读取，
//! 以及启发式类型恢复功能。

use arrow::array::{
    Array, ArrayRef, BooleanArray, DictionaryArray, Float64Array, Int32Array, Int64Array,
    StringArray,
};
use arrow::datatypes::{DataType, Int32Type};
use crosscell::anndata::{read_h5ad, write_h5ad};
use crosscell::ir::{
    DataFrame, DatasetMetadata, ExpressionMatrix, SingleCellData, SparseMatrixCSR,
};
use std::sync::Arc;
use tempfile::NamedTempFile;

/// 辅助函数：创建最小的 SingleCellData 用于测试
fn make_test_data(n_rows: usize, columns: Vec<String>, data: Vec<ArrayRef>) -> SingleCellData {
    let csr = SparseMatrixCSR::new(vec![1.0], vec![0], vec![0, 0, 1], 2, 1).unwrap();
    // 需要 n_rows 行的 cell_metadata
    let expression = if n_rows == 2 {
        ExpressionMatrix::SparseCSR(csr)
    } else {
        // 创建匹配行数的矩阵
        let mut indptr = vec![0usize; n_rows + 1];
        indptr[n_rows] = 1;
        ExpressionMatrix::SparseCSR(SparseMatrixCSR {
            data: vec![1.0],
            indices: vec![0],
            indptr,
            n_rows,
            n_cols: 1,
        })
    };

    let cell_metadata = DataFrame::new(columns, data, n_rows).unwrap();
    let gene_metadata = DataFrame::empty(1);
    let metadata = DatasetMetadata::new(n_rows, 1, "test".to_string());

    SingleCellData::new(expression, cell_metadata, gene_metadata, metadata).unwrap()
}

/// 辅助函数：写入 H5AD 然后读回来
fn roundtrip_h5ad(data: &SingleCellData) -> SingleCellData {
    let temp_file = NamedTempFile::new().unwrap();
    let temp_path = temp_file.path().with_extension("h5ad");

    write_h5ad(data, &temp_path).unwrap();
    let result = read_h5ad(&temp_path).unwrap();

    std::fs::remove_file(&temp_path).ok();
    result
}

// ============================================================================
// Task 3.4: 单元测试
// ============================================================================

#[test]
fn test_int32_roundtrip_with_attribute() {
    // Int32 列写入后应该带有 crosscell_original_dtype="int32" 属性
    // 读回来应该仍然是 Int32（而不是退化为 Float64）
    let int_values: Vec<i32> = vec![1, 2, 3, 42, 100];
    let int_array: ArrayRef = Arc::new(Int32Array::from(int_values.clone()));

    let data = make_test_data(5, vec!["nCount_RNA".to_string()], vec![int_array]);

    let result = roundtrip_h5ad(&data);

    let col = result.cell_metadata.column("nCount_RNA").unwrap();
    assert_eq!(
        *col.data_type(),
        DataType::Int32,
        "Int32 column should survive roundtrip with crosscell_original_dtype attribute"
    );

    let restored = col.as_any().downcast_ref::<Int32Array>().unwrap();
    for (i, &expected) in int_values.iter().enumerate() {
        assert_eq!(restored.value(i), expected);
    }
    println!("✓ Int32 roundtrip with attribute: passed");
}

#[test]
fn test_int64_roundtrip_with_attribute() {
    let int_values: Vec<i64> = vec![100000, 200000, 300000];
    let int_array: ArrayRef = Arc::new(Int64Array::from(int_values.clone()));

    let data = make_test_data(3, vec!["big_count".to_string()], vec![int_array]);

    let result = roundtrip_h5ad(&data);

    let col = result.cell_metadata.column("big_count").unwrap();
    assert_eq!(*col.data_type(), DataType::Int64);

    let restored = col.as_any().downcast_ref::<Int64Array>().unwrap();
    for (i, &expected) in int_values.iter().enumerate() {
        assert_eq!(restored.value(i), expected);
    }
    println!("✓ Int64 roundtrip with attribute: passed");
}

#[test]
fn test_float64_roundtrip_with_attribute() {
    let float_values: Vec<f64> = vec![1.5, 2.7, 3.14];
    let float_array: ArrayRef = Arc::new(Float64Array::from(float_values.clone()));

    let data = make_test_data(3, vec!["score".to_string()], vec![float_array]);

    let result = roundtrip_h5ad(&data);

    let col = result.cell_metadata.column("score").unwrap();
    assert_eq!(*col.data_type(), DataType::Float64);

    let restored = col.as_any().downcast_ref::<Float64Array>().unwrap();
    for (i, &expected) in float_values.iter().enumerate() {
        assert!((restored.value(i) - expected).abs() < 1e-10);
    }
    println!("✓ Float64 roundtrip with attribute: passed");
}

#[test]
fn test_boolean_roundtrip_with_attribute() {
    let bool_values = vec![true, false, true];
    let bool_array: ArrayRef = Arc::new(BooleanArray::from(bool_values.clone()));

    let data = make_test_data(3, vec!["is_primary".to_string()], vec![bool_array]);

    let result = roundtrip_h5ad(&data);

    let col = result.cell_metadata.column("is_primary").unwrap();
    assert_eq!(*col.data_type(), DataType::Boolean);

    let restored = col.as_any().downcast_ref::<BooleanArray>().unwrap();
    for (i, &expected) in bool_values.iter().enumerate() {
        assert_eq!(restored.value(i), expected);
    }
    println!("✓ Boolean roundtrip with attribute: passed");
}

#[test]
fn test_string_roundtrip_with_attribute() {
    let str_values = vec!["hello", "world", "test"];
    let str_array: ArrayRef = Arc::new(StringArray::from(str_values.clone()));

    let data = make_test_data(3, vec!["label".to_string()], vec![str_array]);

    let result = roundtrip_h5ad(&data);

    let col = result.cell_metadata.column("label").unwrap();
    assert_eq!(*col.data_type(), DataType::Utf8);

    let restored = col.as_any().downcast_ref::<StringArray>().unwrap();
    for (i, &expected) in str_values.iter().enumerate() {
        assert_eq!(restored.value(i), expected);
    }
    println!("✓ String roundtrip with attribute: passed");
}

#[test]
fn test_categorical_roundtrip_with_attribute() {
    // Dictionary<Int32, Utf8> 列
    let keys = Int32Array::from(vec![0, 1, 0, 1, 2]);
    let values = Arc::new(StringArray::from(vec!["A", "B", "C"]));
    let dict_array: ArrayRef =
        Arc::new(DictionaryArray::<Int32Type>::try_new(keys, values).unwrap());

    let data = make_test_data(5, vec!["cell_type".to_string()], vec![dict_array]);

    let result = roundtrip_h5ad(&data);

    let col = result.cell_metadata.column("cell_type").unwrap();
    assert!(
        matches!(col.data_type(), DataType::Dictionary(_, _)),
        "Categorical column should survive roundtrip"
    );
    println!("✓ Categorical roundtrip with attribute: passed");
}

#[test]
fn test_mixed_types_roundtrip() {
    // 多种类型混合的 DataFrame
    let int32_col: ArrayRef = Arc::new(Int32Array::from(vec![1, 2, 3]));
    let float64_col: ArrayRef = Arc::new(Float64Array::from(vec![1.5, 2.5, 3.5]));
    let str_col: ArrayRef = Arc::new(StringArray::from(vec!["a", "b", "c"]));

    let data = make_test_data(
        3,
        vec!["count".to_string(), "score".to_string(), "name".to_string()],
        vec![int32_col, float64_col, str_col],
    );

    let result = roundtrip_h5ad(&data);

    // Int32 应该保持
    let count_col = result.cell_metadata.column("count").unwrap();
    assert_eq!(
        *count_col.data_type(),
        DataType::Int32,
        "Int32 should be preserved"
    );

    // Float64 应该保持
    let score_col = result.cell_metadata.column("score").unwrap();
    assert_eq!(
        *score_col.data_type(),
        DataType::Float64,
        "Float64 should be preserved"
    );

    // String 应该保持
    let name_col = result.cell_metadata.column("name").unwrap();
    assert_eq!(
        *name_col.data_type(),
        DataType::Utf8,
        "Utf8 should be preserved"
    );

    println!("✓ Mixed types roundtrip: passed");
}

// ============================================================================
// Task 3.3: 属性测试 - Property 4: DataFrame 类型往返保真
// ============================================================================

#[cfg(test)]
mod prop_tests {
    use super::*;
    use proptest::prelude::*;

    // Feature: reviewer-enhancements, Property 4: DataFrame 类型往返保真
    proptest! {
        #[test]
        fn prop_int32_roundtrip_preserves_type(
            values in proptest::collection::vec(-1000i32..1000, 2..20),
        ) {
            let n = values.len();
            let int_array: ArrayRef = Arc::new(Int32Array::from(values.clone()));

            let data = make_test_data(
                n,
                vec!["test_col".to_string()],
                vec![int_array],
            );

            let result = roundtrip_h5ad(&data);
            let col = result.cell_metadata.column("test_col").unwrap();

            prop_assert_eq!(
                col.data_type().clone(),
                DataType::Int32,
                "Int32 type should be preserved after roundtrip"
            );

            let restored = col.as_any().downcast_ref::<Int32Array>().unwrap();
            for (i, &expected) in values.iter().enumerate() {
                prop_assert_eq!(restored.value(i), expected);
            }
        }

        #[test]
        fn prop_float64_with_decimals_stays_float64(
            values in proptest::collection::vec(-100.0f64..100.0, 2..20),
        ) {
            // 确保至少有一个非整数值，这样启发式不会误转
            let mut values = values;
            values[0] = 1.5; // 保证有小数

            let n = values.len();
            let float_array: ArrayRef = Arc::new(Float64Array::from(values.clone()));

            let data = make_test_data(
                n,
                vec!["test_col".to_string()],
                vec![float_array],
            );

            let result = roundtrip_h5ad(&data);
            let col = result.cell_metadata.column("test_col").unwrap();

            prop_assert_eq!(
                col.data_type().clone(),
                DataType::Float64,
                "Float64 with decimal values should stay Float64"
            );
        }
    }
}
