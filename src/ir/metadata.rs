//! 元数据 IR 类型
//!
//! 使用 Apache Arrow 存储细胞和基因的元数据，支持零拷贝跨语言传输。

use arrow::array::ArrayRef;
use std::fmt;

/// DataFrame 使用 Apache Arrow 存储元数据
///
/// 对应：
/// - AnnData.obs / AnnData.var
/// - Seurat@meta.data / Seurat@assays$RNA@meta.features
#[derive(Debug, Clone)]
pub struct DataFrame {
    /// 列名
    pub columns: Vec<String>,
    /// 列数据（Arrow 数组）
    pub data: Vec<ArrayRef>,
    /// 行数
    pub n_rows: usize,
}

impl DataFrame {
    /// 创建新的 DataFrame
    pub fn new(columns: Vec<String>, data: Vec<ArrayRef>, n_rows: usize) -> Result<Self, String> {
        let df = Self {
            columns,
            data,
            n_rows,
        };
        df.validate()?;
        Ok(df)
    }

    /// 创建空的 DataFrame
    pub fn empty(n_rows: usize) -> Self {
        Self {
            columns: Vec::new(),
            data: Vec::new(),
            n_rows,
        }
    }

    /// 验证 DataFrame 结构
    pub fn validate(&self) -> Result<(), String> {
        // 检查列名和数据长度一致
        if self.columns.len() != self.data.len() {
            return Err(format!(
                "DataFrame: columns length {} != data length {}",
                self.columns.len(),
                self.data.len()
            ));
        }

        // 检查每列的行数
        for (i, array) in self.data.iter().enumerate() {
            if array.len() != self.n_rows {
                return Err(format!(
                    "DataFrame: column '{}' has {} rows, expected {}",
                    self.columns.get(i).unwrap_or(&"unknown".to_string()),
                    array.len(),
                    self.n_rows
                ));
            }
        }

        Ok(())
    }

    /// 获取列数
    pub fn n_cols(&self) -> usize {
        self.columns.len()
    }

    /// 获取列的索引
    pub fn column_index(&self, name: &str) -> Option<usize> {
        self.columns.iter().position(|c| c == name)
    }

    /// 获取列数据
    pub fn column(&self, name: &str) -> Option<&ArrayRef> {
        self.column_index(name).and_then(|i| self.data.get(i))
    }

    /// 获取所有列的迭代器
    pub fn columns(&self) -> impl Iterator<Item = (&String, &ArrayRef)> {
        self.columns.iter().zip(self.data.iter())
    }
}

impl fmt::Display for DataFrame {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "DataFrame({} rows × {} columns)",
            self.n_rows,
            self.columns.len()
        )
    }
}

/// 分类数据映射（用于 Pandas Categorical 和 R Factor）
///
/// 存储整数编码和类别标签的映射关系
#[derive(Debug, Clone, PartialEq)]
pub struct CategoricalMapping {
    /// 类别标签（按编码顺序）
    pub categories: Vec<String>,
    /// 是否有序
    pub ordered: bool,
}

impl CategoricalMapping {
    /// 创建新的分类映射
    pub fn new(categories: Vec<String>, ordered: bool) -> Self {
        Self {
            categories,
            ordered,
        }
    }

    /// 获取类别数量
    pub fn n_categories(&self) -> usize {
        self.categories.len()
    }

    /// 根据编码获取类别标签
    pub fn get_category(&self, code: usize) -> Option<&str> {
        self.categories.get(code).map(|s| s.as_str())
    }

    /// 根据标签获取编码
    pub fn get_code(&self, category: &str) -> Option<usize> {
        self.categories.iter().position(|c| c == category)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow::array::{Int64Array, StringArray};
    use std::sync::Arc;

    #[test]
    fn test_dataframe_creation() {
        let col1 = Arc::new(Int64Array::from(vec![1, 2, 3])) as ArrayRef;
        let col2 = Arc::new(StringArray::from(vec!["a", "b", "c"])) as ArrayRef;

        let df = DataFrame::new(
            vec!["col1".to_string(), "col2".to_string()],
            vec![col1, col2],
            3,
        )
        .unwrap();

        assert_eq!(df.n_rows, 3);
        assert_eq!(df.n_cols(), 2);
    }

    #[test]
    fn test_dataframe_validation() {
        // 列名和数据长度不匹配
        let col1 = Arc::new(Int64Array::from(vec![1, 2, 3])) as ArrayRef;
        let result = DataFrame::new(vec!["col1".to_string(), "col2".to_string()], vec![col1], 3);
        assert!(result.is_err());
    }

    #[test]
    fn test_dataframe_row_mismatch() {
        // 行数不匹配
        let col1 = Arc::new(Int64Array::from(vec![1, 2])) as ArrayRef;
        let result = DataFrame::new(vec!["col1".to_string()], vec![col1], 3);
        assert!(result.is_err());
    }

    #[test]
    fn test_dataframe_column_access() {
        let col1 = Arc::new(Int64Array::from(vec![1, 2, 3])) as ArrayRef;
        let df = DataFrame::new(vec!["col1".to_string()], vec![col1], 3).unwrap();

        assert!(df.column("col1").is_some());
        assert!(df.column("col2").is_none());
        assert_eq!(df.column_index("col1"), Some(0));
    }

    #[test]
    fn test_empty_dataframe() {
        let df = DataFrame::empty(10);
        assert_eq!(df.n_rows, 10);
        assert_eq!(df.n_cols(), 0);
        assert_eq!(df.validate(), Ok(()));
    }

    #[test]
    fn test_categorical_mapping() {
        let mapping = CategoricalMapping::new(
            vec!["A".to_string(), "B".to_string(), "C".to_string()],
            false,
        );

        assert_eq!(mapping.n_categories(), 3);
        assert_eq!(mapping.get_category(0), Some("A"));
        assert_eq!(mapping.get_category(1), Some("B"));
        assert_eq!(mapping.get_code("B"), Some(1));
        assert_eq!(mapping.get_code("D"), None);
    }

    #[test]
    fn test_categorical_ordered() {
        let ordered = CategoricalMapping::new(
            vec!["low".to_string(), "medium".to_string(), "high".to_string()],
            true,
        );
        assert!(ordered.ordered);

        let unordered = CategoricalMapping::new(vec!["A".to_string(), "B".to_string()], false);
        assert!(!unordered.ordered);
    }
}
