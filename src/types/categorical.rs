//! 类别数据映射

/// 类别数据映射
#[derive(Debug, Clone)]
pub struct CategoricalMapping {
    pub categories: Vec<String>,
    pub codes: Vec<i32>,
    pub ordered: bool,
}
