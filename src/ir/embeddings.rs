//! 降维嵌入 IR 类型
//!
//! 存储降维结果（PCA、UMAP、t-SNE 等）

use std::fmt;

/// 降维嵌入
///
/// 对应：
/// - AnnData.obsm['X_pca'], AnnData.obsm['X_umap']
/// - Seurat@reductions$pca, Seurat@reductions$umap
#[derive(Debug, Clone, PartialEq)]
pub struct Embedding {
    /// 嵌入名称（如 "X_pca", "X_umap"）
    pub name: String,
    /// 嵌入数据（行优先，cells × components）
    pub data: Vec<f64>,
    /// 行数（细胞数）
    pub n_rows: usize,
    /// 列数（成分数）
    pub n_cols: usize,
}

impl Embedding {
    /// 创建新的嵌入
    pub fn new(name: String, data: Vec<f64>, n_rows: usize, n_cols: usize) -> Result<Self, String> {
        let embedding = Self {
            name,
            data,
            n_rows,
            n_cols,
        };
        embedding.validate()?;
        Ok(embedding)
    }

    /// 验证嵌入维度一致性
    pub fn validate(&self) -> Result<(), String> {
        let expected_len = self.n_rows * self.n_cols;
        if self.data.len() != expected_len {
            return Err(format!(
                "Embedding '{}': data length {} does not match dimensions {} × {} = {}",
                self.name,
                self.data.len(),
                self.n_rows,
                self.n_cols,
                expected_len
            ));
        }
        Ok(())
    }

    /// 获取指定位置的值
    pub fn get(&self, row: usize, col: usize) -> Option<f64> {
        if row >= self.n_rows || col >= self.n_cols {
            return None;
        }
        Some(self.data[row * self.n_cols + col])
    }
}

impl fmt::Display for Embedding {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Embedding('{}', {} × {})",
            self.name, self.n_rows, self.n_cols
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_embedding_creation() {
        // 3 细胞 × 2 成分
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let emb = Embedding::new("X_pca".to_string(), data, 3, 2).unwrap();

        assert_eq!(emb.name, "X_pca");
        assert_eq!(emb.n_rows, 3);
        assert_eq!(emb.n_cols, 2);
        assert_eq!(emb.get(0, 0), Some(1.0));
        assert_eq!(emb.get(2, 1), Some(6.0));
    }

    #[test]
    fn test_embedding_validation() {
        // 维度不匹配
        let data = vec![1.0, 2.0, 3.0];
        let result = Embedding::new("X_umap".to_string(), data, 3, 2);
        assert!(result.is_err());
    }

    #[test]
    fn test_embedding_get_out_of_bounds() {
        let data = vec![1.0, 2.0, 3.0, 4.0];
        let emb = Embedding::new("X_tsne".to_string(), data, 2, 2).unwrap();

        assert_eq!(emb.get(0, 0), Some(1.0));
        assert_eq!(emb.get(2, 0), None); // 超出行范围
        assert_eq!(emb.get(0, 2), None); // 超出列范围
    }
}
