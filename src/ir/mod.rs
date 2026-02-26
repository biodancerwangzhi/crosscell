//! CrossCell IR (Intermediate Representation)
//!
//! 定义中间表示层，用于在 Python AnnData 和 R Seurat 之间转换数据。

pub mod embeddings;
pub mod expression;
pub mod metadata;
pub mod spatial;
pub mod unstructured;

use std::collections::HashMap;
use std::fmt;

// 重新导出常用类型
pub use embeddings::Embedding;
pub use expression::{
    BackendType, CachePolicy, ChunkedMatrix, DenseMatrix, ExpressionMatrix, LazyMatrix,
    SparseFormat, SparseMatrixCSC, SparseMatrixCSR,
};
pub use metadata::{CategoricalMapping, DataFrame};
pub use spatial::{SpatialData, SpatialImage};
pub use unstructured::UnstructuredValue;

/// 主 IR 结构，表示单细胞数据集
///
/// 对应：
/// - Python: AnnData 对象
/// - R: Seurat 对象
#[derive(Debug, Clone)]
pub struct SingleCellData {
    /// 主表达矩阵（必需）
    /// 对应 AnnData.X 或 Seurat@assays[[active]]@data
    pub expression: ExpressionMatrix,

    /// 多层表达矩阵（可选）
    /// 对应 AnnData.layers 或 Seurat 多个 Assays
    pub layers: Option<HashMap<String, ExpressionMatrix>>,

    /// 细胞元数据（必需）
    /// 对应 AnnData.obs 或 Seurat@meta.data
    pub cell_metadata: DataFrame,

    /// 基因元数据（必需）
    /// 对应 AnnData.var 或 Seurat@assays[[x]]@meta.features
    pub gene_metadata: DataFrame,

    /// 降维嵌入（可选）
    /// 对应 AnnData.obsm 或 Seurat@reductions
    pub embeddings: Option<HashMap<String, Embedding>>,

    /// 细胞-细胞成对关系矩阵（可选）
    /// 对应 AnnData.obsp 或 Seurat@graphs
    pub cell_pairwise: Option<HashMap<String, PairwiseMatrix>>,

    /// 基因-基因成对关系矩阵（可选）
    /// 对应 AnnData.varp
    pub gene_pairwise: Option<HashMap<String, PairwiseMatrix>>,

    /// 空间数据（可选）
    /// 对应 AnnData.obsm['spatial'] 或 Seurat@images
    pub spatial: Option<SpatialData>,

    /// 基因加载矩阵（可选）
    /// 对应 AnnData.varm 或 Seurat@reductions[[x]]@feature.loadings
    /// key = reduction 名称（如 "PCs"），value = Embedding (n_genes × n_components)
    pub gene_loadings: Option<HashMap<String, Embedding>>,

    /// 非结构化数据（可选）
    /// 对应 AnnData.uns 或 Seurat@misc
    pub unstructured: Option<HashMap<String, UnstructuredValue>>,

    /// 数据集元数据
    pub metadata: DatasetMetadata,
}

impl SingleCellData {
    /// 创建新的单细胞数据集
    pub fn new(
        expression: ExpressionMatrix,
        cell_metadata: DataFrame,
        gene_metadata: DataFrame,
        metadata: DatasetMetadata,
    ) -> Result<Self, String> {
        let data = Self {
            expression,
            layers: None,
            cell_metadata,
            gene_metadata,
            embeddings: None,
            cell_pairwise: None,
            gene_pairwise: None,
            spatial: None,
            gene_loadings: None,
            unstructured: None,
            metadata,
        };
        data.validate()?;
        Ok(data)
    }

    /// 验证数据集的一致性
    pub fn validate(&self) -> Result<(), String> {
        let (n_cells, n_genes) = self.expression.shape();

        // 验证元数据维度
        if self.metadata.n_cells != n_cells {
            return Err(format!(
                "SingleCellData: metadata.n_cells {} != expression n_rows {}",
                self.metadata.n_cells, n_cells
            ));
        }

        if self.metadata.n_genes != n_genes {
            return Err(format!(
                "SingleCellData: metadata.n_genes {} != expression n_cols {}",
                self.metadata.n_genes, n_genes
            ));
        }

        // 验证细胞元数据行数
        if self.cell_metadata.n_rows != n_cells {
            return Err(format!(
                "SingleCellData: cell_metadata has {} rows, expected {}",
                self.cell_metadata.n_rows, n_cells
            ));
        }

        // 验证基因元数据行数
        if self.gene_metadata.n_rows != n_genes {
            return Err(format!(
                "SingleCellData: gene_metadata has {} rows, expected {}",
                self.gene_metadata.n_rows, n_genes
            ));
        }

        // 验证 layers 维度
        if let Some(ref layers) = self.layers {
            for (name, layer) in layers {
                let (layer_rows, layer_cols) = layer.shape();
                if layer_rows != n_cells || layer_cols != n_genes {
                    return Err(format!(
                        "SingleCellData: layer '{}' has shape ({}, {}), expected ({}, {})",
                        name, layer_rows, layer_cols, n_cells, n_genes
                    ));
                }
            }
        }

        // 验证嵌入维度
        if let Some(ref embeddings) = self.embeddings {
            for (name, emb) in embeddings {
                if emb.n_rows != n_cells {
                    return Err(format!(
                        "SingleCellData: embedding '{}' has {} rows, expected {}",
                        name, emb.n_rows, n_cells
                    ));
                }
            }
        }

        // 验证细胞成对矩阵维度
        if let Some(ref pairwise) = self.cell_pairwise {
            for (name, pw) in pairwise {
                let (pw_rows, pw_cols) = pw.matrix.shape();
                if pw_rows != n_cells || pw_cols != n_cells {
                    return Err(format!(
                        "SingleCellData: cell_pairwise '{}' has shape ({}, {}), expected ({}, {})",
                        name, pw_rows, pw_cols, n_cells, n_cells
                    ));
                }
            }
        }

        // 验证基因成对矩阵维度
        if let Some(ref pairwise) = self.gene_pairwise {
            for (name, pw) in pairwise {
                let (pw_rows, pw_cols) = pw.matrix.shape();
                if pw_rows != n_genes || pw_cols != n_genes {
                    return Err(format!(
                        "SingleCellData: gene_pairwise '{}' has shape ({}, {}), expected ({}, {})",
                        name, pw_rows, pw_cols, n_genes, n_genes
                    ));
                }
            }
        }

        // 验证基因加载矩阵维度 (varm)
        if let Some(ref gene_loadings) = self.gene_loadings {
            for (name, loading) in gene_loadings {
                if loading.n_rows != n_genes {
                    return Err(format!(
                        "SingleCellData: gene_loading '{}' has {} rows, expected {} (n_genes)",
                        name, loading.n_rows, n_genes
                    ));
                }
            }
        }

        Ok(())
    }
}

impl fmt::Display for SingleCellData {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "SingleCellData({} cells × {} genes, source: {})",
            self.metadata.n_cells, self.metadata.n_genes, self.metadata.source_format
        )
    }
}

/// 成对关系矩阵（用于细胞-细胞或基因-基因关系）
///
/// 对应：
/// - AnnData.obsp (细胞-细胞)
/// - AnnData.varp (基因-基因)
/// - Seurat@graphs (细胞-细胞)
#[derive(Debug, Clone, PartialEq)]
pub struct PairwiseMatrix {
    /// 矩阵名称（如 "connectivities", "distances"）
    pub name: String,
    /// 成对关系矩阵（通常是稀疏方阵）
    pub matrix: ExpressionMatrix,
}

impl PairwiseMatrix {
    /// 创建新的成对矩阵
    pub fn new(name: String, matrix: ExpressionMatrix) -> Result<Self, String> {
        let pw = Self { name, matrix };
        pw.validate()?;
        Ok(pw)
    }

    /// 验证成对矩阵是方阵
    pub fn validate(&self) -> Result<(), String> {
        let (n_rows, n_cols) = self.matrix.shape();
        if n_rows != n_cols {
            return Err(format!(
                "PairwiseMatrix '{}': not a square matrix ({} × {})",
                self.name, n_rows, n_cols
            ));
        }
        Ok(())
    }
}

impl fmt::Display for PairwiseMatrix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (n, _) = self.matrix.shape();
        write!(f, "PairwiseMatrix('{}', {} × {})", self.name, n, n)
    }
}

/// 数据集元数据
#[derive(Debug, Clone, PartialEq)]
pub struct DatasetMetadata {
    /// 细胞数量
    pub n_cells: usize,
    /// 基因数量
    pub n_genes: usize,
    /// 源格式（"anndata" 或 "seurat"）
    pub source_format: String,
    /// 源版本（可选）
    pub source_version: Option<String>,
    /// 活跃 Assay 名称（用于 Seurat 多 Assay 支持）
    pub active_assay: Option<String>,
}

impl DatasetMetadata {
    /// 创建新的数据集元数据
    pub fn new(n_cells: usize, n_genes: usize, source_format: String) -> Self {
        Self {
            n_cells,
            n_genes,
            source_format,
            source_version: None,
            active_assay: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_cell_data_creation() {
        let expr = ExpressionMatrix::Dense(DenseMatrix::zeros(10, 20));
        let cell_meta = DataFrame::empty(10);
        let gene_meta = DataFrame::empty(20);
        let metadata = DatasetMetadata::new(10, 20, "test".to_string());

        let data = SingleCellData::new(expr, cell_meta, gene_meta, metadata).unwrap();
        assert_eq!(data.metadata.n_cells, 10);
        assert_eq!(data.metadata.n_genes, 20);
    }

    #[test]
    fn test_single_cell_data_validation() {
        // 元数据维度不匹配
        let expr = ExpressionMatrix::Dense(DenseMatrix::zeros(10, 20));
        let cell_meta = DataFrame::empty(5); // 错误：应该是 10
        let gene_meta = DataFrame::empty(20);
        let metadata = DatasetMetadata::new(10, 20, "test".to_string());

        let result = SingleCellData::new(expr, cell_meta, gene_meta, metadata);
        assert!(result.is_err());
    }

    #[test]
    fn test_pairwise_matrix_square() {
        // 方阵
        let matrix = ExpressionMatrix::Dense(DenseMatrix::zeros(100, 100));
        let pw = PairwiseMatrix::new("connectivities".to_string(), matrix).unwrap();
        assert_eq!(pw.name, "connectivities");
    }

    #[test]
    fn test_pairwise_matrix_not_square() {
        // 非方阵
        let matrix = ExpressionMatrix::Dense(DenseMatrix::zeros(100, 50));
        let result = PairwiseMatrix::new("distances".to_string(), matrix);
        assert!(result.is_err());
    }

    #[test]
    fn test_single_cell_data_with_layers() {
        let expr = ExpressionMatrix::Dense(DenseMatrix::zeros(10, 20));
        let cell_meta = DataFrame::empty(10);
        let gene_meta = DataFrame::empty(20);
        let metadata = DatasetMetadata::new(10, 20, "test".to_string());

        let mut data = SingleCellData::new(expr, cell_meta, gene_meta, metadata).unwrap();

        // 添加 layers
        let mut layers = HashMap::new();
        layers.insert(
            "counts".to_string(),
            ExpressionMatrix::Dense(DenseMatrix::zeros(10, 20)),
        );
        layers.insert(
            "log1p".to_string(),
            ExpressionMatrix::Dense(DenseMatrix::zeros(10, 20)),
        );
        data.layers = Some(layers);

        assert_eq!(data.validate(), Ok(()));
    }

    #[test]
    fn test_single_cell_data_with_pairwise() {
        let expr = ExpressionMatrix::Dense(DenseMatrix::zeros(10, 20));
        let cell_meta = DataFrame::empty(10);
        let gene_meta = DataFrame::empty(20);
        let metadata = DatasetMetadata::new(10, 20, "test".to_string());

        let mut data = SingleCellData::new(expr, cell_meta, gene_meta, metadata).unwrap();

        // 添加细胞成对矩阵
        let mut cell_pairwise = HashMap::new();
        let pw_matrix = ExpressionMatrix::Dense(DenseMatrix::zeros(10, 10));
        cell_pairwise.insert(
            "connectivities".to_string(),
            PairwiseMatrix::new("connectivities".to_string(), pw_matrix).unwrap(),
        );
        data.cell_pairwise = Some(cell_pairwise);

        assert_eq!(data.validate(), Ok(()));
    }
}
