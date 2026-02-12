//! 必需字段验证
//!
//! 验证 .h5ad 和 .rds 文件包含所有必需的字段

use crate::error::{CrossCellError, Result};
use crate::ir::SingleCellData;

/// 验证 SingleCellData 包含所有必需字段
pub fn validate_required_fields(data: &SingleCellData) -> Result<()> {
    let mut missing_fields = Vec::new();

    // 验证表达矩阵存在
    let (n_cells, n_genes) = data.expression.shape();
    if n_cells == 0 || n_genes == 0 {
        missing_fields.push("expression matrix (empty)".to_string());
    }

    // 验证细胞元数据存在
    if data.cell_metadata.n_rows == 0 {
        missing_fields.push("cell_metadata (empty)".to_string());
    }

    // 验证基因元数据存在
    if data.gene_metadata.n_rows == 0 {
        missing_fields.push("gene_metadata (empty)".to_string());
    }

    // 验证元数据信息
    if data.metadata.n_cells == 0 {
        missing_fields.push("metadata.n_cells (zero)".to_string());
    }

    if data.metadata.n_genes == 0 {
        missing_fields.push("metadata.n_genes (zero)".to_string());
    }

    if data.metadata.source_format.is_empty() {
        missing_fields.push("metadata.source_format (empty)".to_string());
    }

    if !missing_fields.is_empty() {
        return Err(CrossCellError::missing_fields(missing_fields));
    }

    Ok(())
}

/// 验证 AnnData 必需字段
///
/// AnnData 必需字段：
/// - /X (表达矩阵)
/// - /obs (细胞元数据)
/// - /var (基因元数据)
pub fn validate_anndata_required_fields(data: &SingleCellData) -> Result<()> {
    let mut missing_fields = Vec::new();

    // 验证表达矩阵
    let (n_cells, n_genes) = data.expression.shape();
    if n_cells == 0 || n_genes == 0 {
        missing_fields.push("/X (expression matrix)".to_string());
    }

    // 验证细胞元数据
    if data.cell_metadata.n_rows == 0 {
        missing_fields.push("/obs (cell metadata)".to_string());
    }

    // 验证基因元数据
    if data.gene_metadata.n_rows == 0 {
        missing_fields.push("/var (gene metadata)".to_string());
    }

    if !missing_fields.is_empty() {
        return Err(CrossCellError::missing_fields(missing_fields));
    }

    Ok(())
}

/// 验证 Seurat 必需字段
///
/// Seurat 必需字段：
/// - assays (至少一个 Assay)
/// - meta.data (细胞元数据)
/// - 表达矩阵（counts 或 data）
pub fn validate_seurat_required_fields(data: &SingleCellData) -> Result<()> {
    let mut missing_fields = Vec::new();

    // 验证表达矩阵
    let (n_cells, n_genes) = data.expression.shape();
    if n_cells == 0 || n_genes == 0 {
        missing_fields.push("assay expression matrix".to_string());
    }

    // 验证细胞元数据
    if data.cell_metadata.n_rows == 0 {
        missing_fields.push("meta.data (cell metadata)".to_string());
    }

    // 验证基因元数据（Seurat 中基因名存储在 features）
    if data.gene_metadata.n_rows == 0 {
        missing_fields.push("features (gene metadata)".to_string());
    }

    // 验证 active_assay 存在
    if data.metadata.active_assay.is_none() {
        missing_fields.push("active_assay".to_string());
    }

    if !missing_fields.is_empty() {
        return Err(CrossCellError::missing_fields(missing_fields));
    }

    Ok(())
}

/// 生成必需字段报告
pub fn required_fields_report(data: &SingleCellData) -> String {
    let mut report = String::from("Required Fields Report:\n");

    // 检查表达矩阵
    let (n_cells, n_genes) = data.expression.shape();
    if n_cells > 0 && n_genes > 0 {
        report.push_str("✓ Expression matrix present\n");
    } else {
        report.push_str("✗ Expression matrix missing or empty\n");
    }

    // 检查细胞元数据
    if data.cell_metadata.n_rows > 0 {
        report.push_str("✓ Cell metadata present\n");
    } else {
        report.push_str("✗ Cell metadata missing or empty\n");
    }

    // 检查基因元数据
    if data.gene_metadata.n_rows > 0 {
        report.push_str("✓ Gene metadata present\n");
    } else {
        report.push_str("✗ Gene metadata missing or empty\n");
    }

    // 检查元数据信息
    if !data.metadata.source_format.is_empty() {
        report.push_str(&format!(
            "✓ Source format: {}\n",
            data.metadata.source_format
        ));
    } else {
        report.push_str("✗ Source format missing\n");
    }

    if let Some(ref version) = data.metadata.source_version {
        report.push_str(&format!("✓ Source version: {}\n", version));
    }

    if let Some(ref assay) = data.metadata.active_assay {
        report.push_str(&format!("✓ Active assay: {}\n", assay));
    }

    // 检查可选字段
    if data.embeddings.is_some() {
        report.push_str("✓ Embeddings present\n");
    }

    if data.layers.is_some() {
        report.push_str("✓ Layers present\n");
    }

    if data.spatial.is_some() {
        report.push_str("✓ Spatial data present\n");
    }

    report
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::{DatasetMetadata, DenseMatrix, DataFrame, ExpressionMatrix};

    fn create_valid_data() -> SingleCellData {
        SingleCellData {
            expression: ExpressionMatrix::Dense(DenseMatrix {
                data: vec![1.0; 100 * 50],
                n_rows: 100,
                n_cols: 50,
            }),
            cell_metadata: DataFrame::empty(100),
            gene_metadata: DataFrame::empty(50),
            embeddings: None,
            layers: None,
            cell_pairwise: None,
            gene_pairwise: None,
            spatial: None,
            unstructured: None,
            metadata: DatasetMetadata {
                n_cells: 100,
                n_genes: 50,
                source_format: "test".to_string(),
                source_version: Some("1.0".to_string()),
                active_assay: Some("RNA".to_string()),
            },
        }
    }

    #[test]
    fn test_validate_required_fields_success() {
        let data = create_valid_data();
        assert!(validate_required_fields(&data).is_ok());
    }

    #[test]
    fn test_validate_required_fields_empty_expression() {
        let mut data = create_valid_data();
        data.expression = ExpressionMatrix::Dense(DenseMatrix {
            data: vec![],
            n_rows: 0,
            n_cols: 0,
        });

        let result = validate_required_fields(&data);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(matches!(err, CrossCellError::MissingFields { .. }));
    }

    #[test]
    fn test_validate_required_fields_empty_metadata() {
        let mut data = create_valid_data();
        data.cell_metadata = DataFrame::empty(0);

        let result = validate_required_fields(&data);
        assert!(result.is_err());
    }

    #[test]
    fn test_validate_anndata_required_fields() {
        let data = create_valid_data();
        assert!(validate_anndata_required_fields(&data).is_ok());
    }

    #[test]
    fn test_validate_seurat_required_fields() {
        let data = create_valid_data();
        assert!(validate_seurat_required_fields(&data).is_ok());
    }

    #[test]
    fn test_validate_seurat_missing_active_assay() {
        let mut data = create_valid_data();
        data.metadata.active_assay = None;

        let result = validate_seurat_required_fields(&data);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("active_assay"));
    }

    #[test]
    fn test_required_fields_report() {
        let data = create_valid_data();
        let report = required_fields_report(&data);

        assert!(report.contains("Expression matrix present"));
        assert!(report.contains("Cell metadata present"));
        assert!(report.contains("Gene metadata present"));
        assert!(report.contains("Source format: test"));
    }
}
