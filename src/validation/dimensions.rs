//! 维度验证
//!
//! 验证表达矩阵维度与元数据的一致性

use crate::error::{CrossCellError, Result};
use crate::ir::{ExpressionMatrix, SingleCellData};

/// 验证 SingleCellData 的维度一致性
pub fn validate_dimensions(data: &SingleCellData) -> Result<()> {
    let (n_cells, n_genes) = data.expression.shape();

    // 验证细胞元数据
    let meta_rows = data.cell_metadata.n_rows;
    if meta_rows != n_cells {
        return Err(CrossCellError::dimension_mismatch_detailed(
            "cell_metadata",
            format!("{} cells", n_cells),
            format!("{} rows", meta_rows),
        ));
    }

    // 验证基因元数据
    let meta_rows = data.gene_metadata.n_rows;
    if meta_rows != n_genes {
        return Err(CrossCellError::dimension_mismatch_detailed(
            "gene_metadata",
            format!("{} genes", n_genes),
            format!("{} rows", meta_rows),
        ));
    }

    // 验证 layers
    if let Some(ref layers) = data.layers {
        for (layer_name, layer_matrix) in layers {
            let (layer_cells, layer_genes) = layer_matrix.shape();
            if layer_cells != n_cells || layer_genes != n_genes {
                return Err(CrossCellError::dimension_mismatch_detailed(
                    format!("layer '{}'", layer_name),
                    format!("{}x{}", n_cells, n_genes),
                    format!("{}x{}", layer_cells, layer_genes),
                ));
            }
        }
    }

    // 验证 embeddings
    if let Some(ref embeddings) = data.embeddings {
        for (emb_name, embedding) in embeddings {
            let emb_rows = embedding.n_rows;
            if emb_rows != n_cells {
                return Err(CrossCellError::dimension_mismatch_detailed(
                    format!("embedding '{}'", emb_name),
                    format!("{} cells", n_cells),
                    format!("{} rows", emb_rows),
                ));
            }
        }
    }

    // 验证 cell_pairwise 矩阵
    if let Some(ref cell_pairwise) = data.cell_pairwise {
        for (pw_name, pw_matrix) in cell_pairwise {
            let (pw_rows, pw_cols) = pw_matrix.matrix.shape();
            if pw_rows != n_cells || pw_cols != n_cells {
                return Err(CrossCellError::dimension_mismatch_detailed(
                    format!("cell_pairwise '{}'", pw_name),
                    format!("{}x{}", n_cells, n_cells),
                    format!("{}x{}", pw_rows, pw_cols),
                ));
            }
        }
    }

    // 验证 gene_pairwise 矩阵
    if let Some(ref gene_pairwise) = data.gene_pairwise {
        for (pw_name, pw_matrix) in gene_pairwise {
            let (pw_rows, pw_cols) = pw_matrix.matrix.shape();
            if pw_rows != n_genes || pw_cols != n_genes {
                return Err(CrossCellError::dimension_mismatch_detailed(
                    format!("gene_pairwise '{}'", pw_name),
                    format!("{}x{}", n_genes, n_genes),
                    format!("{}x{}", pw_rows, pw_cols),
                ));
            }
        }
    }

    Ok(())
}

/// 验证表达矩阵维度
pub fn validate_expression_matrix_dimensions(
    matrix: &ExpressionMatrix,
    expected_cells: usize,
    expected_genes: usize,
) -> Result<()> {
    let (actual_cells, actual_genes) = matrix.shape();

    if actual_cells != expected_cells {
        return Err(CrossCellError::dimension_mismatch_detailed(
            "expression_matrix (cells)",
            format!("{}", expected_cells),
            format!("{}", actual_cells),
        ));
    }

    if actual_genes != expected_genes {
        return Err(CrossCellError::dimension_mismatch_detailed(
            "expression_matrix (genes)",
            format!("{}", expected_genes),
            format!("{}", actual_genes),
        ));
    }

    Ok(())
}

/// 生成维度报告
pub fn dimension_report(data: &SingleCellData) -> String {
    let (n_cells, n_genes) = data.expression.shape();

    let mut report = format!(
        "Dimension Report:\n\
         - Expression matrix: {} cells × {} genes\n",
        n_cells, n_genes
    );

    report.push_str(&format!("- Cell metadata: {} rows\n", data.cell_metadata.n_rows));
    report.push_str(&format!("- Gene metadata: {} rows\n", data.gene_metadata.n_rows));

    if let Some(ref layers) = data.layers {
        if !layers.is_empty() {
            report.push_str(&format!("- Layers: {}\n", layers.len()));
            for (name, layer) in layers {
                let (lc, lg) = layer.shape();
                report.push_str(&format!("  - {}: {}x{}\n", name, lc, lg));
            }
        }
    }

    if let Some(ref embeddings) = data.embeddings {
        if !embeddings.is_empty() {
            report.push_str(&format!("- Embeddings: {}\n", embeddings.len()));
            for (name, emb) in embeddings {
                report.push_str(&format!(
                    "  - {}: {} cells × {} components\n",
                    name, emb.n_rows, emb.n_cols
                ));
            }
        }
    }

    if let Some(ref cell_pairwise) = data.cell_pairwise {
        if !cell_pairwise.is_empty() {
            report.push_str(&format!("- Cell pairwise: {}\n", cell_pairwise.len()));
        }
    }

    if let Some(ref gene_pairwise) = data.gene_pairwise {
        if !gene_pairwise.is_empty() {
            report.push_str(&format!("- Gene pairwise: {}\n", gene_pairwise.len()));
        }
    }

    if let Some(ref gene_loadings) = data.gene_loadings {
        if !gene_loadings.is_empty() {
            report.push_str(&format!("- Gene loadings (varm): {}\n", gene_loadings.len()));
            for (name, loading) in gene_loadings {
                report.push_str(&format!(
                    "  - {}: {} genes × {} components\n",
                    name, loading.n_rows, loading.n_cols
                ));
            }
        }
    }

    report
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::{DatasetMetadata, DenseMatrix, DataFrame, ExpressionMatrix};

    fn create_test_data(n_cells: usize, n_genes: usize) -> SingleCellData {
        SingleCellData {
            expression: ExpressionMatrix::Dense(DenseMatrix {
                data: vec![0.0; n_cells * n_genes],
                n_rows: n_cells,
                n_cols: n_genes,
            }),
            cell_metadata: DataFrame::empty(n_cells),
            gene_metadata: DataFrame::empty(n_genes),
            embeddings: None,
            layers: None,
            cell_pairwise: None,
            gene_pairwise: None,
            spatial: None,
            gene_loadings: None,
            unstructured: None,
            metadata: DatasetMetadata {
                n_cells,
                n_genes,
                source_format: "test".to_string(),
                source_version: None,
                active_assay: None,
            },
        }
    }

    #[test]
    fn test_validate_dimensions_success() {
        let data = create_test_data(100, 50);
        assert!(validate_dimensions(&data).is_ok());
    }

    #[test]
    fn test_validate_dimensions_cell_metadata_mismatch() {
        let mut data = create_test_data(100, 50);
        data.cell_metadata = DataFrame::empty(90);

        let result = validate_dimensions(&data);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(matches!(
            err,
            CrossCellError::DimensionMismatchDetailed { .. }
        ));
        assert!(err.to_string().contains("cell_metadata"));
    }

    #[test]
    fn test_validate_dimensions_gene_metadata_mismatch() {
        let mut data = create_test_data(100, 50);
        data.gene_metadata = DataFrame::empty(40);

        let result = validate_dimensions(&data);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("gene_metadata"));
    }

    #[test]
    fn test_validate_expression_matrix_dimensions() {
        let matrix = ExpressionMatrix::Dense(DenseMatrix {
            data: vec![0.0; 100 * 50],
            n_rows: 100,
            n_cols: 50,
        });

        assert!(validate_expression_matrix_dimensions(&matrix, 100, 50).is_ok());
        assert!(validate_expression_matrix_dimensions(&matrix, 90, 50).is_err());
        assert!(validate_expression_matrix_dimensions(&matrix, 100, 40).is_err());
    }

    #[test]
    fn test_dimension_report() {
        let data = create_test_data(100, 50);
        let report = dimension_report(&data);

        assert!(report.contains("100 cells"));
        assert!(report.contains("50 genes"));
        assert!(report.contains("Cell metadata: 100 rows"));
        assert!(report.contains("Gene metadata: 50 rows"));
    }
}
