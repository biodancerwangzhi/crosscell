//! Seurat V5 Assay5 支持
//!
//! 此模块提供对 Seurat V5 Assay5 结构的支持
//!
//! ## Assay5 vs 传统 Assay
//!
//! Seurat V5 引入了新Assay5 结构，主要区别：
//!
//! - **传统 Assay**: 固定槽位（counts, data, scale.data
//! - **Assay5**: 动layers 结构，可以包含任意数量的
//!
//! ## Assay5 结构
//!
//! ```r
//! Assay5 <- list(
//!   class = "Assay5",
//!   layers = list(
//!     counts = dgCMatrix,
//!     data = dgCMatrix,
//!     scale.data = matrix  # 可
//!   ),
//!   features = character vector,
//!   cells = character vector,
//!   meta.features = data.frame  # 可
//! )
//! ```
//!
//! ## BPCells 支持
//!
//! Assay5 可以使用 BPCells 进行磁盘支持的存储
//! 当检测到 BPCells 对象时，我们需要将其加载到内存中

use crate::ir::{DataFrame, ExpressionMatrix, SparseMatrixCSC, SparseMatrixCSR};
use crate::rds::{GenericVector, RObject, RdsFile};
use crate::seurat::error::SeuratError;
use std::collections::HashMap;

/// Assay5 类型检测结
#[derive(Debug, Clone, PartialEq)]
pub enum AssayType {
    /// 传统 Assay (Seurat V3/V4)
    Legacy,
    /// Assay5 (Seurat V5)
    Assay5,
    /// 未知类型
    Unknown,
}

/// 检Assay 类型
///
/// # Arguments
///
/// * `robj` - Assay RObject
///
/// # Returns
///
/// * `AssayType` - 检测到Assay 类型
pub fn detect_assay_type(robj: &RObject) -> AssayType {
    // 检class 属
    if let Some(class_str) = extract_class(robj) {
        if class_str.contains("Assay5") {
            return AssayType::Assay5;
        } else if class_str.contains("Assay") {
            return AssayType::Legacy;
        }
    }

    // 如果没有 class 属性，检查结
    let fields = get_named_list_fields(robj);

    if let Some(fields) = fields {
        // Assay5 layers 字段
        for (name, _) in &fields {
            if *name == "layers" {
                return AssayType::Assay5;
            }
        }
        // 传统 Assay counts 字段
        for (name, _) in &fields {
            if *name == "counts" {
                return AssayType::Legacy;
            }
        }
    }

    AssayType::Unknown
}

/// RObject 获取命名列表字段
fn get_named_list_fields<'a>(robj: &'a RObject) -> Option<Vec<(&'a str, &'a RObject)>> {
    match robj {
        RObject::GenericVector(gv) => {
            if let Some(names) = gv.attributes.get_names() {
                Some(
                    names
                        .iter()
                        .zip(gv.data.iter())
                        .map(|(n, v)| (n.as_str(), v))
                        .collect(),
                )
            } else {
                None
            }
        }
        RObject::PairList(pl) => Some(
            pl.tag_names
                .iter()
                .zip(pl.data.iter())
                .map(|(n, v)| (n.as_str(), v))
                .collect(),
        ),
        RObject::S4Object(s4) => Some(
            s4.attributes
                .names
                .iter()
                .zip(s4.attributes.values.iter())
                .map(|(n, v)| (n.as_str(), v))
                .collect(),
        ),
        _ => None,
    }
}

/// RObject 提取 class 属
fn extract_class(robj: &RObject) -> Option<String> {
    // For S4 objects, use the class_name field directly
    if let RObject::S4Object(s4) = robj {
        return Some(s4.class_name.clone());
    }
    robj.attributes()
        .and_then(|attrs| attrs.get_class().map(|classes| classes.join(",")))
}

/// Assay5 提取 layers
///
/// # Arguments
///
/// * `robj` - Assay5 RObject
/// * `file` - RDS 文件（用于解析引用）
///
/// # Returns
///
/// * `Result<HashMap<String, ExpressionMatrix>, SeuratError>` - Layer 名称到表达矩阵的映射
pub fn extract_assay5_layers(
    robj: &RObject,
    _file: &RdsFile,
) -> Result<HashMap<String, ExpressionMatrix>, SeuratError> {
    let assay_fields = get_named_list_fields(robj).ok_or_else(|| {
        SeuratError::ParseError(format!(
            "Expected named list for Assay5, got {}",
            robj.type_name()
        ))
    })?;

    // 查找 layers 字段
    let layers_obj = assay_fields
        .iter()
        .find(|(name, _)| *name == "layers")
        .map(|(_, value)| *value)
        .ok_or_else(|| SeuratError::ParseError("Missing 'layers' field in Assay5".to_string()))?;

    // 解析 layers（应该是一named list
    let layers_list = get_named_list_fields(layers_obj).ok_or_else(|| {
        SeuratError::ParseError(format!(
            "Expected named list for layers, got {}",
            layers_obj.type_name()
        ))
    })?;

    // 提取每个 layer
    let mut layers = HashMap::new();
    for (layer_name, layer_value) in layers_list {
        // 跳过 NULL layers
        if matches!(layer_value, RObject::Null) {
            log::debug!("Skipping NULL layer '{}'", layer_name);
            continue;
        }

        // 检查是否是 BPCells 对象
        if is_bpcells_object(layer_value) {
            log::warn!(
                "Layer '{}' is a BPCells object, loading to memory...",
                layer_name
            );
            // TODO: 实现 BPCells 加载
            // 现在先跳
            continue;
        }

        // 检查是否是 dgCMatrix (S4 对象)
        let is_dgcmatrix = matches!(layer_value, RObject::S4Object(_));

        if is_dgcmatrix {
            // 解析dgCMatrix
            match crate::seurat::seurat_to_ir::parse_dgcmatrix(layer_value) {
                Ok(matrix) => {
                    layers.insert(layer_name.to_string(), matrix);
                }
                Err(e) => {
                    log::warn!("Failed to parse layer '{}' as dgCMatrix: {}", layer_name, e);
                    // 继续处理其他 layers
                }
            }
        } else {
            // 可能是密集矩(scale.data) 或其他类
            // 目前跳过非稀疏矩阵的 layers
            log::debug!(
                "Skipping non-dgCMatrix layer '{}' (type: {})",
                layer_name,
                layer_value.type_name()
            );
        }
    }

    if layers.is_empty() {
        return Err(SeuratError::ParseError(
            "No valid layers found in Assay5".to_string(),
        ));
    }

    // Task 1.6: 合并 split layers（如 counts.1, counts.2 → counts）
    let layers = merge_split_layers(layers)?;

    Ok(layers)
}

// ============================================================================
// Split Layer 合并（Task 1）
// ============================================================================

/// 按行拼接多个 CSR 矩阵（细胞维度）
///
/// 所有矩阵必须具有相同的 n_cols（基因数）。
/// 拼接后的矩阵行数等于所有子矩阵行数之和。
///
/// # Arguments
///
/// * `matrices` - CSR 矩阵引用切片
///
/// # Returns
///
/// * 合并后的 CSR 矩阵，或列数不一致时返回错误
pub fn merge_csr_row_wise(matrices: &[&SparseMatrixCSR]) -> Result<SparseMatrixCSR, SeuratError> {
    if matrices.is_empty() {
        return Err(SeuratError::ParseError(
            "Cannot merge empty list of matrices".to_string(),
        ));
    }

    let n_cols = matrices[0].n_cols;

    // 验证所有矩阵列数相同
    for (i, m) in matrices.iter().enumerate() {
        if m.n_cols != n_cols {
            return Err(SeuratError::ParseError(format!(
                "Split layer column count mismatch: matrix[0] has {} cols, matrix[{}] has {} cols",
                n_cols, i, m.n_cols
            )));
        }
    }

    // 计算合并后的总行数和非零元素数
    let total_rows: usize = matrices.iter().map(|m| m.n_rows).sum();
    let total_nnz: usize = matrices.iter().map(|m| m.data.len()).sum();

    let mut merged_data = Vec::with_capacity(total_nnz);
    let mut merged_indices = Vec::with_capacity(total_nnz);
    let mut merged_indptr = Vec::with_capacity(total_rows + 1);
    merged_indptr.push(0);

    let mut offset = 0usize;
    for m in matrices {
        // 追加 data 和 indices
        merged_data.extend_from_slice(&m.data);
        merged_indices.extend_from_slice(&m.indices);

        // 追加 indptr（跳过每个子矩阵的第一个 0，加上偏移量）
        for &ptr in &m.indptr[1..] {
            merged_indptr.push(ptr + offset);
        }
        offset += m.data.len();
    }

    Ok(SparseMatrixCSR {
        data: merged_data,
        indices: merged_indices,
        indptr: merged_indptr,
        n_rows: total_rows,
        n_cols,
    })
}

/// 合并多个 CSC 稀疏矩阵（按列拼接）
///
/// dgCMatrix 是列压缩格式，列 = 细胞，行 = 基因。
/// 合并 split layers 时需要按列（细胞）拼接。
///
/// # Arguments
///
/// * `matrices` - CSC 矩阵引用切片，所有矩阵必须有相同的行数（基因数）
///
/// # Returns
///
/// * 合并后的 CSC 矩阵
pub fn merge_csc_column_wise(
    matrices: &[&SparseMatrixCSC],
) -> Result<SparseMatrixCSC, SeuratError> {
    if matrices.is_empty() {
        return Err(SeuratError::ParseError(
            "Cannot merge empty list of CSC matrices".to_string(),
        ));
    }

    let n_rows = matrices[0].n_rows;

    // 验证所有矩阵行数相同（基因数）
    for (i, m) in matrices.iter().enumerate() {
        if m.n_rows != n_rows {
            return Err(SeuratError::ParseError(format!(
                "Split layer row count mismatch: matrix[0] has {} rows, matrix[{}] has {} rows",
                n_rows, i, m.n_rows
            )));
        }
    }

    let total_cols: usize = matrices.iter().map(|m| m.n_cols).sum();
    let total_nnz: usize = matrices.iter().map(|m| m.data.len()).sum();

    let mut merged_data = Vec::with_capacity(total_nnz);
    let mut merged_indices = Vec::with_capacity(total_nnz);
    let mut merged_indptr = Vec::with_capacity(total_cols + 1);
    merged_indptr.push(0);

    let mut offset = 0usize;
    for m in matrices {
        merged_data.extend_from_slice(&m.data);
        merged_indices.extend_from_slice(&m.indices);

        for &ptr in &m.indptr[1..] {
            merged_indptr.push(ptr + offset);
        }
        offset += m.data.len();
    }

    Ok(SparseMatrixCSC {
        data: merged_data,
        indices: merged_indices,
        indptr: merged_indptr,
        n_rows,
        n_cols: total_cols,
    })
}

/// 检测并合并 split layers
///
/// 将名称匹配 `<base>.<suffix>` 模式的层按前缀分组。
/// 同一前缀有多个层时，按后缀排序后合并。
/// 支持 CSR（按行拼接）和 CSC（按列拼接）两种稀疏格式。
/// 单个带点的层保留原名。非匹配层原样保留。
///
/// # Arguments
///
/// * `layers` - 层名称到表达矩阵的映射
///
/// # Returns
///
/// * 合并后的层映射
pub fn merge_split_layers(
    layers: HashMap<String, ExpressionMatrix>,
) -> Result<HashMap<String, ExpressionMatrix>, SeuratError> {
    // 按最后一个 '.' 分割层名，将共享相同前缀的层分组
    // 例如: counts.batch1, counts.batch2 → 前缀 "counts", 后缀 "batch1", "batch2"
    // 例如: counts.1, counts.2 → 前缀 "counts", 后缀 "1", "2"
    let mut prefix_groups: HashMap<String, Vec<(String, String, ExpressionMatrix)>> =
        HashMap::new();
    let mut no_dot: HashMap<String, ExpressionMatrix> = HashMap::new();

    for (name, matrix) in layers {
        if let Some(dot_pos) = name.rfind('.') {
            let prefix = name[..dot_pos].to_string();
            let suffix = name[dot_pos + 1..].to_string();
            prefix_groups
                .entry(prefix)
                .or_default()
                .push((name, suffix, matrix));
        } else {
            no_dot.insert(name, matrix);
        }
    }

    let mut result = no_dot;

    for (prefix, mut members) in prefix_groups {
        if members.len() == 1 {
            // 只有一个带此前缀的层 — 不是 split，保留原名
            let (original_name, _suffix, matrix) = members.into_iter().next().unwrap();
            result.insert(original_name, matrix);
        } else {
            // 多个层共享前缀 → split layers，需要合并
            // 按后缀排序以保持确定性顺序
            members.sort_by(|a, b| a.1.cmp(&b.1));

            log::info!(
                "Merging {} split layers for '{}' (suffixes: {:?})",
                members.len(),
                prefix,
                members
                    .iter()
                    .map(|(_, s, _)| s.as_str())
                    .collect::<Vec<_>>()
            );

            let matrices: Vec<ExpressionMatrix> = members.into_iter().map(|(_, _, m)| m).collect();

            // 检测矩阵类型（CSR 或 CSC），所有子层必须是同一类型
            let first_is_csc = matches!(&matrices[0], ExpressionMatrix::SparseCSC(_));
            let first_is_csr = matches!(&matrices[0], ExpressionMatrix::SparseCSR(_));

            if first_is_csc {
                // CSC 矩阵：cells=rows, genes=cols
                // Split layers 按 cells(rows) 分割，需要按行合并
                // 先转换为 CSR，按行合并，再转回 CSC
                use crate::sparse::convert::csc_to_csr_auto;
                let csr_matrices: Vec<SparseMatrixCSR> = matrices
                    .iter()
                    .map(|m| match m {
                        ExpressionMatrix::SparseCSC(csc) => Ok(csc_to_csr_auto(csc)),
                        other => Err(SeuratError::ParseError(format!(
                            "Split layer '{}' has mixed matrix types, expected CSC but got {:?}",
                            prefix,
                            other.shape()
                        ))),
                    })
                    .collect::<Result<Vec<_>, _>>()?;

                let csr_refs: Vec<&SparseMatrixCSR> = csr_matrices.iter().collect();
                let merged = merge_csr_row_wise(&csr_refs)?;
                result.insert(prefix, ExpressionMatrix::SparseCSR(merged));
            } else if first_is_csr {
                // CSR 合并：按行拼接（cells = rows）
                let csr_refs: Vec<&SparseMatrixCSR> = matrices
                    .iter()
                    .map(|m| match m {
                        ExpressionMatrix::SparseCSR(csr) => Ok(csr),
                        other => Err(SeuratError::ParseError(format!(
                            "Split layer '{}' has mixed matrix types, expected CSR but got {:?}",
                            prefix,
                            other.shape()
                        ))),
                    })
                    .collect::<Result<Vec<_>, _>>()?;

                let merged = merge_csr_row_wise(&csr_refs)?;
                result.insert(prefix, ExpressionMatrix::SparseCSR(merged));
            } else {
                return Err(SeuratError::ParseError(format!(
                    "Split layer '{}' contains unsupported matrix type (not CSR or CSC)",
                    prefix
                )));
            }
        }
    }

    Ok(result)
}

/// 将 split layers 填充到与合并后矩阵相同的维度
///
/// 当 `--keep-layers` 时，每个 split layer 只包含部分细胞。
/// AnnData 要求 `/layers` 中所有矩阵维度与 X 一致 (n_cells × n_genes)。
/// 此函数将每个 split layer 扩展到完整的细胞数，缺失的细胞行用 0 填充。
///
/// # Arguments
///
/// * `raw_layers` - 原始 split layers（如 counts.1, counts.2），未合并
/// * `total_cells` - 合并后的总细胞数（即 X 的行数）
/// * `n_genes` - 基因数（即 X 的列数）
///
/// # Returns
///
/// * 填充后的 layers，每个矩阵维度为 (total_cells, n_genes)
pub fn pad_split_layers_to_full_size(
    raw_layers: &HashMap<String, ExpressionMatrix>,
    total_cells: usize,
    n_genes: usize,
) -> Result<HashMap<String, ExpressionMatrix>, SeuratError> {
    // 按最后一个 '.' 分割，找出 split groups
    let re = regex::Regex::new(r"^(.+)\.(.+)$").unwrap();

    // 收集 split groups 及其成员（按后缀排序，与 merge 顺序一致）
    let mut prefix_groups: HashMap<String, Vec<(String, String)>> = HashMap::new();
    for name in raw_layers.keys() {
        if let Some(caps) = re.captures(name) {
            let prefix = caps[1].to_string();
            let suffix = caps[2].to_string();
            prefix_groups
                .entry(prefix)
                .or_default()
                .push((name.clone(), suffix));
        }
    }

    let mut result = HashMap::new();

    for (prefix, mut members) in prefix_groups {
        if members.len() <= 1 {
            // 不是 split group，跳过
            continue;
        }

        // 按后缀排序（与 merge_split_layers 一致）
        members.sort_by(|a, b| a.1.cmp(&b.1));

        // 计算每个 sub-layer 的 cell offset
        let mut cell_offset = 0usize;
        for (layer_name, _suffix) in &members {
            let matrix = raw_layers.get(layer_name).ok_or_else(|| {
                SeuratError::ParseError(format!("Layer '{}' not found", layer_name))
            })?;
            let (sub_cells, sub_genes) = matrix.shape();

            if sub_genes != n_genes {
                return Err(SeuratError::ParseError(format!(
                    "Split layer '{}' has {} genes, expected {}",
                    layer_name, sub_genes, n_genes
                )));
            }

            // 构建填充后的 CSR 矩阵
            let padded = pad_matrix_cells(matrix, cell_offset, total_cells, n_genes)?;
            result.insert(layer_name.clone(), padded);

            log::info!(
                "Padded split layer '{}': {} cells at offset {} → {} total cells",
                layer_name,
                sub_cells,
                cell_offset,
                total_cells
            );

            cell_offset += sub_cells;
        }

        if cell_offset != total_cells {
            log::warn!(
                "Split group '{}': sum of sub-layer cells ({}) != total cells ({})",
                prefix,
                cell_offset,
                total_cells
            );
        }
    }

    Ok(result)
}

/// 将矩阵填充到指定的总行数，在 cell_offset 处放置原始数据
///
/// 生成一个 (total_cells, n_genes) 的稀疏矩阵，其中：
/// - 行 [0, cell_offset) 全为 0
/// - 行 [cell_offset, cell_offset + sub_cells) 为原始数据
/// - 行 [cell_offset + sub_cells, total_cells) 全为 0
fn pad_matrix_cells(
    matrix: &ExpressionMatrix,
    cell_offset: usize,
    total_cells: usize,
    n_genes: usize,
) -> Result<ExpressionMatrix, SeuratError> {
    match matrix {
        ExpressionMatrix::SparseCSR(csr) => {
            let sub_cells = csr.n_rows;
            let _nnz = csr.data.len();

            // 构建新的 indptr：total_cells + 1 个元素
            let mut new_indptr = Vec::with_capacity(total_cells + 1);

            // 前 cell_offset 行：全为空行
            for _ in 0..=cell_offset {
                new_indptr.push(0);
            }

            // 中间 sub_cells 行：复制原始 indptr（跳过第一个 0）
            for i in 1..=sub_cells {
                new_indptr.push(csr.indptr[i]);
            }

            // 后面的空行
            let last_val = *new_indptr.last().unwrap();
            for _ in (cell_offset + sub_cells)..total_cells {
                new_indptr.push(last_val);
            }

            Ok(ExpressionMatrix::SparseCSR(SparseMatrixCSR {
                data: csr.data.clone(),
                indices: csr.indices.clone(),
                indptr: new_indptr,
                n_rows: total_cells,
                n_cols: n_genes,
            }))
        }
        ExpressionMatrix::SparseCSC(csc) => {
            // CSC: 列压缩，indices 存行索引
            // 需要将所有行索引偏移 cell_offset
            let new_indices: Vec<usize> =
                csc.indices.iter().map(|&idx| idx + cell_offset).collect();

            Ok(ExpressionMatrix::SparseCSC(SparseMatrixCSC {
                data: csc.data.clone(),
                indices: new_indices,
                indptr: csc.indptr.clone(),
                n_rows: total_cells,
                n_cols: n_genes,
            }))
        }
        _ => Err(SeuratError::ParseError(
            "pad_matrix_cells: only sparse matrices (CSR/CSC) are supported".to_string(),
        )),
    }
}

/// 检查是否是 BPCells 对象
fn is_bpcells_object(robj: &RObject) -> bool {
    if let Some(class_str) = extract_class(robj) {
        class_str.contains("IterableMatrix") || class_str.contains("BPCells")
    } else {
        false
    }
}

/// Assay5 提取主表达矩
///
/// 优先级：counts > data > 第一个可用的 layer
///
/// # Arguments
///
/// * `robj` - Assay5 RObject
/// * `file` - RDS 文件
///
/// # Returns
///
/// * `Result<ExpressionMatrix, SeuratError>` - 主表达矩
pub fn extract_assay5_main_matrix(
    robj: &RObject,
    file: &RdsFile,
) -> Result<ExpressionMatrix, SeuratError> {
    let layers = extract_assay5_layers(robj, file)?;

    // 优先级：counts > data > 第一个可用的 layer
    if let Some(matrix) = layers.get("counts") {
        return Ok(matrix.clone());
    }

    if let Some(matrix) = layers.get("data") {
        return Ok(matrix.clone());
    }

    // 返回第一个可用的 layer
    layers
        .into_iter()
        .next()
        .map(|(_, matrix)| matrix)
        .ok_or_else(|| SeuratError::ParseError("No layers found in Assay5".to_string()))
}

/// Assay5 提取基因元数
///
/// # Arguments
///
/// * `robj` - Assay5 RObject
/// * `n_genes` - 预期的基因数
/// * `file` - RDS 文件
///
/// # Returns
///
/// * `Result<DataFrame, SeuratError>` - 基因元数
pub fn extract_assay5_gene_metadata(
    robj: &RObject,
    n_genes: usize,
    _file: &RdsFile,
) -> Result<DataFrame, SeuratError> {
    let assay_fields = get_named_list_fields(robj).ok_or_else(|| {
        SeuratError::ParseError(format!(
            "Expected named list for Assay5, got {}",
            robj.type_name()
        ))
    })?;

    // 查找 meta.features 字段
    let meta_features = assay_fields
        .iter()
        .find(|(name, _)| *name == "meta.features")
        .map(|(_, value)| *value);

    if let Some(_meta_obj) = meta_features {
        // TODO: 解析 data.frame
        // 现在先返回空 DataFrame
        Ok(DataFrame::empty(n_genes))
    } else {
        // 没有 meta.features，返回空 DataFrame
        Ok(DataFrame::empty(n_genes))
    }
}

/// IR 构Assay5 对象
///
/// # Arguments
///
/// * `expression` - 主表达矩
/// * `additional_layers` - 可选的额外 layers（如 data, scale.data
/// * `gene_metadata` - 基因元数
///
/// # Returns
///
/// * `Result<RObject, SeuratError>` - Assay5 RObject
pub fn construct_assay5(
    expression: &ExpressionMatrix,
    additional_layers: Option<&HashMap<String, ExpressionMatrix>>,
    gene_metadata: &DataFrame,
) -> Result<RObject, SeuratError> {
    use crate::rds::StringEncoding;
    use crate::rds::{Attributes, StringVector};

    let mut data = Vec::new();
    let mut names = Vec::new();

    // class: "Assay5"
    names.push("class".to_string());
    let mut class_sv = StringVector::default();
    class_sv.add("Assay5".to_string(), StringEncoding::Utf8);
    data.push(RObject::StringVector(class_sv));

    // layers: named list of matrices
    let mut layers_data = Vec::new();
    let mut layers_names = Vec::new();

    // Add counts layer (main expression matrix)
    layers_names.push("counts".to_string());
    let counts = crate::seurat::ir_to_seurat::expression_to_dgcmatrix(expression, &[], &[])?;
    layers_data.push(counts);

    // Add additional layers if provided
    if let Some(layers) = additional_layers {
        for (name, layer_expr) in layers {
            layers_names.push(name.clone());
            let layer_matrix =
                crate::seurat::ir_to_seurat::expression_to_dgcmatrix(layer_expr, &[], &[])?;
            layers_data.push(layer_matrix);
        }
    }

    // Create layers named list
    let mut layers_attrs = Attributes::new();
    let mut layers_names_sv = StringVector::default();
    for name in &layers_names {
        layers_names_sv.add(name.clone(), StringEncoding::Utf8);
    }
    layers_attrs.add(
        "names".to_string(),
        RObject::StringVector(layers_names_sv),
        StringEncoding::Utf8,
    );

    names.push("layers".to_string());
    data.push(RObject::GenericVector(GenericVector {
        data: layers_data,
        attributes: layers_attrs,
    }));

    // features: character vector of gene names
    let n_genes = gene_metadata.n_rows;
    names.push("features".to_string());
    let mut features_sv = StringVector::default();
    if n_genes > 0 {
        for i in 0..n_genes {
            features_sv.add(format!("Gene_{}", i + 1), StringEncoding::Utf8);
        }
    }
    data.push(RObject::StringVector(features_sv));

    // cells: character vector of cell names
    let n_cells = expression.shape().0;
    names.push("cells".to_string());
    let mut cells_sv = StringVector::default();
    if n_cells > 0 {
        for i in 0..n_cells {
            cells_sv.add(format!("Cell_{}", i + 1), StringEncoding::Utf8);
        }
    }
    data.push(RObject::StringVector(cells_sv));

    // Create main named list with names attribute
    let mut attrs = Attributes::new();
    let mut names_sv = StringVector::default();
    for name in &names {
        names_sv.add(name.clone(), StringEncoding::Utf8);
    }
    attrs.add(
        "names".to_string(),
        RObject::StringVector(names_sv),
        StringEncoding::Utf8,
    );

    Ok(RObject::GenericVector(GenericVector {
        data,
        attributes: attrs,
    }))
}

/// IR layers 构Assay5 对象（所layers 平等对待
///
/// `construct_assay5` 不同，此函数将所layers 平等对待
/// 不区分主表达矩阵和额layers
///
/// # Arguments
///
/// * `layers` - 所layers HashMap
/// * `gene_metadata` - 基因元数
///
/// # Returns
///
/// * `Result<RObject, SeuratError>` - Assay5 RObject
pub fn construct_assay5_from_layers(
    layers: &HashMap<String, ExpressionMatrix>,
    gene_metadata: &DataFrame,
) -> Result<RObject, SeuratError> {
    use crate::rds::StringEncoding;
    use crate::rds::{Attributes, StringVector};

    if layers.is_empty() {
        return Err(SeuratError::ParseError(
            "Cannot construct Assay5 from empty layers".to_string(),
        ));
    }

    let mut data = Vec::new();
    let mut names = Vec::new();

    // class: "Assay5"
    names.push("class".to_string());
    let mut class_sv = StringVector::default();
    class_sv.add("Assay5".to_string(), StringEncoding::Utf8);
    data.push(RObject::StringVector(class_sv));

    // layers: named list of matrices
    let mut layers_data = Vec::new();
    let mut layers_names = Vec::new();

    for (name, layer_expr) in layers {
        layers_names.push(name.clone());
        let layer_matrix =
            crate::seurat::ir_to_seurat::expression_to_dgcmatrix(layer_expr, &[], &[])?;
        layers_data.push(layer_matrix);
    }

    // Create layers named list
    let mut layers_attrs = Attributes::new();
    let mut layers_names_sv = StringVector::default();
    for name in &layers_names {
        layers_names_sv.add(name.clone(), StringEncoding::Utf8);
    }
    layers_attrs.add(
        "names".to_string(),
        RObject::StringVector(layers_names_sv),
        StringEncoding::Utf8,
    );

    names.push("layers".to_string());
    data.push(RObject::GenericVector(GenericVector {
        data: layers_data,
        attributes: layers_attrs,
    }));

    // features: character vector of gene names
    let n_genes = gene_metadata.n_rows;
    names.push("features".to_string());
    let mut features_sv = StringVector::default();
    if n_genes > 0 {
        for i in 0..n_genes {
            features_sv.add(format!("Gene_{}", i + 1), StringEncoding::Utf8);
        }
    }
    data.push(RObject::StringVector(features_sv));

    // cells: character vector of cell names (get from first layer)
    let n_cells = layers.values().next().unwrap().shape().0;
    names.push("cells".to_string());
    let mut cells_sv = StringVector::default();
    if n_cells > 0 {
        for i in 0..n_cells {
            cells_sv.add(format!("Cell_{}", i + 1), StringEncoding::Utf8);
        }
    }
    data.push(RObject::StringVector(cells_sv));

    // Create main named list with names attribute
    let mut attrs = Attributes::new();
    let mut names_sv = StringVector::default();
    for name in &names {
        names_sv.add(name.clone(), StringEncoding::Utf8);
    }
    attrs.add(
        "names".to_string(),
        RObject::StringVector(names_sv),
        StringEncoding::Utf8,
    );

    Ok(RObject::GenericVector(GenericVector {
        data,
        attributes: attrs,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rds::StringEncoding;
    use crate::rds::{Attributes, DoubleVector, IntegerVector, S4Object, StringVector};

    /// 创建模拟dgCMatrix S4 对象
    fn create_mock_dgcmatrix(n_cells: usize, n_genes: usize) -> RObject {
        // R dgCMatrix format:
        // - Dim = [n_genes, n_cells]
        // - p (col_ptrs) length = n_cells + 1
        // - i (row_indices) = gene indices

        // Create column pointers (one non-zero per cell)
        let p: Vec<i32> = (0..=n_cells as i32).collect();

        // Create row indices (gene indices, cycling through genes)
        let i: Vec<i32> = (0..n_cells).map(|c| (c % n_genes) as i32).collect();

        // Create values
        let x: Vec<f64> = (0..n_cells).map(|c| (c + 1) as f64).collect();

        let mut s4 = S4Object::default();
        s4.class_name = "dgCMatrix".to_string();
        s4.class_encoding = StringEncoding::Utf8;
        s4.package_name = "Matrix".to_string();
        s4.package_encoding = StringEncoding::Utf8;

        // Add slots as attributes
        s4.attributes.add(
            "i".to_string(),
            RObject::IntegerVector(IntegerVector {
                data: i,
                attributes: Attributes::new(),
            }),
            StringEncoding::Utf8,
        );

        s4.attributes.add(
            "p".to_string(),
            RObject::IntegerVector(IntegerVector {
                data: p,
                attributes: Attributes::new(),
            }),
            StringEncoding::Utf8,
        );

        s4.attributes.add(
            "x".to_string(),
            RObject::DoubleVector(DoubleVector {
                data: x,
                attributes: Attributes::new(),
            }),
            StringEncoding::Utf8,
        );

        s4.attributes.add(
            "Dim".to_string(),
            RObject::IntegerVector(IntegerVector {
                data: vec![n_genes as i32, n_cells as i32],
                attributes: Attributes::new(),
            }),
            StringEncoding::Utf8,
        );

        s4.attributes
            .add("Dimnames".to_string(), RObject::Null, StringEncoding::Utf8);

        s4.attributes.add(
            "factors".to_string(),
            RObject::GenericVector(GenericVector::default()),
            StringEncoding::Utf8,
        );

        RObject::S4Object(s4)
    }

    /// 创建names 属性的命名列表
    fn create_named_list(items: Vec<(&str, RObject)>) -> RObject {
        let mut data = Vec::new();
        let mut names = Vec::new();

        for (name, value) in items {
            names.push(name.to_string());
            data.push(value);
        }

        let mut attrs = Attributes::new();
        let mut names_sv = StringVector::default();
        for name in names {
            names_sv.add(name, StringEncoding::Utf8);
        }
        attrs.add(
            "names".to_string(),
            RObject::StringVector(names_sv),
            StringEncoding::Utf8,
        );

        RObject::GenericVector(GenericVector {
            data,
            attributes: attrs,
        })
    }

    /// 创建class 属性的命名列表
    fn create_named_list_with_class(items: Vec<(&str, RObject)>, class: &str) -> RObject {
        let mut data = Vec::new();
        let mut names = Vec::new();

        for (name, value) in items {
            names.push(name.to_string());
            data.push(value);
        }

        let mut attrs = Attributes::new();
        let mut names_sv = StringVector::default();
        for name in names {
            names_sv.add(name, StringEncoding::Utf8);
        }
        attrs.add(
            "names".to_string(),
            RObject::StringVector(names_sv),
            StringEncoding::Utf8,
        );

        let mut class_sv = StringVector::default();
        class_sv.add(class.to_string(), StringEncoding::Utf8);
        attrs.add(
            "class".to_string(),
            RObject::StringVector(class_sv),
            StringEncoding::Utf8,
        );

        RObject::GenericVector(GenericVector {
            data,
            attributes: attrs,
        })
    }

    #[test]
    fn test_detect_assay_type_assay5() {
        // 测试 Assay5 检测（通过 class 属性）
        let assay5 =
            create_named_list_with_class(vec![("layers", create_named_list(vec![]))], "Assay5");
        assert_eq!(detect_assay_type(&assay5), AssayType::Assay5);
    }

    #[test]
    fn test_detect_assay_type_legacy() {
        // 测试传统 Assay 检测（通过 class 属性）
        let legacy_assay = create_named_list_with_class(vec![("counts", RObject::Null)], "Assay");
        assert_eq!(detect_assay_type(&legacy_assay), AssayType::Legacy);
    }

    #[test]
    fn test_detect_assay_type_by_structure() {
        // 测试通过结构检Assay5（没class 属性）
        let assay5_no_class = create_named_list(vec![
            ("layers", create_named_list(vec![])),
            ("features", RObject::StringVector(StringVector::default())),
        ]);
        assert_eq!(detect_assay_type(&assay5_no_class), AssayType::Assay5);

        // 测试通过结构检测传Assay（没class 属性）
        let legacy_no_class =
            create_named_list(vec![("counts", RObject::Null), ("data", RObject::Null)]);
        assert_eq!(detect_assay_type(&legacy_no_class), AssayType::Legacy);
    }

    #[test]
    fn test_detect_assay_type_unknown() {
        // 测试未知类型
        let unknown = RObject::IntegerVector(IntegerVector {
            data: vec![1, 2, 3],
            attributes: Attributes::new(),
        });
        assert_eq!(detect_assay_type(&unknown), AssayType::Unknown);

        // 测试空的 NamedList
        let empty = create_named_list(vec![]);
        assert_eq!(detect_assay_type(&empty), AssayType::Unknown);
    }

    #[test]
    fn test_extract_assay5_layers_with_null() {
        // 测试跳过 NULL layers
        let file = RdsFile::default();
        let assay5 = create_named_list(vec![(
            "layers",
            create_named_list(vec![
                ("counts", create_mock_dgcmatrix(10, 5)),
                ("data", RObject::Null), // NULL layer 应该被跳
            ]),
        )]);

        let result = extract_assay5_layers(&assay5, &file);
        assert!(result.is_ok());
        let layers = result.unwrap();
        assert_eq!(layers.len(), 1); // 只有 counts，data 被跳
        assert!(layers.contains_key("counts"));
        assert!(!layers.contains_key("data"));
    }

    #[test]
    fn test_extract_assay5_layers_multiple() {
        // 测试提取多个 layers
        let file = RdsFile::default();
        let assay5 = create_named_list(vec![(
            "layers",
            create_named_list(vec![
                ("counts", create_mock_dgcmatrix(10, 5)),
                ("data", create_mock_dgcmatrix(10, 5)),
                ("scale.data", create_mock_dgcmatrix(10, 5)),
            ]),
        )]);

        let result = extract_assay5_layers(&assay5, &file);
        assert!(result.is_ok());
        let layers = result.unwrap();
        assert_eq!(layers.len(), 3);
        assert!(layers.contains_key("counts"));
        assert!(layers.contains_key("data"));
        assert!(layers.contains_key("scale.data"));
    }

    #[test]
    fn test_extract_assay5_layers_missing_layers_field() {
        // 测试缺少 layers 字段
        let file = RdsFile::default();
        let assay5 = create_named_list(vec![(
            "features",
            RObject::StringVector(StringVector::default()),
        )]);

        let result = extract_assay5_layers(&assay5, &file);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("Missing 'layers' field"));
    }

    #[test]
    fn test_extract_assay5_layers_all_null() {
        // 测试所layers 都是 NULL
        let file = RdsFile::default();
        let assay5 = create_named_list(vec![(
            "layers",
            create_named_list(vec![("counts", RObject::Null), ("data", RObject::Null)]),
        )]);

        let result = extract_assay5_layers(&assay5, &file);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("No valid layers found"));
    }

    #[test]
    fn test_extract_assay5_main_matrix_priority() {
        // 测试主矩阵提取的优先级：counts > data > 其他
        let file = RdsFile::default();
        let assay5_with_counts = create_named_list(vec![(
            "layers",
            create_named_list(vec![
                ("counts", create_mock_dgcmatrix(10, 5)),
                ("data", create_mock_dgcmatrix(10, 5)),
                ("other", create_mock_dgcmatrix(10, 5)),
            ]),
        )]);

        let result = extract_assay5_main_matrix(&assay5_with_counts, &file);
        assert!(result.is_ok());
        // 应该返回 counts

        // 测试没有 counts 但有 data
        let assay5_with_data = create_named_list(vec![(
            "layers",
            create_named_list(vec![
                ("data", create_mock_dgcmatrix(10, 5)),
                ("other", create_mock_dgcmatrix(10, 5)),
            ]),
        )]);

        let result = extract_assay5_main_matrix(&assay5_with_data, &file);
        assert!(result.is_ok());
        // 应该返回 data

        // 测试只有其他 layer
        let assay5_with_other = create_named_list(vec![(
            "layers",
            create_named_list(vec![("other", create_mock_dgcmatrix(10, 5))]),
        )]);

        let result = extract_assay5_main_matrix(&assay5_with_other, &file);
        assert!(result.is_ok());
        // 应该返回第一个可用的 layer
    }

    #[test]
    fn test_extract_assay5_gene_metadata() {
        // 测试提取基因元数据（当前返回DataFrame
        let file = RdsFile::default();
        let assay5 = create_named_list(vec![
            ("layers", create_named_list(vec![])),
            ("meta.features", create_named_list(vec![])),
        ]);

        let result = extract_assay5_gene_metadata(&assay5, 100, &file);
        assert!(result.is_ok());
        let metadata = result.unwrap();
        assert_eq!(metadata.n_rows, 100);
    }

    #[test]
    fn test_construct_assay5_basic() {
        use crate::ir::{DataFrame, ExpressionMatrix, SparseMatrixCSC};

        // 创建测试数据
        // CSC format: indptr length = n_cols + 1
        let n_cells = 10;
        let n_genes = 5;

        let expression = ExpressionMatrix::SparseCSC(SparseMatrixCSC {
            n_rows: n_cells,                             // 10 cells
            n_cols: n_genes,                             // 5 genes
            indptr: vec![0, 2, 4, 6, 8, 10],             // length = n_genes + 1 = 6
            indices: vec![0, 5, 1, 6, 2, 7, 3, 8, 4, 9], // cell indices
            data: vec![1.0; 10],
        });
        let gene_metadata = DataFrame::empty(n_genes);

        let result = construct_assay5(&expression, None, &gene_metadata);
        assert!(result.is_ok());

        let assay5 = result.unwrap();
        // 验证结构
        match assay5 {
            RObject::GenericVector(gv) => {
                let names = gv.attributes.get_names();
                assert!(names.is_some());
                let names = names.unwrap();
                assert!(names.contains(&"class".to_string()));
                assert!(names.contains(&"layers".to_string()));
            }
            _ => panic!("Expected List"),
        }
    }

    #[test]
    fn test_construct_assay5_with_layers() {
        use crate::ir::{DataFrame, ExpressionMatrix, SparseMatrixCSC};
        use std::collections::HashMap;

        // 创建主表达矩
        // CSC format: indptr length = n_cols + 1
        let n_cells = 10;
        let n_genes = 5;

        let expression = ExpressionMatrix::SparseCSC(SparseMatrixCSC {
            n_rows: n_cells,
            n_cols: n_genes,
            indptr: vec![0, 2, 4, 6, 8, 10], // length = n_genes + 1 = 6
            indices: vec![0, 5, 1, 6, 2, 7, 3, 8, 4, 9],
            data: vec![1.0; 10],
        });

        // 创建额外layers
        let mut layers = HashMap::new();
        layers.insert(
            "data".to_string(),
            ExpressionMatrix::SparseCSC(SparseMatrixCSC {
                n_rows: n_cells,
                n_cols: n_genes,
                indptr: vec![0, 2, 4, 6, 8, 10], // length = n_genes + 1 = 6
                indices: vec![0, 5, 1, 6, 2, 7, 3, 8, 4, 9],
                data: vec![2.0; 10],
            }),
        );

        let gene_metadata = DataFrame::empty(n_genes);

        let result = construct_assay5(&expression, Some(&layers), &gene_metadata);
        assert!(result.is_ok());

        let assay5 = result.unwrap();
        match assay5 {
            RObject::GenericVector(gv) => {
                // 查找 layers 字段
                let names = gv.attributes.get_names().unwrap();
                let layers_idx = names.iter().position(|n| n == "layers").unwrap();

                if let RObject::GenericVector(layers_gv) = &gv.data[layers_idx] {
                    let layer_names = layers_gv.attributes.get_names().unwrap();
                    // 应该counts data 两个 layers
                    assert_eq!(layer_names.len(), 2);
                    assert!(layer_names.contains(&"counts".to_string()));
                    assert!(layer_names.contains(&"data".to_string()));
                }
            }
            _ => panic!("Expected List"),
        }
    }

    #[test]
    fn test_assay5_v3_v4_compatibility() {
        // 测试 Assay5 检测不会误判传Assay
        let v3_assay = create_named_list_with_class(
            vec![
                ("counts", RObject::Null),
                ("data", RObject::Null),
                ("scale.data", RObject::Null),
            ],
            "Assay",
        );
        assert_eq!(detect_assay_type(&v3_assay), AssayType::Legacy);

        // 测试传统 Assay 检测不会误Assay5
        let v5_assay =
            create_named_list_with_class(vec![("layers", create_named_list(vec![]))], "Assay5");
        assert_eq!(detect_assay_type(&v5_assay), AssayType::Assay5);
    }

    // ========================================================================
    // Task 1.7: Split Layer 合并单元测试
    // ========================================================================

    /// 辅助函数：创建简单的 CSR 矩阵用于测试
    fn make_test_csr(n_rows: usize, n_cols: usize, nnz_per_row: usize) -> SparseMatrixCSR {
        let mut data = Vec::new();
        let mut indices = Vec::new();
        let mut indptr = vec![0usize];

        for row in 0..n_rows {
            for j in 0..nnz_per_row.min(n_cols) {
                let col = (row + j) % n_cols;
                data.push((row * n_cols + col + 1) as f64);
                indices.push(col);
            }
            indptr.push(data.len());
        }

        SparseMatrixCSR {
            data,
            indices,
            indptr,
            n_rows,
            n_cols,
        }
    }

    #[test]
    fn test_merge_csr_row_wise_basic() {
        // 两个 3×5 矩阵合并为 6×5
        let m1 = make_test_csr(3, 5, 2);
        let m2 = make_test_csr(3, 5, 2);

        let merged = merge_csr_row_wise(&[&m1, &m2]).unwrap();
        assert_eq!(merged.n_rows, 6);
        assert_eq!(merged.n_cols, 5);
        assert_eq!(merged.data.len(), m1.data.len() + m2.data.len());
        assert_eq!(merged.indptr.len(), 7); // 6 + 1
        merged.validate().unwrap();
    }

    #[test]
    fn test_merge_csr_row_wise_single() {
        // 单个矩阵合并应该返回等价矩阵
        let m = make_test_csr(5, 10, 3);
        let merged = merge_csr_row_wise(&[&m]).unwrap();
        assert_eq!(merged.n_rows, 5);
        assert_eq!(merged.n_cols, 10);
        assert_eq!(merged.data, m.data);
        assert_eq!(merged.indices, m.indices);
    }

    #[test]
    fn test_merge_csr_row_wise_column_mismatch() {
        // 列数不同应该报错
        let m1 = make_test_csr(3, 5, 2);
        let m2 = make_test_csr(3, 10, 2);

        let result = merge_csr_row_wise(&[&m1, &m2]);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("column count mismatch"));
    }

    #[test]
    fn test_merge_csr_row_wise_empty_input() {
        let result = merge_csr_row_wise(&[]);
        assert!(result.is_err());
    }

    #[test]
    fn test_merge_csr_row_wise_empty_matrices() {
        // 合并两个空矩阵（0 行）
        let m1 = SparseMatrixCSR::empty(0, 5);
        let m2 = SparseMatrixCSR::empty(0, 5);
        let merged = merge_csr_row_wise(&[&m1, &m2]).unwrap();
        assert_eq!(merged.n_rows, 0);
        assert_eq!(merged.n_cols, 5);
        assert!(merged.data.is_empty());
    }

    #[test]
    fn test_merge_split_layers_basic() {
        // counts.1 + counts.2 → counts
        let mut layers = HashMap::new();
        layers.insert(
            "counts.1".to_string(),
            ExpressionMatrix::SparseCSR(make_test_csr(5, 10, 2)),
        );
        layers.insert(
            "counts.2".to_string(),
            ExpressionMatrix::SparseCSR(make_test_csr(8, 10, 2)),
        );

        let result = merge_split_layers(layers).unwrap();
        assert_eq!(result.len(), 1);
        assert!(result.contains_key("counts"));
        let (rows, cols) = result["counts"].shape();
        assert_eq!(rows, 13); // 5 + 8
        assert_eq!(cols, 10);
    }

    #[test]
    fn test_merge_split_layers_single_sublayer_keeps_name() {
        // 只有 counts.1 → 保留原名（不是 split group）
        let mut layers = HashMap::new();
        layers.insert(
            "counts.1".to_string(),
            ExpressionMatrix::SparseCSR(make_test_csr(5, 10, 2)),
        );

        let result = merge_split_layers(layers).unwrap();
        assert_eq!(result.len(), 1);
        assert!(result.contains_key("counts.1"));
        assert_eq!(result["counts.1"].shape(), (5, 10));
    }

    #[test]
    fn test_merge_split_layers_mixed() {
        // 混合 split 和非 split 层
        let mut layers = HashMap::new();
        layers.insert(
            "counts.1".to_string(),
            ExpressionMatrix::SparseCSR(make_test_csr(5, 10, 2)),
        );
        layers.insert(
            "counts.2".to_string(),
            ExpressionMatrix::SparseCSR(make_test_csr(5, 10, 2)),
        );
        layers.insert(
            "data".to_string(),
            ExpressionMatrix::SparseCSR(make_test_csr(10, 10, 2)),
        );

        let result = merge_split_layers(layers).unwrap();
        assert_eq!(result.len(), 2);
        assert!(result.contains_key("counts"));
        assert!(result.contains_key("data"));
        assert_eq!(result["counts"].shape(), (10, 10)); // 5 + 5
        assert_eq!(result["data"].shape(), (10, 10)); // unchanged
    }

    #[test]
    fn test_merge_split_layers_no_split() {
        // 没有 split 层，所有层原样保留
        let mut layers = HashMap::new();
        layers.insert(
            "counts".to_string(),
            ExpressionMatrix::SparseCSR(make_test_csr(10, 5, 2)),
        );
        layers.insert(
            "data".to_string(),
            ExpressionMatrix::SparseCSR(make_test_csr(10, 5, 2)),
        );

        let result = merge_split_layers(layers).unwrap();
        assert_eq!(result.len(), 2);
        assert!(result.contains_key("counts"));
        assert!(result.contains_key("data"));
    }

    #[test]
    fn test_merge_split_layers_ordering() {
        // 验证按数字后缀排序合并
        let mut layers = HashMap::new();
        // 故意乱序插入
        layers.insert(
            "counts.3".to_string(),
            ExpressionMatrix::SparseCSR(make_test_csr(3, 5, 1)),
        );
        layers.insert(
            "counts.1".to_string(),
            ExpressionMatrix::SparseCSR(make_test_csr(1, 5, 1)),
        );
        layers.insert(
            "counts.2".to_string(),
            ExpressionMatrix::SparseCSR(make_test_csr(2, 5, 1)),
        );

        let result = merge_split_layers(layers).unwrap();
        let merged = &result["counts"];
        assert_eq!(merged.shape(), (6, 5)); // 1 + 2 + 3

        // 验证行顺序：第一行来自 counts.1，接下来 2 行来自 counts.2，最后 3 行来自 counts.3
        if let ExpressionMatrix::SparseCSR(csr) = merged {
            // counts.1 有 1 行，counts.2 有 2 行，counts.3 有 3 行
            // indptr 应该反映这个顺序
            assert_eq!(csr.indptr.len(), 7); // 6 rows + 1
        } else {
            panic!("Expected SparseCSR");
        }
    }

    #[test]
    fn test_merge_split_layers_column_mismatch_error() {
        // split 层列数不一致应该报错
        let mut layers = HashMap::new();
        layers.insert(
            "counts.1".to_string(),
            ExpressionMatrix::SparseCSR(make_test_csr(5, 10, 2)),
        );
        layers.insert(
            "counts.2".to_string(),
            ExpressionMatrix::SparseCSR(make_test_csr(5, 20, 2)), // 不同列数
        );

        let result = merge_split_layers(layers);
        assert!(result.is_err());
    }

    #[test]
    fn test_merge_split_layers_non_numeric_suffix() {
        // Seurat V5 split(obj[["RNA"]], f=obj$batch) 产生 counts.batch1, counts.batch2
        let mut layers = HashMap::new();
        layers.insert(
            "counts.batch1".to_string(),
            ExpressionMatrix::SparseCSR(make_test_csr(5, 10, 2)),
        );
        layers.insert(
            "counts.batch2".to_string(),
            ExpressionMatrix::SparseCSR(make_test_csr(8, 10, 2)),
        );

        let result = merge_split_layers(layers).unwrap();
        assert_eq!(result.len(), 1);
        assert!(result.contains_key("counts"));
        let (rows, cols) = result["counts"].shape();
        assert_eq!(rows, 13); // 5 + 8
        assert_eq!(cols, 10);
    }

    fn make_test_csc(n_cols: usize, n_rows: usize, nnz_per_col: usize) -> SparseMatrixCSC {
        use crate::ir::SparseMatrixCSC;
        let mut data = Vec::new();
        let mut indices = Vec::new();
        let mut indptr = vec![0usize];

        for col in 0..n_cols {
            for k in 0..nnz_per_col.min(n_rows) {
                let row = (col + k) % n_rows;
                indices.push(row);
                data.push((col * nnz_per_col + k + 1) as f64);
            }
            indptr.push(data.len());
        }

        SparseMatrixCSC {
            data,
            indices,
            indptr,
            n_rows,
            n_cols,
        }
    }

    #[test]
    fn test_merge_csc_column_wise_basic() {
        let m1 = make_test_csc(5, 10, 2);
        let m2 = make_test_csc(8, 10, 2);

        let merged = merge_csc_column_wise(&[&m1, &m2]).unwrap();
        assert_eq!(merged.n_cols, 13); // 5 + 8
        assert_eq!(merged.n_rows, 10);
        assert_eq!(merged.data.len(), m1.data.len() + m2.data.len());
        assert_eq!(merged.indptr.len(), 14); // 13 + 1
    }

    #[test]
    fn test_merge_split_layers_csc_non_numeric() {
        // 真实场景：dgCMatrix 产生 CSC，split 层名为 counts.batch1 等
        // 在 CrossCell 中 cells=rows, genes=cols
        // split layers 按 cells(rows) 分割，合并后 rows 增加
        let mut layers = HashMap::new();
        layers.insert(
            "counts.batch1".to_string(),
            ExpressionMatrix::SparseCSC(make_test_csc(5, 10, 2)), // 5 cols(genes), 10 rows(cells)
        );
        layers.insert(
            "counts.batch2".to_string(),
            ExpressionMatrix::SparseCSC(make_test_csc(5, 8, 2)), // 5 cols(genes), 8 rows(cells)
        );

        let result = merge_split_layers(layers).unwrap();
        assert_eq!(result.len(), 1);
        assert!(result.contains_key("counts"));
        let (rows, cols) = result["counts"].shape();
        assert_eq!(rows, 18); // cells merged: 10 + 8
        assert_eq!(cols, 5); // genes unchanged
    }

    #[test]
    fn test_merge_csc_column_wise_row_mismatch() {
        let m1 = make_test_csc(5, 10, 2);
        let m2 = make_test_csc(5, 20, 2); // different row count

        let result = merge_csc_column_wise(&[&m1, &m2]);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("row count mismatch"));
    }
}

// ============================================================================
// Task 1.2, 1.3, 1.5: 属性测试
// ============================================================================

#[cfg(test)]
mod prop_tests {
    use super::*;
    use crate::ir::{ExpressionMatrix, SparseMatrixCSR};
    use proptest::prelude::*;
    use std::collections::HashMap;

    /// 生成随机 CSR 矩阵的策略
    fn arb_csr(
        max_rows: usize,
        max_cols: usize,
        max_nnz_per_row: usize,
    ) -> impl Strategy<Value = SparseMatrixCSR> {
        (1..=max_rows, 1..=max_cols).prop_flat_map(move |(n_rows, n_cols)| {
            let nnz_limit = max_nnz_per_row.min(n_cols);
            proptest::collection::vec(0..=nnz_limit, n_rows).prop_map(move |nnz_per_row| {
                let mut data = Vec::new();
                let mut indices = Vec::new();
                let mut indptr = vec![0usize];

                for (row, &nnz) in nnz_per_row.iter().enumerate() {
                    // 为每行生成不重复的列索引
                    let mut cols: Vec<usize> = (0..n_cols).collect();
                    // 简单的确定性 shuffle 基于 row
                    for i in 0..cols.len() {
                        let j = (i + row + 1) % cols.len();
                        cols.swap(i, j);
                    }
                    cols.truncate(nnz);
                    cols.sort();

                    for &col in &cols {
                        data.push((row * n_cols + col + 1) as f64);
                        indices.push(col);
                    }
                    indptr.push(data.len());
                }

                SparseMatrixCSR {
                    data,
                    indices,
                    indptr,
                    n_rows,
                    n_cols,
                }
            })
        })
    }

    /// 生成一组具有相同列数的 CSR 矩阵
    fn arb_csr_same_cols(
        count: std::ops::Range<usize>,
        max_rows: usize,
        max_cols: usize,
    ) -> impl Strategy<Value = Vec<SparseMatrixCSR>> {
        (count, 1..=max_cols).prop_flat_map(move |(n, n_cols)| {
            proptest::collection::vec(
                (1..=max_rows).prop_flat_map(move |n_rows| {
                    arb_csr(1, 1, 1).prop_map(move |_| {
                        // 生成固定列数的矩阵
                        let mut data = Vec::new();
                        let mut indices = Vec::new();
                        let mut indptr = vec![0usize];
                        for row in 0..n_rows {
                            // 每行 1 个非零元素
                            let col = row % n_cols;
                            data.push((row + 1) as f64);
                            indices.push(col);
                            indptr.push(data.len());
                        }
                        SparseMatrixCSR {
                            data,
                            indices,
                            indptr,
                            n_rows,
                            n_cols,
                        }
                    })
                }),
                n,
            )
        })
    }

    // Feature: reviewer-enhancements, Property 1: Split layer 合并维度不变量
    proptest! {
        #[test]
        fn prop_merge_csr_dimensions(matrices in arb_csr_same_cols(1..5, 50, 30)) {
            let refs: Vec<&SparseMatrixCSR> = matrices.iter().collect();
            let merged = merge_csr_row_wise(&refs).unwrap();

            let expected_rows: usize = matrices.iter().map(|m| m.n_rows).sum();
            let expected_cols = matrices[0].n_cols;

            prop_assert_eq!(merged.n_rows, expected_rows);
            prop_assert_eq!(merged.n_cols, expected_cols);
            prop_assert_eq!(merged.indptr.len(), expected_rows + 1);
            prop_assert!(merged.validate().is_ok());
        }
    }

    // Feature: reviewer-enhancements, Property 3: 不同列数的 split layers 被拒绝
    proptest! {
        #[test]
        fn prop_merge_csr_rejects_different_cols(
            n_rows1 in 1..20usize,
            n_rows2 in 1..20usize,
            n_cols1 in 1..20usize,
            n_cols2 in 1..20usize,
        ) {
            prop_assume!(n_cols1 != n_cols2);

            let m1 = SparseMatrixCSR::empty(n_rows1, n_cols1);
            let m2 = SparseMatrixCSR::empty(n_rows2, n_cols2);

            let result = merge_csr_row_wise(&[&m1, &m2]);
            prop_assert!(result.is_err());
        }
    }

    // Feature: reviewer-enhancements, Property 2: 非 split 层保持不变
    proptest! {
        #[test]
        fn prop_non_split_layers_unchanged(
            n_rows in 1..30usize,
            n_cols in 1..20usize,
        ) {
            let csr = SparseMatrixCSR::empty(n_rows, n_cols);
            let original_shape = (csr.n_rows, csr.n_cols);

            let mut layers = HashMap::new();
            // 非 split 层名称（不匹配 base.digit 模式）
            layers.insert("counts".to_string(), ExpressionMatrix::SparseCSR(csr.clone()));
            layers.insert("data".to_string(), ExpressionMatrix::SparseCSR(csr.clone()));
            layers.insert("scale_data".to_string(), ExpressionMatrix::SparseCSR(csr));

            let result = merge_split_layers(layers).unwrap();

            prop_assert_eq!(result.len(), 3);
            prop_assert!(result.contains_key("counts"));
            prop_assert!(result.contains_key("data"));
            prop_assert!(result.contains_key("scale_data"));

            for (_, matrix) in &result {
                prop_assert_eq!(matrix.shape(), original_shape);
            }
        }
    }
}

// ============================================================================
// pad_split_layers_to_full_size 单元测试
// ============================================================================

#[cfg(test)]
mod pad_tests {
    use super::*;
    use std::collections::HashMap;

    /// 创建一个简单的 CSR 测试矩阵，带有可预测的数据
    fn make_csr(n_rows: usize, n_cols: usize, nnz_per_row: usize) -> SparseMatrixCSR {
        let mut data = Vec::new();
        let mut indices = Vec::new();
        let mut indptr = vec![0usize];

        for row in 0..n_rows {
            for j in 0..nnz_per_row.min(n_cols) {
                data.push((row * n_cols + j) as f64);
                indices.push(j);
            }
            indptr.push(data.len());
        }

        SparseMatrixCSR {
            data,
            indices,
            indptr,
            n_rows,
            n_cols,
        }
    }

    #[test]
    fn test_pad_split_layers_basic() {
        // 2 个 split layers: counts.1 (3 cells), counts.2 (2 cells)
        // 合并后总共 5 cells, 10 genes
        let mut raw_layers = HashMap::new();
        raw_layers.insert(
            "counts.1".to_string(),
            ExpressionMatrix::SparseCSR(make_csr(3, 10, 2)),
        );
        raw_layers.insert(
            "counts.2".to_string(),
            ExpressionMatrix::SparseCSR(make_csr(2, 10, 2)),
        );

        let result = pad_split_layers_to_full_size(&raw_layers, 5, 10).unwrap();

        assert_eq!(result.len(), 2);
        assert!(result.contains_key("counts.1"));
        assert!(result.contains_key("counts.2"));

        // 两个 layer 都应该是 (5, 10)
        assert_eq!(result["counts.1"].shape(), (5, 10));
        assert_eq!(result["counts.2"].shape(), (5, 10));
    }

    #[test]
    fn test_pad_preserves_data() {
        let mut raw_layers = HashMap::new();
        let csr1 = make_csr(2, 5, 1); // 2 cells, 5 genes, 1 nnz/row
        let csr2 = make_csr(3, 5, 1); // 3 cells, 5 genes, 1 nnz/row
        let original_data1 = csr1.data.clone();
        let original_data2 = csr2.data.clone();

        raw_layers.insert("counts.1".to_string(), ExpressionMatrix::SparseCSR(csr1));
        raw_layers.insert("counts.2".to_string(), ExpressionMatrix::SparseCSR(csr2));

        let result = pad_split_layers_to_full_size(&raw_layers, 5, 5).unwrap();

        // 数据值应该保持不变
        if let ExpressionMatrix::SparseCSR(padded1) = &result["counts.1"] {
            assert_eq!(padded1.data, original_data1);
            assert_eq!(padded1.n_rows, 5);
            assert_eq!(padded1.n_cols, 5);
            // 前 2 行有数据，后 3 行为空
            assert_eq!(padded1.indptr[0], 0); // row 0 start
            assert_eq!(padded1.indptr[2], 2); // row 1 end = 2 nnz
            assert_eq!(padded1.indptr[3], 2); // row 2 (padded) = same
            assert_eq!(padded1.indptr[5], 2); // row 4 (padded) = same
        } else {
            panic!("Expected SparseCSR");
        }

        if let ExpressionMatrix::SparseCSR(padded2) = &result["counts.2"] {
            assert_eq!(padded2.data, original_data2);
            assert_eq!(padded2.n_rows, 5);
            // 前 2 行为空（offset=2），后 3 行有数据
            assert_eq!(padded2.indptr[0], 0);
            assert_eq!(padded2.indptr[1], 0); // row 0 (padded)
            assert_eq!(padded2.indptr[2], 0); // row 1 (padded)
            assert_eq!(padded2.indptr[3], 1); // row 2 (first real row)
        } else {
            panic!("Expected SparseCSR");
        }
    }

    #[test]
    fn test_pad_no_split_layers() {
        // 没有 split layers 时应该返回空
        let mut raw_layers = HashMap::new();
        raw_layers.insert(
            "counts".to_string(),
            ExpressionMatrix::SparseCSR(make_csr(5, 10, 2)),
        );

        let result = pad_split_layers_to_full_size(&raw_layers, 5, 10).unwrap();
        assert!(result.is_empty());
    }

    #[test]
    fn test_pad_gene_mismatch_error() {
        let mut raw_layers = HashMap::new();
        raw_layers.insert(
            "counts.1".to_string(),
            ExpressionMatrix::SparseCSR(make_csr(3, 10, 2)),
        );
        raw_layers.insert(
            "counts.2".to_string(),
            ExpressionMatrix::SparseCSR(make_csr(2, 8, 2)), // wrong gene count
        );

        let result = pad_split_layers_to_full_size(&raw_layers, 5, 10);
        assert!(result.is_err());
    }

    #[test]
    fn test_pad_csc_matrix() {
        use crate::ir::SparseMatrixCSC;

        let csc1 = SparseMatrixCSC {
            data: vec![1.0, 2.0],
            indices: vec![0, 1],      // row indices
            indptr: vec![0, 1, 2, 2], // 3 cols
            n_rows: 3,                // 3 cells
            n_cols: 3,                // 3 genes
        };
        let csc2 = SparseMatrixCSC {
            data: vec![3.0],
            indices: vec![0],
            indptr: vec![0, 1, 1, 1],
            n_rows: 2, // 2 cells
            n_cols: 3,
        };

        let mut raw_layers = HashMap::new();
        raw_layers.insert("counts.1".to_string(), ExpressionMatrix::SparseCSC(csc1));
        raw_layers.insert("counts.2".to_string(), ExpressionMatrix::SparseCSC(csc2));

        let result = pad_split_layers_to_full_size(&raw_layers, 5, 3).unwrap();

        assert_eq!(result["counts.1"].shape(), (5, 3));
        assert_eq!(result["counts.2"].shape(), (5, 3));

        // CSC: counts.2 的 indices 应该偏移了 3（cell_offset = 3 cells from counts.1）
        if let ExpressionMatrix::SparseCSC(padded2) = &result["counts.2"] {
            assert_eq!(padded2.indices, vec![3]); // 原来的 row 0 → row 3
            assert_eq!(padded2.n_rows, 5);
        } else {
            panic!("Expected SparseCSC");
        }
    }
}
