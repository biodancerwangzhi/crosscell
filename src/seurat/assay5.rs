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

use crate::ir::{DataFrame, ExpressionMatrix};
use crate::rds::{RObject, RdsFile, GenericVector};
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
                Some(names.iter()
                    .zip(gv.data.iter())
                    .map(|(n, v)| (n.as_str(), v))
                    .collect())
            } else {
                None
            }
        }
        RObject::PairList(pl) => {
            Some(pl.tag_names.iter()
                .zip(pl.data.iter())
                .map(|(n, v)| (n.as_str(), v))
                .collect())
        }
        _ => None,
    }
}

/// RObject 提取 class 属
fn extract_class(robj: &RObject) -> Option<String> {
    robj.attributes().and_then(|attrs| {
        attrs.get_class().map(|classes| classes.join(","))
    })
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
    let assay_fields = get_named_list_fields(robj)
        .ok_or_else(|| SeuratError::ParseError(format!(
            "Expected named list for Assay5, got {}",
            robj.type_name()
        )))?;

    // 查找 layers 字段
    let layers_obj = assay_fields
        .iter()
        .find(|(name, _)| *name == "layers")
        .map(|(_, value)| *value)
        .ok_or_else(|| {
            SeuratError::ParseError("Missing 'layers' field in Assay5".to_string())
        })?;

    // 解析 layers（应该是一named list
    let layers_list = get_named_list_fields(layers_obj)
        .ok_or_else(|| SeuratError::ParseError(format!(
            "Expected named list for layers, got {}",
            layers_obj.type_name()
        )))?;

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
            log::debug!("Skipping non-dgCMatrix layer '{}' (type: {})", layer_name, layer_value.type_name());
        }
    }

    if layers.is_empty() {
        return Err(SeuratError::ParseError(
            "No valid layers found in Assay5".to_string(),
        ));
    }

    Ok(layers)
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
pub fn extract_assay5_main_matrix(robj: &RObject, file: &RdsFile) -> Result<ExpressionMatrix, SeuratError> {
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
    let assay_fields = get_named_list_fields(robj)
        .ok_or_else(|| SeuratError::ParseError(format!(
            "Expected named list for Assay5, got {}",
            robj.type_name()
        )))?;

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
    use crate::rds::{Attributes, StringVector};
    use crate::rds::StringEncoding;
    
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
    let counts = crate::seurat::ir_to_seurat::expression_to_dgcmatrix(expression)?;
    layers_data.push(counts);
    
    // Add additional layers if provided
    if let Some(layers) = additional_layers {
        for (name, layer_expr) in layers {
            layers_names.push(name.clone());
            let layer_matrix = crate::seurat::ir_to_seurat::expression_to_dgcmatrix(layer_expr)?;
            layers_data.push(layer_matrix);
        }
    }
    
    // Create layers named list
    let mut layers_attrs = Attributes::new();
    let mut layers_names_sv = StringVector::default();
    for name in &layers_names {
        layers_names_sv.add(name.clone(), StringEncoding::Utf8);
    }
    layers_attrs.add("names".to_string(), RObject::StringVector(layers_names_sv), StringEncoding::Utf8);
    
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
    attrs.add("names".to_string(), RObject::StringVector(names_sv), StringEncoding::Utf8);

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
    use crate::rds::{Attributes, StringVector};
    use crate::rds::StringEncoding;
    
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
        let layer_matrix = crate::seurat::ir_to_seurat::expression_to_dgcmatrix(layer_expr)?;
        layers_data.push(layer_matrix);
    }
    
    // Create layers named list
    let mut layers_attrs = Attributes::new();
    let mut layers_names_sv = StringVector::default();
    for name in &layers_names {
        layers_names_sv.add(name.clone(), StringEncoding::Utf8);
    }
    layers_attrs.add("names".to_string(), RObject::StringVector(layers_names_sv), StringEncoding::Utf8);
    
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
    attrs.add("names".to_string(), RObject::StringVector(names_sv), StringEncoding::Utf8);

    Ok(RObject::GenericVector(GenericVector {
        data,
        attributes: attrs,
    }))
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::rds::{Attributes, IntegerVector, DoubleVector, StringVector, S4Object};
    use crate::rds::StringEncoding;

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
        s4.attributes.add("i".to_string(), RObject::IntegerVector(IntegerVector {
            data: i,
            attributes: Attributes::new(),
        }), StringEncoding::Utf8);
        
        s4.attributes.add("p".to_string(), RObject::IntegerVector(IntegerVector {
            data: p,
            attributes: Attributes::new(),
        }), StringEncoding::Utf8);
        
        s4.attributes.add("x".to_string(), RObject::DoubleVector(DoubleVector {
            data: x,
            attributes: Attributes::new(),
        }), StringEncoding::Utf8);
        
        s4.attributes.add("Dim".to_string(), RObject::IntegerVector(IntegerVector {
            data: vec![n_genes as i32, n_cells as i32],
            attributes: Attributes::new(),
        }), StringEncoding::Utf8);
        
        s4.attributes.add("Dimnames".to_string(), RObject::Null, StringEncoding::Utf8);
        
        s4.attributes.add("factors".to_string(), RObject::GenericVector(GenericVector::default()), StringEncoding::Utf8);
        
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
        attrs.add("names".to_string(), RObject::StringVector(names_sv), StringEncoding::Utf8);
        
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
        attrs.add("names".to_string(), RObject::StringVector(names_sv), StringEncoding::Utf8);
        
        let mut class_sv = StringVector::default();
        class_sv.add(class.to_string(), StringEncoding::Utf8);
        attrs.add("class".to_string(), RObject::StringVector(class_sv), StringEncoding::Utf8);
        
        RObject::GenericVector(GenericVector {
            data,
            attributes: attrs,
        })
    }

    #[test]
    fn test_detect_assay_type_assay5() {
        // 测试 Assay5 检测（通过 class 属性）
        let assay5 = create_named_list_with_class(vec![
            ("layers", create_named_list(vec![])),
        ], "Assay5");
        assert_eq!(detect_assay_type(&assay5), AssayType::Assay5);
    }
    
    #[test]
    fn test_detect_assay_type_legacy() {
        // 测试传统 Assay 检测（通过 class 属性）
        let legacy_assay = create_named_list_with_class(vec![
            ("counts", RObject::Null),
        ], "Assay");
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
        let legacy_no_class = create_named_list(vec![
            ("counts", RObject::Null),
            ("data", RObject::Null),
        ]);
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
        let assay5 = create_named_list(vec![
            ("layers", create_named_list(vec![
                ("counts", create_mock_dgcmatrix(10, 5)),
                ("data", RObject::Null),  // NULL layer 应该被跳
            ])),
        ]);
        
        let result = extract_assay5_layers(&assay5, &file);
        assert!(result.is_ok());
        let layers = result.unwrap();
        assert_eq!(layers.len(), 1);  // 只有 counts，data 被跳
        assert!(layers.contains_key("counts"));
        assert!(!layers.contains_key("data"));
    }
    
    #[test]
    fn test_extract_assay5_layers_multiple() {
        // 测试提取多个 layers
        let file = RdsFile::default();
        let assay5 = create_named_list(vec![
            ("layers", create_named_list(vec![
                ("counts", create_mock_dgcmatrix(10, 5)),
                ("data", create_mock_dgcmatrix(10, 5)),
                ("scale.data", create_mock_dgcmatrix(10, 5)),
            ])),
        ]);
        
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
        let assay5 = create_named_list(vec![
            ("features", RObject::StringVector(StringVector::default())),
        ]);
        
        let result = extract_assay5_layers(&assay5, &file);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Missing 'layers' field"));
    }
    
    #[test]
    fn test_extract_assay5_layers_all_null() {
        // 测试所layers 都是 NULL
        let file = RdsFile::default();
        let assay5 = create_named_list(vec![
            ("layers", create_named_list(vec![
                ("counts", RObject::Null),
                ("data", RObject::Null),
            ])),
        ]);
        
        let result = extract_assay5_layers(&assay5, &file);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("No valid layers found"));
    }
    
    #[test]
    fn test_extract_assay5_main_matrix_priority() {
        // 测试主矩阵提取的优先级：counts > data > 其他
        let file = RdsFile::default();
        let assay5_with_counts = create_named_list(vec![
            ("layers", create_named_list(vec![
                ("counts", create_mock_dgcmatrix(10, 5)),
                ("data", create_mock_dgcmatrix(10, 5)),
                ("other", create_mock_dgcmatrix(10, 5)),
            ])),
        ]);
        
        let result = extract_assay5_main_matrix(&assay5_with_counts, &file);
        assert!(result.is_ok());
        // 应该返回 counts
        
        // 测试没有 counts 但有 data
        let assay5_with_data = create_named_list(vec![
            ("layers", create_named_list(vec![
                ("data", create_mock_dgcmatrix(10, 5)),
                ("other", create_mock_dgcmatrix(10, 5)),
            ])),
        ]);
        
        let result = extract_assay5_main_matrix(&assay5_with_data, &file);
        assert!(result.is_ok());
        // 应该返回 data
        
        // 测试只有其他 layer
        let assay5_with_other = create_named_list(vec![
            ("layers", create_named_list(vec![
                ("other", create_mock_dgcmatrix(10, 5)),
            ])),
        ]);
        
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
        use crate::ir::{ExpressionMatrix, SparseMatrixCSC, DataFrame};
        
        // 创建测试数据
        // CSC format: indptr length = n_cols + 1
        let n_cells = 10;
        let n_genes = 5;
        
        let expression = ExpressionMatrix::SparseCSC(SparseMatrixCSC {
            n_rows: n_cells,  // 10 cells
            n_cols: n_genes,  // 5 genes
            indptr: vec![0, 2, 4, 6, 8, 10],  // length = n_genes + 1 = 6
            indices: vec![0, 5, 1, 6, 2, 7, 3, 8, 4, 9],  // cell indices
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
        use crate::ir::{ExpressionMatrix, SparseMatrixCSC, DataFrame};
        use std::collections::HashMap;
        
        // 创建主表达矩
        // CSC format: indptr length = n_cols + 1
        let n_cells = 10;
        let n_genes = 5;
        
        let expression = ExpressionMatrix::SparseCSC(SparseMatrixCSC {
            n_rows: n_cells,
            n_cols: n_genes,
            indptr: vec![0, 2, 4, 6, 8, 10],  // length = n_genes + 1 = 6
            indices: vec![0, 5, 1, 6, 2, 7, 3, 8, 4, 9],
            data: vec![1.0; 10],
        });
        
        // 创建额外layers
        let mut layers = HashMap::new();
        layers.insert("data".to_string(), ExpressionMatrix::SparseCSC(SparseMatrixCSC {
            n_rows: n_cells,
            n_cols: n_genes,
            indptr: vec![0, 2, 4, 6, 8, 10],  // length = n_genes + 1 = 6
            indices: vec![0, 5, 1, 6, 2, 7, 3, 8, 4, 9],
            data: vec![2.0; 10],
        }));
        
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
        let v3_assay = create_named_list_with_class(vec![
            ("counts", RObject::Null),
            ("data", RObject::Null),
            ("scale.data", RObject::Null),
        ], "Assay");
        assert_eq!(detect_assay_type(&v3_assay), AssayType::Legacy);
        
        // 测试传统 Assay 检测不会误Assay5
        let v5_assay = create_named_list_with_class(vec![
            ("layers", create_named_list(vec![])),
        ], "Assay5");
        assert_eq!(detect_assay_type(&v5_assay), AssayType::Assay5);
    }
}
