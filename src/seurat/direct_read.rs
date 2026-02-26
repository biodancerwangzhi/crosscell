//! Direct read module for Seurat objects
//!
//! Uses the new `crate::rds` parser for correct S4 object parsing.

use crate::ir::{DataFrame, DatasetMetadata, Embedding, ExpressionMatrix, SingleCellData, SpatialData};
use crate::rds::parse::rds::{parse_rds_file, ParseRdsOptions};
use crate::rds::{RObject, RdsFile};
use crate::seurat::error::SeuratError;
use std::collections::HashMap;
use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SeuratVersion { V3, V4, V5, Unknown }

impl std::fmt::Display for SeuratVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SeuratVersion::V3 => write!(f, "V3"),
            SeuratVersion::V4 => write!(f, "V4"),
            SeuratVersion::V5 => write!(f, "V5"),
            SeuratVersion::Unknown => write!(f, "Unknown"),
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct SkippedComponents {
    pub environments: Vec<String>,
    pub closures: Vec<String>,
    pub external_pointers: Vec<String>,
    pub bytecodes: Vec<String>,
    pub languages: Vec<String>,
}

impl SkippedComponents {
    pub fn new() -> Self { Self::default() }
    pub fn has_skipped(&self) -> bool {
        !self.environments.is_empty() || !self.closures.is_empty() ||
        !self.external_pointers.is_empty() || !self.bytecodes.is_empty() || !self.languages.is_empty()
    }
    pub fn total_skipped(&self) -> usize {
        self.environments.len() + self.closures.len() + self.external_pointers.len() +
        self.bytecodes.len() + self.languages.len()
    }
}

#[derive(Debug)]
pub struct DirectReadResult {
    pub data: SingleCellData,
    pub version: SeuratVersion,
    pub skipped: SkippedComponents,
}

// ============================================================
// Public API
// ============================================================

pub fn read_seurat_direct<P: AsRef<Path>>(path: P, keep_layers: bool) -> Result<DirectReadResult, SeuratError> {
    let path = path.as_ref();
    let opts = ParseRdsOptions::default();
    let file = parse_rds_file(path, &opts)
        .map_err(|e| SeuratError::ParseError(format!("Failed to read RDS: {}", e)))?;
    let version = detect_seurat_version(&file.object)?;
    let mut skipped = SkippedComponents::new();
    let data = match version {
        SeuratVersion::V3 | SeuratVersion::V4 => extract_seurat_v3_v4(&file.object, &file, &mut skipped)?,
        SeuratVersion::V5 | SeuratVersion::Unknown => extract_seurat_v5(&file.object, &file, &mut skipped, keep_layers)?,
    };
    Ok(DirectReadResult { data, version, skipped })
}

// ============================================================
// S4 / slot helpers (using new rds types)
// ============================================================

fn extract_slot<'a>(obj: &'a RObject, slot: &str) -> Result<&'a RObject, SeuratError> {
    match obj {
        RObject::S4Object(s4) => s4.attributes.get(slot)
            .ok_or_else(|| SeuratError::MissingSlot(slot.to_string())),
        _ => Err(SeuratError::NotS4Object),
    }
}

fn extract_slot_opt<'a>(obj: &'a RObject, slot: &str) -> Result<Option<&'a RObject>, SeuratError> {
    match obj {
        RObject::S4Object(s4) => Ok(s4.attributes.get(slot)),
        _ => Err(SeuratError::NotS4Object),
    }
}

fn get_named_list_items<'a>(obj: &'a RObject) -> Option<Vec<(&'a str, &'a RObject)>> {
    match obj {
        RObject::GenericVector(gv) => {
            let names: Vec<&str> = gv.attributes.get("names").and_then(|n| {
                if let RObject::StringVector(sv) = n { Some(sv.data.iter().map(|s| s.as_str()).collect()) } else { None }
            })?;
            Some(names.into_iter().zip(gv.data.iter()).collect())
        }
        RObject::PairList(pl) => {
            Some(pl.tag_names.iter().zip(pl.data.iter()).map(|(n, v)| (n.as_str(), v)).collect())
        }
        RObject::S4Object(s4) => {
            if let Some(data_slot) = s4.attributes.get(".Data") {
                return get_named_list_items(data_slot);
            }
            if let Some(RObject::StringVector(sv)) = s4.attributes.get("names") {
                let names = &sv.data;
                for attr_name in &s4.attributes.names {
                    if let Some(RObject::GenericVector(gv)) = s4.attributes.get(attr_name) {
                        if gv.data.len() == names.len() {
                            return Some(names.iter().zip(gv.data.iter()).map(|(n, v)| (n.as_str(), v)).collect());
                        }
                    }
                }
            }
            None
        }
        _ => None,
    }
}

fn extract_string_scalar(obj: &RObject) -> Result<String, SeuratError> {
    match obj {
        RObject::StringVector(sv) if !sv.data.is_empty() => Ok(sv.data[0].clone()),
        _ => Err(SeuratError::ParseError(format!("Expected string, got {:?}", obj.sexp_type()))),
    }
}

fn extract_integer_vector(obj: &RObject) -> Result<Vec<i32>, SeuratError> {
    match obj {
        RObject::IntegerVector(iv) => Ok(iv.data.clone()),
        _ => Err(SeuratError::ParseError(format!("Expected integer vector, got {:?}", obj.sexp_type()))),
    }
}

fn extract_real_vector(obj: &RObject) -> Result<Vec<f64>, SeuratError> {
    match obj {
        RObject::DoubleVector(dv) => Ok(dv.data.clone()),
        _ => Err(SeuratError::ParseError(format!("Expected real vector, got {:?}", obj.sexp_type()))),
    }
}

fn should_skip_object(robj: &RObject, name: &str, skipped: &mut SkippedComponents) -> bool {
    match robj {
        RObject::EnvironmentIndex { .. } => { skipped.environments.push(name.to_string()); true }
        RObject::ExternalPointerIndex(_) => { skipped.external_pointers.push(name.to_string()); true }
        RObject::LanguageObject(_) => { skipped.languages.push(name.to_string()); true }
        RObject::BuiltInFunction(_) => { skipped.closures.push(name.to_string()); true }
        _ => false,
    }
}

// ============================================================
// Version detection
// ============================================================

fn detect_seurat_version(robj: &RObject) -> Result<SeuratVersion, SeuratError> {
    let s4 = match robj {
        RObject::S4Object(s4) => s4,
        _ => return Err(SeuratError::ParseError(format!("Expected S4 for Seurat, got {:?}", robj.sexp_type()))),
    };
    if s4.class_name != "Seurat" {
        return Err(SeuratError::ParseError(format!("Expected Seurat class, got {}", s4.class_name)));
    }
    if let Some(assays) = s4.attributes.get("assays") {
        if let Some(items) = get_named_list_items(assays) {
            if let Some((_, first_assay)) = items.first() {
                return Ok(match detect_assay_type(first_assay) {
                    AssayType::Assay5 => SeuratVersion::V5,
                    AssayType::Legacy => if s4.attributes.get("images").is_some() { SeuratVersion::V4 } else { SeuratVersion::V3 },
                    AssayType::Unknown => SeuratVersion::Unknown,
                });
            }
        }
    }
    Ok(SeuratVersion::Unknown)
}

// ============================================================
// Assay type detection (inline, using new rds types)
// ============================================================

#[derive(Debug, Clone, Copy, PartialEq)]
enum AssayType { Legacy, Assay5, Unknown }

fn detect_assay_type(robj: &RObject) -> AssayType {
    // Check S4 class name
    if let RObject::S4Object(s4) = robj {
        if s4.class_name.contains("Assay5") { return AssayType::Assay5; }
        if s4.class_name.contains("Assay") { return AssayType::Legacy; }
    }
    // Check class attribute on generic vectors
    if let Some(attrs) = robj_attributes(robj) {
        if let Some(RObject::StringVector(sv)) = attrs.get("class") {
            for c in &sv.data {
                if c.contains("Assay5") { return AssayType::Assay5; }
                if c.contains("Assay") { return AssayType::Legacy; }
            }
        }
    }
    // Check structure
    if let Some(fields) = get_named_list_items(robj) {
        for (name, _) in &fields {
            if *name == "layers" { return AssayType::Assay5; }
        }
        for (name, _) in &fields {
            if *name == "counts" { return AssayType::Legacy; }
        }
    }
    AssayType::Unknown
}

fn robj_attributes(robj: &RObject) -> Option<&crate::rds::Attributes> {
    match robj {
        RObject::IntegerVector(v) => Some(&v.attributes),
        RObject::DoubleVector(v) => Some(&v.attributes),
        RObject::StringVector(v) => Some(&v.attributes),
        RObject::LogicalVector(v) => Some(&v.attributes),
        RObject::GenericVector(v) => Some(&v.attributes),
        RObject::PairList(v) => Some(&v.attributes),
        RObject::S4Object(v) => Some(&v.attributes),
        _ => None,
    }
}

// ============================================================
// dgCMatrix parsing (using new rds types)
// ============================================================

fn parse_dgcmatrix(robj: &RObject) -> Result<ExpressionMatrix, SeuratError> {
    let s4 = match robj {
        RObject::S4Object(s4) => s4,
        _ => return Err(SeuratError::ParseError(format!("Expected S4 for dgCMatrix, got {:?}", robj.sexp_type()))),
    };
    let i = s4.attributes.get("i").ok_or_else(|| SeuratError::ParseError("Missing 'i' in dgCMatrix".into()))?;
    let p = s4.attributes.get("p").ok_or_else(|| SeuratError::ParseError("Missing 'p' in dgCMatrix".into()))?;
    let x = s4.attributes.get("x").ok_or_else(|| SeuratError::ParseError("Missing 'x' in dgCMatrix".into()))?;
    let dim = s4.attributes.get("Dim").ok_or_else(|| SeuratError::ParseError("Missing 'Dim' in dgCMatrix".into()))?;

    let r_row_indices = extract_integer_vector(i)?;
    let r_col_ptrs = extract_integer_vector(p)?;
    let data = extract_real_vector(x)?;
    let dimensions = extract_integer_vector(dim)?;

    if dimensions.len() != 2 {
        return Err(SeuratError::ParseError(format!("Expected 2D dimensions, got {}", dimensions.len())));
    }

    let n_genes = dimensions[0] as usize;
    let n_cells = dimensions[1] as usize;

    let r_row_indices: Vec<usize> = r_row_indices.iter().map(|&x| x as usize).collect();
    let r_col_ptrs: Vec<usize> = r_col_ptrs.iter().map(|&x| x as usize).collect();

    use crate::ir::SparseMatrixCSR;
    use crate::sparse::convert::csr_to_csc;

    let csr = SparseMatrixCSR { n_rows: n_cells, n_cols: n_genes, indptr: r_col_ptrs, indices: r_row_indices, data };
    let csc = csr_to_csc(&csr);
    Ok(ExpressionMatrix::SparseCSC(csc))
}

// ============================================================
// Assay5 layer extraction (using new rds types)
// ============================================================

fn extract_assay5_layers(robj: &RObject) -> Result<HashMap<String, ExpressionMatrix>, SeuratError> {
    // Assay5 can be either a named list (simplified) or an S4 object (direct read)
    let layers_obj = if let RObject::S4Object(s4) = robj {
        // S4 Assay5: layers stored as a slot
        s4.attributes.get("layers")
            .ok_or_else(|| SeuratError::ParseError("Missing 'layers' slot in Assay5 S4".into()))?
    } else {
        // Named list: find layers field
        let fields = get_named_list_items(robj)
            .ok_or_else(|| SeuratError::ParseError(format!("Expected named list or S4 for Assay5, got {:?}", robj.sexp_type())))?;
        fields.iter().find(|(n, _)| *n == "layers").map(|(_, v)| *v)
            .ok_or_else(|| SeuratError::ParseError("Missing 'layers' in Assay5".into()))?
    };

    // layers itself can also be S4 or named list
    let layers_list = if let RObject::S4Object(s4) = layers_obj {
        // S4 layers container: look for .Data slot or iterate attributes
        if let Some(data_slot) = s4.attributes.get(".Data") {
            get_named_list_items(data_slot)
                .ok_or_else(|| SeuratError::ParseError("Expected named list in layers .Data".into()))?
        } else {
            // Fallback: treat S4 attributes as layers (excluding class-related ones)
            s4.attributes.names.iter().zip(s4.attributes.values.iter())
                .filter(|(n, _)| !matches!(n.as_str(), "class" | "names" | "package"))
                .map(|(n, v)| (n.as_str(), v))
                .collect()
        }
    } else {
        get_named_list_items(layers_obj)
            .ok_or_else(|| SeuratError::ParseError("Expected named list for layers".into()))?
    };

    let mut layers = HashMap::new();
    for (name, value) in layers_list {
        if matches!(value, RObject::Null) { continue; }
        if matches!(value, RObject::S4Object(_)) {
            if let Ok(matrix) = parse_dgcmatrix(value) {
                layers.insert(name.to_string(), matrix);
            }
        }
    }
    if layers.is_empty() {
        return Err(SeuratError::ParseError("No valid layers in Assay5".into()));
    }
    Ok(layers)
}

fn extract_assay5_main_matrix(robj: &RObject) -> Result<ExpressionMatrix, SeuratError> {
    let layers = extract_assay5_layers(robj)?;
    if let Some(m) = layers.get("counts") { return Ok(m.clone()); }
    if let Some(m) = layers.get("data") { return Ok(m.clone()); }
    layers.into_iter().next().map(|(_, m)| m)
        .ok_or_else(|| SeuratError::ParseError("No layers in Assay5".into()))
}

// ============================================================
// Seurat V3/V4 extraction
// ============================================================

fn extract_seurat_v3_v4(robj: &RObject, _file: &RdsFile, skipped: &mut SkippedComponents) -> Result<SingleCellData, SeuratError> {
    let assays = extract_slot(robj, "assays")?;
    let (expression, gene_metadata, layers) = extract_assays_v3_v4(assays, skipped)?;
    let (n_cells, n_genes) = expression.shape();
    let meta_data = extract_slot(robj, "meta.data")?;
    let cell_metadata = extract_metadata(meta_data, n_cells)?;
    let active_assay = extract_slot_opt(robj, "active.assay")?.and_then(|o| extract_string_scalar(o).ok());
    let embeddings = extract_slot_opt(robj, "reductions")?.and_then(|red| {
        if !should_skip_object(red, "reductions", skipped) { extract_reductions(red, n_cells, skipped).ok() } else { None }
    });
    skip_non_data_slots(robj, skipped)?;
    let spatial = extract_slot_opt(robj, "images")?.and_then(|imgs| extract_spatial_data(imgs).ok());
    let mut metadata = DatasetMetadata::new(n_cells, n_genes, "seurat".to_string());
    metadata.source_version = Some("4.0".to_string());
    metadata.active_assay = active_assay;
    let mut data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata).map_err(|e| SeuratError::ValidationError(e))?;
    data.layers = layers;
    data.embeddings = embeddings;
    data.spatial = spatial;
    Ok(data)
}

// ============================================================
// Seurat V5 extraction
// ============================================================

fn extract_seurat_v5(robj: &RObject, _file: &RdsFile, skipped: &mut SkippedComponents, keep_layers: bool) -> Result<SingleCellData, SeuratError> {
    let assays = extract_slot(robj, "assays")?;
    let (expression, gene_metadata, layers) = extract_assays_v5(assays, keep_layers)?;
    let (n_cells, n_genes) = expression.shape();
    let meta_data = extract_slot(robj, "meta.data")?;
    let cell_metadata = extract_metadata(meta_data, n_cells)?;
    let active_assay = extract_slot_opt(robj, "active.assay")?.and_then(|o| extract_string_scalar(o).ok());
    let embeddings = extract_slot_opt(robj, "reductions")?.and_then(|red| {
        if !should_skip_object(red, "reductions", skipped) { extract_reductions(red, n_cells, skipped).ok() } else { None }
    });
    skip_non_data_slots(robj, skipped)?;
    let spatial = extract_slot_opt(robj, "images")?.and_then(|imgs| extract_spatial_data(imgs).ok());
    let mut metadata = DatasetMetadata::new(n_cells, n_genes, "seurat".to_string());
    metadata.source_version = Some("5.0".to_string());
    metadata.active_assay = active_assay;
    let mut data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata).map_err(|e| SeuratError::ValidationError(e))?;
    data.layers = layers;
    data.embeddings = embeddings;
    data.spatial = spatial;
    Ok(data)
}

fn skip_non_data_slots(robj: &RObject, skipped: &mut SkippedComponents) -> Result<(), SeuratError> {
    for slot in &["graphs", "commands", "tools"] {
        if let Some(obj) = extract_slot_opt(robj, slot)? {
            if let Some(items) = get_named_list_items(obj) {
                for (name, _) in items { skipped.environments.push(format!("{}.{}", slot, name)); }
            }
        }
    }
    Ok(())
}

// ============================================================
// Assay extraction helpers
// ============================================================

fn extract_assays_v3_v4(robj: &RObject, skipped: &mut SkippedComponents)
    -> Result<(ExpressionMatrix, DataFrame, Option<HashMap<String, ExpressionMatrix>>), SeuratError> {
    let assay_list = get_named_list_items(robj).ok_or_else(|| SeuratError::ParseError("Expected named list for assays".into()))?;
    if assay_list.is_empty() { return Err(SeuratError::ParseError("No assays found".into())); }
    let (_, first_assay) = &assay_list[0];
    let (expression, gene_metadata) = extract_single_assay_v3_v4(first_assay, skipped)?;
    let layers = if assay_list.len() > 1 {
        let mut layer_map = HashMap::new();
        for (name, assay) in assay_list.iter().skip(1) {
            if let Ok((layer_expr, _)) = extract_single_assay_v3_v4(assay, skipped) { layer_map.insert(name.to_string(), layer_expr); }
        }
        if layer_map.is_empty() { None } else { Some(layer_map) }
    } else { None };
    Ok((expression, gene_metadata, layers))
}

fn extract_single_assay_v3_v4(robj: &RObject, skipped: &mut SkippedComponents) -> Result<(ExpressionMatrix, DataFrame), SeuratError> {
    let s4 = match robj {
        RObject::S4Object(s4) => s4,
        _ => return Err(SeuratError::ParseError(format!("Expected S4 for Assay, got {:?}", robj.sexp_type()))),
    };
    let counts = s4.attributes.get("counts").or_else(|| s4.attributes.get("data"));
    let expression = if let Some(matrix) = counts {
        if should_skip_object(matrix, "counts", skipped) { return Err(SeuratError::ParseError("Expression matrix is not serializable".into())); }
        parse_dgcmatrix(matrix)?
    } else { return Err(SeuratError::ParseError("Missing counts/data matrix in Assay".into())); };
    let (_, n_genes) = expression.shape();
    Ok((expression, DataFrame::empty(n_genes)))
}

fn extract_assays_v5(robj: &RObject, keep_layers: bool) -> Result<(ExpressionMatrix, DataFrame, Option<HashMap<String, ExpressionMatrix>>), SeuratError> {
    use crate::seurat::assay5::merge_split_layers;

    let assay_list = get_named_list_items(robj).ok_or_else(|| SeuratError::ParseError("Expected named list for assays".into()))?;
    if assay_list.is_empty() { return Err(SeuratError::ParseError("No assays found".into())); }
    let (_, first_assay) = &assay_list[0];

    // Get raw layers from the first assay (before merge)
    let raw_layers = extract_assay5_layers(first_assay)?;

    // Always merge to produce the main expression matrix (X)
    let merged_layers = merge_split_layers(raw_layers.clone())?;
    let expression = if let Some(m) = merged_layers.get("counts") {
        m.clone()
    } else if let Some(m) = merged_layers.get("data") {
        m.clone()
    } else {
        merged_layers.into_iter().next().map(|(_, m)| m)
            .ok_or_else(|| SeuratError::ParseError("No layers in Assay5".into()))?
    };

    let (_, n_genes) = expression.shape();

    // Build the layers map
    let mut layer_map = HashMap::new();

    // If keep_layers is true, preserve the original split layers (e.g. counts.1, counts.2)
    // Pad each sub-layer to full cell count so AnnData dimension requirements are met
    if keep_layers {
        use crate::seurat::assay5::pad_split_layers_to_full_size;

        let re = regex::Regex::new(r"^(.+)\..+$").unwrap();
        let has_split = raw_layers.keys().any(|k| {
            let prefix = re.captures(k).map(|c| c[1].to_string());
            if let Some(p) = prefix {
                raw_layers.keys().filter(|k2| {
                    re.captures(k2).map(|c2| c2[1].to_string()) == Some(p.clone())
                }).count() > 1
            } else {
                false
            }
        });
        if has_split {
            let (total_cells, _) = expression.shape();
            let padded = pad_split_layers_to_full_size(&raw_layers, total_cells, n_genes)?;
            for (name, matrix) in padded {
                layer_map.insert(name, matrix);
            }
            log::info!(
                "Preserved and padded {} split layer(s) with --keep-layers (total_cells={})",
                layer_map.len(),
                total_cells
            );
        }
    }

    // Also add layers from additional assays (e.g. SCT, ADT)
    for (name, assay) in assay_list.iter().skip(1) {
        if let Ok(layer_expr) = extract_assay5_main_matrix(assay) {
            layer_map.insert(name.to_string(), layer_expr);
        }
    }

    let layers = if layer_map.is_empty() { None } else { Some(layer_map) };
    Ok((expression, DataFrame::empty(n_genes), layers))
}

// ============================================================
// Metadata extraction
// ============================================================

fn extract_metadata(robj: &RObject, expected_rows: usize) -> Result<DataFrame, SeuratError> {
    use arrow::array::{ArrayRef, Int32Array, Float64Array, StringArray, BooleanArray, DictionaryArray};
    use arrow::datatypes::Int32Type;
    use std::sync::Arc;
    let columns = get_named_list_items(robj).unwrap_or_default();
    if columns.is_empty() { return Ok(DataFrame::empty(expected_rows)); }
    let mut col_names: Vec<String> = Vec::new();
    let mut col_arrays: Vec<ArrayRef> = Vec::new();
    for (col_name, col_value) in columns {
        if col_name.starts_with('_') { continue; }
        let array: ArrayRef = match col_value {
            RObject::IntegerVector(iv) => {
                if let Some(RObject::StringVector(sv)) = iv.attributes.get("levels") {
                    let keys: Vec<Option<i32>> = iv.data.iter().map(|&v| if v <= 0 || v == i32::MIN { None } else { Some(v - 1) }).collect();
                    let dict = DictionaryArray::<Int32Type>::try_new(Int32Array::from(keys), Arc::new(StringArray::from(sv.data.clone()))).map_err(|e| SeuratError::ParseError(e.to_string()))?;
                    Arc::new(dict)
                } else { Arc::new(Int32Array::from(iv.data.clone())) }
            }
            RObject::DoubleVector(dv) => Arc::new(Float64Array::from(dv.data.clone())),
            RObject::StringVector(sv) => Arc::new(StringArray::from(sv.data.clone())),
            RObject::LogicalVector(lv) => Arc::new(BooleanArray::from(lv.data.iter().map(|&v| v != 0).collect::<Vec<_>>())),
            _ => { let placeholder: Vec<String> = (0..expected_rows).map(|i| format!("Value_{}", i)).collect(); Arc::new(StringArray::from(placeholder)) }
        };
        col_names.push(col_name.to_string());
        col_arrays.push(array);
    }
    if col_names.is_empty() { return Ok(DataFrame::empty(expected_rows)); }
    DataFrame::new(col_names, col_arrays, expected_rows).map_err(|e| SeuratError::ParseError(format!("Failed to create DataFrame: {}", e)))
}

// ============================================================
// Reductions extraction
// ============================================================

fn extract_reductions(robj: &RObject, expected_rows: usize, skipped: &mut SkippedComponents) -> Result<HashMap<String, Embedding>, SeuratError> {
    if matches!(robj, RObject::Null) { return Ok(HashMap::new()); }
    let reduction_list = match get_named_list_items(robj) { Some(list) if !list.is_empty() => list, _ => return Ok(HashMap::new()), };
    let mut embeddings = HashMap::new();
    for (name, reduction) in reduction_list {
        if should_skip_object(reduction, &format!("reductions.{}", name), skipped) { continue; }
        if let Ok(embedding) = extract_single_reduction(reduction, expected_rows) { embeddings.insert(name.to_string(), embedding); }
    }
    Ok(embeddings)
}

fn extract_single_reduction(robj: &RObject, _expected_rows: usize) -> Result<Embedding, SeuratError> {
    let s4 = match robj {
        RObject::S4Object(s4) => s4,
        _ => return Err(SeuratError::ParseError(format!("Expected S4 for DimReduc, got {:?}", robj.sexp_type()))),
    };
    let cell_embeddings = s4.attributes.get("cell.embeddings").ok_or_else(|| SeuratError::ParseError("Missing cell.embeddings in DimReduc".into()))?;
    let key = s4.attributes.get("key").and_then(|o| extract_string_scalar(o).ok()).unwrap_or_else(|| "unknown_".to_string());
    let matrix = parse_dense_matrix(cell_embeddings)?;
    let n_rows = matrix.len();
    let n_cols = if n_rows > 0 { matrix[0].len() } else { 0 };
    let data: Vec<f64> = matrix.into_iter().flatten().collect();
    Embedding::new(key, data, n_rows, n_cols).map_err(|e| SeuratError::ValidationError(e))
}

fn parse_dense_matrix(robj: &RObject) -> Result<Vec<Vec<f64>>, SeuratError> {
    let (data, dim) = match robj {
        RObject::DoubleVector(dv) => {
            if let Some(dim_obj) = dv.attributes.get("dim") { (dv.data.clone(), extract_integer_vector(dim_obj)?) }
            else { return Ok(dv.data.iter().map(|&x| vec![x]).collect()); }
        }
        _ => return Err(SeuratError::ParseError(format!("Expected real matrix, got {:?}", robj.sexp_type()))),
    };
    if dim.len() != 2 { return Err(SeuratError::ParseError(format!("Expected 2D, got {}", dim.len()))); }
    let (n_rows, n_cols) = (dim[0] as usize, dim[1] as usize);
    if data.len() != n_rows * n_cols { return Err(SeuratError::ParseError(format!("Data len {} != {} x {}", data.len(), n_rows, n_cols))); }
    let mut matrix = vec![vec![0.0; n_cols]; n_rows];
    for col in 0..n_cols { for row in 0..n_rows { matrix[row][col] = data[col * n_rows + row]; } }
    Ok(matrix)
}

// ============================================================
// Spatial data extraction
// ============================================================

fn extract_spatial_data(robj: &RObject) -> Result<SpatialData, SeuratError> {
    let image_list = get_named_list_items(robj).ok_or_else(|| SeuratError::ParseError("Expected named list for images".into()))?;
    if image_list.is_empty() { return Err(SeuratError::ParseError("No images found".into())); }
    let (_, image_obj) = &image_list[0];
    extract_single_image(image_obj)
}

fn extract_single_image(robj: &RObject) -> Result<SpatialData, SeuratError> {
    let s4 = match robj {
        RObject::S4Object(s4) => s4,
        _ => return extract_spatial_from_list(robj),
    };
    let coordinates = s4.attributes.get("coordinates").ok_or_else(|| SeuratError::ParseError("Missing coordinates in spatial image".into()))?;
    let coord_matrix = parse_dense_matrix(coordinates)?;
    let n_cells = coord_matrix.len();
    let n_dims = if n_cells > 0 { coord_matrix[0].len() } else { 2 };
    let flat_coords: Vec<f64> = coord_matrix.into_iter().flatten().collect();
    let scale_factors = s4.attributes.get("scale.factors").and_then(|sf| extract_scale_factors(sf).ok());
    SpatialData::new(flat_coords, n_dims, None, scale_factors).map_err(|e| SeuratError::ValidationError(e))
}

fn extract_spatial_from_list(robj: &RObject) -> Result<SpatialData, SeuratError> {
    let items = get_named_list_items(robj).ok_or_else(|| SeuratError::ParseError("Expected named list for spatial".into()))?;
    let mut coordinates = None;
    let mut scale_factors = None;
    for (name, value) in items { match name { "coordinates" => coordinates = Some(value), "scale.factors" => scale_factors = Some(value), _ => {} } }
    let coords = coordinates.ok_or_else(|| SeuratError::ParseError("Missing coordinates".into()))?;
    let coord_matrix = parse_dense_matrix(coords)?;
    let n_cells = coord_matrix.len();
    let n_dims = if n_cells > 0 { coord_matrix[0].len() } else { 2 };
    let flat_coords: Vec<f64> = coord_matrix.into_iter().flatten().collect();
    let sf = scale_factors.and_then(|sf| extract_scale_factors(sf).ok());
    SpatialData::new(flat_coords, n_dims, None, sf).map_err(|e| SeuratError::ValidationError(e))
}

fn extract_scale_factors(robj: &RObject) -> Result<HashMap<String, f64>, SeuratError> {
    let items = get_named_list_items(robj).ok_or_else(|| SeuratError::ParseError("Expected named list for scale.factors".into()))?;
    let mut factors = HashMap::new();
    for (name, value) in items { if let RObject::DoubleVector(dv) = value { if !dv.data.is_empty() { factors.insert(name.to_string(), dv.data[0]); } } }
    Ok(factors)
}

// ============================================================
// Tests
// ============================================================

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_skipped_components() {
        let mut skipped = SkippedComponents::new();
        assert!(!skipped.has_skipped());
        skipped.environments.push("test".to_string());
        assert!(skipped.has_skipped());
        assert_eq!(skipped.total_skipped(), 1);
    }
    #[test]
    fn test_seurat_version_display() { assert_eq!(format!("{}", SeuratVersion::V5), "V5"); }
}
