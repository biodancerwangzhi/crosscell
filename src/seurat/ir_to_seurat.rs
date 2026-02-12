//! IR to Seurat conversion
//!
//! This module provides functions to convert CrossCell IR format to simplified Seurat objects.

use crate::ir::{DataFrame, Embedding, ExpressionMatrix, SingleCellData, SparseMatrixCSC};
use crate::rds::{
    Attributes, DoubleVector, GenericVector, IntegerVector, RObject, RdsFile, S4Object,
    StringVector,
};
use crate::rds::StringEncoding;
use crate::seurat::error::SeuratError;
use arrow::datatypes::Int32Type;
use std::collections::HashMap;

/// Convert IR to simplified Seurat RObject
///
/// This function converts a CrossCell IR representation to a simplified Seurat object
/// that can be written to RDS format.
pub fn ir_to_seurat_rds(ir: &SingleCellData) -> Result<(RObject, RdsFile), SeuratError> {
    let file = RdsFile::default();
    
    // Create simplified Seurat structure as a named list
    let mut data = Vec::new();
    let mut names = Vec::new();

    // 1. project_name
    names.push("project_name".to_string());
    data.push(create_string_scalar("CrossCell"));

    // 2. active_assay
    let active_assay = ir
        .metadata
        .active_assay
        .clone()
        .unwrap_or_else(|| "RNA".to_string());
    names.push("active_assay".to_string());
    data.push(create_string_scalar(&active_assay));

    // 3. assays
    names.push("assays".to_string());
    data.push(construct_assays(ir, &active_assay)?);

    // 4. meta_data
    names.push("meta_data".to_string());
    data.push(construct_metadata(&ir.cell_metadata)?);

    // Get cell count for validation
    let n_cells = ir.metadata.n_cells;

    // 5. reductions (optional)
    if let Some(embeddings) = &ir.embeddings {
        names.push("reductions".to_string());
        data.push(construct_reductions(embeddings)?);
    }

    // 6. images (optional)
    if let Some(spatial) = &ir.spatial {
        names.push("images".to_string());
        data.push(construct_spatial_data(spatial, n_cells)?);
    }

    // 7. class
    names.push("class".to_string());
    data.push(create_string_scalar("SimplifiedSeurat"));

    // Create named list with class attribute
    let mut attrs = Attributes::new();
    let mut names_vec = StringVector::default();
    for name in &names {
        names_vec.add(name.clone(), StringEncoding::Utf8);
    }
    attrs.add("names".to_string(), RObject::StringVector(names_vec), StringEncoding::Utf8);
    
    // Add class attribute
    let mut class_vec = StringVector::default();
    class_vec.add("SimplifiedSeurat".to_string(), StringEncoding::Utf8);
    class_vec.add("list".to_string(), StringEncoding::Utf8);
    attrs.add("class".to_string(), RObject::StringVector(class_vec), StringEncoding::Utf8);

    let seurat_list = RObject::GenericVector(GenericVector {
        data,
        attributes: attrs,
    });

    Ok((seurat_list, file))
}

/// Create a string scalar RObject
fn create_string_scalar(s: &str) -> RObject {
    let mut sv = StringVector::default();
    sv.add(s.to_string(), StringEncoding::Utf8);
    RObject::StringVector(sv)
}

/// Create an integer vector RObject
fn create_integer_vector(data: Vec<i32>) -> RObject {
    RObject::IntegerVector(IntegerVector {
        data,
        attributes: Attributes::new(),
    })
}

/// Create a real vector RObject
fn create_real_vector(data: Vec<f64>) -> RObject {
    RObject::DoubleVector(DoubleVector {
        data,
        attributes: Attributes::new(),
    })
}

/// Create a string vector RObject
fn create_string_vector(data: Vec<String>) -> RObject {
    let mut sv = StringVector::default();
    for s in data {
        sv.add(s, StringEncoding::Utf8);
    }
    RObject::StringVector(sv)
}

/// Create a named list RObject
fn create_named_list(items: Vec<(String, RObject)>) -> RObject {
    let mut data = Vec::new();
    let mut names = Vec::new();
    
    for (name, value) in items {
        names.push(name);
        data.push(value);
    }
    
    let mut attrs = Attributes::new();
    let mut names_vec = StringVector::default();
    for name in names {
        names_vec.add(name, StringEncoding::Utf8);
    }
    attrs.add("names".to_string(), RObject::StringVector(names_vec), StringEncoding::Utf8);
    
    RObject::GenericVector(GenericVector {
        data,
        attributes: attrs,
    })
}

/// Construct assays field (named list of assay objects)
fn construct_assays(ir: &SingleCellData, active_assay: &str) -> Result<RObject, SeuratError> {
    let mut assay_list = Vec::new();

    // Main assay (from primary expression matrix)
    let main_assay = construct_single_assay(&ir.expression, &ir.gene_metadata)?;
    assay_list.push((active_assay.to_string(), main_assay));

    // Additional assays (from layers)
    if let Some(layers) = &ir.layers {
        for (name, layer_expr) in layers {
            let layer_assay = construct_single_assay(layer_expr, &ir.gene_metadata)?;
            assay_list.push((name.clone(), layer_assay));
        }
    }

    Ok(create_named_list(assay_list))
}

/// Construct single assay object
fn construct_single_assay(
    expression: &ExpressionMatrix,
    gene_metadata: &DataFrame,
) -> Result<RObject, SeuratError> {
    let mut assay_fields = Vec::new();

    // class: "Assay5"
    assay_fields.push((
        "class".to_string(),
        create_string_scalar("Assay5"),
    ));

    // counts: dgCMatrix
    let counts = expression_to_dgcmatrix(expression)?;
    assay_fields.push(("counts".to_string(), counts));

    // features: character vector of gene names
    let n_genes = gene_metadata.n_rows;
    let features = if n_genes > 0 {
        let feature_names: Vec<String> = (0..n_genes)
            .map(|i| format!("Gene_{}", i + 1))
            .collect();
        create_string_vector(feature_names)
    } else {
        create_string_vector(vec![])
    };
    assay_fields.push(("features".to_string(), features));

    // cells: character vector of cell names
    let n_cells = expression.shape().0;
    let cells = if n_cells > 0 {
        let cell_names: Vec<String> = (0..n_cells)
            .map(|i| format!("Cell_{}", i + 1))
            .collect();
        create_string_vector(cell_names)
    } else {
        create_string_vector(vec![])
    };
    assay_fields.push(("cells".to_string(), cells));

    Ok(create_named_list(assay_fields))
}

/// Convert ExpressionMatrix to dgCMatrix (R sparse matrix)
pub fn expression_to_dgcmatrix(expression: &ExpressionMatrix) -> Result<RObject, SeuratError> {
    match expression {
        ExpressionMatrix::SparseCSC(sparse) => sparse_csc_to_dgcmatrix(sparse),
        ExpressionMatrix::SparseCSR(sparse) => {
            use crate::sparse::convert::csr_to_csc;
            let csc = csr_to_csc(sparse);
            sparse_csc_to_dgcmatrix(&csc)
        }
        ExpressionMatrix::Dense(dense) => {
            dense_to_dgcmatrix(dense)
        }
        ExpressionMatrix::Lazy(lazy) => {
            if let Some(cached) = lazy.get_cached() {
                expression_to_dgcmatrix(&cached)
            } else {
                Err(SeuratError::ConversionError(
                    "Cannot convert LazyMatrix without loading data first.".to_string()
                ))
            }
        }
    }
}

/// Convert DenseMatrix to dgCMatrix
fn dense_to_dgcmatrix(dense: &crate::ir::DenseMatrix) -> Result<RObject, SeuratError> {
    let n_rows = dense.n_rows;
    let n_cols = dense.n_cols;
    
    let mut data = Vec::new();
    let mut indices = Vec::new();
    let mut indptr = vec![0usize];
    
    for col in 0..n_cols {
        for row in 0..n_rows {
            let val = dense.data[row * n_cols + col];
            if val != 0.0 {
                data.push(val);
                indices.push(row);
            }
        }
        indptr.push(data.len());
    }
    
    let csc = SparseMatrixCSC {
        data,
        indices,
        indptr,
        n_rows,
        n_cols,
    };
    
    sparse_csc_to_dgcmatrix(&csc)
}


/// Convert SparseMatrixCSC to dgCMatrix
fn sparse_csc_to_dgcmatrix(sparse: &SparseMatrixCSC) -> Result<RObject, SeuratError> {
    log::debug!("sparse_csc_to_dgcmatrix: IR input: {} cells × {} genes", sparse.n_rows, sparse.n_cols);
    
    // dgCMatrix in R is CSC format, but stores genes × cells (transposed from our cells × genes)
    // We need to convert our CSC (cells × genes) to CSC (genes × cells)
    // This is equivalent to transposing, which means converting to CSR then treating as CSC
    use crate::sparse::convert::csc_to_csr;
    let csr = csc_to_csr(sparse);
    
    log::debug!("sparse_csc_to_dgcmatrix: CSR result: {} × {}", csr.n_rows, csr.n_cols);
    
    // Create S4 object for dgCMatrix
    let mut s4 = S4Object::default();
    s4.class_name = "dgCMatrix".to_string();
    s4.class_encoding = StringEncoding::Utf8;
    s4.package_name = "Matrix".to_string();
    s4.package_encoding = StringEncoding::Utf8;

    // i: row indices (0-based in R)
    let i_vec: Vec<i32> = csr.indices.iter().map(|&x| x as i32).collect();
    s4.attributes.add("i".to_string(), create_integer_vector(i_vec), StringEncoding::Utf8);

    // p: column pointers
    let p_vec: Vec<i32> = csr.indptr.iter().map(|&x| x as i32).collect();
    s4.attributes.add("p".to_string(), create_integer_vector(p_vec), StringEncoding::Utf8);

    // x: non-zero values
    s4.attributes.add("x".to_string(), create_real_vector(csr.data.clone()), StringEncoding::Utf8);

    // Dim: dimensions [n_genes, n_cells]
    let dim_vec = vec![sparse.n_cols as i32, sparse.n_rows as i32];
    log::debug!("sparse_csc_to_dgcmatrix: Dim = [{}, {}] (genes × cells)", sparse.n_cols, sparse.n_rows);
    s4.attributes.add("Dim".to_string(), create_integer_vector(dim_vec), StringEncoding::Utf8);

    // Dimnames: NULL
    s4.attributes.add("Dimnames".to_string(), RObject::Null, StringEncoding::Utf8);

    // factors: empty list
    s4.attributes.add("factors".to_string(), RObject::GenericVector(GenericVector::default()), StringEncoding::Utf8);

    Ok(RObject::S4Object(s4))
}

/// Construct metadata (data.frame)
fn construct_metadata(df: &DataFrame) -> Result<RObject, SeuratError> {
    if df.n_rows == 0 {
        // Empty data.frame
        let mut attrs = Attributes::new();
        attrs.add("row.names".to_string(), create_integer_vector(vec![]), StringEncoding::Utf8);
        
        let mut class_vec = StringVector::default();
        class_vec.add("data.frame".to_string(), StringEncoding::Utf8);
        attrs.add("class".to_string(), RObject::StringVector(class_vec), StringEncoding::Utf8);

        return Ok(RObject::GenericVector(GenericVector {
            data: vec![],
            attributes: attrs,
        }));
    }

    // Convert DataFrame columns to R vectors
    let mut columns = Vec::new();
    
    for col_name in &df.columns {
        if let Some(array) = df.column(col_name) {
            let r_col = arrow_array_to_robject(array, col_name)?;
            columns.push((col_name.clone(), r_col));
        }
    }
    
    // If no columns were converted, add a dummy column
    if columns.is_empty() {
        let dummy_col: Vec<i32> = (0..df.n_rows as i32).collect();
        columns.push(("_row_id".to_string(), create_integer_vector(dummy_col)));
    }

    // Extract column names
    let col_names: Vec<String> = columns.iter().map(|(n, _)| n.clone()).collect();
    
    // Create data list
    let data: Vec<RObject> = columns.into_iter().map(|(_, v)| v).collect();
    
    // Create attributes
    let mut attrs = Attributes::new();
    
    // row.names attribute (1-based in R)
    let row_names: Vec<i32> = (1..=df.n_rows as i32).collect();
    attrs.add("row.names".to_string(), create_integer_vector(row_names), StringEncoding::Utf8);
    
    // class attribute
    let mut class_vec = StringVector::default();
    class_vec.add("data.frame".to_string(), StringEncoding::Utf8);
    attrs.add("class".to_string(), RObject::StringVector(class_vec), StringEncoding::Utf8);
    
    // names attribute (column names)
    attrs.add("names".to_string(), create_string_vector(col_names), StringEncoding::Utf8);

    Ok(RObject::GenericVector(GenericVector {
        data,
        attributes: attrs,
    }))
}

/// Convert Arrow array to RObject
fn arrow_array_to_robject(array: &dyn arrow::array::Array, col_name: &str) -> Result<RObject, SeuratError> {
    use arrow::array::*;
    use arrow::datatypes::DataType;
    
    match array.data_type() {
        DataType::Int64 => {
            let arr = array.as_any().downcast_ref::<Int64Array>()
                .ok_or_else(|| SeuratError::ConversionError(format!("Failed to downcast {} to Int64Array", col_name)))?;
            let values: Vec<i32> = arr.iter()
                .map(|v| v.map(|x| x as i32).unwrap_or(i32::MIN))
                .collect();
            Ok(create_integer_vector(values))
        }
        DataType::Int32 => {
            let arr = array.as_any().downcast_ref::<Int32Array>()
                .ok_or_else(|| SeuratError::ConversionError(format!("Failed to downcast {} to Int32Array", col_name)))?;
            let values: Vec<i32> = arr.iter()
                .map(|v| v.unwrap_or(i32::MIN))
                .collect();
            Ok(create_integer_vector(values))
        }
        DataType::Float64 => {
            let arr = array.as_any().downcast_ref::<Float64Array>()
                .ok_or_else(|| SeuratError::ConversionError(format!("Failed to downcast {} to Float64Array", col_name)))?;
            let values: Vec<f64> = arr.iter()
                .map(|v| v.unwrap_or(f64::NAN))
                .collect();
            Ok(create_real_vector(values))
        }
        DataType::Float32 => {
            let arr = array.as_any().downcast_ref::<Float32Array>()
                .ok_or_else(|| SeuratError::ConversionError(format!("Failed to downcast {} to Float32Array", col_name)))?;
            let values: Vec<f64> = arr.iter()
                .map(|v| v.map(|x| x as f64).unwrap_or(f64::NAN))
                .collect();
            Ok(create_real_vector(values))
        }
        DataType::Utf8 => {
            let arr = array.as_any().downcast_ref::<StringArray>()
                .ok_or_else(|| SeuratError::ConversionError(format!("Failed to downcast {} to StringArray", col_name)))?;
            let values: Vec<String> = arr.iter()
                .map(|v| v.map(|s| s.to_string()).unwrap_or_else(|| "NA".to_string()))
                .collect();
            Ok(create_string_vector(values))
        }
        DataType::LargeUtf8 => {
            let arr = array.as_any().downcast_ref::<LargeStringArray>()
                .ok_or_else(|| SeuratError::ConversionError(format!("Failed to downcast {} to LargeStringArray", col_name)))?;
            let values: Vec<String> = arr.iter()
                .map(|v| v.map(|s| s.to_string()).unwrap_or_else(|| "NA".to_string()))
                .collect();
            Ok(create_string_vector(values))
        }
        DataType::Boolean => {
            let arr = array.as_any().downcast_ref::<BooleanArray>()
                .ok_or_else(|| SeuratError::ConversionError(format!("Failed to downcast {} to BooleanArray", col_name)))?;
            // R logical is stored as integer with special values
            let values: Vec<i32> = arr.iter()
                .map(|v| match v {
                    Some(true) => 1,
                    Some(false) => 0,
                    None => i32::MIN,  // NA
                })
                .collect();
            Ok(create_integer_vector(values))
        }
        DataType::Dictionary(key_type, _value_type) => {
            // Categorical/Factor data
            match key_type.as_ref() {
                DataType::Int32 => {
                    let dict_arr = array.as_any().downcast_ref::<DictionaryArray<Int32Type>>()
                        .ok_or_else(|| SeuratError::ConversionError(format!("Failed to downcast {} to DictionaryArray", col_name)))?;
                    
                    // Get keys (0-based) and convert to 1-based for R Factor
                    let keys: Vec<i32> = dict_arr.keys().iter()
                        .map(|v: Option<i32>| v.map(|k| k + 1).unwrap_or(i32::MIN))  // 0-based to 1-based, NA as MIN
                        .collect();
                    
                    // Get levels (categories)
                    let values = dict_arr.values();
                    let levels: Vec<String> = if let Some(str_arr) = values.as_any().downcast_ref::<StringArray>() {
                        str_arr.iter().map(|v: Option<&str>| v.unwrap_or("").to_string()).collect()
                    } else {
                        // Fallback: generate level names
                        (0..values.len()).map(|i| format!("Level_{}", i)).collect()
                    };
                    
                    // Create Factor (integer vector with levels and class attributes)
                    let mut int_vec = IntegerVector {
                        data: keys,
                        attributes: Attributes::new(),
                    };
                    
                    // Add levels attribute
                    int_vec.attributes.add("levels".to_string(), create_string_vector(levels), StringEncoding::Utf8);
                    
                    // Add class attribute
                    let mut class_vec = StringVector::default();
                    class_vec.add("factor".to_string(), StringEncoding::Utf8);
                    int_vec.attributes.add("class".to_string(), RObject::StringVector(class_vec), StringEncoding::Utf8);
                    
                    Ok(RObject::IntegerVector(int_vec))
                }
                _ => {
                    // Unsupported dictionary key type, convert to string
                    eprintln!("Warning: Unsupported dictionary key type for column {}, converting to string", col_name);
                    let values: Vec<String> = (0..array.len())
                        .map(|i| format!("Value_{}", i))
                        .collect();
                    Ok(create_string_vector(values))
                }
            }
        }
        _ => {
            // Unsupported type, convert to string representation
            eprintln!("Warning: Unsupported data type {:?} for column {}, converting to string", array.data_type(), col_name);
            let values: Vec<String> = (0..array.len())
                .map(|i| format!("Value_{}", i))
                .collect();
            Ok(create_string_vector(values))
        }
    }
}


/// Construct reductions (named list of reduction objects)
fn construct_reductions(
    embeddings: &HashMap<String, Embedding>,
) -> Result<RObject, SeuratError> {
    let mut reduction_list = Vec::new();

    for (name, embedding) in embeddings {
        let reduction = construct_single_reduction(embedding)?;
        reduction_list.push((name.clone(), reduction));
    }

    Ok(create_named_list(reduction_list))
}

/// Construct single reduction object
fn construct_single_reduction(embedding: &Embedding) -> Result<RObject, SeuratError> {
    let mut reduction_fields = Vec::new();

    // key: reduction name
    reduction_fields.push((
        "key".to_string(),
        create_string_scalar(&embedding.name),
    ));

    // cell_embeddings: matrix
    let cell_embeddings = embedding_to_matrix(embedding)?;
    reduction_fields.push(("cell_embeddings".to_string(), cell_embeddings));

    Ok(create_named_list(reduction_fields))
}

/// Convert Embedding to R matrix
fn embedding_to_matrix(embedding: &Embedding) -> Result<RObject, SeuratError> {
    // R matrix is stored as a real vector with dim attribute
    // R stores matrices in column-major order
    
    let n_rows = embedding.n_rows;
    let n_cols = embedding.n_cols;
    
    // Convert from row-major (our format) to column-major (R format)
    let mut col_major_data = Vec::with_capacity(n_rows * n_cols);
    for col in 0..n_cols {
        for row in 0..n_rows {
            col_major_data.push(embedding.data[row * n_cols + col]);
        }
    }

    let mut double_vec = DoubleVector {
        data: col_major_data,
        attributes: Attributes::new(),
    };

    // Add dim attribute
    let dim_vec = vec![n_rows as i32, n_cols as i32];
    double_vec.attributes.add("dim".to_string(), create_integer_vector(dim_vec), StringEncoding::Utf8);

    Ok(RObject::DoubleVector(double_vec))
}

/// Construct spatial data (images slot)
fn construct_spatial_data(
    spatial: &crate::ir::SpatialData,
    expected_cells: usize,
) -> Result<RObject, SeuratError> {
    // Verify cell count matches
    if spatial.n_cells() != expected_cells {
        return Err(SeuratError::ValidationError(format!(
            "Spatial data cell count {} doesn't match expected {}",
            spatial.n_cells(),
            expected_cells
        )));
    }

    // Create named list of image objects
    // For now, create a single image object named "slice1" (standard for Visium)
    let image_obj = construct_image_object(spatial)?;
    let mut images = Vec::new();
    images.push(("slice1".to_string(), image_obj));

    Ok(create_named_list(images))
}

/// Construct single image object
fn construct_image_object(spatial: &crate::ir::SpatialData) -> Result<RObject, SeuratError> {
    let mut image_fields = Vec::new();

    // 1. coordinates: matrix (cells × 2 or cells × 3)
    let coordinates = construct_spatial_coordinates(spatial)?;
    image_fields.push(("coordinates".to_string(), coordinates));

    // 2. scale.factors: named list (optional)
    if let Some(scale_factors) = &spatial.scale_factors {
        let scale_factors_obj = construct_scale_factors(scale_factors)?;
        image_fields.push(("scale.factors".to_string(), scale_factors_obj));
    }

    // 3. image: array (optional)
    if let Some(images) = &spatial.images {
        if !images.is_empty() {
            // Use first image
            let image_data = construct_image_data(&images[0])?;
            image_fields.push(("image".to_string(), image_data));
        }
    }

    Ok(create_named_list(image_fields))
}

/// Construct spatial coordinates matrix
fn construct_spatial_coordinates(spatial: &crate::ir::SpatialData) -> Result<RObject, SeuratError> {
    let n_cells = spatial.n_cells();
    let n_dims = spatial.n_dims;

    // Convert flat coordinates to matrix (row-major to column-major)
    let mut col_major_data = Vec::with_capacity(n_cells * n_dims);
    for col in 0..n_dims {
        for row in 0..n_cells {
            col_major_data.push(spatial.coordinates[row * n_dims + col]);
        }
    }

    let mut double_vec = DoubleVector {
        data: col_major_data,
        attributes: Attributes::new(),
    };

    // Add dim attribute
    let dim_vec = vec![n_cells as i32, n_dims as i32];
    double_vec.attributes.add("dim".to_string(), create_integer_vector(dim_vec), StringEncoding::Utf8);

    Ok(RObject::DoubleVector(double_vec))
}

/// Construct scale factors (named list)
fn construct_scale_factors(
    scale_factors: &std::collections::HashMap<String, f64>,
) -> Result<RObject, SeuratError> {
    let mut factors = Vec::new();

    for (name, value) in scale_factors {
        factors.push((name.clone(), create_real_vector(vec![*value])));
    }

    Ok(create_named_list(factors))
}

/// Construct image data (array)
fn construct_image_data(image: &crate::ir::SpatialImage) -> Result<RObject, SeuratError> {
    // For now, create a placeholder
    // TODO: Implement actual image data construction
    // This would involve decoding the image bytes and creating an R array
    
    // Return empty array for now
    let mut double_vec = DoubleVector {
        data: Vec::new(),
        attributes: Attributes::new(),
    };
    
    // Add dim attribute (height × width × channels)
    let dim_vec = vec![image.height as i32, image.width as i32, 3]; // Assume RGB
    double_vec.attributes.add("dim".to_string(), create_integer_vector(dim_vec), StringEncoding::Utf8);

    Ok(RObject::DoubleVector(double_vec))
}

/// Write IR to simplified Seurat RDS file
///
/// This is a convenience function that combines IR to Seurat conversion
/// and RDS file writing.
pub fn write_seurat_rds(ir: &SingleCellData, path: &str) -> Result<(), SeuratError> {
    use crate::rds::write_rds;
    use std::path::Path;

    // Convert IR to Seurat RObject
    let (seurat_robj, mut file) = ir_to_seurat_rds(ir)?;

    // Set the object in the file
    file.object = seurat_robj;
    
    write_rds(&file, Path::new(path))
        .map_err(|e| SeuratError::IoError(format!("Failed to write RDS: {}", e)))?;

    Ok(())
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::{DatasetMetadata, ExpressionMatrix, SparseMatrixCSC};

    #[test]
    fn test_ir_to_seurat_minimal() {
        // Create minimal IR
        // CSC format: indptr length = n_cols + 1
        let n_cells = 10;
        let n_genes = 5;
        
        // Create a simple sparse matrix with some non-zero values
        // Each gene has 2 non-zero values
        let sparse = SparseMatrixCSC {
            n_rows: n_cells,  // 10 cells
            n_cols: n_genes,  // 5 genes
            indptr: vec![0, 2, 4, 6, 8, 10],  // length = n_genes + 1 = 6
            indices: vec![0, 5, 1, 6, 2, 7, 3, 8, 4, 9],  // cell indices
            data: vec![1.0; 10],
        };
        let expression = ExpressionMatrix::SparseCSC(sparse);
        let cell_metadata = DataFrame::empty(n_cells);
        let gene_metadata = DataFrame::empty(n_genes);
        let metadata = DatasetMetadata::new(n_cells, n_genes, "test".to_string());

        let ir = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata)
            .expect("Failed to create IR");

        // Convert to Seurat
        let result = ir_to_seurat_rds(&ir);
        assert!(result.is_ok(), "Failed to convert IR to Seurat: {:?}", result.err());

        let (seurat, _file) = result.unwrap();
        match seurat {
            RObject::GenericVector(gv) => {
                // Check that required fields are present
                let field_names = gv.attributes.get_names();
                assert!(field_names.is_some());
                let names = field_names.unwrap();
                assert!(names.contains(&"project_name".to_string()));
                assert!(names.contains(&"active_assay".to_string()));
                assert!(names.contains(&"assays".to_string()));
                assert!(names.contains(&"meta_data".to_string()));
            }
            _ => panic!("Expected List"),
        }
    }

    #[test]
    fn test_write_seurat_rds() {
        // Create minimal IR
        // CSC format: indptr length = n_cols + 1
        let n_cells = 5;
        let n_genes = 10;
        
        let sparse = SparseMatrixCSC {
            n_rows: n_cells,
            n_cols: n_genes,
            indptr: vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],  // length = n_genes + 1 = 11
            indices: vec![0, 1, 2, 3, 4, 0, 1, 2, 3, 4],  // cell indices
            data: vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
        };
        let expression = ExpressionMatrix::SparseCSC(sparse);
        let cell_metadata = DataFrame::empty(n_cells);
        let gene_metadata = DataFrame::empty(n_genes);
        let metadata = DatasetMetadata::new(n_cells, n_genes, "test".to_string());

        let ir = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata)
            .expect("Failed to create IR");

        // Write to RDS file
        let output_path = "tests/data/rust_generated_seurat.rds";
        let result = write_seurat_rds(&ir, output_path);
        assert!(result.is_ok(), "Failed to write Seurat RDS: {:?}", result.err());

        // Verify file exists
        use std::path::Path;
        assert!(Path::new(output_path).exists(), "Output file not created");
    }

    #[test]
    fn test_sparse_csc_to_dgcmatrix() {
        // CSC format: indptr length = n_cols + 1
        let n_cells = 3;
        let n_genes = 5;
        
        let sparse = SparseMatrixCSC {
            n_rows: n_cells,  // 3 cells
            n_cols: n_genes,  // 5 genes
            indptr: vec![0, 1, 2, 3, 4, 5],  // length = n_genes + 1 = 6
            indices: vec![0, 1, 2, 0, 1],  // cell indices
            data: vec![1.0, 2.0, 3.0, 4.0, 5.0],
        };

        let result = sparse_csc_to_dgcmatrix(&sparse);
        assert!(result.is_ok(), "Failed to convert to dgCMatrix: {:?}", result.err());

        match result.unwrap() {
            RObject::S4Object(s4) => {
                assert_eq!(s4.class_name, "dgCMatrix");
                
                // Check dimensions are swapped
                let dim_attr = s4.attributes.get("Dim");
                assert!(dim_attr.is_some());
                
                if let Some(RObject::IntegerVector(int_vec)) = dim_attr {
                    assert_eq!(int_vec.data[0], 5);  // genes (was n_cols)
                    assert_eq!(int_vec.data[1], 3);  // cells (was n_rows)
                }
            }
            _ => panic!("Expected S4 object"),
        }
    }

    #[test]
    fn test_embedding_to_matrix() {
        let embedding = Embedding::new(
            "pca".to_string(),
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            3,  // 3 rows
            2,  // 2 columns
        ).expect("Failed to create embedding");

        let result = embedding_to_matrix(&embedding);
        assert!(result.is_ok(), "Failed to convert embedding to matrix: {:?}", result.err());

        match result.unwrap() {
            RObject::DoubleVector(dv) => {
                // Check data is in column-major order
                // Original row-major: [1, 2, 3, 4, 5, 6]
                // Column-major: [1, 3, 5, 2, 4, 6]
                assert_eq!(dv.data, vec![1.0, 3.0, 5.0, 2.0, 4.0, 6.0]);

                // Check dim attribute
                let dim_attr = dv.attributes.get("dim");
                assert!(dim_attr.is_some());
                
                if let Some(RObject::IntegerVector(int_vec)) = dim_attr {
                    assert_eq!(int_vec.data, vec![3, 2]);
                }
            }
            _ => panic!("Expected Double vector"),
        }
    }

    #[test]
    fn test_construct_spatial_data() {
        use crate::ir::SpatialData;
        use std::collections::HashMap;

        // Create spatial data (3 cells, 2D coordinates)
        let coordinates = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let mut scale_factors = HashMap::new();
        scale_factors.insert("tissue_hires_scalef".to_string(), 0.08);
        scale_factors.insert("spot_diameter_fullres".to_string(), 89.43);
        scale_factors.insert("fiducial_diameter_fullres".to_string(), 144.0);

        let spatial = SpatialData::new(coordinates, 2, None, Some(scale_factors))
            .expect("Failed to create spatial data");

        // Construct Seurat images object
        let result = construct_spatial_data(&spatial, 3);
        assert!(result.is_ok(), "Failed to construct spatial data: {:?}", result.err());

        match result.unwrap() {
            RObject::GenericVector(gv) => {
                // Should have one image named "slice1"
                let names = gv.attributes.get_names();
                assert!(names.is_some());
                assert_eq!(names.unwrap().len(), 1);
                assert_eq!(names.unwrap()[0], "slice1");
            }
            _ => panic!("Expected List for images"),
        }
    }

    #[test]
    fn test_construct_spatial_coordinates() {
        use crate::ir::SpatialData;

        // Create spatial data (3 cells, 2D)
        let coordinates = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let spatial = SpatialData::new(coordinates, 2, None, None)
            .expect("Failed to create spatial data");

        // Construct coordinates matrix
        let result = construct_spatial_coordinates(&spatial);
        assert!(result.is_ok(), "Failed to construct coordinates: {:?}", result.err());

        match result.unwrap() {
            RObject::DoubleVector(dv) => {
                // Check data is in column-major order
                // Original row-major: [1, 2, 3, 4, 5, 6]
                // Column-major: [1, 3, 5, 2, 4, 6]
                assert_eq!(dv.data, vec![1.0, 3.0, 5.0, 2.0, 4.0, 6.0]);

                // Check dim attribute
                let dim_attr = dv.attributes.get("dim");
                assert!(dim_attr.is_some());
                
                if let Some(RObject::IntegerVector(int_vec)) = dim_attr {
                    assert_eq!(int_vec.data, vec![3, 2]); // 3 cells × 2 dims
                }
            }
            _ => panic!("Expected Double vector"),
        }
    }

    #[test]
    fn test_construct_scale_factors() {
        use std::collections::HashMap;

        let mut scale_factors = HashMap::new();
        scale_factors.insert("tissue_hires_scalef".to_string(), 0.08);
        scale_factors.insert("spot_diameter_fullres".to_string(), 89.43);

        let result = construct_scale_factors(&scale_factors);
        assert!(result.is_ok(), "Failed to construct scale factors: {:?}", result.err());

        match result.unwrap() {
            RObject::GenericVector(gv) => {
                assert_eq!(gv.data.len(), 2);
                
                // Check that both factors are present
                let names = gv.attributes.get_names();
                assert!(names.is_some());
                let names = names.unwrap();
                assert!(names.contains(&"tissue_hires_scalef".to_string()) || 
                        names.contains(&"spot_diameter_fullres".to_string()));
            }
            _ => panic!("Expected List"),
        }
    }
}
