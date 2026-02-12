//! User-friendly error message formatting
//!
//! Provides actionable error messages with suggestions for common issues.

use crate::error::CrossCellError;

/// Format an error with user-friendly suggestions
pub fn format_error_with_suggestions(error: &CrossCellError) -> String {
    let base_message = error.to_string();
    let suggestion = get_suggestion(error);
    
    if let Some(suggestion) = suggestion {
        format!("{}\n\n💡 {}", base_message, suggestion)
    } else {
        base_message
    }
}

/// Get a suggestion for a specific error type
fn get_suggestion(error: &CrossCellError) -> Option<String> {
    match error {
        CrossCellError::FileNotFound { path } => {
            Some(format!(
                "Solution: Check that the file exists and the path is correct.\n\
                 Try: ls -la {}\n\
                 Or use an absolute path to the file.",
                path
            ))
        }
        
        CrossCellError::InvalidFormat(msg) => {
            if msg.contains("h5ad") || msg.contains("HDF5") {
                Some(
                    "Solution: Ensure the file is a valid AnnData (.h5ad) file.\n\
                     Try: crosscell inspect -i <file> to check file format.\n\
                     If the file is corrupted, try re-exporting from Python:\n\
                       import anndata\n\
                       adata = anndata.read_h5ad('file.h5ad')\n\
                       adata.write('file_fixed.h5ad')".to_string()
                )
            } else if msg.contains("rds") || msg.contains("RDS") {
                Some(
                    "Solution: Ensure the file is a valid Seurat RDS file.\n\
                     Try: crosscell inspect -i <file> to check file format.\n\
                     If the file is corrupted, try re-saving from R:\n\
                       seurat_obj <- readRDS('file.rds')\n\
                       saveRDS(seurat_obj, 'file_fixed.rds')".to_string()
                )
            } else {
                Some(
                    "Solution: CrossCell supports .h5ad (AnnData) and .rds (Seurat) formats.\n\
                     Check that your file has the correct extension.".to_string()
                )
            }
        }
        
        CrossCellError::UnsupportedFormat { expected, actual } => {
            Some(format!(
                "Solution: Convert your file to a supported format.\n\
                 Expected: {}\n\
                 Got: {}\n\
                 Supported formats: .h5ad (AnnData), .rds (Seurat)",
                expected, actual
            ))
        }
        
        CrossCellError::MissingField { field } => {
            let suggestion = match field.as_str() {
                "X" | "expression" => {
                    "The expression matrix is required for conversion.\n\
                     In Python: adata.X should contain the expression data.\n\
                     In R: seurat_obj@assays$RNA@data should contain the data."
                }
                "obs" | "cell_metadata" => {
                    "Cell metadata is required.\n\
                     In Python: adata.obs should be a DataFrame.\n\
                     In R: seurat_obj@meta.data should contain cell metadata."
                }
                "var" | "gene_metadata" => {
                    "Gene metadata is required.\n\
                     In Python: adata.var should be a DataFrame.\n\
                     In R: seurat_obj@assays$RNA@meta.features should contain gene info."
                }
                _ => {
                    "Check that your data contains all required fields.\n\
                     Run: crosscell inspect -i <file> --detailed to see what's missing."
                }
            };
            Some(format!("Solution: {}", suggestion))
        }
        
        CrossCellError::DimensionMismatch { context } => {
            Some(format!(
                "Solution: Ensure all data components have consistent dimensions.\n\
                 Issue: {}\n\
                 \n\
                 Common fixes:\n\
                 - Check that obs has the same number of rows as X\n\
                 - Check that var has the same number of rows as X columns\n\
                 - Verify embeddings have the correct number of cells\n\
                 \n\
                 Run: crosscell inspect -i <file> --detailed to check dimensions.",
                context
            ))
        }
        
        CrossCellError::OutOfMemory { context } => {
            Some(format!(
                "Solution: Your dataset is too large for available memory.\n\
                 Context: {}\n\
                 \n\
                 Try these options:\n\
                 1. Use lazy loading: crosscell convert -i input -o output -f seurat --lazy\n\
                 2. Increase available memory\n\
                 3. Process a subset of the data first\n\
                 4. Use a machine with more RAM",
                context
            ))
        }
        
        CrossCellError::ConversionFailed { from, to, reason } => {
            Some(format!(
                "Solution: The conversion from {} to {} failed.\n\
                 Reason: {}\n\
                 \n\
                 Try these steps:\n\
                 1. Run: crosscell inspect -i <input_file> to check for issues\n\
                 2. CrossCell will auto-fix common issues during conversion\n\
                 3. If issues persist, check the detailed error message above",
                from, to, reason
            ))
        }
        
        CrossCellError::UnsupportedType { type_name } => {
            Some(format!(
                "Solution: The data type '{}' is not supported.\n\
                 \n\
                 CrossCell supports:\n\
                 - Numeric types: int32, int64, float32, float64\n\
                 - String types: UTF-8 strings\n\
                 - Boolean types\n\
                 - Categorical/Factor types\n\
                 \n\
                 Try converting the problematic column to a supported type before conversion.",
                type_name
            ))
        }
        
        CrossCellError::InvalidSparseMatrix(msg) => {
            Some(format!(
                "Solution: The sparse matrix structure is invalid.\n\
                 Issue: {}\n\
                 \n\
                 This usually means the sparse matrix indices are corrupted.\n\
                 Try re-creating the sparse matrix in Python:\n\
                   import scipy.sparse as sp\n\
                   adata.X = sp.csr_matrix(adata.X)\n\
                   adata.write('fixed.h5ad')",
                msg
            ))
        }
        
        CrossCellError::ValidationFailed(msg) => {
            Some(format!(
                "Solution: Data validation failed.\n\
                 Issue: {}\n\
                 \n\
                 Run: crosscell inspect -i <file> --detailed\n\
                 This will show you exactly what needs to be fixed.",
                msg
            ))
        }
        
        _ => None,
    }
}

/// Error context for enhanced error reporting
#[derive(Debug, Clone)]
pub struct ErrorContext {
    pub input_file: Option<String>,
    pub output_file: Option<String>,
    pub operation: String,
}

impl ErrorContext {
    pub fn new(operation: &str) -> Self {
        Self {
            input_file: None,
            output_file: None,
            operation: operation.to_string(),
        }
    }
    
    pub fn with_input(mut self, path: &str) -> Self {
        self.input_file = Some(path.to_string());
        self
    }
    
    pub fn with_output(mut self, path: &str) -> Self {
        self.output_file = Some(path.to_string());
        self
    }
}

/// Format a complete error report with context
pub fn format_error_report(error: &CrossCellError, context: &ErrorContext) -> String {
    let mut report = String::new();
    
    report.push_str("❌ Error during ");
    report.push_str(&context.operation);
    report.push_str("\n\n");
    
    if let Some(ref input) = context.input_file {
        report.push_str(&format!("Input file: {}\n", input));
    }
    if let Some(ref output) = context.output_file {
        report.push_str(&format!("Output file: {}\n", output));
    }
    
    report.push_str("\n");
    report.push_str(&format_error_with_suggestions(error));
    
    report.push_str("\n\n");
    report.push_str("For more help:\n");
    report.push_str("  - Run: crosscell inspect -i <file> --detailed\n");
    report.push_str("  - Check documentation: https://github.com/your-repo/crosscell\n");
    
    report
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_file_not_found_suggestion() {
        let error = CrossCellError::file_not_found("test.h5ad");
        let formatted = format_error_with_suggestions(&error);
        assert!(formatted.contains("Solution"));
        assert!(formatted.contains("test.h5ad"));
    }
    
    #[test]
    fn test_missing_field_suggestion() {
        let error = CrossCellError::missing_field("X");
        let formatted = format_error_with_suggestions(&error);
        assert!(formatted.contains("Solution"));
        assert!(formatted.contains("expression"));
    }
    
    #[test]
    fn test_error_context() {
        let context = ErrorContext::new("conversion")
            .with_input("input.h5ad")
            .with_output("output.rds");
        
        assert_eq!(context.operation, "conversion");
        assert_eq!(context.input_file, Some("input.h5ad".to_string()));
        assert_eq!(context.output_file, Some("output.rds".to_string()));
    }
    
    #[test]
    fn test_error_report() {
        let error = CrossCellError::file_not_found("test.h5ad");
        let context = ErrorContext::new("conversion")
            .with_input("test.h5ad");
        
        let report = format_error_report(&error, &context);
        assert!(report.contains("Error during conversion"));
        assert!(report.contains("test.h5ad"));
        assert!(report.contains("Solution"));
    }
}
