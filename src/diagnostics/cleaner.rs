//! Data cleaning module
//!
//! Provides automatic data cleaning and fixing capabilities.

use crate::ir::{SingleCellData, DataFrame};
use crate::diagnostics::detector::sanitize_column_name;
use crate::diagnostics::audit::{AuditLog, Modification};
use std::collections::HashMap;

/// Options for data cleaning
#[derive(Debug, Clone)]
pub struct CleaningOptions {
    /// Sanitize column names (replace special characters)
    pub sanitize_names: bool,
    /// Handle duplicate column names
    pub handle_duplicates: bool,
    /// Drop empty columns
    pub drop_empty_columns: bool,
    /// Truncate long column names
    pub truncate_long_names: bool,
    /// Maximum column name length
    pub max_name_length: usize,
}

impl Default for CleaningOptions {
    fn default() -> Self {
        Self {
            sanitize_names: true,
            handle_duplicates: true,
            drop_empty_columns: true,
            truncate_long_names: true,
            max_name_length: 256,
        }
    }
}

/// Report of cleaning operations performed
#[derive(Debug, Clone)]
pub struct CleaningReport {
    pub columns_renamed: Vec<(String, String)>,
    pub columns_dropped: Vec<String>,
    pub duplicates_resolved: Vec<(String, String)>,
    pub names_truncated: Vec<(String, String)>,
    pub total_modifications: usize,
}

impl CleaningReport {
    pub fn new() -> Self {
        Self {
            columns_renamed: Vec::new(),
            columns_dropped: Vec::new(),
            duplicates_resolved: Vec::new(),
            names_truncated: Vec::new(),
            total_modifications: 0,
        }
    }

    pub fn has_modifications(&self) -> bool {
        self.total_modifications > 0
    }

    /// Convert to audit log
    pub fn to_audit_log(&self) -> AuditLog {
        let mut log = AuditLog::new();
        
        for (from, to) in &self.columns_renamed {
            log.add_modification(Modification::ColumnRenamed {
                from: from.clone(),
                to: to.clone(),
            });
        }
        
        for name in &self.columns_dropped {
            log.add_modification(Modification::ColumnDropped {
                name: name.clone(),
                reason: "Empty column".to_string(),
            });
        }
        
        for (from, to) in &self.duplicates_resolved {
            log.add_modification(Modification::ColumnRenamed {
                from: from.clone(),
                to: to.clone(),
            });
        }
        
        for (from, to) in &self.names_truncated {
            log.add_modification(Modification::ColumnRenamed {
                from: from.clone(),
                to: to.clone(),
            });
        }
        
        log
    }
}

impl Default for CleaningReport {
    fn default() -> Self {
        Self::new()
    }
}

/// Data cleaner for automatic issue fixing
pub struct DataCleaner {
    options: CleaningOptions,
}

impl DataCleaner {
    pub fn new(options: CleaningOptions) -> Self {
        Self { options }
    }

    pub fn with_defaults() -> Self {
        Self::new(CleaningOptions::default())
    }

    /// Clean single-cell data and return a report
    pub fn clean(&self, data: &mut SingleCellData) -> CleaningReport {
        let mut report = CleaningReport::new();
        
        // Clean cell metadata
        self.clean_dataframe(&mut data.cell_metadata, "obs", &mut report);
        
        // Clean gene metadata
        self.clean_dataframe(&mut data.gene_metadata, "var", &mut report);
        
        report
    }

    /// Clean a DataFrame
    fn clean_dataframe(&self, df: &mut DataFrame, prefix: &str, report: &mut CleaningReport) {
        // Track name mappings for duplicate handling
        let mut name_counts: HashMap<String, usize> = HashMap::new();
        let mut new_columns = Vec::new();
        
        for col_name in df.columns.iter() {
            let original_name = col_name.clone();
            let mut new_name = original_name.clone();
            let full_original = format!("{}.{}", prefix, original_name);
            
            // Step 1: Sanitize special characters
            if self.options.sanitize_names {
                let sanitized = sanitize_column_name(&new_name);
                if sanitized != new_name {
                    report.columns_renamed.push((full_original.clone(), format!("{}.{}", prefix, sanitized)));
                    report.total_modifications += 1;
                    new_name = sanitized;
                }
            }
            
            // Step 2: Truncate long names
            if self.options.truncate_long_names && new_name.len() > self.options.max_name_length {
                let truncated = new_name[..self.options.max_name_length].to_string();
                report.names_truncated.push((format!("{}.{}", prefix, new_name), format!("{}.{}", prefix, truncated)));
                report.total_modifications += 1;
                new_name = truncated;
            }
            
            // Step 3: Handle duplicates
            if self.options.handle_duplicates {
                let count = name_counts.entry(new_name.clone()).or_insert(0);
                *count += 1;
                
                if *count > 1 {
                    let deduped = format!("{}_dup{}", new_name, count);
                    report.duplicates_resolved.push((format!("{}.{}", prefix, new_name), format!("{}.{}", prefix, deduped)));
                    report.total_modifications += 1;
                    new_name = deduped;
                }
            }
            
            new_columns.push(new_name);
        }
        
        df.columns = new_columns;
    }

    /// Preview cleaning without modifying data
    pub fn preview(&self, data: &SingleCellData) -> CleaningReport {
        let mut data_clone = data.clone();
        self.clean(&mut data_clone)
    }
}

impl Default for DataCleaner {
    fn default() -> Self {
        Self::with_defaults()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::{ExpressionMatrix, DenseMatrix, DatasetMetadata};
    use arrow::array::StringArray;
    use std::sync::Arc;

    fn create_test_data() -> SingleCellData {
        let expression = ExpressionMatrix::Dense(DenseMatrix {
            data: vec![1.0, 2.0, 3.0, 4.0],
            n_rows: 2,
            n_cols: 2,
        });
        
        // DataFrame with columns containing special characters
        let cell_metadata = DataFrame {
            columns: vec!["cell/type".to_string(), "batch".to_string()],
            data: vec![
                Arc::new(StringArray::from(vec!["A", "B"])) as arrow::array::ArrayRef,
                Arc::new(StringArray::from(vec!["1", "2"])) as arrow::array::ArrayRef,
            ],
            n_rows: 2,
        };
        
        let gene_metadata = DataFrame {
            columns: vec!["gene_name".to_string()],
            data: vec![
                Arc::new(StringArray::from(vec!["geneA", "geneB"])) as arrow::array::ArrayRef,
            ],
            n_rows: 2,
        };
        
        let metadata = DatasetMetadata::new(2, 2, "test".to_string());
        
        SingleCellData {
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
        }
    }

    #[test]
    fn test_clean_special_characters() {
        let mut data = create_test_data();
        let cleaner = DataCleaner::with_defaults();
        let report = cleaner.clean(&mut data);
        
        assert!(report.has_modifications());
        assert_eq!(report.columns_renamed.len(), 1);
        assert!(report.columns_renamed[0].0.contains("cell/type"));
        
        // Verify the column was renamed
        assert_eq!(data.cell_metadata.columns[0], "cell_type");
    }

    #[test]
    fn test_preview_does_not_modify() {
        let data = create_test_data();
        let cleaner = DataCleaner::with_defaults();
        let report = cleaner.preview(&data);
        
        assert!(report.has_modifications());
        // Original data should be unchanged
        assert_eq!(data.cell_metadata.columns[0], "cell/type");
    }

    #[test]
    fn test_handle_duplicates() {
        let expression = ExpressionMatrix::Dense(DenseMatrix {
            data: vec![1.0, 2.0, 3.0, 4.0],
            n_rows: 2,
            n_cols: 2,
        });
        
        // DataFrame with duplicate column names
        let cell_metadata = DataFrame {
            columns: vec!["batch".to_string(), "batch".to_string()],
            data: vec![
                Arc::new(StringArray::from(vec!["A", "B"])) as arrow::array::ArrayRef,
                Arc::new(StringArray::from(vec!["1", "2"])) as arrow::array::ArrayRef,
            ],
            n_rows: 2,
        };
        
        let gene_metadata = DataFrame::empty(2);
        let metadata = DatasetMetadata::new(2, 2, "test".to_string());
        
        let mut data = SingleCellData {
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
        
        let cleaner = DataCleaner::with_defaults();
        let report = cleaner.clean(&mut data);
        
        assert!(report.has_modifications());
        assert_eq!(report.duplicates_resolved.len(), 1);
        
        // Verify duplicates were resolved
        assert_eq!(data.cell_metadata.columns[0], "batch");
        assert_eq!(data.cell_metadata.columns[1], "batch_dup2");
    }
}
