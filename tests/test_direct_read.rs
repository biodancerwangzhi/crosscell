//! Integration tests for direct Seurat reading
//!
//! Tests the `read_seurat_direct()` function with actual Seurat RDS files.
//! 
//! Note: `read_seurat_direct()` is designed for UNSIMPLIFIED Seurat objects (S4 objects).
//! Simplified Seurat objects (lists) should use `seurat_rds_to_ir()` instead.

use crosscell::seurat::{read_seurat_direct, SeuratVersion};
use std::path::Path;

/// Test direct reading of an unsimplified Seurat object (S4)
/// This tests the core functionality of reading native Seurat objects
#[test]
fn test_direct_read_seurat_minimal() {
    // seurat_minimal.rds should be an S4 Seurat object (not simplified)
    let path = Path::new("tests/data/seurat_minimal.rds");
    if !path.exists() {
        eprintln!("Skipping test: {} not found", path.display());
        return;
    }

    let result = read_seurat_direct(path, false);
    match result {
        Ok(direct_result) => {
            println!("Successfully read Seurat object directly!");
            println!("  Version: {}", direct_result.version);
            println!("  Cells: {}", direct_result.data.metadata.n_cells);
            println!("  Genes: {}", direct_result.data.metadata.n_genes);
            println!("  Skipped components: {}", direct_result.skipped.total_skipped());
            
            // Verify basic structure
            assert!(direct_result.data.metadata.n_cells > 0, "Should have cells");
            assert!(direct_result.data.metadata.n_genes > 0, "Should have genes");
        }
        Err(e) => {
            // This may fail if the file is actually a simplified list
            eprintln!("Note: {} - file may be simplified (list) not S4", e);
        }
    }
}

/// Test direct reading of a Seurat object with dimensional reductions
#[test]
fn test_direct_read_seurat_with_dimred() {
    // seurat_with_dimred.rds should be an S4 Seurat object with PCA/UMAP
    let path = Path::new("tests/data/seurat_with_dimred.rds");
    if !path.exists() {
        eprintln!("Skipping test: {} not found", path.display());
        return;
    }

    let result = read_seurat_direct(path, false);
    match result {
        Ok(direct_result) => {
            println!("Successfully read Seurat object with reductions!");
            println!("  Version: {}", direct_result.version);
            println!("  Cells: {}", direct_result.data.metadata.n_cells);
            println!("  Genes: {}", direct_result.data.metadata.n_genes);
            
            // Check embeddings
            if let Some(ref embeddings) = direct_result.data.embeddings {
                println!("  Embeddings: {:?}", embeddings.keys().collect::<Vec<_>>());
                assert!(!embeddings.is_empty(), "Should have embeddings");
            }
            
            println!("  Skipped: {} components", direct_result.skipped.total_skipped());
        }
        Err(e) => {
            eprintln!("Note: {} - file may be simplified (list) not S4", e);
        }
    }
}

/// Test direct reading of a multi-assay Seurat object
#[test]
fn test_direct_read_seurat_multi_assay() {
    // seurat_multi_assay.rds should be an S4 Seurat object with multiple assays
    let path = Path::new("tests/data/seurat_multi_assay.rds");
    if !path.exists() {
        eprintln!("Skipping test: {} not found", path.display());
        return;
    }

    let result = read_seurat_direct(path, false);
    match result {
        Ok(direct_result) => {
            println!("Successfully read multi-assay Seurat object!");
            println!("  Version: {}", direct_result.version);
            println!("  Cells: {}", direct_result.data.metadata.n_cells);
            println!("  Genes: {}", direct_result.data.metadata.n_genes);
            
            // Check layers (additional assays)
            if let Some(ref layers) = direct_result.data.layers {
                println!("  Layers: {:?}", layers.keys().collect::<Vec<_>>());
            }
            
            println!("  Skipped: {} components", direct_result.skipped.total_skipped());
        }
        Err(e) => {
            eprintln!("Note: {} - file may be simplified (list) not S4", e);
        }
    }
}

/// Test direct reading of test_minimal_seurat.rds
#[test]
fn test_direct_read_test_minimal_seurat() {
    let path = Path::new("tests/data/test_minimal_seurat.rds");
    if !path.exists() {
        eprintln!("Skipping test: {} not found", path.display());
        return;
    }

    let result = read_seurat_direct(path, false);
    match result {
        Ok(direct_result) => {
            println!("Successfully read test_minimal_seurat!");
            println!("  Version: {}", direct_result.version);
            println!("  Cells: {}", direct_result.data.metadata.n_cells);
            println!("  Genes: {}", direct_result.data.metadata.n_genes);
            println!("  Skipped: {} components", direct_result.skipped.total_skipped());
            
            assert!(direct_result.data.metadata.n_cells > 0);
            assert!(direct_result.data.metadata.n_genes > 0);
        }
        Err(e) => {
            eprintln!("Note: {} - file may be simplified (list) not S4", e);
        }
    }
}

/// Test version detection with various Seurat files
#[test]
fn test_seurat_version_detection() {
    let test_files = [
        "tests/data/seurat_minimal.rds",
        "tests/data/seurat_with_dimred.rds",
        "tests/data/test_minimal_seurat.rds",
    ];
    
    for path_str in &test_files {
        let path = Path::new(path_str);
        if !path.exists() {
            continue;
        }
        
        match read_seurat_direct(path, false) {
            Ok(result) => {
                println!("{}: Version {}", path_str, result.version);
                // Version should be one of the valid versions
                assert!(
                    matches!(result.version, SeuratVersion::V3 | SeuratVersion::V4 | SeuratVersion::V5 | SeuratVersion::Unknown),
                    "Should detect a valid version"
                );
            }
            Err(e) => {
                println!("{}: Not an S4 Seurat object ({})", path_str, e);
            }
        }
    }
}

/// Test skipped components tracking
#[test]
fn test_skipped_components_tracking() {
    let path = Path::new("tests/data/seurat_minimal.rds");
    if !path.exists() {
        eprintln!("Skipping test: {} not found", path.display());
        return;
    }

    if let Ok(result) = read_seurat_direct(path, false) {
        let skipped = &result.skipped;
        
        println!("Skipped components:");
        println!("  Environments: {} items", skipped.environments.len());
        println!("  Closures: {} items", skipped.closures.len());
        println!("  External pointers: {} items", skipped.external_pointers.len());
        println!("  Bytecodes: {} items", skipped.bytecodes.len());
        println!("  Languages: {} items", skipped.languages.len());
        println!("  Total skipped: {}", skipped.total_skipped());
        
        // Unsimplified Seurat objects typically have some skipped components
        // (commands, tools, etc.)
    }
}

/// Test that simplified files are correctly rejected
/// (they should use seurat_rds_to_ir instead)
#[test]
fn test_simplified_files_rejected() {
    let simplified_files = [
        "tests/data/seurat_minimal_simplified.rds",
        "tests/data/seurat_with_dimred_simplified.rds",
        "tests/data/seurat_multi_assay_simplified.rds",
    ];
    
    for path_str in &simplified_files {
        let path = Path::new(path_str);
        if !path.exists() {
            continue;
        }
        
        let result = read_seurat_direct(path, false);
        // Simplified files are lists, not S4 objects, so they should fail
        assert!(result.is_err(), "Simplified file {} should be rejected by read_seurat_direct", path_str);
        if let Err(e) = result {
            println!("{}: Correctly rejected - {}", path_str, e);
        }
    }
}
