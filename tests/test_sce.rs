//! SingleCellExperiment (SCE) format tests
//!
//! Tests for reading and writing SCE objects

use crosscell::sce::{ir_to_sce_rds, sce_rds_to_ir, write_sce_rds};
use std::path::Path;

/// Test reading minimal SCE object
#[test]
fn test_read_sce_minimal() {
    let path = "tests/data/sce_minimal.rds";
    if !Path::new(path).exists() {
        eprintln!("Skipping test: {} not found", path);
        return;
    }

    let result = sce_rds_to_ir(path);
    match result {
        Ok(ir) => {
            println!("Successfully read SCE minimal:");
            println!("  Cells: {}", ir.metadata.n_cells);
            println!("  Genes: {}", ir.metadata.n_genes);
            assert!(ir.metadata.n_cells > 0, "Should have cells");
            assert!(ir.metadata.n_genes > 0, "Should have genes");
        }
        Err(e) => {
            panic!("Failed to read SCE minimal: {}", e);
        }
    }
}

/// Test reading SCE with metadata
#[test]
fn test_read_sce_with_metadata() {
    let path = "tests/data/sce_with_metadata.rds";
    if !Path::new(path).exists() {
        eprintln!("Skipping test: {} not found", path);
        return;
    }

    let result = sce_rds_to_ir(path);
    match result {
        Ok(ir) => {
            println!("Successfully read SCE with metadata:");
            println!("  Cells: {}", ir.metadata.n_cells);
            println!("  Genes: {}", ir.metadata.n_genes);
            println!("  Cell metadata columns: {}", ir.cell_metadata.n_cols());
            println!("  Gene metadata columns: {}", ir.gene_metadata.n_cols());

            assert!(ir.metadata.n_cells > 0, "Should have cells");
            assert!(ir.metadata.n_genes > 0, "Should have genes");
            // SCE with metadata should have cell metadata columns
            assert!(ir.cell_metadata.n_cols() > 0, "Should have cell metadata");
        }
        Err(e) => {
            panic!("Failed to read SCE with metadata: {}", e);
        }
    }
}

/// Test reading SCE with reductions
#[test]
fn test_read_sce_with_reductions() {
    let path = "tests/data/sce_with_reductions.rds";
    if !Path::new(path).exists() {
        eprintln!("Skipping test: {} not found", path);
        return;
    }

    let result = sce_rds_to_ir(path);
    match result {
        Ok(ir) => {
            println!("Successfully read SCE with reductions:");
            println!("  Cells: {}", ir.metadata.n_cells);
            println!("  Genes: {}", ir.metadata.n_genes);

            if let Some(ref embs) = ir.embeddings {
                println!("  Embeddings: {:?}", embs.keys().collect::<Vec<_>>());
                assert!(!embs.is_empty(), "Should have embeddings");
            } else {
                println!("  No embeddings found");
            }

            assert!(ir.metadata.n_cells > 0, "Should have cells");
            assert!(ir.metadata.n_genes > 0, "Should have genes");
        }
        Err(e) => {
            panic!("Failed to read SCE with reductions: {}", e);
        }
    }
}

/// Test reading full SCE object
#[test]
fn test_read_sce_full() {
    let path = "tests/data/sce_full.rds";
    if !Path::new(path).exists() {
        eprintln!("Skipping test: {} not found", path);
        return;
    }

    let result = sce_rds_to_ir(path);
    match result {
        Ok(ir) => {
            println!("Successfully read full SCE:");
            println!("  Cells: {}", ir.metadata.n_cells);
            println!("  Genes: {}", ir.metadata.n_genes);
            println!("  Cell metadata columns: {}", ir.cell_metadata.n_cols());
            println!("  Gene metadata columns: {}", ir.gene_metadata.n_cols());

            if let Some(ref layers) = ir.layers {
                println!("  Layers: {:?}", layers.keys().collect::<Vec<_>>());
            }

            if let Some(ref embs) = ir.embeddings {
                println!("  Embeddings: {:?}", embs.keys().collect::<Vec<_>>());
            }

            assert!(ir.metadata.n_cells > 0, "Should have cells");
            assert!(ir.metadata.n_genes > 0, "Should have genes");
        }
        Err(e) => {
            panic!("Failed to read full SCE: {}", e);
        }
    }
}

/// Test SCE roundtrip (SCE → IR → SCE)
#[test]
fn test_sce_roundtrip() {
    let input_path = "tests/data/sce_minimal.rds";
    let output_path = "tests/data/sce_roundtrip_test.rds";

    if !Path::new(input_path).exists() {
        eprintln!("Skipping test: {} not found", input_path);
        return;
    }

    // Read SCE
    let ir = match sce_rds_to_ir(input_path) {
        Ok(ir) => ir,
        Err(e) => {
            panic!("Failed to read SCE: {}", e);
        }
    };

    println!("Original SCE:");
    println!("  Cells: {}", ir.metadata.n_cells);
    println!("  Genes: {}", ir.metadata.n_genes);

    // Write back to SCE
    match write_sce_rds(&ir, output_path) {
        Ok(_) => {
            println!("Successfully wrote SCE to {}", output_path);
        }
        Err(e) => {
            panic!("Failed to write SCE: {}", e);
        }
    }

    // Read back and verify
    let ir2 = match sce_rds_to_ir(output_path) {
        Ok(ir) => ir,
        Err(e) => {
            panic!("Failed to read roundtrip SCE: {}", e);
        }
    };

    println!("Roundtrip SCE:");
    println!("  Cells: {}", ir2.metadata.n_cells);
    println!("  Genes: {}", ir2.metadata.n_genes);

    // Verify dimensions match
    assert_eq!(
        ir.metadata.n_cells, ir2.metadata.n_cells,
        "Cell count should match"
    );
    assert_eq!(
        ir.metadata.n_genes, ir2.metadata.n_genes,
        "Gene count should match"
    );

    // Clean up
    let _ = std::fs::remove_file(output_path);
}

/// Test SCE roundtrip with metadata
#[test]
fn test_sce_roundtrip_with_metadata() {
    let input_path = "tests/data/sce_with_metadata.rds";
    let output_path = "tests/data/sce_roundtrip_metadata_test.rds";

    if !Path::new(input_path).exists() {
        eprintln!("Skipping test: {} not found", input_path);
        return;
    }

    // Read SCE
    let ir = match sce_rds_to_ir(input_path) {
        Ok(ir) => ir,
        Err(e) => {
            panic!("Failed to read SCE: {}", e);
        }
    };

    println!("Original SCE with metadata:");
    println!("  Cells: {}", ir.metadata.n_cells);
    println!("  Genes: {}", ir.metadata.n_genes);
    println!("  Cell metadata columns: {}", ir.cell_metadata.n_cols());

    // Write back to SCE
    match write_sce_rds(&ir, output_path) {
        Ok(_) => {
            println!("Successfully wrote SCE to {}", output_path);
        }
        Err(e) => {
            panic!("Failed to write SCE: {}", e);
        }
    }

    // Read back and verify
    let ir2 = match sce_rds_to_ir(output_path) {
        Ok(ir) => ir,
        Err(e) => {
            panic!("Failed to read roundtrip SCE: {}", e);
        }
    };

    println!("Roundtrip SCE:");
    println!("  Cells: {}", ir2.metadata.n_cells);
    println!("  Genes: {}", ir2.metadata.n_genes);
    println!("  Cell metadata columns: {}", ir2.cell_metadata.n_cols());

    // Verify dimensions match
    assert_eq!(
        ir.metadata.n_cells, ir2.metadata.n_cells,
        "Cell count should match"
    );
    assert_eq!(
        ir.metadata.n_genes, ir2.metadata.n_genes,
        "Gene count should match"
    );

    // Clean up
    let _ = std::fs::remove_file(output_path);
}

/// Test SCE roundtrip with embeddings
#[test]
fn test_sce_roundtrip_with_embeddings() {
    let input_path = "tests/data/sce_with_reductions.rds";
    let output_path = "tests/data/sce_roundtrip_embeddings_test.rds";

    if !Path::new(input_path).exists() {
        eprintln!("Skipping test: {} not found", input_path);
        return;
    }

    // Read SCE
    let ir = match sce_rds_to_ir(input_path) {
        Ok(ir) => ir,
        Err(e) => {
            panic!("Failed to read SCE: {}", e);
        }
    };

    println!("Original SCE with embeddings:");
    println!("  Cells: {}", ir.metadata.n_cells);
    println!("  Genes: {}", ir.metadata.n_genes);
    if let Some(ref embs) = ir.embeddings {
        println!("  Embeddings: {:?}", embs.keys().collect::<Vec<_>>());
    }

    // Write back to SCE
    match write_sce_rds(&ir, output_path) {
        Ok(_) => {
            println!("Successfully wrote SCE to {}", output_path);
        }
        Err(e) => {
            panic!("Failed to write SCE: {}", e);
        }
    }

    // Read back and verify
    let ir2 = match sce_rds_to_ir(output_path) {
        Ok(ir) => ir,
        Err(e) => {
            panic!("Failed to read roundtrip SCE: {}", e);
        }
    };

    println!("Roundtrip SCE:");
    println!("  Cells: {}", ir2.metadata.n_cells);
    println!("  Genes: {}", ir2.metadata.n_genes);
    if let Some(ref embs) = ir2.embeddings {
        println!("  Embeddings: {:?}", embs.keys().collect::<Vec<_>>());
    }

    // Verify dimensions match
    assert_eq!(
        ir.metadata.n_cells, ir2.metadata.n_cells,
        "Cell count should match"
    );
    assert_eq!(
        ir.metadata.n_genes, ir2.metadata.n_genes,
        "Gene count should match"
    );

    // Verify embeddings exist
    if let Some(ref orig_embs) = ir.embeddings {
        let rt_embs = ir2
            .embeddings
            .as_ref()
            .expect("Roundtrip should have embeddings");
        assert_eq!(
            orig_embs.len(),
            rt_embs.len(),
            "Embedding count should match"
        );
    }

    // Clean up
    let _ = std::fs::remove_file(output_path);
}

/// Test full SCE roundtrip (with layers and embeddings)
#[test]
fn test_sce_roundtrip_full() {
    let input_path = "tests/data/sce_full.rds";
    let output_path = "tests/data/sce_roundtrip_full_test.rds";

    if !Path::new(input_path).exists() {
        eprintln!("Skipping test: {} not found", input_path);
        return;
    }

    // Read SCE
    let ir = match sce_rds_to_ir(input_path) {
        Ok(ir) => ir,
        Err(e) => {
            panic!("Failed to read SCE: {}", e);
        }
    };

    println!("Original full SCE:");
    println!("  Cells: {}", ir.metadata.n_cells);
    println!("  Genes: {}", ir.metadata.n_genes);
    println!("  Cell metadata columns: {}", ir.cell_metadata.n_cols());
    if let Some(ref layers) = ir.layers {
        println!("  Layers: {:?}", layers.keys().collect::<Vec<_>>());
    }
    if let Some(ref embs) = ir.embeddings {
        println!("  Embeddings: {:?}", embs.keys().collect::<Vec<_>>());
    }

    // Write back to SCE
    match write_sce_rds(&ir, output_path) {
        Ok(_) => {
            println!("Successfully wrote SCE to {}", output_path);
        }
        Err(e) => {
            panic!("Failed to write SCE: {}", e);
        }
    }

    // Read back and verify
    let ir2 = match sce_rds_to_ir(output_path) {
        Ok(ir) => ir,
        Err(e) => {
            panic!("Failed to read roundtrip SCE: {}", e);
        }
    };

    println!("Roundtrip full SCE:");
    println!("  Cells: {}", ir2.metadata.n_cells);
    println!("  Genes: {}", ir2.metadata.n_genes);
    println!("  Cell metadata columns: {}", ir2.cell_metadata.n_cols());
    if let Some(ref layers) = ir2.layers {
        println!("  Layers: {:?}", layers.keys().collect::<Vec<_>>());
    }
    if let Some(ref embs) = ir2.embeddings {
        println!("  Embeddings: {:?}", embs.keys().collect::<Vec<_>>());
    }

    // Verify dimensions match
    assert_eq!(
        ir.metadata.n_cells, ir2.metadata.n_cells,
        "Cell count should match"
    );
    assert_eq!(
        ir.metadata.n_genes, ir2.metadata.n_genes,
        "Gene count should match"
    );

    // Verify layers exist
    if let Some(ref orig_layers) = ir.layers {
        let rt_layers = ir2.layers.as_ref().expect("Roundtrip should have layers");
        assert_eq!(
            orig_layers.len(),
            rt_layers.len(),
            "Layer count should match"
        );
    }

    // Verify embeddings exist
    if let Some(ref orig_embs) = ir.embeddings {
        let rt_embs = ir2
            .embeddings
            .as_ref()
            .expect("Roundtrip should have embeddings");
        assert_eq!(
            orig_embs.len(),
            rt_embs.len(),
            "Embedding count should match"
        );
    }

    // Clean up
    let _ = std::fs::remove_file(output_path);
}
