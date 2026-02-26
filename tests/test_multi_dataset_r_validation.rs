//! Integration test: Convert multiple h5ad datasets to Seurat RDS
//! and generate files for R validation.
//!
//! Run: docker-compose run --rm dev cargo test --test test_multi_dataset_r_validation -- --nocapture
//! Then: docker-compose run --rm dev Rscript tests/test_multi_dataset_r_validation.R

use crosscell::read_h5ad;
use crosscell::seurat::write_seurat_rds;
use std::path::Path;

/// Test datasets to convert
const DATASETS: &[(&str, &str)] = &[
    (
        "data/generated/scanpy_pbmc3k.h5ad",
        "tests/data/multi_val/scanpy_pbmc3k.rds",
    ),
    (
        "data/generated/scanpy_pbmc3k_processed.h5ad",
        "tests/data/multi_val/scanpy_pbmc3k_processed.rds",
    ),
    (
        "data/generated/scvelo_pancreas.h5ad",
        "tests/data/multi_val/scvelo_pancreas.rds",
    ),
    (
        "data/generated/scvelo_dentategyrus.h5ad",
        "tests/data/multi_val/scvelo_dentategyrus.rds",
    ),
    (
        "tests/data/small_sparse.h5ad",
        "tests/data/multi_val/small_sparse.rds",
    ),
    (
        "tests/data/small_dense.h5ad",
        "tests/data/multi_val/small_dense.rds",
    ),
    (
        "tests/data/medium_sparse.h5ad",
        "tests/data/multi_val/medium_sparse.rds",
    ),
    (
        "tests/data/single_column.h5ad",
        "tests/data/multi_val/single_column.rds",
    ),
    (
        "tests/data/single_row.h5ad",
        "tests/data/multi_val/single_row.rds",
    ),
];

#[test]
fn test_convert_multiple_h5ad_to_rds() {
    // Create output directory
    std::fs::create_dir_all("tests/data/multi_val").ok();

    let mut success = 0;
    let mut skipped = 0;
    let mut failed = 0;

    for (h5ad_path, rds_path) in DATASETS {
        if !Path::new(h5ad_path).exists() {
            eprintln!("SKIP: {} not found", h5ad_path);
            skipped += 1;
            continue;
        }

        print!("Converting {} ... ", h5ad_path);

        match read_h5ad(h5ad_path) {
            Ok(ir) => {
                println!(
                    "IR: {} cells x {} genes",
                    ir.metadata.n_cells, ir.metadata.n_genes
                );
                match write_seurat_rds(&ir, rds_path) {
                    Ok(_) => {
                        println!("  -> Wrote {}", rds_path);
                        success += 1;
                    }
                    Err(e) => {
                        eprintln!("  -> WRITE FAILED: {}", e);
                        failed += 1;
                    }
                }
            }
            Err(e) => {
                eprintln!("  -> READ FAILED: {}", e);
                failed += 1;
            }
        }
    }

    println!("\n=== Summary ===");
    println!(
        "Success: {}, Skipped: {}, Failed: {}",
        success, skipped, failed
    );
    assert_eq!(failed, 0, "{} datasets failed conversion", failed);
}
