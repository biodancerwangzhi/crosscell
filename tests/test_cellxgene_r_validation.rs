//! Integration test: Convert CellxGENE datasets (larger, real-world) to Seurat RDS
//!
//! Run: docker-compose run --rm dev cargo test --test test_cellxgene_r_validation -- --nocapture --ignored
//! Then: docker-compose run --rm dev Rscript tests/test_multi_dataset_r_validation.R

use crosscell::read_h5ad;
use crosscell::seurat::write_seurat_rds;
use std::path::Path;

/// Large CellxGENE datasets - marked #[ignore] since they take longer
const CELLXGENE_DATASETS: &[(&str, &str)] = &[
    (
        "data/generated/cellxgene_pbmc_15k.h5ad",
        "tests/data/multi_val/cellxgene_pbmc_15k.rds",
    ),
    (
        "data/generated/cellxgene_heart_23k.h5ad",
        "tests/data/multi_val/cellxgene_heart_23k.rds",
    ),
    (
        "data/generated/cellxgene_brain_40k.h5ad",
        "tests/data/multi_val/cellxgene_brain_40k.rds",
    ),
];

#[test]
#[ignore]
fn test_convert_cellxgene_to_rds() {
    std::fs::create_dir_all("tests/data/multi_val").ok();

    let mut success = 0;
    let mut skipped = 0;
    let mut failed = 0;

    for (h5ad_path, rds_path) in CELLXGENE_DATASETS {
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

    println!("\n=== CellxGENE Summary ===");
    println!(
        "Success: {}, Skipped: {}, Failed: {}",
        success, skipped, failed
    );
    assert_eq!(failed, 0, "{} CellxGENE datasets failed conversion", failed);
}
