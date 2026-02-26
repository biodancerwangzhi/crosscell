use crosscell::seurat::seurat_rds_to_ir;

#[test]
fn test_seurat_to_ir_minimal() {
    println!("\n========================================");
    println!("Testing seurat_rds_to_ir");
    println!("========================================\n");

    match seurat_rds_to_ir("tests/data/seurat_minimal_simplified.rds") {
        Ok(data) => {
            println!("✅ Successfully parsed Seurat to IR");
            println!("Cells: {}", data.metadata.n_cells);
            println!("Genes: {}", data.metadata.n_genes);
        }
        Err(e) => {
            eprintln!("❌ Failed to parse: {:?}", e);
            panic!("Seurat to IR conversion failed");
        }
    }
}
