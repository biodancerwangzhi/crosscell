use crosscell::{read_rds, RObject};
use std::path::Path;

#[test]
fn test_simple_dgcmatrix() {
    let path = Path::new("tests/data/simple_dgcmatrix.rds");

    println!("\n========================================");
    println!("Reading simple dgCMatrix (S4 object)");
    println!("========================================\n");

    match read_rds(path) {
        Ok(robj) => {
            println!("✅ Successfully read RDS object");
            println!("Type: {}", robj.type_name());
            println!("{:#?}", robj);
        }
        Err(e) => {
            eprintln!("❌ Failed to read RDS: {:?}", e);
            panic!("RDS read failed");
        }
    }
}
