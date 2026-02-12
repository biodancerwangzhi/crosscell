//! Test reading R-created RDS files

use crosscell::{read_rds, RObject};
use std::path::Path;

#[test]
fn test_read_r_simple_int() {
    let path = Path::new("tests/data/r_simple_int.rds");
    
    if !path.exists() {
        eprintln!("⚠️  File not found, run: Rscript tests/create_simple_dgc.R");
        return;
    }
    
    let result = read_rds(path);
    match result {
        Ok(obj) => {
            println!("✓ Read R-created integer vector");
            println!("  Type: {}", obj.type_name());
        }
        Err(e) => {
            panic!("Failed to read R-created file: {:?}", e);
        }
    }
}

#[test]
fn test_read_r_dgcmatrix() {
    let path = Path::new("tests/data/r_created_dgc.rds");
    
    if !path.exists() {
        eprintln!("⚠️  File not found, run: Rscript tests/create_simple_dgc.R");
        return;
    }
    
    let result = read_rds(path);
    match result {
        Ok(obj) => {
            println!("✓ Read R-created dgCMatrix");
            println!("  Type: {}", obj.type_name());
        }
        Err(e) => {
            panic!("Failed to read R-created dgCMatrix: {:?}", e);
        }
    }
}
