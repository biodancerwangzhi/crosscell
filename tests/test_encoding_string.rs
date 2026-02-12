//! Test encoding string handling in RDS version 2 format

use crosscell::rds::IntegerVector;
use crosscell::{read_rds, write_rds, RObject};
use std::path::Path;

#[test]
fn test_write_and_read_with_encoding_string() {
    let obj = RObject::IntegerVector(IntegerVector { data: vec![1, 2, 3, 4, 5], ..Default::default() });
    let temp_path = Path::new("tests/data/test_encoding_string.rds");

    write_rds(temp_path, &obj).expect("Failed to write RDS");
    let read_obj = read_rds(temp_path).expect("Failed to read RDS");

    match read_obj {
        RObject::IntegerVector(iv) => {
            assert_eq!(iv.data, vec![1, 2, 3, 4, 5]);
        }
        _ => panic!("Expected IntegerVector, got {:?}", read_obj.type_name()),
    }
    let _ = std::fs::remove_file(temp_path);
}
