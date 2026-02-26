//! Test writing dgCMatrix objects

use crosscell::rds::{DoubleVector, GenericVector, IntegerVector, S4Object, StringEncoding};
use crosscell::{read_rds, write_rds, RObject};
use std::path::Path;

#[test]
fn test_write_simple_dgcmatrix() {
    // Create a simple dgCMatrix (3x2 sparse matrix)
    let mut s4 = S4Object::default();
    s4.class_name = "dgCMatrix".to_string();
    s4.class_encoding = StringEncoding::Utf8;
    s4.package_name = "Matrix".to_string();
    s4.package_encoding = StringEncoding::Utf8;

    // Add slots as attributes
    s4.attributes.add(
        "i".to_string(),
        RObject::IntegerVector(IntegerVector {
            data: vec![0, 1, 0],
            ..Default::default()
        }),
        StringEncoding::Utf8,
    );
    s4.attributes.add(
        "p".to_string(),
        RObject::IntegerVector(IntegerVector {
            data: vec![0, 2, 3],
            ..Default::default()
        }),
        StringEncoding::Utf8,
    );
    s4.attributes.add(
        "x".to_string(),
        RObject::DoubleVector(DoubleVector {
            data: vec![1.0, 2.0, 3.0],
            ..Default::default()
        }),
        StringEncoding::Utf8,
    );
    s4.attributes.add(
        "Dim".to_string(),
        RObject::IntegerVector(IntegerVector {
            data: vec![3, 2],
            ..Default::default()
        }),
        StringEncoding::Utf8,
    );
    s4.attributes.add(
        "Dimnames".to_string(),
        RObject::GenericVector(GenericVector {
            data: vec![RObject::Null, RObject::Null],
            ..Default::default()
        }),
        StringEncoding::Utf8,
    );
    s4.attributes.add(
        "factors".to_string(),
        RObject::GenericVector(GenericVector::default()),
        StringEncoding::Utf8,
    );

    let dgc = RObject::S4Object(s4);

    let path = "tests/data/test_dgcmatrix.rds";
    let result = write_rds(Path::new(path), &dgc);
    assert!(
        result.is_ok(),
        "Failed to write dgCMatrix: {:?}",
        result.err()
    );

    let read_result = read_rds(Path::new(path));
    assert!(
        read_result.is_ok(),
        "Failed to read back: {:?}",
        read_result.err()
    );
}
