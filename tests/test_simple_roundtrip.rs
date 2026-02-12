//! Simple RDS roundtrip tests

use crosscell::rds::{IntegerVector, DoubleVector, StringVector, GenericVector, StringEncoding};
use crosscell::{read_rds, write_rds, RObject};
use std::path::Path;

fn make_string_vec(strings: &[&str]) -> RObject {
    let mut sv = StringVector::default();
    for s in strings {
        sv.add(s.to_string(), StringEncoding::Utf8);
    }
    RObject::StringVector(sv)
}

fn make_named_list(items: Vec<(&str, RObject)>) -> RObject {
    let mut gv = GenericVector::default();
    let mut names_sv = StringVector::default();
    for (name, val) in items {
        gv.data.push(val);
        names_sv.add(name.to_string(), StringEncoding::Utf8);
    }
    gv.attributes.add("names".to_string(), RObject::StringVector(names_sv), StringEncoding::Utf8);
    RObject::GenericVector(gv)
}

#[test]
fn test_simple_integer_roundtrip() {
    let original = RObject::IntegerVector(IntegerVector { data: vec![1, 2, 3, 4, 5], ..Default::default() });
    let temp_path = Path::new("tests/data/temp_simple_int.rds");
    write_rds(temp_path, &original).expect("Failed to write");
    let roundtrip = read_rds(temp_path).expect("Failed to read");
    match (&original, &roundtrip) {
        (RObject::IntegerVector(orig), RObject::IntegerVector(rt)) => {
            assert_eq!(orig.data, rt.data);
        }
        _ => panic!("Type mismatch"),
    }
    let _ = std::fs::remove_file(temp_path);
}

#[test]
fn test_simple_real_roundtrip() {
    let original = RObject::DoubleVector(DoubleVector { data: vec![1.1, 2.2, 3.3, 4.4, 5.5], ..Default::default() });
    let temp_path = Path::new("tests/data/temp_simple_real.rds");
    write_rds(temp_path, &original).expect("Failed to write");
    let roundtrip = read_rds(temp_path).expect("Failed to read");
    match (&original, &roundtrip) {
        (RObject::DoubleVector(orig), RObject::DoubleVector(rt)) => {
            assert_eq!(orig.data, rt.data);
        }
        _ => panic!("Type mismatch"),
    }
    let _ = std::fs::remove_file(temp_path);
}

#[test]
fn test_simple_string_roundtrip() {
    let original = make_string_vec(&["hello", "world", "test"]);
    let temp_path = Path::new("tests/data/temp_simple_string.rds");
    write_rds(temp_path, &original).expect("Failed to write");
    let roundtrip = read_rds(temp_path).expect("Failed to read");
    match (&original, &roundtrip) {
        (RObject::StringVector(orig), RObject::StringVector(rt)) => {
            assert_eq!(orig.data, rt.data);
        }
        _ => panic!("Type mismatch"),
    }
    let _ = std::fs::remove_file(temp_path);
}

#[test]
fn test_named_list_roundtrip() {
    let original = make_named_list(vec![
        ("a", RObject::IntegerVector(IntegerVector { data: vec![1, 2, 3], ..Default::default() })),
        ("b", RObject::DoubleVector(DoubleVector { data: vec![1.1, 2.2], ..Default::default() })),
        ("c", make_string_vec(&["test"])),
    ]);
    let temp_path = Path::new("tests/data/temp_named_list.rds");
    write_rds(temp_path, &original).expect("Failed to write");
    let roundtrip = read_rds(temp_path).expect("Failed to read");
    match (&original, &roundtrip) {
        (RObject::GenericVector(orig), RObject::GenericVector(rt)) => {
            assert_eq!(orig.data.len(), rt.data.len());
        }
        _ => panic!("Type mismatch: expected GenericVector"),
    }
    let _ = std::fs::remove_file(temp_path);
}
