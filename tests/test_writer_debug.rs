use crosscell::rds::{IntegerVector, StringVector, GenericVector, StringEncoding};
use crosscell::{read_rds, write_rds, RObject};
use std::path::Path;

#[test]
fn test_write_read_integer() {
    let obj = RObject::IntegerVector(IntegerVector { data: vec![1, 2, 3], ..Default::default() });
    let path = Path::new("tests/data/test_write_int.rds");

    write_rds(path, &obj).expect("Failed to write");
    assert!(path.exists());

    let read_obj = read_rds(path).expect("Failed to read");
    match read_obj {
        RObject::IntegerVector(iv) => {
            assert_eq!(iv.data, vec![1, 2, 3]);
        }
        _ => panic!("Expected IntegerVector, got {:?}", read_obj.type_name()),
    }
}

#[test]
fn test_write_read_string() {
    let mut sv = StringVector::default();
    sv.add("hello".to_string(), StringEncoding::Utf8);
    sv.add("world".to_string(), StringEncoding::Utf8);
    let obj = RObject::StringVector(sv);
    let path = Path::new("tests/data/test_write_str.rds");

    write_rds(path, &obj).expect("Failed to write");
    let read_obj = read_rds(path).expect("Failed to read");

    match read_obj {
        RObject::StringVector(sv) => {
            assert_eq!(sv.data, vec!["hello".to_string(), "world".to_string()]);
        }
        _ => panic!("Expected StringVector, got {:?}", read_obj.type_name()),
    }
}

#[test]
fn test_write_read_list() {
    let mut sv = StringVector::default();
    sv.add("test".to_string(), StringEncoding::Utf8);

    let gv = GenericVector {
        data: vec![
            RObject::IntegerVector(IntegerVector { data: vec![1, 2], ..Default::default() }),
            RObject::StringVector(sv),
        ],
        ..Default::default()
    };
    let obj = RObject::GenericVector(gv);
    let path = Path::new("tests/data/test_write_list.rds");

    write_rds(path, &obj).expect("Failed to write");
    let read_obj = read_rds(path).expect("Failed to read");

    match read_obj {
        RObject::GenericVector(gv) => {
            assert_eq!(gv.data.len(), 2);
        }
        _ => panic!("Expected GenericVector, got {:?}", read_obj.type_name()),
    }
}
