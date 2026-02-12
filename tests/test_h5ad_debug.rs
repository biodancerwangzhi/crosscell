/// Debug test to check if hdf5 crate can read float32 sparse data from old-format h5ad files
#[test]
#[ignore]
fn test_read_old_format_h5ad_sparse() {
    let path = "data/generated/scanpy_pbmc3k.h5ad";
    if !std::path::Path::new(path).exists() {
        eprintln!("Skipping: {} not found", path);
        return;
    }
    
    let file = hdf5::File::open(path).unwrap();
    
    // Check X group
    let x_group = file.group("X").unwrap();
    println!("X group opened successfully");
    
    // Check attributes
    let members = x_group.member_names().unwrap();
    println!("X members: {:?}", members);
    
    // Try reading shape attribute
    if let Ok(attr) = x_group.attr("h5sparse_shape") {
        let shape: Vec<usize> = attr.read_1d().unwrap().to_vec();
        println!("h5sparse_shape: {:?}", shape);
    }
    
    if let Ok(attr) = x_group.attr("h5sparse_format") {
        if let Ok(fmt) = attr.read_scalar::<hdf5::types::VarLenUnicode>() {
            println!("h5sparse_format (VarLenUnicode): {}", fmt);
        } else if let Ok(fmt) = attr.read_scalar::<hdf5::types::FixedAscii<16>>() {
            println!("h5sparse_format (FixedAscii): {}", fmt);
        } else {
            println!("h5sparse_format: could not read");
        }
    }
    
    // Try reading data
    let data_ds = x_group.dataset("data").unwrap();
    println!("data shape: {:?}", data_ds.shape());
    
    match data_ds.read_1d::<f64>() {
        Ok(d) => println!("f64 read OK, len={}", d.len()),
        Err(e) => println!("f64 read FAILED: {}", e),
    }
    
    match data_ds.read_1d::<f32>() {
        Ok(d) => println!("f32 read OK, len={}", d.len()),
        Err(e) => println!("f32 read FAILED: {}", e),
    }
    
    // Try reading indices
    let indices_ds = x_group.dataset("indices").unwrap();
    println!("indices shape: {:?}", indices_ds.shape());
    
    match indices_ds.read_1d::<i32>() {
        Ok(d) => println!("indices i32 read OK, len={}", d.len()),
        Err(e) => println!("indices i32 read FAILED: {}", e),
    }
    
    // Try reading indptr
    let indptr_ds = x_group.dataset("indptr").unwrap();
    println!("indptr shape: {:?}", indptr_ds.shape());
    
    match indptr_ds.read_1d::<i32>() {
        Ok(d) => println!("indptr i32 read OK, len={}", d.len()),
        Err(e) => println!("indptr i32 read FAILED: {}", e),
    }
    
    // Now try the full read_h5ad
    println!("\n--- Testing full read_h5ad ---");
    match crosscell::anndata::read_h5ad(path) {
        Ok(data) => println!("read_h5ad OK: {} cells x {} genes", data.metadata.n_cells, data.metadata.n_genes),
        Err(e) => println!("read_h5ad FAILED: {}", e),
    }
}
