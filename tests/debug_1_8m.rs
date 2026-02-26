#[test]
fn debug_1_8m_x_type() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug")).init();

    let path = "/workspace/data/generated/cellxgene_immune_1.8M.h5ad";

    if !std::path::Path::new(path).exists() {
        println!("File not found, skipping: {}", path);
        return;
    }

    println!("Opening: {}", path);

    // HDF5 版本
    let ver = hdf5::library_version();
    println!("HDF5 library version: {:?}", ver);

    let file = hdf5::File::open(path).expect("Failed to open file");
    println!("File opened successfully");

    let members = file.member_names().expect("Failed to get members");
    println!("Root members: {:?}", members);
    println!("link_exists(\"X\"): {}", file.link_exists("X"));

    // 尝试 group
    println!("\nTrying file.group(\"X\")...");
    match file.group("X") {
        Ok(g) => {
            println!("  SUCCESS: X is a Group");
            let sub = g.member_names().unwrap_or_default();
            println!("  members: {:?}", sub);
        }
        Err(e) => {
            println!("  FAILED: {}", e);
        }
    }

    // 尝试 dataset
    println!("\nTrying file.dataset(\"X\")...");
    match file.dataset("X") {
        Ok(d) => {
            println!("  SUCCESS: X is a Dataset, shape: {:?}", d.shape());
        }
        Err(e) => {
            println!("  FAILED: {}", e);
        }
    }
}
