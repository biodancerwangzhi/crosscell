/// 测试空间转录组学数据的读取
///
/// 验证从 .h5ad 文件中提取空间坐标、图像和缩放因子
use crosscell::anndata::read_h5ad;

#[test]
fn test_read_spatial_with_images() {
    // 读取包含空间坐标、图像和缩放因子的 .h5ad 文件
    let result = read_h5ad("tests/data/spatial_test.h5ad");
    assert!(result.is_ok(), "Failed to read spatial test file");

    let data = result.unwrap();

    // 验证基础数据
    assert_eq!(data.metadata.n_cells, 100);
    assert_eq!(data.metadata.n_genes, 50);

    // 验证空间数据存在
    assert!(data.spatial.is_some(), "Spatial data should be present");

    let spatial = data.spatial.unwrap();

    // 验证空间坐标
    assert_eq!(spatial.n_dims, 2, "Should be 2D coordinates");
    assert_eq!(spatial.n_cells(), 100, "Should have 100 cells");
    assert_eq!(
        spatial.coordinates.len(),
        200,
        "Should have 200 coordinate values (100 cells × 2 dims)"
    );

    // 验证图像数据
    assert!(spatial.images.is_some(), "Images should be present");
    let images = spatial.images.unwrap();
    assert_eq!(images.len(), 2, "Should have 2 images (hires and lowres)");

    // 查找 hires 和 lowres 图像
    let hires = images.iter().find(|img| img.name == "hires");
    let lowres = images.iter().find(|img| img.name == "lowres");

    assert!(hires.is_some(), "Should have hires image");
    assert!(lowres.is_some(), "Should have lowres image");

    // 验证 hires 图像属性
    let hires_img = hires.unwrap();
    assert_eq!(hires_img.width, 200, "Hires width should be 200");
    assert_eq!(hires_img.height, 200, "Hires height should be 200");
    assert_eq!(
        hires_img.data.len(),
        200 * 200 * 3,
        "Hires should have 200×200×3 bytes"
    );

    // 验证 lowres 图像属性
    let lowres_img = lowres.unwrap();
    assert_eq!(lowres_img.width, 100, "Lowres width should be 100");
    assert_eq!(lowres_img.height, 100, "Lowres height should be 100");
    assert_eq!(
        lowres_img.data.len(),
        100 * 100 * 3,
        "Lowres should have 100×100×3 bytes"
    );

    // 验证缩放因子
    assert!(
        spatial.scale_factors.is_some(),
        "Scale factors should be present"
    );
    let scale_factors = spatial.scale_factors.unwrap();
    assert_eq!(scale_factors.len(), 3, "Should have 3 scale factors");

    // 验证特定的缩放因子
    assert!(scale_factors.contains_key("tissue_hires_scalef"));
    assert!(scale_factors.contains_key("tissue_lowres_scalef"));
    assert!(scale_factors.contains_key("spot_diameter_fullres"));

    // 验证缩放因子的值
    assert!((scale_factors["tissue_hires_scalef"] - 0.08).abs() < 1e-6);
    assert!((scale_factors["tissue_lowres_scalef"] - 0.04).abs() < 1e-6);
    assert!((scale_factors["spot_diameter_fullres"] - 89.43).abs() < 1e-6);
}

#[test]
fn test_read_spatial_3d() {
    // 读取包含 3D 空间坐标的文件
    let result = read_h5ad("tests/data/spatial_3d_test.h5ad");
    assert!(result.is_ok(), "Failed to read 3D spatial test file");

    let data = result.unwrap();

    // 验证基础数据
    assert_eq!(data.metadata.n_cells, 50);

    // 验证空间数据
    assert!(data.spatial.is_some(), "Spatial data should be present");
    let spatial = data.spatial.unwrap();

    // 验证 3D 坐标
    assert_eq!(spatial.n_dims, 3, "Should be 3D coordinates");
    assert_eq!(spatial.n_cells(), 50, "Should have 50 cells");
    assert_eq!(
        spatial.coordinates.len(),
        150,
        "Should have 150 coordinate values (50 cells × 3 dims)"
    );

    // 验证没有图像数据（这个测试文件不包含图像）
    assert!(spatial.images.is_none(), "Should not have images");
}

#[test]
fn test_read_spatial_no_images() {
    // 读取只有坐标没有图像的文件
    let result = read_h5ad("tests/data/spatial_no_images_test.h5ad");
    assert!(
        result.is_ok(),
        "Failed to read spatial test file (no images)"
    );

    let data = result.unwrap();

    // 验证基础数据
    assert_eq!(data.metadata.n_cells, 80);

    // 验证空间数据
    assert!(data.spatial.is_some(), "Spatial data should be present");
    let spatial = data.spatial.unwrap();

    // 验证坐标
    assert_eq!(spatial.n_dims, 2, "Should be 2D coordinates");
    assert_eq!(spatial.n_cells(), 80, "Should have 80 cells");

    // 验证没有图像和缩放因子
    assert!(spatial.images.is_none(), "Should not have images");
    assert!(
        spatial.scale_factors.is_none(),
        "Should not have scale factors"
    );
}

#[test]
fn test_read_non_spatial_h5ad() {
    // 读取不包含空间数据的普通 .h5ad 文件
    let result = read_h5ad("tests/data/small_sparse.h5ad");

    if result.is_ok() {
        let data = result.unwrap();
        // 验证没有空间数据
        assert!(
            data.spatial.is_none(),
            "Non-spatial file should not have spatial data"
        );
    }
}

#[test]
fn test_spatial_coordinate_access() {
    // 测试空间坐标的访问方法
    let result = read_h5ad("tests/data/spatial_test.h5ad");
    assert!(result.is_ok());

    let data = result.unwrap();
    let spatial = data.spatial.unwrap();

    // 测试获取特定细胞的坐标
    let cell_0_coords = spatial.get_coordinates(0);
    assert!(cell_0_coords.is_some());
    assert_eq!(cell_0_coords.unwrap().len(), 2);

    let cell_50_coords = spatial.get_coordinates(50);
    assert!(cell_50_coords.is_some());
    assert_eq!(cell_50_coords.unwrap().len(), 2);

    // 测试越界访问
    let invalid_coords = spatial.get_coordinates(1000);
    assert!(
        invalid_coords.is_none(),
        "Out of bounds access should return None"
    );
}

#[test]
fn test_spatial_image_properties() {
    // 测试空间图像的属性和方法
    let result = read_h5ad("tests/data/spatial_test.h5ad");
    assert!(result.is_ok());

    let data = result.unwrap();
    let spatial = data.spatial.unwrap();
    let images = spatial.images.unwrap();

    for img in &images {
        // 验证图像数据不为空
        assert!(!img.data.is_empty(), "Image data should not be empty");

        // 验证图像尺寸合理
        assert!(img.width > 0, "Image width should be positive");
        assert!(img.height > 0, "Image height should be positive");

        // 验证数据大小匹配（RGB 图像）
        assert_eq!(
            img.data.len(),
            img.width * img.height * 3,
            "Image data size should match width × height × 3"
        );

        // 测试 size_bytes 方法
        assert_eq!(img.size_bytes(), img.data.len());
    }
}
