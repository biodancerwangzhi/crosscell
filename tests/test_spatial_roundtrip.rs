/// 测试空间转录组学数据的往返一致性
/// 
/// 验证 AnnData → IR → AnnData 往返后空间数据保持一致

use crosscell::anndata::{read_h5ad, write_h5ad};
use std::path::Path;

#[test]
fn test_spatial_roundtrip_with_images() {
    // 读取包含空间数据的 .h5ad 文件
    let input_path = "tests/data/spatial_test.h5ad";
    let output_path = "tests/data/spatial_roundtrip.h5ad";
    
    // 读取原始文件
    let data = read_h5ad(input_path).expect("Failed to read input file");
    
    // 验证空间数据存在
    assert!(data.spatial.is_some(), "Spatial data should be present");
    let spatial = data.spatial.as_ref().unwrap();
    
    // 记录原始数据
    let orig_n_cells = spatial.n_cells();
    let orig_n_dims = spatial.n_dims;
    let orig_coords = spatial.coordinates.clone();
    let orig_has_images = spatial.images.is_some();
    let orig_has_scale_factors = spatial.scale_factors.is_some();
    
    // 写入新文件
    write_h5ad(&data, output_path).expect("Failed to write output file");
    
    // 读取往返后的文件
    let roundtrip_data = read_h5ad(output_path).expect("Failed to read roundtrip file");
    
    // 验证空间数据仍然存在
    assert!(roundtrip_data.spatial.is_some(), "Spatial data should be preserved");
    let rt_spatial = roundtrip_data.spatial.as_ref().unwrap();
    
    // 验证基本属性
    assert_eq!(rt_spatial.n_cells(), orig_n_cells, "Cell count should match");
    assert_eq!(rt_spatial.n_dims, orig_n_dims, "Dimension count should match");
    
    // 验证坐标数据
    assert_eq!(rt_spatial.coordinates.len(), orig_coords.len(), "Coordinate array length should match");
    for (i, (&orig, &rt)) in orig_coords.iter().zip(rt_spatial.coordinates.iter()).enumerate() {
        assert!(
            (orig - rt).abs() < 1e-10,
            "Coordinate {} mismatch: {} vs {}",
            i, orig, rt
        );
    }
    
    // 验证图像存在性
    assert_eq!(
        rt_spatial.images.is_some(),
        orig_has_images,
        "Image presence should match"
    );
    
    // 如果有图像，验证图像数据
    if let (Some(orig_images), Some(rt_images)) = (&spatial.images, &rt_spatial.images) {
        assert_eq!(orig_images.len(), rt_images.len(), "Number of images should match");
        
        for orig_img in orig_images {
            // 查找对应的图像
            let rt_img = rt_images.iter()
                .find(|img| img.name == orig_img.name)
                .expect(&format!("Image '{}' should be present", orig_img.name));
            
            // 验证图像属性
            assert_eq!(rt_img.width, orig_img.width, "Image '{}' width should match", orig_img.name);
            assert_eq!(rt_img.height, orig_img.height, "Image '{}' height should match", orig_img.name);
            assert_eq!(rt_img.data.len(), orig_img.data.len(), "Image '{}' data length should match", orig_img.name);
            
            // 验证图像数据（逐像素比较）
            for (i, (&orig_pixel, &rt_pixel)) in orig_img.data.iter().zip(rt_img.data.iter()).enumerate() {
                assert_eq!(
                    orig_pixel, rt_pixel,
                    "Image '{}' pixel {} mismatch: {} vs {}",
                    orig_img.name, i, orig_pixel, rt_pixel
                );
            }
        }
    }
    
    // 验证缩放因子存在性
    assert_eq!(
        rt_spatial.scale_factors.is_some(),
        orig_has_scale_factors,
        "Scale factors presence should match"
    );
    
    // 如果有缩放因子，验证数据
    if let (Some(orig_sf), Some(rt_sf)) = (&spatial.scale_factors, &rt_spatial.scale_factors) {
        assert_eq!(orig_sf.len(), rt_sf.len(), "Number of scale factors should match");
        
        for (name, &orig_value) in orig_sf {
            let rt_value = rt_sf.get(name)
                .expect(&format!("Scale factor '{}' should be present", name));
            
            assert!(
                (orig_value - rt_value).abs() < 1e-10,
                "Scale factor '{}' mismatch: {} vs {}",
                name, orig_value, rt_value
            );
        }
    }
    
    // 清理测试文件
    if Path::new(output_path).exists() {
        std::fs::remove_file(output_path).ok();
    }
}

#[test]
fn test_spatial_roundtrip_3d() {
    // 测试 3D 坐标的往返
    let input_path = "tests/data/spatial_3d_test.h5ad";
    let output_path = "tests/data/spatial_3d_roundtrip.h5ad";
    
    // 读取原始文件
    let data = read_h5ad(input_path).expect("Failed to read input file");
    
    // 验证是 3D 坐标
    assert!(data.spatial.is_some());
    let spatial = data.spatial.as_ref().unwrap();
    assert_eq!(spatial.n_dims, 3, "Should be 3D coordinates");
    
    // 写入并读取
    write_h5ad(&data, output_path).expect("Failed to write output file");
    let roundtrip_data = read_h5ad(output_path).expect("Failed to read roundtrip file");
    
    // 验证 3D 坐标保留
    assert!(roundtrip_data.spatial.is_some());
    let rt_spatial = roundtrip_data.spatial.as_ref().unwrap();
    assert_eq!(rt_spatial.n_dims, 3, "Should still be 3D coordinates");
    assert_eq!(rt_spatial.coordinates.len(), spatial.coordinates.len());
    
    // 清理
    if Path::new(output_path).exists() {
        std::fs::remove_file(output_path).ok();
    }
}

#[test]
fn test_spatial_roundtrip_no_images() {
    // 测试只有坐标没有图像的往返
    let input_path = "tests/data/spatial_no_images_test.h5ad";
    let output_path = "tests/data/spatial_no_images_roundtrip.h5ad";
    
    // 读取原始文件
    let data = read_h5ad(input_path).expect("Failed to read input file");
    
    // 验证有坐标但没有图像
    assert!(data.spatial.is_some());
    let spatial = data.spatial.as_ref().unwrap();
    assert!(spatial.images.is_none(), "Should not have images");
    assert!(spatial.scale_factors.is_none(), "Should not have scale factors");
    
    // 写入并读取
    write_h5ad(&data, output_path).expect("Failed to write output file");
    let roundtrip_data = read_h5ad(output_path).expect("Failed to read roundtrip file");
    
    // 验证坐标保留，图像仍然不存在
    assert!(roundtrip_data.spatial.is_some());
    let rt_spatial = roundtrip_data.spatial.as_ref().unwrap();
    assert!(rt_spatial.images.is_none(), "Should still not have images");
    assert!(rt_spatial.scale_factors.is_none(), "Should still not have scale factors");
    
    // 验证坐标数据
    assert_eq!(rt_spatial.coordinates.len(), spatial.coordinates.len());
    
    // 清理
    if Path::new(output_path).exists() {
        std::fs::remove_file(output_path).ok();
    }
}
