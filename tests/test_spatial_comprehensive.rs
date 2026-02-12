/// Comprehensive spatial transcriptomics unit tests
/// 
/// Requirements tested: 8.1, 8.2, 8.3, 8.4, 8.5, 8.6

use crosscell::anndata::{read_h5ad, write_h5ad};
use crosscell::seurat::{seurat_rds_to_ir, ir_to_seurat::write_seurat_rds};
use std::path::Path;

// ============================================================================
// Test 1: AnnData Spatial Coordinates (Requirement 8.1)
// ============================================================================

#[test]
fn test_anndata_spatial_coordinates_extraction() {
    let result = read_h5ad("tests/data/spatial_test.h5ad");
    assert!(result.is_ok(), "Failed to read spatial AnnData");
    
    let data = result.unwrap();
    assert!(data.spatial.is_some(), "Spatial coordinates should be extracted");
    
    let spatial = data.spatial.as_ref().unwrap();
    assert_eq!(spatial.n_dims, 2, "Visium data should have 2D coordinates");
    assert_eq!(spatial.n_cells(), 100, "Should have 100 spots");
    
    println!("✅ Test 1: AnnData spatial coordinates extraction");
}

#[test]
fn test_anndata_spatial_3d_coordinates() {
    let result = read_h5ad("tests/data/spatial_3d_test.h5ad");
    assert!(result.is_ok());
    
    let data = result.unwrap();
    let spatial = data.spatial.as_ref().unwrap();
    assert_eq!(spatial.n_dims, 3, "Should have 3D coordinates");
    
    println!("✅ Test 1b: 3D spatial coordinates");
}


// ============================================================================
// Test 2: AnnData Spatial Images (Requirement 8.2)
// ============================================================================

#[test]
fn test_anndata_spatial_image_extraction() {
    let result = read_h5ad("tests/data/spatial_test.h5ad");
    assert!(result.is_ok());
    
    let data = result.unwrap();
    let spatial = data.spatial.as_ref().unwrap();
    
    assert!(spatial.images.is_some(), "Images should be extracted");
    let images = spatial.images.as_ref().unwrap();
    assert_eq!(images.len(), 2, "Should have 2 images (hires and lowres)");
    
    let hires = images.iter().find(|img| img.name == "hires");
    let lowres = images.iter().find(|img| img.name == "lowres");
    
    assert!(hires.is_some(), "Should have hires image");
    assert!(lowres.is_some(), "Should have lowres image");
    
    let hires_img = hires.unwrap();
    assert_eq!(hires_img.width, 200);
    assert_eq!(hires_img.height, 200);
    assert_eq!(hires_img.data.len(), 200 * 200 * 3, "RGB image");
    
    println!("✅ Test 2: AnnData spatial image extraction");
}

// ============================================================================
// Test 3: AnnData Scale Factors (Requirement 8.3)
// ============================================================================

#[test]
fn test_anndata_scale_factors_extraction() {
    let result = read_h5ad("tests/data/spatial_test.h5ad");
    assert!(result.is_ok());
    
    let data = result.unwrap();
    let spatial = data.spatial.as_ref().unwrap();
    
    assert!(spatial.scale_factors.is_some(), "Scale factors should be extracted");
    let scale_factors = spatial.scale_factors.as_ref().unwrap();
    
    assert!(scale_factors.contains_key("tissue_hires_scalef"));
    assert!(scale_factors.contains_key("tissue_lowres_scalef"));
    assert!(scale_factors.contains_key("spot_diameter_fullres"));
    
    let hires_scale = scale_factors["tissue_hires_scalef"];
    assert!(hires_scale > 0.0 && hires_scale < 1.0);
    
    println!("✅ Test 3: AnnData scale factors extraction");
}

// ============================================================================
// Test 4: AnnData Spatial Roundtrip (Requirements 8.1, 8.2, 8.3)
// ============================================================================

#[test]
fn test_anndata_spatial_roundtrip_coordinates() {
    let input_path = "tests/data/spatial_test.h5ad";
    let output_path = "tests/data/spatial_roundtrip_test.h5ad";
    
    let original = read_h5ad(input_path).expect("Failed to read original");
    let orig_spatial = original.spatial.as_ref().unwrap();
    
    write_h5ad(&original, output_path).expect("Failed to write");
    let roundtrip = read_h5ad(output_path).expect("Failed to read roundtrip");
    let rt_spatial = roundtrip.spatial.as_ref().unwrap();
    
    // Verify coordinate preservation (< 0.01 pixel error as per requirement)
    assert_eq!(rt_spatial.coordinates.len(), orig_spatial.coordinates.len());
    
    for (i, (&orig, &rt)) in orig_spatial.coordinates.iter()
        .zip(rt_spatial.coordinates.iter()).enumerate() {
        let error = (orig - rt).abs();
        assert!(error < 0.01, 
                "Coordinate {} error {} exceeds 0.01 pixel threshold", i, error);
    }
    
    if Path::new(output_path).exists() {
        std::fs::remove_file(output_path).ok();
    }
    
    println!("✅ Test 4: AnnData roundtrip coordinates (error < 0.01 pixels)");
}

#[test]
fn test_anndata_spatial_roundtrip_images() {
    let input_path = "tests/data/spatial_test.h5ad";
    let output_path = "tests/data/spatial_roundtrip_images_test.h5ad";
    
    let original = read_h5ad(input_path).expect("Failed to read original");
    let orig_images = original.spatial.as_ref().unwrap().images.as_ref().unwrap();
    
    write_h5ad(&original, output_path).expect("Failed to write");
    let roundtrip = read_h5ad(output_path).expect("Failed to read roundtrip");
    let rt_images = roundtrip.spatial.as_ref().unwrap().images.as_ref().unwrap();
    
    // Verify pixel-by-pixel preservation
    assert_eq!(rt_images.len(), orig_images.len());
    
    for orig_img in orig_images {
        let rt_img = rt_images.iter()
            .find(|img| img.name == orig_img.name)
            .expect(&format!("Image '{}' not found", orig_img.name));
        
        assert_eq!(rt_img.width, orig_img.width);
        assert_eq!(rt_img.height, orig_img.height);
        assert_eq!(rt_img.data.len(), orig_img.data.len());
        
        for (i, (&orig_pixel, &rt_pixel)) in orig_img.data.iter()
            .zip(rt_img.data.iter()).enumerate() {
            assert_eq!(orig_pixel, rt_pixel, 
                      "Image '{}' pixel {} mismatch", orig_img.name, i);
        }
    }
    
    if Path::new(output_path).exists() {
        std::fs::remove_file(output_path).ok();
    }
    
    println!("✅ Test 4b: AnnData roundtrip images (pixel-perfect)");
}


#[test]
fn test_anndata_spatial_roundtrip_scale_factors() {
    let input_path = "tests/data/spatial_test.h5ad";
    let output_path = "tests/data/spatial_roundtrip_sf_test.h5ad";
    
    let original = read_h5ad(input_path).expect("Failed to read original");
    let orig_sf = original.spatial.as_ref().unwrap().scale_factors.as_ref().unwrap();
    
    write_h5ad(&original, output_path).expect("Failed to write");
    let roundtrip = read_h5ad(output_path).expect("Failed to read roundtrip");
    let rt_sf = roundtrip.spatial.as_ref().unwrap().scale_factors.as_ref().unwrap();
    
    // Verify scale factor preservation (< 1e-6 precision as per requirement)
    assert_eq!(rt_sf.len(), orig_sf.len());
    
    for (name, &orig_value) in orig_sf {
        let rt_value = rt_sf.get(name)
            .expect(&format!("Scale factor '{}' not found", name));
        
        let error = (orig_value - rt_value).abs();
        assert!(error < 1e-6, 
                "Scale factor '{}' error {} exceeds 1e-6 threshold", name, error);
    }
    
    if Path::new(output_path).exists() {
        std::fs::remove_file(output_path).ok();
    }
    
    println!("✅ Test 4c: AnnData roundtrip scale factors (precision < 1e-6)");
}

// ============================================================================
// Test 5: Seurat Spatial Extraction (Requirement 8.4)
// ============================================================================

#[test]
fn test_seurat_spatial_coordinates_extraction() {
    let result = seurat_rds_to_ir("tests/data/seurat_spatial_simplified.rds");
    assert!(result.is_ok(), "Failed to read Seurat spatial");
    
    let ir = result.unwrap();
    assert!(ir.spatial.is_some(), "Spatial data should be extracted");
    
    let spatial = ir.spatial.as_ref().unwrap();
    assert_eq!(spatial.n_cells(), 100);
    assert_eq!(spatial.n_dims, 2);
    
    println!("✅ Test 5: Seurat spatial coordinates extraction");
}

#[test]
fn test_seurat_spatial_scale_factors_extraction() {
    let result = seurat_rds_to_ir("tests/data/seurat_spatial_simplified.rds");
    assert!(result.is_ok());
    
    let ir = result.unwrap();
    let spatial = ir.spatial.as_ref().unwrap();
    
    assert!(spatial.scale_factors.is_some());
    let scale_factors = spatial.scale_factors.as_ref().unwrap();
    
    assert!(scale_factors.contains_key("tissue_hires_scalef"));
    assert!(scale_factors.contains_key("spot_diameter_fullres"));
    assert!(scale_factors.contains_key("fiducial_diameter_fullres"));
    
    println!("✅ Test 5b: Seurat scale factors extraction");
}

// ============================================================================
// Test 6: Seurat Spatial Construction (Requirements 8.5, 8.6)
// ============================================================================

#[test]
fn test_seurat_spatial_roundtrip_coordinates() {
    let input_path = "tests/data/seurat_spatial_simplified.rds";
    let output_path = "tests/data/seurat_spatial_roundtrip_test.rds";
    
    let original = seurat_rds_to_ir(input_path).expect("Failed to read original");
    let orig_spatial = original.spatial.as_ref().unwrap();
    
    write_seurat_rds(&original, output_path).expect("Failed to write");
    let roundtrip = seurat_rds_to_ir(output_path).expect("Failed to read roundtrip");
    let rt_spatial = roundtrip.spatial.as_ref().unwrap();
    
    // Verify coordinate preservation
    assert_eq!(rt_spatial.n_cells(), orig_spatial.n_cells());
    assert_eq!(rt_spatial.n_dims, orig_spatial.n_dims);
    
    for i in 0..orig_spatial.n_cells() {
        let orig_coords = orig_spatial.get_coordinates(i).unwrap();
        let rt_coords = rt_spatial.get_coordinates(i).unwrap();
        
        for j in 0..orig_spatial.n_dims {
            let error = (orig_coords[j] - rt_coords[j]).abs();
            assert!(error < 1e-7, 
                    "Coordinate error at cell {}, dim {}: {}", i, j, error);
        }
    }
    
    if Path::new(output_path).exists() {
        std::fs::remove_file(output_path).ok();
    }
    
    println!("✅ Test 6: Seurat roundtrip coordinates");
}


// ============================================================================
// Test 7: Cross-Format Conversion (AnnData → Seurat)
// ============================================================================

#[test]
fn test_cross_format_anndata_to_seurat_coordinates() {
    // Test AnnData → IR → Seurat conversion preserves coordinates
    let anndata_path = "tests/data/spatial_test.h5ad";
    let seurat_path = "tests/data/spatial_cross_format_test.rds";
    
    // Read AnnData
    let anndata_ir = read_h5ad(anndata_path).expect("Failed to read AnnData");
    let anndata_spatial = anndata_ir.spatial.as_ref().unwrap();
    
    // Convert to Seurat
    write_seurat_rds(&anndata_ir, seurat_path).expect("Failed to write Seurat");
    let seurat_ir = seurat_rds_to_ir(seurat_path).expect("Failed to read Seurat");
    let seurat_spatial = seurat_ir.spatial.as_ref().unwrap();
    
    // Verify coordinate preservation
    assert_eq!(seurat_spatial.n_cells(), anndata_spatial.n_cells());
    assert_eq!(seurat_spatial.n_dims, anndata_spatial.n_dims);
    
    for i in 0..anndata_spatial.n_cells() {
        let anndata_coords = anndata_spatial.get_coordinates(i).unwrap();
        let seurat_coords = seurat_spatial.get_coordinates(i).unwrap();
        
        for j in 0..anndata_spatial.n_dims {
            let error = (anndata_coords[j] - seurat_coords[j]).abs();
            assert!(error < 0.01, 
                    "Cross-format coordinate error at cell {}, dim {}: {}", i, j, error);
        }
    }
    
    if Path::new(seurat_path).exists() {
        std::fs::remove_file(seurat_path).ok();
    }
    
    println!("✅ Test 7: Cross-format AnnData → Seurat coordinates");
}

#[test]
fn test_cross_format_anndata_to_seurat_scale_factors() {
    let anndata_path = "tests/data/spatial_test.h5ad";
    let seurat_path = "tests/data/spatial_cross_format_sf_test.rds";
    
    let anndata_ir = read_h5ad(anndata_path).expect("Failed to read AnnData");
    let anndata_sf = anndata_ir.spatial.as_ref().unwrap()
        .scale_factors.as_ref().unwrap();
    
    write_seurat_rds(&anndata_ir, seurat_path).expect("Failed to write Seurat");
    let seurat_ir = seurat_rds_to_ir(seurat_path).expect("Failed to read Seurat");
    let seurat_sf = seurat_ir.spatial.as_ref().unwrap()
        .scale_factors.as_ref().unwrap();
    
    // Verify scale factors preserved
    for (name, &anndata_value) in anndata_sf {
        let seurat_value = seurat_sf.get(name)
            .expect(&format!("Scale factor '{}' not found in Seurat", name));
        
        let error = (anndata_value - seurat_value).abs();
        assert!(error < 1e-6, 
                "Cross-format scale factor '{}' error: {}", name, error);
    }
    
    if Path::new(seurat_path).exists() {
        std::fs::remove_file(seurat_path).ok();
    }
    
    println!("✅ Test 7b: Cross-format scale factors preservation");
}

// ============================================================================
// Test 8: Coordinate System Alignment (Requirement 8.5)
// ============================================================================

#[test]
fn test_coordinate_system_alignment() {
    // Test that coordinate systems remain aligned across conversions
    let anndata_path = "tests/data/spatial_test.h5ad";
    let seurat_path = "tests/data/spatial_alignment_test.rds";
    let roundtrip_path = "tests/data/spatial_alignment_roundtrip.h5ad";
    
    // Original AnnData
    let original = read_h5ad(anndata_path).expect("Failed to read original");
    let orig_coords = &original.spatial.as_ref().unwrap().coordinates;
    
    // AnnData → Seurat
    write_seurat_rds(&original, seurat_path).expect("Failed to write Seurat");
    let seurat_ir = seurat_rds_to_ir(seurat_path).expect("Failed to read Seurat");
    let seurat_coords = &seurat_ir.spatial.as_ref().unwrap().coordinates;
    
    // Seurat → AnnData
    write_h5ad(&seurat_ir, roundtrip_path).expect("Failed to write roundtrip");
    let roundtrip = read_h5ad(roundtrip_path).expect("Failed to read roundtrip");
    let rt_coords = &roundtrip.spatial.as_ref().unwrap().coordinates;
    
    // Verify alignment through full cycle
    assert_eq!(orig_coords.len(), seurat_coords.len());
    assert_eq!(orig_coords.len(), rt_coords.len());
    
    for i in 0..orig_coords.len() {
        let orig_to_seurat_error = (orig_coords[i] - seurat_coords[i]).abs();
        let orig_to_rt_error = (orig_coords[i] - rt_coords[i]).abs();
        
        assert!(orig_to_seurat_error < 0.01, 
                "Alignment error (AnnData→Seurat) at index {}: {}", i, orig_to_seurat_error);
        assert!(orig_to_rt_error < 0.01, 
                "Alignment error (full cycle) at index {}: {}", i, orig_to_rt_error);
    }
    
    // Cleanup
    if Path::new(seurat_path).exists() {
        std::fs::remove_file(seurat_path).ok();
    }
    if Path::new(roundtrip_path).exists() {
        std::fs::remove_file(roundtrip_path).ok();
    }
    
    println!("✅ Test 8: Coordinate system alignment across formats");
}

// ============================================================================
// Test 9: Visium-Specific Tests
// ============================================================================

#[test]
fn test_visium_spot_coordinates() {
    // Test Visium-specific spot coordinate handling
    let result = read_h5ad("tests/data/spatial_test.h5ad");
    assert!(result.is_ok());
    
    let data = result.unwrap();
    let spatial = data.spatial.as_ref().unwrap();
    
    // Verify Visium has 2D coordinates
    assert_eq!(spatial.n_dims, 2, "Visium should have 2D coordinates");
    
    // Verify coordinates are in reasonable range for Visium
    for i in 0..spatial.n_cells() {
        let coords = spatial.get_coordinates(i).unwrap();
        
        // Visium coordinates should be positive
        assert!(coords[0] >= 0.0, "X coordinate should be non-negative");
        assert!(coords[1] >= 0.0, "Y coordinate should be non-negative");
        
        // Visium coordinates typically in range [0, 10000]
        assert!(coords[0] < 10000.0, "X coordinate should be reasonable");
        assert!(coords[1] < 10000.0, "Y coordinate should be reasonable");
    }
    
    println!("✅ Test 9: Visium spot coordinates validation");
}

#[test]
fn test_visium_scale_factors_values() {
    // Test that Visium scale factors have reasonable values
    let result = read_h5ad("tests/data/spatial_test.h5ad");
    assert!(result.is_ok());
    
    let data = result.unwrap();
    let scale_factors = data.spatial.as_ref().unwrap()
        .scale_factors.as_ref().unwrap();
    
    // Verify scale factor ranges
    let hires_scale = scale_factors.get("tissue_hires_scalef").unwrap();
    assert!(*hires_scale > 0.0 && *hires_scale < 1.0, 
            "Hires scale factor should be between 0 and 1");
    
    let lowres_scale = scale_factors.get("tissue_lowres_scalef").unwrap();
    assert!(*lowres_scale > 0.0 && *lowres_scale < 1.0, 
            "Lowres scale factor should be between 0 and 1");
    
    let spot_diameter = scale_factors.get("spot_diameter_fullres").unwrap();
    assert!(*spot_diameter > 0.0 && *spot_diameter < 200.0, 
            "Spot diameter should be reasonable (typically 50-150 pixels)");
    
    println!("✅ Test 9b: Visium scale factors validation");
}

// ============================================================================
// Test 10: Edge Cases
// ============================================================================

#[test]
fn test_spatial_without_images() {
    // Test spatial data without images (coordinates only)
    let result = read_h5ad("tests/data/spatial_no_images_test.h5ad");
    assert!(result.is_ok());
    
    let data = result.unwrap();
    let spatial = data.spatial.as_ref().unwrap();
    
    assert!(spatial.images.is_none(), "Should not have images");
    assert!(spatial.scale_factors.is_none(), "Should not have scale factors");
    assert!(spatial.coordinates.len() > 0, "Should have coordinates");
    
    println!("✅ Test 10: Spatial data without images");
}

#[test]
fn test_non_spatial_data() {
    // Test that non-spatial data doesn't have spatial field
    let result = read_h5ad("tests/data/small_sparse.h5ad");
    
    if result.is_ok() {
        let data = result.unwrap();
        assert!(data.spatial.is_none(), "Non-spatial data should not have spatial field");
    }
    
    println!("✅ Test 10b: Non-spatial data handling");
}
