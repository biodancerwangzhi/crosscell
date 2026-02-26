use crosscell::seurat::seurat_rds_to_ir;

#[test]
fn test_seurat_spatial_extraction() {
    // Read simplified Seurat with spatial data
    let result = seurat_rds_to_ir("tests/data/seurat_spatial_simplified.rds");

    if let Err(e) = &result {
        eprintln!("Failed to parse Seurat with spatial data: {:?}", e);
    }

    assert!(result.is_ok(), "Failed to parse Seurat with spatial data");

    let ir = result.unwrap();

    // Verify basic structure
    assert_eq!(ir.metadata.n_cells, 100);
    assert_eq!(ir.metadata.n_genes, 50);

    // Verify spatial data exists
    assert!(ir.spatial.is_some(), "Spatial data should be present");

    let spatial = ir.spatial.as_ref().unwrap();

    // Verify spatial coordinates
    assert_eq!(spatial.n_cells(), 100, "Should have 100 cells");
    assert_eq!(spatial.n_dims, 2, "Should be 2D coordinates");
    assert_eq!(
        spatial.coordinates.len(),
        200,
        "Should have 200 coordinate values (100 cells × 2 dims)"
    );

    // Verify scale factors
    assert!(
        spatial.scale_factors.is_some(),
        "Scale factors should be present"
    );
    let scale_factors = spatial.scale_factors.as_ref().unwrap();
    assert_eq!(scale_factors.len(), 3, "Should have 3 scale factors");
    assert!(scale_factors.contains_key("tissue_hires_scalef"));
    assert!(scale_factors.contains_key("spot_diameter_fullres"));
    assert!(scale_factors.contains_key("fiducial_diameter_fullres"));

    // Verify scale factor values
    assert!((scale_factors["tissue_hires_scalef"] - 0.08).abs() < 1e-6);
    assert!((scale_factors["spot_diameter_fullres"] - 89.43).abs() < 1e-6);
    assert!((scale_factors["fiducial_diameter_fullres"] - 144.0).abs() < 1e-6);

    println!("✅ Seurat spatial data extraction test passed!");
    println!("   Cells: {}", spatial.n_cells());
    println!("   Dimensions: {}D", spatial.n_dims);
    println!("   Scale factors: {}", scale_factors.len());
}

#[test]
fn test_seurat_spatial_coordinates() {
    let result = seurat_rds_to_ir("tests/data/seurat_spatial_simplified.rds");
    assert!(result.is_ok());

    let ir = result.unwrap();
    let spatial = ir.spatial.as_ref().unwrap();

    // Test getting individual cell coordinates
    let coord_0 = spatial.get_coordinates(0);
    assert!(
        coord_0.is_some(),
        "Should be able to get coordinates for cell 0"
    );
    assert_eq!(coord_0.unwrap().len(), 2, "Coordinates should be 2D");

    let coord_99 = spatial.get_coordinates(99);
    assert!(
        coord_99.is_some(),
        "Should be able to get coordinates for cell 99"
    );

    let coord_100 = spatial.get_coordinates(100);
    assert!(
        coord_100.is_none(),
        "Should return None for out-of-bounds index"
    );

    println!("✅ Seurat spatial coordinates test passed!");
}

#[test]
fn test_seurat_spatial_roundtrip() {
    use crosscell::seurat::ir_to_seurat::write_seurat_rds;
    use std::path::Path;

    // 1. Read Seurat with spatial data
    let result = seurat_rds_to_ir("tests/data/seurat_spatial_simplified.rds");
    assert!(
        result.is_ok(),
        "Failed to read Seurat spatial: {:?}",
        result.err()
    );

    let ir = result.unwrap();

    // 2. Verify spatial data was extracted
    assert!(ir.spatial.is_some(), "Spatial data not extracted");
    let spatial = ir.spatial.as_ref().unwrap();
    assert_eq!(spatial.n_cells(), 100);
    assert_eq!(spatial.n_dims, 2);

    // 3. Write back to RDS
    let output_path = "tests/data/seurat_spatial_roundtrip.rds";
    let write_result = write_seurat_rds(&ir, output_path);
    assert!(
        write_result.is_ok(),
        "Failed to write Seurat spatial: {:?}",
        write_result.err()
    );

    // 4. Verify file was created
    assert!(Path::new(output_path).exists(), "Output file not created");

    // 5. Read back the roundtrip file
    let roundtrip_result = seurat_rds_to_ir(output_path);
    assert!(
        roundtrip_result.is_ok(),
        "Failed to read roundtrip Seurat: {:?}",
        roundtrip_result.err()
    );

    let roundtrip_ir = roundtrip_result.unwrap();

    // 6. Verify spatial data was preserved
    assert!(
        roundtrip_ir.spatial.is_some(),
        "Spatial data not preserved in roundtrip"
    );
    let roundtrip_spatial = roundtrip_ir.spatial.as_ref().unwrap();

    // Check dimensions
    assert_eq!(roundtrip_spatial.n_cells(), spatial.n_cells());
    assert_eq!(roundtrip_spatial.n_dims, spatial.n_dims);

    // Check coordinates (sample a few cells)
    for i in [0, 49, 99].iter() {
        let orig_coords = spatial.get_coordinates(*i).unwrap();
        let roundtrip_coords = roundtrip_spatial.get_coordinates(*i).unwrap();

        for j in 0..spatial.n_dims {
            let diff = (orig_coords[j] - roundtrip_coords[j]).abs();
            assert!(
                diff < 1e-7,
                "Coordinate mismatch at cell {}, dim {}: {} vs {}",
                i,
                j,
                orig_coords[j],
                roundtrip_coords[j]
            );
        }
    }

    // Check scale factors
    if let Some(orig_sf) = &spatial.scale_factors {
        assert!(
            roundtrip_spatial.scale_factors.is_some(),
            "Scale factors not preserved"
        );
        let roundtrip_sf = roundtrip_spatial.scale_factors.as_ref().unwrap();

        for (key, value) in orig_sf {
            assert!(
                roundtrip_sf.contains_key(key),
                "Scale factor '{}' missing",
                key
            );
            let roundtrip_value = roundtrip_sf.get(key).unwrap();
            let diff = (value - roundtrip_value).abs();
            assert!(
                diff < 1e-7,
                "Scale factor '{}' mismatch: {} vs {}",
                key,
                value,
                roundtrip_value
            );
        }
    }

    println!("✅ Seurat spatial roundtrip test passed!");
    println!("   Original cells: {}", spatial.n_cells());
    println!("   Roundtrip cells: {}", roundtrip_spatial.n_cells());
    println!("   Coordinate preservation: ✓");
    println!("   Scale factors preservation: ✓");
}
