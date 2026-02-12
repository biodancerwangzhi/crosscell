//! 空间转录组学 IR 类型
//!
//! 存储空间坐标、图像和缩放因子

use std::collections::HashMap;
use std::fmt;

/// 空间转录组数据
///
/// 对应：
/// - AnnData.obsm['spatial'], AnnData.uns['spatial']
/// - Seurat@images
#[derive(Debug, Clone, PartialEq)]
pub struct SpatialData {
    /// 空间坐标（cells × 2 或 cells × 3）
    pub coordinates: Vec<f64>,
    /// 空间维度（2D 或 3D）
    pub n_dims: usize,
    /// 组织切片图像（可选）
    pub images: Option<Vec<SpatialImage>>,
    /// 坐标缩放因子（可选）
    pub scale_factors: Option<HashMap<String, f64>>,
}

impl SpatialData {
    /// 创建新的空间数据
    pub fn new(
        coordinates: Vec<f64>,
        n_dims: usize,
        images: Option<Vec<SpatialImage>>,
        scale_factors: Option<HashMap<String, f64>>,
    ) -> Result<Self, String> {
        let spatial = Self {
            coordinates,
            n_dims,
            images,
            scale_factors,
        };
        spatial.validate()?;
        Ok(spatial)
    }

    /// 验证空间数据
    pub fn validate(&self) -> Result<(), String> {
        // 检查维度
        if self.n_dims != 2 && self.n_dims != 3 {
            return Err(format!(
                "SpatialData: n_dims must be 2 or 3, got {}",
                self.n_dims
            ));
        }

        // 检查坐标长度
        if !self.coordinates.len().is_multiple_of(self.n_dims) {
            return Err(format!(
                "SpatialData: coordinates length {} is not divisible by n_dims {}",
                self.coordinates.len(),
                self.n_dims
            ));
        }

        Ok(())
    }

    /// 获取细胞数量
    pub fn n_cells(&self) -> usize {
        self.coordinates.len() / self.n_dims
    }

    /// 获取指定细胞的坐标
    pub fn get_coordinates(&self, cell_idx: usize) -> Option<&[f64]> {
        let start = cell_idx * self.n_dims;
        let end = start + self.n_dims;
        if end <= self.coordinates.len() {
            Some(&self.coordinates[start..end])
        } else {
            None
        }
    }
}

impl fmt::Display for SpatialData {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "SpatialData({} cells, {}D, {} images)",
            self.n_cells(),
            self.n_dims,
            self.images.as_ref().map_or(0, |imgs| imgs.len())
        )
    }
}

/// 空间图像数据
///
/// 对应：
/// - AnnData.uns['spatial']['images']['hires']
/// - Seurat@images$slice1@image
#[derive(Debug, Clone, PartialEq)]
pub struct SpatialImage {
    /// 图像名称/类型（如 "hires", "lowres"）
    pub name: String,
    /// 图像数据（PNG/JPEG 编码）
    pub data: Vec<u8>,
    /// 图像宽度
    pub width: usize,
    /// 图像高度
    pub height: usize,
}

impl SpatialImage {
    /// 创建新的空间图像
    pub fn new(name: String, data: Vec<u8>, width: usize, height: usize) -> Self {
        Self {
            name,
            data,
            width,
            height,
        }
    }

    /// 获取图像大小（字节）
    pub fn size_bytes(&self) -> usize {
        self.data.len()
    }
}

impl fmt::Display for SpatialImage {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "SpatialImage('{}', {} × {}, {} bytes)",
            self.name,
            self.width,
            self.height,
            self.size_bytes()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_spatial_data_2d() {
        // 3 细胞，2D 坐标
        let coords = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let spatial = SpatialData::new(coords, 2, None, None).unwrap();

        assert_eq!(spatial.n_cells(), 3);
        assert_eq!(spatial.n_dims, 2);
        assert_eq!(spatial.get_coordinates(0), Some(&[1.0, 2.0][..]));
        assert_eq!(spatial.get_coordinates(2), Some(&[5.0, 6.0][..]));
    }

    #[test]
    fn test_spatial_data_3d() {
        // 2 细胞，3D 坐标
        let coords = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let spatial = SpatialData::new(coords, 3, None, None).unwrap();

        assert_eq!(spatial.n_cells(), 2);
        assert_eq!(spatial.n_dims, 3);
        assert_eq!(spatial.get_coordinates(0), Some(&[1.0, 2.0, 3.0][..]));
    }

    #[test]
    fn test_spatial_data_invalid_dims() {
        // 无效维度
        let coords = vec![1.0, 2.0];
        let result = SpatialData::new(coords, 4, None, None);
        assert!(result.is_err());
    }

    #[test]
    fn test_spatial_data_mismatched_length() {
        // 坐标长度不能被维度整除
        let coords = vec![1.0, 2.0, 3.0];
        let result = SpatialData::new(coords, 2, None, None);
        assert!(result.is_err());
    }

    #[test]
    fn test_spatial_image() {
        let img = SpatialImage::new("hires".to_string(), vec![0u8; 1000], 100, 100);

        assert_eq!(img.name, "hires");
        assert_eq!(img.width, 100);
        assert_eq!(img.height, 100);
        assert_eq!(img.size_bytes(), 1000);
    }

    #[test]
    fn test_spatial_with_scale_factors() {
        let mut scale_factors = HashMap::new();
        scale_factors.insert("tissue_hires_scalef".to_string(), 0.08);
        scale_factors.insert("spot_diameter_fullres".to_string(), 89.43);

        let coords = vec![1.0, 2.0];
        let spatial = SpatialData::new(coords, 2, None, Some(scale_factors)).unwrap();

        assert!(spatial.scale_factors.is_some());
        assert_eq!(spatial.scale_factors.as_ref().unwrap().len(), 2);
    }
}
