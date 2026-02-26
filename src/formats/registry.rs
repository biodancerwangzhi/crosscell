//! 格式注册表
//!
//! 管理所有已注册的格式转换器，提供查找和列表功能。

use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;

use super::anndata::AnndataConverter;
use super::loom::LoomConverter;
use super::sce::SceConverter;
use super::seurat::SeuratConverter;
use super::{FormatConverter, FormatInfo};

/// 格式注册表
///
/// 存储所有已注册的格式转换器，支持按名称或扩展名查找。
pub struct FormatRegistry {
    /// 按名称索引的转换器
    converters: HashMap<String, Arc<dyn FormatConverter>>,
    /// 按扩展名索引的转换器名称
    extension_map: HashMap<String, String>,
}

impl FormatRegistry {
    /// 创建空的注册表
    pub fn new() -> Self {
        Self {
            converters: HashMap::new(),
            extension_map: HashMap::new(),
        }
    }

    /// 注册格式转换器
    pub fn register(&mut self, converter: Arc<dyn FormatConverter>) {
        let name = converter.name().to_lowercase();

        // 注册扩展名映射
        for ext in converter.extensions() {
            let ext_lower = ext.to_lowercase();
            self.extension_map.insert(ext_lower, name.clone());
        }

        // 注册转换器
        self.converters.insert(name, converter);
    }

    /// 按名称获取转换器
    pub fn get(&self, name: &str) -> Option<Arc<dyn FormatConverter>> {
        self.converters.get(&name.to_lowercase()).cloned()
    }

    /// 按文件路径自动检测格式并获取转换器
    pub fn detect(&self, path: &Path) -> Option<Arc<dyn FormatConverter>> {
        let path_str = path.to_string_lossy().to_lowercase();

        // 首先尝试完整扩展名匹配（如 .sce.rds）
        // 按扩展名长度降序排序，确保更长的扩展名优先匹配
        let mut extensions: Vec<_> = self.extension_map.iter().collect();
        extensions.sort_by(|a, b| b.0.len().cmp(&a.0.len()));

        for (ext, name) in extensions {
            if path_str.ends_with(ext) {
                if let Some(converter) = self.converters.get(name) {
                    return Some(converter.clone());
                }
            }
        }

        // 然后尝试每个转换器的 detect 方法
        for converter in self.converters.values() {
            if converter.detect(path) {
                return Some(converter.clone());
            }
        }

        None
    }

    /// 列出所有已注册的格式名称
    pub fn list(&self) -> Vec<&str> {
        self.converters.keys().map(|s| s.as_str()).collect()
    }

    /// 获取所有格式的详细信息
    pub fn list_info(&self) -> Vec<FormatInfo> {
        self.converters
            .values()
            .map(|c| FormatInfo::from_converter(c.as_ref()))
            .collect()
    }

    /// 检查格式是否已注册
    pub fn contains(&self, name: &str) -> bool {
        self.converters.contains_key(&name.to_lowercase())
    }

    /// 获取支持读取的格式列表
    pub fn list_readable(&self) -> Vec<&str> {
        self.converters
            .iter()
            .filter(|(_, c)| c.can_read())
            .map(|(name, _)| name.as_str())
            .collect()
    }

    /// 获取支持写入的格式列表
    pub fn list_writable(&self) -> Vec<&str> {
        self.converters
            .iter()
            .filter(|(_, c)| c.can_write())
            .map(|(name, _)| name.as_str())
            .collect()
    }
}

impl Default for FormatRegistry {
    fn default() -> Self {
        Self::new()
    }
}

/// 创建包含所有内置格式的默认注册表
pub fn create_default_registry() -> FormatRegistry {
    let mut registry = FormatRegistry::new();

    // 注册内置格式
    // 注意：SCE 需要在 Seurat 之前注册，因为两者都使用 .rds 扩展名
    // SCE 使用 .sce.rds 作为明确扩展名，并通过内容检测区分
    registry.register(Arc::new(AnndataConverter));
    registry.register(Arc::new(SceConverter));
    registry.register(Arc::new(SeuratConverter));
    registry.register(Arc::new(LoomConverter));

    registry
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_registry_creation() {
        let registry = create_default_registry();

        // 验证内置格式已注册
        assert!(registry.contains("anndata"));
        assert!(registry.contains("seurat"));
        assert!(registry.contains("sce"));
        assert!(registry.contains("loom"));
    }

    #[test]
    fn test_registry_get() {
        let registry = create_default_registry();

        // 按名称获取
        let anndata = registry.get("anndata");
        assert!(anndata.is_some());
        assert_eq!(anndata.unwrap().name(), "anndata");

        // 大小写不敏感
        let anndata_upper = registry.get("ANNDATA");
        assert!(anndata_upper.is_some());

        // 获取 SCE
        let sce = registry.get("sce");
        assert!(sce.is_some());
        assert_eq!(sce.unwrap().name(), "sce");
    }

    #[test]
    fn test_registry_detect() {
        let registry = create_default_registry();

        // 按扩展名检测
        let h5ad = registry.detect(Path::new("data.h5ad"));
        assert!(h5ad.is_some());
        assert_eq!(h5ad.unwrap().name(), "anndata");

        // .sce.rds 应该检测为 SCE
        let sce_rds = registry.detect(Path::new("data.sce.rds"));
        assert!(sce_rds.is_some());
        assert_eq!(sce_rds.unwrap().name(), "sce");

        // 普通 .rds 默认检测为 Seurat（除非内容是 SCE）
        let rds = registry.detect(Path::new("data.rds"));
        assert!(rds.is_some());
        assert_eq!(rds.unwrap().name(), "seurat");

        let loom = registry.detect(Path::new("data.loom"));
        assert!(loom.is_some());
        assert_eq!(loom.unwrap().name(), "loom");
    }

    #[test]
    fn test_registry_list() {
        let registry = create_default_registry();

        let formats = registry.list();
        assert!(formats.contains(&"anndata"));
        assert!(formats.contains(&"seurat"));
        assert!(formats.contains(&"sce"));
        assert!(formats.contains(&"loom"));
    }

    #[test]
    fn test_format_info() {
        let registry = create_default_registry();
        let info_list = registry.list_info();

        assert!(!info_list.is_empty());

        for info in &info_list {
            // 验证格式信息完整
            assert!(!info.name.is_empty());
            assert!(!info.extensions.is_empty());
        }
    }
}
