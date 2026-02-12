//! 格式转换插件系统
//!
//! 提供统一的格式转换接口，支持动态注册和查找格式转换器。
//!
//! ## 设计原则
//!
//! - **开闭原则**：添加新格式不需要修改现有代码
//! - **统一接口**：所有格式转换器实现相同的 trait
//! - **可扩展性**：支持运行时注册新格式
//!
//! ## 使用示例
//!
//! ```rust,ignore
//! use crosscell::formats::{FormatRegistry, create_default_registry};
//!
//! let registry = create_default_registry();
//!
//! // 列出所有支持的格式
//! for format in registry.list() {
//!     println!("{}", format);
//! }
//!
//! // 获取转换器
//! let converter = registry.get("anndata").unwrap();
//! let data = converter.read(Path::new("data.h5ad"))?;
//! ```

mod registry;
mod anndata;
mod seurat;
mod sce;
mod loom;

pub use registry::{FormatRegistry, create_default_registry};
pub use anndata::AnndataConverter;
pub use seurat::SeuratConverter;
pub use sce::SceConverter;
pub use loom::LoomConverter;

use crate::error::CrossCellError;
use crate::ir::SingleCellData;
use std::path::Path;

/// 格式转换器 trait
///
/// 所有格式转换器必须实现此 trait，提供统一的读写接口。
pub trait FormatConverter: Send + Sync {
    /// 格式名称（用于 CLI 参数）
    fn name(&self) -> &str;
    
    /// 格式显示名称（用于用户界面）
    fn display_name(&self) -> &str {
        self.name()
    }
    
    /// 支持的文件扩展名
    fn extensions(&self) -> &[&str];
    
    /// 是否支持读取
    fn can_read(&self) -> bool;
    
    /// 是否支持写入
    fn can_write(&self) -> bool;
    
    /// 格式描述
    fn description(&self) -> &str {
        ""
    }
    
    /// 从文件读取数据
    fn read(&self, path: &Path) -> Result<SingleCellData, CrossCellError>;
    
    /// 将数据写入文件
    fn write(&self, data: &SingleCellData, path: &Path) -> Result<(), CrossCellError>;
    
    /// 检测文件是否为此格式
    fn detect(&self, path: &Path) -> bool {
        if let Some(ext) = path.extension() {
            let ext_str = format!(".{}", ext.to_string_lossy().to_lowercase());
            self.extensions().iter().any(|e| e.to_lowercase() == ext_str)
        } else {
            false
        }
    }
}

/// 格式信息（用于显示）
#[derive(Debug, Clone)]
pub struct FormatInfo {
    pub name: String,
    pub display_name: String,
    pub extensions: Vec<String>,
    pub can_read: bool,
    pub can_write: bool,
    pub description: String,
}

impl FormatInfo {
    /// 从转换器创建格式信息
    pub fn from_converter(converter: &dyn FormatConverter) -> Self {
        Self {
            name: converter.name().to_string(),
            display_name: converter.display_name().to_string(),
            extensions: converter.extensions().iter().map(|s| s.to_string()).collect(),
            can_read: converter.can_read(),
            can_write: converter.can_write(),
            description: converter.description().to_string(),
        }
    }
    
    /// 格式化为显示字符串
    pub fn format_display(&self) -> String {
        let exts = self.extensions.join(", ");
        let caps = match (self.can_read, self.can_write) {
            (true, true) => "read/write",
            (true, false) => "read only",
            (false, true) => "write only",
            (false, false) => "none",
        };
        format!("{:<12} ({})  [{}]", self.name, exts, caps)
    }
}
