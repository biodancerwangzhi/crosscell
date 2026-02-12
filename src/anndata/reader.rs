//! AnnData (.h5ad) 文件读取器
//!
//! 从 HDF5 格式的 AnnData 文件读取数据并转换为 IR。
//! 
//! 支持两种读取模式：
//! - 完整读取：读取所有数据（默认）
//! - 部分读取：仅读取元数据，不加载表达矩阵（用于大数据集预览）

use super::{AnnDataError, Result};
use crate::ir::{
    DatasetMetadata, ExpressionMatrix, SparseMatrixCSR, SparseMatrixCSC, DenseMatrix, DataFrame, SingleCellData,
    Embedding, PairwiseMatrix, LazyMatrix, BackendType, SparseFormat,
};
use hdf5::{File, Group, Dataset};
use std::path::{Path, PathBuf};
use std::collections::HashMap;
use arrow::array::{
    ArrayRef, Float64Array, Int64Array, StringArray, BooleanArray, DictionaryArray, Int32Array,
};
use arrow::datatypes::Int32Type;
use std::sync::Arc;

/// 部分加载选项
///
/// 用于控制从 .h5ad 文件加载哪些数据组件。
/// 借鉴 scDIOR 的设计，支持仅加载元数据而不加载表达矩阵。
#[derive(Debug, Clone)]
pub struct PartialLoadOptions {
    /// 是否加载表达矩阵（/X）
    pub load_expression: bool,
    /// 是否加载细胞元数据（/obs）
    pub load_cell_metadata: bool,
    /// 是否加载基因元数据（/var）
    pub load_gene_metadata: bool,
    /// 是否加载降维嵌入（/obsm）
    pub load_embeddings: bool,
    /// 是否加载多层表达矩阵（/layers）
    pub load_layers: bool,
    /// 是否加载成对矩阵（/obsp, /varp）
    pub load_pairwise: bool,
    /// 是否加载空间数据
    pub load_spatial: bool,
    /// 是否使用延迟加载（不立即加载数据，只记录元数据）
    pub lazy_load: bool,
}

impl Default for PartialLoadOptions {
    fn default() -> Self {
        Self {
            load_expression: true,
            load_cell_metadata: true,
            load_gene_metadata: true,
            load_embeddings: true,
            load_layers: true,
            load_pairwise: true,
            load_spatial: true,
            lazy_load: false,
        }
    }
}

impl PartialLoadOptions {
    /// 创建仅加载元数据的选项（不加载表达矩阵）
    pub fn metadata_only() -> Self {
        Self {
            load_expression: false,
            load_cell_metadata: true,
            load_gene_metadata: true,
            load_embeddings: false,
            load_layers: false,
            load_pairwise: false,
            load_spatial: false,
            lazy_load: false,
        }
    }

    /// 创建延迟加载选项（记录元数据但不立即加载数据）
    pub fn lazy() -> Self {
        Self {
            load_expression: true,
            load_cell_metadata: true,
            load_gene_metadata: true,
            load_embeddings: true,
            load_layers: true,
            load_pairwise: true,
            load_spatial: true,
            lazy_load: true,
        }
    }

    /// 创建完整加载选项
    pub fn full() -> Self {
        Self::default()
    }
}

/// H5AD 文件信息（用于快速预览）
#[derive(Debug, Clone)]
pub struct H5adInfo {
    /// 细胞数量
    pub n_cells: usize,
    /// 基因数量
    pub n_genes: usize,
    /// 表达矩阵是否为稀疏格式
    pub is_sparse: bool,
    /// 稀疏格式类型（如果是稀疏矩阵）
    pub sparse_format: Option<String>,
    /// 非零元素数量（如果是稀疏矩阵）
    pub nnz: Option<usize>,
    /// 估算的矩阵内存大小（字节）
    pub estimated_matrix_size: usize,
    /// 细胞元数据列数
    pub n_obs_columns: usize,
    /// 基因元数据列数
    pub n_var_columns: usize,
    /// 降维嵌入名称列表
    pub embedding_names: Vec<String>,
    /// 层名称列表
    pub layer_names: Vec<String>,
    /// 成对矩阵名称列表（obsp）
    pub obsp_names: Vec<String>,
    /// 成对矩阵名称列表（varp）
    pub varp_names: Vec<String>,
    /// 是否包含空间数据
    pub has_spatial: bool,
    /// 文件大小（字节）
    pub file_size: u64,
}

impl H5adInfo {
    /// 估算完整加载所需的内存（字节）
    pub fn estimate_full_memory(&self) -> usize {
        let mut total = self.estimated_matrix_size;
        
        // 元数据估算（每列约 8 字节 * 行数）
        total += self.n_obs_columns * self.n_cells * 8;
        total += self.n_var_columns * self.n_genes * 8;
        
        // 嵌入估算（假设每个嵌入平均 50 维）
        total += self.embedding_names.len() * self.n_cells * 50 * 8;
        
        total
    }

    /// 格式化显示
    pub fn display(&self) -> String {
        let mut output = String::new();
        
        output.push_str(&format!("📊 H5AD File Summary\n"));
        output.push_str(&format!("   Cells: {}\n", self.n_cells));
        output.push_str(&format!("   Genes: {}\n", self.n_genes));
        
        if self.is_sparse {
            let format_str = self.sparse_format.as_deref().unwrap_or("unknown");
            let nnz_str = self.nnz.map(|n| format!("{}", n)).unwrap_or_else(|| "unknown".to_string());
            let sparsity = if let Some(nnz) = self.nnz {
                let total = self.n_cells * self.n_genes;
                if total > 0 {
                    format!("{:.2}%", (1.0 - nnz as f64 / total as f64) * 100.0)
                } else {
                    "N/A".to_string()
                }
            } else {
                "unknown".to_string()
            };
            output.push_str(&format!("   Matrix: Sparse ({}, {} nnz, {} sparse)\n", format_str, nnz_str, sparsity));
        } else {
            output.push_str(&format!("   Matrix: Dense\n"));
        }
        
        output.push_str(&format!("   Estimated matrix size: {:.2} MB\n", self.estimated_matrix_size as f64 / 1024.0 / 1024.0));
        output.push_str(&format!("   Cell metadata: {} columns\n", self.n_obs_columns));
        output.push_str(&format!("   Gene metadata: {} columns\n", self.n_var_columns));
        
        if !self.embedding_names.is_empty() {
            output.push_str(&format!("   Embeddings: {}\n", self.embedding_names.join(", ")));
        }
        
        if !self.layer_names.is_empty() {
            output.push_str(&format!("   Layers: {}\n", self.layer_names.join(", ")));
        }
        
        if self.has_spatial {
            output.push_str(&format!("   Spatial: Yes\n"));
        }
        
        output.push_str(&format!("   File size: {:.2} MB\n", self.file_size as f64 / 1024.0 / 1024.0));
        output.push_str(&format!("   Estimated full load: {:.2} GB\n", self.estimate_full_memory() as f64 / 1024.0 / 1024.0 / 1024.0));
        
        output
    }
}

/// 快速获取 H5AD 文件信息（不加载数据）
///
/// 借鉴 scDIOR 的设计，仅读取元数据以快速预览文件内容。
///
/// # 参数
/// - `path`: .h5ad 文件路径
///
/// # 返回
/// - `Ok(H5adInfo)`: 文件信息
/// - `Err(AnnDataError)`: 读取失败
pub fn inspect_h5ad<P: AsRef<Path>>(path: P) -> Result<H5adInfo> {
    let path_ref = path.as_ref();
    let file = File::open(path_ref)?;
    
    // 获取文件大小
    let file_size = std::fs::metadata(path_ref)
        .map(|m| m.len())
        .unwrap_or(0);
    
    // 读取表达矩阵元数据
    let (n_cells, n_genes, is_sparse, sparse_format, nnz) = read_expression_metadata(&file)?;
    
    // 估算矩阵内存大小
    let estimated_matrix_size = if is_sparse {
        // 稀疏矩阵：data + indices + indptr
        let estimated_nnz = nnz.unwrap_or((n_cells * n_genes) / 20);
        estimated_nnz * (8 + 8) + (n_cells + 1) * 8
    } else {
        n_cells * n_genes * 8
    };
    
    // 读取元数据列数
    let n_obs_columns = if file.link_exists("obs") {
        count_group_members(&file, "obs")?
    } else {
        0
    };
    
    let n_var_columns = if file.link_exists("var") {
        count_group_members(&file, "var")?
    } else {
        0
    };
    
    // 读取嵌入名称
    let embedding_names = if file.link_exists("obsm") {
        list_group_members(&file, "obsm")?
    } else {
        Vec::new()
    };
    
    // 读取层名称
    let layer_names = if file.link_exists("layers") {
        list_group_members(&file, "layers")?
    } else {
        Vec::new()
    };
    
    // 读取成对矩阵名称
    let obsp_names = if file.link_exists("obsp") {
        list_group_members(&file, "obsp")?
    } else {
        Vec::new()
    };
    
    let varp_names = if file.link_exists("varp") {
        list_group_members(&file, "varp")?
    } else {
        Vec::new()
    };
    
    // 检查空间数据
    let has_spatial = embedding_names.contains(&"spatial".to_string());
    
    Ok(H5adInfo {
        n_cells,
        n_genes,
        is_sparse,
        sparse_format,
        nnz,
        estimated_matrix_size,
        n_obs_columns,
        n_var_columns,
        embedding_names,
        layer_names,
        obsp_names,
        varp_names,
        has_spatial,
        file_size,
    })
}

/// 读取表达矩阵元数据（不加载数据）
fn read_expression_metadata(file: &File) -> Result<(usize, usize, bool, Option<String>, Option<usize>)> {
    if !file.link_exists("X") {
        return Err(AnnDataError::MissingField("/X".to_string()));
    }
    
    // 首先尝试作为 Dataset 打开（稠密矩阵或某些旧版本的稀疏矩阵）
    // 这样可以避免在 Group 不存在时的错误
    if let Ok(x_dataset) = file.dataset("X") {
        let shape = x_dataset.shape();
        if shape.len() != 2 {
            return Err(AnnDataError::InvalidFormat(format!(
                "Expected 2D matrix, got {}D",
                shape.len()
            )));
        }
        let (n_rows, n_cols) = (shape[0] as usize, shape[1] as usize);
        return Ok((n_rows, n_cols, false, None, None));
    }
    
    // 尝试作为 Group 打开（稀疏矩阵）
    if let Ok(x_group) = file.group("X") {
        // 尝试读取 shape 属性（新版本格式）
        // 或 h5sparse_shape 属性（旧版本格式）
        let shape: Vec<usize> = if let Ok(shape_attr) = x_group.attr("shape") {
            shape_attr.read_1d()?.to_vec()
        } else if let Ok(shape_attr) = x_group.attr("h5sparse_shape") {
            // 旧版本格式使用 h5sparse_shape
            shape_attr.read_1d()?.to_vec()
        } else {
            return Err(AnnDataError::MissingField("/X shape attribute".to_string()));
        };
        
        if shape.len() != 2 {
            return Err(AnnDataError::InvalidFormat(format!(
                "Expected shape to have 2 dimensions, got {}",
                shape.len()
            )));
        }
        
        let (n_rows, n_cols) = (shape[0], shape[1]);
        
        // 读取 encoding-type 属性来确定稀疏格式（新版本）
        // 或 h5sparse_format 属性（旧版本）
        let sparse_format = if let Ok(encoding_attr) = x_group.attr("encoding-type") {
            if let Ok(encoding) = encoding_attr.read_scalar::<hdf5::types::VarLenUnicode>() {
                Some(encoding.to_string())
            } else {
                Some("csr_matrix".to_string())
            }
        } else if let Ok(format_attr) = x_group.attr("h5sparse_format") {
            // 旧版本格式使用 h5sparse_format
            if let Ok(format_str) = format_attr.read_scalar::<hdf5::types::VarLenUnicode>() {
                Some(format!("{}_matrix", format_str))
            } else {
                Some("csr_matrix".to_string())
            }
        } else {
            Some("csr_matrix".to_string())
        };
        
        // 读取 nnz
        let nnz = if let Ok(data_dataset) = x_group.dataset("data") {
            Some(data_dataset.shape()[0] as usize)
        } else {
            None
        };
        
        return Ok((n_rows, n_cols, true, sparse_format, nnz));
    }
    
    Err(AnnDataError::InvalidFormat(
        "/X is neither a Group (sparse) nor a Dataset (dense)".to_string(),
    ))
}

/// 统计 Group 中的成员数量
fn count_group_members(file: &File, group_name: &str) -> Result<usize> {
    // 首先检查是否是 Group
    if let Ok(group) = file.group(group_name) {
        let members = group.member_names()?;
        // 排除特殊成员（如 __categories）
        Ok(members.iter().filter(|n| !n.starts_with("__")).count())
    } else {
        // 如果不是 Group，可能是旧版本格式的 Dataset
        // 在这种情况下返回 0
        Ok(0)
    }
}

/// 列出 Group 中的成员名称
fn list_group_members(file: &File, group_name: &str) -> Result<Vec<String>> {
    // 首先检查是否是 Group
    if let Ok(group) = file.group(group_name) {
        let members = group.member_names()?;
        // 排除特殊成员（如 __categories）
        Ok(members.into_iter().filter(|n| !n.starts_with("__")).collect())
    } else {
        // 如果不是 Group，可能是旧版本格式的 Dataset
        // 在这种情况下返回空列表
        Ok(Vec::new())
    }
}

/// 使用部分加载选项读取 .h5ad 文件
///
/// # 参数
/// - `path`: .h5ad 文件路径
/// - `options`: 部分加载选项
///
/// # 返回
/// - `Ok(SingleCellData)`: 成功读取的 IR 数据
/// - `Err(AnnDataError)`: 读取失败
pub fn read_h5ad_partial<P: AsRef<Path>>(path: P, options: &PartialLoadOptions) -> Result<SingleCellData> {
    let path_ref = path.as_ref();
    let file = File::open(path_ref)?;
    
    // 首先读取表达矩阵元数据以获取维度
    let (n_cells, n_genes, is_sparse, sparse_format_str, nnz) = read_expression_metadata(&file)?;
    
    // 根据选项读取表达矩阵
    let expression = if options.load_expression {
        if options.lazy_load {
            // 延迟加载：只记录元数据
            let sparse_format = if is_sparse {
                if sparse_format_str.as_deref() == Some("csc_matrix") {
                    Some(SparseFormat::CSC)
                } else {
                    Some(SparseFormat::CSR)
                }
            } else {
                None
            };
            
            let mut lazy = LazyMatrix::new(
                BackendType::HDF5,
                PathBuf::from(path_ref),
                "/X".to_string(),
                (n_cells, n_genes),
                is_sparse,
                sparse_format,
            );
            
            if let Some(n) = nnz {
                lazy = lazy.with_nnz(n);
            }
            
            ExpressionMatrix::Lazy(lazy)
        } else {
            // 完整加载
            read_expression_matrix(&file)?
        }
    } else {
        // 不加载表达矩阵，创建空的延迟加载矩阵
        let sparse_format = if is_sparse {
            if sparse_format_str.as_deref() == Some("csc_matrix") {
                Some(SparseFormat::CSC)
            } else {
                Some(SparseFormat::CSR)
            }
        } else {
            None
        };
        
        let mut lazy = LazyMatrix::new(
            BackendType::HDF5,
            PathBuf::from(path_ref),
            "/X".to_string(),
            (n_cells, n_genes),
            is_sparse,
            sparse_format,
        );
        
        if let Some(n) = nnz {
            lazy = lazy.with_nnz(n);
        }
        
        ExpressionMatrix::Lazy(lazy)
    };
    
    // 读取细胞元数据
    let cell_metadata = if options.load_cell_metadata && file.link_exists("obs") {
        read_metadata_dataframe(&file, "obs", n_cells)?
    } else {
        DataFrame::empty(n_cells)
    };
    
    // 读取基因元数据
    let gene_metadata = if options.load_gene_metadata && file.link_exists("var") {
        read_metadata_dataframe(&file, "var", n_genes)?
    } else {
        DataFrame::empty(n_genes)
    };
    
    // 读取降维嵌入
    let embeddings = if options.load_embeddings && file.link_exists("obsm") {
        Some(read_embeddings(&file, n_cells)?)
    } else {
        None
    };
    
    // 读取 layers
    let layers = if options.load_layers && file.link_exists("layers") {
        Some(read_layers(&file, n_cells, n_genes)?)
    } else {
        None
    };
    
    // 读取成对矩阵
    let cell_pairwise = if options.load_pairwise && file.link_exists("obsp") {
        Some(read_pairwise_matrices(&file, "obsp", n_cells)?)
    } else {
        None
    };
    
    let gene_pairwise = if options.load_pairwise && file.link_exists("varp") {
        Some(read_pairwise_matrices(&file, "varp", n_genes)?)
    } else {
        None
    };
    
    // 读取空间数据
    let spatial = if options.load_spatial {
        read_spatial_data(&file, n_cells)?
    } else {
        None
    };
    
    // 创建数据集元数据
    let metadata = DatasetMetadata::new(n_cells, n_genes, "anndata".to_string());
    
    // 构造 SingleCellData
    let mut data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata)
        .map_err(|e| AnnDataError::InvalidFormat(e))?;
    
    data.embeddings = embeddings;
    data.layers = layers;
    data.cell_pairwise = cell_pairwise;
    data.gene_pairwise = gene_pairwise;
    data.spatial = spatial;
    
    Ok(data)
}

/// 从 .h5ad 文件读取 AnnData 对象
///
/// # 参数
/// - `path`: .h5ad 文件路径
///
/// # 返回
/// - `Ok(SingleCellData)`: 成功读取的 IR 数据
/// - `Err(AnnDataError)`: 读取失败
///
/// # 示例
/// ```no_run
/// use crosscell::anndata::read_h5ad;
///
/// let data = read_h5ad("pbmc3k.h5ad").unwrap();
/// println!("Loaded {} cells × {} genes", data.metadata.n_cells, data.metadata.n_genes);
/// ```
pub fn read_h5ad<P: AsRef<Path>>(path: P) -> Result<SingleCellData> {
    let file = File::open(path.as_ref())?;

    // 读取主表达矩阵 /X
    let expression = read_expression_matrix(&file)?;
    let (n_cells, n_genes) = expression.shape();

    // 读取细胞元数据 /obs
    let cell_metadata = if file.link_exists("obs") {
        read_metadata_dataframe(&file, "obs", n_cells)?
    } else {
        DataFrame::empty(n_cells)
    };

    // 读取基因元数据 /var
    let gene_metadata = if file.link_exists("var") {
        read_metadata_dataframe(&file, "var", n_genes)?
    } else {
        DataFrame::empty(n_genes)
    };

    // 读取降维嵌入 /obsm
    let embeddings = if file.link_exists("obsm") {
        Some(read_embeddings(&file, n_cells)?)
    } else {
        None
    };

    // 读取 layers /layers
    let layers = if file.link_exists("layers") {
        Some(read_layers(&file, n_cells, n_genes)?)
    } else {
        None
    };

    // 读取细胞-细胞成对矩阵 /obsp
    let cell_pairwise = if file.link_exists("obsp") {
        Some(read_pairwise_matrices(&file, "obsp", n_cells)?)
    } else {
        None
    };

    // 读取基因-基因成对矩阵 /varp
    let gene_pairwise = if file.link_exists("varp") {
        Some(read_pairwise_matrices(&file, "varp", n_genes)?)
    } else {
        None
    };

    // 读取空间数据 /obsm['spatial'] 和 /uns['spatial']
    let spatial = read_spatial_data(&file, n_cells)?;

    // 创建数据集元数据
    let metadata = DatasetMetadata::new(n_cells, n_genes, "anndata".to_string());

    // 构造 SingleCellData
    let mut data = SingleCellData::new(expression, cell_metadata, gene_metadata, metadata)
        .map_err(|e| AnnDataError::InvalidFormat(e))?;
    
    // 设置可选字段
    data.embeddings = embeddings;
    data.layers = layers;
    data.cell_pairwise = cell_pairwise;
    data.gene_pairwise = gene_pairwise;
    data.spatial = spatial;
    
    Ok(data)
}

/// 从 HDF5 文件读取表达矩阵 (/X)
///
/// 支持两种格式：
/// 1. 稀疏矩阵（CSR 格式）：/X 是一个 Group，包含 data, indices, indptr
/// 2. 稠密矩阵：/X 是一个 Dataset
fn read_expression_matrix(file: &File) -> Result<ExpressionMatrix> {
    // 检查 /X 是否存在
    if !file.link_exists("X") {
        return Err(AnnDataError::MissingField("/X".to_string()));
    }

    // 尝试作为 Group 打开（稀疏矩阵）
    if let Ok(x_group) = file.group("X") {
        return read_sparse_csr_matrix(&x_group);
    }

    // 尝试作为 Dataset 打开（稠密矩阵）
    if let Ok(x_dataset) = file.dataset("X") {
        return read_dense_matrix(&x_dataset);
    }

    Err(AnnDataError::InvalidFormat(
        "/X is neither a Group (sparse) nor a Dataset (dense)".to_string(),
    ))
}

/// 读取稀疏矩阵（CSR 或 CSC）
///
/// 支持两种 HDF5 属性格式：
/// - 新版本 (anndata >= 0.8): `shape`, `encoding-type` (csr_matrix/csc_matrix)
/// - 旧版本 (anndata < 0.8): `h5sparse_shape`, `h5sparse_format` (csr/csc)
///
/// HDF5 结构：
/// - /X/data: 非零元素值 (float32/float64)
/// - /X/indices: 索引 (int32 或 int64)
/// - /X/indptr: 指针 (int32 或 int64)
fn read_sparse_csr_matrix(group: &Group) -> Result<ExpressionMatrix> {
    // 读取 shape 属性（兼容新旧格式）
    let shape: Vec<usize> = if let Ok(shape_attr) = group.attr("shape") {
        shape_attr.read_1d()?.to_vec()
    } else if let Ok(shape_attr) = group.attr("h5sparse_shape") {
        // 旧版本格式 (anndata < 0.8)
        shape_attr.read_1d()?.to_vec()
    } else {
        return Err(AnnDataError::MissingField("/X shape attribute".to_string()));
    };
    
    if shape.len() != 2 {
        return Err(AnnDataError::InvalidFormat(format!(
            "Expected shape to have 2 dimensions, got {}",
            shape.len()
        )));
    }
    let (n_rows, n_cols) = (shape[0], shape[1]);

    // 判断稀疏格式（CSR 或 CSC）
    // 新版本: encoding-type = "csr_matrix" / "csc_matrix"
    // 旧版本: h5sparse_format = "csr" / "csc"
    let is_csc = if let Ok(encoding_attr) = group.attr("encoding-type") {
        if let Ok(encoding) = encoding_attr.read_scalar::<hdf5::types::VarLenUnicode>() {
            encoding.to_string().contains("csc")
        } else {
            false
        }
    } else if let Ok(format_attr) = group.attr("h5sparse_format") {
        if let Ok(format_str) = format_attr.read_scalar::<hdf5::types::VarLenUnicode>() {
            format_str.to_string() == "csc"
        } else if let Ok(format_bytes) = format_attr.read_scalar::<hdf5::types::FixedAscii<16>>() {
            format_bytes.to_string().trim() == "csc"
        } else {
            false
        }
    } else {
        false // 默认 CSR
    };

    // 读取 data (非零元素值) — 支持 float32 和 float64
    let data_dataset = group
        .dataset("data")
        .map_err(|_| AnnDataError::MissingField("/X/data".to_string()))?;
    let data: Vec<f64> = if let Ok(d) = data_dataset.read_1d::<f64>() {
        d.to_vec()
    } else if let Ok(d) = data_dataset.read_1d::<f32>() {
        d.to_vec().into_iter().map(|x| x as f64).collect()
    } else {
        return Err(AnnDataError::UnsupportedType(
            "X/data must be float32 or float64".to_string(),
        ));
    };

    // 读取 indices
    let indices_dataset = group
        .dataset("indices")
        .map_err(|_| AnnDataError::MissingField("/X/indices".to_string()))?;
    
    let indices: Vec<usize> = if let Ok(indices_i32) = indices_dataset.read_1d::<i32>() {
        indices_i32.to_vec().into_iter().map(|x| x as usize).collect()
    } else if let Ok(indices_i64) = indices_dataset.read_1d::<i64>() {
        indices_i64.to_vec().into_iter().map(|x| x as usize).collect()
    } else {
        return Err(AnnDataError::UnsupportedType(
            "indices must be int32 or int64".to_string(),
        ));
    };

    // 读取 indptr
    let indptr_dataset = group
        .dataset("indptr")
        .map_err(|_| AnnDataError::MissingField("/X/indptr".to_string()))?;
    
    let indptr: Vec<usize> = if let Ok(indptr_i32) = indptr_dataset.read_1d::<i32>() {
        indptr_i32.to_vec().into_iter().map(|x| x as usize).collect()
    } else if let Ok(indptr_i64) = indptr_dataset.read_1d::<i64>() {
        indptr_i64.to_vec().into_iter().map(|x| x as usize).collect()
    } else {
        return Err(AnnDataError::UnsupportedType(
            "indptr must be int32 or int64".to_string(),
        ));
    };

    if data.len() != indices.len() {
        return Err(AnnDataError::DimensionMismatch(format!(
            "data length {} != indices length {}",
            data.len(),
            indices.len()
        )));
    }

    if is_csc {
        // CSC: indptr 长度 = n_cols + 1, indices 是行索引
        if indptr.len() != n_cols + 1 {
            return Err(AnnDataError::DimensionMismatch(format!(
                "CSC indptr length {} != n_cols + 1 = {}",
                indptr.len(),
                n_cols + 1
            )));
        }
        let csc = SparseMatrixCSC::new(data, indices, indptr, n_rows, n_cols)
            .map_err(|e| AnnDataError::InvalidFormat(e))?;
        Ok(ExpressionMatrix::SparseCSC(csc))
    } else {
        // CSR: indptr 长度 = n_rows + 1, indices 是列索引
        if indptr.len() != n_rows + 1 {
            return Err(AnnDataError::DimensionMismatch(format!(
                "CSR indptr length {} != n_rows + 1 = {}",
                indptr.len(),
                n_rows + 1
            )));
        }
        let csr = SparseMatrixCSR::new(data, indices, indptr, n_rows, n_cols)
            .map_err(|e| AnnDataError::InvalidFormat(e))?;
        Ok(ExpressionMatrix::SparseCSR(csr))
    }
}

/// 读取稠密矩阵
///
/// HDF5 结构：
/// - /X: Dataset (n_cells × n_genes)
fn read_dense_matrix(dataset: &Dataset) -> Result<ExpressionMatrix> {
    // 获取矩阵形状
    let shape = dataset.shape();
    if shape.len() != 2 {
        return Err(AnnDataError::InvalidFormat(format!(
            "Expected 2D matrix, got {}D",
            shape.len()
        )));
    }
    let (n_rows, n_cols) = (shape[0] as usize, shape[1] as usize);

    // 读取数据（行优先）— 支持 float64 和 float32
    let data: Vec<f64> = if let Ok(arr) = dataset.read_2d::<f64>() {
        arr.into_raw_vec()
    } else if let Ok(arr) = dataset.read_2d::<f32>() {
        arr.into_raw_vec().into_iter().map(|x| x as f64).collect()
    } else {
        return Err(AnnDataError::UnsupportedType(
            "Dense matrix X must be float32 or float64".to_string(),
        ));
    };

    // 创建稠密矩阵
    let dense = DenseMatrix::new(data, n_rows, n_cols)
        .map_err(|e| AnnDataError::InvalidFormat(e))?;

    Ok(ExpressionMatrix::Dense(dense))
}

/// 读取元数据 DataFrame (/obs 或 /var)
///
/// HDF5 结构：
/// - /obs 或 /var: Group，包含多个列
/// - 每列是一个 Dataset
/// - Categorical 列：
///   - 列本身存储整数编码
///   - /obs/__categories/{column_name}: 存储类别标签
fn read_metadata_dataframe(file: &File, group_name: &str, expected_rows: usize) -> Result<DataFrame> {
    // 先尝试作为 Group 打开（新版本 anndata >= 0.7）
    // 如果失败，尝试作为 Dataset 打开（旧版本 compound Dataset）
    let group = match file.group(group_name) {
        Ok(g) => g,
        Err(_) => {
            // 可能是旧版本 compound Dataset 格式 (anndata < 0.7)
            if let Ok(dataset) = file.dataset(group_name) {
                log::info!("Reading legacy compound Dataset for /{}", group_name);
                return read_metadata_dataframe_legacy(&dataset, group_name, expected_rows);
            }
            return Err(AnnDataError::MissingField(format!("/{}", group_name)));
        }
    };
    
    // 获取所有列名（排除 __categories 和其他非 Dataset 成员）
    let all_members = group.member_names()?;
    let mut column_names: Vec<String> = Vec::new();
    
    for name in all_members {
        // 跳过特殊组（旧版本格式）
        if name == "__categories" {
            continue;
        }
        
        // 检查是否是 Dataset 或 Group（新版本 categorical）
        if group.dataset(&name).is_ok() {
            column_names.push(name);
        } else if let Ok(subgroup) = group.group(&name) {
            // 检查是否是 categorical 列（新版本格式）
            if subgroup.link_exists("codes") && subgroup.link_exists("categories") {
                column_names.push(name);
            }
            // 检查是否是 nullable 列（nullable-integer, nullable-boolean）
            else if subgroup.link_exists("values") && subgroup.link_exists("mask") {
                column_names.push(name);
            }
        }
    }
    
    // 按字母顺序排序以保证一致性
    column_names.sort();
    
    let mut columns = Vec::new();
    let mut data: Vec<ArrayRef> = Vec::new();
    
    for col_name in &column_names {
        // 尝试作为 Dataset 打开
        if let Ok(dataset) = group.dataset(col_name) {
            // 检查是否是 Categorical 列（旧版本格式）
            let is_categorical_old = if let Ok(categories_group) = group.group("__categories") {
                categories_group.link_exists(col_name)
            } else {
                false
            };
            
            if is_categorical_old {
                // 旧版本 Categorical 列
                let array = read_categorical_column(&group, col_name, expected_rows)?;
                columns.push(col_name.clone());
                data.push(array);
            } else {
                // 普通列
                let array = read_column(&dataset, expected_rows)?;
                columns.push(col_name.clone());
                data.push(array);
            }
        } else if let Ok(subgroup) = group.group(col_name) {
            // 新版本 Categorical 列（Group 格式）
            if subgroup.link_exists("codes") && subgroup.link_exists("categories") {
                let array = read_categorical_column_v2(&subgroup, col_name, expected_rows)?;
                columns.push(col_name.clone());
                data.push(array);
            }
            // Nullable 列（nullable-integer, nullable-boolean）
            else if subgroup.link_exists("values") && subgroup.link_exists("mask") {
                let array = read_nullable_column(&subgroup, col_name, expected_rows)?;
                columns.push(col_name.clone());
                data.push(array);
            }
        }
    }
    
    DataFrame::new(columns, data, expected_rows)
        .map_err(|e| AnnDataError::InvalidFormat(e))
}

/// 读取旧版本 compound Dataset 格式的元数据 (anndata < 0.7)
///
/// 旧版本 anndata 将 obs/var 存储为 HDF5 compound Dataset，
/// 每个字段对应一列。字段类型可能是：
/// - 固定长度字节字符串 (S-type) → 需要转换为 UTF-8
/// - float64, int64, int32, bool 等数值类型
///
/// 参考: anndata Python 实现的 read_dataframe_legacy()
fn read_metadata_dataframe_legacy(dataset: &Dataset, group_name: &str, expected_rows: usize) -> Result<DataFrame> {
    use hdf5_sys::h5t::{
        H5Tget_nmembers, H5Tget_member_name, H5Tget_member_type, H5Tget_member_offset,
        H5Tget_class, H5Tget_size, H5Tclose,
        H5T_class_t,
    };
    use hdf5_sys::h5d::H5Dread;
    use hdf5_sys::h5s::H5S_ALL;
    use hdf5_sys::h5p::H5P_DEFAULT;
    use std::ffi::CStr;

    let shape = dataset.shape();
    if shape.is_empty() {
        return Ok(DataFrame::empty(0));
    }
    let n_rows = shape[0];
    
    if n_rows != expected_rows {
        return Err(AnnDataError::DimensionMismatch(format!(
            "/{} has {} rows, expected {}", group_name, n_rows, expected_rows
        )));
    }

    let ds_id = dataset.id();
    let dtype_id = unsafe { hdf5_sys::h5d::H5Dget_type(ds_id) };
    if dtype_id < 0 {
        return Err(AnnDataError::InvalidFormat(format!(
            "Failed to get dtype for /{}", group_name
        )));
    }

    let n_fields = unsafe { H5Tget_nmembers(dtype_id) } as usize;
    let total_size = unsafe { H5Tget_size(dtype_id) };

    // 读取整个 compound dataset 为原始字节
    let buf_size = n_rows * total_size;
    let mut raw_buf: Vec<u8> = vec![0u8; buf_size];
    
    let read_status = unsafe {
        H5Dread(
            ds_id,
            dtype_id,
            H5S_ALL,
            H5S_ALL,
            H5P_DEFAULT,
            raw_buf.as_mut_ptr() as *mut std::ffi::c_void,
        )
    };
    
    if read_status < 0 {
        unsafe { H5Tclose(dtype_id); }
        return Err(AnnDataError::InvalidFormat(format!(
            "Failed to read compound dataset /{}", group_name
        )));
    }

    let mut columns = Vec::new();
    let mut data: Vec<ArrayRef> = Vec::new();
    
    // 跳过第一个字段（通常是 index）
    let start_field = if n_fields > 1 { 1 } else { 0 };

    for i in start_field..n_fields {
        let field_name_ptr = unsafe { H5Tget_member_name(dtype_id, i as u32) };
        if field_name_ptr.is_null() {
            continue;
        }
        let field_name = unsafe { CStr::from_ptr(field_name_ptr) }
            .to_string_lossy()
            .to_string();
        // Free the HDF5-allocated string
        unsafe { hdf5_sys::h5::H5free_memory(field_name_ptr as *mut std::ffi::c_void); }

        let member_type_id = unsafe { H5Tget_member_type(dtype_id, i as u32) };
        let member_offset = unsafe { H5Tget_member_offset(dtype_id, i as u32) };
        let member_class = unsafe { H5Tget_class(member_type_id) };
        let member_size = unsafe { H5Tget_size(member_type_id) };

        let array: Option<ArrayRef> = match member_class {
            // H5T_FLOAT = 1
            H5T_class_t::H5T_FLOAT => {
                if member_size == 8 {
                    // float64
                    let values: Vec<f64> = (0..n_rows)
                        .map(|row| {
                            let offset = row * total_size + member_offset;
                            let bytes: [u8; 8] = raw_buf[offset..offset + 8].try_into().unwrap();
                            f64::from_ne_bytes(bytes)
                        })
                        .collect();
                    Some(Arc::new(Float64Array::from(values)) as ArrayRef)
                } else if member_size == 4 {
                    // float32 → promote to float64
                    let values: Vec<f64> = (0..n_rows)
                        .map(|row| {
                            let offset = row * total_size + member_offset;
                            let bytes: [u8; 4] = raw_buf[offset..offset + 4].try_into().unwrap();
                            f32::from_ne_bytes(bytes) as f64
                        })
                        .collect();
                    Some(Arc::new(Float64Array::from(values)) as ArrayRef)
                } else {
                    log::warn!("Unsupported float size {} for field '{}' in /{}", member_size, field_name, group_name);
                    None
                }
            }
            // H5T_INTEGER = 0
            H5T_class_t::H5T_INTEGER => {
                if member_size == 8 {
                    let values: Vec<i64> = (0..n_rows)
                        .map(|row| {
                            let offset = row * total_size + member_offset;
                            let bytes: [u8; 8] = raw_buf[offset..offset + 8].try_into().unwrap();
                            i64::from_ne_bytes(bytes)
                        })
                        .collect();
                    Some(Arc::new(Int64Array::from(values)) as ArrayRef)
                } else if member_size == 4 {
                    let values: Vec<i64> = (0..n_rows)
                        .map(|row| {
                            let offset = row * total_size + member_offset;
                            let bytes: [u8; 4] = raw_buf[offset..offset + 4].try_into().unwrap();
                            i32::from_ne_bytes(bytes) as i64
                        })
                        .collect();
                    Some(Arc::new(Int64Array::from(values)) as ArrayRef)
                } else if member_size == 1 {
                    // int8 or bool-like
                    let values: Vec<i64> = (0..n_rows)
                        .map(|row| {
                            let offset = row * total_size + member_offset;
                            raw_buf[offset] as i8 as i64
                        })
                        .collect();
                    Some(Arc::new(Int64Array::from(values)) as ArrayRef)
                } else {
                    log::warn!("Unsupported integer size {} for field '{}' in /{}", member_size, field_name, group_name);
                    None
                }
            }
            // H5T_STRING = 3
            H5T_class_t::H5T_STRING => {
                // Fixed-length byte strings (S-type) — convert to UTF-8
                let values: Vec<String> = (0..n_rows)
                    .map(|row| {
                        let offset = row * total_size + member_offset;
                        let end = offset + member_size;
                        let bytes = &raw_buf[offset..end];
                        // Trim trailing null bytes and convert to UTF-8
                        let trimmed = match bytes.iter().position(|&b| b == 0) {
                            Some(pos) => &bytes[..pos],
                            None => bytes,
                        };
                        String::from_utf8_lossy(trimmed).to_string()
                    })
                    .collect();
                Some(Arc::new(StringArray::from(values)) as ArrayRef)
            }
            // H5T_ENUM = 8 (boolean in HDF5)
            H5T_class_t::H5T_ENUM => {
                if member_size == 1 {
                    let values: Vec<bool> = (0..n_rows)
                        .map(|row| {
                            let offset = row * total_size + member_offset;
                            raw_buf[offset] != 0
                        })
                        .collect();
                    Some(Arc::new(BooleanArray::from(values)) as ArrayRef)
                } else {
                    log::warn!("Unsupported enum size {} for field '{}' in /{}", member_size, field_name, group_name);
                    None
                }
            }
            _ => {
                log::warn!("Unsupported HDF5 type class {:?} for field '{}' in /{}", member_class, field_name, group_name);
                None
            }
        };

        unsafe { H5Tclose(member_type_id); }

        if let Some(arr) = array {
            columns.push(field_name);
            data.push(arr);
        }
    }

    unsafe { H5Tclose(dtype_id); }

    DataFrame::new(columns, data, expected_rows)
        .map_err(|e| AnnDataError::InvalidFormat(e))
}

/// 读取普通列（非 Categorical）
fn read_column(dataset: &Dataset, expected_rows: usize) -> Result<ArrayRef> {
    let shape = dataset.shape();
    let n_rows = shape[0] as usize;
    
    if n_rows != expected_rows {
        return Err(AnnDataError::DimensionMismatch(format!(
            "Column has {} rows, expected {}",
            n_rows, expected_rows
        )));
    }
    
    // 尝试不同的数据类型
    // 1. Float64
    if let Ok(values) = dataset.read_1d::<f64>() {
        let array = Float64Array::from(values.to_vec());
        return Ok(Arc::new(array) as ArrayRef);
    }
    
    // 2. Int64
    if let Ok(values) = dataset.read_1d::<i64>() {
        let array = Int64Array::from(values.to_vec());
        return Ok(Arc::new(array) as ArrayRef);
    }
    
    // 3. Int32
    if let Ok(values) = dataset.read_1d::<i32>() {
        let array = Int32Array::from(values.to_vec());
        return Ok(Arc::new(array) as ArrayRef);
    }
    
    // 4. Boolean
    if let Ok(values) = dataset.read_1d::<bool>() {
        let array = BooleanArray::from(values.to_vec());
        return Ok(Arc::new(array) as ArrayRef);
    }
    
    // 5. String (variable length unicode)
    if let Ok(values) = dataset.read_1d::<hdf5::types::VarLenUnicode>() {
        let strings: Vec<String> = values.iter().map(|s| s.to_string()).collect();
        let array = StringArray::from(strings);
        return Ok(Arc::new(array) as ArrayRef);
    }
    
    // 5b. String (variable length ascii)
    if let Ok(values) = dataset.read_1d::<hdf5::types::VarLenAscii>() {
        let strings: Vec<String> = values.iter().map(|s| s.to_string()).collect();
        let array = StringArray::from(strings);
        return Ok(Arc::new(array) as ArrayRef);
    }
    
    // 6. String (fixed length unicode)
    if let Ok(values) = dataset.read_1d::<hdf5::types::FixedUnicode<64>>() {
        let strings: Vec<String> = values.iter().map(|s| s.to_string()).collect();
        let array = StringArray::from(strings);
        return Ok(Arc::new(array) as ArrayRef);
    }
    
    // 6b. String (fixed length ascii)
    if let Ok(values) = dataset.read_1d::<hdf5::types::FixedAscii<64>>() {
        let strings: Vec<String> = values.iter().map(|s| s.to_string()).collect();
        let array = StringArray::from(strings);
        return Ok(Arc::new(array) as ArrayRef);
    }
    
    Err(AnnDataError::UnsupportedType(format!(
        "Unsupported column data type for dataset"
    )))
}

/// 读取 Nullable 列（nullable-integer 或 nullable-boolean）
///
/// HDF5 结构（anndata encoding-type: "nullable-integer" 或 "nullable-boolean"）：
/// - /obs/{column_name}/values: 数值数组（int64, int32, float64, bool 等）
/// - /obs/{column_name}/mask: 布尔数组（true 表示 NA）
///
/// 对于 NA 位置，values 中存储 0（整数）或 false（布尔），mask 为 true。
/// 我们将 NA 位置的值替换为 null（使用 Arrow 的 null bitmap）。
fn read_nullable_column(subgroup: &Group, col_name: &str, expected_rows: usize) -> Result<ArrayRef> {
    let values_dataset = subgroup.dataset("values")?;
    let mask_dataset = subgroup.dataset("mask")?;
    
    let shape = values_dataset.shape();
    let n_rows = if shape.is_empty() { 0 } else { shape[0] as usize };
    
    if n_rows != expected_rows {
        return Err(AnnDataError::DimensionMismatch(format!(
            "Nullable column '{}' has {} rows, expected {}",
            col_name, n_rows, expected_rows
        )));
    }
    
    // 读取 mask（true = NA）
    let mask_values: Vec<bool> = mask_dataset.read_1d::<bool>()
        .map(|v| v.to_vec())
        .unwrap_or_else(|_| {
            // 如果 mask 读取失败，假设没有 NA
            vec![false; n_rows]
        });
    
    // validity bitmap: Arrow 中 true = valid, anndata mask 中 true = NA，所以取反
    let validity: Vec<bool> = mask_values.iter().map(|&m| !m).collect();
    
    // 尝试不同的数值类型
    // 1. Int64
    if let Ok(values) = values_dataset.read_1d::<i64>() {
        let vals: Vec<i64> = values.to_vec();
        let array = Int64Array::from(
            vals.into_iter()
                .zip(validity.iter())
                .map(|(v, &valid)| if valid { Some(v) } else { None })
                .collect::<Vec<Option<i64>>>()
        );
        return Ok(Arc::new(array) as ArrayRef);
    }
    
    // 2. Int32
    if let Ok(values) = values_dataset.read_1d::<i32>() {
        let vals: Vec<i32> = values.to_vec();
        let array = Int32Array::from(
            vals.into_iter()
                .zip(validity.iter())
                .map(|(v, &valid)| if valid { Some(v) } else { None })
                .collect::<Vec<Option<i32>>>()
        );
        return Ok(Arc::new(array) as ArrayRef);
    }
    
    // 3. Float64
    if let Ok(values) = values_dataset.read_1d::<f64>() {
        let vals: Vec<f64> = values.to_vec();
        let array = Float64Array::from(
            vals.into_iter()
                .zip(validity.iter())
                .map(|(v, &valid)| if valid { Some(v) } else { None })
                .collect::<Vec<Option<f64>>>()
        );
        return Ok(Arc::new(array) as ArrayRef);
    }
    
    // 4. Boolean (nullable-boolean)
    if let Ok(values) = values_dataset.read_1d::<bool>() {
        let vals: Vec<bool> = values.to_vec();
        let array = BooleanArray::from(
            vals.into_iter()
                .zip(validity.iter())
                .map(|(v, &valid)| if valid { Some(v) } else { None })
                .collect::<Vec<Option<bool>>>()
        );
        return Ok(Arc::new(array) as ArrayRef);
    }
    
    // 5. Int8
    if let Ok(values) = values_dataset.read_1d::<i8>() {
        let vals: Vec<i64> = values.to_vec().into_iter().map(|x| x as i64).collect();
        let array = Int64Array::from(
            vals.into_iter()
                .zip(validity.iter())
                .map(|(v, &valid)| if valid { Some(v) } else { None })
                .collect::<Vec<Option<i64>>>()
        );
        return Ok(Arc::new(array) as ArrayRef);
    }
    
    // 6. Int16
    if let Ok(values) = values_dataset.read_1d::<i16>() {
        let vals: Vec<i64> = values.to_vec().into_iter().map(|x| x as i64).collect();
        let array = Int64Array::from(
            vals.into_iter()
                .zip(validity.iter())
                .map(|(v, &valid)| if valid { Some(v) } else { None })
                .collect::<Vec<Option<i64>>>()
        );
        return Ok(Arc::new(array) as ArrayRef);
    }
    
    Err(AnnDataError::UnsupportedType(format!(
        "Unsupported nullable integer type for column '{}'", col_name
    )))
}

/// 读取 Categorical 列
///
/// HDF5 结构：
/// - /obs/{column_name}: 整数编码 (int32 或 int64)
/// - /obs/__categories/{column_name}: 类别标签 (string array)
///
/// 注意：-1 表示缺失值（NA），需要特殊处理
fn read_categorical_column(group: &Group, col_name: &str, expected_rows: usize) -> Result<ArrayRef> {
    // 读取整数编码
    let codes_dataset = group.dataset(col_name)?;
    let shape = codes_dataset.shape();
    let n_rows = shape[0] as usize;
    
    if n_rows != expected_rows {
        return Err(AnnDataError::DimensionMismatch(format!(
            "Categorical column '{}' has {} rows, expected {}",
            col_name, n_rows, expected_rows
        )));
    }
    
    // 读取编码（可能是 int8, int16, int32 或 int64）
    // pandas 根据类别数量自动选择最小整数类型：<128→int8, <32768→int16, 否则 int32
    let codes: Vec<i32> = if let Ok(codes_i8) = codes_dataset.read_1d::<i8>() {
        codes_i8.to_vec().into_iter().map(|x| x as i32).collect()
    } else if let Ok(codes_i16) = codes_dataset.read_1d::<i16>() {
        codes_i16.to_vec().into_iter().map(|x| x as i32).collect()
    } else if let Ok(codes_i32) = codes_dataset.read_1d::<i32>() {
        codes_i32.to_vec()
    } else if let Ok(codes_i64) = codes_dataset.read_1d::<i64>() {
        codes_i64.to_vec().into_iter().map(|x| x as i32).collect()
    } else {
        return Err(AnnDataError::UnsupportedType(
            "Categorical codes must be int8, int16, int32 or int64".to_string(),
        ));
    };
    
    // 读取类别标签
    let categories_group = group.group("__categories")?;
    let categories_dataset = categories_group.dataset(col_name)?;
    
    // 读取类别标签（支持字符串和数值类型）
    let mut categories: Vec<String> = if let Ok(cat_varlen) = categories_dataset.read_1d::<hdf5::types::VarLenUnicode>() {
        cat_varlen.iter().map(|s| s.to_string()).collect()
    } else if let Ok(cat_ascii) = categories_dataset.read_1d::<hdf5::types::VarLenAscii>() {
        cat_ascii.iter().map(|s| s.to_string()).collect()
    } else if let Ok(cat_fixed) = categories_dataset.read_1d::<hdf5::types::FixedUnicode<64>>() {
        cat_fixed.iter().map(|s| s.to_string()).collect()
    } else if let Ok(cat_fixed_ascii) = categories_dataset.read_1d::<hdf5::types::FixedAscii<64>>() {
        cat_fixed_ascii.iter().map(|s| s.to_string()).collect()
    } else if let Ok(cat_i64) = categories_dataset.read_1d::<i64>() {
        cat_i64.iter().map(|v| v.to_string()).collect()
    } else if let Ok(cat_f64) = categories_dataset.read_1d::<f64>() {
        cat_f64.iter().map(|v| v.to_string()).collect()
    } else if let Ok(cat_i32) = categories_dataset.read_1d::<i32>() {
        cat_i32.iter().map(|v| v.to_string()).collect()
    } else {
        return Err(AnnDataError::UnsupportedType(
            format!("Categorical labels for '{}' must be string type (v1)", col_name),
        ));
    };
    
    // 处理 -1（缺失值）：添加一个特殊的 NA 类别，并将 -1 映射到它
    // 检查是否有 -1 值
    let has_na = codes.iter().any(|&c| c < 0);
    
    let (adjusted_codes, final_categories) = if has_na {
        // 添加 NA 类别到末尾
        let na_index = categories.len() as i32;
        categories.push("NA".to_string());
        
        // 将 -1 映射到 NA 类别索引
        let adjusted: Vec<i32> = codes.iter().map(|&c| {
            if c < 0 { na_index } else { c }
        }).collect();
        
        (adjusted, categories)
    } else {
        (codes, categories)
    };
    
    // 创建 Arrow Dictionary Array
    let keys = Int32Array::from(adjusted_codes);
    let values = StringArray::from(final_categories);
    let dict_array = DictionaryArray::<Int32Type>::try_new(keys, Arc::new(values))
        .map_err(|e| AnnDataError::InvalidFormat(format!("Failed to create dictionary array: {}", e)))?;
    
    Ok(Arc::new(dict_array) as ArrayRef)
}

/// 读取降维嵌入 (/obsm)
///
/// HDF5 结构：
/// - /obsm: Group，包含多个嵌入
/// - 每个嵌入是一个 Dataset (n_cells × n_components)
fn read_embeddings(file: &File, expected_cells: usize) -> Result<HashMap<String, Embedding>> {
    // 先尝试作为 Group 打开（新版本格式）
    // 如果失败，尝试作为 compound Dataset 打开（旧版本格式）
    let obsm_group = match file.group("obsm") {
        Ok(g) => g,
        Err(_) => {
            // 旧版本格式：obsm 是 compound Dataset，每个字段是一个嵌入
            if let Ok(dataset) = file.dataset("obsm") {
                log::info!("Reading legacy compound Dataset for /obsm");
                return read_embeddings_legacy(&dataset, expected_cells);
            }
            return Ok(HashMap::new());
        }
    };
    let member_names = obsm_group.member_names()?;
    
    let mut embeddings = HashMap::new();
    
    for name in member_names {
        // 尝试作为 Dataset 打开
        if let Ok(dataset) = obsm_group.dataset(&name) {
            let embedding = read_embedding_dataset(&dataset, &name, expected_cells)?;
            embeddings.insert(name.clone(), embedding);
        }
    }
    
    Ok(embeddings)
}

/// 读取单个嵌入 Dataset
fn read_embedding_dataset(dataset: &Dataset, name: &str, expected_cells: usize) -> Result<Embedding> {
    let shape = dataset.shape();
    
    if shape.len() != 2 {
        return Err(AnnDataError::InvalidFormat(format!(
            "Embedding '{}' must be 2D, got {}D",
            name,
            shape.len()
        )));
    }
    
    let n_rows = shape[0] as usize;
    let n_cols = shape[1] as usize;
    
    if n_rows != expected_cells {
        return Err(AnnDataError::DimensionMismatch(format!(
            "Embedding '{}' has {} rows, expected {}",
            name, n_rows, expected_cells
        )));
    }
    
    // 读取数据（行优先）— 支持 float64 和 float32
    let data: Vec<f64> = if let Ok(arr) = dataset.read_2d::<f64>() {
        arr.into_raw_vec()
    } else if let Ok(arr) = dataset.read_2d::<f32>() {
        arr.into_raw_vec().into_iter().map(|x| x as f64).collect()
    } else {
        return Err(AnnDataError::UnsupportedType(
            format!("Embedding '{}' must be float32 or float64", name),
        ));
    };
    
    Embedding::new(name.to_string(), data, n_rows, n_cols)
        .map_err(|e| AnnDataError::InvalidFormat(e))
}

/// 读取旧版本 compound Dataset 格式的嵌入 (anndata < 0.7)
///
/// 旧版本 anndata 将 obsm 存储为 compound Dataset，
/// 每个字段是一个嵌入（可能是多维数组类型）。
fn read_embeddings_legacy(dataset: &Dataset, _expected_cells: usize) -> Result<HashMap<String, Embedding>> {
    use hdf5_sys::h5t::{
        H5Tget_nmembers, H5Tget_member_name, H5Tget_member_type, H5Tget_member_offset,
        H5Tget_class, H5Tget_size, H5Tclose, H5Tget_array_dims,
        H5T_class_t,
    };
    use hdf5_sys::h5d::H5Dread;
    use hdf5_sys::h5s::H5S_ALL;
    use hdf5_sys::h5p::H5P_DEFAULT;
    use std::ffi::CStr;

    let shape = dataset.shape();
    if shape.is_empty() {
        return Ok(HashMap::new());
    }
    let n_rows = shape[0];

    let ds_id = dataset.id();
    let dtype_id = unsafe { hdf5_sys::h5d::H5Dget_type(ds_id) };
    if dtype_id < 0 {
        return Ok(HashMap::new());
    }

    let n_fields = unsafe { H5Tget_nmembers(dtype_id) } as usize;
    let total_size = unsafe { H5Tget_size(dtype_id) };

    // 读取整个 compound dataset 为原始字节
    let buf_size = n_rows * total_size;
    let mut raw_buf: Vec<u8> = vec![0u8; buf_size];

    let read_status = unsafe {
        H5Dread(
            ds_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            raw_buf.as_mut_ptr() as *mut std::ffi::c_void,
        )
    };

    if read_status < 0 {
        unsafe { H5Tclose(dtype_id); }
        return Ok(HashMap::new());
    }

    let mut embeddings = HashMap::new();

    for i in 0..n_fields {
        let field_name_ptr = unsafe { H5Tget_member_name(dtype_id, i as u32) };
        if field_name_ptr.is_null() { continue; }
        let field_name = unsafe { CStr::from_ptr(field_name_ptr) }
            .to_string_lossy().to_string();
        unsafe { hdf5_sys::h5::H5free_memory(field_name_ptr as *mut std::ffi::c_void); }

        let member_type_id = unsafe { H5Tget_member_type(dtype_id, i as u32) };
        let member_offset = unsafe { H5Tget_member_offset(dtype_id, i as u32) };
        let member_class = unsafe { H5Tget_class(member_type_id) };

        // 嵌入字段通常是 H5T_ARRAY 类型（固定大小数组）
        if member_class == H5T_class_t::H5T_ARRAY {
            // 获取数组维度
            let ndims = unsafe { hdf5_sys::h5t::H5Tget_array_ndims(member_type_id) };
            if ndims == 1 {
                let mut dims = [0u64; 1];
                unsafe { H5Tget_array_dims(member_type_id, dims.as_mut_ptr()); }
                let n_cols = dims[0] as usize;

                // 获取基础类型
                let base_type_id = unsafe { hdf5_sys::h5t::H5Tget_super(member_type_id) };
                let base_class = unsafe { H5Tget_class(base_type_id) };
                let base_size = unsafe { H5Tget_size(base_type_id) };

                let data: Vec<f64> = if base_class == H5T_class_t::H5T_FLOAT {
                    if base_size == 8 {
                        // float64 array
                        let mut vals = Vec::with_capacity(n_rows * n_cols);
                        for row in 0..n_rows {
                            let offset = row * total_size + member_offset;
                            for col in 0..n_cols {
                                let pos = offset + col * 8;
                                let bytes: [u8; 8] = raw_buf[pos..pos + 8].try_into().unwrap();
                                vals.push(f64::from_ne_bytes(bytes));
                            }
                        }
                        vals
                    } else if base_size == 4 {
                        // float32 array → promote to f64
                        let mut vals = Vec::with_capacity(n_rows * n_cols);
                        for row in 0..n_rows {
                            let offset = row * total_size + member_offset;
                            for col in 0..n_cols {
                                let pos = offset + col * 4;
                                let bytes: [u8; 4] = raw_buf[pos..pos + 4].try_into().unwrap();
                                vals.push(f32::from_ne_bytes(bytes) as f64);
                            }
                        }
                        vals
                    } else {
                        unsafe { H5Tclose(base_type_id); H5Tclose(member_type_id); }
                        continue;
                    }
                } else {
                    unsafe { H5Tclose(base_type_id); H5Tclose(member_type_id); }
                    continue;
                };

                unsafe { H5Tclose(base_type_id); }

                if let Ok(embedding) = Embedding::new(field_name.clone(), data, n_rows, n_cols) {
                    embeddings.insert(field_name, embedding);
                }
            }
        }

        unsafe { H5Tclose(member_type_id); }
    }

    unsafe { H5Tclose(dtype_id); }
    Ok(embeddings)
}

/// 读取 layers (/layers)
///
/// HDF5 结构：
/// - /layers: Group，包含多个层
/// - 每层是一个表达矩阵（稀疏或稠密）
fn read_layers(file: &File, expected_cells: usize, expected_genes: usize) -> Result<HashMap<String, ExpressionMatrix>> {
    let layers_group = file.group("layers")?;
    let member_names = layers_group.member_names()?;
    
    let mut layers = HashMap::new();
    
    for name in member_names {
        // 尝试作为 Group 打开（稀疏矩阵）
        if let Ok(layer_group) = layers_group.group(&name) {
            let matrix = read_sparse_csr_matrix(&layer_group)?;
            let (n_rows, n_cols) = matrix.shape();
            
            if n_rows != expected_cells || n_cols != expected_genes {
                return Err(AnnDataError::DimensionMismatch(format!(
                    "Layer '{}' has shape ({}, {}), expected ({}, {})",
                    name, n_rows, n_cols, expected_cells, expected_genes
                )));
            }
            
            layers.insert(name.clone(), matrix);
        }
        // 尝试作为 Dataset 打开（稠密矩阵）
        else if let Ok(layer_dataset) = layers_group.dataset(&name) {
            let matrix = read_dense_matrix(&layer_dataset)?;
            let (n_rows, n_cols) = matrix.shape();
            
            if n_rows != expected_cells || n_cols != expected_genes {
                return Err(AnnDataError::DimensionMismatch(format!(
                    "Layer '{}' has shape ({}, {}), expected ({}, {})",
                    name, n_rows, n_cols, expected_cells, expected_genes
                )));
            }
            
            layers.insert(name.clone(), matrix);
        }
    }
    
    Ok(layers)
}

/// 读取成对矩阵 (/obsp 或 /varp)
///
/// HDF5 结构：
/// - /obsp 或 /varp: Group，包含多个成对矩阵
/// - 每个矩阵是一个方阵（稀疏或稠密）
fn read_pairwise_matrices(file: &File, group_name: &str, expected_size: usize) -> Result<HashMap<String, PairwiseMatrix>> {
    let group = file.group(group_name)?;
    let member_names = group.member_names()?;
    
    let mut pairwise_matrices = HashMap::new();
    
    for name in member_names {
        // 尝试作为 Group 打开（稀疏矩阵）
        if let Ok(pw_group) = group.group(&name) {
            let matrix = read_sparse_csr_matrix(&pw_group)?;
            let (n_rows, n_cols) = matrix.shape();
            
            if n_rows != expected_size || n_cols != expected_size {
                return Err(AnnDataError::DimensionMismatch(format!(
                    "Pairwise matrix '{}' in {} has shape ({}, {}), expected ({}, {})",
                    name, group_name, n_rows, n_cols, expected_size, expected_size
                )));
            }
            
            let pw = PairwiseMatrix::new(name.clone(), matrix)
                .map_err(|e| AnnDataError::InvalidFormat(e))?;
            pairwise_matrices.insert(name.clone(), pw);
        }
        // 尝试作为 Dataset 打开（稠密矩阵）
        else if let Ok(pw_dataset) = group.dataset(&name) {
            let matrix = read_dense_matrix(&pw_dataset)?;
            let (n_rows, n_cols) = matrix.shape();
            
            if n_rows != expected_size || n_cols != expected_size {
                return Err(AnnDataError::DimensionMismatch(format!(
                    "Pairwise matrix '{}' in {} has shape ({}, {}), expected ({}, {})",
                    name, group_name, n_rows, n_cols, expected_size, expected_size
                )));
            }
            
            let pw = PairwiseMatrix::new(name.clone(), matrix)
                .map_err(|e| AnnDataError::InvalidFormat(e))?;
            pairwise_matrices.insert(name.clone(), pw);
        }
    }
    
    Ok(pairwise_matrices)
}

/// 读取空间数据 (/obsm['spatial'] 和 /uns['spatial'])
///
/// AnnData 空间数据结构：
/// - /obsm/spatial: 空间坐标矩阵 (n_cells × 2 或 n_cells × 3)
/// - /uns/spatial: Group，包含图像和缩放因子
///   - /uns/spatial/{library_id}/images: Group，包含多分辨率图像
///   - /uns/spatial/{library_id}/scalefactors: 缩放因子
fn read_spatial_data(file: &File, n_cells: usize) -> Result<Option<crate::ir::SpatialData>> {
    use crate::ir::{SpatialData, SpatialImage};
    
    // 1. 读取空间坐标 from /obsm/spatial
    let coordinates = if file.link_exists("obsm") {
        // obsm may be a compound Dataset in old anndata (<0.7) instead of a Group
        let obsm_group = match file.group("obsm") {
            Ok(g) => g,
            Err(_) => {
                // Old format: obsm is a compound Dataset, skip spatial coordinate reading
                log::info!("obsm is not a Group (likely old anndata format), skipping spatial data");
                return Ok(None);
            }
        };
        if obsm_group.link_exists("spatial") {
            let spatial_dataset = obsm_group.dataset("spatial")?;
            let shape = spatial_dataset.shape();
            
            if shape.len() != 2 {
                return Err(AnnDataError::InvalidFormat(
                    "Spatial coordinates must be 2D array".to_string(),
                ));
            }
            
            let n_rows = shape[0] as usize;
            let n_cols = shape[1] as usize;
            
            if n_rows != n_cells {
                return Err(AnnDataError::DimensionMismatch(format!(
                    "Spatial coordinates have {} rows, expected {}",
                    n_rows, n_cells
                )));
            }
            
            if n_cols != 2 && n_cols != 3 {
                return Err(AnnDataError::InvalidFormat(format!(
                    "Spatial coordinates must have 2 or 3 columns, got {}",
                    n_cols
                )));
            }
            
            // 读取坐标数据
            let coords: Vec<f64> = spatial_dataset.read_2d::<f64>()?.into_raw_vec();
            Some((coords, n_cols))
        } else {
            None
        }
    } else {
        None
    };
    
    // 如果没有空间坐标，返回 None
    let (coordinates, n_dims) = match coordinates {
        Some((coords, dims)) => (coords, dims),
        None => return Ok(None),
    };
    
    // 2. 读取图像和缩放因子 from /uns/spatial
    let (images, scale_factors) = if file.link_exists("uns") {
        let uns_group = file.group("uns")?;
        if uns_group.link_exists("spatial") {
            let spatial_group = uns_group.group("spatial")?;
            
            // 获取第一个 library_id（通常是唯一的）
            let library_ids = spatial_group.member_names()?;
            if library_ids.is_empty() {
                (None, None)
            } else {
                let library_id = &library_ids[0];
                let library_group = spatial_group.group(library_id)?;
                
                // 读取图像
                let images = if library_group.link_exists("images") {
                    let images_group = library_group.group("images")?;
                    let image_names = images_group.member_names()?;
                    
                    let mut spatial_images = Vec::new();
                    for img_name in image_names {
                        if let Ok(img_dataset) = images_group.dataset(&img_name) {
                            // 读取图像数据（通常是 uint8 数组）
                            let shape = img_dataset.shape();
                            if shape.len() == 3 {
                                let height = shape[0] as usize;
                                let width = shape[1] as usize;
                                let _channels = shape[2] as usize;
                                
                                // 读取图像数据（HDF5 存储为 C 顺序，即行优先）
                                // 使用 read_raw 读取原始字节数据
                                let img_data: Vec<u8> = img_dataset.read_raw()?;
                                
                                // 创建 SpatialImage
                                let spatial_img = SpatialImage::new(
                                    img_name.clone(),
                                    img_data,
                                    width,
                                    height,
                                );
                                spatial_images.push(spatial_img);
                            }
                        }
                    }
                    
                    if spatial_images.is_empty() {
                        None
                    } else {
                        Some(spatial_images)
                    }
                } else {
                    None
                };
                
                // 读取缩放因子
                let scale_factors = if library_group.link_exists("scalefactors") {
                    let scalefactors_group = library_group.group("scalefactors")?;
                    let factor_names = scalefactors_group.member_names()?;
                    
                    let mut factors = std::collections::HashMap::new();
                    for factor_name in factor_names {
                        if let Ok(factor_dataset) = scalefactors_group.dataset(&factor_name) {
                            // 读取单个浮点数
                            if let Ok(value) = factor_dataset.read_scalar::<f64>() {
                                factors.insert(factor_name.clone(), value);
                            }
                        }
                    }
                    
                    if factors.is_empty() {
                        None
                    } else {
                        Some(factors)
                    }
                } else {
                    None
                };
                
                (images, scale_factors)
            }
        } else {
            (None, None)
        }
    } else {
        (None, None)
    };
    
    // 3. 构造 SpatialData
    let spatial_data = SpatialData::new(coordinates, n_dims, images, scale_factors)
        .map_err(|e| AnnDataError::InvalidFormat(e))?;
    
    Ok(Some(spatial_data))
}

#[cfg(test)]
mod tests {
    use super::*;

    // 注意：这些测试需要实际的 .h5ad 文件
    // 在实际测试中，我们会使用 Python 创建测试文件

    #[test]
    fn test_read_h5ad_missing_file() {
        let result = read_h5ad("nonexistent.h5ad");
        assert!(result.is_err());
    }

    // TODO: 添加更多测试，使用 Python 创建的测试 .h5ad 文件
}

/// 读取 Categorical 列（新版本格式 0.2.0）
///
/// HDF5 结构：
/// - /obs/{column_name}: Group
/// - /obs/{column_name}/codes: 整数编码 (int32)
/// - /obs/{column_name}/categories: 类别标签 (string array)
///
/// 注意：-1 表示缺失值（NA），需要特殊处理
fn read_categorical_column_v2(subgroup: &Group, col_name: &str, expected_rows: usize) -> Result<ArrayRef> {
    // 读取编码
    let codes_dataset = subgroup.dataset("codes")?;
    let shape = codes_dataset.shape();
    let n_rows = shape[0] as usize;
    
    if n_rows != expected_rows {
        return Err(AnnDataError::DimensionMismatch(format!(
            "Categorical column '{}' has {} rows, expected {}",
            col_name, n_rows, expected_rows
        )));
    }
    
    // 读取编码（可能是 int8, int16, int32 或 int64）
    // pandas 根据类别数量自动选择最小整数类型：<128→int8, <32768→int16, 否则 int32
    let codes: Vec<i32> = if let Ok(codes_i8) = codes_dataset.read_1d::<i8>() {
        codes_i8.to_vec().into_iter().map(|x| x as i32).collect()
    } else if let Ok(codes_i16) = codes_dataset.read_1d::<i16>() {
        codes_i16.to_vec().into_iter().map(|x| x as i32).collect()
    } else if let Ok(codes_i32) = codes_dataset.read_1d::<i32>() {
        codes_i32.to_vec()
    } else if let Ok(codes_i64) = codes_dataset.read_1d::<i64>() {
        codes_i64.to_vec().into_iter().map(|x| x as i32).collect()
    } else {
        return Err(AnnDataError::UnsupportedType(
            "Categorical codes must be int8, int16, int32 or int64".to_string(),
        ));
    };
    
    // 读取类别标签（支持字符串和数值类型）
    let categories_dataset = subgroup.dataset("categories")?;
    
    let mut categories: Vec<String> = if let Ok(cat_varlen) = categories_dataset.read_1d::<hdf5::types::VarLenUnicode>() {
        cat_varlen.iter().map(|s| s.to_string()).collect()
    } else if let Ok(cat_ascii) = categories_dataset.read_1d::<hdf5::types::VarLenAscii>() {
        cat_ascii.iter().map(|s| s.to_string()).collect()
    } else if let Ok(cat_fixed) = categories_dataset.read_1d::<hdf5::types::FixedUnicode<64>>() {
        cat_fixed.iter().map(|s| s.to_string()).collect()
    } else if let Ok(cat_fixed_ascii) = categories_dataset.read_1d::<hdf5::types::FixedAscii<64>>() {
        cat_fixed_ascii.iter().map(|s| s.to_string()).collect()
    } else if let Ok(cat_i64) = categories_dataset.read_1d::<i64>() {
        // 数值类型的 categorical（如基因长度）— 转为字符串
        cat_i64.iter().map(|v| v.to_string()).collect()
    } else if let Ok(cat_f64) = categories_dataset.read_1d::<f64>() {
        cat_f64.iter().map(|v| v.to_string()).collect()
    } else if let Ok(cat_i32) = categories_dataset.read_1d::<i32>() {
        cat_i32.iter().map(|v| v.to_string()).collect()
    } else {
        return Err(AnnDataError::UnsupportedType(
            format!("Categorical labels for '{}' must be string type (v2)", col_name),
        ));
    };
    
    // 处理 -1（缺失值）：添加一个特殊的 NA 类别，并将 -1 映射到它
    let has_na = codes.iter().any(|&c| c < 0);
    
    let (adjusted_codes, final_categories) = if has_na {
        // 添加 NA 类别到末尾
        let na_index = categories.len() as i32;
        categories.push("NA".to_string());
        
        // 将 -1 映射到 NA 类别索引
        let adjusted: Vec<i32> = codes.iter().map(|&c| {
            if c < 0 { na_index } else { c }
        }).collect();
        
        (adjusted, categories)
    } else {
        (codes, categories)
    };
    
    // 创建 Arrow Dictionary Array
    let keys = Int32Array::from(adjusted_codes);
    let values = StringArray::from(final_categories);
    let dict_array = DictionaryArray::<Int32Type>::try_new(keys, Arc::new(values))
        .map_err(|e| AnnDataError::InvalidFormat(format!("Failed to create dictionary array: {}", e)))?;
    
    Ok(Arc::new(dict_array) as ArrayRef)
}
