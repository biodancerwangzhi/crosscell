//! 内存映射文件支持
//!
//! 提供零拷贝的文件读取，适用于大文件的随机访问。

use std::fs::File;
use std::path::Path;

use memmap2::{Mmap, MmapOptions};

use crate::error::{CrossCellError, Result};

/// 内存映射读取器
///
/// 使用 mmap 实现零拷贝文件读取，由操作系统管理页面缓存。
#[derive(Debug)]
pub struct MmapReader {
    /// 内存映射
    mmap: Mmap,
    /// 文件大小
    size: usize,
}

impl MmapReader {
    /// 打开文件并创建内存映射
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path.as_ref())?;

        let metadata = file.metadata()?;
        let size = metadata.len() as usize;

        // 创建只读内存映射
        let mmap = unsafe {
            MmapOptions::new()
                .map(&file)
                .map_err(|e| CrossCellError::InvalidFormat(format!("Failed to mmap file: {}", e)))?
        };

        Ok(Self { mmap, size })
    }

    /// 获取文件大小
    pub fn size(&self) -> usize {
        self.size
    }

    /// 获取整个文件的字节切片
    pub fn as_slice(&self) -> &[u8] {
        &self.mmap
    }

    /// 读取指定范围的字节
    pub fn read_range(&self, start: usize, end: usize) -> Result<&[u8]> {
        if end > self.size {
            return Err(CrossCellError::InvalidFormat(format!(
                "Read range {}..{} exceeds file size {}",
                start, end, self.size
            )));
        }
        if start > end {
            return Err(CrossCellError::InvalidFormat(format!(
                "Invalid range: start {} > end {}",
                start, end
            )));
        }
        Ok(&self.mmap[start..end])
    }

    /// 读取 u64（小端序）
    pub fn read_u64(&self, offset: usize) -> Result<u64> {
        let bytes = self.read_range(offset, offset + 8)?;
        Ok(u64::from_le_bytes(bytes.try_into().unwrap()))
    }

    /// 读取 f64（小端序）
    pub fn read_f64(&self, offset: usize) -> Result<f64> {
        let bytes = self.read_range(offset, offset + 8)?;
        Ok(f64::from_le_bytes(bytes.try_into().unwrap()))
    }

    /// 读取 f64 数组
    pub fn read_f64_array(&self, offset: usize, count: usize) -> Result<Vec<f64>> {
        let bytes = self.read_range(offset, offset + count * 8)?;
        let mut result = Vec::with_capacity(count);
        for i in 0..count {
            let start = i * 8;
            let val = f64::from_le_bytes(bytes[start..start + 8].try_into().unwrap());
            result.push(val);
        }
        Ok(result)
    }

    /// 读取 usize 数组
    pub fn read_usize_array(&self, offset: usize, count: usize) -> Result<Vec<usize>> {
        let bytes = self.read_range(offset, offset + count * 8)?;
        let mut result = Vec::with_capacity(count);
        for i in 0..count {
            let start = i * 8;
            let val = u64::from_le_bytes(bytes[start..start + 8].try_into().unwrap()) as usize;
            result.push(val);
        }
        Ok(result)
    }
}

/// 预取提示（告诉操作系统即将访问的区域）
#[cfg(unix)]
pub fn prefetch_range(_mmap: &Mmap, _offset: usize, _length: usize) {
    // 在 Unix 系统上可以使用 madvise 提示
    // 这里简化处理，实际可以调用 libc::madvise
}

#[cfg(not(unix))]
pub fn prefetch_range(_mmap: &Mmap, _offset: usize, _length: usize) {
    // Windows 上暂不实现预取
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_mmap_reader() {
        // 创建测试文件
        let mut file = NamedTempFile::new().unwrap();
        let data: Vec<u8> = (0..=255u8).collect();
        file.write_all(&data).unwrap();
        file.flush().unwrap();

        // 打开并读取
        let reader = MmapReader::open(file.path()).unwrap();
        assert_eq!(reader.size(), 256);

        // 读取范围
        let slice = reader.read_range(10, 20).unwrap();
        assert_eq!(slice, &data[10..20]);
    }

    #[test]
    fn test_mmap_read_numbers() {
        let mut file = NamedTempFile::new().unwrap();

        // 写入测试数据
        let val: u64 = 12345678;
        file.write_all(&val.to_le_bytes()).unwrap();
        let fval: f64 = 3.14159;
        file.write_all(&fval.to_le_bytes()).unwrap();
        file.flush().unwrap();

        let reader = MmapReader::open(file.path()).unwrap();
        assert_eq!(reader.read_u64(0).unwrap(), 12345678);
        assert!((reader.read_f64(8).unwrap() - 3.14159).abs() < 1e-10);
    }
}
