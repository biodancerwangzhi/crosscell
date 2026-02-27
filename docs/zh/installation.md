# 安装

CrossCell 提供三种接口：CLI、Python、R。

## CLI

### 预编译二进制（推荐）

从 [GitHub Releases](https://github.com/biodancerwangzhi/crosscell/releases) 下载：

```bash
# Linux x86_64
curl -L https://github.com/biodancerwangzhi/crosscell/releases/latest/download/crosscell-linux-x64.tar.gz | tar xz
sudo mv crosscell-linux-x64/crosscell /usr/local/bin/

# macOS Intel
curl -L https://github.com/biodancerwangzhi/crosscell/releases/latest/download/crosscell-macos-x64.tar.gz | tar xz
sudo mv crosscell-macos-x64/crosscell /usr/local/bin/

# macOS Apple Silicon
curl -L https://github.com/biodancerwangzhi/crosscell/releases/latest/download/crosscell-macos-arm64.tar.gz | tar xz
sudo mv crosscell-macos-arm64/crosscell /usr/local/bin/
```

Windows：从 Release 页面下载 `crosscell-windows-x64.zip`，解压后添加到 PATH。

### 从源码编译

需要 [Rust 工具链](https://rustup.rs) 和 CMake ≥ 3.26。

```bash
git clone https://github.com/biodancerwangzhi/crosscell.git
cd crosscell
cargo build --release
```

### 验证

```bash
crosscell --version
crosscell --help
```

## Python

### 从 PyPI

```bash
pip install crosscell
```

### 从源码

```bash
pip install maturin
git clone https://github.com/biodancerwangzhi/crosscell.git
cd crosscell/crosscell-py
maturin develop --release
```

### 验证

```python
import crosscell
print(crosscell.__version__)
```

## R

### 从 R-universe

```r
install.packages("crosscell", repos = "https://biodancerwangzhi.r-universe.dev")
```

### 从源码

需要 Rust 工具链和 CMake ≥ 3.26。

```r
install.packages("path/to/crosscell-r", repos = NULL, type = "source")
```

### 验证

```r
library(crosscell)
info <- inspect_file("test.h5ad")
```

## 依赖

| 接口 | 要求 |
|------|------|
| CLI | 无（静态链接） |
| Python | Python 3.9+, anndata, scipy, numpy, pandas |
| R | R 4.0+, Matrix；可选：Seurat, SingleCellExperiment |
