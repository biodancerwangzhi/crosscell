# 安装

## Docker（推荐）

```bash
# 克隆仓库
git clone https://github.com/yourusername/crosscell.git
cd crosscell

# 构建 Docker 镜像
docker-compose build dev

# 编译 release 版本（静态链接 HDF5）
docker-compose run --rm dev cargo build --release

# 二进制文件生成在 target/release/crosscell
# 在 Linux/WSL 中运行
chmod +x target/release/crosscell
./target/release/crosscell --version
```

## 从源码编译

### 环境要求

- Rust 1.75+
- HDF5 库
- C 编译器 (gcc/clang)

### 编译

```bash
# 克隆仓库
git clone https://github.com/yourusername/crosscell.git
cd crosscell

# 编译 release 版本
cargo build --release

# 二进制文件位置
./target/release/crosscell --version
```

### 添加到 PATH

```bash
# Linux/macOS
sudo cp target/release/crosscell /usr/local/bin/

# 或添加到 PATH
export PATH="$PATH:$(pwd)/target/release"
```

## 验证安装

```bash
# 检查版本
crosscell --version

# 显示帮助
crosscell --help

# 测试示例文件
crosscell inspect -i sample.h5ad
```

## Docker 命令

本文档中的所有命令都可以通过 Docker 运行：

```bash
# 格式
docker-compose run --rm dev crosscell <命令>

# 示例
docker-compose run --rm dev crosscell convert -i input.h5ad -o output.rds -f seurat
docker-compose run --rm dev crosscell inspect -i data.h5ad
```
