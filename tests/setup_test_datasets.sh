#!/bin/bash
# 自动下载和设置真实测试数据集
# Task 20.1 - 准备真实测试数据集

set -e  # 遇到错误立即退出

echo "======================================================================"
echo "CrossCell 真实测试数据集设置脚本"
echo "Task 20.1 - 准备真实测试数据集"
echo "======================================================================"
echo ""

# 检查是否在 Docker 容器中
if [ -f /.dockerenv ]; then
    echo "✓ 在 Docker 容器中运行"
else
    echo "⚠️  建议在 Docker 容器中运行此脚本"
    echo ""
    echo "运行命令:"
    echo "  docker run -it --rm -v \${PWD}:/workspace crosscell-dev bash"
    echo "  cd /workspace"
    echo "  bash tests/setup_test_datasets.sh"
    echo ""
    read -p "是否继续在本地环境运行? (y/N) " -n 1 -r
    echo ""
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# 创建数据目录
echo "创建数据目录..."
mkdir -p tests/data/real_datasets
echo "✓ 数据目录: tests/data/real_datasets"
echo ""

# 下载 Python/AnnData 格式数据
echo "======================================================================"
echo "步骤 1/2: 下载 AnnData 格式数据 (Python)"
echo "======================================================================"
echo ""

if command -v python3 &> /dev/null; then
    echo "✓ Python3 已安装"
    
    # 检查 Python 依赖
    if python3 -c "import scanpy, anndata" 2>/dev/null; then
        echo "✓ Python 依赖已安装 (scanpy, anndata)"
    else
        echo "⚠️  缺少 Python 依赖"
        echo "安装命令: pip install scanpy anndata"
        echo ""
        read -p "是否现在安装? (y/N) " -n 1 -r
        echo ""
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            pip install scanpy anndata
        else
            echo "跳过 Python 数据下载"
        fi
    fi
    
    # 运行 Python 下载脚本
    if python3 -c "import scanpy, anndata" 2>/dev/null; then
        echo ""
        echo "运行 Python 下载脚本..."
        python3 tests/download_test_datasets.py
        echo ""
    fi
else
    echo "❌ Python3 未安装，跳过 AnnData 格式数据下载"
    echo ""
fi

# 下载 R/Seurat 格式数据
echo "======================================================================"
echo "步骤 2/2: 下载 Seurat 格式数据 (R)"
echo "======================================================================"
echo ""

if command -v Rscript &> /dev/null; then
    echo "✓ R 已安装"
    
    # 检查 R 依赖
    if Rscript -e "library(Seurat)" 2>/dev/null; then
        echo "✓ R 依赖已安装 (Seurat)"
    else
        echo "⚠️  缺少 R 依赖"
        echo "安装命令:"
        echo "  R -e \"install.packages('Seurat')\""
        echo "  R -e \"install.packages('remotes')\""
        echo "  R -e \"remotes::install_github('satijalab/seurat-data')\""
        echo ""
        echo "跳过 R 数据下载"
    fi
    
    # 运行 R 下载脚本
    if Rscript -e "library(Seurat)" 2>/dev/null; then
        echo ""
        echo "运行 R 下载脚本..."
        Rscript tests/download_test_datasets.R
        echo ""
    fi
else
    echo "❌ R 未安装，跳过 Seurat 格式数据下载"
    echo ""
fi

# 验证下载结果
echo "======================================================================"
echo "验证下载结果"
echo "======================================================================"
echo ""

downloaded_files=0
total_size=0

for file in tests/data/real_datasets/*.h5ad tests/data/real_datasets/*.rds; do
    if [ -f "$file" ]; then
        size=$(du -h "$file" | cut -f1)
        echo "✓ $(basename $file) - $size"
        downloaded_files=$((downloaded_files + 1))
        
        # 累加文件大小（转换为字节）
        size_bytes=$(stat -f%z "$file" 2>/dev/null || stat -c%s "$file" 2>/dev/null || echo 0)
        total_size=$((total_size + size_bytes))
    fi
done

if [ $downloaded_files -eq 0 ]; then
    echo "❌ 没有下载任何数据集"
    echo ""
    echo "请检查:"
    echo "  1. Python 和 R 是否正确安装"
    echo "  2. 必需的包是否已安装"
    echo "  3. 网络连接是否正常"
    exit 1
fi

echo ""
echo "======================================================================"
echo "✅ 数据集设置完成！"
echo "======================================================================"
echo ""
echo "下载统计:"
echo "  文件数: $downloaded_files"
echo "  总大小: $(echo "scale=2; $total_size / 1024 / 1024" | bc) MB"
echo ""
echo "数据位置: $(pwd)/tests/data/real_datasets/"
echo ""
echo "下一步:"
echo "  1. 查看数据集信息: cat tests/data/real_datasets/README.md"
echo "  2. 运行准确性测试: cargo test --test test_real_datasets"
echo "  3. 手动验证数据: python3 或 R"
echo ""
