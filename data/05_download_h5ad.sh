#!/bin/bash
# ============================================
# 下载 H5AD 基准测试数据集
# ============================================
# 在 benchmark 容器中运行:
#   docker-compose -f benchmark/docker-compose.benchmark.yml run --rm benchmark bash /benchmark/data/download_h5ad.sh
#
# 分步运行:
#   docker-compose -f benchmark/docker-compose.benchmark.yml run --rm benchmark bash /benchmark/data/download_h5ad.sh scanpy
#   docker-compose -f benchmark/docker-compose.benchmark.yml run --rm benchmark bash /benchmark/data/download_h5ad.sh scvelo
#   docker-compose -f benchmark/docker-compose.benchmark.yml run --rm benchmark bash /benchmark/data/download_h5ad.sh squidpy
#   docker-compose -f benchmark/docker-compose.benchmark.yml run --rm benchmark bash /benchmark/data/download_h5ad.sh cellxgene-large
# ============================================

set -e

OUTPUT_DIR="/benchmark/data/generated"
SOURCE="${1:-all}"

echo "============================================"
echo "CrossCell Benchmark - H5AD Dataset Download"
echo "============================================"
echo "Output: $OUTPUT_DIR"
echo "Source: $SOURCE"
echo ""

mkdir -p "$OUTPUT_DIR"

python3 /benchmark/data/05_download_h5ad_datasets.py \
    --output-dir "$OUTPUT_DIR" \
    --source "$SOURCE"

echo ""
echo "============================================"
echo "下载完成！文件列表:"
echo "============================================"
ls -lh "$OUTPUT_DIR"/*.h5ad 2>/dev/null || echo "（无 H5AD 文件）"
