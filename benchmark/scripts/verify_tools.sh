#!/bin/bash
# ============================================
# 验证所有工具是否正确安装
# ============================================

echo "============================================"
echo "CrossCell Benchmark - Tool Verification"
echo "============================================"
echo ""

# 颜色定义
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

check_tool() {
    local name=$1
    local cmd=$2
    echo -n "Checking $name... "
    if eval "$cmd" > /dev/null 2>&1; then
        echo -e "${GREEN}✓${NC}"
        return 0
    else
        echo -e "${RED}✗${NC}"
        return 1
    fi
}

echo "--- CLI Tools ---"
check_tool "CrossCell" "crosscell --version"

echo ""
echo "--- R Packages ---"
check_tool "Zellkonverter" "Rscript -e \"library(zellkonverter)\""
check_tool "anndataR" "Rscript -e \"library(anndataR)\""
check_tool "convert2anndata" "Rscript -e \"library(convert2anndata)\""
check_tool "easySCF (R)" "Rscript -e \"library(easySCFr)\""
check_tool "Seurat" "Rscript -e \"library(Seurat)\""

echo ""
echo "--- Python Packages ---"
check_tool "easySCFpy" "python3 -c \"import easySCFpy\""
check_tool "anndata" "python3 -c \"import anndata\""
check_tool "scanpy" "python3 -c \"import scanpy\""

echo ""
echo "--- Test Data ---"
if [ -f "/benchmark/data/pbmc3k_raw.h5ad" ]; then
    echo -e "pbmc3k_raw.h5ad: ${GREEN}✓${NC}"
else
    echo -e "pbmc3k_raw.h5ad: ${RED}✗ (not found)${NC}"
fi

echo ""
echo "============================================"
echo "Verification complete!"
echo "============================================"
