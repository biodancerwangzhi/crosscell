#!/bin/bash
# ============================================
# CrossCell RDS 读取正确性测试脚本
# ============================================
# 测试 CrossCell 读取 SeuratData 生成的 RDS 文件
#
# 使用方法:
#   docker-compose -f benchmark/docker-compose.benchmark.yml run --rm benchmark bash /benchmark/scripts/test_rds_read.sh
#
# 或者在开发容器中:
#   docker-compose run --rm dev bash scripts/test_rds_read.sh
# ============================================

set -e

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 数据目录（根据运行环境自动检测）
if [ -d "/benchmark/data/generated" ]; then
    DATA_DIR="/benchmark/data/generated"
elif [ -d "data/generated" ]; then
    DATA_DIR="data/generated"
else
    echo -e "${RED}错误: 找不到数据目录${NC}"
    exit 1
fi

echo "============================================"
echo "CrossCell RDS 读取正确性测试"
echo "============================================"
echo "数据目录: $DATA_DIR"
echo ""

# 统计变量
TOTAL=0
PASSED=0
FAILED=0
SKIPPED=0

# 测试函数
test_rds_file() {
    local file=$1
    local expected_cells=$2
    local expected_genes=$3
    local description=$4
    
    TOTAL=$((TOTAL + 1))
    
    if [ ! -f "$file" ]; then
        echo -e "${YELLOW}⏭️  跳过: $description${NC}"
        echo "   文件不存在: $file"
        SKIPPED=$((SKIPPED + 1))
        return
    fi
    
    echo -e "${BLUE}测试: $description${NC}"
    echo "   文件: $(basename $file)"
    
    # 运行 crosscell inspect
    output=$(crosscell inspect --input "$file" 2>&1) || {
        echo -e "${RED}❌ 失败: 读取错误${NC}"
        echo "   错误信息: $output"
        FAILED=$((FAILED + 1))
        return
    }
    
    # 提取细胞数和基因数
    cells=$(echo "$output" | grep -oP 'Cells: \K[0-9]+' || echo "0")
    genes=$(echo "$output" | grep -oP 'Genes: \K[0-9]+' || echo "0")
    sparsity=$(echo "$output" | grep -oP 'Sparsity: \K[0-9.]+' || echo "0")
    
    # 验证
    if [ "$cells" -eq "$expected_cells" ] && [ "$genes" -eq "$expected_genes" ]; then
        echo -e "${GREEN}✅ 通过${NC}"
        echo "   细胞: $cells | 基因: $genes | 稀疏度: ${sparsity}%"
        PASSED=$((PASSED + 1))
    else
        echo -e "${RED}❌ 失败: 维度不匹配${NC}"
        echo "   期望: $expected_cells cells × $expected_genes genes"
        echo "   实际: $cells cells × $genes genes"
        FAILED=$((FAILED + 1))
    fi
    echo ""
}

echo "============================================"
echo "Seurat V4 数据集测试"
echo "============================================"
echo ""

# Seurat V4 数据集（基于 seurat_v4_summary.csv）
test_rds_file "$DATA_DIR/seurat_v4_pbmc3k_raw.rds" 2700 13714 "pbmc3k (小规模)"
test_rds_file "$DATA_DIR/seurat_v4_cbmc_raw.rds" 8617 20501 "cbmc (CITE-seq)"
test_rds_file "$DATA_DIR/seurat_v4_panc8_raw.rds" 14890 34363 "panc8 (中规模)"
test_rds_file "$DATA_DIR/seurat_v4_ifnb_raw.rds" 13999 14053 "ifnb (刺激响应)"
test_rds_file "$DATA_DIR/seurat_v4_pbmcsca_raw.rds" 31021 33694 "pbmcsca (大规模)"
test_rds_file "$DATA_DIR/seurat_v4_hcabm40k_raw.rds" 40000 17369 "hcabm40k (40k cells)"
test_rds_file "$DATA_DIR/seurat_v4_ssHippo_raw.rds" 53173 23264 "ssHippo (空间数据)"
test_rds_file "$DATA_DIR/seurat_v4_bmcite_raw.rds" 30672 17009 "bmcite (多模态, 含R闭包)"
test_rds_file "$DATA_DIR/seurat_v4_thp1.eccite_raw.rds" 20729 18649 "thp1.eccite (ECCITE-seq)"
test_rds_file "$DATA_DIR/seurat_v4_celegans.embryo_raw.rds" 6188 20222 "celegans.embryo"
test_rds_file "$DATA_DIR/seurat_v4_stxKidney_raw.rds" 1438 31053 "stxKidney (空间数据)"

echo "============================================"
echo "Seurat V5 数据集测试"
echo "============================================"
echo ""

# Seurat V5 数据集（基于 seurat_v5_summary.csv）
# 注意: V5 使用 Assay5 格式，当前 CrossCell 尚未完全支持
test_rds_file "$DATA_DIR/seurat_v5_pbmc3k_raw.rds" 2700 13714 "pbmc3k V5 (Assay5)"
test_rds_file "$DATA_DIR/seurat_v5_cbmc_raw.rds" 8617 20501 "cbmc V5 (Assay5)"
test_rds_file "$DATA_DIR/seurat_v5_panc8_raw.rds" 14890 34363 "panc8 V5 (Assay5)"
test_rds_file "$DATA_DIR/seurat_v5_ifnb_raw.rds" 13999 14053 "ifnb V5 (Assay5)"
test_rds_file "$DATA_DIR/seurat_v5_pbmcsca_raw.rds" 31021 33694 "pbmcsca V5 (Assay5)"
test_rds_file "$DATA_DIR/seurat_v5_hcabm40k_raw.rds" 40000 17369 "hcabm40k V5 (Assay5)"
test_rds_file "$DATA_DIR/seurat_v5_ssHippo_raw.rds" 53173 23264 "ssHippo V5 (Assay5)"
test_rds_file "$DATA_DIR/seurat_v5_bmcite_raw.rds" 30672 17009 "bmcite V5 (Assay5)"
test_rds_file "$DATA_DIR/seurat_v5_thp1.eccite_raw.rds" 20729 18649 "thp1.eccite V5 (Assay5)"
test_rds_file "$DATA_DIR/seurat_v5_celegans.embryo_raw.rds" 6188 20222 "celegans.embryo V5 (Assay5)"

echo "============================================"
echo "测试汇总"
echo "============================================"
echo ""
echo "总计: $TOTAL"
echo -e "${GREEN}通过: $PASSED${NC}"
echo -e "${RED}失败: $FAILED${NC}"
echo -e "${YELLOW}跳过: $SKIPPED${NC}"
echo ""

if [ $FAILED -eq 0 ]; then
    echo -e "${GREEN}✅ 所有测试通过！${NC}"
    exit 0
else
    echo -e "${RED}❌ 有 $FAILED 个测试失败${NC}"
    exit 1
fi
