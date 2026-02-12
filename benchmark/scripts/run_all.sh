#!/bin/bash
# ============================================
# CrossCell Benchmark — 一键运行所有测试
# ============================================
#
# 用法 (在 benchmark 容器内):
#   bash /benchmark/scripts/run_all.sh
#   bash /benchmark/scripts/run_all.sh --scale small
#   bash /benchmark/scripts/run_all.sh --tool crosscell
#   bash /benchmark/scripts/run_all.sh --tool crosscell --direction rds_to_h5ad
#
# 从宿主机运行:
#   docker-compose -f benchmark/docker-compose.benchmark.yml run --rm benchmark \
#     bash /benchmark/scripts/run_all.sh
# ============================================

set -e

# 颜色
GREEN='\033[0;32m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

# 默认参数
TOOL="all"
DIRECTION="all"
SCALE=""
RUNS=3
TIMEOUT=600

# 解析参数
while [[ $# -gt 0 ]]; do
    case $1 in
        --tool) TOOL="$2"; shift 2 ;;
        --direction) DIRECTION="$2"; shift 2 ;;
        --scale) SCALE="$2"; shift 2 ;;
        --runs) RUNS="$2"; shift 2 ;;
        --timeout) TIMEOUT="$2"; shift 2 ;;
        -h|--help)
            echo "用法: bash run_all.sh [选项]"
            echo ""
            echo "选项:"
            echo "  --tool NAME       测试指定工具 (crosscell/zellkonverter/anndataR/convert2anndata/easySCF/all)"
            echo "  --direction DIR   转换方向 (rds_to_h5ad/h5ad_to_rds/all)"
            echo "  --scale SIZE      只测试指定规模 (small/medium/large)"
            echo "  --runs N          每个测试运行次数 (默认: 3)"
            echo "  --timeout SEC     超时秒数 (默认: 600)"
            exit 0
            ;;
        *) echo "未知参数: $1"; exit 1 ;;
    esac
done

echo "============================================"
echo "  CrossCell Performance Benchmark"
echo "============================================"
echo ""
echo "  工具:     $TOOL"
echo "  方向:     $DIRECTION"
echo "  规模:     ${SCALE:-all}"
echo "  运行次数: $RUNS"
echo "  超时:     ${TIMEOUT}s"
echo ""

# 构建 python 命令参数
BENCH_ARGS="--tool $TOOL --direction $DIRECTION --runs $RUNS --timeout $TIMEOUT"
if [ -n "$SCALE" ]; then
    BENCH_ARGS="$BENCH_ARGS --scale $SCALE"
fi

# Step 1: 运行基准测试
echo -e "${BLUE}[1/2] 运行基准测试...${NC}"
echo ""
python3 /benchmark/scripts/run_benchmark.py $BENCH_ARGS

# Step 2: 验证输出
echo ""
echo -e "${BLUE}[2/2] 验证转换输出...${NC}"
echo ""
python3 /benchmark/scripts/validate_output.py --tool "$TOOL"

# 完成
echo ""
echo "============================================"
echo -e "${GREEN}✅ 全部完成!${NC}"
echo "============================================"
echo ""
echo "结果文件:"
ls -la /benchmark/results/*.json 2>/dev/null || echo "  (无结果文件)"
