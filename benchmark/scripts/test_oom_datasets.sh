#!/bin/bash
# CrossCell 大数据集测试脚本 — 在云服务器上运行
# 这些数据集在本地 16GB 机器上 OOM，需要更大内存
#
# 使用方法（两种模式）：
#
#   模式 A：直接调用静态二进制（推荐）
#     export CROSSCELL_BIN=~/crosscell
#     export DATA_DIR=./data/generated
#     export RESULTS_DIR=./results/cellxgene_11
#     ./test_oom_datasets.sh
#
#   模式 B：通过 Docker 运行
#     export USE_DOCKER=1
#     ./test_oom_datasets.sh
#
# 前提：
#   模式 A：~/crosscell 静态二进制已 scp 到服务器
#   模式 B：crosscell-benchmark:latest 镜像已 load

# 配置（可通过环境变量覆盖）
CROSSCELL_BIN="${CROSSCELL_BIN:-$HOME/crosscell}"
DATA_DIR="${DATA_DIR:-./data/generated}"
RESULTS_DIR="${RESULTS_DIR:-./results/cellxgene_11}"
USE_DOCKER="${USE_DOCKER:-0}"

# Docker 模式的路径映射
DOCKER_DATA_DIR="./crosscell-benchmark/data"
DOCKER_RESULTS_DIR="./crosscell-benchmark/results/cellxgene_11"

mkdir -p "$RESULTS_DIR"

# 6 个 OOM 数据集，按规模从小到大
# 注意：retina_244k 已排除（异常 dense 结构，18GB 文件）
DATASETS=(
    "cellxgene_brain_dlpfc_172k"
    "cellxgene_gut_428k"
    "cellxgene_heart_486k"
    "cellxgene_hlca_core_585k"
    "cellxgene_combat_pbmc_836k"
    "cellxgene_eqtl_autoimmune_1.2M"
)

# 运行命令的封装函数
run_crosscell() {
    local args="$@"
    if [ "$USE_DOCKER" = "1" ]; then
        docker run --rm \
            -v "$DOCKER_DATA_DIR":/benchmark/data \
            -v "$DOCKER_RESULTS_DIR":/benchmark/results/cellxgene_11 \
            crosscell-benchmark:latest \
            /usr/bin/time -v crosscell $args
    else
        /usr/bin/time -v $CROSSCELL_BIN $args
    fi
}

echo "========================================"
echo "CrossCell 大数据集测试 (OOM on 16GB)"
echo "========================================"
if [ "$USE_DOCKER" = "1" ]; then
    echo "模式: Docker (crosscell-benchmark:latest)"
else
    echo "模式: 直接二进制 ($CROSSCELL_BIN)"
    echo "数据目录: $DATA_DIR"
fi
echo "结果目录: $RESULTS_DIR"
echo ""

for ds in "${DATASETS[@]}"; do
    echo "──────────────────────────────────────"
    echo "📊 ${ds}"
    echo "──────────────────────────────────────"

    if [ "$USE_DOCKER" = "1" ]; then
        H5AD="/benchmark/data/generated/${ds}.h5ad"
        RDS="/benchmark/results/cellxgene_11/${ds}.rds"
        H5AD_RT="/benchmark/results/cellxgene_11/${ds}_roundtrip.h5ad"
    else
        H5AD="${DATA_DIR}/${ds}.h5ad"
        RDS="${RESULTS_DIR}/${ds}.rds"
        H5AD_RT="${RESULTS_DIR}/${ds}_roundtrip.h5ad"
    fi

    # 1. inspect
    echo "[inspect]"
    run_crosscell inspect -i "$H5AD" --detailed \
        2>&1 | tee "${RESULTS_DIR}/${ds}_inspect.log"
    echo ""

    # 2. H5AD → RDS
    echo "[h5ad→rds]"
    run_crosscell convert -i "$H5AD" -o "$RDS" -f seurat -v \
        2>&1 | tee "${RESULTS_DIR}/${ds}_h2r.log"
    echo ""

    # 3. RDS → H5AD (roundtrip) — 仅在上一步成功时执行
    if [ -f "${RESULTS_DIR}/${ds}.rds" ]; then
        echo "[rds→h5ad]"
        run_crosscell convert -i "$RDS" -o "$H5AD_RT" -f anndata -v \
            2>&1 | tee "${RESULTS_DIR}/${ds}_r2h.log"

        # 清理 roundtrip 文件节省磁盘
        rm -f "${RESULTS_DIR}/${ds}_roundtrip.h5ad"
        echo "  🗑️ 已清理 roundtrip 文件"
    else
        echo "[rds→h5ad] ⏭️ 跳过（H5AD→RDS 未成功）"
    fi

    echo ""
done

echo "========================================"
echo "✅ 全部完成，日志保存在 ${RESULTS_DIR}/"
echo "========================================"
