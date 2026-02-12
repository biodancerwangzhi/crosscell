#!/bin/bash
# 运行 .h5ad 文件读取测试

set -e

echo "=========================================="
echo "Task 7.1: .h5ad 文件读取测试"
echo "=========================================="
echo ""

# Step 1: 生成测试数据
echo "[1/3] 生成测试数据..."
python3 tests/create_test_h5ad.py
echo ""

# Step 2: 编译项目
echo "[2/3] 编译项目..."
cargo build --release
echo ""

# Step 3: 运行测试
echo "[3/3] 运行测试..."
cargo test --test test_h5ad_read -- --nocapture
echo ""

echo "=========================================="
echo "✅ 所有测试完成！"
echo "=========================================="
