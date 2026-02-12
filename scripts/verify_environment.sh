#!/bin/bash
# ============================================
# CrossCell 开发环境验证脚本
# ============================================
# 用途：验证 Docker 容器中的开发环境是否正确配置
# 运行：bash scripts/verify_environment.sh

# 注释掉 set -e，让脚本在遇到错误时继续运行
# set -e  # 遇到错误立即退出

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 计数器
TOTAL=0
PASSED=0
FAILED=0

# 打印函数
print_header() {
    echo -e "\n${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n"
}

print_test() {
    echo -e "${YELLOW}▶ 测试: $1${NC}"
}

print_pass() {
    echo -e "${GREEN}  ✓ $1${NC}"
    ((PASSED++))
    ((TOTAL++))
}

print_fail() {
    echo -e "${RED}  ✗ $1${NC}"
    ((FAILED++))
    ((TOTAL++))
}

print_info() {
    echo -e "  ℹ $1"
}

# ============================================
# 1. 系统基础工具
# ============================================
print_header "1️⃣  系统基础工具"

print_test "检查 bash"
if command -v bash &> /dev/null; then
    VERSION=$(bash --version | head -n 1)
    print_pass "bash: $VERSION"
else
    print_fail "bash 未安装"
fi

print_test "检查 git"
if command -v git &> /dev/null; then
    VERSION=$(git --version)
    print_pass "$VERSION"
else
    print_fail "git 未安装"
fi

print_test "检查 curl"
if command -v curl &> /dev/null; then
    VERSION=$(curl --version | head -n 1)
    print_pass "$VERSION"
else
    print_fail "curl 未安装"
fi

print_test "检查 pkg-config"
if command -v pkg-config &> /dev/null; then
    VERSION=$(pkg-config --version)
    print_pass "pkg-config: $VERSION"
else
    print_fail "pkg-config 未安装"
fi

# ============================================
# 2. Rust 工具链
# ============================================
print_header "2️⃣  Rust 工具链"

print_test "检查 rustc"
if command -v rustc &> /dev/null; then
    VERSION=$(rustc --version)
    print_pass "$VERSION"
else
    print_fail "rustc 未安装"
fi

print_test "检查 cargo"
if command -v cargo &> /dev/null; then
    VERSION=$(cargo --version)
    print_pass "$VERSION"
else
    print_fail "cargo 未安装"
fi

print_test "检查 rustfmt"
if command -v rustfmt &> /dev/null; then
    VERSION=$(rustfmt --version)
    print_pass "$VERSION"
else
    print_fail "rustfmt 未安装"
fi

print_test "检查 clippy"
if cargo clippy --version &> /dev/null; then
    VERSION=$(cargo clippy --version)
    print_pass "$VERSION"
else
    print_fail "clippy 未安装"
fi

print_test "检查 maturin"
if command -v maturin &> /dev/null; then
    VERSION=$(maturin --version)
    print_pass "$VERSION"
else
    print_fail "maturin 未安装（如果不需要 Python 绑定可忽略）"
fi

# ============================================
# 3. Python 环境
# ============================================
print_header "3️⃣  Python 环境"

print_test "检查 Python"
if command -v python3 &> /dev/null; then
    VERSION=$(python3 --version)
    print_pass "$VERSION"
    
    # 检查 Python 路径
    PYTHON_PATH=$(which python3)
    print_info "Python 路径: $PYTHON_PATH"
else
    print_fail "Python 未安装"
fi

print_test "检查 Python 包: numpy"
if python3 -c "import numpy; print(numpy.__version__)" &> /dev/null; then
    VERSION=$(python3 -c "import numpy; print(numpy.__version__)")
    print_pass "numpy: $VERSION"
else
    print_fail "numpy 未安装"
fi

print_test "检查 Python 包: scipy"
if python3 -c "import scipy; print(scipy.__version__)" &> /dev/null; then
    VERSION=$(python3 -c "import scipy; print(scipy.__version__)")
    print_pass "scipy: $VERSION"
else
    print_fail "scipy 未安装"
fi

print_test "检查 Python 包: pandas"
if python3 -c "import pandas; print(pandas.__version__)" &> /dev/null; then
    VERSION=$(python3 -c "import pandas; print(pandas.__version__)")
    print_pass "pandas: $VERSION"
else
    print_fail "pandas 未安装"
fi

print_test "检查 Python 包: anndata"
if python3 -c "import anndata; print(anndata.__version__)" &> /dev/null; then
    VERSION=$(python3 -c "import anndata; print(anndata.__version__)")
    print_pass "anndata: $VERSION"
else
    print_fail "anndata 未安装（如果不需要 Python 互操作可忽略）"
fi

# ============================================
# 4. R 环境
# ============================================
print_header "4️⃣  R 环境"

print_test "检查 R"
if command -v R &> /dev/null; then
    VERSION=$(R --version | head -n 1)
    print_pass "$VERSION"
    
    # 检查 R 路径
    R_PATH=$(which R)
    print_info "R 路径: $R_PATH"
else
    print_fail "R 未安装"
fi

print_test "检查 Rscript"
if command -v Rscript &> /dev/null; then
    VERSION=$(Rscript --version 2>&1 | head -n 1)
    print_pass "$VERSION"
else
    print_fail "Rscript 未安装"
fi

print_test "检查 R 包: Matrix"
if Rscript -e "library(Matrix)" &> /dev/null; then
    VERSION=$(Rscript -e "cat(as.character(packageVersion('Matrix')))" 2>/dev/null)
    print_pass "Matrix: $VERSION"
else
    print_fail "Matrix 未安装"
fi

print_test "检查 R 包: Seurat"
if Rscript -e "library(Seurat)" &> /dev/null; then
    VERSION=$(Rscript -e "cat(as.character(packageVersion('Seurat')))" 2>/dev/null)
    print_pass "Seurat: $VERSION"
else
    print_fail "Seurat 未安装（需要安装以生成测试数据）"
fi

# ============================================
# 5. 系统库
# ============================================
print_header "5️⃣  系统库"

print_test "检查 HDF5"
if pkg-config --exists hdf5 2>/dev/null; then
    VERSION=$(pkg-config --modversion hdf5)
    print_pass "HDF5: $VERSION"
else
    print_fail "HDF5 未安装"
fi

print_test "检查 OpenBLAS"
if pkg-config --exists openblas 2>/dev/null; then
    VERSION=$(pkg-config --modversion openblas)
    print_pass "OpenBLAS: $VERSION"
else
    print_fail "OpenBLAS 未安装（可能使用其他 BLAS 实现）"
fi

print_test "检查 libcurl"
if pkg-config --exists libcurl 2>/dev/null; then
    VERSION=$(pkg-config --modversion libcurl)
    print_pass "libcurl: $VERSION"
else
    print_fail "libcurl 未安装"
fi

# ============================================
# 6. 项目编译测试
# ============================================
print_header "6️⃣  项目编译测试"

print_test "检查 Cargo.toml"
if [ -f "Cargo.toml" ]; then
    print_pass "Cargo.toml 存在"
    
    print_test "尝试编译项目（cargo check）"
    if cargo check 2>&1 | grep -q "Finished"; then
        print_pass "项目编译检查通过"
    else
        print_fail "项目编译检查失败（可能有依赖问题）"
    fi
else
    print_fail "Cargo.toml 不存在（请在项目根目录运行）"
fi

# ============================================
# 7. 网络连接测试
# ============================================
print_header "7️⃣  网络连接测试"

print_test "测试 crates.io 连接"
if curl -s --connect-timeout 5 https://crates.io &> /dev/null; then
    print_pass "crates.io 可访问"
else
    print_fail "crates.io 不可访问（可能影响 cargo 下载）"
fi

print_test "测试 CRAN 镜像连接"
if curl -s --connect-timeout 5 https://mirrors.nju.edu.cn/CRAN/ &> /dev/null; then
    print_pass "CRAN 镜像可访问"
else
    print_fail "CRAN 镜像不可访问（可能影响 R 包安装）"
fi

print_test "测试 conda 镜像连接"
if curl -s --connect-timeout 5 https://mirror.nju.edu.cn/anaconda/ &> /dev/null; then
    print_pass "conda 镜像可访问"
else
    print_fail "conda 镜像不可访问（可能影响 conda 包安装）"
fi

# ============================================
# 8. 总结
# ============================================
print_header "📊 测试总结"

echo -e "总测试数: ${BLUE}$TOTAL${NC}"
echo -e "通过: ${GREEN}$PASSED${NC}"
echo -e "失败: ${RED}$FAILED${NC}"

if [ $FAILED -eq 0 ]; then
    echo -e "\n${GREEN}✅ 所有测试通过！开发环境配置正确。${NC}\n"
    exit 0
else
    echo -e "\n${YELLOW}⚠️  有 $FAILED 个测试失败，请检查上述错误信息。${NC}\n"
    exit 1
fi
