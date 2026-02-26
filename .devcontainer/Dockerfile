# ============================================
# CrossCell Development Container
# ============================================
# FROM ubuntu:22.04
FROM satijalab/seurat:5.0.0


ENV DEBIAN_FRONTEND=noninteractive

# 1. 安装基本证书支持
RUN apt-get update && apt-get install -y ca-certificates wget gnupg && \
    update-ca-certificates

# 2. 替换为 CERNET 镜像（官方推荐方式）
RUN cat > /etc/apt/sources.list <<EOF
deb https://mirrors.cernet.edu.cn/ubuntu/ jammy main restricted universe multiverse
deb https://mirrors.cernet.edu.cn/ubuntu/ jammy-updates main restricted universe multiverse
deb https://mirrors.cernet.edu.cn/ubuntu/ jammy-backports main restricted universe multiverse
deb http://security.ubuntu.com/ubuntu/ jammy-security main restricted universe multiverse
EOF



# ============================================
# 安装基础系统工具 + R 环境
# ============================================
RUN apt-get update && apt-get install -y \
    build-essential \
    pkg-config \
    cmake \
    git \
    curl \
    wget \
    vim \
    ca-certificates \
    software-properties-common \
    dirmngr \
    gnupg \
    libclang-dev \
    llvm-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    hdf5-tools \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libglpk-dev \
    libopenblas-dev \
    && rm -rf /var/lib/apt/lists/*

# ============================================
# 安装 CMake >= 3.26（hdf5-metno-sys 静态编译 HDF5 2.0.0 需要）
# HDF5 完全静态链接，不需要系统 libhdf5-dev
# ============================================
RUN wget -q https://github.com/Kitware/CMake/releases/download/v3.31.4/cmake-3.31.4-linux-x86_64.tar.gz -O /tmp/cmake.tar.gz && \
    tar xzf /tmp/cmake.tar.gz -C /opt/ && \
    rm /tmp/cmake.tar.gz
ENV PATH=/opt/cmake-3.31.4-linux-x86_64/bin:$PATH

# ============================================
# 安装 Rust（使用 rustup）
# ============================================
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain stable
ENV PATH=/root/.cargo/bin:$PATH

# 安装 Rust 工具
RUN rustup component add rustfmt clippy && \
    cargo install cargo-watch

# ============================================
# 安装 Miniforge（Python + R）
# ============================================
# 复制本地的 Miniforge 安装包
# COPY miniforge.sh /tmp/miniforge.sh

# 安装 Miniforge
# RUN bash /tmp/miniforge.sh -b -p /opt/conda && \
#     rm /tmp/miniforge.sh

# ENV PATH=/opt/conda/bin:$PATH


# 彻底禁用官方 conda-forge
# =========================
# RUN conda config --remove-key channels || true && \
#     conda config --remove-key default_channels || true && \
#     conda config --set channel_priority strict && \
#     conda config --set auto_activate_base false

# =========================
# 强制使用国内 Conda 镜像
# =========================
# RUN conda config --remove-key channels || true && \
#     conda config --set channel_priority strict && \
#     conda config --add channels https://mirrors.nju.edu.cn/anaconda/cloud/conda-forge && \
#     conda config --add channels https://mirrors.nju.edu.cn/anaconda/cloud/bioconda && \
#     conda config --add channels https://mirrors.nju.edu.cn/anaconda/pkgs/main && \
#     conda config --add channels https://mirrors.nju.edu.cn/anaconda/pkgs/r && \
#     conda config --add channels https://mirrors.nju.edu.cn/anaconda/pkgs/msys2 && \
#     conda config --set show_channel_urls yes

# ============================================
# 安装 Python 环境和科学计算库（通过 conda）
# ============================================
# RUN conda install -y mamba -n base -c https://mirrors.nju.edu.cn/anaconda/cloud/conda-forge

# RUN mamba install -y \
#     python \
#      pyarrow \
#     && conda clean -afy

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-dev \
    build-essential \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# ----------------------------
# Python 科学计算
# ----------------------------
RUN pip install numpy 
RUN pip install scipy 
RUN pip install pandas 
RUN pip install anndata 
RUN pip install scanpy 

RUN pip install igraph 
RUN pip install leidenalg
RUN pip install scvelo
RUN pip install diopy
RUN pip install easySCFpy

# ----------------------------
# CrossCell Reader 开发依赖
# ----------------------------
RUN pip install pytest hypothesis
# ============================================
# 使用 pip 安装单细胞分析库和开发工具
# ============================================
# RUN pip install --no-cache-dir \
#     numpy scipy pandas anndata \
#     --index-url https://pypi.org/simple \
#     --trusted-host pypi.org \
#     --trusted-host files.pythonhosted.org


# ============================================
# 安装 R ≥ 4.3（CRAN 官方 apt 源，关键修改）
# ============================================
# RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | apt-key add -

# RUN add-apt-repository \
#   "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/"

# RUN apt-get update && apt-get install -y \
#     r-base \
#     r-base-dev \
#     && rm -rf /var/lib/apt/lists/*

# ============================================
# 配置 CRAN 国内镜像（全局，极其重要）
# ============================================
# RUN echo "options( \
#   repos = c(CRAN = 'https://mirrors.tuna.tsinghua.edu.cn/CRAN'), \
#   download.file.method = 'libcurl', \
#   Ncpus = 2 \
# )" >> /usr/lib/R/etc/Rprofile.site

# 清空旧 R 包（防止 Matrix 幽灵版本）
# RUN rm -rf /usr/local/lib/R/site-library/*
# ============================================

# ============================================
# R 包安装（CRAN，正确顺序）
# ============================================
# RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"

# RUN R -e "BiocManager::install( \
#     c('SingleCellExperiment','SummarizedExperiment','Biobase'), \
#     ask = FALSE, update = FALSE \
# )"

# RUN R -e "install.packages(c('arrow','Matrix'), repos='https://cloud.r-project.org')"
# RUN R -e "install.packages('Seurat', dependencies = TRUE, repos='https://cloud.r-project.org')"

RUN apt-get update && apt-get install -y \
    build-essential \
    pkg-config \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libwebp-dev

RUN R -e "install.packages('devtools')"
RUN R -e "devtools::install_github('xleizi/easySCF/r')"
RUN R -e "devtools::install_github('JiekaiLab/dior')"

# 设置环境变量
# ============================================
ENV RUST_BACKTRACE=1
ENV CARGO_INCREMENTAL=1
ENV PYTHONPATH=/workspace/python
# ENV R_LIBS_USER=/workspace/R/library

# ============================================
# 初始化 conda
# ============================================
# RUN /opt/conda/bin/conda init bash

WORKDIR /workspace

CMD ["/bin/bash"]
