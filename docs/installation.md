# Installation

CrossCell provides three interfaces: CLI, Python, and R.

## CLI

### Prebuilt Binaries (Recommended)

Download from [GitHub Releases](https://github.com/biodancerwangzhi/crosscell/releases):

```bash
# Linux x86_64
curl -L https://github.com/biodancerwangzhi/crosscell/releases/latest/download/crosscell-linux-x64.tar.gz | tar xz
sudo mv crosscell-linux-x64/crosscell /usr/local/bin/

# macOS Intel
curl -L https://github.com/biodancerwangzhi/crosscell/releases/latest/download/crosscell-macos-x64.tar.gz | tar xz
sudo mv crosscell-macos-x64/crosscell /usr/local/bin/

# macOS Apple Silicon
curl -L https://github.com/biodancerwangzhi/crosscell/releases/latest/download/crosscell-macos-arm64.tar.gz | tar xz
sudo mv crosscell-macos-arm64/crosscell /usr/local/bin/
```

Windows: download `crosscell-windows-x64.zip` from the release page and add to PATH.

### From Source

Requires [Rust toolchain](https://rustup.rs) and CMake ≥ 3.26.

```bash
git clone https://github.com/biodancerwangzhi/crosscell.git
cd crosscell
cargo build --release
# Binary at target/release/crosscell
```

### Verify

```bash
crosscell --version
crosscell --help
```

## Python

### From PyPI

```bash
pip install crosscell
```

### From Source

```bash
pip install maturin
git clone https://github.com/biodancerwangzhi/crosscell.git
cd crosscell/crosscell-py
maturin develop --release
```

### Verify

```python
import crosscell
print(crosscell.__version__)
```

## R

### From R-universe

```r
install.packages("crosscell", repos = "https://biodancerwangzhi.r-universe.dev")
```

### From Source

Requires Rust toolchain and CMake ≥ 3.26.

```r
# Install from local clone
install.packages("path/to/crosscell-r", repos = NULL, type = "source")
```

### Verify

```r
library(crosscell)
info <- inspect_file("test.h5ad")
```

## Dependencies

| Interface | Requirements |
|-----------|-------------|
| CLI | None (statically linked) |
| Python | Python 3.9+, anndata, scipy, numpy, pandas |
| R | R 4.0+, Matrix; optional: Seurat, SingleCellExperiment |
