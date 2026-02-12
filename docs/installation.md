# Installation

## Docker (Recommended)

```bash
# Clone repository
git clone https://github.com/yourusername/crosscell.git
cd crosscell

# Build Docker image
docker-compose build dev

# Build release binary (static HDF5 linking)
docker-compose run --rm dev cargo build --release

# Binary is generated at target/release/crosscell
# Run on Linux/WSL
chmod +x target/release/crosscell
./target/release/crosscell --version
```

## From Source

### Requirements

- Rust 1.75+
- HDF5 library
- C compiler (gcc/clang)

### Build

```bash
# Clone repository
git clone https://github.com/yourusername/crosscell.git
cd crosscell

# Build release binary
cargo build --release

# Binary location
./target/release/crosscell --version
```

### Install to PATH

```bash
# Linux/macOS
sudo cp target/release/crosscell /usr/local/bin/

# Or add to PATH
export PATH="$PATH:$(pwd)/target/release"
```

## Verify Installation

```bash
# Check version
crosscell --version

# Show help
crosscell --help

# Test with sample file
crosscell inspect -i sample.h5ad
```

## Docker Commands

All commands in this documentation can be run with Docker:

```bash
# Pattern
docker-compose run --rm dev crosscell <command>

# Examples
docker-compose run --rm dev crosscell convert -i input.h5ad -o output.rds -f seurat
docker-compose run --rm dev crosscell inspect -i data.h5ad
```
