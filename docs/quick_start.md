# Quick Start

## Basic Conversion

### AnnData to Seurat

```bash
crosscell convert -i input.h5ad -o output.rds -f seurat
```

### Seurat to AnnData

```bash
crosscell convert -i input.rds -o output.h5ad -f anndata
```

### With Validation

```bash
crosscell convert -i input.h5ad -o output.rds -f seurat --validate
```

## Inspect Files

```bash
# Basic info
crosscell inspect -i data.h5ad

# Detailed info
crosscell inspect -i data.h5ad --detailed
```

Output:
```
📊 File Information
Format: AnnData (h5ad)
Cells: 2,700
Genes: 32,738
Sparsity: 93.4%

✅ Data Components
  ├─ Expression Matrix (X): Sparse CSR
  ├─ Cell Metadata (obs): 5 columns
  ├─ Embeddings (obsm): PCA, UMAP
  └─ Layers: counts, log1p
```

## Large Datasets

### Streaming Mode (Recommended for 1M+ cells)

```bash
# Enable streaming with custom chunk size
crosscell convert -i huge.h5ad -o huge.rds -f seurat --streaming --chunk-size 10000
```

Output:
```
🔄 CrossCell Streaming Conversion
   huge.h5ad → huge.rds
   Format: ANNDATA → SEURAT
   Mode: Streaming (constant memory)
   Chunk size: 10000 cells

📊 File Information:
   File size: 15.2 GB
   Estimated cells: ~5,000,000
   Estimated memory: ~2.4 GB

   Progress: 100% (5,000,000/5,000,000 cells)

✅ Streaming conversion completed!
   Time: 245.32s
   Output size: 12.8 GB
```

### Auto-Detection

For files > 5GB, CrossCell automatically suggests streaming mode:

```bash
crosscell convert -i huge.h5ad -o huge.rds -f seurat
# ⚠ Large file detected (15 GB, ~5M cells)
# 💡 Using streaming mode with chunk_size=10000
```

## Validate Roundtrip

```bash
# Convert both ways
crosscell convert -i original.h5ad -o converted.rds -f seurat
crosscell convert -i converted.rds -o roundtrip.h5ad -f anndata

# Validate
crosscell validate --original original.h5ad --converted roundtrip.h5ad
```

## Simplify Seurat Objects

Complex Seurat objects may need simplification before conversion:

```bash
# Simplify
crosscell simplify-seurat -i complex.rds -o simplified.rds

# Then convert
crosscell convert -i simplified.rds -o output.h5ad -f anndata
```

## List Supported Formats

```bash
crosscell formats
```

Output:
```
📋 Supported Formats
─────────────────────────────────────────

  anndata      (.h5ad)  [read/write]
               AnnData format for Python/Scanpy
  seurat       (.rds)   [read/write]
               Seurat format for R
  loom         (.loom)  [read/write]
               Loom format
```

## Commands Reference

| Command | Description |
|---------|-------------|
| `convert` | Convert between formats |
| `inspect` | Show file information |
| `validate` | Validate conversion accuracy |
| `simplify-seurat` | Simplify Seurat objects |
| `formats` | List supported formats |

## Global Options

| Option | Description |
|--------|-------------|
| `-v, --verbose` | Verbose output |
| `-q, --quiet` | Errors only |
| `-t, --threads` | Thread count |
| `-h, --help` | Show help |
| `-V, --version` | Show version |
