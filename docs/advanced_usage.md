# Advanced Usage

## Large Datasets

CrossCell provides multiple modes for handling large datasets:

| Mode | Cells | Memory | Speed | Use Case |
|------|-------|--------|-------|----------|
| Standard | < 50k | High | Fastest | Small datasets |
| Lazy | 50k - 1M | Medium | Fast | Medium datasets |
| Streaming | 1M+ | Constant | Good | Large datasets |
| Disk-backed | 10M+ | Minimal | Slower | Ultra-large datasets |

### Memory Estimation

```bash
crosscell convert -i large.h5ad -o large.rds -f seurat --estimate-memory
```

Output:
```
📊 Memory Estimation:
   File size: 5.2 GB
   Full load: ~15.6 GB
   Lazy load: ~2.6 GB (peak)
   Streaming: ~2.4 GB (constant)

💡 Tip: Use --streaming flag for large datasets
```

### Streaming Mode

For datasets with 1M+ cells, streaming mode processes data in chunks with constant memory:

```bash
# Basic streaming
crosscell convert -i huge.h5ad -o huge.rds -f seurat --streaming

# Custom chunk size (default: 10000)
crosscell convert -i huge.h5ad -o huge.rds -f seurat --streaming --chunk-size 20000

# Dry run to preview
crosscell convert -i huge.h5ad -o huge.rds -f seurat --streaming --dry-run
```

Streaming supports:
- H5AD → H5AD (re-chunking)
- H5AD → RDS
- RDS → H5AD

### Lazy Loading

For datasets with 50k-1M cells:

```bash
crosscell convert -i large.h5ad -o large.rds -f seurat --lazy --chunk-size 5000
```

### Disk-Backed Mode

For ultra-large datasets (10M+ cells):

```bash
# Use disk cache
crosscell convert -i huge.h5ad -o huge.rds -f seurat --disk-backed

# Custom temp directory
crosscell convert -i huge.h5ad -o huge.rds -f seurat --disk-backed --temp-dir /fast/ssd
```

### Compression Options

```bash
# Fast (LZ4)
crosscell convert -i input.h5ad -o output.rds -f seurat -c lz4

# Maximum compression
crosscell convert -i input.h5ad -o output.rds -f seurat -c gzip --compression-level 9

# No compression (fastest)
crosscell convert -i input.h5ad -o output.rds -f seurat -c none
```

## Batch Processing

```bash
#!/bin/bash
for file in data/*.h5ad; do
    name=$(basename "$file" .h5ad)
    crosscell convert -i "$file" -o "output/${name}.rds" -f seurat --validate
done
```

## Python Integration

```python
import subprocess
import scanpy as sc

# Process in Python
adata = sc.read_h5ad('raw.h5ad')
sc.pp.normalize_total(adata)
sc.tl.pca(adata)
adata.write('processed.h5ad')

# Convert to Seurat
subprocess.run(['crosscell', 'convert', '-i', 'processed.h5ad', '-o', 'for_r.rds', '-f', 'seurat'])
```

## R Integration

```r
# After conversion from Python
seurat_obj <- readRDS('from_python.rds')
seurat_obj <- FindClusters(seurat_obj)
saveRDS(seurat_obj, 'analyzed.rds')

# Convert back (run in shell)
# crosscell simplify-seurat -i analyzed.rds -o simplified.rds
# crosscell convert -i simplified.rds -o back_to_python.h5ad -f anndata
```

## Spatial Transcriptomics

```bash
# Convert Visium data
crosscell convert -i visium.h5ad -o visium.rds -f seurat

# Validate spatial coordinates
crosscell validate --original visium.h5ad --converted visium_roundtrip.h5ad
```

## Convert Command Options

| Option | Default | Description |
|--------|---------|-------------|
| `-i, --input` | required | Input file |
| `-o, --output` | required | Output file |
| `-f, --format` | required | Target format (anndata/seurat/loom) |
| `--validate` | false | Validate after conversion |
| `--sparse-format` | auto | Sparse format (csr/csc/auto) |
| `-c, --compression` | gzip | Compression (none/gzip/lz4) |
| `--compression-level` | 6 | Compression level (0-9) |
| `--lazy` | false | Enable lazy loading |
| `--streaming` | false | Enable streaming mode |
| `--chunk-size` | 10000 | Chunk size for lazy/streaming |
| `--disk-backed` | false | Enable disk-backed mode |
| `--temp-dir` | system | Temp directory for disk-backed |
| `--estimate-memory` | false | Estimate memory usage |
| `--dry-run` | false | Preview without converting |
| `--direct` | true | Direct read for Seurat (default) |
| `--simplify-first` | false | Legacy mode requiring R preprocessing |

## Validate Command Options

| Option | Default | Description |
|--------|---------|-------------|
| `--original` | required | Original file |
| `--converted` | required | Converted file |
| `--tolerance` | 1e-7 | Numerical tolerance |
| `-o, --output` | | Save report to JSON |
| `--strict` | false | Fail on any difference |

## Environment Variables

| Variable | Description |
|----------|-------------|
| `CROSSCELL_THREADS` | Default thread count |
| `CROSSCELL_LOG_LEVEL` | Log level (error/warn/info/debug) |
| `RAYON_NUM_THREADS` | Rayon thread pool size |

## Exit Codes

| Code | Description |
|------|-------------|
| 0 | Success |
| 1 | General error |
| 2 | Argument error |
| 3 | File not found |
| 4 | Unsupported format |
| 5 | Out of memory |
| 6 | Validation failed |
| 7 | Conversion failed |

## Performance Tips

1. **Use streaming for large files**: Files > 5GB automatically suggest streaming mode
2. **Adjust chunk size**: Larger chunks = faster but more memory
3. **Use SSD for disk-backed mode**: Significantly faster than HDD
4. **Enable validation selectively**: Skip `--validate` for trusted conversions
5. **Use LZ4 compression**: Best speed/size tradeoff for most cases
