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

The Python API now supports the same AI-Ready transforms as the CLI:

```python
import crosscell

# AI-Ready pipeline in Python (no subprocess needed)
adata = crosscell.read_rds(
    "seurat.rds",
    normalize=True,
    top_genes=2000,
    gene_id_column="ensembl_id",
)

# Validate roundtrip
result = crosscell.validate("original.h5ad", "roundtrip.h5ad", cluster_column="cell_type")
print(f"Passed: {result['passed']}, ARI: {result.get('ari', 'N/A')}")
```

See [Python API Reference](python_api.md) for full details.

## R Integration

The R API also supports AI-Ready transforms:

```r
library(crosscell)

# AI-Ready export from R
write_as_h5ad(seurat_obj, "ai_ready.h5ad",
              normalize = TRUE, top_genes = 2000L, gene_id_column = "ensembl_id")

# Validate roundtrip
result <- validate_conversion("original.rds", "converted.h5ad", cluster_column = "cell_type")
cat(sprintf("Passed: %s, ARI: %.4f\n", result$passed, result$ari))
```

See [R API Reference](r_api.md) for full details.

## Spatial Transcriptomics

```bash
# Convert Visium data
crosscell convert -i visium.h5ad -o visium.rds -f seurat

# Validate spatial coordinates
crosscell validate --original visium.h5ad --converted visium_roundtrip.h5ad
```

## AI-Ready Export

CrossCell can prepare data for AI foundation models (Geneformer, scGPT, scBERT, etc.) in a single conversion step. Three transform flags can be combined freely and are applied in order: normalize → top-genes → gene-id-column.

### Library-Size Normalization

Applies per-cell library-size normalization (scale to 10,000) followed by log1p transform:

```bash
crosscell convert -i raw.h5ad -o normalized.h5ad -f anndata --normalize
```

All-zero cells (empty droplets) are left as zero to avoid division errors.

### Highly Variable Gene Selection

Keeps only the top N genes ranked by variance:

```bash
crosscell convert -i raw.h5ad -o hvg.h5ad -f anndata --top-genes 2000
```

If N exceeds the total gene count, all genes are kept and a warning is emitted.

### Gene Identifier Replacement

Replaces the gene index with values from a specified column in `var`:

```bash
crosscell convert -i raw.h5ad -o ensembl.h5ad -f anndata --gene-id-column ensembl_id
```

### Combined Example

A typical AI-ready pipeline in one command:

```bash
crosscell convert \
  -i seurat_obj.rds \
  -o ai_ready.h5ad \
  -f anndata \
  --normalize \
  --top-genes 2000 \
  --gene-id-column ensembl_id
```

## CELLxGENE Schema Validation

Validate that a dataset meets CELLxGENE Discover Schema 5.0 requirements before submission:

```bash
crosscell validate \
  --original data.h5ad \
  --converted data.h5ad \
  --validate-schema cellxgene-v5
```

This checks for 12 required `obs` fields (e.g. `assay_ontology_term_id`, `cell_type_ontology_term_id`, `donor_id`, `organism_ontology_term_id`, etc.) and 1 required `var` field (`feature_is_filtered`).

## Cluster-Aware Validation (ARI / NMI)

When validating roundtrip conversions, you can measure the impact on clustering by specifying a cluster label column:

```bash
crosscell validate \
  --original original.h5ad \
  --converted roundtrip.h5ad \
  --cluster-column cell_type
```

This adds Adjusted Rand Index (ARI) and Normalized Mutual Information (NMI) to the validation report alongside the standard metrics (Pearson correlation, MSE, exact match ratio, etc.).

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
| `--normalize` | false | Library-size normalization + log1p |
| `--top-genes` | | Select top N variable genes by variance |
| `--gene-id-column` | | Use specified var column as gene identifiers |
| `--keep-layers` | false | Preserve Seurat V5 split layers as separate AnnData layers |

## Validate Command Options

| Option | Default | Description |
|--------|---------|-------------|
| `--original` | required | Original file |
| `--converted` | required | Converted file |
| `--tolerance` | 1e-7 | Numerical tolerance |
| `-o, --output` | | Save report to JSON |
| `--strict` | false | Fail on any difference |
| `--cluster-column` | | obs column for cluster labels (enables ARI/NMI) |
| `--validate-schema` | | Schema to validate against (e.g. `cellxgene-v5`) |

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
