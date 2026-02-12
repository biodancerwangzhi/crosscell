# FAQ

## General

### What formats are supported?

- **AnnData** (.h5ad) - Python/Scanpy format
- **Seurat** (.rds) - R/Seurat format (V3, V4, V5)

### Is Python/R required?

No. CrossCell is a standalone Rust binary with no runtime dependencies.

### What data is preserved?

Expression matrices, cell/gene metadata, embeddings (PCA, UMAP, tSNE), layers, categorical variables, and spatial coordinates are fully preserved.

---

## Conversion Issues

### "Out of memory" error

Use lazy loading for large datasets:

```bash
crosscell convert -i large.h5ad -o large.rds -f seurat --lazy --chunk-size 5000
```

Or limit threads:

```bash
crosscell convert -i large.h5ad -o large.rds -f seurat -t 4
```

### "Unsupported format" error

Ensure the file extension matches the content:
- `.h5ad` for AnnData
- `.rds` for Seurat

### Seurat conversion fails

Complex Seurat objects may contain runtime objects (closures, environments). Simplify first:

```bash
crosscell simplify-seurat -i complex.rds -o simplified.rds
crosscell convert -i simplified.rds -o output.h5ad -f anndata
```

### Validation fails

Try relaxing the tolerance:

```bash
crosscell validate --original a.h5ad --converted b.h5ad --tolerance 1e-5
```

---

## Performance

### How to speed up conversion?

1. Use LZ4 compression: `-c lz4`
2. Disable compression: `-c none`
3. Use all CPU cores (default)

### How to reduce memory usage?

1. Enable lazy loading: `--lazy`
2. Reduce chunk size: `--chunk-size 2000`
3. Limit threads: `-t 4`

### Estimated conversion time?

| Dataset Size | Time |
|--------------|------|
| 3k cells | ~5 seconds |
| 20k cells | ~30 seconds |
| 100k cells | ~3 minutes |

---

## Data Integrity

### Is data loss possible?

Some `.uns` (unstructured) data may not be preserved. Core data (expression, metadata, embeddings) is fully preserved.

### How to verify conversion?

```bash
crosscell convert -i input.h5ad -o output.rds -f seurat --validate
```

Or manually:

```bash
crosscell validate --original input.h5ad --converted roundtrip.h5ad
```

### What tolerance is acceptable?

Default tolerance is `1e-7`. For most use cases, `1e-5` is acceptable.

---

## Docker

### Why use Docker?

Docker ensures consistent environment across platforms and avoids dependency issues.

### How to run commands?

```bash
docker-compose run --rm dev crosscell <command>
```

### How to access files?

The workspace is mounted at `/workspace`. Use relative paths:

```bash
docker-compose run --rm dev crosscell convert -i data/input.h5ad -o data/output.rds -f seurat
```
