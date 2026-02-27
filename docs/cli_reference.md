# CLI Reference

## Global Options

| Option | Short | Description |
|--------|-------|-------------|
| `--verbose` | `-v` | Enable debug-level output |
| `--quiet` | `-q` | Errors only (conflicts with `-v`) |
| `--threads` | `-t` | Number of threads (default: all cores) |
| `--log-file` | | Write logs to file |
| `--version` | `-V` | Print version |
| `--help` | `-h` | Print help |

## Commands

### convert

Convert between H5AD and RDS formats.

```bash
crosscell convert -i <INPUT> -o <OUTPUT> -f <FORMAT> [OPTIONS]
```

#### Required

| Option | Short | Description |
|--------|-------|-------------|
| `--input` | `-i` | Input file path (.h5ad or .rds) |
| `--output` | `-o` | Output file path |
| `--format` | `-f` | Target format: `anndata`, `seurat` |

#### AI-Ready Transforms

Applied in order: normalize → top-genes → gene-id-column.

| Option | Default | Description |
|--------|---------|-------------|
| `--normalize` | false | Library-size normalization (scale to 10k) + log1p |
| `--top-genes` | — | Select top N genes by variance |
| `--gene-id-column` | — | Replace gene index with values from this var column |

#### Seurat V5

| Option | Default | Description |
|--------|---------|-------------|
| `--keep-layers` | false | Preserve split layers as separate AnnData layers instead of merging |
| `--direct` | true | Direct Rust-native RDS parsing (default) |
| `--simplify-first` | false | Legacy mode: requires R-preprocessed simplified RDS |

#### Quality

| Option | Default | Description |
|--------|---------|-------------|
| `--validate` | false | Run roundtrip validation after conversion |
| `--dry-run` | false | Preview conversion without writing output |
| `--estimate-memory` | false | Print memory usage estimate |

#### Performance

| Option | Default | Description |
|--------|---------|-------------|
| `--streaming` | false | Constant-memory streaming mode (1M+ cells) |
| `--lazy` | false | Lazy loading mode (50k–1M cells) |
| `--chunk-size` | 10000 | Rows per chunk for streaming/lazy |
| `--disk-backed` | false | Disk-cached mode for 10M+ cells |
| `--temp-dir` | system | Temp directory for disk-backed mode |

#### Compression

| Option | Default | Description |
|--------|---------|-------------|
| `--compression` / `-c` | gzip | Algorithm: `none`, `gzip`, `lz4` |
| `--compression-level` | 6 | Level 0–9 (gzip only) |

#### Sparse Matrix

| Option | Default | Description |
|--------|---------|-------------|
| `--sparse-format` | auto | Output sparse format: `csr`, `csc`, `auto` |

#### Examples

```bash
# Basic conversion
crosscell convert -i pbmc3k.h5ad -o pbmc3k.rds -f seurat

# AI-Ready pipeline
crosscell convert -i seurat.rds -o ai_ready.h5ad -f anndata \
  --normalize --top-genes 2000 --gene-id-column ensembl_id

# Keep Seurat V5 split layers
crosscell convert -i integrated.rds -o integrated.h5ad -f anndata --keep-layers

# Streaming for large datasets
crosscell convert -i big.h5ad -o big.rds -f seurat --streaming --chunk-size 20000

# Conversion with validation
crosscell convert -i input.h5ad -o output.rds -f seurat --validate
```

---

### inspect

Display file metadata.

```bash
crosscell inspect -i <INPUT> [OPTIONS]
```

| Option | Short | Description |
|--------|-------|-------------|
| `--input` | `-i` | Input file path (.h5ad or .rds) |
| `--detailed` | `-d` | Show detailed information |
| `--output` | `-o` | Save report as JSON |

#### Example

```bash
crosscell inspect -i pbmc3k.h5ad --detailed
```

---

### validate

Compare original and converted files for fidelity.

```bash
crosscell validate --original <FILE> --converted <FILE> [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--original` | required | Original file path |
| `--converted` | required | Converted file path |
| `--tolerance` | 1e-7 | Numerical tolerance for float comparison |
| `--output` / `-o` | — | Save report as JSON |
| `--strict` | false | Exit non-zero on any difference |
| `--cluster-column` | — | obs column for cluster labels (enables ARI/NMI) |
| `--validate-schema` | — | Validate against schema (e.g. `cellxgene-v5`) |

#### Metrics reported

Pearson R, Spearman ρ, MSE, RMSE, Max Diff, Mean Diff, Exact Match %, and optionally ARI + NMI.

#### Example

```bash
crosscell validate \
  --original original.h5ad \
  --converted roundtrip.h5ad \
  --cluster-column cell_type \
  --validate-schema cellxgene-v5
```

---

### simplify-seurat

Remove runtime objects (closures, environments, bytecode) from a Seurat RDS file.

```bash
crosscell simplify-seurat -i <INPUT> -o <OUTPUT> [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--input` / `-i` | required | Input Seurat RDS file |
| `--output` / `-o` | required | Output simplified RDS file |
| `--keep-graphs` | false | Preserve graph objects |
| `--keep-commands` | false | Preserve command history |
| `--keep-tools` | false | Preserve tool functions |

---

### formats

List all supported formats.

```bash
crosscell formats
```

## Environment Variables

| Variable | Description |
|----------|-------------|
| `CROSSCELL_THREADS` | Default thread count |
| `CROSSCELL_LOG_LEVEL` | Log level: error, warn, info, debug, trace |
| `RAYON_NUM_THREADS` | Rayon thread pool size |

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | General error |
| 2 | Argument error |
| 3 | File not found |
| 4 | Unsupported format |
| 5 | Out of memory |
| 6 | Validation failed |
| 7 | Conversion failed |
| 8 | User interrupt (Ctrl+C) |
