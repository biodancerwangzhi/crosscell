#!/usr/bin/env python3
import anndata

ad = anndata.read_h5ad('tests/data/rust_written_metadata.h5ad')
print('✅ Successfully read Rust-written .h5ad file!')
print(f'Shape: {ad.shape}')
print(f'obs columns: {list(ad.obs.columns)}')
print(f'var columns: {list(ad.var.columns)}')
print(f'cell_type dtype: {ad.obs["cell_type"].dtype}')
print(f'cell_type categories: {list(ad.obs["cell_type"].cat.categories)}')
print(f'n_genes values: {ad.obs["n_genes"].values}')
print(f'cell_id values: {list(ad.obs["cell_id"].values)}')
