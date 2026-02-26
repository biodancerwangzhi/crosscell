//! High-performance sparse matrix subset algorithms

use crate::ir::expression::SparseMatrixCSR;
use std::collections::HashMap;
use std::ptr;

pub fn is_contiguous(indices: &[usize]) -> bool {
    if indices.len() <= 1 {
        return true;
    }
    indices.windows(2).all(|w| w[1] == w[0] + 1)
}

fn create_empty_csr_matrix(n_rows: usize, n_cols: usize) -> SparseMatrixCSR {
    SparseMatrixCSR {
        data: Vec::new(),
        indices: Vec::new(),
        indptr: vec![0; n_rows + 1],
        n_rows,
        n_cols,
    }
}

pub fn subset_rows_only_in_place(
    row_offsets: Vec<usize>,
    mut col_indices: Vec<usize>,
    mut values: Vec<f64>,
    row_indices: &[usize],
    n_cols: usize,
) -> SparseMatrixCSR {
    if row_indices.is_empty() {
        return create_empty_csr_matrix(0, n_cols);
    }

    let total_nnz: usize = row_indices
        .iter()
        .map(|&row_idx| row_offsets[row_idx + 1] - row_offsets[row_idx])
        .sum();

    let mut new_row_offsets = Vec::with_capacity(row_indices.len() + 1);
    new_row_offsets.push(0);

    let mut write_pos = 0;

    for &row_idx in row_indices {
        let start = row_offsets[row_idx];
        let end = row_offsets[row_idx + 1];
        let row_nnz = end - start;

        if write_pos != start {
            unsafe {
                ptr::copy(
                    col_indices.as_ptr().add(start),
                    col_indices.as_mut_ptr().add(write_pos),
                    row_nnz,
                );
                ptr::copy(
                    values.as_ptr().add(start),
                    values.as_mut_ptr().add(write_pos),
                    row_nnz,
                );
            }
        }

        write_pos += row_nnz;
        new_row_offsets.push(write_pos);
    }

    col_indices.truncate(total_nnz);
    values.truncate(total_nnz);

    if col_indices.capacity() > total_nnz * 2 {
        col_indices.shrink_to_fit();
        values.shrink_to_fit();
    }

    SparseMatrixCSR {
        data: values,
        indices: col_indices,
        indptr: new_row_offsets,
        n_rows: row_indices.len(),
        n_cols,
    }
}

pub fn subset_contiguous_columns_in_place(
    row_offsets: Vec<usize>,
    mut col_indices: Vec<usize>,
    mut values: Vec<f64>,
    row_indices: &[usize],
    col_range: &[usize],
) -> SparseMatrixCSR {
    if row_indices.is_empty() || col_range.is_empty() {
        return create_empty_csr_matrix(row_indices.len(), col_range.len());
    }

    let col_start = col_range[0];
    let col_end = col_range[col_range.len() - 1] + 1;

    let mut new_row_offsets = Vec::with_capacity(row_indices.len() + 1);
    new_row_offsets.push(0);

    let mut write_pos = 0;

    for &row_idx in row_indices {
        let start = row_offsets[row_idx];
        let end = row_offsets[row_idx + 1];

        let mut row_write_pos = write_pos;
        for read_pos in start..end {
            let col = col_indices[read_pos];
            if col >= col_start && col < col_end {
                let new_col = col - col_start;
                col_indices[row_write_pos] = new_col;
                if row_write_pos != read_pos {
                    values[row_write_pos] = unsafe { ptr::read(&values[read_pos]) };
                }
                row_write_pos += 1;
            }
        }
        write_pos = row_write_pos;
        new_row_offsets.push(write_pos);
    }

    let final_nnz = write_pos;
    col_indices.truncate(final_nnz);
    values.truncate(final_nnz);

    if col_indices.capacity() > final_nnz * 2 {
        col_indices.shrink_to_fit();
        values.shrink_to_fit();
    }

    SparseMatrixCSR {
        data: values,
        indices: col_indices,
        indptr: new_row_offsets,
        n_rows: row_indices.len(),
        n_cols: col_range.len(),
    }
}

pub fn subset_sparse_columns_in_place(
    row_offsets: Vec<usize>,
    col_indices: Vec<usize>,
    values: Vec<f64>,
    row_indices: &[usize],
    target_cols: &[usize],
) -> SparseMatrixCSR {
    if row_indices.is_empty() || target_cols.is_empty() {
        return create_empty_csr_matrix(row_indices.len(), target_cols.len());
    }

    let max_col = target_cols.iter().max().copied().unwrap_or(0);
    let density = target_cols.len() as f64 / (max_col + 1) as f64;

    if target_cols.len() < 512 && max_col < 8192 && density < 0.1 {
        subset_sparse_with_dense_lookup(
            row_offsets,
            col_indices,
            values,
            row_indices,
            target_cols,
            max_col,
        )
    } else if target_cols.len() < 2048 {
        subset_sparse_with_binary_search(row_offsets, col_indices, values, row_indices, target_cols)
    } else {
        subset_sparse_with_hash_lookup(row_offsets, col_indices, values, row_indices, target_cols)
    }
}

fn subset_sparse_with_dense_lookup(
    row_offsets: Vec<usize>,
    mut col_indices: Vec<usize>,
    mut values: Vec<f64>,
    row_indices: &[usize],
    target_cols: &[usize],
    max_col: usize,
) -> SparseMatrixCSR {
    let mut lookup = vec![None; max_col + 1];
    for (new_idx, &old_idx) in target_cols.iter().enumerate() {
        lookup[old_idx] = Some(new_idx);
    }

    let mut new_row_offsets = Vec::with_capacity(row_indices.len() + 1);
    new_row_offsets.push(0);
    let mut write_pos = 0;

    for &row_idx in row_indices {
        let start = row_offsets[row_idx];
        let end = row_offsets[row_idx + 1];

        let mut row_write_pos = write_pos;
        for read_pos in start..end {
            let old_col = col_indices[read_pos];
            if let Some(new_col) = lookup.get(old_col).and_then(|&x| x) {
                col_indices[row_write_pos] = new_col;
                if row_write_pos != read_pos {
                    values[row_write_pos] = unsafe { ptr::read(&values[read_pos]) };
                }
                row_write_pos += 1;
            }
        }
        write_pos = row_write_pos;
        new_row_offsets.push(write_pos);
    }

    col_indices.truncate(write_pos);
    values.truncate(write_pos);

    SparseMatrixCSR {
        data: values,
        indices: col_indices,
        indptr: new_row_offsets,
        n_rows: row_indices.len(),
        n_cols: target_cols.len(),
    }
}

fn subset_sparse_with_binary_search(
    row_offsets: Vec<usize>,
    mut col_indices: Vec<usize>,
    mut values: Vec<f64>,
    row_indices: &[usize],
    target_cols: &[usize],
) -> SparseMatrixCSR {
    let mut sorted_mapping: Vec<(usize, usize)> = target_cols
        .iter()
        .enumerate()
        .map(|(new_idx, &old_idx)| (old_idx, new_idx))
        .collect();
    sorted_mapping.sort_unstable_by_key(|&(old_idx, _)| old_idx);

    let mut new_row_offsets = Vec::with_capacity(row_indices.len() + 1);
    new_row_offsets.push(0);
    let mut write_pos = 0;

    for &row_idx in row_indices {
        let start = row_offsets[row_idx];
        let end = row_offsets[row_idx + 1];

        let mut row_write_pos = write_pos;
        for read_pos in start..end {
            let old_col = col_indices[read_pos];
            if let Ok(found_idx) =
                sorted_mapping.binary_search_by_key(&old_col, |&(old_idx, _)| old_idx)
            {
                let new_col = sorted_mapping[found_idx].1;
                col_indices[row_write_pos] = new_col;
                if row_write_pos != read_pos {
                    values[row_write_pos] = unsafe { ptr::read(&values[read_pos]) };
                }
                row_write_pos += 1;
            }
        }
        write_pos = row_write_pos;
        new_row_offsets.push(write_pos);
    }

    col_indices.truncate(write_pos);
    values.truncate(write_pos);

    SparseMatrixCSR {
        data: values,
        indices: col_indices,
        indptr: new_row_offsets,
        n_rows: row_indices.len(),
        n_cols: target_cols.len(),
    }
}

fn subset_sparse_with_hash_lookup(
    row_offsets: Vec<usize>,
    mut col_indices: Vec<usize>,
    mut values: Vec<f64>,
    row_indices: &[usize],
    target_cols: &[usize],
) -> SparseMatrixCSR {
    let col_map: HashMap<usize, usize> = target_cols
        .iter()
        .enumerate()
        .map(|(new_idx, &old_idx)| (old_idx, new_idx))
        .collect();

    let mut new_row_offsets = Vec::with_capacity(row_indices.len() + 1);
    new_row_offsets.push(0);
    let mut write_pos = 0;

    for &row_idx in row_indices {
        let start = row_offsets[row_idx];
        let end = row_offsets[row_idx + 1];

        let mut row_write_pos = write_pos;
        for read_pos in start..end {
            let old_col = col_indices[read_pos];
            if let Some(&new_col) = col_map.get(&old_col) {
                col_indices[row_write_pos] = new_col;
                if row_write_pos != read_pos {
                    values[row_write_pos] = unsafe { ptr::read(&values[read_pos]) };
                }
                row_write_pos += 1;
            }
        }
        write_pos = row_write_pos;
        new_row_offsets.push(write_pos);
    }

    col_indices.truncate(write_pos);
    values.truncate(write_pos);

    SparseMatrixCSR {
        data: values,
        indices: col_indices,
        indptr: new_row_offsets,
        n_rows: row_indices.len(),
        n_cols: target_cols.len(),
    }
}
