#!/bin/bash
# Test old-format H5AD files that were previously failing
set -e

FILES="scanpy_pbmc3k scanpy_pbmc3k_processed scvelo_pancreas scvelo_dentategyrus"

for f in $FILES; do
    echo "=== Testing $f ==="
    INPUT="/workspace/data/generated/${f}.h5ad"
    OUTPUT="/tmp/test_${f}.rds"
    
    if [ ! -f "$INPUT" ]; then
        echo "  SKIP: $INPUT not found"
        echo
        continue
    fi
    
    ./benchmark/crosscell convert -i "$INPUT" -o "$OUTPUT" -f seurat 2>&1 | tail -8
    echo
done

echo "=== All tests complete ==="
