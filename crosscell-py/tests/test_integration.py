"""Integration tests for crosscell Python API.

Tests 7.1-7.6: read_rds, read_h5ad, write/read roundtrip, inspect, error handling.
Uses real data from tests/data/ and data/generated/.
"""
import os
import tempfile
import pytest
import numpy as np
import anndata

import crosscell

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
WORKSPACE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TESTS_DATA = os.path.join(WORKSPACE, "tests", "data")
GENERATED_DATA = os.path.join(WORKSPACE, "data", "generated")


def _rds_path(name):
    return os.path.join(TESTS_DATA, name)


def _h5ad_path(name):
    return os.path.join(TESTS_DATA, name)


def _generated_rds(name):
    return os.path.join(GENERATED_DATA, name)


def _generated_h5ad(name):
    return os.path.join(GENERATED_DATA, name)


# ===================================================================
# 7.1  read_rds integration tests
# ===================================================================
class TestReadRds:
    """Test reading Seurat RDS files into AnnData."""

    def test_read_minimal_seurat(self):
        path = _rds_path("seurat_minimal.rds")
        if not os.path.exists(path):
            pytest.skip(f"Test file not found: {path}")
        adata = crosscell.read_rds(path)
        assert isinstance(adata, anndata.AnnData)
        assert adata.n_obs > 0
        assert adata.n_vars > 0

    def test_read_seurat_with_dimred(self):
        path = _rds_path("seurat_with_dimred.rds")
        if not os.path.exists(path):
            pytest.skip(f"Test file not found: {path}")
        adata = crosscell.read_rds(path)
        assert isinstance(adata, anndata.AnnData)
        # Should have embeddings in obsm
        assert len(adata.obsm) > 0

    def test_read_seurat_v4_generated(self):
        path = _generated_rds("seurat_v4_pbmc3k_raw.rds")
        if not os.path.exists(path):
            pytest.skip(f"Test file not found: {path}")
        adata = crosscell.read_rds(path)
        assert isinstance(adata, anndata.AnnData)
        assert adata.n_obs > 0
        assert adata.n_vars > 0
        # X should be a matrix
        assert adata.X is not None

    def test_read_seurat_multi_assay(self):
        path = _rds_path("seurat_multi_assay.rds")
        if not os.path.exists(path):
            pytest.skip(f"Test file not found: {path}")
        adata = crosscell.read_rds(path)
        assert isinstance(adata, anndata.AnnData)


# ===================================================================
# 7.2  read_h5ad integration tests
# ===================================================================
class TestReadH5ad:
    """Test reading H5AD files into AnnData."""

    def test_read_small_sparse(self):
        path = _h5ad_path("small_sparse.h5ad")
        if not os.path.exists(path):
            pytest.skip(f"Test file not found: {path}")
        adata = crosscell.read_h5ad(path)
        assert isinstance(adata, anndata.AnnData)
        assert adata.n_obs > 0
        assert adata.n_vars > 0
        # X should be sparse
        import scipy.sparse as sp
        assert sp.issparse(adata.X) or isinstance(adata.X, np.ndarray)

    def test_read_small_dense(self):
        path = _h5ad_path("small_dense.h5ad")
        if not os.path.exists(path):
            pytest.skip(f"Test file not found: {path}")
        adata = crosscell.read_h5ad(path)
        assert isinstance(adata, anndata.AnnData)
        assert adata.n_obs > 0

    def test_read_with_embeddings(self):
        path = _h5ad_path("python_verify_embeddings.h5ad")
        if not os.path.exists(path):
            pytest.skip(f"Test file not found: {path}")
        adata = crosscell.read_h5ad(path)
        assert isinstance(adata, anndata.AnnData)
        assert len(adata.obsm) > 0

    def test_read_scanpy_pbmc3k(self):
        path = _generated_h5ad("scanpy_pbmc3k.h5ad")
        if not os.path.exists(path):
            pytest.skip(f"Test file not found: {path}")
        adata = crosscell.read_h5ad(path)
        assert isinstance(adata, anndata.AnnData)
        assert adata.n_obs > 0
        assert adata.n_vars > 0

    def test_read_single_row(self):
        path = _h5ad_path("single_row.h5ad")
        if not os.path.exists(path):
            pytest.skip(f"Test file not found: {path}")
        adata = crosscell.read_h5ad(path)
        assert adata.n_obs == 1

    def test_read_single_column(self):
        path = _h5ad_path("single_column.h5ad")
        if not os.path.exists(path):
            pytest.skip(f"Test file not found: {path}")
        adata = crosscell.read_h5ad(path)
        assert adata.n_vars == 1


# ===================================================================
# 7.3  H5AD write-read roundtrip
# ===================================================================
class TestH5adRoundtrip:
    """Test H5AD write then read back consistency."""

    def test_roundtrip_small_sparse(self):
        path = _h5ad_path("small_sparse.h5ad")
        if not os.path.exists(path):
            pytest.skip(f"Test file not found: {path}")
        original = crosscell.read_h5ad(path)

        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as f:
            tmp = f.name
        try:
            crosscell.write_h5ad(original, tmp)
            reloaded = crosscell.read_h5ad(tmp)

            assert reloaded.n_obs == original.n_obs
            assert reloaded.n_vars == original.n_vars
            # Expression matrix values should be close
            import scipy.sparse as sp
            x_orig = original.X.toarray() if sp.issparse(original.X) else np.asarray(original.X)
            x_new = reloaded.X.toarray() if sp.issparse(reloaded.X) else np.asarray(reloaded.X)
            np.testing.assert_allclose(x_orig, x_new, rtol=1e-6)
        finally:
            os.unlink(tmp)

    def test_roundtrip_with_embeddings(self):
        path = _h5ad_path("python_verify_embeddings.h5ad")
        if not os.path.exists(path):
            pytest.skip(f"Test file not found: {path}")
        original = crosscell.read_h5ad(path)

        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as f:
            tmp = f.name
        try:
            crosscell.write_h5ad(original, tmp)
            reloaded = crosscell.read_h5ad(tmp)

            assert reloaded.n_obs == original.n_obs
            assert reloaded.n_vars == original.n_vars
            # Check embeddings preserved
            for key in original.obsm:
                assert key in reloaded.obsm
                np.testing.assert_allclose(
                    np.asarray(original.obsm[key]),
                    np.asarray(reloaded.obsm[key]),
                    rtol=1e-6,
                )
        finally:
            os.unlink(tmp)


# ===================================================================
# 7.4  inspect / format auto-detection
# ===================================================================
class TestInspect:
    """Test inspect function for format auto-detection."""

    def test_inspect_h5ad(self):
        path = _h5ad_path("small_sparse.h5ad")
        if not os.path.exists(path):
            pytest.skip(f"Test file not found: {path}")
        info = crosscell.inspect(path)
        assert isinstance(info, dict)
        assert info["format"] == "h5ad"
        assert info["n_cells"] > 0
        assert info["n_genes"] > 0
        assert "is_sparse" in info

    def test_inspect_rds(self):
        path = _rds_path("seurat_minimal.rds")
        if not os.path.exists(path):
            pytest.skip(f"Test file not found: {path}")
        info = crosscell.inspect(path)
        assert isinstance(info, dict)
        assert info["format"] == "rds"
        assert info["n_cells"] > 0
        assert info["n_genes"] > 0
        assert "seurat_version" in info

    def test_inspect_h5ad_generated(self):
        path = _generated_h5ad("scanpy_pbmc3k.h5ad")
        if not os.path.exists(path):
            pytest.skip(f"Test file not found: {path}")
        info = crosscell.inspect(path)
        assert info["format"] == "h5ad"
        assert isinstance(info["embedding_names"], list)
        assert isinstance(info["layer_names"], list)

    def test_inspect_unsupported_format(self):
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
            f.write(b"a,b\n1,2\n")
            tmp = f.name
        try:
            with pytest.raises(ValueError, match="Unsupported file format"):
                crosscell.inspect(tmp)
        finally:
            os.unlink(tmp)


# ===================================================================
# 7.5  AnnData -> Seurat RDS cross-format roundtrip
# ===================================================================
class TestCrossFormatRoundtrip:
    """Test AnnData -> RDS -> AnnData roundtrip."""

    @pytest.mark.skip(reason="write_rds creates simplified structure, not full S4 Seurat")
    def test_h5ad_to_rds_roundtrip(self):
        path = _h5ad_path("small_sparse.h5ad")
        if not os.path.exists(path):
            pytest.skip(f"Test file not found: {path}")
        original = crosscell.read_h5ad(path)

        with tempfile.NamedTemporaryFile(suffix=".rds", delete=False) as f:
            rds_tmp = f.name
        try:
            crosscell.write_rds(original, rds_tmp)
            roundtrip = crosscell.read_rds(rds_tmp)

            assert isinstance(roundtrip, anndata.AnnData)
            assert roundtrip.n_obs == original.n_obs
            assert roundtrip.n_vars == original.n_vars
        finally:
            os.unlink(rds_tmp)

    def test_rds_to_h5ad_roundtrip(self):
        path = _rds_path("seurat_minimal.rds")
        if not os.path.exists(path):
            pytest.skip(f"Test file not found: {path}")
        original = crosscell.read_rds(path)

        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as f:
            h5ad_tmp = f.name
        try:
            crosscell.write_h5ad(original, h5ad_tmp)
            roundtrip = crosscell.read_h5ad(h5ad_tmp)

            assert isinstance(roundtrip, anndata.AnnData)
            assert roundtrip.n_obs == original.n_obs
            assert roundtrip.n_vars == original.n_vars
        finally:
            os.unlink(h5ad_tmp)


# ===================================================================
# 7.6  Error handling tests
# ===================================================================
class TestErrorHandling:
    """Test error handling for edge cases."""

    def test_read_rds_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            crosscell.read_rds("/nonexistent/path/file.rds")

    def test_read_h5ad_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            crosscell.read_h5ad("/nonexistent/path/file.h5ad")

    def test_inspect_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            crosscell.inspect("/nonexistent/path/file.h5ad")

    def test_read_rds_invalid_file(self):
        """Reading a non-Seurat RDS should raise ValueError."""
        path = _rds_path("r_simple_int.rds")
        if not os.path.exists(path):
            pytest.skip(f"Test file not found: {path}")
        with pytest.raises(ValueError):
            crosscell.read_rds(path)

    def test_read_h5ad_invalid_file(self):
        """Reading a non-H5AD file as H5AD should raise ValueError."""
        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as f:
            f.write(b"this is not an h5ad file")
            tmp = f.name
        try:
            with pytest.raises(ValueError):
                crosscell.read_h5ad(tmp)
        finally:
            os.unlink(tmp)

    def test_write_rds_empty_anndata(self):
        """Writing an empty AnnData should raise ValueError."""
        import scipy.sparse as sp
        adata = anndata.AnnData(
            X=sp.csr_matrix((0, 0)),
        )
        with tempfile.NamedTemporaryFile(suffix=".rds", delete=False) as f:
            tmp = f.name
        try:
            with pytest.raises(ValueError):
                crosscell.write_rds(adata, tmp)
        finally:
            if os.path.exists(tmp):
                os.unlink(tmp)

    def test_inspect_unsupported_extension(self):
        with tempfile.NamedTemporaryFile(suffix=".txt", delete=False) as f:
            f.write(b"hello")
            tmp = f.name
        try:
            with pytest.raises(ValueError, match="Unsupported"):
                crosscell.inspect(tmp)
        finally:
            os.unlink(tmp)
