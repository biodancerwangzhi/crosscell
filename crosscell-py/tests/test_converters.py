"""Tests for IR ↔ AnnData converters."""
import numpy as np
import scipy.sparse as sp
import pandas as pd
import anndata


def test_import():
    """Test that the module can be imported."""
    import crosscell
    assert hasattr(crosscell, "__version__")
    print(f"crosscell version: {crosscell.__version__}")
