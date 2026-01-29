"""
pytest configuration for SparseSolv tests
"""

import pytest
import sys


def pytest_configure(config):
    """Configure pytest"""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )


def pytest_collection_modifyitems(config, items):
    """Skip tests if sparsesolv module is not available"""
    try:
        import sparsesolv
    except ImportError:
        skip_no_module = pytest.mark.skip(reason="sparsesolv module not built")
        for item in items:
            item.add_marker(skip_no_module)
