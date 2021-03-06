"""
Unit and regression test for the parsnip package.
"""

# Import package, test suite, and other packages as needed
import parsnip
import pytest
import sys

def test_parsnip_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "parsnip" in sys.modules
