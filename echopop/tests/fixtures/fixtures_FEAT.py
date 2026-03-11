"""Test fixtures for FEAT workflow."""

import pytest


@pytest.fixture
def tmp_excel(tmp_path):
    """Return a temporary Excel file path."""
    return tmp_path / "test.xlsx"
