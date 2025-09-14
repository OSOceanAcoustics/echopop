import pytest


@pytest.fixture
def tmp_excel(tmp_path):
    return tmp_path / "test.xlsx"
