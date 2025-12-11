import contextlib
import os
import signal
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd
import psycopg
import pytest
import testing.postgresql

HERE = Path(__file__).parent.absolute()
TEST_DATA_ROOT = HERE.parent / "test_data"
TEST_SQL_FILE = TEST_DATA_ROOT / "Biological" / "test_bio_data.sql"


def _create_creds_dict(db_url):
    """Helper function to parse the DB URL into a dict."""
    url = urlparse(db_url)
    return {
        "dbname": url.path.lstrip("/"),
        "user": url.username,
        "password": url.password,
        "host": url.hostname,
        "port": url.port,
    }


@contextlib.contextmanager
def get_postgresql_context(original_sigint):
    """
    Context manager that either connects to an existing CI server
    or starts a local testing.postgresql instance.
    """

    is_github_action = os.environ.get("GITHUB_ACTIONS")

    if is_github_action:

        class CI_Postgres_Wrapper:
            def url(self):
                return "postgresql://test_user:postgres@localhost:5432/test"

        yield CI_Postgres_Wrapper()

    else:
        if os.name == "nt":
            signal.SIGINT = signal.SIGTERM

        try:
            with testing.postgresql.Postgresql() as postgresql:
                yield postgresql

        finally:
            if os.name == "nt":
                signal.SIGINT = original_sigint


@pytest.fixture(scope="session")
def database_credentials():
    """
    Session-scoped fixture to:
    1. Connect to CI DB or Start a local Postgres database.
    2. Load 'test_bio_data.sql' into it.
    3. Yield the *credentials dictionary* for this test DB.
    """
    original_sigint = signal.SIGINT

    with get_postgresql_context(original_sigint) as postgresql:

        try:
            conn = psycopg.connect(postgresql.url())
            conn.autocommit = True
            with conn.cursor() as cur:
                with open(TEST_SQL_FILE, "r") as f:
                    cur.execute(f.read())
            conn.close()
        except Exception as e:
            pytest.fail(f"Failed to load {TEST_SQL_FILE}: {e}")

        creds_dict = _create_creds_dict(postgresql.url())

        yield creds_dict  # Wait for tests to run


@pytest.fixture
def biological_data():
    """Create sample biological data for testing."""
    # Create sample data for different sheets
    catch_data = pd.DataFrame(
        {
            "haul": [1, 2, 3, 4],
            "ship_id": [160, 160, 584, 584],
            "survey": [201906, 201906, 2019097, 2019097],
            "species_code": [22500, 22500, 22500, 30420],
            "weight_in_haul": [150.5, 200.3, 175.8, 90.2],
        }
    )

    length_data = pd.DataFrame(
        {
            "haul": [1, 1, 2, 3, 4],
            "ship_id": [160, 160, 160, 584, 584],
            "survey": [201906, 201906, 201906, 2019097, 2019097],
            "species_code": [22500, 22500, 22500, 22500, 30420],
            "sex": [1, 2, 1, 3, 2],
            "length": [45.2, 48.7, 42.5, 51.3, 38.9],
            "frequency": [10, 15, 8, 12, 5],
        }
    )

    specimen_data = pd.DataFrame(
        {
            "haul": [1, 2, 3, 4],
            "ship_id": [160, 160, 584, 584],
            "survey": [201906, 201906, 2019097, 2019097],
            "species_code": [22500, 22500, 22500, 30420],
            "sex": [1, 1, 2, 2],
            "length": [45.2, 42.5, 51.3, 38.9],
            "age": [3, 2, 4, 1],
        }
    )

    return {"catch": catch_data, "length": length_data, "specimen": specimen_data}


@pytest.fixture
def bio_excel_file(biological_data, tmp_path):
    """Create a temporary Excel file with biological data sheets."""
    file_path = tmp_path / "test_biodata.xlsx"

    with pd.ExcelWriter(file_path) as writer:
        biological_data["catch"].to_excel(writer, sheet_name="biodata_catch", index=False)
        biological_data["length"].to_excel(writer, sheet_name="biodata_length", index=False)
        biological_data["specimen"].to_excel(writer, sheet_name="biodata_specimen", index=False)

    return file_path


@pytest.fixture
def bio_sheet_map():
    """Create sheet mapping for biological data."""
    return {"catch": "biodata_catch", "length": "biodata_length", "specimen": "biodata_specimen"}


@pytest.fixture
def bio_column_map():
    """Create column mapping for biological data."""
    return {"frequency": "length_count", "haul": "haul_num", "weight_in_haul": "haul_weight"}


@pytest.fixture
def subset_dict():
    """Create subset dictionary for filtering biological data."""
    return {
        "ships": {160: {"survey": 201906}, 584: {"survey": 2019097, "haul_offset": 200}},
        "species_code": [22500],
    }


@pytest.fixture
def pg_subset_dict():
    """Create subset dictionary for filtering biological data."""
    return {
        "ship": {101: {"survey": 2025}},
        "species_code": [22500],
    }


@pytest.fixture
def label_map():
    """Create label mapping dictionary for biological data."""
    return {"sex": {1: "male", 2: "female", 3: "unsexed"}}
