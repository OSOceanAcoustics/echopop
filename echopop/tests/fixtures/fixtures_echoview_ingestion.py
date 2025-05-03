import pytest
from pathlib import Path
from unittest.mock import MagicMock
import pandas as pd
import numpy as np
import tempfile
import os

# ==================================================================================================
# Mock Paths
# ----------
@pytest.fixture
def mock_paths():
    """Create mock Path objects for testing."""
    paths = {
        "analysis": [
            Path("T1_survey_data_(analysis).csv"),
            Path("T2_survey_data_(analysis).csv"),
            Path("T3_survey_data_(analysis).csv"),
        ],
        "cells": [
            Path("T1_survey_data_(cells).csv"),
            Path("T2_survey_data_(cells).csv"),
            Path("T3_survey_data_(cells).csv"),
        ],
        "intervals": [
            Path("T1_survey_data_(intervals).csv"),
            Path("T2_survey_data_(intervals).csv"),
            Path("T3_survey_data_(intervals).csv"),
        ],
        "layers": [
            Path("T1_survey_data_(layers).csv"),
            Path("T2_survey_data_(layers).csv"),
            Path("T3_survey_data_(layers).csv"),
        ],
    }
    
    # Convert lists to mock glob iterators
    mock_dict = {}
    for key, path_list in paths.items():
        mock_generator = MagicMock()
        mock_generator.__iter__.return_value = iter(path_list)
        mock_dict[key] = mock_generator
    
    return mock_dict

# ==================================================================================================
# Coordinate imputation/repair
# ----------------------------
@pytest.fixture
def bad_coords_at_start():    
    return pd.DataFrame({
        'latitude': [999.0, 999.0, 45.2, 45.3, 45.4],
        'longitude': [-125.1, -125.2, -125.3, -125.4, -125.5]
    })

@pytest.fixture
def bad_coords_in_middle():    
    return pd.DataFrame({
        'latitude': [45.1, 45.2, 999.0, 999.0, 45.5],        
        'longitude': [-125.1, -125.2, -125.3, -125.4, -125.5]
    })

@pytest.fixture
def bad_coords_at_end():    
    return pd.DataFrame({
        'latitude': [45.1, 45.2, 45.3, 999.0, 999.0],
        'longitude': [-125.1, -125.2, -125.3, -125.4, -125.5]
    })

@pytest.fixture
def multiple_bad_coord_groups():
    return pd.DataFrame({
        'latitude': [999.0, 45.1, 45.2, 999.0, 999.0, 45.5, 999.0],
        'longitude': [-125.1, -125.2, -125.3, -125.4, -125.5, -125.6, -125.7]
    })

@pytest.fixture
def no_bad_coords():
    return pd.DataFrame({
        'latitude': [45.1, 45.2, 45.3, 45.4, 45.5],
        'longitude': [-125.1, -125.2, -125.3, -125.4, -125.5]
    })

# ==================================================================================================
# Reading Echoview data files
# ---------------------------
@pytest.fixture
def sample_echoview_data():
    """Create sample Echoview export data for testing."""
    # Create data with original Echoview column names
    return pd.DataFrame({
        'date_s': ['2023-01-01', '2023-01-01', '2023-01-01'],
        'exclude_below_line_depth_mean': [200.5, 198.2, 201.3],
        'lat_s': [45.1, 45.2, 45.3],
        'lon_s': [-125.1, -125.2, -125.3],
        'prc_nasc': [120.5, 135.8, 118.2],
        'time_s': ['10:15:30', '10:16:00', '10:16:30'],
        'vl_end': [10.5, 11.0, 11.5],
        'vl_start': [10.0, 10.5, 11.0],
        'other_column': ['a', 'b', 'c']  # Column not in mapping
    })

@pytest.fixture
def echoview_temp_csv(sample_echoview_data):
    """Create a temporary CSV file with sample Echoview data."""
    with tempfile.NamedTemporaryFile(suffix='.csv', delete=False, mode='w+') as f:
        sample_echoview_data.to_csv(f.name, index=False)
        temp_path = f.name
    
    yield Path(temp_path)
    
    # Clean up the temporary file
    os.unlink(temp_path)

@pytest.fixture
def empty_echoview_data():
    """Create empty dataframe with Echoview headers."""
    return pd.DataFrame(columns=[
        'date_s', 'exclude_below_line_depth_mean', 'lat_s', 'lon_s', 
        'prc_nasc', 'time_s', 'vl_end', 'vl_start'
    ])

@pytest.fixture
def missing_columns_data():
    """Create dataframe with some columns missing."""
    return pd.DataFrame({
        'date_s': ['2023-01-01', '2023-01-01'],
        # Missing exclude_below_line_depth_mean
        'lat_s': [45.1, 45.2],
        'lon_s': [-125.1, -125.2],
        # Missing prc_nasc
        'time_s': ['10:15:30', '10:16:00'],
        # Missing vl_end
        'vl_start': [10.0, 10.5],
    })

@pytest.fixture
def extreme_values_data():
    """Create dataframe with extreme/edge values."""
    return pd.DataFrame({
        'date_s': ['9999-12-31', '0001-01-01', '2023-01-01'],
        'exclude_below_line_depth_mean': [np.inf, -np.inf, np.nan],
        'lat_s': [90.0, -90.0, 0.0],  # Max, min, and zero latitude
        'lon_s': [180.0, -180.0, 0.0],  # Max, min, and zero longitude
        'prc_nasc': [1e10, -1e10, 0.0],  # Very large, very negative, zero
        'time_s': ['23:59:59', '00:00:00', '12:00:00'],
        'vl_end': [1e6, -1e6, 0.0],
        'vl_start': [1e6, -1e6, 0.0],
    })


@pytest.fixture
def duplicate_columns_data():
    """Create dataframe with duplicate column names (after reading)."""
    # This will simulate a CSV with duplicate columns after reading
    df = pd.DataFrame()
    df['date_s'] = ['2023-01-01', '2023-01-01']
    df['lat_s'] = [45.1, 45.2]
    df['lat_s.1'] = [45.3, 45.4]  # Duplicate column that pandas would create
    df['lon_s'] = [-125.1, -125.2]
    return df

@pytest.fixture
def bad_coords_middle():
    """Create dataframe with bad coordinates in middle."""
    return pd.DataFrame({
        'date_s': ['2023-01-01', '2023-01-01', '2023-01-01', '2023-01-01', '2023-01-01'],
        'lat_s': [45.1, 45.2, 999.0, 999.0, 45.5],
        'lon_s': [-125.1, -125.2, -125.3, -125.4, -125.5]
    })


@pytest.fixture
def bad_coords_end():
    """Create dataframe with bad coordinates at end."""
    return pd.DataFrame({
        'date_s': ['2023-01-01', '2023-01-01', '2023-01-01', '2023-01-01', '2023-01-01'],
        'lat_s': [45.1, 45.2, 45.3, 999.0, 999.0],
        'lon_s': [-125.1, -125.2, -125.3, -125.4, -125.5]
    })


@pytest.fixture
def bad_coords_mixed():
    """Create dataframe with mixed patterns of bad coordinates."""
    return pd.DataFrame({
        'date_s': ['2023-01-01', '2023-01-01', '2023-01-01', '2023-01-01', '2023-01-01', 
                   '2023-01-01', '2023-01-01'],
        'lat_s': [999.0, 45.1, 45.2, 999.0, 999.0, 45.5, 999.0],
        'lon_s': [-125.1, -125.2, -125.3, -125.4, -125.5, -125.6, -125.7]
    })

@pytest.fixture
def bad_coords_all():
    """Create dataframe with all bad coordinates."""
    return pd.DataFrame({
        'date_s': ['2023-01-01', '2023-01-01', '2023-01-01'],
        'lat_s': [999.0, 999.0, 999.0],
        'lon_s': [999.0, 999.0, 999.0],
        'prc_nasc': [100, 200, 300]
    })


@pytest.fixture
def bad_coords_start():
    """Create dataframe with bad coordinates at start."""
    return pd.DataFrame({
        'date_s': ['2023-01-01', '2023-01-01', '2023-01-01', '2023-01-01', '2023-01-01'],
        'lat_s': [999.0, 999.0, 45.2, 45.3, 45.4],
        'lon_s': [-125.1, -125.2, -125.3, -125.4, -125.5]
    })