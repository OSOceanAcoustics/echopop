import pandas as pd
import pytest
from pathlib import Path

from echopop.nwfsc_feat.load_data import (
    load_biological_data,
    preprocess_biological_data,
    apply_ship_survey_filters,
)

def test_load_biological_data_basic(bio_excel_file, bio_sheet_map):
    """Test basic loading of biological data without optional parameters."""
    # Pass an empty dict instead of None for column_name_map
    result = load_biological_data(bio_excel_file, bio_sheet_map, column_name_map={})
    
    assert isinstance(result, dict)
    assert set(result.keys()) == set(bio_sheet_map.keys())
    
    for df in result.values():
        assert isinstance(df, pd.DataFrame)
        assert not df.empty

def test_load_biological_data_with_column_map(bio_excel_file, bio_sheet_map, bio_column_map):
    """Test loading with column name mapping."""
    result = load_biological_data(bio_excel_file, bio_sheet_map, column_name_map=bio_column_map)
    
    if "length" in result:
        assert "length_count" in result["length"].columns
        assert "haul_num" in result["length"].columns
        assert "haul" not in result["length"].columns
    
    if "catch" in result:
        assert "haul_weight" in result["catch"].columns
        assert "haul_num" in result["catch"].columns
        assert "weight_in_haul" not in result["catch"].columns

def test_load_biological_data_with_subset(bio_excel_file, bio_sheet_map, subset_dict):
    """Test loading with subset filtering."""
    # Pass an empty dict for column_name_map
    result = load_biological_data(bio_excel_file, bio_sheet_map, column_name_map={}, subset_dict=subset_dict)
    
    for df in result.values():
        if "species_code" in df.columns:
            assert df["species_code"].unique() == [22500]
            
        if "ship_id" in df.columns:
            assert set(df["ship_id"].unique()).issubset({160, 584})

def test_load_biological_data_with_label_map(bio_excel_file, bio_sheet_map, label_map):
    """Test loading with label mapping."""
    result = load_biological_data(bio_excel_file, bio_sheet_map, column_name_map={}, biodata_label_map=label_map)
    
    for df in result.values():
        if "sex" in df.columns:
            assert set(df["sex"].unique()).issubset({"male", "female", "unsexed"})
            assert not any(pd.api.types.is_numeric_dtype(val) for val in df["sex"].unique())

def test_load_biological_data_file_not_found(bio_sheet_map):
    """Test error handling for non-existent file."""
    non_existent_file = Path("non_existent_file.xlsx")
    
    with pytest.raises(FileNotFoundError):
        load_biological_data(non_existent_file, bio_sheet_map)

# Biological preprocessing tests
def test_preprocess_biological_data(biological_data, label_map):
    """Test biological data preprocessing."""
    result = preprocess_biological_data(biological_data, label_map)
    
    assert set(result["length"]["sex"].unique()).issubset({"male", "female", "unsexed"})
    assert set(result["specimen"]["sex"].unique()).issubset({"male", "female"})
    assert all(result["catch"] == biological_data["catch"])

def test_preprocess_biological_data_no_label_map(biological_data):
    """Test preprocessing without label map."""
    result = preprocess_biological_data(biological_data, None)
    
    for key in biological_data:
        pd.testing.assert_frame_equal(result[key], biological_data[key])

def test_preprocess_biological_data_empty_dict():
    """Test preprocessing with empty dictionary."""
    result = preprocess_biological_data({}, {"sex": {1: "male"}})
    
    assert result == {}

# Ship/survey filter tests
def test_apply_ship_survey_filters(biological_data, subset_dict):
    """Test ship/survey filtering."""
    df = biological_data["length"]
    result = apply_ship_survey_filters(df, subset_dict)
    
    assert result["species_code"].unique() == [22500]
    assert set(result["ship_id"].unique()).issubset({160, 584})

def test_apply_ship_survey_filters_no_subset(biological_data):
    """Test with no subset dict."""
    df = biological_data["length"]
    result = apply_ship_survey_filters(df, None)
    
    assert result is not df  # Not the same object
    pd.testing.assert_frame_equal(result, df)  # But same content