import pytest
import pandas as pd
import os
import tempfile

# ==================================================================================================
# Mock *.csv files
# ----------------
@pytest.fixture
def sample_csv_content():
    """
    Provide sample CSV content with mixed-case column names for testing.
    
    Returns
    -------
        str: Sample CSV data with 3 columns and 2 rows
    """
    return "Column1,COLUMN2,CoLuMn3\n1,2,3\n4,5,6"

@pytest.fixture
def sample_csv_file(sample_csv_content):
    """
    Create a temporary CSV file for testing.
    
    Parameters
    ----------
        sample_csv_content: Content to write to the CSV file
        
    Yields
    ------
        str: Path to the temporary CSV file
        
    Notes
    -----
        File is automatically deleted after the test
    """
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(sample_csv_content)  # Now properly indented under 'with'
        filename = f.name
    yield filename
    os.unlink(filename)

# ==================================================================================================
# Mock *.xlsx files
# -----------------
@pytest.fixture
def sample_excel_file():
    """Create a temporary Excel file with multiple sheets for testing."""
    # Create a temporary file
    with tempfile.NamedTemporaryFile(suffix='.xlsx', delete=False) as temp_file:
        temp_filename = temp_file.name
    
    # Create DataFrames for multiple sheets
    df1 = pd.DataFrame({
        'Column1': [1, 2, 3],
        'COLUMN2': ['a', 'b', 'c'],
        'MixedCase_Column': [1.1, 2.2, 3.3]
    })
    
    df2 = pd.DataFrame({
        'SHEET2_COL1': [4, 5, 6],
        'sheet2_col2': ['d', 'e', 'f']
    })
    
    # Write to Excel file with multiple sheets
    with pd.ExcelWriter(temp_filename) as writer:
        df1.to_excel(writer, sheet_name='Sheet1', index=False)
        df2.to_excel(writer, sheet_name='Sheet2', index=False)
    
    yield temp_filename
    
    # Clean up the temporary file
    os.unlink(temp_filename)