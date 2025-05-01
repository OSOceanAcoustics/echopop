import pandas as pd

# TODO: Add validation settings as a kwarg
def read_csv_file(filename: str) -> pd.DataFrame:
    """
    Read a *.csv file and convert column names to lowercase.
    
    Parameters
    ----------
    filename: str
        Path to the *.csvfile
        
    Returns
    -------
    pandas.DataFrame
        DataFrame with lowercase column names
    """

    # Read in the CSV file
    export_file = pd.read_csv(filename, index_col=None, header=0, skipinitialspace=True)

    # Set column names to lowercase
    export_file.columns = export_file.columns.str.lower()

    # Export the resulting `pandas.DataFrame`
    return export_file

# TODO: Add validation settings as a kwarg
def read_xlsx_file(filename: str, sheetname: str) -> pd.DataFrame:
    """
    Read an *.xlsx file and convert column names to lowercase.
    
    Parameters
    ----------
    filename: str
        Path to the *.xlsx file
    sheetname: str
        Name of the sheet to read
        
    Returns
    -------
    pandas.DataFrame
        DataFrame with lowercase column names
    """    
    pass