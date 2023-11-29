from pathlib import Path

import pandas as pd


def check_column_names(
    df: pd.DataFrame, expected_names: set, path_for_df: Path
) -> None:
    """
    Ensures that the input DataFrame has the expected column names.

    Parameters
    ----------
    df: pd.DataFrame
        The DataFrame for which the check should be applied
    expected_names: set
        All column names that should be within ``df``
    path_for_df: Path
        The full path to the Excel file used to create ``df``

    Raises
    ------
    RuntimeError
        If ``df`` does not have the expected columns
    """

    # get those expected column names that df has
    df_expected_columns = set(df.columns).intersection(expected_names)

    if len(df_expected_columns) != len(expected_names):
        raise RuntimeError(
            f"The file '{str(path_for_df.absolute())}' does not contain the "
            f"following columns: {expected_names - df_expected_columns}"
        )


def check_existence_of_file(file_path: Path) -> None:
    """
    Makes sure that the input ``file_path`` exists.

    Parameters
    ----------
    file_path: Path
        The path we want to make sure exists

    Raises
    ------
    FileNotFoundError
        If the path does not exist
    """

    if not file_path.exists():
        raise FileNotFoundError(
            f"The file '{str(file_path.absolute())}' does not exist."
        )


def check_and_read(
        param_filename: str,
        param_sheetname: str,
        cols_types: dict,
        params: dict
) -> pd.DataFrame:
    """
    For input data table, check file path and column names, read into DataFrame,
    retain only the target columns, and set the target column data types.

    Parameters
    ----------
    param_filename : str
        Name of configuration parameter specifying the file path
    param_sheetname : str
        Name of configuration parameter specifying the Excel file sheet name
    cols_types : dict
        Dictionary specifying the column names and types for the target input data
    params: dict
        Survey parameter dictionary

    Returns
    -------
    pd.DataFrame
        Read and checked Dataframe
    """
    # check existence of the files
    file_path = params["data_root_dir"] / params[param_filename]
    check_existence_of_file(file_path)

    # read in and check Excel file
    df = pd.read_excel(file_path, sheet_name=params[param_sheetname])
    check_column_names(df=df, expected_names=set(cols_types.keys()), path_for_df=file_path)

    # obtaining those columns that are required
    df = df[list(cols_types.keys())].copy()

    # set data types of dataframe
    df = df.astype(cols_types)

    return df
