import pandas as pd
from pathlib import Path


def check_column_names(df: pd.DataFrame, expected_names: set, path_for_df: Path) -> None:
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
        raise RuntimeError(f"The file '{str(path_for_df.absolute())}' does not contain the "
                           f"following columns: {expected_names - df_expected_columns}")


def check_existence_of_file(file_path: Path):
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
        raise FileNotFoundError(f"The file '{str(file_path.absolute())}' does not exist.")





