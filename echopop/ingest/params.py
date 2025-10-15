from pathlib import Path
from typing import Dict, Tuple, Union

import pandas as pd


def load_kriging_variogram_params(
    geostatistic_params_filepath: Union[str, Path],
    sheet_name: str,
    column_name_map: Dict[str, str] = {},
) -> Tuple[Dict, Dict]:
    """
    Load kriging and variogram parameters from Excel file.

    Parameters
    ----------
    geostatistic_params_filepath : str or Path
        Path to the Excel file containing kriging parameters
    sheet_name : str
        Name of the sheet to load
    column_name_map : dict, optional
        Dictionary mapping original column names to new column names

    Returns
    -------
    tuple
        Tuple containing (kriging_params_dict, variogram_params_dict)
    """

    if not geostatistic_params_filepath.exists():
        raise FileNotFoundError(
            f"Kriging parameter file not found: {geostatistic_params_filepath.as_posix()}"
        )

    # Read Excel file into memory and then transpose
    df_initial = pd.read_excel(
        geostatistic_params_filepath, sheet_name=sheet_name, index_col=None, header=None
    ).T

    # Take the values from the first row and redefine them as the column headers
    df_initial.columns = df_initial.iloc[0]
    df_initial = df_initial.drop(0)

    # Extract kriging parameters
    kriging_params = (
        df_initial.filter(regex="krig[.]")
        .rename(columns=lambda x: x.replace("krig.", ""))
        .rename(columns=column_name_map)
        .to_dict(orient="records")[0]
    )

    # Extract variogram parameters
    variogram_params = (
        df_initial.filter(regex="vario[.]")
        .rename(columns=lambda x: x.replace("vario.", ""))
        .rename(columns=column_name_map)
        .to_dict(orient="records")[0]
    )

    return kriging_params, variogram_params
