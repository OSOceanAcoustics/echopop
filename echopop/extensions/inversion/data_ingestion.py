import glob
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

import pandas as pd
from lmfit import Parameters
from pandera import DataFrameModel
from pydantic import BaseModel

####################################################################################################
# Validation / preparation
####################################################################################################


class inversion_configuration_validator(BaseModel):
    """
    Pydantic model for validating configuration parameters
    """

    # RETURNS: Dict[str, Any]
    pass


class dataset_validator(DataFrameModel):
    """
    Pandera model for validating dataset values
    """

    # RETURNS: pd.DataFrame
    pass


def validate_inversion_export_directories(configuration_dict: dict) -> Tuple[str, str, list]:

    # Get the data root directory
    root_directory = Path(configuration_dict["data_root_dir"])

    # Get NASC export settings
    export_settings = configuration_dict["inversion_exports"]

    # Construct the directorypaths: Save files
    # ---- Save file directory
    save_file_directory = export_settings["inversion_save_directory"]
    # ---- Drop prepended "/" or "\\" to avoid using an absolute path
    if save_file_directory.startswith("/") or save_file_directory.startswith("\\"):
        save_file_directory = save_file_directory[1:]
    # ---- Create directorypath
    save_folder = str(root_directory / save_file_directory)
    # ---- Validate existence
    if not Path(save_folder).exists():
        raise FileNotFoundError(
            f"Save directory for inversion files ({save_folder}) does not exist."
        )

    # Construct the directorypaths: Export files
    # ---- Export file directory
    export_file_directory = export_settings["inversion_file_directory"]
    # ---- Drop prepended "/" or "\\" to avoid using an absolute path
    if export_file_directory.startswith("/") or export_file_directory.startswith("\\"):
        export_file_directory = export_file_directory[1:]
    # ---- Create directorypath
    file_folder = str(root_directory / export_file_directory)
    # ---- Validate existence
    if not Path(file_folder).exists():
        raise FileNotFoundError(f"The export file directory {{{file_folder}}} not found!")

    # Validate export files existence
    # ---- Check whether files exist at all
    if not any(Path(file_folder).iterdir()):
        raise FileNotFoundError(f"The export file directory {{{file_folder}}} contains no files!")
    # ---- Get raw inversion files
    inversion_files = glob.glob(file_folder + "/*/*", recursive=True)

    # Return
    return save_folder, file_folder, inversion_files


def prepare_scattering_model_inputs(scattering_config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Prepare scattering model parameter inputs
    """
    # == functions/set_para.m
    # == functions/inversion_para.m

    # PHASE 1) INGEST VALUES FROM CONFIGURATION FILE
    # PHASE 2) VALIDATE USING `inversion_configuration_validator`
    # PHASE 3) COMPUTE INTERMEDIATE VARIABLES (e.g. acoustic wavenumber, position matrix)
    # PHASE 4) PASS TO SCATTERER CLASS
    #          --> EXTERNAL TO THIS FUNCTION

    # RETURNS: Validated scattering model inputs
    pass


def prepare_dataset(dataset: pd.DataFrame) -> Dict[str, Any]:
    """
    Prepare dataset inputs
    """

    # PHASE 1) INGEST DATASET (*.xlsx)
    # PHASE 2) VALIDATE USING `dataset_validator`
    # PHASE 3) PARTITION DATASET BASED ON DIFFERENT ECHOMETRICS (e.g. mean Sv, median Sv)

    # RETURNS: Validated dataset DataFrame objects used for inversion
    pass


def prepare_inversion_settings(inversion_config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Prepare inversion configuration and parameterization
    """

    # PHASE 1) INGEST VALUES FROM CONFIGURATION FILE
    # PHASE 2) VALIDATE USING `inversion_configuration_validator`
    # PHASE 3) COMPUTE INTERMEDIATE VARIABLES (e.g. acoustic wavenumber, position matrix)

    # RETURNS: Validated inversion and optimization parameters
    pass


####################################################################################################
# Data ingestion
####################################################################################################
def read_echoview_inversion_file(
    data_file: Union[str, Path], transect_number: Union[int, float]
) -> pd.DataFrame:

    # Set of valid columns required for the inversion analysis
    valid_cols = [
        "Interval",
        "Layer",
        "Sv_mean",
        "NASC",
        "Height_mean",
        "Depth_mean",
        "Lon_M",
        "Lat_M",
        "Frequency",
    ]

    # Read in the correct file and subsetted columns
    export_file = pd.read_csv(
        data_file,
        index_col=None,
        header=0,
        skipinitialspace=True,
        usecols=lambda x: x in valid_cols,
    )

    # Change the case of the column names to lower case
    export_file.columns = map(str.lower, export_file.columns)

    # Rename certain columns
    export_file.rename(columns={"lat_m": "latitude", "lon_m": "longitude"}, inplace=True)

    # Append transect number
    export_file["transect_num"] = transect_number

    # Assign datatypes
    export_file = export_file.astype(
        {
            "transect_num": float,
            "interval": int,
            "layer": int,
            "sv_mean": float,
            "nasc": float,
            "height_mean": float,
            "depth_mean": float,
            "longitude": float,
            "latitude": float,
            "frequency": float,
        }
    )

    # Return the export file
    return export_file


def yaml_configuration_reader(
    config_file: Union[str, Path],
) -> Dict[str, Union[float, int, Parameters, pd.DataFrame, str]]:
    """
    Read and validate the input parameterization YAML configuration
    """
    # == functions/load_para_data.m
    # == functions/load_geo_phy_para.m
    # == functions/get_simu_para.m

    # PHASE 1) READ CONFIGURATION FILE

    # RETURNS: Raw Dict[str, Any]
    pass


def dataset_reader(
    inversion_files: Union[str, Path], transect_reference: Dict[str, Union[int, float]]
) -> pd.DataFrame:
    """
    Read aggregate acoustic backscatter measurements
    """

    # Helper generator function for producing a single dataframe from entire file directory
    def generate_dataframes(files, transect_reference=transect_reference):
        # ---- Iterate through directory
        for filename in files:
            # ---- Get transect number
            transect_num = transect_reference.get(filename, None)
            # ---- Read in file and impute, if necessary plus validate columns and dtypes
            yield read_echoview_inversion_file(filename, transect_num)

    # Read in the files
    export_df = pd.concat(
        generate_dataframes(inversion_files, transect_reference),
        axis=0,
        ignore_index=True,
    )

    # Return the files
    return export_df
