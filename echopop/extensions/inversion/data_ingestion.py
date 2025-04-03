import glob
from pathlib import Path
from typing import Any, Dict, List, Literal, Tuple, Union

import numpy as np
import pandas as pd
import yaml
from lmfit import Parameters
from pandera import DataFrameModel
from pydantic import BaseModel

from ...utils.load_nasc import get_transect_numbers
from ...utils.validate_df import IsobathData, KrigedMesh
from ...utils.validate_dict import Geospatial, KrigingParameters, VariogramBase, VariogramEmpirical

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


def validate_inversion_export_directories(
    data_root_dir: Union[str, Path], file_directory: Union[str, Path]
) -> Tuple[str, list]:

    # Get the data root directory
    root_directory = Path(data_root_dir)

    # Construct the directorypaths: Export files
    # ---- Drop prepended "/" or "\\" to avoid using an absolute path
    if file_directory.startswith("/") or file_directory.startswith("\\"):
        file_directory = file_directory[1:]
    # ---- Create directorypath
    file_folder = str(root_directory / file_directory)
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
    return file_folder, inversion_files


####################################################################################################
# Data ingestion
####################################################################################################
def read_echoview_inversion_file(
    data_file: Union[str, Path], transect_number: Union[int, float]
) -> pd.DataFrame:
    """
    Read and validate Echoview export files
    """

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
        "VL_end",
        "VL_start",
        "Dist_E",
        "Dist_S",
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

    # Compute distance
    export_file["distance"] = export_file["vl_end"] - export_file["vl_start"]

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
    if not export_file.empty:
        return export_file.dropna()


def load_inversion_configuration(inversion_config_path: Union[str, Path]):
    """
    Loads the inversion configuration settings.

    Parameters
    ----------
    inversion_config_path : Union[str, Path]
        A string or Path specifying the path to the inversion YAML file

    """

    # Validate the configuration file
    if not isinstance(inversion_config_path, Path):
        inversion_config_path = Path(inversion_config_path)

    # Validate configuration files
    config_existence = inversion_config_path.exists()

    # Error evaluation and print message (if applicable)
    if not config_existence:
        # ---- Raise Error
        raise FileNotFoundError(
            f"The inversion configuration file [{inversion_config_path.as_posix()}] does not exist!"
        ).with_traceback(None)

    # Read configuration files
    # ---- Initialization
    inversion_params = yaml.safe_load(inversion_config_path.read_text())

    # Convert center frequencies list to a NumPy array
    inversion_params["processing_parameters"]["acoustic_files"]["center_frequencies"] = np.array(
        inversion_params["processing_parameters"]["acoustic_files"]["center_frequencies"]
    )

    # Update the inversion parameter configuration
    inversion_params["processing_parameters"]["inversion"].update(
        {
            "index_columns": (
                ["transect_num"]
                if inversion_params["processing_parameters"]["acoustic_files"]["aggregate"]
                == "transect"
                else ["transect_num", "interval"]
            )
        }
    )

    # Update the Monte Carlo argument if it is not enabled
    if not inversion_params["processing_parameters"]["simulation"]["monte_carlo"]:
        inversion_params["processing_parameters"]["simulation"]["n_realizations"] = 1

    # Parse the bounds for all parameters and update the configuration
    parameter_limits = get_parameter_limits(inversion_params["scattering_parameters"])

    # Further update the parameterization
    inversion_params["processing_parameters"]["acoustic_models"].update(
        {
            "sv_threshold": (
                inversion_params["processing_parameters"]["acoustic_files"]["sv_threshold"]
            ),
            "parameter_bounds": parameter_limits,
            "scale_parameters": (
                inversion_params["processing_parameters"]["simulation"]["scale_parameters"]
            ),
        },
    )

    # Pass 'inversion_params' to the class instance
    return inversion_params


def get_parameter_limits(
    scattering_parameters: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Extract the lower and upper bounds for scattering model parameters
    """

    # Get the upper and lower values in case of rescaling
    parameter_limits = {
        key: {"low": value["low"], "high": value["high"]}
        for key, value in scattering_parameters.items()
    }

    return parameter_limits


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


def aggregate_interval(
    data_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Aggregate integration over intervals along each transect
    """

    # Drop empty cells
    data_reduced_df = data_df.copy()[data_df["sv_mean"] > -999.0]

    # Multiply `sv` by the cell height
    data_reduced_df["sv_h"] = data_reduced_df["sv"] * data_reduced_df["height_mean"]

    # Vertically integrate `sv_h`
    sv_interval = (
        data_reduced_df.groupby(["frequency", "transect_num", "interval"])["sv_h"]
        .sum()
        .to_frame(name="sv")
    )

    # Convert to the logarithmic domain
    sv_interval = 10 * np.log10(sv_interval)

    # Reset the index
    sv_interval.reset_index(inplace=True)

    # Pivot the DataFrame and return the dataset
    return sv_interval.pivot_table(
        index=["transect_num", "interval"], columns=["frequency"], values="sv"
    ).fillna(-999.0)


def aggregate_transect(
    data_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Aggregate integration over transect
    """

    # Drop empty cells
    data_reduced_df = data_df.copy()[data_df["sv_mean"] > -999.0]

    # Compute cell area (distance * thickness)
    data_reduced_df["area"] = data_reduced_df["distance"] * data_reduced_df["height_mean"]

    # Compute the total cell area per transect
    total_area = data_reduced_df.groupby(["frequency", "transect_num"])["area"].sum()

    # Compute the total NASC per transect
    total_nasc = data_reduced_df.groupby(["frequency", "transect_num"])["nasc"].sum()

    # Set the index
    data_reduced_df.set_index(["frequency", "transect_num"], inplace=True)

    # Compute the areal normalization weights
    data_reduced_df["area_weight"] = data_reduced_df["area"] / total_area.reindex(
        data_reduced_df.index
    )

    # Compute the area-weight
    data_reduced_df["sv_area"] = data_reduced_df["sv"] * data_reduced_df["area_weight"]

    # Compute the NASC-weight
    data_reduced_df["nasc_weight"] = data_reduced_df["nasc"] / total_nasc.reindex(
        data_reduced_df.index
    )

    # Compute the NASC-based weights for all coordinates
    # ---- Longitude
    data_reduced_df["longitude_nasc"] = (
        data_reduced_df["longitude"] * data_reduced_df["nasc_weight"]
    )
    # ---- Latitude
    data_reduced_df["latitude_nasc"] = data_reduced_df["latitude"] * data_reduced_df["nasc_weight"]

    # Reset the index
    data_reduced_df.reset_index(inplace=True)

    # Vertically integrate `sv_a`
    sv_transect = (
        data_reduced_df.groupby(["frequency", "transect_num"])["sv_area"].sum().to_frame(name="sv")
    )

    # Convert to the logarithmic domain
    sv_transect = 10 * np.log10(sv_transect)

    # Compute the NASC-weighted coordinates
    # ---- Longitude
    sv_transect["longitude"] = data_reduced_df.groupby(["frequency", "transect_num"])[
        "longitude_nasc"
    ].sum()
    # ---- Latitude
    sv_transect["latitude"] = data_reduced_df.groupby(["frequency", "transect_num"])[
        "latitude_nasc"
    ].sum()

    # Reset the index
    sv_transect.reset_index(inplace=True)

    # Pivot the DataFrame and return the dataset
    return sv_transect.pivot_table(
        index=["transect_num"], columns=["frequency"], values="sv"
    ).fillna(-999.0)


def integrate_measurements(
    measurement_data: pd.DataFrame,
    aggregate: Literal["interval", "transect"],
    sv_threshold: float,
) -> pd.DataFrame:
    """
    Integrate measurements based on transect and interval, or just transect
    """

    # Initialize copy
    measurement_data = measurement_data.copy()

    # Apply minimum threshold
    measurement_data.loc[measurement_data["sv_mean"] < sv_threshold, "sv_mean"] = -999.0

    # Linearize Sv to sv
    measurement_data["sv"] = 10.0 ** (measurement_data["sv_mean"] / 10.0)

    # Aggregation methods
    if aggregate == "interval":
        sv_indexed = aggregate_interval(measurement_data)
    elif aggregate == "transect":
        sv_indexed = aggregate_transect(measurement_data)
    else:
        raise ValueError(
            f"The defined aggregation method ['{aggregate}'] is invalid. This method must either be "
            f"'interval' or 'transect'."
        ).with_traceback(None)

    # Return the dataset
    return sv_indexed


def extract_metadata(
    measurement_data: pd.DataFrame, sv_threshold: float, **kwargs
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Extract metadata
    """

    # Sum up NASC per frequency for later reference in case of weighted apportionment
    nasc_table = measurement_data.pivot_table(
        index=["transect_num", "interval"], columns=["frequency"], values="nasc", aggfunc="sum"
    )

    # Sum up total layer height per interval
    metadata_layer_df = measurement_data.copy()
    # ---- Zero out empty cells
    metadata_layer_df.loc[metadata_layer_df["sv_mean"] < sv_threshold, "height_mean"] = 0.0
    # ---- Pivot table
    interval_height_table = metadata_layer_df.pivot_table(
        index=["transect_num", "interval", "layer"],
        columns=["frequency"],
        values="height_mean",
        aggfunc="sum",
    )

    # Drop duplicate values
    metadata_df = measurement_data.drop_duplicates(
        ["transect_num", "interval", "longitude", "latitude"]
    )

    # Set index
    metadata_df.set_index(["transect_num", "interval"], inplace=True)

    # Extract coordinates and return
    return metadata_df.filter(["longitude", "latitude"]), nasc_table, interval_height_table


def ingest_inversion_files(
    data_root_dir: Union[str, Path],
    file_directory: Union[str, Path],
    transect_pattern: str,
    aggregate: Literal["interval", "transect"],
    center_frequencies: np.ndarray[float],
    sv_threshold: float,
    **kwargs,
) -> Tuple[Dict[str, pd.DataFrame], pd.DataFrame]:
    """
    Ingest the inversion input files
    """

    # Validate relevant directories and file existence
    file_directory, inversion_files = validate_inversion_export_directories(
        data_root_dir, file_directory
    )

    # Get the transect numbers
    transect_reference = get_transect_numbers(inversion_files, transect_pattern, file_directory)

    # Read in the data
    input_df = dataset_reader(inversion_files, transect_reference)

    # Remove extraneous frequencies
    input_filtered_df = input_df[input_df["frequency"].isin(center_frequencies)]

    # Retrieve the metadata DataFrame
    metadata_df, nasc_table, interval_height_table = extract_metadata(
        input_filtered_df, sv_threshold
    )

    # Integrate Sv across the defined aggregation index
    input_integrated_df = integrate_measurements(input_filtered_df, aggregate, sv_threshold)

    # Organize the metadata DataFrames into a single dictionary
    metadata_dict = {
        "coordinates_df": metadata_df,
        "nasc_df": nasc_table,
        "layer_height_df": interval_height_table,
    }

    # Return the datasets
    return metadata_dict, input_integrated_df


def load_analysis_files(
    input_dict: Dict[str, Any],
    configuration_dict: Dict[str, Any],
):
    """
    Read files and additional analysis parameters
    """

    # Check for geostatistic files
    if "geostatistic_file_settings" in configuration_dict:
        # ---- File settings
        file_config = configuration_dict["geostatistic_file_settings"]
        # ---- Get root directory and format as Path
        data_root_directory = Path(file_config["data_root_dir"])
        # ---- Read and validate kriging mesh
        if "mesh" in file_config:
            # ---- Read in temporary DataFrame
            df_tmp = pd.read_excel(
                data_root_directory / file_config["mesh"]["filename"],
                sheet_name=file_config["mesh"]["sheetname"],
            )
            ## ---- Rename, if needed
            df_tmp.rename(
                columns={"centroid_longitude": "longitude", "centroid_latitude": "latitude"},
                inplace=True,
            )
            # ---- Validate
            input_dict.update({"mesh": KrigedMesh.validate_df(df_tmp)})
        # ---- Read and validate isobath reference coordinates
        if "isobath_reference" in file_config:
            # ---- Read in temporary DataFrame
            df_tmp = pd.read_excel(
                data_root_directory / file_config["isobath_reference"]["filename"],
                sheet_name=file_config["isobath_reference"]["sheetname"],
            )
            # ---- Validate
            input_dict.update({"isobath_reference": IsobathData.validate_df(df_tmp)})

    # Check for geospatial settings/configuration
    if "geospatial" in configuration_dict:
        # ---- Validate geospatial settings
        configuration_dict.update(
            {"geospatial": Geospatial.create(**configuration_dict["geospatial"])}
        )

    # Check for variogram parameterization/configuration
    if "variogram_parameters" in configuration_dict:
        # ---- Create temporary dictionary
        tmp_dict = {
            k: v["initial"] if isinstance(v, dict) and "initial" in v else v
            for k, v in configuration_dict["variogram_parameters"].items()
        }
        # ---- Validate the values
        tmp_validated_dict = VariogramBase.create(**tmp_dict)
        # ---- Validate the remaining keys
        tmp_validated_dict.update(**{**VariogramEmpirical.create(**tmp_dict)})
        # ---- If validators all succeeded, update the original values
        configuration_dict["variogram_parameters"].update(
            {
                k: (
                    {**v, "initial": tmp_validated_dict[k]}
                    if isinstance(v, dict)
                    else tmp_validated_dict.get(k, v)
                )
                for k, v in configuration_dict["variogram_parameters"].items()
            }
        )
