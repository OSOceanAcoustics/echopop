from pathlib import Path
from typing import Dict, Generator, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import xarray as xr
import glob
import re

from echopop.kriging import Kriging
from echopop.nwfsc_feat import get_proportions, ingest_nasc, load_data
from echopop.nwfsc_feat.ingest_nasc import read_echoview_nasc, map_transect_num, validate_transect_exports, clean_echoview_cells_df
from echopop.ingest.common import read_csv_file
from echopop.core.echopop_columns import ECHOVIEW_TO_ECHOPOP
import functools
# ===========================================
# Organize NASC file
nasc_path: Path = Path("C:/Users/Brandyn/Documents/GitHub/EchoPro_data/echopop_2019/raw_nasc")
filename_transect_pattern: str = r"T(\d+)"
default_transect_spacing = 10. # nmi
default_transect_spacing_latitude = 60. # deg N

ev_export_paths: dict = {
    "analysis": nasc_path.glob("*(analysis).csv"),  # Removed leading / only
    "cells": nasc_path.glob("*(cells).csv"),
    "intervals": nasc_path.glob("*(intervals).csv"),
    "layers": nasc_path.glob("*(layers).csv"),
}

transect_num_df = map_transect_num(ev_export_paths, filename_transect_pattern)
valid_transect_num_df = validate_transect_exports(transect_num_df)

df_intervals = pd.concat(
    echoview_nasc_to_df(valid_transect_num_df[valid_transect_num_df["file_type"] == "intervals"])
)
df_cells: pd.DataFrame = pd.concat(
    echoview_nasc_to_df(valid_transect_num_df[valid_transect_num_df["file_type"] == "cells"])
)
df_layers: pd.DataFrame = pd.concat(
    echoview_nasc_to_df(valid_transect_num_df[valid_transect_num_df["file_type"] == "layers"])
)

clean_echoview_cells_df(df_cells, inplace=True)

# Sort/reindex each fileset
sort_echoview_export_df(df_intervals, inplace=True)
sort_echoview_export_df(df_cells, inplace=True)
sort_echoview_export_df(df_layers, inplace=True)

def update_transect_spacing(
    transect_data: pd.DataFrame, 
    default_transect_spacing: float,
    latitude_threshold: float = 60.0,
    inplace: bool = False
) -> Optional[pd.DataFrame]:
    """
    Calculate and update the maximum spacing between transects.
    
    Parameters
    ----------
    transect_data : pd.DataFrame
        DataFrame containing transect data with 'transect_num', 'latitude' columns
    default_transect_spacing : float
        Default spacing to use (in nautical miles)
    latitude_threshold : float, default 60.0
        Maximum latitude to consider for spacing calculations
    inplace : bool, default False
        If True, modify the input DataFrame in-place and return None
    
    Returns
    -------
    pd.DataFrame or None
        Updated DataFrame if inplace=False, None otherwise
    """
    # Work on a copy or the original based on inplace parameter
    df = transect_data if inplace else transect_data.copy()
    
    # Get unique transect numbers
    transect_number = np.unique(df["transect_num"])
    
    # Initialize max transect spacing column
    df["transect_spacing"] = default_transect_spacing
    
    # Iterate through the transects to determine the maximum spacing
    for i in range(len(transect_number)):
        if i >= 2:
            # ---- For 2 transects prior to the current transect
            lag_2_index = df.index[
                (df["transect_num"] == transect_number[i - 2])
                & (df["latitude"] < latitude_threshold)
            ]
            # ---- For 1 transect prior to the current transect
            lag_1_index = df.index[
                (df["transect_num"] == transect_number[i - 1])
            ]
            # ---- Current transect
            current_index = df.index[
                (df["transect_num"] == transect_number[i])
                & (df["latitude"] < latitude_threshold)
            ]
            
            # Check if we have data for all three transects
            if len(lag_2_index) > 0 and len(lag_1_index) > 0 and len(current_index) > 0:
                # ---- Calculate the mean transect latitude (lag-2)
                lag_2_latitude = df.loc[lag_2_index, "latitude"].mean()
                # ---- Calculate the mean transect latitude (current)
                current_latitude = df.loc[current_index, "latitude"].mean()
                # ---- Compute the difference in the latitudes of adjacent transects
                delta_latitude = np.abs(current_latitude - lag_2_latitude)
                # ---- Get latitude range for the lag-2 transect
                latitude_2_range = (
                    df.loc[lag_2_index, "latitude"].max()
                    - df.loc[lag_2_index, "latitude"].min()
                )
                # ---- Get latitude range for current transect
                latitude_range = (
                    df.loc[current_index, "latitude"].max()
                    - df.loc[current_index, "latitude"].min()
                )
                # ---- Assign maximum spacing
                if (
                    (delta_latitude <= 2.0 * default_transect_spacing * 1.1 / 30.0)
                    & (latitude_2_range < 1 / 6)
                    & (latitude_range < 1 / 6)
                ):
                    df.loc[lag_1_index, "transect_spacing"] = delta_latitude * 30.0
    
    # Return the updated dataframe or None if inplace
    return None if inplace else df

out1 = update_transect_spacing(df_intervals, 10.)
out1[out1.transect_num == 15]
out1 = export_transect_spacing(df_intervals, default_transect_spacing=10.)
update_transect_spacing(df_intervals, 10., inplace=True)
df_intervals.loc[df_intervals["transect_num"] == 15]


transect_data = df_intervals.copy()
i = 14

# ---- For 2 transects prior to the current transect
lag_2_index = transect_data.index[
    (transect_data["transect_num"] == transect_number[i - 2])
    & (transect_data["latitude"] < latitude_threshold)
]
# ---- For 1 transect prior to the current transect
lag_1_index = transect_data.index[
    (transect_data["transect_num"] == transect_number[i - 1])
]
# ---- Current transect
current_index = transect_data.index[
    (transect_data["transect_num"] == transect_number[i])
    & (transect_data["latitude"] < latitude_threshold)
]
# ---- Calculate the mean transect latitude (lag-2)
lag_2_latitude = transect_data.loc[lag_2_index, "latitude"].mean()
# ---- Calculate the mean transect latitude (lag-2)
current_latitude = transect_data.loc[current_index, "latitude"].mean()
# ---- Compute the difference in the latitudes of adjacent transects
delta_latitude = np.abs(current_latitude - lag_2_latitude)
# ---- Get latitude range for the lag-2 transect
latitude_2_range = (
    transect_data.loc[lag_2_index, "latitude"].max()
    - transect_data.loc[lag_2_index, "latitude"].min()
)
# ---- Get latitude range for current transect
latitude_range = (
    transect_data.loc[current_index, "latitude"].max()
    - transect_data.loc[current_index, "latitude"].min()
)
# ---- Assign maximum spacing
if (
    (delta_latitude <= 2.0 * default_transect_spacing * 1.1 / 30.0)
    & (latitude_2_range < 1 / 6)
    & (latitude_range < 1 / 6)
):
    transect_data.loc[lag_1_index, "transect_spacing"] = delta_latitude * 30.0

transect_data[transect_data.transect_num == 15]
    
transect_data = out1.copy()
transect_number = np.unique(transect_data["transect_num"])
transect_data["transect_spacing"] = default_transect_spacing


out1.transect_spacing.max()
out1[out1.transect_num == 15]
out2.transect_spacing.max()

ev_export_paths: dict = {
    "analysis": nasc_path.glob("*(analysis).csv"),  # Removed leading / only
    "cells": nasc_path.glob("*(cells).csv"),
    "intervals": nasc_path.glob("*(intervals).csv"),
    "layers": nasc_path.glob("*(layers).csv"),
}

transect_num_df = map_transect_num(ev_export_paths, filename_transect_pattern)
valid_transect_num_df = validate_transect_exports(transect_num_df)

df_intervals = pd.concat(
    echoview_nasc_to_df(valid_transect_num_df[valid_transect_num_df["file_type"] == "intervals"])
)
df_cells: pd.DataFrame = pd.concat(
    echoview_nasc_to_df(valid_transect_num_df[valid_transect_num_df["file_type"] == "cells"])
)
df_layers: pd.DataFrame = pd.concat(
    echoview_nasc_to_df(valid_transect_num_df[valid_transect_num_df["file_type"] == "layers"])
)

clean_echoview_cells_df(df_cells, inplace=True)





# # ONLY do data ingestion and organization, not writing out anything
# def merge_echoview_nasc(
#     nasc_path: Path, nasc_filename_pattern: str, default_transect_spacing: float
# ) -> pd.DataFrame:
#     """
#     Ingest and merge all Echoview NASC files (intervals, cells, layers).

#     Parameters
#     ----------
#     nasc_path : Path
#         Directory containing Echoview export files (*.csv).
#     default_transect_spacing : float
#         Default spacing to impute where missing.

#     Returns
#     -------
#     pd.DataFrame
#         Merged dataframe from intervals, cells, and layers.
#     """
    
#     # Get all echoview NASC files: analysis, cells, intervals, layers
#     # -- do not need to validate, since non-existing folders/files will error out automatically
#     # -- use .glob to get all NASC-related files
#     ev_nasc_files: dict = {
#         "analysis": nasc_path.glob("*\.csv"),  # PATTERN built from nasc_filename_pattern
#         "cells": nasc_path.glob("PATTERN"),
#         "intervals": nasc_path.glob("PATTERN"),
#         "layers": nasc_path.glob("PATTERN"),
#     }

#     # Get all transect numbers
#     # -- for each transect number, there should be 4 files
#     # -- store only transect numbers with a complete set (4) csv files
#     # -- raise warning for transect numbers with an incomplete set (<4) csv files
#     transect_num: list = get_transect_num(ev_nasc_files)

#     # Read and concat intervals, cells, and layers dataframes
#     # -- use current code in consolidate_exports
#     # -- but do not worry about validator at this time
    # df_intervals: pd.DataFrame = pd.concat(
    #     echoview_nasc_to_df(ev_nasc_files["intervals"], transect_num), ...
    # )
    # df_cells: pd.DataFrame = pd.concat(
    #     echoview_nasc_to_df(ev_nasc_files["cells"], transect_num), ...
    # )
    # df_layers: pd.DataFrame = pd.concat(
    #     echoview_nasc_to_df(ev_nasc_files["layers"], transect_num), ...
    # )

#     # Wrangle cells_df columns

#     # Update transect spacing
#     # TODO: what does update_transect_spacing do?
#     df_intervals = update_transect_spacing(df_intervals, default_transect_spacing)

#     # Explicitly merge the 3 dataframes
#     # -- do not need group_merge as a separate method
#     df_merged: pd.DataFrame
#     return df_merged

# Get the transect numbers
transect_pattern: str = r"T(\d+)"
# transect_reference = ingest_nasc.map_transect_num(ev_export_paths, transect_pattern)

filename = transect_reference["intervals"][0][-1]
transect_number = transect_reference["intervals"][0][0]
df = read_csv_file(filename)
df.rename(columns=ECHOVIEW_TO_ECHOPOP, inplace=True)
df["transect_num"] = transect_number

# # Get all echoview NASC files: analysis, cells, intervals, layers
# # -- do not need to validate, since non-existing folders/files will error out automatically
# # -- use .glob to get all NASC-related files
# ev_export_paths: dict = {
#     "analysis": nasc_path.glob("(analysis).csv"),  # PATTERN built from nasc_filename_pattern
#     "cells": nasc_path.glob("(cells).csv"),
#     "intervals": nasc_path.glob("(intervals).csv"),
#     "layers": nasc_path.glob("(layers).csv"),
# }

# ev_export_paths: dict = {
#     "analysis": nasc_path.glob("*(analysis).csv"),  # Removed leading / only
#     "cells": nasc_path.glob("*(cells).csv"),
#     "intervals": nasc_path.glob("*(intervals).csv"),
#     "layers": nasc_path.glob("*(layers).csv"),
# }

# path_data = functools.reduce(
#     lambda acc, item: acc + [(
#         float(compiled_pattern.search(str(path)).group(1)), 
#         item[0], 
#         path
#     ) for path in list(item[1]) if compiled_pattern.search(str(path))],
#     ev_export_paths.items(),
#     []
# )


# {
#     t_num: dict([(ft, p) for tn, ft, p in path_data if tn == t_num])
#     for t_num in sorted(set(t[0] for t in path_data))
# }

# def map_transect_num(
#     ev_export_paths: Dict[str, Path], transect_pattern: str = r"T(\d+)" 
# ) -> Dict[float, Dict[str, Path]]:
#     """
#     Extract and map transect numbers to complete sets of Echoview export files.

#     Parameters
#     ----------
#     ev_export_paths : Dict[str, pathlib.Path]
#         Dictionary with keys "analysis", "cells", "intervals", "layers", with values being 
#         `pathlib.Path` glob generators pointing torward Echoview export *.csv files.
#     transect_pattern : str
#         Regex pattern to extract transect numbers from file names.

#     Returns
#     -------
#     Dict[float, Dict[str, pathlib.Path]]:
#         Dictionary with transect numbers as keys and dictionaries of file types ('analysis', 
#         'cells', 'intervals', and 'layers') and filepaths as values.
#     """

#     # Compite the transect pattern regex
#     compiled_pattern = re.compile(transect_pattern) 

#     # Reduce dictitonary into a list of tuples with transect number, file type, and path
#     path_mapping = functools.reduce(
#         lambda acc, item: acc + [(
#             float(compiled_pattern.search(str(path)).group(1)), 
#             item[0], 
#             path
#         ) for path in list(item[1]) if compiled_pattern.search(str(path))],
#         ev_export_paths.items(),
#         []
#     )

#     # Dictionary comprehension to create a mapping of transect numbers to file types and paths
#     transect_file_mapping = {
#         t_num: dict([(ft, p) for tn, ft, p in path_mapping if tn == t_num])
#         for t_num in sorted(set(t[0] for t in path_mapping))
#     }
    
#     # Return the transect_file_mapping
#     return transect_file_mapping

    

# # Get unique transect numbers
# transect_numbers = {t_num for t_num, _, _ in path_data}

# files = {
#     t_num: {
#         file_type: path
#         for t, file_type, path in path_data
#         if t == t_num
#     }
#     for t_num in transect_numbers
# }

# ##
# transect_number = 1

# # read_echoview_export
# filename = files[transect_number]["intervals"]

# # TODO: Add validation settings as a kwarg
# def read_csv_file(filename: str) -> pd.DataFrame:
#     """
#     Read a CSV file and convert column names to lowercase.
    
#     Parameters
#     ----------
#     filename: str
#         Path to the CSV file
        
#     Returns
#     -------
#     pandas.DataFrame
#         DataFrame with lowercase column names
#     """

#     # Read in the CSV file
#     export_file = pd.read_csv(filename, index_col=None, header=0, skipinitialspace=True)

#     # Set column names to lowercase
#     export_file.columns = export_file.columns.str.lower()

#     # Export the resulting `pandas.DataFrame`
#     return export_file


# # Read csv file
# export_file = pd.read_csv(filename, index_col=None, header=0, skipinitialspace=True)
# # ---- Set the column names to lowercase
# export_file.columns = export_file.columns.str.lower()

# # TODO: Add validation step
# # Validate
# # export_valid = validator.validate_df(export_file, filename)
# export_valid = export_file.copy()

# # Assign transect number column
# export_valid["transect_num"] = transect_number

# ##
# from echopop.core.echopop_columns import ECHOVIEW_TO_ECHOPOP

# # Assign the names
# export_valid.rename(columns=ECHOVIEW_TO_ECHOPOP, inplace=True)



# path_data = [
#     (int(match.group(1)), file_type, path)
#     for file_type, paths in file_dict.items()
#     for path in paths
#     if (match := compiled_pattern.search(str(path)))
# ]

# def example(data: dict) -> dict:
#     data["hello"] = "world"

#     return data

# def example(data: dict) -> dict:
#     return data | {"hello": "world"}

# example({})

# type_patterns = {
#     "analysis": "(analysis).csv",
#     "cells": "(cells).csv", 
#     "intervals": "(intervals).csv",
#     "layers": "(layers).csv"
# }
# {
#     key: [(p, int(match.group(1)) if (match := compiled_pattern.search(str(p))) else None)
#             for p in nasc_path.rglob(f"*{pattern}")]
#     for key, pattern in type_patterns.items()
# }

# pattern = r"T(\d+)"
# file_dict = ev_nasc_files

# def extract_transect_numbers(file_dict: Dict[str, Path], 
#                            pattern: str = r"T(\d+)") -> Dict[str, List[Tuple[Path, Optional[int]]]]:
#     """
#     Extract transect numbers from Path objects using regex pattern without loops.
    
#     Args:
#         file_dict: Dictionary with generators of Path objects
#         pattern: Regex pattern to extract transect number
    
#     Returns:
#         Dictionary with lists of (path, transect_number) tuples
#     """

#     files_as_lists = {
#         key: list(generator) 
#         for key, generator in file_dict.items()
#     }
    
#     compiled_pattern = re.compile(pattern)
#     {
#         key: [(path, int(match.group(1)) if (match := compiled_pattern.search(str(path))) else None)
#               for path in paths]
#         for key, paths in file_dict.items()
#     }
    
#     return {
#         key: [(path, int(match.group(1)) if (match := compiled_pattern.search(str(path))) else None)
#               for path in list(paths)]
#         for key, paths in file_dict.items()
#     }

# df_merged = ingest_nasc.merge_echoview_nasc(nasc_path, nasc_filename_pattern)

# # Optional: only use for years needing this as external resources
# df_transect_region_key = ingest_nasc.load_transect_region_key(region_class_filepath)

# # Use df.to_csv to save df_transect_region_key, in place of the specialized transect_region_key file
# # Keep read_transect_region_file and make sure its output is the same as construct_transect_region_key


# # Age-1+
# df_nasc_all_ages = ingest_nasc.consolidate_echoview_nasc(
#     df_merged,
#     region_names=["Age-1 Hake", "Age-1 Hake Mix", "Hake", "Hake Mix"],
#     survey_identifier=survey_identifier,
# )

# # Age-2+ (no age 1)
# df_nasc_no_age1 = ingest_nasc.consolidate_echoview_nasc(
#     df_merged,
#     region_names=["Hake", "Hake Mix"],
#     survey_identifier=survey_identifier,
# )

# # Use df.to_csv to save df_nasc_all_ages and df_nasc_no_age1 if needed

# # Use regular pd.read_csv to read df_nasc_*, effectively break up the current load_data()
# # -- there is no need to have a one-size-fits-all load_data function
# # -- just read them in without validation is fine: these are all files under our control


# # ===========================================
# # Execute what's in Survey.load_survey_data()
# # All *_dict below are a subdict from the original config yaml

# root_path = "WHERE_ALL_DATA_ARE"
# species_code = "SPECIES_CODE"
# df_nasc_no_age1: pd.DataFrame  # extracted nasc data from above, can also be df_nasc_all_ages

# bio_path_dict: dict  # the "biological" section of year_config.yml
# # this will be simplified now that we read from the master spreadsheet
# strata_path_dict: dict  # the "stratification" section of year_config.yml

# df_bio_dict = load_data.load_biological_data(root_path, bio_path_dict, species_code)
# df_strata_dict = load_data.load_stratification(root_path, strata_path_dict)

# # Consolidate all input data into df_acoustic_dict
# df_nasc_no_age1 = load_data.consolidate_all_data(
#     df_nasc=df_nasc_no_age1, df_bio_dict=df_bio_dict, df_strata_dict=df_strata_dict
# )


# # ===========================================
# # Compute biological composition based on stratum
# length_bins: np.array  # length bin specification
# df_length_weight, df_regression = get_proportions.length_weight_regression(
#     df_bio_dict["specimen"], length_bins
# )
# # df_regression seems unused afterwards -- good as a record?

# # Get counts ----------------
# df_aged_counts = get_proportions.fish_count(  # previously "aged_number_distribution"
#     df_specimen=df_bio_dict["specimen"],
#     df_length=df_bio_dict["length"],
#     aged=True,
#     sexed=True,
# )
# df_unaged_counts = get_proportions.fish_count(  # previously "unaged_number_distribution"
#     df_specimen=df_bio_dict["specimen"],
#     df_length=df_bio_dict["length"],
#     aged=False,
#     sexed=True,
# )
# # Previously there was also "aged_number_distribution_filtered"
# # but it is simply df_aged_counts with unsexed fish removed,
# # I think it is better to have that explicitly in the code,
# # so removed OUTSIDE of the get_fish_count function


# # Get number proportions ----------------
# # TODO: DISCUSS THIS!
# # TODO: what does _overall stand for?
# da_number_proportion = get_proportions.number_proportion()


# # Get weight proportions ----------------
# # aged fish - weight distribution
# da_sex_length_age: xr.DataArray = get_proportions.weight_distributions(
#     df_specimen=df_bio_dict["specimen"],
#     df_length=df_bio_dict["length"],
#     df_length_weight=df_length_weight,
#     aged=True,
# )

# # unaged fish - weight distribution
# da_sex_length: xr.DataArray = get_proportions.weight_distributions(
#     df_specimen=df_bio_dict["specimen"],
#     df_length=df_bio_dict["length"],
#     df_length_weight=df_length_weight,
#     aged=False,
# )

# # Get stratum averaged weight for all sex, male, female
# df_averaged_weight = get_proportions.stratum_averaged_weight()


# # Get weight proportions ----------------
# # TODO: DISCUSS THIS!
# # TODO: what does _overall stand for?
# da_weight_proportion = get_proportions.weight_proportion()


# # ===========================================
# # NASC to number density


# # ===========================================
# # Perform kriging using class Kriging
# # TODO:
# # put back FEAT-specific kriging files

# # Load kriging-related params
# kriging_const: dict  # from initalization_config.yaml:
# # A0, longitude_reference, longitude/latitude_offset
# kriging_path_dict: dict  # the "kriging" section of year_config.yml
# # combined with the "kriging" section of init_config.yml
# kriging_param_dict, variogram_param_dict = load_data.load_kriging_variogram_params(
#     root_path=root_path,
#     file_path_dict=kriging_path_dict,
#     kriging_const=kriging_const,
# )

# kriging = Kriging(
#     kriging_param_dict=kriging_param_dict,
#     variogram_param_dict=variogram_param_dict,
#     mesh_template="PATH_TO_MESH_TEMPLATE",
#     isobath_template="PATH_TO_ISOBATH_REFERENCE",
# )

# # Create kriging mesh including cropping based on transects
# # Created mesh is stored in kriging.df_mesh
# kriging.create_mesh()

# # Perform coordinate transformation based on isobath if needed
# # This adds columns x/y to kriging.df_mesh
# kriging.latlon_to_xy()

# # Perform kriging
# # This adds kriging result columns to df_in
# df_out = kriging.krige(df_in=df_nasc_no_age1, variables="biomass")
