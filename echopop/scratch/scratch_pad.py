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

## Merge the exports
merged_exports = merge_exports(df_intervals, df_cells, df_layers)

# Resort and reindex
sort_echoview_export_df(merged_exports, inplace=TR)

# Sort/reindex each fileset
sort_echoview_export_df(df_intervals, inplace=True)
update_transect_spacing(df_intervals, 10., inplace=True)

out1 = group_merge(df_cells, [df_intervals, df_layers])
out2 = merge_exports(df_intervals, df_cells, df_layers)

out1.reset_index(drop=True).equals(out2.sort_values(["transect_num", "interval", "layer"]).reset_index(drop=True))


out1.transect_num
sort_echoview_export_df(out2).equals(out1)
out1.equals(sort_echoview_export_df(out2))
out2.sort_values(["transect_num", "interval", "layer"]).reset_index(drop=True).equals(out1)

out2.transect_num.unique()

out2.equals(out1)
out11 = out1.sort_values(["transect_num", "interval", "layer"]).reset_index(drop=True)
out22 = out2.sort_values(["transect_num", "interval", "layer"]).reset_index(drop=True)
out22.equals(out1)
# Store original datatypes using dictionary comprehensions
cells_dtypes = df_cells.dtypes.to_dict()
intervals_dtypes = df_intervals.dtypes.to_dict()
layers_dtypes = df_layers.dtypes.to_dict()

# Combine all dtypes, prioritizing integer types
all_dtypes = {**layers_dtypes, **intervals_dtypes, **cells_dtypes}

# Get the common column names
interval_cells_cols = list(set(df_intervals.columns).intersection(df_cells.columns))

interval_cells_df = df_intervals.merge(df_cells, on=interval_cells_cols, how="outer")

interval_cells_layers_cols = list(set(df_layers.columns).intersection(interval_cells_df.columns))
merged_df = interval_cells_df.merge(df_layers, on=interval_cells_layers_cols, how="outer")

# Drop NA's
merged_df.dropna(inplace=True)

# Separate integer and non-integer columns
integer_cols = {col: dtype for col, dtype in all_dtypes.items() 
                if col in merged_df.columns and pd.api.types.is_integer_dtype(dtype)}

non_integer_cols = {col: dtype for col, dtype in all_dtypes.items() 
                    if col in merged_df.columns and not pd.api.types.is_integer_dtype(dtype)}

# Convert non-integer columns
merged_df = merged_df.astype(non_integer_cols)

merged_df.astype(all_dtypes)
# Explicit merges with hard-coded keys
merged_df = (
    df_intervals
    .merge(df_cells, how="outer")
).drop_duplicates()


    .merge(df_layers)
    df_cells
    .merge(df_intervals, on=["transect_num", "interval"], how="outer")
    .merge(df_layers, how="outer")
)

# Separate integer and non-integer columns
integer_cols = {col: dtype for col, dtype in all_dtypes.items() 
                if col in merged_df.columns and pd.api.types.is_integer_dtype(dtype)}

non_integer_cols = {col: dtype for col, dtype in all_dtypes.items() 
                    if col in merged_df.columns and not pd.api.types.is_integer_dtype(dtype)}

# Convert non-integer columns
merged_df = merged_df.astype(non_integer_cols)

# Handle integer columns without loops
if integer_cols:
    # Get column names that should be integers
    int_col_names = list(integer_cols.keys())
    
    # Create mask of columns with no NaN values
    no_na_mask = merged_df[int_col_names].notna().all()
    
    # Filter to only columns without NaN values
    valid_int_cols = no_na_mask[no_na_mask].index.tolist()
    
    # Create conversion dictionary and convert in one operation
    if valid_int_cols:
        int_conversion_dict = {col: integer_cols[col] for col in valid_int_cols}
        merged_df.loc[:, valid_int_cols] = merged_df.loc[:, valid_int_cols].astype(int_conversion_dict)

# Restore datatypes without loops using dict comprehension and masking
integer_cols = {col: dtype for col, dtype in all_dtypes.items() 
                if col in merged_df.columns and pd.api.types.is_integer_dtype(dtype)}

non_integer_cols = {col: dtype for col, dtype in all_dtypes.items() 
                    if col in merged_df.columns and not pd.api.types.is_integer_dtype(dtype)}

# Convert non-integer columns directly
merged_df = merged_df.astype(non_integer_cols)

# Convert integer columns only if they don't contain NaN
for col, dtype in integer_cols.items():
    if merged_df[col].notna().all():
        merged_df[col] = merged_df[col].astype(dtype)

merged_df.dropna()

        
df_intervals.merge(df_cells, on=["transect_num", "interval", "process_id"], how="outer")

df_cells.merge(df_intervals, how="")

# Explicit merges with hard-coded keys
merged_df = (
    df_cells
    .merge(df_intervals, on=["transect_num", "interval", "process_id"], how="outer")
    .merge(df_layers, on=["transect_num", "interval", "process_id"], how="outer")
)

# Restore datatypes without loops using dict comprehension and masking
integer_cols = {col: dtype for col, dtype in all_dtypes.items() 
                if col in merged_df.columns and pd.api.types.is_integer_dtype(dtype)}

non_integer_cols = {col: dtype for col, dtype in all_dtypes.items() 
                    if col in merged_df.columns and not pd.api.types.is_integer_dtype(dtype)}


# Validate inputs
if not isinstance(dataframe, pd.DataFrame):
    raise ValueError("dataframe must be a pandas DataFrame")
if not isinstance(dataframes_to_add, list) or not all(
    isinstance(df, pd.DataFrame) for df in dataframes_to_add
):
    raise ValueError("dataframes_to_add must be a list of pandas DataFrames")

# Store original data types of columns in the parent dataframe
original_dtypes = {id(dataframe): dataframe.dtypes}
original_dtypes.update({id(df): df.dtypes for df in dataframes_to_add})

# Merge all dictionaries into a single dictionary with unique keys
unique_dtypes = {}
for dtypes_dict in original_dtypes.values():
    unique_dtypes.update(dtypes_dict.to_dict())

# Ensure inner_on and outer_on are lists
inner_on = inner_on if inner_on is not None else []
outer_on = outer_on if outer_on is not None else []

inner_on_lst = [inner_on] if isinstance(inner_on, str) else inner_on
outer_on_lst = [outer_on] if isinstance(outer_on, str) else outer_on

# If inner_on is None, find common columns across all dataframes_to_add
if not inner_on:
    common_columns_inner = set.intersection(*(set(df.columns) for df in dataframes_to_add))
    inner_on_lst = list(common_columns_inner)

# Merge dataframes within dataframes_to_add on inner_on
if inner_on_lst:
    merged_inner_frames = reduce(
        lambda left, right: pd.merge(left, right, on=inner_on_lst, how=how), dataframes_to_add
    )
else:
    merged_inner_frames = pd.DataFrame()

# If outer_on is None, find common columns between dataframe and merged_inner_frames
if not outer_on:
    common_columns_outer = set(dataframe.columns).intersection(set(merged_inner_frames.columns))
    outer_on_lst = list(common_columns_outer)

# Merge dataframe with merged_inner_frames
merged_frame = dataframe.merge(merged_inner_frames, on=outer_on_lst, how=how)

# Perform merge with drop_na option
if drop_na:
    merged_frame = merged_frame.dropna()

# Restore original dtypes of columns in merged_frame
for col in merged_frame.columns:
    dtype_to_restore = unique_dtypes.get(col)
    if dtype_to_restore is not None:
        if pd.api.types.is_integer_dtype(dtype_to_restore):
            # Check if column contains NaN or inf values
            if merged_frame[col].isnull().any() or not merged_frame[col].notna().all():
                continue  # Skip conversion if NaN or inf present
            merged_frame[col] = merged_frame[col].astype(dtype_to_restore)

return merged_frame

(
    df_cells
    .merge(df_intervals, on=["transect_num", "interval"], how="outer")
    .merge(df_layers, on=["transect_num", "interval"], how="outer")
)


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

# Explicitly merge the 3 dataframes
# -- do not need group_merge as a separate method
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
