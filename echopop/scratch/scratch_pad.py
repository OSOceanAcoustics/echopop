from pathlib import Path
from typing import Dict, Generator, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import xarray as xr
import glob
import re

from echopop.kriging import Kriging
from echopop.nwfsc_feat import get_proportions, ingest_nasc, load_data
from echopop.ingest.common import read_csv_file
from echopop.core.echopop_columns import ECHOVIEW_TO_ECHOPOP
import functools
# ===========================================
# Organize NASC file
nasc_path = Path("C:/Users/Brandyn/Documents/GitHub/EchoPro_data/echopop_2019/raw_nasc")
type(nasc_path.glob("*.csv"))  # <class 'generator'>
ev_export_paths: dict = {
    "analysis": nasc_path.glob("*(analysis).csv"),  # Removed leading / only
    "cells": nasc_path.glob("*(cells).csv"),
    "intervals": nasc_path.glob("*(intervals).csv"),
    "layers": nasc_path.glob("*(layers).csv"),
}

# Get the transect numbers
transect_pattern: str = r"T(\d+)"
# transect_reference = ingest_nasc.map_transect_num(ev_export_paths, transect_pattern)

filename = transect_reference["intervals"][0][-1]
transect_number = transect_reference["intervals"][0][0]
df = read_csv_file(filename)
df.rename(columns=ECHOVIEW_TO_ECHOPOP, inplace=True)
df["transect_num"] = transect_number

def read_echoview_export(filename: Path, 
                         validator: Optional[Any] = None) -> pd.DataFrame:
    """
    Generic reader for Echoview export CSVs. Used for files like analysis, cells, layers, 
    intervals.

    Parameters
    ----------
    filename : pathlib.Path
        Full path to the NASC CSV file.
    validator : Any
        File-specific validator, if defined.

    Returns
    -------
    pd.DataFrame
        Cleaned and formatted data.
    """

    # Read the CSV file
    df = read_csv_file(filename)

    # Rename columns used by Echopop
    df.rename(columns=ECHOVIEW_TO_ECHOPOP, inplace=True)

    # TODO: Validation step would be here

    return df


def read_echoview_nasc(filename: Path,
                       transect_num: float,
                       validator: Optional[Any] = None) -> pd.DataFrame:
    """
    Generic reader for Echoview export CSVs. Used for files like analysis, cells, layers, 
    intervals.

    Parameters
    ----------
    filename : pathlib.Path
        Full path to the NASC CSV file.
    transect_num : float
        Transect number to use for filtering or labeling.

    Returns
    -------
    pd.DataFrame
        Cleaned and formatted data.
    """

    # Read in the defined CSV file
    nasc_df = read_echoview_export(filename, validator)
    
    # Add transect number
    nasc_df["transect_num"] = transect_num
    
    # Fix latitude and longitude
    # ---- Latitude
    impute_bad_coordinates(nasc_df, "latitude")
    # ---- Longitude
    impute_bad_coordinates(nasc_df, "longitude")

    # Return the cleaned DataFrame
    return nasc_df

[read_echoview_nasc(path, t_num) for t_num, path in transect_reference["intervals"]]
[read_echoview_nasc(f, tnum) for f, tnum in transect_reference["intervals"]]
transect_reference["intervals"]
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
