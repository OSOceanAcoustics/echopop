# from pathlib import Path
# import pandas as pd
# import numpy as np
# from typing import Optional, Union
# from pandera.api.base.model import MetaModel

# # filepath: echopop/utils/load_nasc.py


# def read_csv_file(filename: str) -> pd.DataFrame:
#     """
#     Read a CSV file and convert column names to lowercase.

#     Parameters
#     ----------
#     filename: str
#         Path to the CSV file

#     Returns
#     -------
#     pd.DataFrame
#         DataFrame with lowercase column names
#     """
#     export_file = pd.read_csv(filename, index_col=None, header=0, skipinitialspace=True)
#     export_file.columns = export_file.columns.str.lower()
#     return export_file

# def impute_bad_coordinates(data: pd.DataFrame, column: str) -> None:
#     """
#     Impute bad or missing coordinates in a DataFrame.

#     Parameters
#     ----------
#     data: pd.DataFrame
#         DataFrame containing coordinate data
#     column: str
#         Name of the column containing coordinates to impute

#     Returns
#     -------
#     None
#         The function modifies the DataFrame in-place
#     """
#     # Find indices where values are "999.0" or `NaN`
#     invalid_coords = data.index[(data[column] == 999.0) | (data[column].isna())].to_numpy()

#     # Break from helper if there are no invalid coordinates
#     if len(invalid_coords) == 0:
#         return

#     # Evaluate cases where there are multiple steps between indices (i.e. multiple groups)
#     delta_invalid_coords = np.where(np.diff(invalid_coords) > 1)[0]

#     # Digitize into discrete groups of "bad" data
#     idx_groups = np.digitize(invalid_coords, invalid_coords[delta_invalid_coords], right=True) \
#                  if len(delta_invalid_coords) > 0 else np.zeros(len(invalid_coords), dtype=int)

#     # Iterate across the bad data groups
#     for group in np.unique(idx_groups):
#         # Get the matching group index
#         in_group_idx = np.where(idx_groups == group)[0]
#         # Parse the index
#         in_group = invalid_coords[in_group_idx]

#         # Case 1: Single or consecutive bad coordinates at head of transect
#         if any(in_group == 0):
#             # Rolling imputation
#             for i in range(in_group.size - 1, -1, -1):
#                 data.loc[i, column] = 2 * data.loc[i + 1, column] - data.loc[i + 2, column]

#         # Case 2: Single or consecutive bad coordinates at tail of transect
#         elif any(in_group + 2 > len(data)):
#             # Get starting index
#             start_index = in_group.min()
#             # Get ending index
#             end_index = in_group.max()
#             # Rolling imputation
#             for i in range(start_index, end_index + 1, 1):
#                 data.loc[i, column] = 2 * data.loc[i - 1, column] - data.loc[i - 2, column]

#         # Case 3: Single or consecutive bad coordinates in the middle of transect
#         else:
#             # Get starting index
#             start_index = in_group.min() - 1
#             # Get ending index
#             end_index = in_group.max() + 1
#             # Calculate the mean difference between "good/valid" coordinates
#             step_mean = (data.loc[end_index, column] - data.loc[start_index, column]) / (
#                 in_group.size + 1
#             )
#             # Impute over "bad" points
#             data.loc[in_group, column] = data.loc[start_index, column] + step_mean * np.arange(
#                 1, in_group.size + 1, 1
#             )

# def read_echoview_export(filename: str, transect_number: int, validator: MetaModel,
#                          column_mapping: Optional[dict] = None) -> pd.DataFrame:
#     """
#     Read Echoview export file, validate it, and prepare it for further processing.

#     Parameters
#     ----------
#     filename: str
#         Export file path
#     transect_number: int
#         Transect number associated with the export file
#     validator: MetaModel
#         A ``pandera`` ``DataFrameModel``-class validator
#     column_mapping: Optional[dict]
#         Dictionary mapping Echoview column names to EchoPop column names

#     Returns
#     -------
#     pd.DataFrame
#         Validated and processed DataFrame

#     See Also
#     --------
#     echopop.load.ingest_echoview_exports
#         A more detailed description
#     """
#     # Read and prepare the CSV file
#     export_file = read_csv_file(filename)

#     # Validate
#     export_valid = validator.validate_df(export_file, filename)

#     # Assign transect number column
#     export_valid["transect_num"] = transect_number

#     # Assign the names if column_mapping is provided
#     if column_mapping:
#         export_valid.rename(columns=column_mapping, inplace=True)

#     # Get coordinate column names, if present
#     lon_col = next((col for col in export_valid.columns if "lon" in col.lower()), None)
#     lat_col = next((col for col in export_valid.columns if "lat" in col.lower()), None)

#     # Impute latitude
#     if lat_col is not None:
#         impute_bad_coordinates(export_valid, lat_col)

#     # Impute longitude
#     if lon_col is not None:
#         impute_bad_coordinates(export_valid, lon_col)

#     return export_valid
