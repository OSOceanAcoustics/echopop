import glob
import os
import re
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

import numpy as np
import pandas as pd
from pandera.api.base.model import MetaModel

from ..core import ECHOVIEW_TO_ECHOPOP_NAMES, NAME_CONFIG, REGION_EXPORT_MAP
from ..spatial.transect import export_transect_layers, export_transect_spacing
from .operations import compile_patterns, extract_parts_and_labels, group_merge
from .validate_df import ECHOVIEW_DF_MODEL, KSStrata


def read_echoview_export(filename: str, transect_number: int, validator: MetaModel):
    """
    Read Echoview export file

    Parameters
    ----------
    filename: str
        Export file path
    transect_number: int
        Transect number associated with the export file
    validator: MetaModel
        A ``pandera`` ``DataFrameModel``-class validator

    See Also
    ----------
    echopop.load.ingest_echoview_exports
        A more detailed description
    """

    # Read csv file
    export_file = pd.read_csv(filename, index_col=None, header=0, skipinitialspace=True)
    # ---- Set the column names to lowercase
    export_file.columns = export_file.columns.str.lower()

    # Validate
    export_valid = validator.validate_df(export_file, filename)

    # Assign transect number column
    export_valid["transect_num"] = transect_number

    # Assign the names
    export_valid.rename(columns=ECHOVIEW_TO_ECHOPOP_NAMES, inplace=True)

    # Get coordinate column names, if present
    # ---- Extract longitude column name
    lon_col = next((col for col in export_valid.columns if "lon" in col.lower()), None)
    # ---- Extract latitude column name
    lat_col = next((col for col in export_valid.columns if "lat" in col.lower()), None)

    # Automatic "bad data" removal/imputation for georeferenced (Latitude/Longitude) data
    # ---- Helper imputation function
    def impute_bad_data(data, column):
        # ---- Find indices where values are "999.0" or `NaN`
        invalid_coords = data.index[(data[column] == 999.0) | (data[column].isna())].to_numpy()
        # ---- Break from helper if there are no invalid coordinates
        if len(invalid_coords) == 0:
            return
        # ---- Evaluate cases where there are multiple steps between indices (i.e. multiple groups)
        delta_invalid_coords = np.where(np.diff(invalid_coords) > 1)[0]
        # ---- Digitize into discrete groups of "bad" data
        idx_groups = np.digitize(invalid_coords, invalid_coords[delta_invalid_coords], right=True)
        # ---- Iterate across the bad data groups
        for group in np.unique(idx_groups):
            # ---- Get the matching group index
            in_group_idx = np.where(idx_groups == group)[0]
            # ---- Parse the index
            in_group = invalid_coords[in_group_idx]
            # ---- Case 1: Single or consecutive bad coordinates at head of transect
            if any(in_group == 0):
                # ---- Rolling imputation
                for i in range(in_group.size - 1, -1, -1):
                    data.loc[i, column] = 2 * data.loc[i + 1, column] - data.loc[i + 2, column]
            # ---- Case 2: Single or consecutive bad coordinates at tail of transect
            elif any(in_group + 2 > len(data)):
                # ---- Get starting index
                start_index = in_group.min()
                # ---- Get ending index
                end_index = in_group.max()
                # ---- Rolling imputation
                for i in range(start_index, end_index + 1, 1):
                    data.loc[i, column] = 2 * data.loc[i - 1, column] - data.loc[i - 2, column]
            # ---- Case 3: Single or consecutive bad coordinates in the middle of transect
            else:
                # ---- Get starting index
                start_index = in_group.min() - 1
                # ---- Get ending index
                end_index = in_group.max() + 1
                # ---- Calculate the mean difference between "good/valid" coordinates
                step_mean = (data.loc[end_index, column] - data.loc[start_index, column]) / (
                    in_group.size + 1
                )
                # ---- Impute over "bad" points
                data.loc[in_group, column] = data.loc[start_index, column] + step_mean * np.arange(
                    1, in_group.size + 1, 1
                )

    # ---- Impute latitude
    if lat_col is not None:
        impute_bad_data(export_valid, lat_col)
    # ---- Impute longitude
    if lon_col is not None:
        impute_bad_data(export_valid, lon_col)

    # Return dataframe
    return export_valid


def get_transect_numbers(export_files: list, transect_pattern: str, file_directory: str):
    """
    Extract transect numbers from export file strings.
    """

    # Validate `transect_pattern`
    if not isinstance(transect_pattern, str):
        raise TypeError(
            f"Input value for `transect_pattern` ({transect_pattern}) must be a regular expression "
            f"or `str`."
        )

    # Helper function for creating a filename-transect association reference dictionary
    def extract_transect_number(file, file_pattern=transect_pattern):
        # ---- Get base filename
        filename = os.path.basename(file)
        # ---- Parse string
        parsed_transect = re.search(file_pattern, filename)
        # ---- Extract transect number, if present
        if parsed_transect:
            try:
                return float(parsed_transect.group(1)), file
            # ---- Return None if the transect number is not parseable
            except ValueError:
                return None, file
        else:
            return None, file

    # ---- Generate filename-transect association reference dictionary
    transect_reference = {file: num for num, file in map(extract_transect_number, export_files)}

    # Evaluate whether any files failed to be parsed
    unparsed_transects = [file for file, num in transect_reference.items() if num is None]
    # ---- If some are missing
    if unparsed_transects:
        # ---- Raise error
        raise ValueError(
            f"Transect numbers could not be parsed from the following filename(s) in "
            f"{file_directory}:\n {unparsed_transects}"
        )

    # Return output
    return transect_reference


def consolidate_exports(transect_reference: dict, default_transect_spacing: float):
    """
    Read in the intervals, layers, and cells datasets associated with Echoview exports.
    """

    # Helper generator function for producing a single dataframe from entire file directory
    def generate_dataframes(files, validator, transect_reference=transect_reference):
        # ---- Iterate through directory
        for filename in files:
            # ---- Get transect number
            transect_num = transect_reference.get(filename, None)
            # ---- Read in file and impute, if necessary plus validate columns and dtypes
            yield read_echoview_export(filename, transect_num, validator)

    # Read in and concatenate all dataframes
    # ---- Intervals
    # -------- Interval files
    interval_files = [file for file in transect_reference.keys() if "(intervals)" in file]
    # -------- Get validator
    validator = ECHOVIEW_DF_MODEL["intervals"]
    # -------- Read in files
    intervals_df = pd.concat(
        generate_dataframes(interval_files, validator), axis=0, ignore_index=True
    )
    # ---- Cells
    # -------- Cells files
    cells_files = [file for file in transect_reference.keys() if "(cells)" in file]
    # -------- Get validator
    validator = ECHOVIEW_DF_MODEL["cells"]
    # -------- Read in files
    cells_df = pd.concat(generate_dataframes(cells_files, validator), axis=0, ignore_index=True)
    # ---- Layer
    # -------- Layer files
    layer_files = [file for file in transect_reference.keys() if "(layers)" in file]
    # -------- Get validator
    validator = ECHOVIEW_DF_MODEL["layers"]
    # -------- Read in files
    layers_df = pd.concat(generate_dataframes(layer_files, validator), axis=0, ignore_index=True)

    # Consolidate all three datasets
    # ---- Adjust the strings to avoid odd formatting
    cells_df["region_class"] = (
        cells_df["region_class"]
        .replace('"', "", regex=True)
        .replace(r"^\s*(.*?)\s*$", r"\1", regex=True)
    )
    # ---- Adjust the strings to avoid odd formatting for region names as well
    cells_df["region_name"] = (
        cells_df["region_name"]
        .replace('"', "", regex=True)
        .replace(r"^\s*(.*?)\s*$", r"\1", regex=True)
    )
    # ---- Convert to lower case
    cells_df["region_class"] = cells_df["region_class"].str.lower()

    # Get the transect spacing
    updated_intervals = export_transect_spacing(intervals_df, default_transect_spacing)

    # Merge the datasets together
    return group_merge(cells_df, [updated_intervals, layers_df]), updated_intervals


def load_export_regions(region_files: dict, region_names: dict, root_directory: str):
    """
    Load region definitions.
    """

    # Create helper function
    def read_export_region(
        root_directory: Union[str, Path],
        filename: str,
        sheetname: str,
        group: str,
        region_names: Union[List[str], str],
    ):

        # Validate `region_name`
        if not isinstance(region_names, (list, str)):
            raise TypeError("Argument 'region_name' must be either a `str` or `list`.")
        elif not isinstance(region_names, list):
            region_names = [region_names]

        # Import the `REGION_EXPORT_MAP`
        region_settings = REGION_EXPORT_MAP[group]

        # Generate full filepath
        full_path = Path(root_directory) / filename

        # Validate existence
        if not full_path.exists():
            raise FileNotFoundError(f"Region definition file '{str(full_path)}' not found!")

        # Read file
        if filename.lower().endswith(".csv"):
            # ---- Initially read in the file
            new_region_df = pd.read_csv(full_path)
            # ---- Adjust column name cases
            new_region_df.columns = new_region_df.columns.str.lower()
            # ---- Adjust column names, if needed
            new_region_df.rename(
                columns={
                    "transect": "transect_num",
                    "region id": "region_id",
                    "assigned haul": "haul_num",
                },
                inplace=True,
            )
        else:
            new_region_df = pd.read_excel(full_path, sheet_name=sheetname)

        # Fix strings if required
        new_region_df = new_region_df.replace('"', "", regex=True).replace(
            r"^\s*(.*?)\s*$", r"\1", regex=True
        )

        # Update column names if needed
        # ---- Get new column names
        new_column_names = {col: info["name"] for col, info in region_settings.items()}
        # ---- Update
        new_region_df = new_region_df.rename(columns=new_column_names)

        # Validate column presence/structure
        # ---- Expected columns
        # expected_columns = list(new_column_names.values())
        expected_columns = ["region_id", "haul_num", "transect_num"]
        # ---- Missing columns
        missing_columns = [col for col in expected_columns if col not in new_region_df.columns]
        # ---- Raise Error
        if missing_columns:
            raise ValueError(f"Missing columns in DataFrame for group '{group}': {missing_columns}")

        # Validate column datatypes
        # ---- Drop NaN values
        new_region_df = new_region_df.dropna(subset=["transect_num", "haul_num", "region_id"])
        # ---- Assign datatypes
        # new_region_df = new_region_df.astype(
        #     {info["name"]: info["type"] for info in region_settings.values()
        #      if info["name"] in new_region_df.columns}
        # )

        # Add group ID key
        new_region_df["group"] = group

        # Convert the region class column to lowercase
        if "region_class" in new_region_df.columns:
            new_region_df["region_class"] = new_region_df["region_class"].str.lower()

        # Filter the region dataframe based on regular expression filter
        if "region_class" in new_region_df.columns:
            # ---- Create filter
            pattern = "|".join(
                [r"(?<!-)\b{}\b".format(re.escape(name.lower())) for name in region_names]
            )
            # ---- Apply filter and return
            return new_region_df[
                new_region_df["region_class"].str.contains(pattern, case=False, regex=True)
            ]
        else:
            return new_region_df

    # Return the vertically concatenated DataFrame
    # ---- If there are not multiple groups, mutate the region files dictionary to match the groups
    if set(["filename"]).issubset(region_files):
        # ---- Transform dictionary
        region_files = {group: region_files for group in region_names}
        # ---- Add sheetname entry (None) if csv
        region_files = {
            k: {**v, "sheetname": None} if "sheetname" not in v else v
            for k, v in region_files.items()
        }

    return pd.concat(
        [
            read_export_region(
                root_directory,
                region_files[key]["filename"],
                region_files[key]["sheetname"],
                key,
                region_names[key],
            )
            for key in region_files
        ],
        ignore_index=True,
    )


def get_haul_strata_key(
    configuration_dict: Dict[str, Any], root_dir: str, transect_region_key: pd.DataFrame
) -> pd.Series:
    """
    Get the INPFC stratum and haul mapping

    """

    # Create copy
    transect_regions = transect_region_key.copy()

    # Flatten the configuration settings
    flat_table = pd.json_normalize(configuration_dict)

    # Get the configuration specific to KS strata
    strata_config = flat_table.filter(regex=r"\bstrata\b")

    # Get the target file and sheet
    # ---- Get filename
    file = strata_config.filter(regex="filename").values.flatten()
    # ---- Get sheetname
    sheet = strata_config.filter(regex="sheetname").values.flatten()[0]
    # ---- Find the sheetname with INPFC in it
    sheetname = [s for s in sheet if "inpfc" in s.lower()][0]

    # Create filepath
    # ---- Create filepath
    inpfc_strata_path = Path(root_dir) / file[0]
    # ---- Validate
    if not inpfc_strata_path.exists():
        raise FileExistsError(
            f"The INPFC latitude-based stratification file '{inpfc_strata_path.as_posix()}' does "
            "not exist!"
        )

    # Read in data
    strata_key = pd.read_excel(inpfc_strata_path, sheet_name=sheetname)
    # ---- Make column names lowercase
    strata_key.columns = strata_key.columns.str.lower()
    # ---- Rename the columns
    strata_key.rename(columns=NAME_CONFIG, inplace=True)
    # ---- Validate
    strata_key_filtered = KSStrata.validate_df(strata_key).set_index(["haul_num"])
    # ---- Set index for `transect_regions`
    trans_regions_copy = transect_regions.copy().set_index(["haul_num"])
    # ---- Cut the transect-region key
    trans_regions_copy["stratum_inpfc"] = strata_key_filtered["stratum_num"]
    # ---- Reset index
    trans_regions_copy.reset_index(inplace=True)
    # ---- Set new index
    trans_regions_copy.set_index(["stratum_inpfc"], inplace=True)

    # Get the haul-strata-region maps
    strata_map = configuration_dict["transect_region_mapping"]["inpfc_strata_region"]

    # Pre-allocate the column
    trans_regions_copy["country"] = None

    # Iterate through
    for region in strata_map:
        # ---- Get strata
        region_strata = strata_map[region]
        # ---- Assign country codes
        trans_regions_copy.loc[region_strata, "country"] = region

    # Get columns present that can be used for indexing
    cols_idx = list(
        set(trans_regions_copy.columns).intersection(
            set(["haul_num", "transect_num", "region_id", "region_name", "region_class"])
        )
    )

    # Return the country column Series
    return trans_regions_copy.reset_index().set_index(cols_idx)["country"]


def filter_export_regions(
    transect_data: pd.DataFrame,
    region_filter: Union[str, list[str]],
    index_variable: Union[str, List[str]] = ["transect_num", "interval"],
    unique_region_id: str = "region_id",
    region_class_column: str = "region_class",
    impute_regions: bool = True,
):

    # Check argument datatypes
    # ---- Convert into validation dictionary
    validation_dict = dict(region_filter=region_filter, index_variable=index_variable)
    # ---- `region_filter` and `index_variable`
    for keys, vars in validation_dict.items():
        if not isinstance(vars, (list, str)):
            raise TypeError(f"Argument '{keys}' must be either a `str` or `list`.")
        elif not isinstance(vars, list):
            validation_dict[keys] = [vars]

    # Validate that `unique_region_id`, `region_class_column`, and `index_variable` all exist
    for col in [unique_region_id, region_class_column] + index_variable:
        if col not in transect_data.columns:
            raise ValueError(f"Expected column '{col}' is missing from transect data DataFrame.")

    # Pass back parameters
    region_filter, index_variable = (
        validation_dict.get("region_filter"),
        validation_dict.get("index_variable"),
    )

    # Format the region ID pattern (as a regular expression)
    region_pattern = rf"^(?:{'|'.join([re.escape(name.lower()) for name in region_filter])})"

    # Apply the filter to only include the regions-of-interest
    transect_data_regions = transect_data[
        transect_data[region_class_column].str.contains(region_pattern, case=False, regex=True)
    ]

    # Impute region IDs for cases where multiple regions overlap within the same interval
    if impute_regions:
        # ---- Sort values
        transect_data_regions = transect_data_regions.sort_values(
            ["transect_num", "interval", unique_region_id], ignore_index=True
        )
        # ---- Re-assign the region ID based on the first element
        transect_data_regions.loc[:, unique_region_id] = transect_data_regions.groupby(
            ["interval", "transect_num"]
        )[unique_region_id].transform("first")
    # Return the output
    return transect_data_regions


def ingest_echoview_exports(
    configuration_dict: dict,
    transect_pattern: str,
    index_variable: Union[str, List[str]],
    read_transect_region_file: bool,
    region_class_column: str,
    unique_region_id: str,
    verbose: bool,
    write_transect_region_file: bool,
):
    """
    Import a directory comprising Echoview exports via batch processing

    Parameters
    ----------
    configuration_dict: dict
        Dictionary containing file information, directories, and other values stored within the
        configuration YAML files.
    transect_pattern: str
        See :func:`echopop.survey.load_acoustic_data`.
    index_variable: Union[str, List[str]]
        See :func:`echopop.survey.load_acoustic_data`.
    read_transect_region_file: bool
        See :func:`echopop.survey.load_acoustic_data`.
    region_class_column: str
        See :func:`echopop.survey.load_acoustic_data`.
    unique_region_id: str
        See :func:`echopop.survey.load_acoustic_data`.
    verbose: bool
        See :func:`echopop.survey.load_acoustic_data`.
    write_transect_region_file: bool
        See :func:`echopop.survey.load_acoustic_data`.

    Returns
    ----------
    pd.DataFrame:
        A `pandas.DataFrame` that includes the following columns:
        - `transect_num` (int): Transect number.
        - `distance_s` (float): The starting interval distance for each transect interval.
        - `distance_e` (float): The ending interval distance for each transect interval.
        - `latitude` (float): Latitude (degrees).
        - `longitude` (float): Longitude (degrees).
        - `transect_spacing` (float): Distance between transects (nmi).
        - `nasc_no_age` (float): Age-2+ NASC.
        - `nasc_all_ages` (float): Age-1+ NASC.
        - `haul_num` (int): The trawl/haul number associated with each transect interval.
        - `stratum_num` (Optional[int]): Length-based (KS clustered) stratum numbers that are added
        when `include_stratum = True` and `strata = 'ks'`.
        - `stratum_inpfc` (Optional[int]): Latitude-based (INPFC) stratum numbers that are added
        when `include_stratum = True` and `strata = 'inpfc'`.

    See Also
    ----------
    :func:`echopop.survey.load_acoustic_data`
    """

    # Get NASC export settings
    export_settings = configuration_dict["nasc_exports"]

    # Validate relevant directories and file existence
    save_folder, file_directory, export_files = validate_export_directories(configuration_dict)

    # Get the transect numbers
    transect_reference = get_transect_numbers(export_files, transect_pattern, file_directory)

    # Read in and consolidate all export files (intervals, layers, cells)
    transect_data, interval_template = consolidate_exports(
        transect_reference, export_settings["max_transect_spacing"]
    )
    # ---- Rename a column
    interval_template = interval_template.rename(columns={"max_depth": "bottom_depth"})

    # Process and construct the transect-region map key
    # ---- Ingest a pre-defined transect-region key
    if read_transect_region_file:
        # ---- Transect-region mapping file metadata
        region_files = configuration_dict["export_regions"]
        # ---- Region names
        region_names = export_settings["regions"]
        # ---- Root directory
        root_dir = configuration_dict["data_root_dir"]
        # ---- Read in the file
        transect_region_key = load_export_regions(region_files, region_names, root_dir)
        # ---- Get the country-code information, if parameterized
        if "inpfc_strata_region" in configuration_dict["transect_region_mapping"]:
            # ---- Assign country
            country_codes = get_haul_strata_key(configuration_dict, root_dir, transect_region_key)
            # ---- Get columns present that can be used for indexing
            cols_idx = list(
                set(transect_region_key.columns).intersection(
                    set(["haul_num", "transect_num", "region_id", "region_name", "region_class"])
                )
            )
            # ---- Set index
            transect_region_key.set_index(cols_idx, inplace=True)
            # ---- Assign
            transect_region_key["country"] = country_codes
            # ---- Reset index
            transect_region_key.reset_index(inplace=True)
    # ---- Construct the transect-region-key
    else:
        transect_region_key = construct_transect_region_key(transect_data, configuration_dict)

    # Write the file, if configured
    if write_transect_region_file:
        # ---- Root directory
        root_dir = configuration_dict["data_root_dir"]
        # --- Write the files
        write_transect_region_key(configuration_dict, root_dir, transect_region_key, verbose)

    # Get region info
    # ---- Region names
    region_names = export_settings["regions"]

    # Vertical/areal integration (over groups)
    for key, values in region_names.items():
        # ---- Get valid region names
        region_filter = values
        # ---- Apply the region filter
        transect_regions = filter_export_regions(
            transect_data,
            index_variable=index_variable,
            region_filter=region_filter,
            unique_region_id=unique_region_id,
            region_class_column=region_class_column,
            impute_regions=True,
        )
        # ---- Calculate transect spatial metrics
        transect_layer_summary = export_transect_layers(transect_regions)
        # ---- Vertically sum backscatter
        nasc_intervals = (
            transect_regions.groupby(list(index_variable + [unique_region_id]))["nasc"]
            .sum()
            .reset_index()
        )
        # ---- Partition specific grouped regions
        grouped_region = transect_region_key[transect_region_key["group"] == key]
        # ---- Create interval copy
        interval_copy = (
            interval_template.copy().set_index(["interval", "transect_num"]).sort_index()
        )
        # ---- Merge stratum and haul information into the integrated NASC dataframe
        nasc_hauls = nasc_intervals.merge(grouped_region).set_index(["interval", "transect_num"])
        # ---- Append the haul numbers
        interval_copy.loc[:, nasc_hauls.columns] = nasc_hauls
        # ---- Reset the index
        interval_copy.reset_index(inplace=True)
        # ---- Combine with the transect layer summaries
        full_interval_strata_df = interval_copy.merge(transect_layer_summary, how="left")
        # ---- Sort
        full_interval_strata_df = full_interval_strata_df.sort_values(
            ["transect_num", "distance_s", "distance_e"]
        )
        # ---- Fill NaN region id column with 999
        full_interval_strata_df[unique_region_id] = (
            full_interval_strata_df[unique_region_id].fillna(999).astype(int)
        )
        # ---- Fill haul with 0's
        full_interval_strata_df["haul_num"] = (
            full_interval_strata_df["haul_num"].fillna(0).astype(int)
        )
        # ---- Fill float/continuous columns
        full_interval_strata_df[["nasc", "layer_mean_depth", "layer_height", "bottom_depth"]] = (
            full_interval_strata_df[["nasc", "layer_mean_depth", "layer_height", "bottom_depth"]]
            .fillna(0)
            .astype(float)
        )
        # ---- Drop unused columns
        output_nasc = full_interval_strata_df.filter(
            [e.lower() for e in export_settings["file_columns"]]
        )
        # ---- Add the region group key
        output_nasc["group"] = key
        # ---- Format the save filename
        save_filename = export_settings["save_file_template"]
        # ---- Parse unique regions
        unique_country = [i for i in grouped_region["country"].dropna().unique() if i != "None"]
        # ---- REGION replacement
        save_filename = save_filename.replace("{REGION}", "_".join(unique_country))
        # ---- YEAR replacement
        save_filename = save_filename.replace("{YEAR}", str(configuration_dict["survey_year"]))
        # --- GROUP replacement
        save_filename = save_filename.replace("{GROUP}", key)
        # --- Generate output filepath
        save_filepath = "/".join([save_folder.replace("\\", "/"), save_filename])
        # ---- Update the configuration
        configuration_dict["NASC"][key]["filename"] = save_filepath
        # ---- Update the sheetname
        configuration_dict["NASC"][key]["sheetname"] = export_settings["save_file_sheetname"]
        # ---- Save xlsx file
        output_nasc.to_excel(
            excel_writer=str(Path(save_folder) / save_filename),
            sheet_name=export_settings["save_file_sheetname"],
            index=False,
        )
        # ---- Print out message
        if verbose:
            print(
                f"Updated NASC export file for group '{key}' saved at "
                f"'{str(Path(save_folder) / save_filename)}'."
            )


def construct_transect_region_key(transect_data: pd.DataFrame, configuration_dict: dict):

    # Get pattern configuration for filtering region names
    pattern_config = configuration_dict["transect_region_mapping"]["parts"]

    # Get unique region names
    unique_regions = pd.DataFrame({"region_name": transect_data["region_name"].unique()})

    # Compile string label patterns
    compiled_patterns = compile_patterns(pattern_config)

    # Helper function for extracting
    def extract_parts(row):
        """Extract parts and labels from the region_name and return as a Series."""
        extracted_labels = extract_parts_and_labels(
            row["region_name"], compiled_patterns, pattern_config
        )
        # Initialize column dictionary
        column_data = {"region_name": row["region_name"]}
        # Dynamically create the column names based on the pattern configuration
        column_data.update(
            {
                f"{part_name.lower()}": extracted_labels.get(part_name, None)
                for part_name in pattern_config.keys()
            }
        )
        # Return the column dictionary as a series
        return pd.Series(column_data)

    # Extract the various components from the configuration string
    unique_regions_coded = unique_regions.apply(extract_parts, axis=1)
    # ---- Set the index
    unique_regions_coded.set_index("region_name", inplace=True)

    # Apply valid types
    valid_dtypes = {
        "region_class": str,
        "haul_num": int,
        "country": str,
    }
    # ---- Replace missing `haul_num` with 0
    unique_regions_coded["haul_num"] = unique_regions_coded["haul_num"].fillna(0)
    # ---- Apply conversion
    unique_regions_filtered = unique_regions_coded.apply(
        lambda col: col.astype(valid_dtypes.get(col.name, type(col.iloc[0])))
    )
    # ---- Apply haul number offset if so defined
    if "CAN_haul_offset" in configuration_dict:
        unique_regions_filtered.loc[unique_regions_filtered["country"] == "CAN", "haul_num"] = (
            unique_regions_filtered.loc[unique_regions_filtered["country"] == "CAN", "haul_num"]
            + configuration_dict["CAN_haul_offset"]
        )

    # Map the regions-hauls
    # ---- Index the data based on region name
    transect_data_filtered = transect_data.set_index("region_name")
    # ---- Map the values
    transect_data_filtered.loc[:, unique_regions_filtered.columns] = unique_regions_filtered
    # ---- Reset the index
    transect_data_filtered.reset_index(inplace=True)

    # Pull out NASC export region settings
    nasc_regions = configuration_dict["nasc_exports"]["regions"]

    # Initialize a list to contain the grouped transect-region keys
    full_transect_region_map = []
    # ---- Iterate through
    for region in nasc_regions.keys():
        # ---- Format region pattern
        region_pattern = (
            rf"^(?:{'|'.join([re.escape(name.lower()) for name in nasc_regions[region]])})"
        )
        # ---- Filter the dataframe
        transect_region = transect_data_filtered[
            transect_data_filtered["region_class"].str.contains(
                region_pattern, case=False, regex=True
            )
        ]
        # ---- Ensure that `region_class` is lowercase
        transect_region.loc[:, "region_class"] = transect_region.loc[:, "region_class"].str.lower()
        # ---- Find the unique keys
        unique_regions_map = (
            transect_region.groupby(["haul_num", "transect_num", "country", "region_id"])[
                "region_name"
            ]
            .first()
            .reset_index()
            .sort_values(["haul_num"])
        )
        # ---- Add `group` name
        unique_regions_map["group"] = region
        # ---- Append
        full_transect_region_map.append(unique_regions_map)

    # Return the values
    return pd.concat(full_transect_region_map).sort_values(["haul_num"])


def write_transect_region_key(
    configuration_dict: Dict[str, Any],
    root_dir: str,
    transect_region_key: pd.DataFrame,
    verbose: bool,
) -> None:

    # Get the configuration information
    meta = configuration_dict["transect_region_mapping"]

    # Get the file template
    file_template = meta["save_file_template"]
    # ---- Parse unique regions
    unique_country = transect_region_key["country"].dropna().unique()
    # ---- REGION replacement
    file_template = file_template.replace("{COUNTRY}", "_".join(unique_country.tolist()))
    # ---- YEAR replacement
    file_template = file_template.replace("{YEAR}", str(configuration_dict["survey_year"]))

    # Get the sheetname
    sheetname = meta["save_file_sheetname"]

    # Get save directory
    save_dir = meta["save_file_directory"]

    # Iterate through the groups
    for group, ids in configuration_dict["nasc_exports"]["regions"].items():
        # --- GROUP replacement
        file_template_group = file_template.replace("{GROUP}", group)
        # --- Convert to lowercase
        ids_proc = [i.lower() for i in ids]
        # ---- Filter out any erroneous groups
        transect_filtered = transect_region_key.loc[transect_region_key.region_class.isin(ids_proc)]
        # ---- Assign group
        transect_filtered.loc[:, "group"] = group
        # ---- Format the path
        file_path = Path(root_dir + save_dir) / file_template_group
        # ---- Write file
        transect_filtered.to_excel(excel_writer=file_path, sheet_name=sheetname, index=False)
        # ---- Print out message
        if verbose:
            print(
                f"Updated export region, transect, and haul key map for '{group}' saved at "
                f"'{file_path.as_posix()}'."
            )


def validate_export_directories(configuration_dict: dict) -> Tuple[str, str, list]:

    # Get the data root directory
    root_directory = Path(configuration_dict["data_root_dir"])

    # Get NASC export settings
    export_settings = configuration_dict["nasc_exports"]

    # Construct the directorypaths: Save files
    # ---- Save file directory
    save_file_directory = export_settings["nasc_export_directory"]
    # ---- Drop prepended "/" or "\\" to avoid using an absolute path
    if save_file_directory.startswith("/") or save_file_directory.startswith("\\"):
        save_file_directory = save_file_directory[1:]
    # ---- Create directorypath
    save_folder = str(root_directory / save_file_directory)
    # ---- Validate existence
    if not Path(save_folder).exists():
        raise FileNotFoundError(f"Save directory for NASC file ({save_folder}) does not exist.")

    # Construct the directorypaths: Export files
    # ---- Export file directory
    export_file_directory = export_settings["export_file_directory"]
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
    # ---- Get export files
    export_files = glob.glob(file_folder + "/*.csv")

    # Return
    return save_folder, file_folder, export_files
