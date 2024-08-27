import glob
import os
import re
from pathlib import Path
from typing import List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from ..core import BIODATA_HAUL_MAP, CONFIG_MAP, ECHOVIEW_EXPORT_MAP, REGION_EXPORT_MAP
from ..spatial.transect import export_transect_layers, export_transect_spacing
from .operations import group_merge


def read_echoview_export(filename: str, transect_number: int):
    """
    Read Echoview export file

    Parameters
    ----------
    filename: str
        Export file path.
    transect_number: int
        Transect number associated with the export file.

    See Also
    ----------
    echopop.load.batch_read_echoview_exports
        A more detailed description.
    """

    # Read csv file
    export_file = pd.read_csv(filename, index_col=None, header=0, skipinitialspace=True)

    # Assign transect number column
    export_file["transect_num"] = transect_number

    # Get coordinate column names, if present
    # ---- Extract longitude column name
    lon_col = next((col for col in export_file.columns if "lon" in col.lower()), None)
    # ---- Extract latitude column name
    lat_col = next((col for col in export_file.columns if "lat" in col.lower()), None)

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
        impute_bad_data(export_file, lat_col)
    # ---- Impute longitude
    if lon_col is not None:
        impute_bad_data(export_file, lon_col)

    # Return dataframe
    return export_file


def validate_echoview_exports(export_dataframe: pd.DataFrame):
    """
    Validate column names and data types of the Echoview export dataframe.
    """

    # Get the export dataframe columns present in the `ECHOVIEW_EXPORT_MAP` API
    valid_columns = {
        old: new_info
        for old, new_info in ECHOVIEW_EXPORT_MAP.items()
        if old in export_dataframe.columns
    }

    # Drop unnecessary columns
    export_dataframe_filtered = export_dataframe.copy().filter(valid_columns.keys())

    # Validate the datatypes
    for old, new_info in valid_columns.items():
        # ---- Get expected datatype for the matched key name
        expected_dtype = new_info["type"]
        # ---- Transform the datatype, if needed
        export_dataframe_filtered[old] = export_dataframe_filtered[old].astype(expected_dtype)

    # Transform the column names
    # ---- Get the new names
    rename_dict = {old: new_info["name"] for old, new_info in valid_columns.items()}
    # ---- Rename the columns
    return export_dataframe_filtered.rename(columns=rename_dict)


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
                return int(parsed_transect.group(1)), file
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
    def generate_dataframes(files, transect_reference=transect_reference):
        # ---- Iterate through directory
        for filename in files:
            # ---- Get transect number
            transect_num = transect_reference.get(filename, None)
            # ---- Read in file and impute, if necessary plus validate columns and dtypes
            yield validate_echoview_exports(read_echoview_export(filename, transect_num))

    # Read in and concatenate all dataframes
    # ---- Intervals
    # -------- Interval files
    interval_files = [file for file in transect_reference.keys() if "(intervals)" in file]
    # -------- Read in files
    intervals_df = pd.concat(
        generate_dataframes(interval_files, transect_reference), axis=0, ignore_index=True
    )
    # ---- Cells
    # -------- Cells files
    cells_files = [file for file in transect_reference.keys() if "(cells)" in file]
    # -------- Read in files
    cells_df = pd.concat(
        generate_dataframes(cells_files, transect_reference), axis=0, ignore_index=True
    )
    # ---- Layer
    # -------- Layer files
    layer_files = [file for file in transect_reference.keys() if "(layers)" in file]
    # -------- Read in files
    layers_df = pd.concat(
        generate_dataframes(layer_files, transect_reference), axis=0, ignore_index=True
    )

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
        expected_columns = list(new_column_names.values())
        # ---- Missing columns
        missing_columns = [col for col in expected_columns if col not in new_region_df.columns]
        # ---- Raise Error
        if missing_columns:
            raise ValueError(f"Missing columns in DataFrame for group '{group}': {missing_columns}")

        # Validate column datatypes
        # ---- Drop NaN values
        new_region_df = new_region_df.dropna(subset=["transect_num", "haul_num", "region_class"])
        # ---- Assign datatypes
        new_region_df = new_region_df.astype(
            {info["name"]: info["type"] for info in region_settings.values()}
        )

        # Add group ID key
        new_region_df["group"] = group

        # Convert the region class column to lowercase
        new_region_df["region_class"] = new_region_df["region_class"].str.lower()

        # Filter the region dataframe based on regular expression filter
        # ---- Create filter
        pattern = "|".join(
            [r"(?<!-)\b{}\b".format(re.escape(name.lower())) for name in region_names]
        )
        # ---- Apply filter and return
        return new_region_df[
            new_region_df["region_class"].str.contains(pattern, case=False, regex=True)
        ]

    # Return the vertically concatenated DataFrame
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


def get_haul_transect_key(configuration_dict: dict, root_dir: str):
    """
    Get the haul-to-transect KS stratification key.
    """

    # Flatten the configuration settings
    flat_table = pd.json_normalize(configuration_dict)

    # Get the configuration specific to KS strata
    strata_config = flat_table.filter(regex=r"\bstrata\b")

    # Get the target file and sheet
    # ---- Get filename
    file = strata_config.filter(regex="filename").values.flatten()
    # ---- Get sheetname
    sheet = strata_config.filter(regex="sheetname").values.flatten()

    # Get configuration map settings specific to the KS strata
    # ---- Dictionary to DataFrame
    CONFIG_DF = pd.DataFrame(CONFIG_MAP)
    # ---- Find 'strata'
    strata_config = CONFIG_DF.loc["strata"].dropna()
    # ---- Get the first valid index
    strata_reference = strata_config[strata_config.first_valid_index()]

    # Create filepath
    # ---- Create filepath
    haul_key_path = Path(root_dir) / file[0]
    # ---- Validate
    if not haul_key_path.exists():
        raise FileExistsError(f"The haul-to-transect file {{{str(haul_key_path)}}} does not exist!")

    # Read in data
    haul_key = pd.read_excel(haul_key_path, sheet_name=sheet[0])
    # ---- Retain only the necessary columns
    haul_key_filtered = haul_key.filter(strata_reference.keys())
    # ---- Assign datatypes
    haul_key_filtered = haul_key_filtered.astype(strata_reference, errors="ignore")

    # Return output
    return haul_key_filtered


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
        # ---- Re-assign the region ID based on the first element
        transect_data_regions.loc[:, unique_region_id] = transect_data_regions.groupby(
            ["interval", "transect_num"]
        )[unique_region_id].transform("first")
    # Return the output
    return transect_data_regions


def batch_read_echoview_exports(
    configuration_dict: dict,
    transect_pattern: str = r"T(\d+)",
    index_variable: Union[str, List[str]] = ["transect_num", "interval"],
    unique_region_id: str = "region_id",
    region_class_column: str = "region_class",
    create_haul_transect_key: bool = True,
    create_transect_region_key: bool = True,
):
    """
    Import a directory comprising Echoview exports via batch processing

    Parameters
    ----------
    file_directory: Union[str, Path]
        Filepath/directory where acoustic transect exports are located. It is assumed that the files
        in this directory include file quartets for each export with the following labels:
        'cells', 'layers', 'analysis', and 'intervals'.
    transect_pattern: str
        A (raw) string that corresponds to the transect number embedded within the base name of the
        file path associated with each export file. Defaults to ``r'T(\\d+)'``. See a further
        description below for more details.
    configuration_dict: dict
        Dictionary containing file information, directories, and other values stored within the
        configuration YAML files.

    Returns
    ----------
    pd.DataFrame:
        A `pandas.DataFrame` that includes the following columns:
        - `transect_num` (int): Transect number.
        - `vessel_log_start` (float): The starting vessel log distance for each transect interval.
        - `vessel_log_end` (float): The ending vessel log distance for each transect interval.
        - `latitude` (float): Latitude (degrees).
        - `longitude` (float): Longitude (degrees).
        - `transect_spacing` (float): Distance between transects (nmi).
        - `NASC_no_age` (float): Age-2+ NASC.
        - `NASC_all_ages` (float): Age-1+ NASC.
        - `haul_num` (int): The trawl/haul number associated with each transect interval.
        - `stratum_num` (Optional[int]): Length-based (KS clustered) stratum numbers that are added
        when `include_stratum = True` and `strata = 'ks'`.
        - `stratum_inpfc` (Optional[int]): Latitude-based (INPFC) stratum numbers that are added
        when `include_stratum = True` and `strata = 'inpfc'`.

    Notes
    ----------
    The string pattern for `transect_pattern` requires a consistent filename format that can be
    readily parsed by `read_echoview_exports`. The default value for `transect_pattern`
    (``r'T(\\d+)'``) enables the function to parse all numbers that trail the letter "T" in the base
    filename. For example, an example file of "C:/Path/User/Data/random_12_V34-T56-A78_G9.csv" would
    yield a transect number of '56' since it trails the "T" and does not accidentally use other
    numbers in the string. However, a filename like
    "C:/Path/User/Data/randomT_12_T34-T56-T78_G9.csv" would detect multiple transect numbers ('34',
    '56', and '78') since numbers trail the letter "T" in three places. Therefore, it is important
    to ensure that embedded transect numbers are differentiable from other digits that may appear
    in the base filename.
    """

    # Get NASC export settings
    export_settings = configuration_dict["nasc_exports"]

    # Write haul-to-transect key, if the appropriate data inputs are available
    if create_transect_region_key:
        write_haul_to_transect_key(configuration_dict)

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

    # Write the transect-region-haul key xlsx file
    if create_haul_transect_key:
        write_transect_region_key(transect_data, configuration_dict)

    # Get region info
    # ---- Region filenames
    region_files = configuration_dict["export_regions"]
    # ---- Region names
    region_names = export_settings["regions"]
    # ---- Root directory
    root_dir = configuration_dict["data_root_dir"]
    # ---- Read in region definitions
    regions_df = load_export_regions(region_files, region_names, root_dir)

    # Read in the haul-to-transect KS strata key
    haul_transect = get_haul_transect_key(configuration_dict, root_dir)

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
            .to_frame("NASC")
            .reset_index()
        )
        # ---- Partition specific grouped regions
        grouped_region = regions_df[regions_df["group"] == key]
        # ---- Merge stratum and haul information into the integrated NASC dataframe
        nasc_haul_strata = group_merge(nasc_intervals, [grouped_region, haul_transect])
        #
        full_interval_df = interval_template.merge(nasc_haul_strata, how="left")
        #
        full_interval_strata_df = full_interval_df.merge(transect_layer_summary, how="left")
        # ---- Sort
        full_interval_strata_df = full_interval_strata_df.sort_values(
            ["transect_num", "vessel_log_start", "vessel_log_end"]
        )
        # ---- Fill NaN region id column with 999
        full_interval_strata_df[unique_region_id] = (
            full_interval_strata_df[unique_region_id].fillna(999).astype(int)
        )
        # ---- Fill haul with 0's
        full_interval_strata_df["haul_num"] = (
            full_interval_strata_df["haul_num"].fillna(0).astype(int)
        )
        # ---- Fill stratum with 1's
        # TODO: Add clarification that these don't make an impact or something to that effect
        # ! Or can this be avoided (e.g. just replace with '0' instead of '1'?)
        full_interval_strata_df["stratum_num"] = (
            full_interval_strata_df["stratum_num"].fillna(1).astype(int)
        )
        # ---- Fill float/continuous columns
        full_interval_strata_df[["NASC", "layer_mean_depth", "layer_height", "bottom_depth"]] = (
            full_interval_strata_df[["NASC", "layer_mean_depth", "layer_height", "bottom_depth"]]
            .fillna(0)
            .astype(float)
        )
        # ---- Drop unused columns
        output_nasc = full_interval_strata_df.filter(export_settings["file_columns"])
        # ---- Format the save filename
        save_filename = export_settings["save_file_template"]
        # ---- REGION replacement
        save_filename = save_filename.replace("{REGION}", "US_CAN")
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
        print(
            f"Updated NASC export file for group '{key}' saved at "
            f"'{str(Path(save_folder) / save_filename)}'."
        )


def compile_patterns(pattern_config: dict):
    """Compile patterns for each part in the configuration."""
    compiled_patterns = {}
    for part_name, part_patterns in pattern_config.items():
        compiled_patterns[part_name] = [
            re.compile(pattern, re.IGNORECASE) for pattern in [p["pattern"] for p in part_patterns]
        ]
    return compiled_patterns


def extract_parts_and_labels(region_name: str, compiled_patterns: re.Pattern, pattern_config: dict):
    """Extract corresponding labels from a region name based on compiled patterns."""
    # extracted_parts = {}
    labels = {}
    remaining_name = region_name

    for part_name, patterns in compiled_patterns.items():
        for pattern in patterns:
            match = pattern.search(remaining_name)
            if match:
                matched_value = match.group(0)
                label = next(
                    (
                        p["label"]
                        for p in pattern_config[part_name]
                        if p["pattern"] == pattern.pattern
                    ),
                    matched_value,
                )
                labels[part_name] = label if label != "None" else matched_value
                # BELOW FOR DEBUGGING
                # --------
                # extracted_parts[part_name] = matched_value
                # --------
                remaining_name = remaining_name.replace(matched_value, "", 1)
                break

    # return extracted_parts, labels
    return labels


def write_transect_region_key(transect_data: pd.DataFrame, configuration_dict: dict):

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
                f"{part_name.lower()}": extracted_labels.get(part_name, "")
                for part_name in pattern_config.keys()
            }
        )

        # Return the column dictionary as a series
        return pd.Series(column_data)

    # Extract the various components from the configuration string
    unique_regions_coded = unique_regions.apply(extract_parts, axis=1)
    # ---- Set the index
    unique_regions_coded.set_index("region_name", inplace=True)

    # Index the data based on region name
    transect_data_filtered = transect_data.set_index("region_name")

    # Assign the haul numbers and country IDs
    transect_data_filtered.loc[:, unique_regions_coded.columns] = unique_regions_coded

    # Pull out transect region mapping
    transect_region_settings = configuration_dict["transect_region_mapping"]
    # ---- Get the save directory
    save_dir = transect_region_settings["save_file_directory"]
    # ---- Get sheetname
    sheet = transect_region_settings["save_file_sheetname"]
    # ---- Get root directory
    root_dir = configuration_dict["data_root_dir"]

    # Pull out NASC export region settings
    nasc_regions = configuration_dict["nasc_exports"]["regions"]

    # Extract save file template and format
    # ---- Template
    template = transect_region_settings["save_file_template"]
    # ---- Get unique country codes
    country_codes = "_".join(np.unique(transect_data_filtered["country"]).tolist()).strip("_")
    # ---- Get survey year
    survey_year = configuration_dict["survey_year"]
    # ---- Swap out elements in template
    file_template = template.replace("{COUNTRY}", country_codes).replace("{YEAR}", str(survey_year))

    # Iterate through each expected NASC region groups
    for region in nasc_regions.keys():
        # ---- Format region pattern
        region_pattern = (
            rf"^(?:{'|'.join([re.escape(name.lower()) for name in nasc_regions[region]])})"
        )
        # ---- Swap out the {GROUP} component of the filename string
        filename = file_template.replace("{GROUP}", region)
        # ---- Filter the dataframe
        transect_region = transect_data_filtered[
            transect_data_filtered["region_class"].str.contains(
                region_pattern, case=False, regex=True
            )
        ]
        # ---- Filter out specific column names
        transect_out = (
            transect_region.reset_index()
            .filter(["transect_num", "region_id", "haul_num", "region_name", "region_class"])
            .sort_values(["transect_num", "haul_num"])
        )
        # ---- Reduce the duplicates
        transect_out = transect_out.drop_duplicates()
        # ---- Write xlsx file
        transect_out.to_excel(
            excel_writer=str(Path(root_dir + save_dir) / filename), sheet_name=sheet, index=False
        )
        # ---- Print out message
        print(
            f"Update export region, transect, and haul key map for '{region}' saved at "
            f"'{str(Path(root_dir + save_dir) / filename)}'."
        )


def write_haul_to_transect_key(configuration_dict: dict):

    # Get the haul-to-transect mapping settings
    haul_to_transect_settings = configuration_dict["haul_to_transect_mapping"]

    # Get root directory
    root_directory = configuration_dict["data_root_dir"]

    if "gear_data" in configuration_dict:
        # ---- Get filename template
        name_template = haul_to_transect_settings["save_file_template"]
        # ---- Define datatypes
        dtypes = {key: values["type"] for key, values in BIODATA_HAUL_MAP.items()}
        # ---- Define column name conversion
        col_names = {key: values["name"] for key, values in BIODATA_HAUL_MAP.items()}
        # ---- Get gear data dictionary
        gear_data = configuration_dict["gear_data"]
        # ---- Get survey year
        survey_year = configuration_dict["survey_year"]
        # ---- Swap out the {YEAR} component
        name_template = name_template.replace("{YEAR}", str(survey_year))
        # ---- Iterate through regions
        for region in gear_data.keys():
            # ---- Update {COUNTRY_CODE} component
            save_file = name_template.replace("{COUNTRY}", region + ".xlsx")
            # ---- Get directory settings
            dir_settings = haul_to_transect_settings["file_settings"][region]
            # ---- Get directory
            save_dir = dir_settings["directory"]
            # ---- Create filepath
            filepath = Path(root_directory) / gear_data[region]["filename"]
            # ---- Get sheetname
            sheet = gear_data[region]["sheetname"]
            # ---- Read file
            data = pd.read_excel(
                filepath, sheet_name=sheet, usecols=BIODATA_HAUL_MAP.keys()
            ).dropna(subset=BIODATA_HAUL_MAP.keys())
            # ---- Update datatypes
            data = data.astype(dtypes)
            # ---- Change column names
            data = data.rename(columns=col_names).reset_index(drop=True)

            # ---- Write the file
            data.to_excel(
                excel_writer=str(Path(root_directory + save_dir) / save_file),
                sheet_name=dir_settings["sheetname"],
                index=False,
            )
            # ---- Print out message
            print(
                f"Haul-to-transect mapping file for '{region}' saved at "
                f"'{Path(root_directory + save_dir) / save_file}'."
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
        raise FileNotFoundError(f"The export file directory {{{file_folder}}} contains no files!")
    # ---- Get export files
    export_files = glob.glob(file_folder + "/*")

    # Return
    return save_folder, file_folder, export_files
