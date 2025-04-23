import copy
from pathlib import Path
from typing import List, Literal, Optional, Union

import numpy as np
import pandas as pd
import yaml

from ..core import BIODATA_HAUL_MAP, DATA_STRUCTURE, LAYER_NAME_MAP, NAME_CONFIG
from .data_structure_utils import map_imported_datasets
from .validate_df import DATASET_DF_MODEL
from .validate_dict import CONFIG_DATA_MODEL, CONFIG_INIT_MODEL


def load_configuration(init_config_path: Path, survey_year_config_path: Path):
    """
    Loads the biological, NASC, and stratification
    data using parameters obtained from the configuration
    files.

    Parameters
    ----------
    init_config_path : pathlib.Path
        A string specifying the path to the initialization YAML file
    survey_year_config_path : pathlib.Path
        A string specifying the path to the survey year YAML file

    Notes
    -----
    This function parses the configuration files and incorporates them into
    the Survey class object. This initializes the `config` attribute that
    becomes available for future reference and functions.
    """
    # Validate configuration files
    # Retrieve the module directory to begin mapping the configuration file location
    # current_directory = os.path.dirname(os.path.abspath(__file__))

    # Build the full configuration file paths and verify they exist
    config_files = [init_config_path, survey_year_config_path]
    config_existence = [init_config_path.exists(), survey_year_config_path.exists()]

    # Error evaluation and print message (if applicable)
    if not all(config_existence):
        # ---- Get missing config filenames
        missing_config = [
            files for files, exists in zip(config_files, config_existence) if not exists
        ]
        # ---- Join the strings
        missing_str = ",\n   ".join(config.as_posix() for config in missing_config)
        # ---- Raise Error
        raise FileNotFoundError(
            f"The following configuration files do not exist:\n   {missing_str}."
        ).with_traceback(None)

    # Read configuration files
    # ---- Initialization
    init_config_params = yaml.safe_load(init_config_path.read_text())
    # -------- Validate
    valid_init_config_params = CONFIG_INIT_MODEL(
        init_config_path.as_posix(), **init_config_params
    ).model_dump(exclude_none=True)
    # ---- Survey year data
    survey_year_config_params = yaml.safe_load(survey_year_config_path.read_text())
    # -------- Validate
    valid_survey_year_config_params = CONFIG_DATA_MODEL(
        survey_year_config_path.as_posix(), **survey_year_config_params
    ).model_dump(exclude_none=True)

    # Validate that initialization and survey year configuration parameters do not intersect
    config_intersect = set(valid_init_config_params.keys()).intersection(
        set(valid_survey_year_config_params.keys())
    )

    # Error evaluation, if applicable
    if config_intersect:
        raise RuntimeError(
            f"The initialization and survey year configuration files comprise the following"
            f"intersecting variables: {', '.join(config_intersect)}"
        )

    # Format dictionary that will parameterize the `config` class attribute
    # Join the initialization and survey year parameters into a single dictionary
    config_to_add = {**valid_init_config_params, **valid_survey_year_config_params}

    # Amend length/age distribution locations within the configuration attribute
    config_to_add["biometrics"] = {
        "bio_hake_len_bin": valid_init_config_params["bio_hake_len_bin"],
        "bio_hake_age_bin": valid_init_config_params["bio_hake_age_bin"],
    }

    del config_to_add["bio_hake_len_bin"], config_to_add["bio_hake_age_bin"]

    # Pass 'full_params' to the class instance
    return config_to_add


def load_dataset(
    input_dict: dict, configuration_dict: dict, dataset_type: Optional[Union[str, List[str]]] = None
):
    """
    Loads the biological, NASC, and stratification
    data using parameters obtained from the configuration
    files. This will generate data attributes associated with the tags
    defined in both the configuration yml files and the reference CONFIG_MAP
    and LAYER_NAME_MAP dictionaries.

    Parameters
    ----------
    input_dict: dict
        A dictionary containing the loaded survey data.
    configuration_dict: dict
        Dictionary that contains all of the `Survey`-object configurations found within
        the `config` attribute.
    dataset_type: Union[str, List[str]], optional
        A string (or list of strings) corresponding to the named datasets defined in both
        'CONFIG_MAP' and the dataset definitions located in the file configuration .yaml.
    """

    # Re-initialize the input keys, if needed
    # ---- Get name of the proposed input dictionary keys
    input_keys = [
        (
            LAYER_NAME_MAP[key]["superlayer"][0]
            if LAYER_NAME_MAP[key]["superlayer"]
            else LAYER_NAME_MAP[key]["name"]
        )
        # for key in CONFIG_MAP.keys()
        for key in DATASET_DF_MODEL.keys()
        if key in list(dataset_type)
    ]
    # ---- Map the complete datasets
    imported_data = map_imported_datasets(input_dict)
    # ---- Re-initialize if data already loaded to avoid duplication issues
    if set(input_keys).issubset(imported_data):
        # ---- Reset the relevant keys
        input_dict.update({key: copy.deepcopy(DATA_STRUCTURE["input"][key]) for key in input_keys})

    # Check whether data files defined from the configuration file exists
    # ---- Generate flat JSON table comprising all configuration parameter names
    flat_configuration_table = pd.json_normalize(configuration_dict).filter(regex="filename")

    # Coerce `dataset_type` into List[str], if needed
    dataset_type = [dataset_type] if isinstance(dataset_type, str) else dataset_type

    # Get the subset table if specific `dataset_type` is defined
    if dataset_type:
        # ---- Get the outermost dictionary keys
        outer_keys = flat_configuration_table.columns.str.split(".").str[0]
        # ---- Get the associated column names
        matching_columns = flat_configuration_table.columns[outer_keys.isin(dataset_type)]
        # ---- Filter the columns
        flat_configuration_table = flat_configuration_table.filter(matching_columns)
    # ---- Default to `CONFIG_MAP` keys otherwise
    else:
        # dataset_type = list(CONFIG_MAP.keys())
        dataset_type = list(DATASET_DF_MODEL.keys())
    # ---- Parse the flattened configuration table to identify data file names and paths
    parsed_filenames = flat_configuration_table.values.flatten()
    # ---- Evaluate whether either file is missing
    data_existence = [
        (Path(configuration_dict["data_root_dir"]) / file).exists() for file in parsed_filenames
    ]

    # Assign the existence status to each configuration file for error evaluation
    # ---- Error evaluation and print message (if applicable)
    if not all(data_existence):
        missing_data = parsed_filenames[~np.array(data_existence)]
        raise FileNotFoundError(f"The following data files do not exist: {missing_data}")

    # Get the applicable `CONFIG_MAP` keys for the defined datasets
    expected_datasets = set(DATASET_DF_MODEL.keys()).intersection(dataset_type)

    # Define root data directory
    data_root_directory = Path(configuration_dict["data_root_dir"])

    # Data validation and import
    # ---- Iterate through known datasets and datalayers
    for dataset in list(expected_datasets):
        # ---- Define datalayers
        if dataset == "biological":
            if "filename" and "sheetname" in configuration_dict[dataset]:
                datalayers = {
                    k: {"filename": configuration_dict[dataset]["filename"], "sheetname": v}
                    for k, v in configuration_dict[dataset]["sheetname"].items()
                }
        else:
            datalayers = configuration_dict[dataset]
        # ---- Iterate through
        for datalayer, config_settings in datalayers.items():

            # Define validation settings from CONFIG_MAP
            validation_settings = DATASET_DF_MODEL[dataset][datalayer]

            # Create reference index of the dictionary path
            config_map = [dataset, datalayer]

            # Create the file and sheetnames for the associated datasets
            file_name = data_root_directory / config_settings["filename"]
            sheet_name = config_settings["sheetname"]

            # Ensure sheet_name is a list (to handle cases with multiple lists)
            sheet_name = [sheet_name] if isinstance(sheet_name, str) else sheet_name

            # Validate data for each sheet
            for sheets in sheet_name:
                read_validated_data(
                    input_dict,
                    configuration_dict,
                    file_name,
                    sheets,
                    config_map,
                    validation_settings,
                )

    # If all files could be successfully read in, update the `imported_data` list
    imported_data = map_imported_datasets(input_dict)

    # Update the data format of various inputs within `Survey`
    prepare_input_data(input_dict, configuration_dict, imported_data)


def read_validated_data(
    input_dict: dict,
    configuration_dict: dict,
    file_name: Path,
    sheet_name: str,
    config_map: list,
    validation_settings: dict,
):
    """
    Reads in data and validates the data type of each column/variable

    Parameters
    ----------
    input_dict: dict
        Dictionary represent the `input` attribute for the `Survey`-class object
    file_name: Path
        The file name without the prepended file path
    sheet_name: str
        The Excel sheet name containing the target data
    config_map: list
        A list parsed from the file name that indicates how data attributes
        within `self` are organized
    validation_settings: dict
        The subset CONFIG_MAP settings that contain the target column names
    """

    # Based on the configuration settings, read the Excel files into memory. A format
    # exception is made for 'kriging.vario_krig_para' since it requires additional
    # data wrangling (i.e. transposing) to resemble the same dataframe format applied
    # to all other data attributes.
    if "vario_krig_para" in config_map:
        # Read Excel file into memory and then transpose
        df_initial = pd.read_excel(file_name, header=None).T

        # Take the values from the first row and redfine them as the column headers
        df_initial.columns = df_initial.iloc[0]
        df_initial = df_initial.drop(0)

        # Slice only the columns that are relevant to the echopop module functionality
        # df_filtered = df_initial.filter(validation_settings)
        df = validation_settings.validate_df(df_initial)
    else:
        # Read Excel file into memory -- this only reads in the required columns
        df_initial = pd.read_excel(file_name, sheet_name=sheet_name)
        # ---- Force the column names to be lower case
        df_initial.columns = df_initial.columns.str.lower()
        # ---- Rename the columns
        df_initial.rename(columns=NAME_CONFIG, inplace=True)
        # ---- Filter the dataset accordingly, if required
        if "ship_id" in configuration_dict and config_map[0] == "biological":
            # ---- Collect ship configuration
            ship_config = configuration_dict["ship_id"]
            # ---- Get ship IDs
            ship_ids = [*ship_config.keys()]
            # ---- Apply ship-based filter
            if "ship_id" in df_initial.columns:
                df_filtered = df_initial.loc[df_initial["ship_id"].isin(ship_ids)]
            # ---- Collect survey IDs, if present
            survey_ids = [
                v["survey"] for v in configuration_dict["ship_id"].values() if "survey" in v
            ]
            # ---- Apply survey-based filter
            if survey_ids and "survey" in df_initial.columns:
                df_filtered = df_filtered.loc[df_filtered["survey"].isin(survey_ids)]
            # ---- Collect haul offsets, if any are present
            haul_offsets = {
                k: v["haul_offset"]
                for k, v in configuration_dict["ship_id"].items()
                if "haul_offset" in v
            }
            # ---- Apply haul number offsets, if defined
            if haul_offsets:
                df_filtered["haul_num"] = df_filtered["haul_num"] + df_filtered["ship_id"].map(
                    haul_offsets
                ).fillna(0)
        else:
            df_filtered = df_initial.copy()
        # ---- Validate the dataframes
        df = validation_settings.validate_df(df_filtered)

    # Assign the data to their correct data attributes/keys
    if LAYER_NAME_MAP[config_map[0]]["superlayer"] == []:
        sub_attribute = LAYER_NAME_MAP[config_map[0]]["name"]
    else:
        sub_attribute = LAYER_NAME_MAP[config_map[0]]["superlayer"][0]

    # Step 2: Determine whether the dataframe already exists
    if sub_attribute in ["biology", "statistics", "spatial"]:
        # ---- A single dataframe per entry is expected, so no other fancy operations are needed
        if sheet_name.lower() == "inpfc":
            df_list = [input_dict[sub_attribute]["inpfc_strata_df"], df]
            input_dict[sub_attribute]["inpfc_strata_df"] = pd.concat(df_list)
        else:
            if config_map[0] == "kriging" and config_map[1] == "vario_krig_para":
                df_list = [input_dict[sub_attribute]["kriging"][config_map[1] + "_df"], df]
                input_dict[sub_attribute]["kriging"][config_map[1] + "_df"] = pd.concat(
                    df_list
                ).tail(1)
            elif config_map[0] == "kriging":
                df_list = [input_dict[sub_attribute]["kriging"][config_map[1] + "_df"], df]
                input_dict[sub_attribute]["kriging"][config_map[1] + "_df"] = pd.concat(df_list)
            else:
                df_list = [input_dict[sub_attribute][config_map[1] + "_df"], df]
                input_dict[sub_attribute][config_map[1] + "_df"] = pd.concat(df_list)
    # TODO: This can be refactored out
    elif sub_attribute == "acoustics":
        # ---- Toggle through including and excluding age-1
        if config_map[1] == "no_age1":
            df = df.rename(
                columns={
                    "nasc": "NASC_no_age1",
                    "haul_num": "haul_no_age1",
                    "stratum_num": "stratum_no_age1",
                }
            )
        else:
            df = df.rename(
                columns={
                    "nasc": "NASC_all_ages",
                    "haul_num": "haul_all_ages",
                    "stratum_num": "stratum_all_ages",
                }
            )

        if input_dict["acoustics"]["nasc_df"].shape == (0, 0):
            input_dict["acoustics"]["nasc_df"] = df
        else:
            column_to_add = df.columns.difference(
                input_dict["acoustics"]["nasc_df"].columns
            ).tolist()
            input_dict["acoustics"]["nasc_df"][column_to_add] = df[column_to_add]
    else:
        raise ValueError(
            "Unexpected data attribute structure. Check the settings in "
            "the configuration YAML and core.py."
        )


def write_haul_to_transect_key(configuration_dict: dict, verbose: bool):
    """
    Function for writing the haul-transect mapping key .xlsx file.

    Parameters
    ----------
    configuration_dict: dict
        Dictionary containing file information, directories, and other values stored within the
        configuration YAML files.
    verbose: bool
        Console messages that will print various messages, updates, etc. when set to True.
    """

    # Get the haul-to-transect mapping settings
    haul_to_transect_settings = configuration_dict["haul_to_transect_mapping"]

    # Get root directory
    root_directory = configuration_dict["data_root_dir"]

    # Validate that `gear_data` exists in configuration settings
    if "gear_data" not in configuration_dict:
        raise KeyError(
            "Directory information for 'gear_data' missing from file configuration settings."
        )

    # Generate haul-to-transect map key
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
        data = pd.read_excel(filepath, sheet_name=sheet, usecols=BIODATA_HAUL_MAP.keys()).dropna(
            subset=BIODATA_HAUL_MAP.keys()
        )
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
        if verbose:
            print(
                f"Haul-to-transect mapping file for '{region}' saved at "
                f"'{Path(root_directory + save_dir) / save_file}'."
            )


def preprocess_biodata(input_dict: dict, configuration_dict: dict) -> None:
    """
    Preprocess just biological data

    Parameters
    ----------
    input_dict: dict
        Dictionary corresponding to the `input` attribute belonging to `Survey`-class
    configuration_dict: dict
        Dictionary corresponding to the `config` attribute belonging to `Survey`-class
    """

    # Generate length and age vectors
    # ---- Length vector
    length_bins = np.linspace(
        configuration_dict["biometrics"]["bio_hake_len_bin"][0],
        configuration_dict["biometrics"]["bio_hake_len_bin"][1],
        configuration_dict["biometrics"]["bio_hake_len_bin"][2],
        dtype=np.float64,
    )
    # ---- Age vector
    age_bins = np.linspace(
        configuration_dict["biometrics"]["bio_hake_age_bin"][0],
        configuration_dict["biometrics"]["bio_hake_age_bin"][1],
        configuration_dict["biometrics"]["bio_hake_age_bin"][2],
        dtype=np.float64,
    )

    # Discretize these values into discrete intervals
    # ---- Calculate binwidths
    # -------- Length
    length_binwidth = np.mean(np.diff(length_bins / 2.0))
    # -------- Age
    age_binwidth = np.mean(np.diff(age_bins / 2.0))
    # ---- Center the bins within the binwidths
    # -------- Length
    length_centered_bins = np.concatenate(
        ([length_bins[0] - length_binwidth], length_bins + length_binwidth)
    )
    # -------- Age
    age_centered_bins = np.concatenate(([age_bins[0] - age_binwidth], age_bins + age_binwidth))

    # Merge the vector and centered bins into dataframes that will be added into the `input`
    # attribute
    # ---- Generate DataFrame for length
    length_bins_df = pd.DataFrame({"length_bins": length_bins})
    # -------- Discretize the bins as categorical intervals
    length_bins_df["length_intervals"] = pd.cut(length_bins_df["length_bins"], length_centered_bins)
    # ---- Generate DataFrame for age
    age_bins_df = pd.DataFrame({"age_bins": age_bins})
    # -------- Discretize the bins as categorical intervals
    age_bins_df["age_intervals"] = pd.cut(age_bins_df["age_bins"], age_centered_bins)
    # ---- Update `input` attribute
    # -------- Length
    input_dict["biology"]["distributions"]["length_bins_df"] = length_bins_df
    # -------- Age
    input_dict["biology"]["distributions"]["age_bins_df"] = age_bins_df
    # -------- Delete the duplicate configuration keys
    # del configuration_dict["biometrics"]

    # Relabel sex to literal words among biological data
    # ---- Specimen
    input_dict["biology"]["specimen_df"]["sex"] = np.where(
        input_dict["biology"]["specimen_df"]["sex"] == int(1),
        "male",
        np.where(
            input_dict["biology"]["specimen_df"]["sex"] == int(2),
            "female",
            np.where(
                input_dict["biology"]["specimen_df"]["sex"].isin(["male", "female"]),
                input_dict["biology"]["specimen_df"]["sex"],
                "unsexed",
            ),
        ),
    )
    # -------- Sex group
    input_dict["biology"]["specimen_df"]["group_sex"] = np.where(
        input_dict["biology"]["specimen_df"]["sex"] != "unsexed", "sexed", "unsexed"
    )
    # ---- Length
    input_dict["biology"]["length_df"]["sex"] = np.where(
        input_dict["biology"]["length_df"]["sex"] == int(1),
        "male",
        np.where(
            input_dict["biology"]["length_df"]["sex"] == int(2),
            "female",
            np.where(
                input_dict["biology"]["length_df"]["sex"].isin(["male", "female"]),
                input_dict["biology"]["length_df"]["sex"],
                "unsexed",
            ),
        ),
    )
    # -------- Sex group
    input_dict["biology"]["length_df"]["group_sex"] = np.where(
        input_dict["biology"]["length_df"]["sex"] != "unsexed", "sexed", "unsexed"
    )

    # Discretize the age and length bins of appropriate biological data
    # ---- Specimen
    input_dict["biology"]["specimen_df"] = input_dict["biology"]["specimen_df"].bin_variable(
        [length_centered_bins, age_centered_bins], ["length", "age"]
    )
    # ---- Length
    input_dict["biology"]["length_df"] = input_dict["biology"]["length_df"].bin_variable(
        length_centered_bins, "length"
    )


def preprocess_spatial(input_dict: dict) -> None:
    """
    Preprocess just spatial data

    Parameters
    ----------
    input_dict: dict
        Dictionary corresponding to the `input` attribute belonging to `Survey`-class
    """

    # Update column names
    # ---- `geo_strata`
    input_dict["spatial"]["geo_strata_df"].columns = input_dict["spatial"][
        "geo_strata_df"
    ].columns.str.replace(" ", "_")
    # ---- `inpfc_strata`
    input_dict["spatial"]["inpfc_strata_df"].columns = input_dict["spatial"][
        "inpfc_strata_df"
    ].columns.str.replace(" ", "_")
    # ---- `inpfc_strata`: rename stratum column name to avoid conflicts
    input_dict["spatial"]["inpfc_strata_df"].rename(
        columns={"stratum_num": "stratum_inpfc"}, inplace=True
    )

    # Bin data
    # ---- Create latitude intervals to bin the strata
    latitude_bins = np.concatenate(
        [[-90], input_dict["spatial"]["inpfc_strata_df"]["northlimit_latitude"], [90]]
    )
    # ---- Add categorical intervals
    input_dict["spatial"]["inpfc_strata_df"]["latitude_interval"] = pd.cut(
        input_dict["spatial"]["inpfc_strata_df"]["northlimit_latitude"] * 0.99, latitude_bins
    )


def preprocess_acoustic_spatial(input_dict: dict) -> None:
    """
    Preprocess joint acoustic and spatial data

    Parameters
    ----------
    input_dict: dict
        Dictionary corresponding to the `input` attribute belonging to `Survey`-class
    """

    # Bin data
    # ---- Create latitude intervals to bin the strata
    latitude_bins = np.concatenate(
        [[-90], input_dict["spatial"]["inpfc_strata_df"]["northlimit_latitude"], [90]]
    )
    # ---- Bin NASC transects into appropriate INPFC strata
    input_dict["acoustics"]["nasc_df"]["stratum_inpfc"] = (
        pd.cut(
            input_dict["acoustics"]["nasc_df"]["latitude"],
            latitude_bins,
            right=True,
            labels=range(len(latitude_bins) - 1),
        )
    ).astype(int) + 1

    # KS strata
    # ---- Map hauls to `all_ages`
    input_dict["acoustics"]["nasc_df"].set_index("haul_all_ages", inplace=True)
    input_dict["acoustics"]["nasc_df"]["stratum_all_ages"] = (
        input_dict["spatial"]["strata_df"]
        .rename(columns={"haul_num": "haul_all_ages"})
        .set_index("haul_all_ages")["stratum_num"]
    )
    input_dict["acoustics"]["nasc_df"]["stratum_all_ages"] = input_dict["acoustics"]["nasc_df"][
        "stratum_all_ages"
    ].fillna(1)
    input_dict["acoustics"]["nasc_df"] = input_dict["acoustics"]["nasc_df"].reset_index()
    # ---- Map hauls to `no_age1`
    input_dict["acoustics"]["nasc_df"].set_index("haul_no_age1", inplace=True)
    input_dict["acoustics"]["nasc_df"]["stratum_no_age1"] = (
        input_dict["spatial"]["strata_df"]
        .rename(columns={"haul_num": "haul_no_age1"})
        .set_index("haul_no_age1")["stratum_num"]
    )
    input_dict["acoustics"]["nasc_df"]["stratum_no_age1"] = input_dict["acoustics"]["nasc_df"][
        "stratum_no_age1"
    ].fillna(1)
    input_dict["acoustics"]["nasc_df"] = input_dict["acoustics"]["nasc_df"].reset_index()


def preprocess_biology_spatial(input_dict: dict) -> None:
    """
    Preprocess joint biological and spatial data

    Parameters
    ----------
    input_dict: dict
        Dictionary corresponding to the `input` attribute belonging to `Survey`-class
    """

    # Merge haul numbers and spatial information across biological variables
    # ---- Create interval key for haul numbers to assign INPFC stratum
    haul_bins = np.sort(
        np.unique(
            np.concatenate(
                [
                    input_dict["spatial"]["inpfc_strata_df"]["haul_start"] - int(1),
                    input_dict["spatial"]["inpfc_strata_df"]["haul_end"],
                ]
            )
        )
    )
    # ---- Quantize the INPFC dataframe hauls based on strata
    input_dict["spatial"]["inpfc_strata_df"]["haul_bin"] = pd.cut(
        (
            input_dict["spatial"]["inpfc_strata_df"]["haul_start"]
            + input_dict["spatial"]["inpfc_strata_df"]["haul_end"]
        )
        / 2,
        haul_bins,
    )
    # ---- Rename `stratum_num` column
    input_dict["spatial"]["inpfc_strata_df"].rename(
        columns={"stratum_num": "stratum_inpfc"}, inplace=True
    )
    # ---- Set the index to `haul_bins`
    inpfc_df = input_dict["spatial"]["inpfc_strata_df"].copy().set_index(["haul_bin"])

    # Get the KS-strata
    strata_df = input_dict["spatial"]["strata_df"].copy().set_index(["haul_num"])

    # Loop through the KS-strata to map the correct strata values
    for keys, values in input_dict["biology"].items():
        if isinstance(values, pd.DataFrame) and "haul_num" in values.columns:
            # ---- Index based on `haul_num`
            input_dict["biology"][keys].set_index(["haul_num"], inplace=True)
            # ---- Map the correct `stratum_num` value
            input_dict["biology"][keys]["stratum_num"] = strata_df["stratum_num"]
            # ---- Correct the value datatype
            input_dict["biology"][keys]["stratum_num"] = (
                input_dict["biology"][keys]["stratum_num"].fillna(0.0).astype(int)
            )
            # ---- Reset the index
            input_dict["biology"][keys].reset_index(inplace=True)
            # ---- Bin for `stratum_inpfc`
            input_dict["biology"][keys]["haul_bin"] = pd.cut(
                input_dict["biology"][keys]["haul_num"], haul_bins
            )
            # ---- Set index to `haul_bins`
            input_dict["biology"][keys].set_index(["haul_bin"], inplace=True)
            # ---- Merge
            input_dict["biology"][keys]["stratum_inpfc"] = inpfc_df["stratum_inpfc"]
            # ---- Reset indices
            input_dict["biology"][keys].reset_index(inplace=True)
            # ---- Drop `haul_bin`
            input_dict["biology"][keys].drop(columns=["haul_bin"], inplace=True)


def preprocess_acoustic_biology_spatial(input_dict: dict, configuration_dict: dict) -> None:
    """
    Preprocess joint acoustic, biological, and spatial data

    Parameters
    ----------
    input_dict: dict
        Dictionary corresponding to the `input` attribute belonging to `Survey`-class
    configuration_dict: dict
        Dictionary corresponding to the `config` attribute belonging to `Survey`-class
    """

    # Identify haul sub-groups, if they exist
    groups = [
        col.replace("haul_", "")
        for col in input_dict["acoustics"]["nasc_df"].columns
        if "haul" in col
    ]
    # ---- Get unique species to be processed
    spp = configuration_dict["species"]["number_code"]
    # ---- Change to list, if needed
    spp = [spp] if not isinstance(spp, list) else spp
    # ---- Get the unique values in the table
    nonzero_trawls = list(
        set(
            val
            for df in input_dict["biology"].values()
            if isinstance(df, pd.DataFrame)
            and all(col in df.columns for col in ["haul_num", "species_id"])
            for val in df.loc[df["species_id"].isin(spp), "haul_num"].unique()
        )
    )
    # ---- Assign the correct KS-strata to NASC region groups
    for grp in groups:
        # ---- Index `nasc_df`
        input_dict["acoustics"]["nasc_df"].set_index(f"haul_{grp}", inplace=True)
        # ---- Create copy of `strata_df` and index
        strata_df = input_dict["spatial"]["strata_df"].copy()
        # ---- Remove zero-trawls
        strata_df = (
            strata_df[strata_df["haul_num"].isin(nonzero_trawls)]
            .filter(["stratum_num", "haul_num"])
            .rename(columns={"haul_num": f"haul_{grp}"})
            .set_index([f"haul_{grp}"])
        )
        # ---- Merge the strata information based on haul
        input_dict["acoustics"]["nasc_df"][f"stratum_{grp}"] = strata_df["stratum_num"]
        # ---- Fill empty with 0.0's
        input_dict["acoustics"]["nasc_df"][f"stratum_{grp}"] = (
            input_dict["acoustics"]["nasc_df"][f"stratum_{grp}"].fillna(0).astype(int)
        )
        # ---- Reset index
        input_dict["acoustics"]["nasc_df"].reset_index(inplace=True)


def preprocess_statistics(input_dict: dict, configuration_dict: dict) -> None:
    """
    Preprocess statistical data and settings

    Parameters
    ----------
    input_dict: dict
        Dictionary corresponding to the `input` attribute belonging to `Survey`-class
    configuration_dict: dict
        Dictionary corresponding to the `config` attribute belonging to `Survey`-class
    """

    # Reorganize kriging/variogram parameters
    # ---- Kriging
    # -------- Generate dictionary comprising kriging model configuration
    kriging_params = (
        input_dict["statistics"]["kriging"]["vario_krig_para_df"]
        .filter(regex="krig[.]")
        .rename(columns=lambda x: x.replace("krig.", ""))
        .rename(columns={"ratio": "anisotropy", "srad": "search_radius"})
        .to_dict(orient="records")[0]
    )
    # -------- Concatenate configuration settings for kriging
    kriging_params.update(configuration_dict["kriging_parameters"])
    # ---- Variogram
    # -------- Generate dictionary comprising variogram model configuration
    variogram_params = (
        input_dict["statistics"]["kriging"]["vario_krig_para_df"]
        .filter(regex="vario[.]")
        .rename(columns=lambda x: x.replace("vario.", ""))
        .rename(
            columns={
                "lscl": "correlation_range",
                "powr": "decay_power",
                "hole": "hole_effect_range",
                "res": "lag_resolution",
                "nugt": "nugget",
            }
        )
        .to_dict(orient="records")[0]
    )
    # ---- Update the input attribute with the reorganized parameters
    input_dict["statistics"]["variogram"].update({"model_config": variogram_params})
    input_dict["statistics"]["kriging"].update({"model_config": kriging_params})
    # -------- Delete the duplicate dataframe
    del input_dict["statistics"]["kriging"]["vario_krig_para_df"]


def prepare_input_data(input_dict: dict, configuration_dict: dict, imported_data: List[str]):
    """
    Rearranges and organizes data formats of the initial file inputs

    Parameters
    ----------
    input_dict: dict
        Dictionary corresponding to the `input` attribute belonging to `Survey`-class
    configuration_dict: dict
        Dictionary corresponding to the `config` attribute belonging to `Survey`-class
    imported_data: List[str]
        A list of datasets for determining which pre-processing steps to run
    """

    if "biology" in imported_data:
        preprocess_biodata(input_dict, configuration_dict)

    # SPATIAL
    if "spatial" in imported_data:
        preprocess_spatial(input_dict)

    # ACOUSTICS + SPATIAL
    if set(["acoustics", "spatial"]).issubset(imported_data):
        preprocess_acoustic_spatial(input_dict)

    # BIOLOGY + SPATIAL
    if set(["biology", "spatial"]).issubset(imported_data):
        preprocess_biology_spatial(input_dict)

    # ACOUSTICS + SPATIAL + BIOLOGICAL
    if set(["acoustics", "biology", "spatial"]).issubset(imported_data):
        preprocess_acoustic_biology_spatial(input_dict, configuration_dict)

    # STATISTICS
    if (
        set(["statistics"]).issubset(imported_data)
        and "vario_krig_para_df" in input_dict["statistics"]["kriging"]
    ):
        preprocess_statistics(input_dict, configuration_dict)


def validate_config_structure(yaml_data, config_spec):
    """
    Validate configuration YAML dictionary structure and entry types
    """

    # Helper function for validating dictionary key names/structure
    def validate_dict(data, spec, branch=""):
        for key, value in spec.items():
            current_branch = f"{branch}.{key}" if branch else key
            if key == "ANY":
                # If the spec key is "ANY", check all keys in the data against the "ANY" spec
                if not isinstance(data, dict):
                    raise ValueError(
                        f"Expected a dictionary at {current_branch}, got {type(data).__name__}"
                    )
                for sub_key in data:
                    validate_value(data[sub_key], spec[key], f"{current_branch}.{sub_key}")
            else:
                if key not in data:
                    raise KeyError(f"Missing key: {current_branch}")
                validate_value(data[key], value, current_branch)

    # Helper function for parsing values entered at different points throughout the dictionary
    def validate_value(data, spec, branch):
        if isinstance(spec, dict):
            validate_dict(data, spec, branch)
        elif isinstance(spec, type):
            if not is_valid_type(data, spec):
                raise TypeError(
                    f"Expected type {spec.__name__} at {branch}, got {type(data).__name__}"
                )
        elif isinstance(spec, list):
            if not isinstance(data, list):
                raise TypeError(f"Expected a list at {branch}, got {type(data).__name__}")
            if len(spec) != 1:
                raise ValueError(
                    f"Spec list at {branch} should contain exactly one type element, got {spec}"
                )
            for index, item in enumerate(data):
                validate_value(item, spec[0], f"{branch}[{index}]")
        else:
            raise ValueError(f"Unknown spec type at {branch}: {spec}")

    # Helper function for assessing the data type of each nested value
    def is_valid_type(data, spec_type):
        if spec_type is float:
            return isinstance(data, (int, float))
        return isinstance(data, spec_type)

    try:
        validate_dict(yaml_data, config_spec)
    except (KeyError, TypeError, ValueError) as e:
        print(f"Validation error: {e}")


def dataset_integrity(
    input_dict: dict,
    analysis: Literal[
        "transect", "stratified:kriging", "stratified:transect", "variogram", "kriging"
    ],
):
    """
    Determine whether all of the necessary datasets are contained within the `Survey`-class object
    for each analysis

    Parameters
    ----------
    input_dict: dict
        Dictionary corresponding to the `input` attribute belonging to `Survey`-class
    analysis: Literal["transect", "stratified", "variogram", "kriging"]
        The name of the analysis to be performed
    """

    # Map the imported datasets
    imported_data = map_imported_datasets(input_dict)

    # Initialize `missing`
    missing = None

    # Transect analysis
    if analysis == "transect":
        # ---- Create expected datasets list
        expected_datasets = ["acoustics", "biology", "spatial"]
        # ---- Missing analysis str
        missing_analysis_method = "'transect_analysis'"

    # Stratified analysis (transect)
    if analysis == "stratified:transect":
        # ---- Create expected datasets list
        expected_datasets = ["acoustics", "spatial"]
        # ---- Missing analysis str
        missing_analysis_method = "'stratified_analysis' (for transect data)"

    # Stratified analysis (kriging)
    if analysis == "stratified:kriging":
        # ---- Create expected datasets list
        expected_datasets = ["acoustics", "spatial", "statistics"]
        # ---- Missing analysis str
        missing_analysis_method = "'stratified_analysis' (for kriged data)"

    # Kriging analysis
    if analysis == "kriging":
        # ---- Create expected datasets list
        expected_datasets = ["acoustics", "biology", "spatial", "statistics"]
        # ---- Missing analysis str
        missing_analysis_method = "'fit_variogram'/'variogram_gui'"

    # Variogram analysis
    if analysis == "variogram":
        # ---- Create expected datasets list
        expected_datasets = ["acoustics", "spatial", "statistics"]
        # ---- Missing analysis str
        missing_analysis_method = "'fit_variogram'/'variogram_gui'"

    # Raise Error, if appropriate
    if not set(expected_datasets).issubset(imported_data):
        # ---- Collect missing values
        missing = set(expected_datasets) - set(imported_data)
        # ---- Format string
        missing_str = ", ".join([f"'{i}'" for i in missing])
        # ---- Raise error
        raise ValueError(
            f"The following input datasets are missing for {missing_analysis_method}: "
            f"{missing_str}."
        ).with_traceback(None)
