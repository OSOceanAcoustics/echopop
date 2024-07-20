import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
from geopy.distance import distance
from shapely.geometry import Polygon, Point, box
import geopandas as gpd
from shapely.ops import unary_union
import pyproj
import geopy
from echopop.spatial.projection import wgs84_to_utm, utm_string_generator
import shapely.geometry
from echopop.survey import Survey
survey = Survey( init_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config.yml" ,
                 survey_year_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_config.yml" )


grid_settings = file_configuration["geospatial"]["griddify"]
# lat_min = grid_settings["bounds"]["latitude"][0]
lat_min = 33.75
# lat_max = grid_settings["bounds"]["latitude"][1]
lat_max = 55.50
# lon_min = grid_settings["bounds"]["longitude"][0]
lon_min = -134.25
lon_max = grid_settings["bounds"]["longitude"][1]

projection = file_configuration["geospatial"]["projection"]

utm_code = utm_string_generator((lon_max + lon_min)/2, (lat_max + lat_min)/2)
utm_num = int(utm_code)
utm_str = f"epsg:{utm_num}"

biology_data = filtered_biology_output

from sqlalchemy import create_engine, text, Engine, inspect
root_dir = file_configuration["data_root_dir"]
db_directory = Path(root_dir) / "database"
db_directory.mkdir(parents=True, exist_ok=True)
db_file = db_directory / "biology.db"
# Create the engine with the full path
engine = create_engine(f'sqlite:///{db_file}')

SQL_COMMANDS = {
    "create": "CREATE TABLE IF NOT EXISTS {table_name} ({column_definitions});",
    "check": "SELECT name FROM sqlite_master WHERE type='table' AND name='{table_name}';",
    "drop": "DROP TABLE IF EXISTS {table_name};",
    "select": "SELECT {columns} FROM {table_name};",
    "index": "CREATE UNIQUE INDEX IF NOT EXISTS {index_name} ON {table_name} ({columns})",
    # "insert": "INSERT INTO {table_name} ({columns});",
    "insert": """
        INSERT INTO {table_name} ({columns}) 
        SELECT {columns} 
        FROM (SELECT VALUES {values} FROM (VALUES {value_placeholder})) AS source ({columns}) 
        {filter_clause};
        """,
    "inspect": None,
}

SQL_DTYPES = {
    'int32': 'INTEGER',
    'int64': 'INTEGER',
    'float64': 'FLOAT',
    'bool': 'BOOLEAN',
    'datetime64[ns]': 'DATETIME',
    'object': 'TEXT'
}

def SQL(db_file: str, command: str, **kwargs):

    # Create engine from `db_file` string
    engine = create_engine(f"sqlite:///{db_file}")

    # Format `columns`, if there are any and more than 1
    if "columns" in kwargs.keys():
        if isinstance(kwargs["columns"], list):
            kwargs["columns"] = ", ".join(kwargs["columns"])
    else:
        kwargs["columns"] = "*"

    # Format `columns`, if there are any and more than 1
    # if "filter_columns" in kwargs.keys():
    #     # ---- Store the value for later
    #     kwargs["filter_columns_store"] = kwargs["filter_columns"]
    #     if isinstance(kwargs["filter_columns"], list):
    #         kwargs["filter_columns"] = ", ".join(kwargs["filter_columns"])

    # Run the command
    try:
        with engine.connect() as connection:
            # ---- SELECT
            if command == "select":
                return pd.read_sql(text(SQL_COMMANDS[command].format(**kwargs)), con=connection)
            # ---- CREATE
            elif command == "create":
                # ---- Extract dataframe
                df_to_add = kwargs["dataframe"]
                # ---- Check whether the table already exists or not
                table_exists = (
                    connection.execute(text(SQL_COMMANDS["check"].format(**kwargs))).fetchone()
                )
                # ---- If it doesn't, pre-allocate the table 
                if table_exists is None:
                    # ---- Get column definitions as a string
                    column_def_dict = {
                        col: SQL_DTYPES.get(str(dtype), 'TEXT') 
                        for col, dtype in zip(df_to_add.columns, df_to_add.dtypes)
                    }
                    # ---- Convert to a single string
                    kwargs["column_definitions"] = (
                        ", ".join([f"{col} {dtype}" for col, dtype in column_def_dict.items()])
                    )
                    # ---- Create table
                    connection.execute(text(SQL_COMMANDS["create"].format(**kwargs)))
            # ---- REPLACE
            elif command == "replace":
                # ---- Extract dataframe
                df_to_add = kwargs["dataframe"]
                # ---- Replace current
                df_to_add.to_sql(name=kwargs["table_name"], 
                                 con=connection, 
                                 if_exists="replace", index=False)

            # ---- INSERT
            elif command == "insert": 
                # ---- Extract dataframe
                df_to_add = kwargs["dataframe"]
                # ---- Check if 
                # table_exists = (
                #     connection.execute(text(SQL_COMMANDS["check"].format(**kwargs))).fetchone()
                # )
                # tables = SQL(db_file, "inspect")
                # ---- If it doesn't, pre-allocate the table 
                # if kwargs["table_name"] not in tables and "filter_columns" in kwargs.keys():
                df_to_add.to_sql(name=kwargs["table_name"], 
                                    con=connection, 
                                    if_exists="append", index=False)
                # else:
                    #     # ---- Format `filter_columns` command if present
                    # if "filter_columns" in kwargs.keys():
                    #     # ---- Fetch table
                    #     fetch_table = (
                    #         connection.execute(text(
                    #             ("SELECT DISTINCT {filter_columns} FROM {table_name}")
                    #             .format(**kwargs))
                    #         )
                    #     )
                    #     # ---- Format the SQL data into a DataFrame
                    #     fetched_df = pd.DataFrame(fetch_table.fetchall(), columns=fetch_table.keys())              
                    #     # ---- Create an index tuples
                    #     index_tuples = (
                    #         set(fetched_df[kwargs["filter_columns_store"]]
                    #             .itertuples(index=False, name=None))
                    #     )
                    #     # ---- Filter the dataframe
                    #     filtered_df = (
                    #         df_to_add[
                    #             ~df_to_add[fetched_df.columns].apply(tuple, axis=1)
                    #             .isin(index_tuples)
                    #             ]
                    #     )
                    #     # ---- Insert the data
                    #     filtered_df.to_sql(name=kwargs["table_name"], 
                    #                         con=connection, 
                    #                         if_exists="append", index=False)
                    # else:
                    # df_to_add.to_sql(name=kwargs["table_name"], 
                    #                 con=connection, 
                    #                 if_exists="append", index=False)
            # ---- INSPECT
            elif command == "inspect":
                return inspect(engine).get_table_names()
            else: 
                connection.execute(text(SQL_COMMANDS[command].format(**kwargs)))
    finally: 
        # ---- Dispose of the engine to release any resources being pooled/used
        engine.dispose()

_ = SQL(db_file, "drop", table_name="catch_df")
_ = SQL(db_file, "drop", table_name="specimen_df")
_ = SQL(db_file, "drop", table_name="length_df")
_ = SQL(db_file, "drop", table_name="files_read")

_ = SQL(db_file, "insert", table_name="files_read", dataframe=current_files)
current = SQL(db_file, "select", table_name="files_read", columns="filepath")
current


# Get acoustic directory and initialization settings
# ---- Files
biology_file_settings = file_configuration["input_directories"]["biological"]
# ---- General settings
biology_analysis_settings = file_configuration["biology"]

# Get the file-specific settings, datatypes, columns, etc.
# ---- Get defined columns and datatypes from `LIVE_INPUT_FILE_CONFIG_MAP`
biology_config_map = LIVE_INPUT_FILE_CONFIG_MAP["biology"]
# ---- Extract the expected file name ID's
biology_file_ids = biology_file_settings["file_name_formats"]
# ---- Extract all of the file ids
biology_config_ids = list(biology_file_ids.keys())
# ---- Initialize the dictionary that will define this key in the `input` attribute
biology_output = {f"{key}_df": pd.DataFrame() for key in biology_config_ids}
# ---- Initialize the SQL dictionary
sql_biology_output = {f"{key}_df": pd.DataFrame() for key in biology_config_ids}

# Create full filepath
biology_directory_path = (
    Path(file_configuration["data_root_dir"]) / biology_file_settings["directory"]
)
# ---- Directory check
directory_existence = biology_directory_path.exists()
# ---- Error evaluation (if applicable)
if not directory_existence:
    raise FileNotFoundError(
        f"The acoustic data directory [{biology_directory_path}] does not exist."
    )
# ---- Get the defined file extension
file_extension = biology_file_settings["extension"]
# ---- Create Path.glob generator object
file_path_obj = biology_directory_path.glob(f"*{'.'+file_extension}")
#---- Create list of `*.csv`` files
csv_files = list(file_path_obj)
# ---- Ensure files exist or raise error otherwise
if len(csv_files) < 1:
    raise FileNotFoundError(
        f"No `*.csv` files found in [{biology_directory_path}]!"
    )
else: 
    # ---- Create Path to SQL database file
    db_directory = Path(file_configuration["data_root_dir"]) / "database"
    # ---- Create the directory if it does not already exist
    db_directory.mkdir(parents=True, exist_ok=True)
    # ---- Complete path to `biology.db`
    db_file = db_directory / "biology.db"
    # ---- Query the external SQL database to see if the file tracking table exists
    tables = SQL(db_file, "inspect")
    # ---- Create a list of string-formatted Path names
    csv_files_str = [str(file) for file in csv_files]
    # ---- Create DataFrame
    current_files = pd.DataFrame(csv_files_str, columns=["filepath"])
    # ---- Create if it is missing and then advance `csv_files`
    if "files_read" not in tables:
        # ---- Insert into the SQL database file
        _ = SQL(db_file, "insert", table_name="files_read", columns="filepath",
                    dataframe=current_files)        
        # ---- Create empty list for later comparison
        new_files = []
    else:
        # ---- Pull already processed filenames
        previous_files = SQL(db_file, "select", table_name="files_read")
        # ---- Compare against the current filelist 
        new_files = (
            [file for file in csv_files_str if file not in set(previous_files["filepath"])]
        )  
        # ---- Create a DataFrame for the new files
        new_files_df = pd.DataFrame(new_files, columns=["filepath"])
        # ---- Insert into the SQL database file
        _ = SQL(db_file, "insert", table_name="files_read", dataframe=new_files_df) 

# Iterate through each of the file ids and read in the data 
for id in list(biology_file_ids.keys()): 
    # ---- Extract the specific config mapping for this tag/id
    sub_config_map = biology_config_map[id]
    # ---- Drop the `{FIELD_ID}` tag identifier
    file_id_format = re.sub(r'\{FILE_ID:([^}]+)\}', r'\1', biology_file_ids[id])
    # ---- Replace all other tags with `*` placeholders
    file_id_format = re.sub(r"\{[^{}]+\}", "*", file_id_format)
    # ---- Create Path object with the generalized format
    subfile_path_obj = biology_directory_path.glob(f"{file_id_format}.{file_extension}")
    # ---- List all files that match this pattern
    subcsv_files_str = [str(file) for file in list(subfile_path_obj)]
    # ---- Filter for only new files
    subset_files = set(subcsv_files_str).intersection(set(new_files))
    # ---- Pull from SQL database, if applicable
    if f"{id}_df" in tables:
        # ---- SELECT
        sql_df = SQL(db_file, "select", table_name=f"{id}_df", columns="*")
        # ---- Concatenate to the dictionary
        sql_biology_output[f"{id}_df"] = pd.concat([biology_output[f"{id}_df"], sql_df])
    # ---- Add data files not stored in SQL database
    if len(subset_files) > 0 or len(subset_files)== 0 and f"{id}_df" not in tables:
        if len(subset_files) > 0:
            file_list = subset_files
        else:
            file_list = subcsv_files_str
        # ---- Create a list of relevant dataframes
        sub_df_lst = [read_biology_csv(Path(file), biology_file_ids[id], sub_config_map) 
                        for file in file_list]
        # ---- Concatenate into a single DataFrame
        sub_df = pd.concat(sub_df_lst, ignore_index=True)
        # ---- Concatenate to the dictionary DataFrame
        biology_output[f"{id}_df"] = pd.concat([biology_output[f"{id}_df"], sub_df])

# Get contrasts used for filtering the dataset
# ---- Species
species_filter = file_configuration["species"]["number_code"]
# ---- Trawl partition information
trawl_filter = biology_analysis_settings["catch"]["partition"]
# ---- Apply the filter
filtered_biology_output = {
    key: df[
        (df['species_id'] == species_filter if 'species_id' in df.columns else True) &
        (df['trawl_partition'].str.lower() == trawl_filter if 'trawl_partition' in df.columns else True)
    ]
    for key, df in biology_output.items() if isinstance(df, pd.DataFrame) and not df.empty
}

# Update the SQL database
for table_name, df in filtered_biology_output.items():
    # ---- Update        
    _ = SQL(db_file, "insert", table_name=table_name, columns="*", 
            dataframe=df)
    
# Combine the two datasets 
merged_output = {
    key: pd.concat([
        sql_biology_output.get(key, pd.DataFrame()), 
        filtered_biology_output.get(key, pd.DataFrame())
    ]).drop_duplicates().reset_index(drop=True)
    for key in set(sql_biology_output) | set(filtered_biology_output)
}
# ---- Return output
merged_output

coordinate_metadata.attrs[]

SQL(biology_db, command="drop", table_name="catch_df")
SQL(biology_db, command="drop", table_name="specimen_df")
SQL(biology_db, command="drop", table_name="length_df")
SQL(biology_db, command="drop", table_name="files_read")
_ = SQL(db_file=db_file, command="create", table_name="files_read", columns="filepath")
tables = SQL(db_file, "inspect")
tables
current = SQL(db_file, "select", table_name="files_read", columns=["filepath"])
current

SQL(db_file, "select", table_name="catch_df", columns="*")
new_files_df = pd.DataFrame(csv_files_str, columns=['file_path'])
_ = SQL("insert", engine, table_name="files_read",dataframe=new_files_df)
current = SQL("select", engine, table_name="csv_files_read", columns="file_path")
current
for table_name, df in biology_data.items():
    df.to_sql(table_name, con=engine, if_exists='append', index=False)
command = "read"
engine = create_engine(f'sqlite:///{db_file}')
table_name = "files_read"
columns = "file_path"

kwargs = {
    "table_name": table_name,
    "columns": columns, 
}

zarr_data_ds["depth"].diff(dim="depth")

prc_nasc_df.groupby(["longitude", "latitude"])

from pandas.core.groupby import DataFrameGroupBy

def estimate_echometrics(acoustic_data_df: pd.DataFrame):

    # Create copy
    acoustic_df = acoustic_data_df.copy().reset_index(drop=True)

    # Pre-compute the change in depth
    acoustic_df["dz"] = acoustic_df["depth"].diff()

    # Initialize echometrics dictionary
    echometrics = {}

    # Compute the metrics center-of-mass
    if acoustic_df["NASC"].sum() == 0.0:
        echometrics.update({
            "n_layers": 0,
            "mean_Sv": -999,
            "max_Sv": -999,
            "nasc_db": np.nan,
            "center_of_mass": np.nan,
            "dispersion": np.nan,
            "evenness": np.nan,
            "aggregation": np.nan,    
            "occupied_area": 0.0,        
        })
    else:
        
        # Compute the number of layers
        echometrics.update({
            "n_layers": acoustic_df["depth"][acoustic_df["NASC"] > 0.0].size
        })

        # Compute ABC
        # ---- Convert NASC to ABC
        acoustic_df["ABC"] = acoustic_df["NASC"] / (4 * np.pi * 1852 ** 2)
        # ---- Estimate mean Sv
        echometrics.update({
            "mean_Sv": 10.0 * np.log10(acoustic_df["ABC"].sum() / acoustic_df["depth"].max())
        })
        # --- Estimate max Sv (i.e. )
        echometrics.update({
            "max_Sv": 10 * np.log10(acoustic_df["ABC"].max() 
                                    / acoustic_df.loc[np.argmax(acoustic_df["ABC"]), "dz"])
        })

        # Compute (acoustic) abundance
        echometrics.update({
            "nasc_db": 10 * np.log10(acoustic_df["ABC"].sum())
        })

        # Compute center of mass
        echometrics.update({
            "center_of_mass": (
                (acoustic_df["depth"] * acoustic_df["NASC"]).sum()
                / (acoustic_df["NASC"]).sum()
            )
        })

        # Compute the dispersion
        echometrics.update({
            "dispersion": (
                ((acoustic_df["depth"] - echometrics["center_of_mass"]) ** 2 
                * acoustic_df["NASC"]).sum() / (acoustic_df["NASC"]).sum()                
            )
        })

        # Compute the evenness
        echometrics.update({
            "evenness": (acoustic_df["NASC"] **2).sum() / ((acoustic_df["NASC"]).sum()) ** 2
        })

        # Compute the index of aggregation
        echometrics.update({
            "aggregation": 1 / echometrics["evenness"]
        })

        # Get the occupied area
        echometrics.update({
            "occupied_area": (
                acoustic_df["dz"][acoustic_df["ABC"] > 0.0].sum() / acoustic_df["depth"].max()
            )
        })

    # Return the dictionary
    return echometrics

def integrate_nasc(acoustic_data_df: pd.DataFrame, echometrics: bool = True):

    # Vertically integrate PRC NASC
    nasc_dict = {"nasc": acoustic_data_df["NASC"].sum()}
    
    # Horizontally concatenate `echometrics`, if `True`
    if echometrics:
        # ---- Compute values
        # NOTE: This uses NASC instead of linear `sv`
        echometrics_dict = estimate_echometrics(acoustic_data_df)
        # ---- Merge
        nasc_dict.update(echometrics_dict)

    # Convert `nasc_dict` to a DataFrame and return the output
    return pd.Series(nasc_dict)

def process_group(group):
    result = integrate_nasc(group, echometrics=True)
    result = result.reset_index(drop=True)
    # Concatenate the result back to the original group for alignment
    group = group.reset_index(drop=True)
    combined = pd.concat([group, result], axis=1)
    return combined

acoustic_data_df = acoustic_data["prc_nasc_df"]


rc_nasc_df[prc_nasc_df["distance"] == 0.0]
acoustic_data_df = mek[mek["distance"] == 0.0]
pd.DataFrame(nasc_dict, index=[0]).reset_index(drop=True).unstack()
nasc_data_df = (
    prc_nasc_df.groupby(["longitude", "latitude", "ping_time"])
    .apply(lambda group: integrate_nasc(group, echometrics=False), include_groups=False)
    .reset_index()
)




kwargs = {
    "table_name": "csv_files_read",
    "columns": "file_path",
    "dataframe": new_files_df
}

current_process = psutil.Process()
import logging 

# Create a session
Session = sessionmaker(bind=engine)
session = Session()

# Perform database operations
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
logger.info("Performing database operations")

# Create a session
Session = sessionmaker(bind=engine)
session = Session()

# Perform database operations
logger.info("Performing database operations")

# Close the session
session.close()
logger.info("Session closed")

# Dispose the engine
engine.dispose()
logger.info("Engine disposed")

# Force garbage collection
import gc
gc.collect()
logger.info("Garbage collection performed")

import psutil

pid = psutil.Process().pid
process = psutil.Process(pid)
open_files = process.open_files()
db_path = r'C:\Users\Brandyn\Documents\GitHub\EchoPro_data\live_2019_files\database\biology.db'

# Check if the file is still in use
for file in open_files:
    if db_path in file.path:
        logger.info(f"File {db_path} is still in use.")
    else:
        logger.info(f"File {db_path} is not in use.")

# Define the SQL to drop the table
drop_table_sql = "DROP TABLE IF EXISTS csv_files_read;"
# Execute the drop table SQL
with engine.connect() as connection:
    _ = connection.execute(text(drop_table_sql))

import sqlite3
if os.path.exists(db_path):
    conn = sqlite3.connect(db_path)
    conn.close()
    # Force the file to be removed
    try:
        os.remove(db_path)
        print(f"Database file {db_path} has been deleted.")
    except PermissionError:
        print(f"Failed to delete {db_path}. The file is still in use.")
        
create_table_sql = """
CREATE TABLE IF NOT EXISTS csv_files_read (
    file_path TEXT UNIQUE
);
"""
# Execute the create table SQL
with engine.connect() as connection:
    _ = connection.execute(text(create_table_sql))

root_directory =  Path(root_dir)
dataset = "biology"

# Convert to strings
csv_files_str = [str(file) for file in csv_files]

existing_files_df = pd.read_sql('SELECT file_path FROM csv_files_read', con=engine)
existing_files_set = set(existing_files_df['file_path'])
# Filter out duplicates from the csv_files list
new_files = [file for file in csv_files_str if file not in existing_files_set]
# Insert only new file paths into the SQL table
if new_files:
    new_files_df = pd.DataFrame(new_files, columns=['file_path'])
    _ = new_files_df.to_sql('csv_files_read', con=engine, if_exists='append', index=False)


with engine.connect() as conn:
    conn.execute("""
        CREATE TABLE IF NOT EXISTS csv_files_read (
            file_path TEXT UNIQUE
        )
    """)

csv_files
files_df.to_sql('csv_files_read', con=engine, if_exists='append', index=False)
file_name_format = biology_file_ids[id]
def compile_filename_format(file_name_format: str):

    # Create a copy of `file_name_format`
    regex_pattern = file_name_format
    
    # Iterate through the keys from `LIVE_FILE_FORMAT_MAP` to format a regex pattern
    for key, value in LIVE_FILE_FORMAT_MAP.items():
        regex_pattern = regex_pattern.replace(f"{{{key}}}", value["expression"])
    # ---- Replace the `FILE_ID` tag
    regex_pattern = re.sub(r'\{FILE_ID:(.+?)\}', r'(?P<FILE_ID>\1)', regex_pattern)

    # Compile the regex pattern and return the output
    return re.compile(regex_pattern)

from sqlalchemy.orm import sessionmaker
Session = sessionmaker(bind=engine)
session = Session()
session.close()
engine.pool.status()
# Dispose the engine to close all connections
engine.dispose()
import gc
gc.collect()
import psutil
dbapi_conn = engine.raw_connection()
dbapi_conn.close()
# Get the process ID of the current process
pid = psutil.Process().pid

# List all open files for the current process
process = psutil.Process(pid)
open_files = process.open_files()

for file in open_files:
    print(file.path)


pattern = filename_format
config_settings = sub_config_map
regex_pattern = pattern

# Replace patterns based on LIVE_FILE_FORMAT_MAP
for key, value in LIVE_FILE_FORMAT_MAP.items():
    regex_pattern = regex_pattern.replace(f'{{{key}}}', value['expression'])
regex_pattern = re.sub(r'\{FILE_ID:(.+?)\}', r'(?P<FILE_ID>\1)', regex_pattern)
new_pattern = compile_filename_format(regex_pattern)
match_obj = new_pattern.search(file.name)
# Get substring components as a list
filename_substrings = re.findall(r'\{([^:}]+)(?::[^}]+)?}', pattern)
valid_tags = list(set(["HAUL", "SPECIES_CODE"]).intersection(set(filename_substrings)))

for i in valid_tags: 
    matched_key = LIVE_FILE_FORMAT_MAP[i]
    df[matched_key["name"]] = matched_key["dtype"](match_obj.group(i))



# Assign the data as new columns to the DataFrame
for key, value in data_to_add.items():
    df[key] = value

for i in valid_tags: 
    matched_key = LIVE_FILE_FORMAT_MAP[i]
    df[matched_key["name"]] = matched_key["dtype"](match_obj.group(i))
biology_analysis_settings
species_id_value = 22500
trawl_partition_value = 'Codend'  # Adjust as needed
{
    key: df[
        (('species_id' not in df.columns) or (df['species_id'] == species_id_value)) &
        (('trawl_partition' not in df.columns) or (df['trawl_partition'] == trawl_partition_value))
    ]
    for key, df in biology_output.items() if isinstance(df, pd.DataFrame)
}

(match_obj.group(i)).astype(matched_key["dtype"])
pattern = '{DATE:YYYYMM}_{HAUL}_{FILE_ID:catch_perc}'
modified_pattern = re.sub(r'\{FILE_ID:(.+?)\}', r'\1', pattern)
# Create the regex pattern
regex_pattern = modified_pattern.replace('{', '(?P<').replace('}', '>.+?)')
re.compile(regex_pattern)

modified_pattern = re.sub(r'\{FILE_ID:(.+?)\}', r'\1', pattern)
    
# Create the regex pattern
regex_pattern = modified_pattern.replace('{', '(?P<').replace('}', '>.+?)')
compile_filename_format(regex_pattern)
# Regular expression to capture values inside the curly braces
regex = r'\{([^:}]+):([^}]+)\}'

# Find all matches
matches = re.findall(regex, modified_pattern)

# Get substring components as a list
filename_substrings = re.findall(r'\{([^:}]+)(?::[^}]+)?}', pattern)

pattern_changed = pattern.replace("FILE_ID:", "")

# Compilte the filename regular expression format
compiled_regex = compile_filename_format(pattern_changed)

file_id_tag = pattern.split('{FILE_ID:')[1].split('}')[0]

 # Get the file name and produce a `re.Match` object
match_obj = compiled_regex.search(file.name)


def read_biology_csv(file: Path, pattern: re.Pattern, config_settings: dict):

    # Get the file name and produce a `re.Match` object
    match_obj = pattern.search(file.name)

    # Read in the `*.csv` file
    df = pd.read_csv(file, usecols=list(config_settings["dtypes"].keys()))

    # Validate the dataframe
    # ---- Check for any missing columns
    missing_columns = (
        [key for key in config_settings["dtypes"].keys() if key not in df.columns]
    )
    # ---- Raise Error, if needed
    if missing_columns: 
        raise ValueError(
            f"The following columns are missing from [{file}]: {', '.join(missing_columns)}!"
        )
    # ---- Ensure the correct datatypes
    df_validated = df.astype(config_settings["dtypes"])

    # Replace column names and drop 
    df_validated = df_validated.rename(columns=config_settings["names"])

    # Get the haul number and add the the dataframe
    # ---- Extract the haul number and convert to an integer
    haul_num = int(match_obj.group("HAUL"))
    # ---- Add the column
    df_validated["haul_num"] = haul_num

    # Return the resulting DataFrame
    return df_validated

##
grid_settings["grid_resolution"]["x"] = 50
grid_settings["grid_resolution"]["y"] = 50
lat_step = distance(nautical=grid_settings["grid_resolution"]["x"]).meters
lon_step = distance(nautical=grid_settings["grid_resolution"]["y"]).meters

# CREATE BOUNDING
bound_df = pd.DataFrame({
    "lon": np.array([lon_min, lon_max, lon_max, lon_min, lon_min]),
    "lat": np.array([lat_min, lat_min, lat_max, lat_max, lat_min])
})

bound_gdf = gpd.GeoDataFrame(
    data=bound_df,
    geometry=gpd.points_from_xy(bound_df["lon"], bound_df["lat"]),
    crs = projection
)

utm_string_generator(-117.0, 33.75)
bound_gdf.total_bounds
# Convert to UTM
bound_utm = bound_gdf.to_crs(utm_num)
bound_utm.total_bounds
y_step = lat_step
x_step = lon_step
# bound_utm = bound_gdf
# y_step = grid_settings["grid_resolution"]["y"] * 1852 / 110574
# x_step = grid_settings["grid_resolution"]["x"] * 1852 / 60.0

xmin, ymin, xmax, ymax = bound_utm.total_bounds

# Get number of cells
n_x_cells = int(np.ceil((xmax - xmin) / x_step))
n_y_cells = int(np.ceil((ymax - ymin) / y_step))

import pyproj
# create the cells in a loop
# grid_cells = []
# for x0 in np.arange(xmin, xmax, x_step):
#     for y0 in np.arange(ymin, ymax, y_step):
#         # bounds
#         utm_zone = utm_string_generator(x0, y0)
#         proj = pyproj.Proj(f"epsg:{utm_code}")
#         x1 = x0-x_step
#         y1 = y0+y_step
#         grid_cells.append(shapely.geometry.box(x0, y0, x1, y1))

grid_cells = []
for y0 in np.arange(ymin, ymax, y_step):

    # x_step = grid_settings["grid_resolution"]["x"] * 1852 / (1852 * 60 * np.cos(np.radians(y0)))

    for x0 in np.arange(xmin, xmax, x_step):
        # bounds
        # utm_zone = utm_string_generator(x0, y0)
        # proj = pyproj.Proj(f"epsg:{utm_code}")
        # x1, y1 = proj(x0, y0)
        # x2, y2 = proj(x0 - x_step, y0 + y_step)
        # grid_cells.append(box(x1, y1, x2, y2))
        x1 = x0-x_step
        y1 = y0+y_step
        grid_cells.append(shapely.geometry.box(x0, y0, x1, y1))

cells_gdf = gpd.GeoDataFrame(grid_cells, columns=["geometry"], crs=utm_code)
cells_gdf.shape
n_x_cells * n_y_cells
# cells_gdf = gpd.GeoDataFrame(grid_cells, columns=["geometry"])
cells_gdf.total_bounds
cells_gdf.to_crs(projection).total_bounds
from shapely.validation import make_valid
from shapely.geometry import mapping
########
world = gpd.read_file("C:/Users/Brandyn/Documents/GitHub/EchoPro_data/live_2019_files/coastline/ne_10m_land/ne_10m_land.shp")
bb_orig = box(lon_min, lat_min, lon_max, lat_max)
boundary_box = box(lon_min - 5, lat_min - 5, lon_max + 5, lat_max + 5)
world_orig = gpd.clip(world, box(lon_min-1, lat_min-1, lon_max+1, lat_max+1))
world_clipped_latlon = gpd.clip(world, boundary_box)
world_clipped = gpd.clip(world, boundary_box).to_crs(utm_code)

world_utm = world.to_crs(utm_code)
world_utm = world_utm[~world_utm.is_empty]

bbox_latlon = box(lon_min, lat_min, lon_max, lat_max)

gpd.GeoDataFrame(geometry=[bbox_latlon], crs=projection).to_crs(utm_code)

bbox_utm = bound_utm.total_bounds

buffer = [-lon_step * 1.01, -lat_step * 1.01, lon_step * 1.01, lat_step * 1.01]
array_buffer = bbox_utm + buffer
array_names = ["minx", "miny", "maxx", "maxy"]
buffered = dict(zip(array_names, array_buffer))
buffer_boundary = box(**buffered)
# box(array_buffer[0], array_buffer[1], array_buffer[2], array_buffer[3])
# buffer_boundary = buffer_boundary.to_crs(world_utm.crs)

buffer_boundary_gdf = gpd.GeoDataFrame(geometry=[buffer_boundary], crs=world_utm.crs)  # Replace with the correct EPSG code
bb_orig_gdf = gpd.GeoDataFrame(geometry=[bb_orig], crs=projection) 
# sub_clipped = gpd.clip(world_utm, buffer_boundary)
# sub_clipped = gpd.clip(world_utm, bbox_utm) 

# fig, ax = plt.subplots(figsize=(10, 10))
# # Plot the buffer_boundary
# world.plot(ax=ax, linewidth=2, color='gray')
# buffer_boundary_gdf.to_crs(projection).plot(ax=ax, facecolor='none', edgecolor='blue')
# bb_orig_gdf.plot(ax=ax, facecolor='none', edgecolor='red')
# plt.xlim(lon_min-3, lon_max+3)
# plt.ylim(lat_min-3, lat_max+3)
# plt.show()

len(bbox_latlon.exterior.coords)
len(buffer_boundary.exterior.coords)

# world_clipped_latlon = gpd.clip(world_utm, buffer_boundary).to_crs(projection)
world_clipped_latlon
########
cells_clipped = cells_gdf["geometry"].difference(world_clipped.geometry.union_all()).to_frame("geometry")
# cells_clipped = cells_gdf["geometry"].difference(world_clipped_latlon.geometry.union_all()).to_frame("geometry")
cell_colors = cells_clipped.area / (lat_step * lon_step)
# cell_colors = cells_clipped.to_crs({"proj": "cea"}).area / 46300.00000000001**2
cells_clipped['cell_colors'] = cell_colors
# ---> back to epsg lat/long
cells_latlon = cells_clipped.to_crs(projection)
cells_latlon_clipped = gpd.clip(cells_latlon, bb_orig_gdf)
cell_colors_clipped = cells_latlon_clipped.to_crs(utm_code).area / (lat_step * lon_step)
# cell_colors = cells_clipped.to_crs({"proj": "cea"}).area / 46300.00000000001**2
cells_latlon_clipped['cell_colors'] = cell_colors_clipped
########
from shapely.geometry import Point, LineString, shape
nasc_df = survey.input["acoustics"]["nasc_df"]
nasc_gdf = gpd.GeoDataFrame(data=nasc_df, geometry=gpd.points_from_xy(nasc_df["longitude"], nasc_df["latitude"]), crs=projection)
geo_df = nasc_gdf.groupby(["transect_num"])['geometry'].apply(lambda x: LineString(x.tolist())).to_frame("geometry").set_crs(projection)
custom_crs = '+proj=epsg:4326 +lat_ts=0 +lat_0=0 +lon_0=-180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs'
cells_latlon_clipped.to_crs(custom_crs).crs
########
import matplotlib.colors as colors
import matplotlib.cm as cm
cells_transformed = cells_latlon.to_crs(utm_code)
lims = cells_transformed.total_bounds

fig, ax = plt.subplots(figsize=(10, 10))
# cells_clipped.plot(ax=ax, column="cell_colors", edgecolor="black", cmap="viridis", legend=True)
# cells_clipped.plot.hexbin()
cells_latlon.to_crs(utm_code).plot(ax=ax, column="cell_colors", edgecolor="black", cmap="viridis", legend=False)
# cells_latlon.plot(ax=ax, column="cell_colors", edgecolor="black", cmap="viridis", legend=False)
# cells_latlon_clipped.plot(ax=ax, column="cell_colors", edgecolor="black", cmap="viridis", legend=False)
# cells_clipped.plot(ax=ax, column="cell_colors", edgecolor="black", cmap="viridis", legend=True)
# cells_gdf.plot(ax=ax, facecolor="none", edgecolor="black")
norm = colors.Normalize(vmin=cells_latlon["cell_colors"].min(), vmax=cells_latlon["cell_colors"].max())
cbar = plt.colorbar(cm.ScalarMappable(norm=norm, cmap="viridis"), ax=ax, orientation="horizontal", shrink=0.5)
cbar.set_label("Normalized grid area (50x50 nmi)", fontsize=12, labelpad=10, loc='center')  
cbar.ax.xaxis.set_label_position('top')
cbar.ax.xaxis.set_ticks_position('top')
geo_df.reset_index().to_crs(utm_code).plot(ax=ax, color="red")
# geo_df.reset_index().plot(ax=ax, color="red")
# plt.plot(ax=ax, nasc_df["longitude"], nasc_df["latitude"], color="red")
ax.margins(0.00, 0.00)
world_orig.to_crs(utm_code).plot(ax=ax, linewidth=1.2, color='gray', edgecolor="black")
# world_orig.plot(ax=ax, linewidth=1.2, color='gray', edgecolor="black")
# bb_orig_gdf.to_crs(utm_code).plot(ax=ax, facecolor='none', edgecolor='red')
plt.xlim(lims[0]*1.02, lims[2]*1.01)
# ax.set_yticks([4e6, 5e6, 6e6])
# ax.set_yticklabels(["4000", "5000", "6000"], fontsize=10)
plt.ylim(lims[1]*0.98, lims[3]*1.005)
ax.set_yticks([4e6, 5e6, 6e6])
ax.set_yticklabels(["4000", "5000", "6000"], fontsize=10)
plt.xlabel("Eastings (km)")
plt.ylabel("Northings (km)")
# plt.xlabel("Longitude (°E)")
# ax.set_xticks([-135, -130, -125, -120])
# plt.ylabel("Latitude (°N)")
ax.set_xticks([-600e3, -400e3, -200e3, 0, 200e3, 400e3, 600e3, 800e3])
ax.set_xticklabels(["-600", "-400", "-200", "0", "200", "400", "600", "800"], fontsize=10)
# Adding the colorbar title
# cax = fig.get_axes()[1]  # Assuming the colorbar is the second axis
# cax.set_ylabel("Normalized grid area (25x25 nmi)")  # Setting the title of the colorbar
plt.tight_layout()
plt.show()