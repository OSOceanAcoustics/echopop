# from echopop.live.live_survey import LiveSurvey
# from echopop.live.sql_methods import SQL
# import echopop.live.live_visualizer as elv
# from pathlib import Path
# from echopop.live import live_data_processing as eldp
# from echopop.live import live_data_loading as eldl
# from echopop.live.live_core import(
#     LIVE_DATA_STRUCTURE, LIVE_INPUT_FILE_CONFIG_MAP
# )
# import boto3
# from botocore.exceptions import NoCredentialsError, ClientError
# import pandas as pd
# import numpy as np
# from echopop.live.sql_methods import SQL, sql_data_exchange, get_table_key_names,
# sql_group_update, query_processed_files, sql_update_strata_summary
# from echopop.live.live_spatial_methods import apply_spatial_definitions
# from echopop.live.live_acoustics import average_sigma_bs, compute_nasc
# from echopop.live.live_biology import compute_sigma_bs
# from echopop.acoustics import ts_length_regression, to_dB, to_linear
# from echopop.utils.operations import group_interpolator_creator
# from functools import reduce
# from echopop.live.live_data_loading import filter_filenames, read_biology_csv

# ##################################################################################################
# # TEST: Set up `LiveSurvey` object
# # NOTE: General initialization parameter configuration
# live_init_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/live_initializat
# ion_config.yml"
# # NOTE: File configuration
# live_file_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/live_survey_yea
# r_2019_config.yml"
# # NOTE: Create object
# realtime_survey = LiveSurvey(live_init_config_path, live_file_config_path, verbose=True)
# realtime_survey = LiveSurvey(live_file_config_path, live_init_config_path, verbose=True)

# # NOTE: String-representation via `LiveSurvey.__repr__`:
# # NOTE: Lists current files being processed and linked databases (WIP)
# self = realtime_survey
# file_configuration = self.config

# input_filenames = ["202407_003_operation_info.csv", "202407_22500_003_lf.csv",
# "202407_22500_003_spec.csv", "202407_003_catch_perc.csv"]
# realtime_survey.config["input_directories"]["biology"]["directory"] =
# "s3://sh2407-upload/data/Echopop-biology"

# survey_data = SQL("C:/Users/Brandyn/Downloads/acoustics.db", "select",
# table_name="survey_data_df")


# del realtime_survey.config["data_root_dir"]
# self = realtime_survey

# # realtime_survey.config["storage_options"] = aws_credentials
# realtime_survey = LiveSurvey(live_init_config_path, live_file_config_path, verbose=True)
# realtime_survey.load_biology_data(input_filenames=input_filenames)
# realtime_survey.input["biology"]
# def is_s3_path(path):
#     """Check if a path is an S3 path."""
#     return path.startswith("s3://")

# dataset_directory = realtime_survey.config["input_directories"]["biology"]["directory"]
# s3_path = dataset_directory
# is_s3_path(dataset_directory)

# cloud_credentials = aws_credentials
# cloud_credentials = {}
# def validate_s3_path(s3_path: str, cloud_credentials: dict):
#     """Check if (parts of) S3 path exists."""

#     # Redundant validation that S3 object validation is appropriate
#     if not is_s3_path(s3_path):
#         raise ValueError("The path is not an S3 path.")

#     # Validate credentials
#     if not all([True if param in cloud_credentials.keys() else False
#                 for param in ["key", "secret"]]):
#         # ---- Find missing credentials
#         missing_creds = set(["key", "secret"]) - set(cloud_credentials)
#         # ---- Format into string
#         missing_creds_str = ", ".join(["'{}'".format(x.replace("'", "''")) for x in
# missing_creds])
#         # ---- Raise Error
#         raise PermissionError(
#             f"Required S3 credentials missing: {missing_creds_str}."
#         )

#     # Remove the s3:// prefix
#     s3_path_reduced = s3_path[len("s3://"):]

#     # Split into bucket and key
#     parts = s3_path_reduced.split("/", 1)
#     if len(parts) < 2:
#         raise ValueError(f"Invalid S3 path format for '{s3_path}'.")

#     # Get bucket name and directory keys
#     bucket_name, directory = parts

#     # Initialize the S3 client
#     s3_client = boto3.client("s3",
#                              aws_access_key_id=cloud_credentials["key"],
#                              aws_secret_access_key=cloud_credentials["secret"])

#     # Check if the bucket exists
#     try:
#         s3_client.head_bucket(Bucket=bucket_name)
#     except ClientError as e:
#         raise FileNotFoundError(
#             f"S3 bucket '{bucket_name}' does not exist or you do not have access."
#         )

#     # Check if the S3 directory exists
#     try:
#         # ---- Ping a response from the bucket
#         response = s3_client.list_objects_v2(Bucket=bucket_name, Prefix=directory, MaxKeys=1)
#         # ---- Check for `Contents`
#         if "Contents" not in response:
#             raise FileNotFoundError(f"S3 path '{s3_path}' does not exist.")
#     except ClientError as e:
#         # --- Raise Error and propagate it upwards
#         raise e

# validate_s3_path(s3_path, cloud_credentials)

# import pandas as pd

# self = realtime_survey
# biology_files = self.meta["provenance"]["biology_files_read"]
# file_configuration = self.config
# dataset = "biology"

# # Get the dataset file settings
# file_settings = file_configuration["input_directories"][dataset]

# def construct_directorypath(file_configuration: dict, file_settings: dict):
#     """Construct the root directory path."""

#     # Get the general root_directory, if present
#     if "data_root_dir" in file_configuration:
#         root_directory = file_configuration["data_root_dir"]
#     else:
#         root_directory = ""

#     # Get the local directory (or this may be the root directory depending on the config)
#     data_directory = file_settings["directory"]

#     # Return the directory path
#     if root_directory != "":
#         return "/".join([root_directory, data_directory])
#     else:
#         return data_directory

# directory_path = construct_directorypath(file_configuration, file_settings)

# def validate_local_path(directory_path: str):

#     # Validate filepath
#     # ---- Error evaluation (if applicable)
#     if not Path(directory_path).exists():
#         raise FileNotFoundError(
#             f"The acoustic data directory [{directory_path}] does not exist."
#         )

#     # Validate that files even exist
#     # ---- List available files of target extension
#     data_files = list(directory_path.glob(f"*{'.'+file_settings['extension']}"))
#     # ---- Error evaluation (if applicable)
#     if not data_files:
#         raise FileNotFoundError(
#             f"No `*.{file_settings['extension']}` files found in [{directory_path}]!"
#         )


# # Get the biology data file settings
# file_settings = file_configuration["input_directories"]["biology"]

# # Get the file-specific settings, datatypes, columns, etc.
# # ---- Get defined columns and datatypes from `LIVE_INPUT_FILE_CONFIG_MAP`
# biology_config_map = LIVE_INPUT_FILE_CONFIG_MAP["biology"]
# # ---- Extract the expected file name ID's
# biology_file_ids = file_settings["file_name_formats"]
# # ---- Extract all of the file ids
# biology_config_ids = list(biology_file_ids.keys())
# # ---- Initialize the dictionary that will define this key in the `input` attribute
# biology_output = {f"{key}_df": pd.DataFrame() for key in biology_config_ids}


# # Initialize a session with AWS credentials
# s3_client = boto3.client(
#     's3',
#     aws_access_key_id=aws_credentials["key"],
#     aws_secret_access_key=aws_credentials["secret"]
# )
# response = s3_client.list_buckets()
# buckets = response.get('Buckets', [])
# for bucket in buckets:
#     print(f"Bucket Name: {bucket['Name']}")
# s3_client.head_bucket(Bucket="sh2407-upload")
# realtime_survey.load_biology_data(pandas_kwargs=aws_credentials, input_filenames=input_filenames)
# realtime_survey.config["ship_id"]
# grid_data = SQL(realtime_survey.config["database"]["grid"], "select", table_name="grid_df")
# grid_data[grid_data.abundance > 0]
# bucket = boto3.client("s3", region_name=None)
# bucket.head_bucket(Bucket=realtime_survey.config["input_directories"]["biology"]["directory"]
# +"/")
# bucket.list_objects_v2(Bucket=realtime_survey.config["input_directories"]["biology"]["directory"],
# Prefix=path, MaxKeys=1)
# #################################################################################################
# # TEST: TRIGGER --> NEW ACOUSTIC DATA
# # NOTE: Load new acoustic data (Either glob file search or `input_filenames Optional[List[str]]`)
# realtime_survey.load_acoustic_data()
# # NOTE: Process new acoustic data
# # NOTE: This will update linked database tables
# realtime_survey.process_acoustic_data()
# # NOTE: Generate population estimates (or pass if there are no biological data)
# # NOTE: `working_dataset = Literal["acoustic", "biology"]`
# realtime_survey.estimate_population(working_dataset="acoustic")
# # NOTE: String-representation via `LiveSurvey.__repr__`:
# # NOTE: Lists current files being processed and linked databases (WIP)
# realtime_survey.input["acoustics"]
# ##################################################################################################
# # TEST: TRIGGER --> NEW BIOLOGY DATA
# # NOTE: Load new biological data (Either glob file search or `input_filenames Optional[List[str]]`
# realtime_survey.load_biology_data()
# len(realtime_survey.meta["provenance"]["biology_files_checkpoint1"])
# realtime_survey.meta["provenance"]["biology_files_checkpoint3"]
# # NOTE: Process new biological data
# # NOTE: This will update linked database tables
# realtime_survey.process_biology_data()
# # NOTE: Generate population estimates (or pass if there are no acoustic data)
# # NOTE: `working_dataset = Literal["acoustic", "biology"]`
# realtime_survey.estimate_population(working_dataset="biology")
# # NOTE: String-representation via `LiveSurvey.__repr__`:
# # NOTE: Lists current files being processed and linked databases (WIP)
# realtime_survey
# ##################################################################################################
# # TEST: `LiveSurvey` --[`files_processed`]--> `Echodataflow`
# # NOTE: `LiveSurvey.meta` attribute
# # ---- ACOUSTIC
# realtime_survey.meta["provenance"]["acoustic_files"]
# # ---- BIOLOGICAL
# realtime_survey.meta["provenance"]["biology_files"]
# # NOTE: SQL function query from database file [cumulative list]
# # ---- ACOUSTIC
# SQL(db_file=realtime_survey.config["database"]["acoustics"],
#     command="select", table_name="files_processed")
# dat = SQL(db_file=realtime_survey.config["database"]["acoustics"],command="select",
# table_name="files_processed")
# # ---- BIOLOGICAL
# SQL(db_file=realtime_survey.config["database"]["biology"],command="select",
# table_name="files_processed")
# dat.loc[0:, "filepath"][105]
# ##################################################################################################
# # TEST: `LiveSurvey` --[(key) SQL tables]--> Users
# # !!! The SQL functions will fail if the tables have not yet been created/initialized
# # ---- ACOUSTICS
# # NOTE: Mean linear backscatter coefficient (`sigma_bs`) keyed for each haul and stratum
# SQL(realtime_survey.config["database"]["biology"], "select", table_name="sigma_bs_mean_df")
# SQL(realtime_survey.config["database"]["biology"], "select", table_name="specimen_df")
# .latitude.max()
# realtime_survey.input["spatial"]["strata"]
# # NOTE: Along-track acoustically-derived number/biomass densities and NASC
# SQL(realtime_survey.config["database"]["acoustics"], "select", table_name="survey_data_df")
# # ---- BIOLOGICAL
# # NOTE: Fitted (discretized) length-weight relationship
# SQL(realtime_survey.config["database"]["biology"], "select", table_name="weight_fitted_df")
# # NOTE: Quantized length-binned weights (summed)
# SQL(realtime_survey.config["database"]["biology"], "select", table_name="length_weight_df")
# # NOTE: Average weights per stratum
# SQL(realtime_survey.config["database"]["biology"], "select", table_name="weight_stratum_df")
# # NOTE: Stratum summary tables
# SQL(realtime_survey.config["database"]["biology"], "select", table_name="strata_summary_df")
# ##################################################################################################
# # FROM THE `LiveSurvey` object !
# # ---- Convert to a Panel
# import panel as pn
# # ---- Either have the db file already called in as a `pandas.DataFrame`, or query the table
# survey_data_db = Path(realtime_survey.config["database"]["acoustics"])
# # grid_db = Path(realtime_survey.config["database"]["grid"])
# grid_db = Path("C:/Users/Brandyn/Downloads/grid.db")
# dat = SQL(realtime_survey.config["database"]["acoustics"], "select", table_name="survey_data_df")
# dat
# dat1 = SQL(grid_db, "select", table_name="grid_df")
# SQL(realtime_survey.config["database"]["biology"], "select", table_name="sigma_bs_mean_df")

# sql_cmd = "SELECT * FROM sigma_bs_mean_df ORDER BY stratum, haul_num, species_id"
# # Create the engine
# engine = create_engine(f"sqlite:///{"C:/Users/Brandyn/Downloads/biology.db"}")
# # Create the SQL database connection and send the script
# with engine.connect() as connection:
#     table = connection.execute(text(sql_cmd))

# data = table.fetchall()
# dd = pd.DataFrame(data, columns=table.keys()).loc[0:1, :]
# dd = dd[["stratum", "haul_num", "species_id", "sigma_bs", "sigma_bs_count", "sigma_bs_sum", "id"]]
# dd.loc[:, "id"] = pd.Series([f"{(4,4,4)}", f"{(5,5,5)}"])
# SQL("C:/Users/Brandyn/Downloads/biology.db", "insert", table_name="sigma_bs_mean_df",
# dataframe=dd)
# SQL("C:/Users/Brandyn/Downloads/biology.db", "map")
# SQL(biology_db, "drop", table_name="sigma_bs_mean_df")
# SQL(biology_db, "select", table_name="sigma_bs_mean_df")
# dd.loc[:, "haul_num"] = pd.Series([101, 103])
# dd = dd[["species_id", "haul_num", "id", "stratum", "sigma_bs", "sigma_bs_count", "sigma_bs_sum"]]
# SQL(biology_db, "insert", table_name="sigma_bs_mean_df", dataframe=dd, id_columns=key_list+["id"])
# SQL(biology_db, "select", table_name="sigma_bs_mean_df")
# import numpy as np; import pandas as pd
# SQL("C:/Users/Brandyn/Downloads/biology.db", "select", table_name="length_weight_df")
# sigma_bs_df = SQL("C:/Users/Brandyn/Downloads/biology.db", "select",
# table_name="sigma_bs_mean_df")
# table_df = SQL(realtime_survey.config["database"]["biology"], "select",
# table_name="sigma_bs_mean_df")
# sigma_bs_df = table_df
# # ---- Check the table keys
# table_keys = np.unique(table_df["id"]).tolist()
# # ---- Get unique values
# current_keys = np.unique(sigma_bs_df["id"]).tolist()
# # ---- Get INSERTION keys
# insertion_keys = list(set(current_keys).difference(set(table_keys)))
# # ---- Get UPDATE keys
# update_keys = list(set(current_keys).intersection(set(table_keys)))
# insertion_df = sigma_bs_df[sigma_bs_df["id"].isin(insertion_keys)]
# insertion_df.loc[0, "species_id"] = 22500
# insertion_df.loc[0, "stratum"] = 5
# insertion_df.loc[0, "haul_num"] = 100
# insertion_df.loc[0, "sigma_bs"] = 1e-10
# insertion_df.loc[0, "sigma_bs_count"] = 100
# insertion_df.loc[0, "sigma_bs_sum"] = 1e10 * 100
# insertion_df.loc[0, "id"] = f"{(1,1,1)}"
# SQL(realtime_survey.config["database"]["biology"], "insert", table_name="sigma_bs_mean_df",
#     dataframe=insertion_df)
# SQL(realtime_survey.config["database"]["acoustics"], "select", table_name="sigma_bs_mean_df")
# survey_data = SQL(realtime_survey.config["database"]["acoustics"], "select",
# table_name="survey_data_df")
# dat1[dat1.abundance > 0]
# dat[dat.number_density > 0]
# coast_db = grid_db
# biology_db = Path(realtime_survey.config["database"]["biology"])
# projection = realtime_survey.config["geospatial"]["projection"]
# # NOTE: PLOTS
# # Ensure Panel is initialized
# pn.extension()
# # ---- Helper function
# def plt_to_pn(fig):
#     # Convert to a panel object
#     panel = pn.panel(fig)
#     # Display
#     panel.show() # OR panel.servable() if you want to serve it in a Panel server
# # ---- PLOT GRID
# fig = elv.plot_livesurvey_grid(grid_db, projection, coast_db)
# fig.show()
# plt_to_pn(fig)
# # ---- PLOT TRACK
# from echopop.live.live_visualizer import plot_livesurvey_track
# fig1 = plot_livesurvey_track(survey_data, projection, coast_db)
# fig1.show()
# plt_to_pn(fig1)
# # ---- PLOT DISTRIBUTIONS
# weight_table = SQL(biology_db, "select",
#                    table_name="length_weight_df")
# stratum_table = SQL(biology_db, "select",
#                     table_name="strata_summary_df")
# specimen_table = SQL(biology_db, "select",
#                      table_name="specimen_data_df")
# length_table = SQL(biology_db, "select",
#                    table_name="length_df")
# fig2 = elv.plot_livesurvey_distributions(weight_table, stratum_table, specimen_table,
# length_table)
# plt_to_pn(fig2)
# ### MULTIPANEL
# panel0 = pn.panel(fig, name='Gridded population estimates')
# panel1 = pn.panel(fig1, name='Alongtrack population estimates')
# panel2 = pn.panel(fig2, name='Length and weight distributions')

# def serve_panels():
#     # Create links to each panel
#     home = pn.Column(
#         pn.pane.Markdown("# Main Page"),
#         pn.pane.Markdown("[Gridded population estimates](gridded_population_estimates)",
# sizing_mode="stretch_width"),
#         pn.pane.Markdown("[Alongtrack population estimates](alongtrack_population_estimates)",
# sizing_mode="stretch_width"),
#         pn.pane.Markdown("[Length and weight distributions](length_weight_distributions)",
# sizing_mode="stretch_width")
#     )

#     # Serve the home page and individual panels
#     pn.serve({
#         'Main Page': home,
#         'gridded_population_estimates': panel0,
#         'alongtrack_population_estimates': panel1,
#         'length_weight_distributions': panel2
#     },  show=True)
# # Run the function to serve panels
# serve_panels()
