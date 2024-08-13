from echopop.live.live_survey import LiveSurvey
from echopop.live.sql_methods import SQL
####################################################################################################
# TEST: Set up `LiveSurvey` object
# NOTE: General initialization parameter configuration
live_init_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/live_initialization_config.yml"
# NOTE: File configuration
live_file_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/live_survey_year_2019_config.yml"
# NOTE: Create object
realtime_survey = LiveSurvey(live_file_config_path, live_init_config_path, verbose=True)
# NOTE: String-representation via `LiveSurvey.__repr__`: 
# NOTE: Lists current files being processed and linked databases (WIP)
realtime_survey
####################################################################################################
# TEST: TRIGGER --> NEW ACOUSTIC DATA
# NOTE: Load new acoustic data (Either glob file search or `input_filenames Optional[List[str]]`)
realtime_survey.load_acoustic_data()
# NOTE: Process new acoustic data
# NOTE: This will update linked database tables
realtime_survey.process_acoustic_data()
# NOTE: Generate population estimates (or pass if there are no biological data)
# NOTE: `working_dataset = Literal["acoustic", "biology"]`
realtime_survey.estimate_population(working_dataset="acoustic")
# NOTE: String-representation via `LiveSurvey.__repr__`: 
# NOTE: Lists current files being processed and linked databases (WIP)
realtime_survey
####################################################################################################
# TEST: TRIGGER --> NEW BIOLOGY DATA
# NOTE: Load new biological data (Either glob file search or `input_filenames Optional[List[str]]`)
realtime_survey.load_biology_data()
# NOTE: Process new biological data
# NOTE: This will update linked database tables
realtime_survey.process_biology_data()
# NOTE: Generate population estimates (or pass if there are no acoustic data)
# NOTE: `working_dataset = Literal["acoustic", "biology"]`
realtime_survey.estimate_population(working_dataset="biology")
# NOTE: String-representation via `LiveSurvey.__repr__`: 
# NOTE: Lists current files being processed and linked databases (WIP)
realtime_survey
####################################################################################################
# TEST: `LiveSurvey` --[`files_processed`]--> `Echodataflow`
# NOTE: `LiveSurvey.meta` attribute
# ---- ACOUSTIC
realtime_survey.meta["provenance"]["acoustic_files"]
# ---- BIOLOGICAL
realtime_survey.meta["provenance"]["biology_files"]
# NOTE: SQL function query from database file [cumulative list]
# ---- ACOUSTIC
SQL(db_file=realtime_survey.config["database"]["acoustics"],
    command="select", table_name="files_processed")
# ---- BIOLOGICAL
SQL(db_file=realtime_survey.config["database"]["biology"],
    command="select", table_name="files_processed")
####################################################################################################
# TEST: `LiveSurvey` --[(key) SQL tables]--> Users
# !!! The SQL functions will fail if the tables have not yet been created/initialized
# ---- ACOUSTICS
# NOTE: Mean linear backscatter coefficient (`sigma_bs`) keyed for each haul and stratum
SQL(realtime_survey.config["database"]["acoustics"], "select", table_name="sigma_bs_mean_df")
# NOTE: Along-track acoustically-derived number/biomass densities and NASC 
SQL(realtime_survey.config["database"]["acoustics"], "select", table_name="survey_data_df")
# ---- BIOLOGICAL
# NOTE: Fitted (discretized) length-weight relationship
SQL(realtime_survey.config["database"]["biology"], "select", table_name="weight_fitted_df")
# NOTE: Quantized length-binned weights (summed)
SQL(realtime_survey.config["database"]["biology"], "select", table_name="length_weight_df")
# NOTE: Average weights per stratum
SQL(realtime_survey.config["database"]["biology"], "select", table_name="weight_stratum_df")