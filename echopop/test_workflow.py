from echopop.live.live_survey import LiveSurvey
from echopop.live.sql_methods import SQL

# Set up `LiveSurvey` object
live_init_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/live_initialization_config.yml"
live_file_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/live_survey_year_2019_config.yml"
realtime_survey = LiveSurvey(live_file_config_path, live_init_config_path, verbose=True)
realtime_survey
####################################################################################################
# TEST: ACOUSTICS
####################################################################################################
# NOTE: LOAD DATA
realtime_survey.load_acoustic_data()
realtime_survey
SQL(realtime_survey.config["database"]["acoustics"], "select", table_name="files_read")
SQL(realtime_survey.config["database"]["acoustics"], "select", table_name="files_processed")
SQL(realtime_survey.config["database"]["acoustics"], "map")
realtime_survey.config["database"]
realtime_survey.meta["provenance"]
# NOTE: INITIAL PROCESSING [JUST ACOUSTIC]
# ! ERRORS OUT WHEN NUMBER OF FILES == 1
realtime_survey.process_acoustic_data()
realtime_survey.estimate_population(working_dataset="acoustic")
self = realtime_survey
SQL(realtime_survey.config["database"]["acoustics"], "select", table_name="survey_data_df")
####################################################################################################
# TEST: BIOLOGY
####################################################################################################
# NOTE: LOAD DATA
realtime_survey.load_biology_data()
# NOTE: INITIAL PROCESSING [JUST BIOLOGY]
realtime_survey.process_biology_data()
realtime_survey.estimate_population(working_dataset="biology")
SQL(realtime_survey.config["database"]["biology"], "select", table_name="files_read")
SQL(realtime_survey.config["database"]["biology"], "select", table_name="files_processed")
SQL(realtime_survey.config["database"]["biology"], "map")
####################################################################################################
# TEST: POPULATION ESTIMATES
####################################################################################################
# NOTE: Acoustic / biological data converge here to derive population estimates 
# TODO: Add argument that indicates what the new datasets and what data need to be pulled in
# TODO: ARGUMENT {working_dataset: Literal["acoustic", "biology"]}
# ! SQL ARGUMENT STRINGS FAIL ON > 1000 ENTRIES (250 ROWS)
realtime_survey.estimate_population(working_dataset="biology")
realtime_survey.estimate_population(working_dataset="acoustic")
####################################################################################################
# TEST: GET DATA
####################################################################################################
SQL(realtime_survey.config["database"]["acoustics"], "select", table_name="survey_data_df")
SQL(realtime_survey.config["database"]["acoustics"], "select", table_name="sigma_bs_mean_df")
SQL(realtime_survey.config["database"]["biology"], "select", table_name="strata_summary_df")