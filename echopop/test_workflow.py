from echopop.live.live_survey import LiveSurvey
from echopop.live.sql_methods import reset_db_files
from echopop.live.sql_methods import query_processed_files
from echopop.live.live_acoustics import preprocess_acoustic_data, compute_nasc
from echopop.live.live_biology import preprocess_biology_data
from echopop.live.sql_methods import SQL, SQL_COMMANDS, query_processed_files, format_sql_columns, sql_group_update, sql_data_exchange

from echopop.live.live_core import(
    LIVE_DATA_STRUCTURE,
)
from echopop.live.live_biology import (
    bin_length_data,
    compute_average_weights,
    compute_sigma_bs,
    length_bin_counts,
    length_weight_regression,
    number_proportions,
    length_bin_weights,
    preprocess_biology_data,
    weight_proportions
)
from echopop.live import live_data_processing as eldp
from echopop.live import live_data_loading as eldl

live_init_config_path = "C:/Users/15052/Documents/GitHub/echopop/config_files/live_initialization_config.yml"
live_file_config_path = "C:/Users/15052/Documents/GitHub/echopop/config_files/live_survey_year_2019_config.yml"

realtime_survey = LiveSurvey(live_file_config_path, live_init_config_path)

####################################################################################################
# TEST: ACOUSTICS
####################################################################################################
# NOTE: Reset database file for utility purposes
reset_db_files(realtime_survey.config)

# NOTE: LOAD DATA
realtime_survey.load_acoustic_data()
# NOTE: INITIAL PROCESSING [JUST ACOUSTIC]
realtime_survey.process_acoustic_data()
realtime_survey.input
####################################################################################################
# TEST: BIOLOGY
####################################################################################################
# NOTE: Reset database file for utility purposes
reset_db_files(realtime_survey.config)

# NOTE: LOAD DATA
realtime_survey.load_biology_data()
realtime_survey.input
# NOTE: INITIAL PROCESSING [JUST BIOLOGY]
realtime_survey.process_biology_data()
realtime_survey.input
####################################################################################################
# TEST: POPULATION ESTIMATES
####################################################################################################
# NOTE: Acoustic / biological data converge here to derive population estimates 
# TODO: Add argument that indicates what the new datasets and what data need to be pulled in
# TODO: ARGUMENT {working_dataset: Literal["acoustic", "biology"]}
realtime_survey.estimate_population()