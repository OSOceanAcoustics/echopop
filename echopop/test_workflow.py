from echopop.live.live_survey import LiveSurvey
from echopop.live.sql_methods import reset_db_files

live_init_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/live_initialization_config.yml"
live_file_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/live_survey_year_2019_config.yml"

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
realtime_survey.estimate_population()