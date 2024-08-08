from echopop.live.live_survey import LiveSurvey
from echopop.live.sql_methods import SQL
from echopop.live.live_biology import (
    bin_length_data,
    compute_average_weights,
    compute_sigma_bs,
    length_bin_counts,
    length_bin_weights,
    length_weight_regression,
    number_proportions,    
    preprocess_biology_data,
    weight_proportions
)
# Set up `LiveSurvey` object
live_init_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/live_initialization_config.yml"
live_file_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/live_survey_year_2019_config.yml"
realtime_survey = LiveSurvey(live_file_config_path, live_init_config_path, verbose=True)
realtime_survey
####################################################################################################
# TEST: ACOUSTICS
# Actual flow:
realtime_survey.load_acoustic_data() #`input_filenames` = Optional[List[str]]
realtime_survey.process_acoustic_data()
realtime_survey.estimate_population(working_dataset="acoustic")
amo = SQL(realtime_survey.config["database"]["acoustics"], "select", table_name="survey_data_df")
amo[amo.nasc > 0]
SQL(realtime_survey.config["database"]["biology"], "select", table_name="weight_fitted_df")

realtime_survey.load_biology_data() #`input_filenames` = Optional[List[str]]
realtime_survey.process_biology_data()
realtime_survey.estimate_population(working_dataset="biology")
SQL(realtime_survey.config["database"]["acoustics"], "select", table_name="sigma_bs_mean_df")
SQL(realtime_survey.config["database"]["biology"], "select", table_name="strata_summary_df")
tbl = SQL(realtime_survey.config["database"]["biology"], "select", table_name="length_weight_df")
tbl[tbl.weight > 0]
tbl.weight.sum()
# NOTE: Pulling successfully processed filenames
# ! This dictionary key name will change
realtime_survey.meta["provenance"][f"{working_dataset}_files"]
SQL(realtime_survey.config["database"]["acoustics"], "select", table_name="files_processed")
####################################################################################################
# NOTE: LOAD DATA
table_df[table_df.weight > 0]
realtime_survey.load_acoustic_data()
realtime_survey
SQL(realtime_survey.config["database"]["acoustics"], "select", table_name="files_read")
SQL(realtime_survey.config["database"]["acoustics"], "select", table_name="files_processed")
out = SQL(realtime_survey.config["database"]["acoustics"], "map")
SQL(realtime_survey.config["database"]["biology"], "select", table_name="length_weight_df")
_ = SQL(realtime_survey.config["database"]["biology"], "drop", table_name="length_weight_df")
realtime_survey.config["database"]
realtime_survey.meta["provenance"]
# NOTE: INITIAL PROCESSING [JUST ACOUSTIC]
# ! ERRORS OUT WHEN NUMBER OF FILES == 1
realtime_survey.process_acoustic_data()
# ! sqlalchemy.exc.OperationalError: (sqlite3.OperationalError) near "18": syntax error
realtime_survey.estimate_population(working_dataset="acoustic")
self = realtime_survey
SQL(realtime_survey.config["database"]["acoustics"], "select", table_name="survey_data_df")
####################################################################################################
# TEST: BIOLOGY
# Actual flow
realtime_survey.load_biology_data() #`input_filenames` = Optional[List[str]]
realtime_survey.process_biology_data()
realtime_survey.estimate_population(working_dataset="biology")
self = realtime_survey
biology_unprocessed = self.input["biology"]
specimen_data = biology_unprocessed["specimen_df"]
length_data = biology_unprocessed["length_df"]
biology_dict = self.input["biology"]
file_configuration = self.config
strata_df = self.input["spatial"]["strata"]
from echopop.live.live_acoustics import average_sigma_bs, compute_nasc, estimate_echometrics, integrate_nasc
from echopop.live.sql_methods import sql_group_update
from echopop.live.live_biology import summarize_strata
import numpy as np; import pandas as pd
echometrics: bool = True
acoustic_data_df = self.input["acoustics"]["prc_nasc_df"].copy()
spatial_column = ["stratum"]
# acoustic_data_df_copy = acoustic_data_df.copy()
acoustic_data_df_copy.groupby(["longitude", "latitude", "ping_time"] + spatial_column).apply(integrate_nasc, echometrics, include_groups=False)

acoustic_data_df.groupby(["longitude", "latitude", "ping_time"] + spatial_column).apply(integrate_nasc, echometrics, include_groups=False).droplevel(-1).reset_index()
acoustic_data_df_copy.groupby(["longitude", "latitude", "ping_time"] + spatial_column).apply(integrate_nasc, echometrics, include_groups=False)
acoustic_data_df.groupby(["longitude", "latitude", "ping_time"] + spatial_column).apply(integrate_nasc, echometrics, include_groups=False)

acoustic_data_df_copy.groupby(["longitude", "latitude", "ping_time"] + spatial_column, as_index=True, group_keys=True).apply(integrate_nasc, echometrics, include_groups=False)
acoustic_data_df.groupby(["longitude", "latitude", "ping_time"] + spatial_column, as_index=True, group_keys=True).apply(integrate_nasc, echometrics, include_groups=False)
acoustic_data_df.groupby(["longitude", "latitude", "ping_time"] + spatial_column).apply(lambda g: integrate_nasc(g, echometrics)).reset_index(drop=True)
dd.index.get_level_values(-1)
cc.index.get_level_values(-1)
(
    acoustic_data_df
    .groupby(['longitude', 'latitude', 'ping_time', 'source'] + spatial_column)
    .apply(lambda g: integrate_nasc(g, echometrics=True), include_groups=False)
    .reset_index()
    # .rename_axis(None, axis=0)  # Remove any unwanted hierarchical index
)
(acoustic_data_df.groupby(['longitude', 'latitude', 'ping_time', 'source', 'stratum'])
              .apply(integrate_nasc, echometrics=True)
              .reset_index())
acoustic_data_df = acoustic_data_df[acoustic_data_df.distance == 0.0]
acoustic_data_df = acoustic_data_df_copy[acoustic_data_df_copy.distance==0.0]
pd.Series(nasc_dict).index
pd.DataFrame.from_dict(nasc_dict, orient="columns") 
pd.DataFrame(nasc_dict, index=[0])
SQL(realtime_survey.config["database"]["acoustics"], "select", table_name="sigma_bs_mean_df")
(
    acoustic_data_df.groupby(["longitude", "latitude", "ping_time", "source"] + spatial_column, 
                             observed=False)
                             .apply(lambda df: integrate_nasc(df, echometrics)).reset_index()
)

(
    acoustic_data_df.groupby(["longitude", "latitude", "ping_time", "source", "stratum"], observed=False)
    .apply(lambda df: integrate_nasc(df, echometrics=True), include_groups=False)
    .reset_index()
)

result_df = acoustic_data_df.groupby(['longitude', 'latitude', 'ping_time', 'source', 'stratum']) \
    .apply(lambda g: integrate_nasc(g, echometrics=True), include_groups=False) \
    .reset_index()

print(acoustic_data_df.columns)
print(acoustic_data_df_copy.columns)
print(acoustic_data_df.dtypes)
print(acoustic_data_df_copy.dtypes)
# Inspect DataFrame before groupby
print(acoustic_data_df.head())
print(acoustic_data_df_copy.head())

print(acoustic_data_df["longitude"].unique())
print(acoustic_data_df_copy["longitude"].unique())

print(acoustic_data_df["latitude"].unique())
print(acoustic_data_df_copy["latitude"].unique())

print("Grouped original index levels:", grouped_original.size().index.names)
print("Grouped reset index levels:", grouped_reset.size().index.names)

print(acoustic_data_df["ping_time"].unique())
print(acoustic_data_df_copy["ping_time"].unique())

print(acoustic_data_df["source"].unique())
print(acoustic_data_df_copy["source"].unique())

grouped_original = acoustic_data_df_copy.groupby(["longitude", "latitude", "ping_time", "source"] + spatial_column, observed=False)
grouped_reset = acoustic_data_df.groupby(["longitude", "latitude", "ping_time", "source"] + spatial_column, observed=False)
print(acoustic_data_df.index)
print(acoustic_data_df_copy.index)
grouped_original.index
print(grouped_original.size())
print(grouped_reset.size())
sql_group_update(acoustic_db, sigma_bs_df, table_name="sigma_bs_mean_df")
####################################################################################################
# NOTE: LOAD DATA
realtime_survey.load_biology_data()
# NOTE: INITIAL PROCESSING [JUST BIOLOGY]
realtime_survey.process_biology_data()
realtime_survey.estimate_population(working_dataset="biology")
SQL(realtime_survey.config["database"]["biology"], "select", table_name="files_read")
SQL(realtime_survey.config["database"]["biology"], "select", table_name="files_processed")
SQL(realtime_survey.config["database"]["biology"], "map")
SQL(realtime_survey.config["database"]["acoustics"], "select", table_name="sigma_bs_mean_df")
SQL(realtime_survey.config["database"]["acoustics"], "drop", table_name="sigma_bs_mean_df")
####################################################################################################
# TEST: POPULATION ESTIMATES
####################################################################################################
# NOTE: Acoustic / biological data converge here to derive population estimates 
# TODO: Add argument that indicates what the new datasets and what data need to be pulled in
# TODO: ARGUMENT {working_dataset: Literal["acoustic", "biology"]}
# ! SQL ARGUMENT STRINGS FAIL ON > 1000 ENTRIES (250 ROWS)
realtime_survey.estimate_population(working_dataset="biology")
realtime_survey.estimate_population(working_dataset="acoustic")
self = realtime_survey
acoustic_dict = self.input["acoustics"]
strata_df = self.input["spatial"]["strata"]
file_configuration = self.config
from echopop.live.sql_methods import SQL, sql_group_update
from echopop.live.live_biology import summarize_strata

db_file = acoustic_db
dataframe=nasc_biology
table_name="survey_data_df"
columns=["number_density", "biomass_density"]
unique_columns = ["stratum", "longitude", "latitude", "ping_time"]
####################################################################################################
# TEST: GET DATA
####################################################################################################
SQL(realtime_survey.config["database"]["acoustics"], "select", table_name="survey_data_df")
SQL(realtime_survey.config["database"]["acoustics"], "select", table_name="sigma_bs_mean_df")
SQL(realtime_survey.config["database"]["biology"], "select", table_name="strata_summary_df")

SQL(acoustic_db, "drop", table_name="sigma_bs_mean_df")

#####
# NOTE: Below are hypothetical visualizations
# 
survey_data = SQL(realtime_survey.config["database"]["acoustics"], "select", 
                  table_name="survey_data_df")

import matplotlib.pyplot as plt
import numpy as np

survey_data.loc[0, "nasc"] = 1e3

plt.plot(survey_data["longitude"], survey_data["latitude"])
plt.scatter(survey_data["longitude"], survey_data["latitude"], s=survey_data["nasc"])
plt.show()

SQL(realtime_survey.config["database"]["biology"], "map")
# ! NEED TO ENSURE THAT TABLE FOR LENGTH/WEIGHT HISTOGRAM IS AVAILABLE
SQL(realtime_survey.config["database"]["biology"], "select", table_name="length_weight_df")
SQL(realtime_survey.config["database"]["biology"], "select", table_name="weight_fitted_df")
SQL(realtime_survey.config["database"]["biology"], "select", table_name="weight_stratum_df")
realtime_survey.input["spatial"]["strata"]
#
SQL(realtime_survey.config["database"]["acoustics"], "select", table_name="sigma_bs_mean_df")
SQL(realtime_survey.config["database"]["biology"], "select", table_name="strata_summary_df")