####################################################################################################
# TEST BACKWARDS COMPATIBILITY WITH OTHER SURVEY YEARS
# LINKED ISSUE: https://github.com/OSOceanAcoustics/echopop/issues/307
####################################################################################################
# SUCCESSFUL TESTS
# 2017, 2019, 2021, 2023
####################################################################################################
# FAILURE NOTES
# 2015: Strings with spaces in some entries produce errors
####################################################################################################
# CURRENT SURVEY YEAR BEING TESTED: 2017
SURVEY_YEAR = 2013
####################################################################################################
# Initialize year-specific survey parameters
from echopop.survey import Survey
from echopop.extensions import generate_reports
import pandas as pd
import glob

# Initialization configuration
init_config_path = f"C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config_\
{SURVEY_YEAR}.yml"

# Filepath/dataset configuration
survey_year_config_path = f"C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_\
{SURVEY_YEAR}_config.yml"
####################################################################################################
# Run
# ---- Create object
survey = Survey(init_config_path, survey_year_config_path)
# ---- Load survey data
survey.load_survey_data()
# ---- Load acoustic data
survey.load_acoustic_data(ingest_exports="echoview")
# survey.load_acoustic_data()
# ---- Transect analysis
survey.transect_analysis()
# ---- Stratified analysis
survey.stratified_analysis(bootstrap_ci_method="percentile")
# ---- Fit variogram
survey.fit_variogram()
# ---- Kriging analysis
survey.kriging_analysis()
# ---- Stratified analysis
survey.stratified_analysis(dataset="kriging", bootstrap_ci_method="percentile")
# ---- Report generation
survey.generate_reports()
####################################################################################################
# Get report outputs
# ---- Define directory
output_dir = f"C:/Users/Brandyn/Documents/GitHub/EchoPro_data/echopo_backtest_years/Outputs/\
Historical Outputs (KS Stratification with aged data)/without extrapolation/{SURVEY_YEAR}"
# ---- Get files
files = glob.glob(output_dir + "/EchoPro_un-kriged_output*")
# ---- Filter files
filtered_files = [f for f in files if f.endswith("_1.xlsx")]
# ---- Read
echopro_data = pd.read_excel(filtered_files[0])
####################################################################################################
# Compare
echopop_data = survey.analysis["transect"]["acoustics"]["adult_transect_df"].copy()
# ---- Calculate: NASC
echopro_data["NASC"].sum() - echopop_data["nasc"].sum()
# ---- Calculate: biomass
(echopro_data["wgt_total"].sum() - echopop_data["biomass"].sum()) * 1e-6
