####################################################################################################
# TEST BACKWARDS COMPATIBILITY WITH OTHER SURVEY YEARS
# LINKED ISSUE: https://github.com/OSOceanAcoustics/echopop/issues/307
####################################################################################################
# SUCCESSFUL TESTS
# 2015, 2017, 2019, 2021, 2023
####################################################################################################
# FAILURE NOTES
# 2012: Transect-region mapping
# Initialize year-specific survey parameters
####################################################################################################
# Import libraries
from echopop.survey import Survey
from echopop.extensions import generate_reports
import pandas as pd
from pathlib import Path
import glob
import json
import os

####################################################################################################
# CURRENT SURVEY YEAR BEING TESTED: 2019
####################################################################################################
# Define current survey year
SURVEY_YEAR = 2017

# Initialization configuration
init_config_path = f"C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_\
config_{SURVEY_YEAR}.yml"

# Filepath/dataset configuration
survey_year_config_path = f"C:/Users/Brandyn/Documents/GitHub/echopop/config_files\
/survey_year_{SURVEY_YEAR}_config.yml"

# Load json settings
# ---- File open
with open(Path(os.getcwd() + "\\echopop\\compatibility_parameters_test.json").as_posix()) as f:
    json_dict = json.load(f)
# ---- Load
parameters = json_dict[f"{SURVEY_YEAR}"]

####################################################################################################
# Run
# ---- Create object
survey = Survey(init_config_path, survey_year_config_path)
# ---- Load acoustic data (Echoview export ingestion)
survey.load_acoustic_data(ingest_exports="echoview", 
                          read_transect_region_file=parameters["read_transect_region_file"], 
                          write_transect_region_file=parameters["write_transect_region_file"])
# ---- Load acoustic data (Already-defined file)
if parameters["default_acoustics"]:
    survey =  Survey(init_config_path, survey_year_config_path)
    survey.load_acoustic_data()
else:
    Survey(init_config_path, survey_year_config_path).load_acoustic_data()
# ---- Load survey data
survey.load_survey_data()
# ---- Initial transect analysis test
survey.transect_analysis()
# ---- Counter
counter = 1
# ---- Iterate across multiple strata types
for stratum in parameters["strata_types"]:
    # ---- Iterate across different age-1 exclusion definitions
    for excl in parameters["exclude_age1"]:           
            # ---- Transect analysis
            survey.transect_analysis(exclude_age1=excl, stratum=stratum, verbose=False)
            # ---- Stratified analysis (transect analysis)
            survey.stratified_analysis(bootstrap_ci_method=parameters["bootstrap_ci_method"],
                                       transect_replicates=parameters["transect_replicates"], 
                                       verbose=False)
            # ---- Fit variogram
            survey.fit_variogram(verbose=False)
            # ---- Iterate across different extrapolation schema
            for extrap in parameters["extrapolate"]: 
                # ---- Kriging analysis (no variogram fitting)
                survey.kriging_analysis(extrapolate=extrap, verbose=False)
                # ---- Apply best-fit variogram
                survey.kriging_analysis(best_fit_variogram=True, extrapolate=extrap, verbose=False)
                # ---- Stratified analysis (kriging analysis)
                survey.stratified_analysis(dataset="kriging",
                                           bootstrap_ci_method=parameters["bootstrap_ci_method"],
                                           transect_replicates=parameters["transect_replicates"],
                                           verbose=False)
                # ---- Test reports  
                survey.generate_reports(reports=["aged_length_haul_counts",
                                                "kriging_input",
                                                "kriged_length_age_abundance",
                                                "kriged_length_age_biomass",
                                                "kriged_mesh_results",
                                                "total_length_haul_counts",
                                                "transect_length_age_abundance",
                                                "transect_length_age_biomass",
                                                "transect_population_results"])
                # ---- Print out success
                print(
                    f"""            
                    Year: {SURVEY_YEAR} success [{counter}/8 configurations]!
                        Stratum: {stratum}
                        Age-1 fish excluded: {excl}
                        Extrapolated: {extrap}            
                    """
                )
                # ---- Advance counter
                counter += 1
            
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
