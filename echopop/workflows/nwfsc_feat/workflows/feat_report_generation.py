import os
import pickle
from pathlib import Path

import pandas as pd

from echopop.workflows.nwfsc_feat.reporter import Reporter

# ==================================================================================================
# ==================================================================================================
# DEFINE DATA ROOT DIRECTORY
# --------------------------
DATA_ROOT = Path("C:/Data/EchopopData/echopop_2019")
# --------------------------
# SAVE DIRECTORY FOR REPORTS
# --------------------------
SAVE_DIRECTORY = DATA_ROOT / "reports_test"

# ==================================================================================================
# ==================================================================================================
# DEFINE DEMO DATA ROOT DIRECTORY
# -------------------------------

# Estabalish workflow directory
WORKFLOW_DIR = Path(os.getcwd()) / "echopop/workflow"
# ---- Validate existence
WORKFLOW_DIR.exists()

# Demo folder
DEMO_DIR = WORKFLOW_DIR / "demo"
# ---- Validate existence
DEMO_DIR.exists()

# Assign sub-folder for files
FILES_DIR = DEMO_DIR / "files"

# Assign sub-folder for output figures
FIGURES_DIR = DEMO_DIR / "figures"

# ==================================================================================================
# Load the pickled dataframes
# ---------------------------
try:
    # NASC - transect data
    df_nasc_noage1_prt = pd.read_pickle(FILES_DIR / "df_nasc_no_age1_prt.pkl")
    # Mesh - kriged data
    df_kriged_results = pd.read_pickle(FILES_DIR / "df_kriged_results.pkl")
    # Abundance table - kriged data
    df_kriged_abundance_table = pd.read_pickle(FILES_DIR / "df_kriged_abundance_table.pkl")
    # Biomass table - kriged data
    df_kriged_biomass_table = pd.read_pickle(FILES_DIR / "df_kriged_biomass_table.pkl")
    # Abundance table - transect data
    with open(FILES_DIR / "dict_transect_abundance_table.pkl", "rb") as f:
        dict_transect_abundance_table = pickle.load(f)
    # Biomass table - transect data
    with open(FILES_DIR / "dict_transect_biomass_table.pkl", "rb") as f:
        dict_transect_biomass_table = pickle.load(f)
    # Biomass table - transect aged-only data
    df_transect_aged_biomass_table = pd.read_pickle(
        FILES_DIR / "df_transect_aged_biomass_table.pkl"
    )
    # Abundance tables - kriged data
    with open(FILES_DIR / "dict_kriged_abundance_table.pkl", "rb") as f:
        dict_kriged_abundance_table = pickle.load(f)
    # Biomass tables - kriged data
    with open(FILES_DIR / "dict_kriged_biomass_table.pkl", "rb") as f:
        dict_kriged_biomass_table = pickle.load(f)
    # Biohaul data
    with open(FILES_DIR / "biohaul_data.pkl", "rb") as f:
        biohaul_data = pickle.load(f)
    # Binned weights
    with open(FILES_DIR / "dict_df_weight_distr.pkl", "rb") as f:
        dict_df_weight_distr = pickle.load(f)
    # Linear scattering coefficient
    sigma_bs_strata = pd.read_pickle(FILES_DIR / "stratum_sigma_bs.pkl")
    # Stratified weights
    df_averaged_weight = pd.read_pickle(FILES_DIR / "df_averaged_weight.pkl")
    # Verbose validation upon success
    print("Pickled demo DataFrames and Dictionaries successfully 'unpickled'.")
except Exception as e:
    raise e from None

# ==================================================================================================
# ==================================================================================================
# REPORT GENERATION
# -----------------

# Initialize report generator
reporter = Reporter(SAVE_DIRECTORY)

####################################################################################################
# Aged-length haul counts report
# ------------------------------

reporter.aged_length_haul_counts_report(
    filename="aged_length_haul_counts.xlsx",
    sheetnames={"male": "Sheet1", "female": "Sheet2", "all": "Sheet3"},
    bio_data=biohaul_data["specimen"].dropna(subset=["age", "length", "weight"]),
)

####################################################################################################
# Total length haul counts report
# -------------------------------

reporter.total_length_haul_counts_report(
    filename="total_length_haul_counts.xlsx",
    sheetnames={"male": "Sheet1", "female": "Sheet2", "all": "Sheet3"},
    bio_data=biohaul_data,
)

####################################################################################################
# Kriged aged biomass mesh report
# -------------------------------

# All values
reporter.kriged_aged_biomass_mesh_report(
    filename="kriged_aged_biomass_mesh_full.xlsx",
    sheetnames={"all": "Sheet1", "male": "Sheet2", "female": "Sheet3"},
    kriged_data=df_kriged_results,
    weight_data=dict_df_weight_distr["aged"],
    kriged_stratum_link={"geostratum_ks": "stratum_ks"},
)

# Nonzero values
reporter.kriged_aged_biomass_mesh_report(
    filename="kriged_aged_biomass_mesh_nonzero.xlsx",
    sheetnames={"all": "Sheet1", "male": "Sheet2", "female": "Sheet3"},
    kriged_data=df_kriged_results[df_kriged_results["biomass"] > 0.0],
    weight_data=dict_df_weight_distr["aged"],
    kriged_stratum_link={"geostratum_ks": "stratum_ks"},
)

####################################################################################################
# Kriged mesh results report
# --------------------------

# All values
reporter.kriged_mesh_results_report(
    filename="kriged_biomass_mesh_full.xlsx",
    sheetname="Sheet1",
    kriged_data=df_kriged_results,
    kriged_stratum="geostratum_ks",
    kriged_variable="biomass",
    sigma_bs_data=sigma_bs_strata,
    sigma_bs_stratum="stratum_ks",
)

# Nonzero values
reporter.kriged_mesh_results_report(
    filename="kriged_biomass_mesh_nonzero.xlsx",
    sheetname="Sheet1",
    kriged_data=df_kriged_results[df_kriged_results["biomass"] > 0.0],
    kriged_stratum="geostratum_ks",
    kriged_variable="biomass",
    sigma_bs_data=sigma_bs_strata,
    sigma_bs_stratum="stratum_ks",
)

####################################################################################################
# Kriged length-age abundance report
# ----------------------------------

reporter.kriged_length_age_abundance_report(
    filename="kriged_length_age_abundance_report.xlsx",
    sheetnames={"male": "Sheet1", "female": "Sheet2", "all": "Sheet3"},
    datatables=dict_kriged_abundance_table,
)

####################################################################################################
# Kriged length-age biomass report
# --------------------------------

reporter.kriged_length_age_biomass_report(
    filename="kriged_length_age_biomass_report.xlsx",
    sheetnames={"male": "Sheet1", "female": "Sheet2", "all": "Sheet3"},
    datatables=dict_kriged_biomass_table,
)

####################################################################################################
# Kriging input report
# --------------------

reporter.kriging_input_report(
    filename="kriging_input_report.xlsx",
    sheetname="Sheet1",
    transect_data=df_nasc_noage1_prt,
)

####################################################################################################
# Transect length-age abundance report
# ------------------------------------

reporter.transect_length_age_abundance_report(
    filename="transect_length_age_abundance_report.xlsx",
    sheetnames={"male": "Sheet1", "female": "Sheet2", "all": "Sheet3"},
    datatables=dict_transect_abundance_table,
)

####################################################################################################
# Kriged length-age biomass report
# --------------------------------

reporter.transect_length_age_biomass_report(
    filename="kriged_length_age_biomass_report.xlsx",
    sheetnames={"male": "Sheet1", "female": "Sheet2", "all": "Sheet3"},
    datatable=dict_transect_biomass_table["aged"],
)

####################################################################################################
# Transect aged biomass report
# ----------------------------

# Full values
reporter.transect_aged_biomass_report(
    filename="transect_aged_biomass_report_full.xlsx",
    sheetnames={"all": "Sheet1", "male": "Sheet2", "female": "Sheet3"},
    transect_data=df_nasc_noage1_prt,
    weight_data=dict_df_weight_distr["aged"],
)

# Nonzero values
reporter.transect_aged_biomass_report(
    filename="transect_aged_biomass_report_nonzero.xlsx",
    sheetnames={"all": "Sheet1", "male": "Sheet2", "female": "Sheet3"},
    transect_data=df_nasc_noage1_prt[df_nasc_noage1_prt["biomass"] > 0.0],
    weight_data=dict_df_weight_distr["aged"],
)

####################################################################################################
# Transect population results report
# ----------------------------------

# Full values
reporter.transect_population_results_report(
    filename="transect_population_results_full.xlsx",
    sheetname="Sheet1",
    transect_data=df_nasc_noage1_prt,
    weight_strata_data=df_averaged_weight,
    sigma_bs_stratum=sigma_bs_strata,
    stratum_name="stratum_ks",
)

# Nonzero values
reporter.transect_population_results_report(
    filename="transect_population_results_nonzero.xlsx",
    sheetname="Sheet1",
    transect_data=df_nasc_noage1_prt[df_nasc_noage1_prt["nasc"] > 0.0],
    weight_strata_data=df_averaged_weight,
    sigma_bs_stratum=sigma_bs_strata,
    stratum_name="stratum_ks",
)
