from pathlib import Path

from echopop import utils
from echopop.ingest import sv
from echopop.inversion import InversionMatrix, InvParameters, estimate_population

# ==================================================================================================
# ==================================================================================================
# DEFINE DATA ROOT DIRECTORY
# --------------------------
DATA_ROOT = Path("C:/Data/EchopopData/echopop_inversion/")

# ==================================================================================================
# ==================================================================================================
# DATA INGESTION
# ==================================================================================================

SV_PATH = DATA_ROOT / "raw_files/2023/"
TRANSECT_PATTERN = r"x(\d+)"
IMPUTE_BAD_COORDINATES = True
CENTER_FREQUENCIES = {
    18e3: {"min": -90.0, "max": -50.0},
    38e3: {"min": -90.0, "max": -50.0},
    70e3: {"min": -90.0, "max": -50.0},
    120e3: {"min": -90.0, "max": -50.0},
    200e3: {"min": -90.0, "max": -50.0},
}
PROCESSING_METHOD = "transect"

sv_data, nasc_coordinates = sv.ingest_echoview_sv(
    sv_path=SV_PATH,
    center_frequencies=CENTER_FREQUENCIES,
    transect_pattern=TRANSECT_PATTERN,
    aggregate_method=PROCESSING_METHOD,
    impute_coordinates=IMPUTE_BAD_COORDINATES,
)

# ==================================================================================================
# ==================================================================================================
# DATA SUBSET
# ==================================================================================================

sv_data_sub = utils.apply_filters(sv_data, include_filter={"transect_num": [1, 2, 3]})

# ==================================================================================================
# ==================================================================================================
# DEFINE PARAMETERS
# ==================================================================================================

MODEL_PARAMETERS = {
    "number_density": {"value": 500.0, "min": 10.0, "max": 10000.0, "vary": True},
    "theta_mean": {"value": 10.0, "min": 0.0, "max": 90.0, "vary": True},
    "theta_sd": {"value": 20.0, "min": 0.0, "max": 90.0, "vary": False},
    "length_mean": {"value": 0.030, "min": 0.008, "max": 0.040, "vary": True},
    "length_sd_norm": {"value": 0.15, "min": 0.05, "max": 0.15, "vary": False},
    "g": {"value": 1.015, "min": 1.015, "max": 1.060, "vary": False},
    "h": {"value": 1.020, "min": 1.015, "max": 1.060, "vary": False},
    "radius_of_curvature_ratio": {"value": 3.0, "min": 0.5, "max": 100.0, "vary": False},
    "length_radius_ratio": {"value": 18.2, "min": 14.0, "max": 20.0, "vary": False},
}

# Model-specific settings, including distributions
MODEL_SETTINGS = {
    "type": "pcdwba",
    "taper_order": 10.0,
    "frequency_interval": 2000.0,
    "n_integration": 50,
    "n_wavelength": 10,
    "orientation_distribution": {"family": "gaussian", "bins": 60},
    "length_distribution": {"family": "gaussian", "bins": 100},
}

ENVIRONMENT = {
    "sound_speed_sw": 1500.0,
    "density_sw": 1026.9,
}

SIMULATION_SETTINGS = {
    "environment": ENVIRONMENT,
    "monte_carlo": True,
    "mc_realizations": 20,
    "minimum_frequency_count": 2,
}

OPTIMIZATION_KWARGS = {
    "max_nfev": 2000,
    # "method": "lbfgsb",
    "method": "least_squares",
    "loss": "cauchy",
    "xtol": 1e-8,
    # "tol": 1e-8,
    # "options": {"maxiter": 100, "disp": False, "eps": 1e-6},
    "ftol": 1e-6,
    "gtol": 1e-6,
    "diff_step": 1e-6,
    "verbose": 1,
    "restart_strategy": {"max_attempts": 2, "Q_threshold": 3.0, "scale": 0.2},
    "burnin": {
        "method": "nelder",
        "max_nfev": 200,
        "tol": 1e-4,
    },
}


def run_krill_inversion_workflow():
    # ==================================================================================================
    # ==================================================================================================
    # FORMAT SCATTERING MODEL PARAMETERS
    # ==================================================================================================

    scattering_parameters = InvParameters(MODEL_PARAMETERS)
    # ---- SCALE THE PARAMETERS FOR IMPROVED NUMERICAL STABILITY
    scattering_parameters.scale()

    # ==================================================================================================
    # ==================================================================================================
    # INITIALIZE
    # ==================================================================================================

    inversion_krill = InversionMatrix(sv_data_sub, SIMULATION_SETTINGS)

    # ==================================================================================================
    # ==================================================================================================
    # BUILD AND FORMAT SCATTERING MODEL OPTIMIZERS
    # ==================================================================================================

    inversion_krill.build_scattering_model(scattering_parameters, MODEL_SETTINGS)

    # ==================================================================================================
    # ==================================================================================================
    # BUILD AND FORMAT SCATTERING MODEL OPTIMIZERS
    # ==================================================================================================
    df_inversion_results = inversion_krill.invert(optimization_kwargs=OPTIMIZATION_KWARGS)

    # ==================================================================================================
    # ==================================================================================================
    # CONVERT TO POPULATION ESTIMATES
    # ==================================================================================================

    return estimate_population(
        df_inversion_results,
        nasc_coordinates,
        aggregate_method="transect",
        density_sw=ENVIRONMENT["density_sw"],
        reference_frequency=120e3,
    )
