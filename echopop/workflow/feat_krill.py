from pathlib import Path
from echopop.nwfsc_feat import ingest_sv
from echopop.inversion.inversion_matrix import InversionMatrix, estimate_population
from echopop.typing import InvParameters
from echopop.nwfsc_feat import utils

# ==================================================================================================
# ==================================================================================================
# DEFINE DATA ROOT DIRECTORY
# --------------------------
DATA_ROOT = Path("C:/Users/Brandyn/Documents/GitHub/EchoPro_data/echopop_inversion/")

# ==================================================================================================
# ==================================================================================================
# DATA INGESTION
# ==================================================================================================

sv_path = DATA_ROOT / "raw_files/2023/"
transect_pattern = r"x(\d+)"
impute_coordinates = True
center_frequencies = {18e3: {"min": -90., "max": -50.},
                      38e3: {"min": -90., "max": -50.},
                      70e3: {"min": -90., "max": -50.},
                      120e3: {"min": -90., "max": -50.},
                      200e3: {"min": -90., "max": -50.}}
method="transect"

sv_data, nasc_coordinates = ingest_sv.ingest_echoview_sv(sv_path=sv_path, 
                                                         center_frequencies=center_frequencies, 
                                                         transect_pattern=transect_pattern, 
                                                         aggregate_method=method,
                                                         impute_coordinates=True)

# ==================================================================================================
# ==================================================================================================
# DATA SUBSET
# ==================================================================================================

sv_data_sub = utils.apply_filters(sv_data, include_filter={"transect_num":[1, 2, 3]})

# ==================================================================================================
# ==================================================================================================
# DEFINE PARAMETERS
# ==================================================================================================

MODEL_PARAMETERS = {
    "number_density": {"value": 500., "min": 10., "max": 10000., "vary": True},
    "theta_mean": {"value": 10., "min": 0., "max": 90., "vary": True},
    "theta_sd": {"value": 20., "min": 0., "max": 90., "vary": False},
    "length_mean": {"value": 0.030, "min": 0.008, "max": 0.040, "vary": True},
    "length_sd_norm": {"value": 0.15, "min": 0.05, "max": 0.15, "vary": False},
    "g": {"value": 1.015, "min": 1.015, "max": 1.060, "vary": False},
    "h": {"value": 1.020, "min": 1.015, "max": 1.060, "vary":False},
    "radius_of_curvature_ratio": {"value": 3.0, "min": 0.5, "max": 100.0, "vary": False},
    "length_radius_ratio": {"value": 18.2, "min": 14.0, "max": 20.0, "vary": False},
}

# Model-specific settings, including distributions
MODEL_SETTINGS = {
    "type": "pcdwba",
    "taper_order": 10.,
    "frequency_interval": 2000.,
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
    "mc_realizations": 2,
    "scale_parameters": True, 
    "minimum_frequency_count": 2,
    "reference_frequency": 120e3,
}

OPTIMIZATION_KWARGS = {
    "max_nfev": 1000,
    "method": "least_squares",
    "loss": "huber",
    "xtol": 1e-8,
    "ftol": 1e-8,
    "gtol": 1e-8,
    "diff_step": 1e-7,
    "verbose": 2,
    "restart_strategy": {
        "max_attempts": 3,
        "Q_threshold": 3.,
        "scale": 0.10
    },
    "burnin": {
        "method": "nelder",
        "max_nfev": 100,
        "tol": 1e-3,
    }
}

def run_krill_inversion_workflow():
# ==================================================================================================
# ==================================================================================================
# FORMAT SCATTERING MODEL PARAMETERS
# ==================================================================================================

    scattering_parameters = InvParameters(MODEL_PARAMETERS)

# ==================================================================================================
# ==================================================================================================
# INITIALIZE
# ==================================================================================================

    INVERSION = InversionMatrix(sv_data_sub, SIMULATION_SETTINGS)

# ==================================================================================================
# ==================================================================================================
# BUILD AND FORMAT SCATTERING MODEL OPTIMIZERS
# ==================================================================================================

    INVERSION.build_scattering_model(scattering_parameters, MODEL_SETTINGS)

# ==================================================================================================
# ==================================================================================================
# BUILD AND FORMAT SCATTERING MODEL OPTIMIZERS
# ==================================================================================================
    inversion_results = INVERSION.invert(optimization_kwargs=OPTIMIZATION_KWARGS)

# ==================================================================================================
# ==================================================================================================
# CONVERT TO POPULATION ESTIMATES
# ==================================================================================================

    return estimate_population(inversion_results, 
                               nasc_coordinates, 
                               aggregate_method="transect",
                               density_sw=ENVIRONMENT["density_sw"],
                               reference_frequency=120e3)