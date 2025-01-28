import numpy as np
import pandas as pd
from scipy.special import j1

from echopop.extensions.inversion.math import (
    generate_frequency_interval,
    wavenumber,
)
from echopop.extensions.inversion.scatterer import compute_Sv, compute_ts
from echopop.extensions.inversion.optimize import (
    invert_population,
    normalize_parameters,
    prepare_optimization,
    optimize_scattering_model, 
    simulate_Sv,
    simulate_Sv_fit
)

####################################################################################################
# PARAMETERIZE DATASETS
# ---------------------
# Metadata
# ---------------------
metadata_df = pd.DataFrame(
    {
        "interval": np.repeat([1, 2, 3, 4, 5, 6], 3),
        "layer": np.tile([1, 2, 3], 6),
        "layer_thickness": 70.,
        "layer_volume": 3970.,
        "transect_num": np.repeat([1, 2], 9),
        "longitude": np.repeat([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], 3),
        "latitude": np.repeat([3.0, 5.0], 9),
    },
)
# ---------------------
# Data
# ---------------------
data_df = pd.DataFrame(
    {
        "interval": np.repeat([1, 2, 3, 4, 5, 6], 3),
        "layer": np.tile([1, 2, 3], 6),
        "Sv_mean_18": [
            -999.9,
            -999.9,
            -999.9,
            -999.9,
            -99.0,
            -999.9,
            -100.0,
            -102.8,
            -97.0,
            -999.9,
            -98.0,
            -999.9,
            -999.9,
            -999.9,
            -999.9,
            -999.9,
            -99.0,
            -999.9,
        ],
        "Sv_mean_38": [
            -999.9,
            -999.9,
            -999.9,
            -999.9,
            -99.0,
            -999.9,
            -92.0,
            -86.2,
            -90.0,
            -999.9,
            -88.0,
            -999.9,
            -999.9,
            -999.9,
            -999.9,
            -999.9,
            -90.0,
            -999.9,
        ],
        "Sv_mean_70": [
            -999.9,
            -999.9,
            -999.9,
            -999.9,
            -89.0,
            -999.9,
            -77.7,
            -83.0,
            -86.0,
            -999.9,
            -84.0,
            -999.9,
            -999.9,
            -999.9,
            -999.9,
            -999.9,
            -83.0,
            -999.9,
        ],
        "Sv_mean_120": [
            -999.9,
            -999.9,
            -999.9,
            -999.9,
            -77.0,
            -999.9,
            -76.0,
            -72.1,
            -77.0,
            -999.9,
            -76.0,
            -999.9,
            -999.9,
            -999.9,
            -999.9,
            -999.9,
            -74.0,
            -999.9,
        ],
        "Sv_mean_200": [
            -999.9,
            -999.9,
            -999.9,
            -999.9,
            -75.0,
            -999.9,
            -70.7,
            -73.0,
            -78.0,
            -999.9,
            -75.0,
            -999.9,
            -999.9,
            -999.9,
            -999.9,
            -999.9,
            -74.0,
            -999.9,
        ],
    }
)
# ---------------------------
# Scattering model parameters
# ---------------------------
scattering_parameters = {
    "number_density": {"distribution": "uniform", "initial": 3.0, "low": 1.0, "high": 1000.0, 
                       "vary": True},
    "theta_mean": {
        "distribution": "uniform",
        "initial": 10.0,
        "low": 0.0,
        "high": 90.0,
        "vary": True,
    },
    "theta_sd": {
        "distribution": "uniform",
        "initial": 20.0,
        "low": 0.0,
        "high": 90.0,
        "vary": False,
    },
    "length_mean": {
        "distribution": "uniform",
        "initial": 30e-3,
        "low": 8e-3,
        "high": 30e-3,
        "vary": True,
    },
    "length_sd_norm": {
        "distribution": "uniform",
        "initial": 0.15,
        "low": 0.05,
        "high": 0.15,
        "vary": False,
    },
    "g": {
        "distribution": "uniform",
        "initial": 1.015,
        "low": 1.015,
        "high": 1.060,
        "vary": False,
    },
    "h": {
        "distribution": "uniform",
        "initial": 1.020,
        "low": 1.015,
        "high": 1.060,
        "vary": False,
    },
    "radius_of_curvature_ratio": {
        "distribution": "uniform",
        "initial": 3.0,
        "low": 0.5,
        "high": 100.0,
        "vary": False,
    },
    "length_radius_ratio": {
        "distribution": "uniform",
        "initial": 18.2,
        "low": 14.0,
        "high": 20.0,
        "vary": False,
    },
}
# -----------------------
# Optimization parameters
# -----------------------
optimization_parameters = {
    # "max_iterations": 30,  # maximum number of iterations; == MaxIter
    # "max_fun_evaluations": 200,  # maximum number of function evaluations == MaxFunEvals
    # # "fdgradient_max": 1.0,  # maximum change in variables for finite-differences gradients == DiffMaxChange
    # # "fgradient_min": 1e-3,  # minimum change in variables for finite-differences gradients == DiffMinChange
    # "cost_fun_tolerance": 1e-3,  # cost-function tolerance for termination == TolFun
    # "gradient_tolerance": 1e-3,  # gradient step tolerance for termination == TolX
    "gtol": 1e-3,
    "ftol": 1e-3,
    "max_nfev": 30,
    # "jac": "3-point",
    # "maxiter": 30,
    "diff_step": 1e-6,
}
# --------------------------------
# Simulation/processing parameters
# --------------------------------
processing_parameters = {
    "data_source": "measured",  # options: "measured", "simulated"
    "ts_model": "pcdwba",
    "inversion_source": "acoustics",  # options: "acoustics", "biological"
    "center_frequencies": np.array([18e3, 38e3, 70e3, 120e3, 200e3]),  # center frequencies (array), Hz
    "taper_order": 10, # shape tapering order
    "frequency_interval": 2e3, # Hz
    "n_integration": 50, # number of integration points for the TS model
    "water_sound_speed": 1500, # seawater sound speed, m s^-1
    "water_density": 1.0279, # seawater density, kg m^-3
    "n_theta": 60, # number of orientation values to use
    "ni_wavelen": 10, # number of sample points per wave length
    "theta_distribution": "gaussian", # distribution type for theta distribution
    "length_distribution": "gaussian", # distribution type for length distribution
    "length_bin_count": 100, # number of length bins for averaging
}
####################################################################################################
# Test params
L = 16.351316286953450e-3  # SL2 length, mm
L_std = 0.0900  # SL2 length standard deviation [ L_std / L]
length_radius_ratio = 18.2  # length-to-radius ratio
# theta = np.array([-0.9022, -0.8380, 2.8210, 2.8852]) # incident angle (broadside incidence = pi/2)
theta_mean = 43.9199  # mean orientation
g = 1.015  # density contrast
h = 1.020  # sound speed contrast
radius_of_curvature_ratio = 3.0  # radius of curvature ratio
n_theta = 60  # number of orientation values to use
theta_sd = 35  # number of degrees offset for the incidence angle
length_mean = L
length_deviation = 0.09 * L  # SL2 length standard deviation, mm
number_density = 4.720963086533185e3  # animal number density (animals m^-3)
Sv_measured = data_df.iloc[7, 2:].to_numpy()
####################################################################################################
# Normalize scattering parameters
scattering_params_norm = normalize_parameters(scattering_parameters)

# Prepare parameters and minimization procedure
parameters, minimizer = prepare_optimization(scattering_parameters, 
                                             processing_parameters, 
                                             Sv_measured)

# Run inversion
# ---- Compute the initial fit (for comparison -- pre-optimized values)
initial_fit = simulate_Sv_fit(parameters, Sv_measured, processing_parameters)

# Optimize
best_fit_params, best_fit_Sv, fit_error = optimize_scattering_model(minimizer, 
                                                                    Sv_measured,
                                                                    optimization_parameters, 
                                                                    processing_parameters)
# ----- Compute delta
delta = initial_fit - fit_error


# Extract the best-fit parameter values
optimized_scattering_params = best_fit_params.params.valuesdict()
# from lmfit import fit_report
# print(fit_report(optimized_scattering_params))

# Generate population estimate array
population_estimate = invert_population(optimized_scattering_params, 
                                        metadata_df.loc[7], 
                                        processing_parameters)