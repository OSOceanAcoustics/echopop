import numpy as np
import pandas as pd
from scipy.special import j1

from echopop.extensions.inversion.math import (
    generate_frequency_interval,
    length_average,
    orientation_average,
    wavenumber,
)
from echopop.extensions.inversion.scatterer import compute_Sv, compute_ts, uniformly_bent_cylinder
from echopop.extensions.inversion.scattering_models import pcdwba

####################################################################################################
# PARAMETERIZE DATASETS
# ---------------------
# Metadata
# ---------------------
metadata_df = pd.DataFrame(
    {
        "interval": np.repeat([1, 2, 3, 4, 5, 6], 3),
        "layer": np.tile([1, 2, 3], 6),
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
            -95.0,
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
            -89.0,
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
            -85.0,
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
            -72.0,
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
            -74.0,
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
    "number_density": {"distribution": "uniform", "initial": 3.0, "low": 1.0, "high": 1000.0},
    "theta_mean": {
        "distribution": "uniform",
        "initial": 10.0,
        "low": 0.0,
        "high": 90.0,
    },
    "theta_sd": {
        "distribution": "uniform",
        "initial": 20.0,
        "low": 0.0,
        "high": 90.0,
    },
    "length_mean": {
        "distribution": "uniform",
        "initial": 30e-3,
        "low": 8e-3,
        "high": 30e-3,
    },
    "length_sd_norm": {
        "distribution": "uniform",
        "initial": 0.15,
        "low": 0.05,
        "high": 0.15,
    },
    "g": {
        "distribution": "uniform",
        "initial": 1.015,
        "low": 1.015,
        "high": 1.060,
    },
    "h": {
        "distribution": "uniform",
        "initial": 1.020,
        "low": 1.015,
        "high": 1.060,
    },
    "radius_of_curvature_ratio": {
        "distribution": "uniform",
        "initial": 3.0,
        "low": 0.5,
        "high": 100.0,
    },
    "length_radius_ratio": {
        "distribution": "uniform",
        "initial": 18.2,
        "low": 14.0,
        "high": 20.0,
    },
    "taper_order": 10.0,
}
# -----------------------
# Optimization parameters
# -----------------------
optimization_parameters = {
    "max_iterations": 30,  # maximum number of iterations; == MaxIter
    "max_fun_evaluations": 200,  # maximum number of function evaluations == MaxFunEvals
    "fdgradient_max": 1.0,  # maximum change in variables for finite-differences gradients == DiffMaxChange
    "fgradient_min": 1e-3,  # minimum change in variables for finite-differences gradients == DiffMinChange
    "cost_fun_tolerance": 1e-3,  # cost-function tolerance for termination == TolFun
    "gradient_tolerance": 1e-3,  # gradient step tolerance for termination == TolX
}
# --------------------------------
# Simulation/processing parameters
# --------------------------------
processing_parameters = {
    "data_source": "measured", # options: "measured", "simulated"
    "inversion_source": "acoustics", # options: "acoustics", "biological"
    "center_frequencies": np.array([18e3, 38e3, 120e3]), # center frequencies (array), Hz
}
####################################################################################################
# Test params
water_sound_speed = 1500  # seawater sound speed, m s^-1
water_density = 1.0279  # seawater density, kg m^-3
L = 16.351316286953450e-3  # SL2 length, mm
L_std = 0.0900  # SL2 length standard deviation [ L_std / L]
length_radius_ratio = 18.2  # length-to-radius ratio
frequency_interval = 2e3  # Hz
n = 7  # number of frequencies/ka values
n_integration = 50  # minimum number of integration points
# theta = np.array([-0.9022, -0.8380, 2.8210, 2.8852]) # incident angle (broadside incidence = pi/2)
theta_mean = 43.9199  # mean orientation
g = 1.015  # density contrast
h = 1.020  # sound speed contrast
taper_order = 10  # shape tapering order
radius_of_curvature_ratio = 3.0  # radius of curvature ratio
center_frequencies = np.array([18e3, 38e3, 120e3])
# frequencies = np.array([12.9780e3, 14.9780e3, 16.9780e3, 18.9780e3, 20.9780e3, 22.9780e3, 24.9780e3]) # transmit frequencies, Hz
n_theta = 60  # number of orientation values to use
theta_sd = 35  # number of degrees offset for the incidence angle
ni_wavelen = 10  # number of sample points per wave length
theta_distribution = "gaussian"
length_mean = L
length_deviation = 0.09 * L  # SL2 length standard deviation, mm
length_bin_count = 100  # number of length bins for averaging
number_density = 4.720963086533185e3  # animal number density (animals m^-3)
####################################################################################################
# Compute acoustic property metrics
# ---------------------------------

# Normalize the length standard deviation
length_sd_norm = length_deviation / length_mean

# Generate frequency intervals centered on the central frequencies
frequencies = generate_frequency_interval(center_frequencies, length_sd_norm, frequency_interval)

# Compute the acoustic wavenumbers weighted by target size
# ---- Center frequencies
ka_center = wavenumber(center_frequencies, water_sound_speed) * length_mean / length_radius_ratio
# ---- Frequency intervals
# -------- Just wavenumber (`k`)
k = wavenumber(frequencies, water_sound_speed)
# -------- Now `ka`
ka = k * length_mean / length_radius_ratio

####################################################################################################
# Compute over a vector of angles (centered on 90 degrees)
theta_values = np.linspace(theta_mean - 3.1 * theta_sd, theta_mean + 3.1 * theta_sd, n_theta)
theta_radians = theta_values * np.pi / 180.0

length_values = np.linspace(
    length_mean - 3 * (L_std * length_mean),
    length_mean + 3 * (L_std * length_mean),
    length_bin_count,
)
####################################################################################################
# PCDWBA (TS modeling step)
# ------
fbs = compute_ts(
    taper_order,
    length_sd_norm,
    length_mean,
    length_radius_ratio,
    radius_of_curvature_ratio,
    theta_radians,
    k,
    ka,
    g,
    h,
    n_integration,
    ni_wavelen,
    model="pcdwba",
)
####################################################################################################
# Compute S_V
Sv_prediction = compute_Sv(
    number_density,
    theta_values,
    theta_mean,
    theta_sd,
    length_values,
    length_mean,
    length_deviation,
    fbs,
    ka,
    ka_center,
)
