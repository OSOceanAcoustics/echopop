import numpy as np
import pandas as pd
from scipy.special import j1
from echopop.extensions.inversion.scattering_models import pcdwba
from echopop.extensions.inversion.scatterer import uniformly_bent_cylinder, compute_Sv
from echopop.extensions.inversion.scatterer import compute_Sv
from echopop.extensions.inversion.math import wavenumber, reflection_coefficient, orientation_average, length_average
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
# ---------------------
# Scattering model parameters
# ---------------------
{
    "g": {
        "distribution": "uniform",
        "low": 1.015,
        "high": 1.060,
    },
    "h": {
        "distribution": "uniform",
        "low": 1.015,
        "high": 1.060,
    },
}
####################################################################################################
# Test params
water_sound_speed = 1500 # seawater sound speed, m s^-1
water_density = 1.0279 # seawater density, kg m^-3
L = 16.351316286953450e-3 # SL2 length, mm
L_std = 0.0900 # SL2 length standard deviation [ L_std / L]
L_a = 18.2 # length-to-radius ratio
frequency_interval = 2e3 # Hz
n = 7 # number of frequencies/ka values
n_integration = 50 # minimum number of integration points
# theta = np.array([-0.9022, -0.8380, 2.8210, 2.8852]) # incident angle (broadside incidence = pi/2)
theta_mean = 43.9199 # mean orientation
g = 1.015 # density contrast
h = 1.020 # sound speed contrast
taper_order = 10 # shape tapering order
rho_L = 3.0 # radius of curvature ratio
center_frequency = np.array([18e3, 38e3, 120e3])
# frequencies = np.array([12.9780e3, 14.9780e3, 16.9780e3, 18.9780e3, 20.9780e3, 22.9780e3, 24.9780e3]) # transmit frequencies, Hz
n_theta = 60 # number of orientation values to use
theta_sd = 35 # number of degrees offset for the incidence angle
ni_wavelen = 10 # number of sample points per wave length
theta_distribution = "gaussian"
length_mean = L
length_deviation = 0.09 * L # SL2 length standard deviation, mm
length_bin_count = 100 # number of length bins for averaging
number_density = 4.720963086533185e3 # animal number density (animals m^-3)
####################################################################################################
# Compute acoustic property metrics
# ---------------------------------
# Compute frequency values that will be averaged together
# ---- Normalize the length standard deviation
length_sd_norm = length_deviation / length_mean
# ---- Generate new frequency values
def generate_frequency_interval(frequency: np.ndarray[float], 
                                length_sd_norm: float, 
                                frequency_interval: float):
    frequency_lst = [
        np.arange(freq * (1-3.1*length_sd_norm), 
                  freq * (1+3.1*length_sd_norm) + frequency_interval,
                  frequency_interval)
        for freq in frequency
    ]

    # Find the maximum length of the generated arrays
    max_length = max(len(arr) for arr in frequency_lst)
    
    # Create a padded 2D array with NaN for shorter arrays
    padded_results = np.full((len(frequency_lst), max_length), np.nan)
    for i, arr in enumerate(frequency_lst):
        padded_results[i, :len(arr)] = arr
    
    return padded_results

frequencies = generate_frequency_interval(center_frequency, length_sd_norm, frequency_interval)

# Compute acoustic wavenumber
k = wavenumber(frequencies, water_sound_speed)
# Convert length to m
# L_m = L * 1e-3
# Initial ka
ka = k*length_mean/L_a
ka_center = wavenumber(center_frequency, water_sound_speed)*length_mean/L_a

kLmax = np.nanmax(k*length_mean, axis=1)*(1+3.1*length_sd_norm)
# Calculate the number of integration points
n_int = np.where(kLmax < n_integration, n_integration, np.ceil(kLmax*ni_wavelen/(2*np.pi))).astype(int)

n_int = np.array([50, 30, 40])

####################################################################################################
# Build position vector
taper, gamma_tilt, beta_tilt, r_pos, dr_pos = uniformly_bent_cylinder(n_int, rho_L, taper_order)
####################################################################################################
# Compute over a vector of angles (centered on 90 degrees)
theta_values = np.linspace(theta_mean - 3.1*theta_sd, theta_mean + 3.1*theta_sd, n_theta) 
theta_radians = theta_values * np.pi / 180.0

length_values = np.linspace(length_mean - 3*(L_std*length_mean), length_mean + 3*(L_std*length_mean), length_bin_count)
####################################################################################################
# PCDWBA
# ------
length_radius_ratio = L_a
theta = theta_radians
fbs = pcdwba(taper, gamma_tilt, beta_tilt, r_pos, dr_pos, length_mean, L_a, g, h, ka, theta_radians)
####################################################################################################
# Compute S_V
Sv_prediction = compute_Sv(number_density, theta_values, theta_mean, theta_sd, length_values, 
                           length_mean, length_deviation, fbs, ka, ka_center)

np.arctan(dz/dx).shape
np.concatenate([np.arctan(dz/dx), np.arctan(dz[-1]/dx[-1])], axis=1)
A = np.arctan(dz/dx)

dx[:, 0]
np.concatenate([np.sqrt(dx[:, 0]*dx[:, 0] + dz[:, 0]*dz[:, 0]).reshape(-1, 1), np.sqrt(dx*dx + dz*dz)], axis=1)

np.concatenate([A, np.ones((A.shape[0],1),dtype=A.dtype)], axis=1)

[ka[i, :n_k[i]].reshape(-1, 1) * taper[i] / h for i in range(len(n_k))]

ka[0, :n_k[0]].reshape(-1, 1) * taper[0] / h

(ka[0].reshape(-1, 1) * taper[0] / h).reshape
beta_tilt.shape
gamma_tilt.shape
(gamma_tilt.reshape(-1, 1) - theta).shape
def process_and_multiply(ka, taper):
    # Initialize an empty 3D array to store the results
    result = np.full((ka.shape[0], taper.shape[1], taper.shape[1]), np.nan)
    
    # Iterate over each row in ka to copy valid values into result
    for i, row in enumerate(ka):
        valid_values = row[~np.isnan(row)]  # Extract valid (non-NaN) values
        
        # Fill the result with valid values, padded to match taper's shape
        result[i, :len(valid_values), :len(valid_values)] = valid_values[:, np.newaxis]  # Reshape to 2D
        
    # Multiply result with taper across all valid elements
    result *= taper[:, np.newaxis]
    
    return result

# Example usage with ka and taper of different shapes
process_and_multiply(ka, taper)
np.concatenate([alpha_tilt, new_column.reshape(-1, 1)], axis=1)
np.ndindex(ka.shape)
ka.reshape(-1, 1).squeeze() * taper.reshape(1, -1).squeeze()