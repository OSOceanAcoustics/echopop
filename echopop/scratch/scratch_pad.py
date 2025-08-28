import pandas as pd
import numpy as np
from typing import Any, Dict, List, Union, Optional, Literal
from pydantic import ValidationError
from pathlib import Path
from typing import Any, Dict
import numpy as np
import numpy.typing as npt
import pandas as pd
from lmfit import Parameters
from echopop.inversion.inversion_base import InversionBase
from echopop.inversion.operations import impute_missing_sigma_bs
from echopop import validators as val
from echopop.nwfsc_feat import ingest_nasc, utils
from echopop import acoustics
from echopop.core.echoview import (
    ECHOVIEW_DATABASE_EXPORT_FILESET,
    ECHOVIEW_EXPORT_ROW_SORT,
    ECHOVIEW_TO_ECHOPOP,
)
from echopop.scratch.scratch_tests import *
from scipy.special import j1
import hashlib
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
transect_pattern: Optional[str] = r"x(\d+)"
impute_coordinates: bool = True
center_frequencies = {18e3: {"min": -90., "max": -50.},
                      38e3: {"min": -90., "max": -50.},
                      70e3: {"min": -90., "max": -50.},
                      120e3: {"min": -90., "max": -50.},
                      200e3: {"min": -90., "max": -50.}}
method="transects"

sv_data = ingest_echoview_sv(sv_path=sv_path, 
                             center_frequencies=center_frequencies, 
                             transect_pattern=transect_pattern, 
                             aggregate_method=method,
                             impute_coordinates=True)

# !!! EQUIVALENT TO `reduce_dataset` from inversion extension PR
utils.apply_filters(sv_data, include_filter={"transect_num": np.linspace(1, 30, 30)})

####################################################################################################
MODEL_PARAMETERS = {
    "number_density": {"value": 500., "min": 10., "max": 10000., "vary": True},
    "theta_mean": {"value": 10., "min": 0., "max": 90., "vary": True},
    "theta_sd": {"value": 20., "min": 0., "max": 90., "vary": False},
    "length_mean": {"value": 0.030, "min": 0.008, "max": 0.040, "vary": True},
    "length_sd_norm": {"value": 0.15, "min": 0.05, "max": 0.15, "vary": False},
    "g": {"value": 1.015, "min": 1.015, "max": 1.060, "vary": False},
    "h": {"value": 1.020, "min": 1.015, "max": 1.060, "vary":False},
    "radius_of_curvature_ratio": {"value": 3.0, "min": 0.5, "max": 100.0, "vary": True},
    "length_radius_ratio": {"value": 18.2, "min": 14.0, "max": 20.0, "vary": True},
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
    "n_realizations": 10,
    "scale_parameters": True, 
    "minimum_frequency_count": 2,
    "reference_frequency": 120e3,
}

OPTIMIZATION_KWARGS = {
    "max_nfev": 500,
    "gtol": 1e-3,
    "ftol": 1e-3,
    "xtol": 1e-3,
    "diff_step": 1e-4,
}

####################################################################################################
from echopop.inversion import InversionBase
from echopop.validators.base import BaseDataFrame, BaseDictionary
from pydantic import ConfigDict, field_validator, Field, model_validator, RootModel, SerializeAsAny
import pandera as pa
import pandera.typing as pat
import echopop.validators as val
import awkward as awk

    
MODEL_PARAMETERS = {
    "number_density": {"value": 500., "min": 10., "max": 10000., "vary": True},
    "theta_mean": {"value": 10., "min": 0., "max": 90., "vary": True},
    "theta_sd": {"value": 20., "min": 0., "max": 90., "vary": False},
    "length_mean": {"value": 0.030, "min": 0.008, "max": 0.040, "vary": True},
    "length_sd_norm": {"value": 0.15, "min": 0.05, "max": 0.15, "vary": False},
    "g": {"value": 1.015, "min": 1.015, "max": 1.060, "vary": False},
    "h": {"value": 1.020, "min": 1.015, "max": 1.060, "vary":False},
    "radius_of_curvature_ratio": {"value": 3.0, "min": 0.5, "max": 100.0, "vary": True},
    "length_radius_ratio": {"value": 18.2, "min": 14.0, "max": 20.0, "vary": True},
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




MOD = InversionMatrix(sv_data, SIMULATION_SETTINGS)
MOD.build_scattering_model(
    MODEL_PARAMETERS, MODEL_SETTINGS
)
self = MOD
self.optimization_kwargs = OPTIMIZATION_KWARGS


def parameterize(
    Sv_measured: pd.Series,
    scattering_parameters: Dict[str, Any],
    simulation_parameters: Dict[str, Any],
    model_parameters: Dict[str, Any],
    **kwargs,
) -> List[Minimizer]:
    """
    Generate `lmfit.Minimizer` object used for optimizing theroetical Sv (dB re. m^-1) based on
    measured values
    """

    # Prepare parameters and minimization procedure
    parameters = prepare_optimization(
        scattering_parameters,
        simulation_parameters=simulation_parameters,
    )

    # Find which values are below the defined threshold
    valid_idx = np.argwhere(Sv_measured > -999.).flatten()

    # Convert to a `numpy.ndarray`
    if not isinstance(Sv_measured, np.ndarray):
        Sv_measured = Sv_measured.to_numpy()[valid_idx]
    else:
        Sv_measured = Sv_measured[valid_idx]

    # Generate `Minimizer` function class required for bounded optimization
    # ---- This will generate per realization within a List
    return [
        Minimizer(
            simulate_Sv_fit,
            parameters[idx],
            fcn_args=(
                Sv_measured,
                {**processing_parameters, **{"center_frequencies": center_frequencies[valid_idx]}},
            ),
            nan_policy="omit",
        )
        for idx in parameters
    ]


# Initialize the Minimizer class for each index
sv_data["minimizer"] = np.nan

# Convert to an `object` for flexible datatyping
sv_data["minimizer"] = sv_data["minimizer"].astype(object)

# Create list of `lmfit.Minimizer` objects for all realizations
sv_data["minimizer"] = sv_data["sv_mean"].apply(
    parameterize,
    axis=1,
    args=(
        scattering_parameters,
        SIMULATION_PARAMETERS,
        SCATTERING_MODEL,
    )
)

# Otherwise, produce just the single-set
if SIMULATION_PARAMETERS["monte_carlo"]:   
    scattering_parameters=monte_carlo_initialization(
        scattering_parameters, 
        n_realizations=SIMULATION_PARAMETERS["n_realizations"]
    )

# Convert to a `lmfit.Parameters` object
scattering_lmfit = {
    k: dict_to_Parameters(v)
    for k, v in scattering_parameters.items()
}

# Define label for verbosity
sv_data["label"] = [
    "; ".join(
        f"{name if name is not None else 'index'}: {val}"
        for name, val in zip(sv_data.index.names, (idx if isinstance(idx, tuple) else (idx,)))
    )
    for idx in sv_data.index
]

Sv_measured = sv_data["sv_mean"].iloc[0, :]


from typing import Tuple


length_mean = scattering_parameters["length_mean"]["value"]
length_sd_norm=0.15
theta_mean=10.0
theta_sd=20.0
g=1.02
h=1.025
radius_of_curvature_ratio=3.0
length_radius_ratio=18.2
taper_order=10.
sound_speed_sw=1500.0
frequency=120e3
theta_mean=10.





# Generate frequency intervals centered on the central frequencies
frequencies = generate_frequency_interval(
    np.array(list(center_frequencies)),
    length_sd_norm,
    SCATTERING_MODEL["frequency_interval"],
)

# Compute the acoustic wavenumbers weighted by target size
# ---- Center frequencies
k_center = wavenumber(np.array(list(center_frequencies)), sound_speed_sw)
# ---- Compute ka (center frequencies)
ka_center = k_center * length_mean / length_radius_ratio
# ---- Frequency intervals
# -------- Just wavenumber (`k`)
k = wavenumber(frequencies, sound_speed_sw)
# -------- Now `ka`
ka = k * length_mean / length_radius_ratio

# Compute over a vector of angles (centered on 90 degrees)
theta_values = np.linspace(
    theta_mean - 3.1 * theta_sd,
    theta_mean + 3.1 * theta_sd,
    SCATTERING_MODEL["orientation_distribution"]["bins"],
)
theta_radians = theta_values * np.pi / 180.0

# Compute over vector lengths
length_values = np.linspace(
    length_mean - 3 * (length_sd_norm * length_mean),
    length_mean + 3 * (length_sd_norm * length_mean),
    SCATTERING_MODEL["length_distribution"]["bins"],
)

n_wavelength = 10
n_integration = 50

# Calculate the appropriate number of integration points
# ---- Compute threshold
kL_max = np.nanmax(k * length_mean, axis=1) * (1 + 3.1 * length_sd_norm)
# ---- Adjust number of integration points based on `kL_max`, if needed
n_int
n_int = awk.values_astype(
    np.where(
        kL_max < n_integration, 
        n_integration, 
        np.ceil(kL_max * n_wavelength / (2 * np.pi))
    ),
    int
)
n_int

# Create shape to build position vector and other required arrays
awk.mean(taper)
awk.mean(gamma_tilt)
awk.mean(beta_tilt)
awk.mean(r_pos)
awk.mean(dr_pos)
taper, gamma_tilt, beta_tilt, r_pos, dr_pos = uniformly_bent_cylinder(
    n_segments=n_int, 
    radius_of_curvature_ratio=radius_of_curvature_ratio, 
    taper_order=taper_order
)
awk.mean(taper)
awk.mean(gamma_tilt)
awk.mean(beta_tilt)
awk.mean(r_pos)
awk.mean(dr_pos)

# taper: float,
# gamma_tilt: np.ndarray[float],
# beta_tilt: np.ndarray[float],
# r_pos: np.ndarray[float],
# dr_pos: np.ndarray[float],
# length_radius_ratio: float,
# g: Union[np.ndarray[float], float],
# h: Union[np.ndarray[float], float],
# ka: Union[np.ndarray[float], float],
# theta: Union[np.ndarray[float], float],

def reflection_coefficient(
    g: Union[np.ndarray, float],
    h: Union[np.ndarray, float],
) -> np.ndarray[float]:
    """
    Compute the reflection coefficient based on material properties
    """

    return (1 - g * h * h) / (g * h * h) - (g - 1) / g

from scipy.special import j1

# Get the array sizes for later broadcasting
# ---- Number of frequencies/wavenumbers
n_k = awk.num(ka, axis=-1)
# ---- Number of segments and minimum number of integration points
n_segments = awk.num(r_pos, axis=-1)
# ---- Number of orientation values
n_theta = len(theta_values)

# Compute the reflection coefficient, `C_b`
C_b = reflection_coefficient(g, h)

# Vectorized computation without for loop
# ---- Adjust `ka` to account for body shape tapering (vectorized across all frequencies)
awk.mean(ka_tapered)
ka_tapered = awk.Array([
    np.array(ka[i, :n_k[i]]).reshape(-1, 1) * taper[i, :n_segments[i]] / h
    for i in range(len(n_k))
])
awk.mean(ka_tapered)
# ---- Adjust along-axis tilt angles and slopes to be relative to the incident planar wave (vectorized)
# -------- Along-axis curvature slopes
awk.mean(delta_gamma_cos)
delta_gamma_cos = awk.Array([
    np.cos(np.array(gamma_tilt[i, :n_segments[i]]).reshape(-1, 1) - theta_radians)
    for i in range(len(n_k))
])
awk.mean(delta_gamma_cos)
# -------- Along-axis tilt angles between segments
awk.mean(delta_theta_cos)
delta_theta_cos = awk.Array([
    np.abs(np.cos(np.array(beta_tilt[i, :n_segments[i]]).reshape(-1, 1) - theta_radians))
    for i in range(len(n_k))
])
awk.mean(delta_theta_cos)

# ---- Generate matrices (vectorized)
awk.mean(M1)
M1 = awk.Array([
    length_radius_ratio * np.array(ka[i, :n_k[i]]).reshape(-1, 1) * (np.array(r_pos[i, :n_segments[i]]) / h)
    for i in range(len(n_k))
])
awk.mean(M1)

awk.mean(M2)
M2 = awk.Array([
    h**2 * C_b * np.array(dr_pos[i, :n_segments[i]]) / 4
    for i in range(len(n_k))
])
awk.mean(M2)

# ---- Vectorized Bessel function computation
awk.mean(J1_vals)
ka_last = awk.Array([np.array(ka[i, n_k[i]-1]) for i in range(len(n_k))])
ka_norm = awk.Array([
    np.linspace(-2 * ka_last[i], 2 * ka_last[i], 2 * n_segments[i])
    for i in range(len(n_k))
])
J1_vals = awk.Array([j1(ka_norm[i]) for i in range(len(ka_norm))])
awk.mean(J1_vals)



# ---- Vectorized form function computation
awk.mean(f_bs)
f_bs = awk.Array([
    _compute_single_frequency(
        ka_tapered[i], 
        delta_gamma_cos[i], delta_theta_cos[i],
        M1[i], M2[i], 
        ka_norm[i], J1_vals[i],
        n_k[i], n_segments[i], n_theta
    )
    for i in range(len(n_k))
])
awk.mean(f_bs)

# This below is EQUIVALENT to the argument `angle`
angle=theta_values
# This below is EQUIVALENT to the argument `f_bs`
form_function=f_bs
theta_mean
theta_sd
# This below is EQUIVALENT to the argument `distribution`
distribution=SCATTERING_MODEL["orientation_distribution"]["family"]

awk.mean(PDF)
# ---- Get interval
orientation_interval = np.diff(angle).mean()
# ---- Compute the PDF
PDF = (
    orientation_interval
    * np.exp(-0.5 * (angle - theta_mean) ** 2 / theta_sd**2)
    / (np.sqrt(2 * np.pi) * theta_sd)
)
awk.mean(PDF)

f_bs_orientation = awk.Array([
    np.sqrt(np.matmul((np.array(form_function[i][0]).real ** 2 + 
                        np.array(form_function[i][0]).imag ** 2), PDF))
    for i in range(len(form_function))
])

length_deviation = length_mean * length_sd_norm

length=length_values
ka
ka_center
f_bs_orientation
length_mean
length_deviation

# Normalize the length values, if needed
length_norm = length / length_mean
# ---- Also normalize the standard deviation
length_sd_norm = length_deviation / length_mean

# ---- Get the interval
length_interval = np.diff(length_norm).mean()
# ---- Compute the PDF
PDF = (
    length_interval
    * np.exp(-0.5 * (length_norm - 1) ** 2 / length_sd_norm**2)
    / (np.sqrt(2 * np.pi) * length_sd_norm)
)

# Vectorized computation - compute length-weighted ka for all frequencies
np.mean(ka_weighted)
ka_weighted = length_norm * ka_center.reshape(-1, 1)
np.mean(ka_weighted)
#
# Vectorized computation using awkward arrays
awk.mean(sigma_bs_length)
sigma_bs_length = awk.Array([
    (
        length_norm**2
        * PDF
        * np.interp(ka_weighted[i], np.array(ka[i]), np.array(f_bs_orientation[i]) ** 2)
    ).sum()
    for i in range(len(f_bs_orientation))
])
awk.mean(sigma_bs_length)

# ---- Convert to sigma_bs (linear backscattering cross-section)
awk.mean(sigma_bs)
sigma_bs = sigma_bs_length * (length_mean) ** 2
awk.mean(sigma_bs)

# Switch to logarithmic domain to compute S_V (volumetric backscattering strength)
Sv_prediction
Sv_prediction = 10 * np.log10(number_density * sigma_bs)
Sv_prediction

TS
TS = 10 * np.log10(awk.to_numpy(sigma_bs))
TS
####################
scattering_params = scattering_lmfit[0]
model_params = SCATTERING_MODEL
model_params["frequency_interval"] = model_params["frequency_interval"] * 1e3
simulation_params = SIMULATION_PARAMETERS
simulation_params["parameter_bounds"] = parameter_bounds
center_frequencies
####################
# Extract parameter values from dictionary for parsing
scattering_params = scattering_params.valuesdict()

# Rescale parameters to their original scales
if simulation_params["scale_parameters"]:
    scattering_params = minmax_normalize(
        parameter_sets=scattering_params, 
        inverse=True, 
        inverse_reference=simulation_params["parameter_bounds"],
    )
    
# Compute the acoustic wavenumber for the center frequencies (by default)
k_center = wavenumber(np.array(list(center_frequencies)), model_params["sound_speed_sw"])

def pcdwba(
    k_center: np.ndarray[float],
    scattering_params: Dict[str, Any],
    model_params: Dict[str, Any],
) -> List[np.ndarray[complex]]:
    """
    Phase-compensated distorted wave Born approximation (DWBA)

    Defined as [1]_,

    ..math::
        f(\phi)=(kaT)^2\frac{J_1(2kaT \cos{\theta})}{2kaT \cos{\theta}}
        e^{i\varepsilon ka(\vec{r}_{pos}/h)\cos{\beta}}

    where :math:`f` is the scattering amplitude, :math:`\phi` is the scatterer orientation angle
    relative to the incident sound wave, :math:`k` is the acoustic wavenumber, :math:`T` is the
    taper coefficient, :math:`a` is the radius at a specific point along the body, :math:`J_1` is
    the cylindrical Bessel function of the first kind, :math:`\theta` is the orientation angle
    relative to the incident sound wave at a specific point along the body, :math:`\varepsilon` is
    the length-to-radius ratio, :math:`\vec{r}_{pos}` is the positional vector, :math:`h` is the
    soundspeed contrast, and :math:`\beta` is the orientation of a specific point along the body.

    References
    ----------
    ..[1] Chu, D., and Ye, Z. (1999). A phase-compensated distorted wave Born approximation
    representation of the bistatic scattering by weakly scattering objects: Application to
    zooplankton
    """

    # Compute the radius-weighted wavenumber
    ka_center = (
        k_center * scattering_params["length_mean"] / scattering_params["length_radius_ratio"]
    )

    # Normalize the length standard deviation
    if "length_sd_norm" not in scattering_params and "length_deviation" in scattering_params:
        length_sd_norm = scattering_params["length_deviation"] / scattering_params["length_mean"]
    else:
        length_sd_norm = scattering_params["length_sd_norm"]

    # Generate frequency intervals centered on the central frequencies
    frequencies = generate_frequency_interval(
        np.array(list(center_frequencies)) * 1e3,
        length_sd_norm,
        model_params["frequency_interval"],
    )

    # Compute the wavenumbers and radius-weighted wavenumbers for the frequency intervals
    # ---- wavenumber
    k = wavenumber(frequencies, model_params["sound_speed_sw"])
    # ---- ka
    ka = k * scattering_params["length_mean"] / scattering_params["length_radius_ratio"]

    # Compute over a vector of angles (centered on 90 degrees)
    theta_values = np.linspace(
        theta_mean - 3.1 * theta_sd,
        theta_mean + 3.1 * theta_sd,
        model_params["orientation_distribution"]["bins"],
    )
    # ---- Convert to radians
    theta_radians = theta_values * np.pi / 180.0


# ---- Compute ka (center frequencies)
ka_center = k_center * parameters_dict["length_mean"] / parameters_dict["length_radius_ratio"]
# ---- Frequency intervals
# -------- Just wavenumber (`k`)
k = wavenumber(frequencies, model_params["sound_speed_sw"])
# -------- Now `ka`
ka = k * parameters_dict["length_mean"] / parameters_dict["length_radius_ratio"]



    
# Normalize the length standard deviation
if "length_sd_norm" not in parameters_dict and "length_deviation" in parameters_dict:
    length_sd_norm = parameters_dict["length_deviation"] / parameters_dict["length_mean"]
else:
    length_sd_norm = parameters_dict["length_sd_norm"]






length_values = np.linspace(
    length_mean - 3 * (length_sd_norm * length_mean),
    length_mean + 3 * (length_sd_norm * length_mean),
    model_params["length_distribution"]["bins"],
)

# Compute over a vector of angles (centered on 0 degrees)
theta_values = np.linspace(
    theta_mean - 3.1 * theta_sd,
    theta_mean + 3.1 * theta_sd,
    model_parameters["orientation_distribution"]["bins"],
)
# ---- Convert to radians
theta_radians = theta_values * np.pi / 180.0

def pcdwba(
    length_mean: float,
    length_radius_ratio: float,
    length_sd_norm: float,
    theta_mean: float,
    theta_sd: float,
    center_frequencies: Dict[str, Any],
    frequency_interval: float,
    sound_speed_sw: float,
    density_sw: float,
    **kwargs,
) -> List[np.ndarray[complex]]:
    """
    Phase-compensated distorted wave Born approximation (DWBA)

    Defined as [1]_,

    ..math::
        f(\phi)=(kaT)^2\frac{J_1(2kaT \cos{\theta})}{2kaT \cos{\theta}}
        e^{i\varepsilon ka(\vec{r}_{pos}/h)\cos{\beta}}

    where :math:`f` is the scattering amplitude, :math:`\phi` is the scatterer orientation angle
    relative to the incident sound wave, :math:`k` is the acoustic wavenumber, :math:`T` is the
    taper coefficient, :math:`a` is the radius at a specific point along the body, :math:`J_1` is
    the cylindrical Bessel function of the first kind, :math:`\theta` is the orientation angle
    relative to the incident sound wave at a specific point along the body, :math:`\varepsilon` is
    the length-to-radius ratio, :math:`\vec{r}_{pos}` is the positional vector, :math:`h` is the
    soundspeed contrast, and :math:`\beta` is the orientation of a specific point along the body.

    References
    ----------
    ..[1] Chu, D., and Ye, Z. (1999). A phase-compensated distorted wave Born approximation
    representation of the bistatic scattering by weakly scattering objects: Application to
    zooplankton
    """

    # Generate the frequency interval
    frequencies = generate_frequency_interval(
        np.array(list(center_frequencies)) * 1e3,
        length_sd_norm,
        frequency_interval,
    )

    # Compute the acoustic wavenumbers weighted by target size
    # ---- Center frequencies
    k_center = wavenumber(np.array(list(center_frequencies)), sound_speed_sw)
    # ---- Compute ka (center frequencies)
    ka_center = k_center * length_mean / length_radius_ratio
    # ---- Frequency intervals
    # -------- Just wavenumber (`k`)
    k = wavenumber(frequencies, sound_speed_sw)
    # -------- Now `ka`
    ka = k * length_mean / length_radius_ratio

    # Compute over a vector of angles (centered on 0 degrees)
    theta_values = np.linspace(
        theta_mean - 3.1 * theta_sd,
        theta_mean + 3.1 * theta_sd,
        SCATTERING_MODEL["orientation_distribution"]["bins"],
    )
    # ---- Convert to radians
    theta_radians = theta_values * np.pi / 180.0

