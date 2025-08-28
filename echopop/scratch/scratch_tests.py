from typing import Any, Dict, Literal, Optional, Tuple, Union
from lmfit import Minimizer, Parameters
import awkward as awk
import numpy as np
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from echopop.nwfsc_feat import ingest_nasc
from echopop.validators.base import BaseDataFrame, BaseDictionary
from pydantic import ConfigDict, field_validator, Field, model_validator, RootModel, SerializeAsAny, ValidationError
import pandas as pd
import numpy as np
from scipy.special import j1
import time
from echopop import validators as val
import warnings
from echopop.inversion.inversion_base import InversionBase
# !!! ADD WARNING IF NUMBER OF PARAMETERS EXCEEDS NUMBER OF POINTS
# !!! ADD RNG SEED SETTER 
# # Initialize the random number generator
# self.rng = np.random.default_rng(resample_seed)

class ScatterDF(BaseDataFrame):
    """
    Validate hierarchical columns in a MultiIndex DataFrame.

    Ensures the top-level column names include 'sv_mean', 'nasc', and 'thickness_mean',
    and that the second-level column names ('frequency') are numeric. All values
    under each frequency must be floats.
    """
    class Config:
        multiindex_strict = True
        multiindex_coerce = True

    @classmethod
    def pre_validate(cls, df: pd.DataFrame) -> pd.DataFrame:
        # Raise Error if not a MultiIndex DataFrame
        if not isinstance(df.columns, pd.MultiIndex):
            raise TypeError("Expected MultiIndex columns.")
        
        # Check for required top-level columns
        # ---- Get the top level indices
        top_level = set(df.columns.get_level_values(0))
        # ---- Identify missing columns
        missing_columns = {"sv_mean", "nasc", "thickness_mean"} - top_level
        # ---- Raise Error, if needed
        if missing_columns:
            # ---- Format
            missing_str = ", ".join(f"'{col}'" for col in missing_columns)
            raise KeyError(
                f"Missing required top-level columns: {missing_str}."
            )
            
        # Check for nested column 'frequency'
        if "frequency" not in df.columns.names:
            # ---- Get current column names
            current_names = list(df.columns.names)
            # ---- Format 
            current_str = ", ".join(f"'{col}'" for col in current_names if col)
            # ---- Raise Error
            raise KeyError(
                f"Missing required nested column index 'frequency'. Current column indices are: "
                f"{current_str}."
            )
            
        return df

class EnvironmentParameters(
    BaseDictionary,
    arbitrary_types_allowed=True,
    title="scattering model environment/medium parameters",
):
    sound_speed_sw: float = Field(default=1500., gt=0., allow_inf_nan=False)
    density_sw: float = Field(default=1025., gt=0., allow_inf_nan=False)

class SimulationParameters(
    BaseDictionary,
    arbitrary_types_allowed=True,
    title="inversion simulation parameters",
):
    environment: EnvironmentParameters = Field(
        default_factory=EnvironmentParameters
    )
    monte_carlo: Optional[bool] = Field(default=None)
    mc_realizations: Optional[int] = Field(default=None)
    mc_seed: Optional[int] = Field(default=None) 
    scale_parameters: bool = Field(default=True)
    minimum_frequency_count: int = Field(default=2)    

    @model_validator(mode="after")
    def validate_monte_carlo_realizations(self):
        # Get `monte_carlo`
        monte_carlo = self.monte_carlo or False

        # Default n_realizations if monte_carlo
        if monte_carlo and self.mc_realizations is None:
            self.mc_realizations = 100

        # Always update monte_carlo
        self.monte_carlo = monte_carlo
        
        return self

class ValidateInversionMatrix(
    BaseDictionary,
    arbitrary_types_allowed=True,
    title="scattering model inversion analysis",
):
    data: pd.DataFrame
    simulation_settings: SimulationParameters = Field(
        default_factory=SimulationParameters
    )
    
    @field_validator("data", mode="after")
    def validate_data(cls, v):
        # Validate with pandera
        return ScatterDF.validate(v)

ValidateInversionMatrix.create(**{"data": sv_data})

class ValidatePCDWBAParams(
    BaseDictionary,
    arbitrary_types_allowed=True,
    title="PCDWBA model parameters",
):
    g: float = Field(gt=0., le=np.inf, allow_inf_nan=True)
    h: float = Field(gt=0., le=np.inf, allow_inf_nan=True)
    length_mean: float = Field(gt=0., le=np.inf, allow_inf_nan=True)
    length_radius_ratio: float = Field(gt=0., le=np.inf, allow_inf_nan=True)
    length_sd_norm: float = Field(ge=0., le=np.inf, allow_inf_nan=True)
    number_density: float = Field(gt=0., le=np.inf, allow_inf_nan=True)
    radius_of_curvature_ratio: float = Field(gt=0., le=np.inf, allow_inf_nan=True)
    theta_mean: float = Field(ge=-180., le=180., allow_inf_nan=False)
    theta_sd: float = Field(ge=0., le=np.inf, allow_inf_nan=True)

class DistributionParameters(
    BaseDictionary,
    arbitrary_types_allowed=True,
    title="distribution parameters",    
):
    bins: int = Field(default=30)
    family: Literal["gaussian", "uniform"] = Field(default="gaussian")
class ValidatePCDWBASettings(
    BaseDictionary,
    arbitrary_types_allowed=True,
    title="PCDWBA model settings and additional arguments",
):
    frequency_interval: float = Field(gt=0., allow_inf_nan=False)
    length_distribution: DistributionParameters = Field(
        default_factory=DistributionParameters
    )
    n_integration: int = Field(default=50, gt=0)
    n_wavelength: int = Field(default=10, gt=0)
    orientation_distribution: DistributionParameters = Field(
        default_factory=DistributionParameters
    )
    taper_order: float = Field(default=10., gt=0., allow_inf_nan=False)
    type: str

class SingleParameter(
    BaseDictionary,
    arbitrary_types_allowed=True,
    title="values for lmfit.Parameters class required for optimization",
    ):    
    min: Optional[float] = Field(default=-np.inf, allow_inf_nan=True)
    value: float = Field(allow_inf_nan=False)
    max: Optional[float] = Field(default=np.inf, allow_inf_nan=True)
    vary: bool = Field(default=False)

    @model_validator(mode="after")
    def check_bounds(self):
        if self.min > self.max:
            raise ValueError(f"min ({self.min}) cannot be greater than max ({self.max}).")
        if not (self.min <= self.value <= self.max):
            raise ValueError(
                f"value ({self.value}) must be between min ({self.min}) and max ({self.max})."
            )
        return self

class ModelInputParameters(
    RootModel[Dict[str, "SingleParameter"]],
    ):

    model_config = ConfigDict(
        arbitrary_types_allowed=True,
        title="scattering model parameters",
    )

class ModelSettingsParameters(
    BaseDictionary,
    arbitrary_types_allowed=True,
    title="general model settings", 
):
    type: str  
    model_config = ConfigDict(extra="allow")

class ValidateBuildModelArgs(
    BaseDictionary,
    arbitrary_types_allowed=True,
    title="scattering model preparation and assembly",
):
    model_parameters: ModelInputParameters
    model_settings: ModelSettingsParameters
    model_config = ConfigDict(protected_namespaces=())

    @staticmethod
    def _check_variable(params: Dict[str, Any], variable: str, validator: Any):

        # Get the schema JSON
        properties = validator.model_json_schema()["properties"]

        # Get fill variable if needed
        if variable in ["min", "max"]:
            # ---- Format variable
            var = f"{variable}imum"
            var_ex = f"exclusive{variable.capitalize()}imum"
            # ---- Get variable defaults for boundaries
            var_defaults = {
                k: (v[var] if var in v else v[var_ex] + np.finfo(float).eps)
                for k, v in properties.items()
                if var in v or var_ex in v
            }
        else:
            var_defaults = {k: None for k in properties}
            
        # Set up the values
        param_slice = {k: v[variable] if variable in v else var_defaults[k] 
                       for k, v in params.items()}

        # Get the validation scheme
        return {k: {variable: v} for k, v in validator.create(**param_slice).items()}

    @model_validator(mode="after")
    def validate_model_parameterization(self):

        # Check for model-type
        # ---- Dump the model
        model_settings = self.model_settings.model_dump()
        # ---- Check against reference
        if model_settings["type"] not in SCATTERING_MODEL_PARAMETERS:
            raise LookupError(
                f"Scattering model '{model_settings["type"]}' could not be found."
            )

        # Get the model-specific validator
        model_validators = SCATTERING_MODEL_PARAMETERS[model_settings["type"]]

        # Validate `model_parameters`
        # ---- Dump the model
        model_parameters = self.model_parameters.model_dump()
        # ---- Check minima
        min_params = self._check_variable(model_parameters, "min", model_validators["parameters"])
        # ---- Set up values
        value_params = self._check_variable(model_parameters, "value", 
                                            model_validators["parameters"])
        # --- Set up maximum values
        max_params = self._check_variable(model_parameters, "max", model_validators["parameters"])
        # ---- Get the `vary` definitions
        vary_params = {k: {"vary": v["vary"]} for k, v in model_parameters.items()}
        # ---- Rebuild the model parameters
        self.model_parameters = ModelInputParameters(
            {
                k: {**min_params.get(k, {}), **value_params.get(k, {}), 
                    **max_params.get(k, {}), **vary_params.get(k, {})} 
                for k in set(min_params) | set(value_params) | set(max_params) | set(vary_params)
            }
        ) 
        
        # Validate `model_settings`
        self.model_settings = ModelSettingsParameters.model_validate(
            model_validators["settings"].create(**self.model_settings.model_dump())
        )
        
        return self

# Interpolate the Bessel function output over all orientations
def _fast_bessel_interp(ARG_flat, ka_norm_np, J1_np, n_theta):
    """Pure numpy interpolation - no compilation overhead"""
    result = np.zeros((ARG_flat.shape[0], n_theta))
    
    # Use numpy's vectorized operations instead of explicit loops
    for th in range(n_theta):
        result[:, th] = np.interp(ARG_flat[:, th], ka_norm_np, J1_np)
    
    return result

def pcdwba(
    center_frequencies: np.ndarray[float],
    length_mean: float,
    length_sd_norm: float,
    length_radius_ratio: float,
    taper_order: float,
    radius_of_curvature_ratio: float,
    theta_mean: float,  
    theta_sd: float,
    orientation_distribution: Dict[str, Any],
    g: float,
    h: float,
    sound_speed_sw: float,
    frequency_interval: float,
    n_integration: int,
    n_wavelength: int, 
    number_density: Optional[float] = None,
    length_distribution: Optional[Dict[str, Any]] = None,
    **kwargs
):
    
    # Use pre-computed constants
    eps = np.finfo(float).eps

    # Pre-allocate arrays based on input size
    n_theta = orientation_distribution["bins"]
    n_length = length_distribution["bins"]

    # Generate frequency intervals centered on the central frequencies
    frequencies = generate_frequency_interval(
        center_frequencies,
        length_sd_norm,
        frequency_interval
    )
    
    # Compute the acoustic wavenumbers weighted by target size
    # ---- Center frequencies
    k_c = wavenumber(center_frequencies, sound_speed_sw)
    # ---- Compute ka (center frequencies)
    ka_c = k_c * length_mean / length_radius_ratio
    # ---- Frequency intervals
    # -------- Just wavenumber (`k`)
    k_f = wavenumber(frequencies, sound_speed_sw)
    # -------- Now `ka`
    ka_f = k_f * length_mean / length_radius_ratio
    
    # Compute over a vector of angles (centered on 90 degrees)
    theta_values = np.linspace(
        theta_mean - 3.1 * theta_sd,
        theta_mean + 3.1 * theta_sd,
        n_theta,
    )
    # ---- Convert to radians
    theta_radians = theta_values * np.pi / 180.0
    
    # Calculate the appropriate number of integration points
    # ---- Compute threshold
    kL_max = np.nanmax(k_f * length_mean, axis=1) * (1 + 3.1 * length_sd_norm)
    # ---- Adjust number of integration points based on `kL_max`, if needed
    n_int = awk.values_astype(
        np.where(
            kL_max < n_integration, 
            n_integration, 
            np.ceil(kL_max * n_wavelength / (2 * np.pi))
        ),
        int
    )
    
    # Create shape/position matrix
    taper, gamma_tilt, beta_tilt, r_pos, dr_pos = uniformly_bent_cylinder(
        n_int, radius_of_curvature_ratio, taper_order
    )

    # Get the array sizes for later broadcasting
    # ---- Number of frequencies/wavenumbers
    n_k = awk.num(ka_f, axis=-1)
    # ---- Number of segments and minimum number of integration points
    n_segments = awk.num(r_pos, axis=-1)
    # ---- Number of orientation values
    n_theta = len(theta_radians)

    # Compute the reflection coefficient, `C_b`
    C_b = reflection_coefficient(g, h)

    # Adjust `ka` to account for body shape tapering
    ka_f_tapered = awk.Array([
        np.outer(np.array(ka_f[i, :n_k[i]]), taper[i, :n_segments[i]]) / h
        for i in range(len(n_k))
    ])

    # Adjust the along-axis tilt angles and slopes so they are relative to the incident planar 
    # wave
    # ---- Curvature slopes
    delta_gamma_cos = awk.Array([
        np.cos(np.array(gamma_tilt[i, :n_segments[i]]).reshape(-1, 1) - theta_radians)
        for i in range(len(n_k))
    ])
    # ---- Along-axis intersegment tilt angles
    delta_theta_cos = awk.Array([
        np.abs(np.cos(np.array(beta_tilt[i, :n_segments[i]]).reshape(-1, 1) - theta_radians))
        for i in range(len(n_k))
    ])

    # Calculate the phase term
    # ---- Compute the bulk of the exponentiated matrix
    M1 = awk.Array([
        length_radius_ratio * np.array(ka_f[i, :n_k[i]]).reshape(-1, 1) * 
        (np.array(r_pos[i, :n_segments[i]]) / h)
        for i in range(len(n_k))
    ])
    # ---- Compute phase
    phase = awk.Array([
        np.exp(1j * np.array(M1[i][:, :, np.newaxis]) * 
               np.array(delta_gamma_cos[i][np.newaxis, :, :]))
        for i in range(len(n_k))
    ])

    # Calculate the effect of the material properties with respect to the position matrix derivative
    M2 = awk.Array([
        h**2 * C_b * np.array(dr_pos[i, :n_segments[i]]) / 4
        for i in range(len(n_k))
    ])

    # Prepare the ka_f values for the Bessel function
    ka_f_last = awk.Array([np.array(ka_f[i, n_k[i]-1]) for i in range(len(n_k))])
    ka_f_norm = awk.Array([
        np.linspace(-2 * ka_f_last[i], 2 * ka_f_last[i], 2 * n_segments[i])
        for i in range(len(n_k))
    ])
    # ---- Compute the cylindrical Bessel function(s) of the first kind
    J1 = awk.Array([j1(ka_f_norm[i]) for i in range(len(ka_f_norm))])

    # Calculate product of tapered ka_f and along-axis tilts
    ARG = 2 * awk.Array([
        np.array(ka_f_tapered[i][:, :, np.newaxis]) * np.array(delta_theta_cos[i][np.newaxis, :, :])
        for i in range(len(n_k))
    ]) + eps
    
    J1_interp = awk.Array([
        _fast_bessel_interp(
            np.array(ARG[i]).ravel(order="F").reshape((n_k[i] * n_segments[i], n_theta), order="F"),
            np.array(ka_f_norm[i]), np.array(J1[i]), n_theta
        )
        for i in range(len(n_k))
    ])

    # Normalize the interpolated Bessel function
    J1_norm = awk.Array([
        (
            np.array(J1_interp[i]) / 
            np.array(ARG[i]).ravel(order="F").reshape((n_k[i] * n_segments[i], n_theta), order="F")
        ).reshape((n_k[i], n_segments[i], n_theta), order="F")
        for i in range(len(n_k))
    ])

    # Calculate the primary DWBA equation
    f_j = ka_f_tapered ** 2 * J1_norm * phase

    # Matrix multiplication via Einstein summation to get the final complex form function result
    f_bs = awk.Array([np.einsum("ijk, j->ik", f_j[i], M2[i]) for i in range(len(n_k))])

    # Orientation averaging
    f_bs_orientation = orientation_average(
        theta_values, f_bs, theta_mean, theta_sd, orientation_distribution["family"]
    )

    # Compute over vector lengths
    length_values = np.linspace(
        length_mean - 3 * (length_sd_norm * length_mean),
        length_mean + 3 * (length_sd_norm * length_mean),
        n_length,
    )

    # Length-averaged sigma_bs (normalized to length)
    sigma_bs_length = length_average(
        length_values, ka_f, ka_c, f_bs_orientation, length_mean, length_mean * length_sd_norm, 
        length_distribution["family"]
    )

    # Convert to sigma_bs (linear backscattering cross-section)
    sigma_bs = sigma_bs_length * (length_mean) ** 2

    # Switch to logarithmic domain to compute S_V (volumetric backscattering strength)
    Sv_prediction = 10 * np.log10(number_density * sigma_bs)
    
    return Sv_prediction

SCATTERING_MODEL_PARAMETERS = {
    "pcdwba": {
        "function": pcdwba,
        "parameters": ValidatePCDWBAParams,
        "settings": ValidatePCDWBASettings,
    }
}


def read_echoview_sv(filename: Path, 
                     impute_coordinates: bool = True, 
                     transect_num: Optional[float] = None, 
                     validator: Optional[Any]=None,):

    # Read in the defined CSV file
    sv_df = ingest_nasc.read_echoview_export(filename, validator)
    
    # Don't read in if the file contents are empty
    if sv_df.empty or sv_df.dropna(axis=1, how="all").empty:
        return None

    # Add transect number, if defined
    if transect_num is not None:
        sv_df["transect_num"] = transect_num

    # Append filename
    sv_df["filename"] = filename.as_posix()

    # Fix latitude and longitude
    # ---- Latitude
    if "latitude" in sv_df.columns and impute_coordinates:
        ingest_nasc.impute_bad_coordinates(sv_df, "latitude")
    # ---- Longitude
    if "longitude" in sv_df.columns and impute_coordinates:
        ingest_nasc.impute_bad_coordinates(sv_df, "longitude")

    # Return the cleaned DataFrame
    return sv_df

def echoview_sv_to_df(
    filtered_df: pd.DataFrame, impute_coordinates: bool = True
):
    return [
        read_echoview_sv(row["file_path"], impute_coordinates, row["transect_num"])
        for _, row in filtered_df.iterrows()
    ]    

def apply_Sv_thresholds(data: pd.DataFrame, thresholds: Dict[str, Any]):

    # Create copy
    data = data.copy()
    
    # Create mapping series for minimum and maximum Sv thresholds
    freq_min_map = pd.Series({freq: vals["min"] for freq, vals in thresholds.items()})
    freq_max_map = pd.Series({freq: vals["max"] for freq, vals in thresholds.items()})

    # Apply frequency-specific thresholds
    min_thresh = data["frequency"].map(freq_min_map)
    max_thresh = data["frequency"].map(freq_max_map)

    # Create mask
    if "sv_mean" in data.columns:
        mask = (data["sv_mean"] < min_thresh) | (data["sv_mean"] > max_thresh)
    else:
        raise KeyError(
            "Could not apply Sv thresholds. Column 'sv_mean' could not be found within the "
            "ingested acoustic DataFrame."
        )

    # Apply the thresholds
    data.loc[mask, "sv_mean"] = -999.

    return data

def sv_to_nasc(sv_linear, thickness_mean):
    """
    Convert volume backscattering coefficient (sv, linear) to NASC (sA)
    
    Parameters:
    -----------
    sv_linear : array-like
        Volume backscattering coefficient in linear units (m^-1)
    thickness_mean : array-like  
        Mean thickness of integration layer (m)
    
    Returns:
    --------
    nasc : array-like
        Nautical area scattering coefficient (m^2 nmi^-2)
    """
    # From the equations in the image:
    # sa = ∫ sv dz  (area backscattering coefficient)
    # sA = 4π (1852)² sa  (NASC conversion)
    
    # Constants
    STERADIANS_SPHERE = 4 * np.pi  # 4π steradians
    METERS_PER_NMILE_SQUARED = 1852**2  # (1852 m/nmi)²
    
    # Calculate area backscattering coefficient (sa) by integrating sv over thickness
    sa = sv_linear * thickness_mean  # ∫ sv dz approximated as sv × thickness
    
    # Convert to NASC using sA = 4π(1852)² × sa
    nasc = STERADIANS_SPHERE * METERS_PER_NMILE_SQUARED * sa
    
    return nasc

def aggregate_cells(data):
    
    # Find the overlapping columns
    valid_idx_cols = [
        col for col in data.columns 
        if col in ["transect_num", "longitude", "latitude", "interval", "layer"]
    ]
    
    # Compute pivot table
    return data.pivot_table(
        index=valid_idx_cols,
        columns=["frequency"],
        values=["sv_mean", "nasc", "thickness_mean"]
    ).fillna({"nasc": 0., "sv_mean": -999., "thickness_mean": 0.})
        

def aggregate_intervals(
    data: pd.DataFrame,
) -> pd.DataFrame:
    """
    Aggregate integration over intervals along each transect
    """

    # Create copy
    data = data.copy()

    # Find the overlapping columns
    valid_idx_cols = [
        col for col in data.columns 
        if col in ["transect_num", "longitude", "latitude", "interval"]
    ]

    # 'Interval' must be present
    if "interval" not in valid_idx_cols:
        raise KeyError(
            "Integration over intervals requires column 'interval' to be present within the "
            "ingested acoustic DataFrame."
        )

    # Weight the linear Sv values by thickness
    data["sv_t"] = data["sv_mean_linear"] * data["thickness_mean"]

    data.groupby()

    # Aggregate the values over each interval
    data_pvt = data.groupby(valid_idx_cols + ["frequency"]).agg({
        "sv_t": lambda x: 10 * np.log10(x.sum()),
        "nasc": "sum", 
        "thickness_mean": "sum"
    })

    # Rename the Sv column
    data_pvt.rename(columns={"sv_t": "sv_mean"}, inplace=True)

    # Unstack for 'frequency'
    return data_pvt.unstack("frequency").fillna(
        {"nasc": 0., "sv_mean": -999., "thickness_mean": 0.}
    )

def aggregate_transects(
    data: pd.DataFrame,
):

    # Check for 'transect_num'
    if "transect_num" not in data.columns:
        raise KeyError(
            "Integration over transects requires column 'transect_num' to be present within the "
            "ingested acoustic DataFrame."
        )
    
    # Create copy
    data = data.copy()

    # Compute distance, if missing
    if "distance" not in data.columns:
        data["distance"] = data["distance_e"] - data["distance_s"]
    
    # Calculate 2D cell areas (needed for calculating the 'line backscattering coefficient', S_L)
    data["cell_area"] = data["distance"] * data["thickness_mean"]

    # Sum the total cell areas and NASC for each transect
    # ---- Cell area
    data["total_cell_area"] = (
        data.groupby(["frequency", "transect_num"])["cell_area"].transform("sum")
    )
    # ---- NASC
    data["total_nasc"] = data.groupby(["frequency", "transect_num"])["nasc"].transform("sum")

    # Calculate weights
    # ---- Cell areas
    data["cell_area_weight"] = data["cell_area"] / data["total_cell_area"]
    # ---- NASC
    data["nasc_weight"] = data["nasc"] / data["total_nasc"]

    # Calculate sv(L)
    data["sv_L"] = data["sv_mean_linear"] * data["cell_area_weight"]

    # Weight the coordinates
    # ---- Longitude
    data["longitude_weight"] = data["longitude"] * data["nasc_weight"]
    # ---- Latitude
    data["latitude_weight"] = data["latitude"] * data["nasc_weight"]

    # Sum the thicknesses per interval
    data["thickness_interval"] = (
        data.groupby(["frequency", "transect_num", "interval"])["thickness_mean"].transform("sum")
    )
    
    # Aggregate the values over each interval
    data_pvt = data.groupby(["frequency", "transect_num"]).agg({
        "longitude_weight": "sum",
        "latitude_weight": "sum",
        "sv_L": lambda x: 10 * np.log10(x.sum()),
        "nasc": "sum", 
        "thickness_interval": "mean"
    })

    # Rename the columns
    data_pvt.rename(columns={"sv_L": "sv_mean", 
                             "longitude_weight": "longitude_weighted", 
                             "latitude_weight": "latitude_weighted", 
                             "thickness_interval": "thickness_mean"},
                    inplace=True)

    # Add the coordinates to the index
    return data_pvt.unstack("frequency").fillna(
        {"nasc": 0., "sv_mean": -999., "thickness_mean": 0.}
    )

def integrate_measurements(
    data: pd.DataFrame,
    method: Literal["cells", "interval", "transect"],
    sv_thresholds: Dict[str, float],
) -> pd.DataFrame:
    """
    Integrate measurements based on transect and interval, or just transect
    """

    # Create copy
    data = data.copy()

    # Apply minimum and maximum Sv thresholds
    sv_thresholded = apply_Sv_thresholds(data, sv_thresholds)

    # Compute the linear volumetric scattering coefficient (sv)
    sv_thresholded["sv_mean_linear"] = 10. ** (sv_thresholded["sv_mean"] / 10.) 

    # Drop empty cells 
    sv_reduced = sv_thresholded.loc[sv_thresholded["sv_mean"] > -999.]

    # Transform frequency 
    sv_reduced.loc[:, "frequency"] = sv_reduced.loc[:, "frequency"] * 1e3

    # If NASC is not a column
    if "nasc" not in sv_reduced.columns:
        # ---- Calculate NASC
        sv_reduced.loc[:, "nasc"] = sv_to_nasc(sv_reduced["sv_mean_linear"], 
                                               sv_reduced["thickness_mean"])

    # Aggregation methods
    if method == "cells":
        sv_indexed = aggregate_cells(sv_reduced)
    elif method == "intervals":
        sv_indexed = aggregate_intervals(sv_reduced)
    elif method == "transects":
        sv_indexed = aggregate_transects(sv_reduced)
    else:
        raise ValueError(
            f"The defined aggregation method ('{method}') is invalid. This method must either be "
            f"one of the folmining: 'cells', 'intervals', 'transects'."
        ).with_traceback(None)

    # Return the dataset
    return sv_indexed

def ingest_echoview_sv(
    sv_path: Path,
    center_frequencies: Optional[Dict[str, float]] = None,
    transect_pattern: Optional[str] = None,    
    aggregate_method: Literal["cells", "interval", "transect"] = "cells",
    impute_coordinates: bool = True
):
    # Validate directory existence
    if not sv_path.exists():
        raise FileNotFoundError(f"The export file directory ({sv_path.as_posix()}) not found!")

    # Validate files existence
    if not sv_path.iterdir():
        raise FileNotFoundError(
            f"The export file directory ({sv_path.as_posix()}) contains no files!"
        )

    # Update the units for `center_frequencies` to match expected values from Echoview
    # ---- Hz -> kHz
    center_frequencies = {freq * 1e-3: value for freq, value in center_frequencies.items()}

    # Get the target inversion files
    sv_filepaths = {"cells": [p for p in sv_path.rglob("*.csv") if p.is_file()]}

    # Get the trasnect numbers
    if transect_pattern:
        transect_num_df = ingest_nasc.map_transect_num(sv_filepaths, transect_pattern)    
    else:
        # ---- valueize DataFrame
        transect_num_df = pd.DataFrame({
            "file_type": ["cells"] * len(sv_filepaths["cells"]),
            "file_path": sv_filepaths["cells"],
            "transect_num": None
        })

    # Concatenate the files
    sv = pd.concat(
        echoview_sv_to_df(
            transect_num_df, impute_coordinates,
        )
    )

    # Sort
    ingest_nasc.sort_echoview_export_df(sv, inplace=True) 
    
    # Keep target transmit frequencies
    if center_frequencies:
        sv_subset = sv[sv["frequency"].isin(center_frequencies)]
    else:
        sv_subset = sv.copy()
        center_frequencies = {
            freq: {"min": -999., "max": 999.} for freq in sv_subset["frequency"].unique()
        }

    # Integrate the backscatter based on the defined aggregation method
    sv_integrated = integrate_measurements(data=sv_subset, 
                                           method=aggregate_method, 
                                           sv_thresholds=center_frequencies)

    # Return
    return sv_integrated


def dict_to_Parameters(params_dict: Dict[str, Any]):
    params = Parameters()
    for name, kwargs in params_dict.items():
        params.add(
            name,
            value=kwargs["value"],
            min=kwargs.get("min", 0.),
            max=kwargs.get("max", np.inf),
            vary=kwargs.get("vary", True)
        )
    return params

def get_parameter_limits(
    scattering_parameters: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Extract the lower and upper bounds for scattering model parameters
    """

    # Get the upper and lower values in case of rescaling
    parameter_limits = {
        key: {"min": value["min"], "max": value["max"]}
        for key, value in scattering_parameters.items()
    }

    return parameter_limits


def minmax_normalize(
    parameter_sets: Dict[str, Any],
    inverse: bool = False,
    inverse_reference: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Normalize the optimization parameters

    Parameters
    ----------
    parameter_sets: Dict[str, Any]
        Dictionary comprising acoustic scattering model parameters that are to be optimized.
    inverse: bool, default: False
        Boolean flag to indicate whether the parameters should be normalized via min-max
        normalization (`inverse=True`) or denormalized (i.e. inverse normalization,
        `inverse=False`).
    """

    # Min-max normalization
    if inverse is False:
        return {
            realization: {
                key: (
                    {
                        **kwargs,
                        "value": (
                            (kwargs["value"] - kwargs["min"]) / (kwargs["max"] - kwargs["min"])
                            if kwargs["max"] != kwargs["min"]
                            else 0.0
                        ),
                        "min": 0.0,
                        "max": 1.0,
                    }
                    if isinstance(kwargs, dict)
                    else kwargs
                )  # Skip scalar entries
                for key, kwargs in sets.items()
            }
            for realization, sets in parameter_sets.items()
        }
    # Min-max inverse normalization
    else:
        return {
            key: (
                kwargs * (inverse_reference[key]["max"] - inverse_reference[key]["min"])
                + inverse_reference[key]["min"]
            )
            for key, kwargs in parameter_sets.items()
        }

def monte_carlo_initialization(
    scattering_parameters: Dict[str, Any], 
    n_realizations: int, 
) -> Dict[str, Any]:
    """
    Monte Carlo simulation of initial values for scattering parameters
    """

    # Create parameter sets for the defined number of realizations
    parameter_sets = dict(
        map(
            lambda i: (
                i,
                {
                    key: (
                        {
                            "value": (
                                np.random.uniform(values["min"], values["max"])
                                if values["vary"] else values["value"]
                            ),
                            "min": values["min"],
                            "max": values["max"],
                        }

                    )
                    for key, values in scattering_parameters[0].items()
                },
            ),
            range(1, n_realizations + 1),
        )
    )
    
    # Return the parameter sets
    return {**scattering_parameters, **parameter_sets}

def inverse_normalize_series(series: pd.Series, parameter_bounds: Dict[str, Any]):
    """
    Apply inverse min-max normalization to a `pandas.Series`
    """

    # Create copy
    series_copy = series.copy()

    # Get the columns that will be inverse transformed
    transform_cols = list(parameter_bounds.keys())

    # Convert `ranges_dict` to a `pd.DataFrame`
    ranges_df = pd.DataFrame(parameter_bounds).T

    # Apply inverse transformation
    series_copy[transform_cols] = (
        series_copy[transform_cols] * (ranges_df["max"] - ranges_df["min"]) + ranges_df["min"]
    )

    # Return the inverse min-max normalized series
    return series_copy

def uniformly_bent_cylinder(
    n_segments: Union[int, awk.Array],
    radius_of_curvature_ratio: float,
    taper_order: float,
) -> Tuple[awk.Array, awk.Array, awk.Array, awk.Array, awk.Array]:
    """
    Generates the normalized position matrix for an uniformly bent cylinder
    """

    # Curvature coefficients (for shape normalization)
    # ---- Gamma
    gamma = 0.5 / radius_of_curvature_ratio
    # ---- Normalization
    norm_ratio = radius_of_curvature_ratio * 2

    # Create normalized horizontal increments along the anterioposterior (z) axis of the body shape
    # ---- If array of values
    if isinstance(n_segments, awk.Array):
        z = awk.Array([np.linspace(-1.0, 1.0, n) for n in n_segments])
    # ---- If only a single value
    else:
        z = awk.Array(np.linspace(-1.0, 1.0, n_segments))
        
    # Compute the taper vector
    taper = np.sqrt(1 - z**taper_order)
    
    # Bend the cylinder
    # ---- z-axis
    z_curved = np.sin(gamma) * z
    # ---- Dorsoventral axis (x-axis)
    x_curved = 1 - np.sqrt(1 - z_curved**2)
    
    # Normalize the curvature
    # ---- z-axis
    z_norm = z_curved * norm_ratio
    # ---- x-axis
    x_norm = x_curved * norm_ratio
    
    # Calculate the slope between curved segments
    gamma_tilt = np.arctan2(z_norm, x_norm)

    # Calculate the orientation angles between curved segments
    # ---- z-axis differences
    dz = z_norm[:, 1:] - z_norm[:, :-1]
    # ---- x-axis differences
    dx = x_norm[:, 1:] - x_norm[:, :-1] + np.finfo(float).eps
    # ---- alpha tilt angles
    alpha_tilt = np.arctan(dz / dx)
    # ---- Compute final element
    alpha_tilt_end = awk.Array([[v] for v in np.arctan(dz[:, -1] / dx[:, -1])])
    # ---- Concatenate 
    alpha_tilts = awk.concatenate([alpha_tilt, alpha_tilt_end], axis=1)
    # ---- beta tilt angles
    beta_tilt = awk.where(alpha_tilts >= 0.0, alpha_tilts - np.pi / 2, alpha_tilts + np.pi / 2)

    # Compute the along-axis Euclidean distances to construct the position vector
    # ---- Position vector
    r_pos = np.sqrt(x_norm**2 + z_norm**2)
    # ---- Compute the first derivative
    if r_pos.ndim > 1:
        # ---- Get the first derivative of the initial index
        pos1 = awk.Array([[v] for v in np.sqrt(dx[:, 0] * dx[:, 0] + dz[:, 0] * dz[:, 0])])
        # ---- Prepend to the rest of the derivatives
        dr_pos = awk.concatenate([pos1, np.sqrt(dx * dx + dz * dz)], axis=1)
    else:
        # ---- Get the first derivative of the initial index
        pos1 = awk.Array([np.sqrt(dx[0] * dx[0] + dz[0] * dz[0])])
        ## ---- Prepend to the remaining derivatives
        dr_pos = awk.concatenate([pos1, np.sqrt(dx * dx + dz * dz)], axis=0)
        
    # Return the relevant parameters
    return taper, gamma_tilt, beta_tilt, r_pos, dr_pos

def generate_frequency_interval(
    frequency: np.ndarray[float], length_sd_norm: float, frequency_interval: float
) -> awk.Array:
    """
    Generate frequency interval 2D array centered on an array input of center frequencies.
    """

    frequency_lst = [
        np.arange(
            freq * (1 - 3.1 * length_sd_norm),
            freq * (1 + 3.1 * length_sd_norm) + frequency_interval,
            frequency_interval,
        )
        for freq in frequency
    ]
    
    # Return the frequency list
    return awk.Array(frequency_lst)

def wavenumber(
    frequency: Union[np.ndarray, float],
    sound_speed_sw: float,
) -> np.ndarray[float]:
    """
    Compute the acoustic wavenumber
    """

    return 2 * np.pi * frequency / sound_speed_sw

def reflection_coefficient(
    g: Union[np.ndarray, float],
    h: Union[np.ndarray, float],
) -> np.ndarray[float]:
    """
    Compute the reflection coefficient based on material properties
    """

    return (1 - g * h * h) / (g * h * h) - (g - 1) / g

def _compute_single_frequency(ka_tapered_i, delta_gamma_cos_i, delta_theta_cos_i,
                                       M1_i, M2_i, ka_norm_i, J1_i, n_k_i, n_segments_i, n_theta):
    """Vectorized computation for single frequency"""
    
    # Convert awkward arrays to numpy arrays
    ka_tapered_np = np.array(ka_tapered_i)
    delta_gamma_cos_np = np.array(delta_gamma_cos_i)
    delta_theta_cos_np = np.array(delta_theta_cos_i)
    M1_np = np.array(M1_i)
    M2_np = np.array(M2_i)
    ka_norm_np = np.array(ka_norm_i)
    J1_np = np.array(J1_i)
    
    # ---- Broadcast `ka_tapered` for multiplication with `delta_theta_cos`
    ARG = (
        2 * ka_tapered_np[:, :, np.newaxis] * delta_theta_cos_np[np.newaxis, :, :]
        + np.finfo(float).eps
    )
    
    # -------- Flatten for subsequent interpolation
    ARG_flat = ARG.ravel(order="F").reshape((n_k_i * n_segments_i, n_theta), order="F")
    
    # ---- Vectorized interpolation
    J1_interp = np.column_stack([
        np.interp(ARG_flat[:, m], ka_norm_np, J1_np) 
        for m in range(n_theta)
    ])
    
    # -------- Reshape and normalize
    J1_norm = (J1_interp / ARG_flat).reshape((n_k_i, n_segments_i, n_theta), order="F")
    
    # ---- Exponentiate the terms to compute the phase
    phase = np.exp(1j * M1_np[:, :, np.newaxis] * delta_gamma_cos_np[np.newaxis, :, :])
    
    # ---- Combine the adjusted `ka`, interpolated Bessel function output, and phase
    M3 = ka_tapered_np[:, :, np.newaxis] ** 2 * J1_norm * phase
    
    # ---- Compute the complex form function (matrix multiplication via Einstein summation)
    return [np.einsum("ijk, j->ik", M3, M2_np) + np.finfo(float).eps]

def orientation_average(
    angle: np.ndarray[float],
    form_function: np.ndarray[complex],
    theta_mean: float,
    theta_sd: float,
    distribution: Literal["gaussian", "uniform"] = "gaussian",
) -> np.ndarray[float]:
    """
    Compute the orientation-averaged linear backscattering cross-section :math:`\sigma_{bs}(\theta)`
    """

    # Weight based on distribution input
    # ---- Gaussian (Normal)
    if distribution == "gaussian":
        # ---- Get interval
        orientation_interval = np.diff(angle).mean()
        # ---- Compute the PDF
        PDF = (
            orientation_interval
            * np.exp(-0.5 * (angle - theta_mean) ** 2 / theta_sd**2)
            / (np.sqrt(2 * np.pi) * theta_sd)
        )
    # ---- Uniform
    elif distribution == "uniform":
        # ---- Compute the PDF
        PDF = np.ones(len(form_function)) / len(form_function)
    else:
        raise ValueError("Invalid distribution type. Choose 'gaussian' or 'uniform'.")

    # Return the weighted form function
    # ---- If complex
    return awk.Array([
        np.sqrt(np.matmul((np.array(form_function[i]).real ** 2 + 
                            np.array(form_function[i]).imag ** 2), PDF))
        for i in range(len(form_function))
    ])

def length_average(
    length_values: np.ndarray[float],
    ka_f: np.ndarray[float],
    ka_c: np.ndarray[float],
    form_function: np.ndarray[complex],
    length_mean: float,
    length_deviation: float,
    distribution: Literal["gaussian", "uniform"] = "gaussian",
) -> np.ndarray[float]:
    """
    Compute the length-averaged linear backscattering cross-section (:math:`\sigma_{bs}(L)`)
    """
    
    # Normalize the length values, if needed
    length_norm = length_values / length_mean
    # ---- Also normalize the standard deviation
    length_sd_norm = length_deviation / length_mean

    # Weight based on distribution input
    # ---- Gaussian (Normal)
    if distribution == "gaussian":
        # ---- Get the interval
        length_interval = np.diff(length_norm).mean()
        # ---- Compute the PDF
        PDF = (
            length_interval
            * np.exp(-0.5 * (length_norm - 1) ** 2 / length_sd_norm**2)
            / (np.sqrt(2 * np.pi) * length_sd_norm)
        )
    # ---- Uniform
    elif distribution == "uniform":
        # ---- Compute the PDF
        PDF = np.ones(len(form_function)) / len(form_function)
    else:
        raise ValueError("Invalid distribution type. Choose 'gaussian' or 'uniform'.")

    # Vectorized computation - compute length-weighted ka for all frequencies
    ka_weighted = length_norm * ka_c.reshape(-1, 1)

    # Vectorized computation using awkward arrays
    sigma_bs_length = awk.Array([
        (
            length_norm**2
            * PDF
            * np.interp(ka_weighted[i], np.array(ka_f[i]), np.array(form_function[i]) ** 2)
        ).sum()
        for i in range(len(form_function))
    ])

    # Return the weighted average
    return sigma_bs_length

class InversionMatrix(InversionBase):
    """
    !!! DOCSTRING
    """
    def __new__(
        cls,      
        data: pd.DataFrame,   
        simulation_settings: Dict[str, Any],
    ):
        # Validate
        try:
            # ---- Check
            valid_args = ValidateInversionMatrix.create(
                **dict(data=data, simulation_settings=simulation_settings)
            )
        # Break creation
        except ValidationError as e:
            raise val.EchopopValidationError(str(e)) from None

        # Create instance
        self = super().__new__(cls)

        # Update attributes
        self.measurements = valid_args["data"].copy()
        self.simulation_settings = valid_args["simulation_settings"]
        
        # Generate
        return self
    
    def __init__(
        self, 
        data: pd.DataFrame,   
        simulation_settings: Dict[str, Any],
    ):
        
        # Set inversion method
        self.inversion_method = "scattering_model"

        # Initialize attributes
        self._sv_cache = {}  # Initialize cache
        self.parameter_bounds = {}
        self.model = None
        self.model_params = {}
        self.model_settings = {}
        self.rng = None

        # Store random number generator, if required
        self._set_rng(**self.simulation_settings)

    def _set_rng(
        self,
        monte_carlo: bool,
        mc_seed: Optional[int] = None,
        **kwargs
    ):
        if monte_carlo:
            self.rng = np.random.default_rng(mc_seed)

    def _set_monte_carlo(
        self,
        initial_parameters: Dict[str, Any],
        mc_realizations: int,
        **kwargs,
    ) -> Dict[str, Any]:

        # Create parameter sets for the defined number of realizations
        parameter_sets = dict(
            map(
                lambda i: (
                    i,
                    {
                        key: (
                            {
                                "value": (
                                    self.rng.uniform(values["min"], values["max"])
                                    if values["vary"] else values["value"]
                                ),
                                "min": values["min"],
                                "max": values["max"],
                            }

                        )
                        for key, values in initial_parameters[0].items()
                    },
                ),
                range(1, mc_realizations + 1),
            )
        )
        
        # Return the parameter sets
        return {**initial_parameters, **parameter_sets}
    
    def _set_minimizer(
        self,
        Sv_measured: pd.Series,
        parameter_set: Dict[int, Any],
    ):
        
        # Extract the frequencies from the index
        center_frequencies = np.array(Sv_measured.index.values, dtype=float)
        
        # Generate 'n' realizations from Monte Carlo method
        if self.simulation_settings["monte_carlo"]:
            realizations = self._set_monte_carlo(
                initial_parameters=parameter_set, 
                **self.simulation_settings
            )
            
        # Convert model parameters to the correct `lmfit.Parameters` class object
        parameter_realizations = {r: dict_to_Parameters(p) for r, p in realizations.items()}
        
        # Find which values are below the defined threshold
        valid_idx = np.argwhere(Sv_measured > -999.).flatten()
        
        # Convert to a `numpy.ndarray`
        Sv_measured = Sv_measured.to_numpy()[valid_idx]
        
        # Generate `Minimizer` function class required for bounded optimization
        # ---- This will generate per realization within a List
        return [
            Minimizer(
                self._objective,
                parameter_realizations[realization],
                fcn_args=(
                    Sv_measured,
                    center_frequencies[valid_idx],
                ),
                nan_policy="omit",
            )
            for realization in parameter_realizations
        ]
        
    def _objective(self, 
                   parameter_set: Parameters, 
                   Sv_measured: np.ndarray[float],
                   center_frequencies: np.ndarray[float],
                   **kwargs) -> float:
        
        # Compute Sv fit using the incoming parameter set
        Sv_prediction = self._predict_Sv(parameter_set=parameter_set, 
                                         center_frequencies=center_frequencies,
                                         **kwargs)

        # Pre-allocate weight array
        wd = np.ones(len(Sv_measured))

        # Compute deviation
        deviation = Sv_prediction - Sv_measured

        # Return the summed absolute deviation (Q)
        return np.sum(np.abs(deviation) * wd)

    def _predict_Sv(
        self,
        parameter_set: Parameters,
        center_frequencies: np.ndarray[float],
        **model_kwargs,
    ):
        
        # Extract parameter values from dictionary for parsing
        parameters_dict = parameter_set.valuesdict()
        
        # Inverse normalization, if needed
        if self.simulation_settings.get("scale_parameters", False):
            parameters_dict = minmax_normalize(
                parameters_dict, 
                inverse=True,
                inverse_reference=self.parameter_bounds
            )
            
        # Run scattering model with the new parameterization to get the linear backscattering 
        # coefficient (f_bs)
        Sv_prediction = self.model(center_frequencies=center_frequencies,
                                   **parameters_dict, 
                                   **self.model_settings, 
                                   **self.simulation_settings["environment"])   
        
        return Sv_prediction     
        
    def _set_minimizers(
        self,
        parameter_set: Dict[int, Any],      
    ):
        
        # Definine a new column for the `lmfit.Minimizer` class
        self.measurements["minimizer"] = np.array(np.nan).astype(object)
        
        # Create list of `lmfit.Minimizer` objects for all realizations
        self.measurements["minimizer"] = self.measurements["sv_mean"].apply(
            self._set_minimizer,
            axis=1,
            args=(
                parameter_set,
            )
        )

        # Define label for verbosity
        self.measurements["label"] = [
            "; ".join(
                f"{name if name is not None else 'index'}: {val}"
                for name, val in zip(self.measurements.index.names, 
                                     (idx if isinstance(idx, tuple) else (idx,)))
            )
            for idx in self.measurements.index
        ]
        
    def build_scattering_model(
        self,
        model_parameters: Dict[str, Any],
        model_settings: Dict[str, Any],
    ) -> None:
        
        # Validate
        try:
            # ---- Check
            valid_args = ValidateBuildModelArgs.create(
                **dict(
                    model_parameters=model_parameters,
                    model_settings=model_settings
                )
            )
        # Break creation
        except (ValidationError, Exception) as e:
            raise val.EchopopValidationError(str(e)) from None
        
        # Update attributes
        self.model_params.update(valid_args["model_parameters"])
        self.model_settings.update(valid_args["model_settings"])
        
        # Retrieve the scattering model
        self.model = SCATTERING_MODEL_PARAMETERS[self.model_settings["type"]].get("function")
        
        # Get the parameter boundaries (required for inverting minmax normalization if applied)
        self.parameter_bounds.update(get_parameter_limits(valid_args["model_parameters"]))
        
        # Initialize the parameter sets dictionary
        parameter_set = {0: self.model_params}
        
        # Apply minmax normalization, if set
        if self.simulation_settings["scale_parameters"]:
            parameter_set = minmax_normalize(
                parameter_sets=parameter_set, 
            )

        # Build the minimizers and labels
        self._set_minimizers(parameter_set)

    def _optim(
        self,
        Sv_measured: pd.Series,
        verbose: bool = True,
        **kwargs
    ):

        # Catch start time in case `verbose=True`
        start_time = time.time()

        # Find which values are below the defined threshold
        valid_idx = np.argwhere(Sv_measured["sv_mean"] > -999.).flatten()
        center_frequencies = np.array(Sv_measured["sv_mean"].index.values[valid_idx], dtype=float)

        # Only run if the correct number of frequencies are valid
        if len(valid_idx) >= self.simulation_settings["minimum_frequency_count"]:

            # Assign message string
            frequency_msg = ""

            # Convert to a numpy array
            Sv = Sv_measured["sv_mean"].to_numpy()[valid_idx].astype(float)

            # Optimize over all realizations [initialize]
            # ---- Fit errors
            fit_errors = np.ones(self.simulation_settings["mc_realizations"]) * np.nan
            # ---- Parameter sets (column names/structure)
            parameter_fits = (
                pd.DataFrame(
                    Sv_measured.minimizer.iloc[0][0].params.valuesdict(),
                    index=range(self.simulation_settings["mc_realizations"]),
                )
                * np.nan
            )
            
            # Iterate through the realizations
            for realization in range(self.simulation_settings["mc_realizations"]):
                with warnings.catch_warnings():
                    warnings.filterwarnings(action="ignore", category=RuntimeWarning)
                    # ---- Optimize
                    minimizer = Sv_measured["minimizer"].iloc[0][realization]
                    parameters_optimized = minimizer.minimize(
                        method="least_squares",
                        **self.optimization_kwargs
                    )
                    
                    # best_fit_Sv = self._predict_Sv(
                    #     parameter_set=parameters_optimized.params,
                    #     center_frequencies=center_frequencies
                    # )
                    fit_errors[realization] = self._objective(
                        parameters_optimized.params, 
                        Sv, 
                        center_frequencies
                    )
                    parameter_fits.loc[realization] = pd.Series(parameters_optimized.params.valuesdict())
                    print(realization)

            # Find the parameter set with the lowest Q
            best_fit_set = parameter_fits.loc[np.nanargmin(fit_errors)]

            # Add error, Q, to series
            best_fit_set["Q"] = np.nanmin(fit_errors)
        else:
            # Assign message string
            frequency_msg = (
                f"\nWARNING: The number of frequencies with valid Sv [{len(valid_idx)}] "
                f"was fewer than the minimum frequency count "
                f"[{self.simulation_settings["minimum_frequency_count"]}]. Values were not optimized."
            )

            # Create `pandas.Series`
            best_fit_set = pd.Series(Sv_measured.minimizer.iloc[0][0].params.valuesdict()) * np.nan

            # Add error, Q, to series
            best_fit_set["Q"] = np.nan

        # Inverse transformation if values are scaled
        if self.simulation_settings["scale_parameters"]:
            best_fit_set = inverse_normalize_series(
                best_fit_set, self.parameter_bounds
            )

        # Catch end time in case `verbose=True`
        end_time = time.time()

        # Print results
        if verbose:
            # ---- Row label
            row = f"{Sv_measured['label']}"
            # ---- Get error value
            error_value = f" Sv error (Q): {np.round(best_fit_set['Q'], 3)} dB "
            # ---- Parameter values
            parameter_values = f"\n{best_fit_set[:-1].to_frame().T}"
            # ---- Get elapsed time (s)
            elapsed_time = f" Elapsed time: {np.round(end_time - start_time, 2)} s;"
            # ---- Number of frequencies
            valid_freq = (
                "["
                + "/".join(f"{freq * 1e-3}" for freq in center_frequencies[valid_idx])
                + " kHz"
                + "]"
            )

            # Print out
            print(row + elapsed_time + error_value + valid_freq + frequency_msg + parameter_values)

    def invert():
        pass