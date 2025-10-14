from echopop.scratch.scratch_tests import *
import numpy as np
import pandas as pd
import time
import warnings
import hashlib
from scipy.special import j1
import awkward as awk
from concurrent.futures import ThreadPoolExecutor

try:
    from numba import njit
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    def njit(*args, **kwargs):
        def decorator(func):
            return func
        return decorator

# ==============================================================================
# OPTIMIZED HELPER FUNCTIONS
# ==============================================================================

def ragged_rowwise_apply(arr, fn):
    """Apply function to each row of ragged array efficiently"""
    return awk.Array([fn(row) for row in awk.to_list(arr)])

def vectorized_bessel_interp(ARG_ragged, ka_norm_ragged, J1_ragged, n_theta, n_k_array, n_segments_array):
    """Optimized Bessel interpolation with correct dimension handling"""
    def _fast_bessel_interp(ARG_flat, ka_norm_np, J1_np, n_theta):
        """Vectorized interpolation with minimal memory allocations"""
        n_points = ARG_flat.shape[0]
        
        # Pre-allocate result array
        result = np.empty((n_points, n_theta), dtype=np.float64)
        
        # For small n_theta, unroll loop for better cache performance  
        if n_theta == 1:
            result[:, 0] = np.interp(ARG_flat[:, 0], ka_norm_np, J1_np)
        elif n_theta == 2:
            result[:, 0] = np.interp(ARG_flat[:, 0], ka_norm_np, J1_np)
            result[:, 1] = np.interp(ARG_flat[:, 1], ka_norm_np, J1_np)
        elif n_theta <= 8:
            # Small theta: vectorized loop is optimal
            for th in range(n_theta):
                result[:, th] = np.interp(ARG_flat[:, th], ka_norm_np, J1_np)
        else:
            # Large theta: use numpy broadcasting tricks
            result = np.array([np.interp(ARG_flat[:, th], ka_norm_np, J1_np) for th in range(n_theta)]).T
        
        return result
    
    return awk.Array([
        _fast_bessel_interp(
            np.array(ARG_ragged[i]).ravel(order="F").reshape((n_k_array[i] * n_segments_array[i], n_theta), order="F"),
            np.array(ka_norm_ragged[i]), 
            np.array(J1_ragged[i]), 
            n_theta
        )
        for i in range(len(ARG_ragged))
    ])

def fast_einsum_reduction(arr_3d_ragged, weights_ragged):
    """Fast Einstein summation for ragged 3D arrays"""
    return ragged_rowwise_apply(
        awk.zip(arr_3d_ragged, weights_ragged),
        lambda args: np.einsum('ijk,j->ik', args[0], args[1]) + np.finfo(float).eps
    )

def batch_bessel_j1(ka_norm_ragged):
    """Batch compute Bessel functions for ragged arrays"""
    return awk.Array([j1(ka_norm_array) for ka_norm_array in ka_norm_ragged])

def _original_bessel_interp(ARG_flat, ka_norm_np, J1_np, n_theta):
    """Original working Bessel interpolation - DO NOT CHANGE"""
    result = np.zeros((ARG_flat.shape[0], n_theta))
    
    # Use the exact original logic
    for th in range(n_theta):
        result[:, th] = np.interp(ARG_flat[:, th], ka_norm_np, J1_np)
    
    return result

def _normalize_bessel_result(J1_interp, ARG, shape):
    """Efficiently normalize a single Bessel interpolation result"""
    n_k, n_segments, n_theta = shape
    
    # Ensure we're working with numpy arrays
    ARG_np = np.array(ARG) if not isinstance(ARG, np.ndarray) else ARG
    J1_interp_np = np.array(J1_interp) if not isinstance(J1_interp, np.ndarray) else J1_interp
    
    # Flatten and reshape ARG
    ARG_flat = ARG_np.ravel(order="F").reshape((n_k * n_segments, n_theta), order="F")
    
    # Normalize
    normalized = J1_interp_np / ARG_flat
    return normalized.reshape(shape, order="F")

def compute_phase_efficiently(M1_ragged, delta_gamma_cos_ragged):
    """Optimized phase computation avoiding redundant array conversions"""
    return awk.Array([
        np.exp(1j * np.array(M1_ragged[i][:, :, np.newaxis]) * 
               np.array(delta_gamma_cos_ragged[i][np.newaxis, :, :]))
        for i in range(len(M1_ragged))
    ])

# ==============================================================================
# OPTIMIZED EXTERNAL PCDWBA FUNCTION
# ==============================================================================
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
    """
    Optimized pcdwba with ONLY the fast Bessel interpolation improvement
    """
    # Use pre-computed constants
    sqrt_2pi = np.sqrt(2 * np.pi)
    eps = np.finfo(float).eps

    # Pre-allocate arrays based on input size
    n_freqs = len(center_frequencies)
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

    # Calculate product of tapered ka_f and along-axis tilts - OPTIMIZED
    ARG = 2 * awk.Array([
        np.array(ka_f_tapered[i][:, :, np.newaxis]) * 
        np.array(delta_theta_cos[i][np.newaxis, :, :])
        for i in range(len(n_k))
    ]) + eps

    # Use the original working Bessel interpolation approach
    J1_interp = awk.Array([
        _original_bessel_interp(
            np.array(ARG[i]).ravel(order="F").reshape((n_k[i] * n_segments[i], n_theta), order="F"),
            np.array(ka_f_norm[i]), np.array(J1[i]), n_theta
        )
        for i in range(len(n_k))
    ])

    # Normalize the interpolated Bessel function - ORIGINAL METHOD
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
    # ---- Get interval
    orientation_interval = np.diff(theta_values).mean()
    # ---- Compute the PDF
    PDF = (
        orientation_interval
        * np.exp(-0.5 * (theta_values - theta_mean) ** 2 / theta_sd**2)
        / (sqrt_2pi * theta_sd)
    )

    f_bs_orientation = awk.Array([
        np.sqrt(np.matmul((np.array(f_bs[i]).real ** 2 + 
                           np.array(f_bs[i]).imag ** 2), PDF))
        for i in range(len(f_bs))
    ])

    # Compute over vector lengths
    length_values = np.linspace(
        length_mean - 3 * (length_sd_norm * length_mean),
        length_mean + 3 * (length_sd_norm * length_mean),
        n_length,
    )

    # Generate length deviation statistics
    length_deviation = length_mean * length_sd_norm

    # Normalize by mean length
    length_norm = length_values / length_mean
    length_sd_norm_calc = length_deviation / length_mean

    # ---- Get interval
    length_interval = np.diff(length_norm).mean()
    # ---- Compute the PDF
    PDF_length = (
        length_interval
        * np.exp(-0.5 * (length_norm - 1) ** 2 / length_sd_norm_calc**2)
        / (sqrt_2pi * length_sd_norm_calc)
    )

    # Integrate length vector over the form function
    sigma_bs_length = awk.Array([
        np.sum(length_norm**2 * PDF_length * 
               np.interp(length_norm * ka_c[i], np.array(ka_f[i]), 
                         np.array(f_bs_orientation[i]) ** 2))
        for i in range(len(f_bs_orientation))
    ])

    # Final length and normalization
    sigma_bs = sigma_bs_length * length_mean**2
    
    # Convert to Sv
    Sv_prediction = 10 * np.log10(number_density * sigma_bs)
    return Sv_prediction


# ==============================================================================
# OPTIMIZED HELPER FUNCTIONS  
# ==============================================================================

def monte_carlo_initialization(
    scattering_parameters: Dict[str, Any], 
    n_realizations: int, 
) -> Dict[str, Any]:
    """
    Optimized Monte Carlo simulation of initial values for scattering parameters
    """
    # Get the base parameters
    base_params = scattering_parameters[0]
    
    # Pre-allocate result dictionary
    parameter_sets = {}
    
    # Vectorized parameter generation
    for i in range(1, n_realizations + 1):
        param_dict = {}
        for key, values in base_params.items():
            if values["vary"]:
                param_dict[key] = {
                    "value": np.random.uniform(values["min"], values["max"]),
                    "min": values["min"],
                    "max": values["max"],
                }
            else:
                param_dict[key] = {
                    "value": values["value"],
                    "min": values["min"],
                    "max": values["max"],
                }
        parameter_sets[i] = param_dict
    
    # Return the parameter sets combined with original
    return {**scattering_parameters, **parameter_sets}


# ==============================================================================
# OPTIMIZED INVERSION MATRIX CLASS
# ==============================================================================

class InversionMatrix(InversionBase):
    """
    Highly optimized inversion class with significant performance improvements:
    - Batch processing of minimizers
    - Vectorized operations where possible
    - Efficient caching strategies
    - Reduced memory allocations
    - Parallel processing support
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

        # Initialize attributes with performance optimizations
        self._prediction_cache = {}  # Cache for model predictions
        self._minimizer_cache = {}   # Cache for minimizer objects
        self.parameter_bounds = {}
        self.model = pcdwba  # Use the original pcdwba function with proven optimizations
        self.model_params = {}
        self.model_settings = {}
        self.optimization_kwargs = {}
        self.rng = None

        # Pre-allocate arrays for better memory management
        self._preallocated_arrays = {}

        # Store random number generator, if required
        self._set_rng(**self.simulation_settings)

    def _create_prediction_cache_key(self, parameter_set: Parameters, center_frequencies: np.ndarray) -> str:
        """Create efficient cache key for predictions"""
        param_hash = hashlib.md5(str(sorted(parameter_set.valuesdict().items())).encode()).hexdigest()
        freq_hash = hashlib.md5(center_frequencies.tobytes()).hexdigest()
        return f"{param_hash}_{freq_hash}"

    def _batch_create_minimizers(self, parameter_realizations: Dict[int, Any], Sv_measured: np.ndarray, center_frequencies: np.ndarray) -> list:
        """Create all minimizers at once for better efficiency"""
        return [
            Minimizer(
                self._objective,
                dict_to_Parameters(parameter_realizations[realization]),
                fcn_args=(Sv_measured, center_frequencies),
                nan_policy="omit",
            )
            for realization in parameter_realizations
        ]
        self._set_rng(**self.simulation_settings)

    def _cache_key(self, parameters_dict, center_frequencies):
        """Generate cache key without aggressive rounding to preserve precision"""
        key_data = tuple(sorted(parameters_dict.items())) + tuple(center_frequencies)
        return hashlib.md5(str(key_data).encode()).hexdigest()
    
    def _predict_Sv(
        self,
        parameter_set: Parameters,
        center_frequencies: np.ndarray[float],
        **model_kwargs,
    ):
        """Highly optimized Sv prediction with intelligent caching"""
        # Extract parameter values from dictionary for parsing
        parameters_dict = parameter_set.valuesdict()
        
        # Create efficient cache key
        cache_key = self._create_prediction_cache_key(parameter_set, center_frequencies)
        
        # Check cache first - this can provide massive speedups
        if cache_key in self._prediction_cache:
            return self._prediction_cache[cache_key]
        
        # Inverse normalization, if needed
        if self.simulation_settings.get("scale_parameters", False):
            parameters_dict = minmax_normalize(
                parameters_dict, 
                inverse=True,
                inverse_reference=self.parameter_bounds
            )
            
        # Run optimized scattering model
        Sv_prediction = self.model(
            center_frequencies=center_frequencies,
            **parameters_dict, 
            **self.model_settings, 
            **self.simulation_settings["environment"]
        )   
        
        # Cache the result for future use
        self._prediction_cache[cache_key] = Sv_prediction
        
        # Prevent cache from growing too large (memory management)
        if len(self._prediction_cache) > 1000:
            # Remove oldest entries (simple FIFO)
            oldest_keys = list(self._prediction_cache.keys())[:100]
            for old_key in oldest_keys:
                del self._prediction_cache[old_key]
        
        return Sv_prediction
        
        # Cache the result for future use
        self._sv_cache[cache_key] = Sv_prediction
        
        return Sv_prediction     
    
    def _objective(self, 
                   parameter_set: Parameters, 
                   Sv_measured: np.ndarray[float],
                   center_frequencies: np.ndarray[float],
                   **kwargs) -> float:
        """Objective function using L1 norm (sum of absolute deviations)"""
        
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
        """Create parameter sets for Monte Carlo realizations"""

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
        """Set up single minimizer for optimization"""
        
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
        
    def _set_minimizers(
        self,
        parameter_set: Dict[int, Any],      
    ):
        """Set up minimizers for all realizations"""
        
        # Define a new column for the `lmfit.Minimizer` class
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
        """Build the scattering model configuration with validation"""
        
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
        
        # Retrieve the optimized scattering model
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
        """Highly optimized version with parallel processing and advanced caching"""

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

            # Pre-allocate arrays for better memory performance
            n_realizations = self.simulation_settings["mc_realizations"]
            fit_errors = np.full(n_realizations, np.nan)
            
            # Get parameter structure once
            param_template = Sv_measured.minimizer.iloc[0][0].params.valuesdict()
            parameter_fits = pd.DataFrame(
                index=range(n_realizations), 
                columns=param_template.keys(),
                dtype=float
            )
            
            # Get minimizers for all realizations
            minimizers = Sv_measured["minimizer"].iloc[0]
            
            for realization in range(n_realizations):
                with warnings.catch_warnings():
                    warnings.filterwarnings(action="ignore", category=RuntimeWarning)
                    
                    # Get minimizer for this realization
                    minimizer = minimizers[realization]
                    
                    # Optimize with pre-configured settings
                    try:
                        parameters_optimized = minimizer.minimize(
                            method="least_squares",
                            **self.optimization_kwargs
                        )
                        
                        # Store results efficiently
                        fit_errors[realization] = self._objective(
                            parameters_optimized.params, 
                            Sv, 
                            center_frequencies
                        )
                        
                        # Use vectorized assignment
                        parameter_fits.iloc[realization] = list(parameters_optimized.params.valuesdict().values())
                        
                    except Exception as e:
                        if verbose:
                            print(f"Realization {realization} failed: {str(e)}")
                        continue
                        
                    if verbose and realization % 10 == 0:  # Reduce print frequency for performance
                        print(f"Completed {realization + 1}/{n_realizations} realizations")

            # Find best parameters efficiently
            best_idx = np.nanargmin(fit_errors)
            best_fit_set = parameter_fits.iloc[best_idx].copy()
            best_fit_set["Q"] = fit_errors[best_idx]
            
        else:
            # Assign message string
            frequency_msg = (
                f"\nWARNING: The number of frequencies with valid Sv [{len(valid_idx)}] "
                f"was fewer than the minimum frequency count "
                f"[{self.simulation_settings['minimum_frequency_count']}]. Values were not optimized."
            )

            # Create `pandas.Series` efficiently
            param_template = Sv_measured.minimizer.iloc[0][0].params.valuesdict()
            best_fit_set = pd.Series(param_template) * np.nan
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
                + "/".join(f"{freq * 1e-3}" for freq in center_frequencies)
                + " kHz"
                + "]"
            )

            # Print out
            print(row + elapsed_time + error_value + valid_freq + frequency_msg + parameter_values)

        return best_fit_set

    def _parallel_optim_realization(self, args):
        """Single realization optimization for parallel processing"""
        realization, minimizer, Sv, center_frequencies = args
        
        with warnings.catch_warnings():
            warnings.filterwarnings(action="ignore", category=RuntimeWarning)
            
            try:
                parameters_optimized = minimizer.minimize(
                    method="least_squares",
                    **self.optimization_kwargs
                )
                
                fit_error = self._objective(
                    parameters_optimized.params, 
                    Sv, 
                    center_frequencies
                )
                
                param_values = list(parameters_optimized.params.valuesdict().values())
                
                return realization, fit_error, param_values
                
            except Exception as e:
                return realization, np.nan, [np.nan] * len(minimizer.params)

    def _optim_parallel(
        self,
        Sv_measured: pd.Series,
        verbose: bool = True,
        n_workers: int = 4,
        **kwargs
    ):
        """Parallel version of optimization for significant speedup"""
        
        start_time = time.time()
        
        # Find valid indices
        valid_idx = np.argwhere(Sv_measured["sv_mean"] > -999.).flatten()
        center_frequencies = np.array(Sv_measured["sv_mean"].index.values[valid_idx], dtype=float)
        
        if len(valid_idx) >= self.simulation_settings["minimum_frequency_count"]:
            
            Sv = Sv_measured["sv_mean"].to_numpy()[valid_idx].astype(float)
            n_realizations = self.simulation_settings["mc_realizations"]
            minimizers = Sv_measured["minimizer"].iloc[0]
            
            # Prepare arguments for parallel processing
            args_list = [
                (i, minimizers[i], Sv, center_frequencies)
                for i in range(n_realizations)
            ]
            
            # Use ThreadPoolExecutor for I/O bound operations
            with ThreadPoolExecutor(max_workers=n_workers) as executor:
                results = list(executor.map(self._parallel_optim_realization, args_list))
            
            # Process results
            fit_errors = np.full(n_realizations, np.nan)
            param_template = minimizers[0].params.valuesdict()
            parameter_fits = pd.DataFrame(
                index=range(n_realizations), 
                columns=param_template.keys(),
                dtype=float
            )
            
            for realization, fit_error, param_values in results:
                fit_errors[realization] = fit_error
                parameter_fits.iloc[realization] = param_values
            
            # Find best parameters
            best_idx = np.nanargmin(fit_errors)
            best_fit_set = parameter_fits.iloc[best_idx].copy()
            best_fit_set["Q"] = fit_errors[best_idx]
            
            frequency_msg = ""
            
        else:
            frequency_msg = (
                f"\nWARNING: The number of frequencies with valid Sv [{len(valid_idx)}] "
                f"was fewer than the minimum frequency count "
                f"[{self.simulation_settings['minimum_frequency_count']}]. Values were not optimized."
            )
            param_template = Sv_measured.minimizer.iloc[0][0].params.valuesdict()
            best_fit_set = pd.Series(param_template) * np.nan
            best_fit_set["Q"] = np.nan
        
        # Inverse transformation if needed
        if self.simulation_settings["scale_parameters"]:
            best_fit_set = inverse_normalize_series(
                best_fit_set, self.parameter_bounds
            )
        
        # Print results if verbose
        if verbose:
            end_time = time.time()
            elapsed_time = f" Elapsed time: {np.round(end_time - start_time, 2)} s;"
            row = f"{Sv_measured['label']}"
            error_value = f" Sv error (Q): {np.round(best_fit_set['Q'], 3)} dB "
            print(row + elapsed_time + error_value + frequency_msg)
        
        return best_fit_set

    def invert(self, optimization_kwargs: Dict[str, Any], verbose: bool = True, parallel: bool = False, n_workers: int = 4):
        """
        Main inversion method with optional parallel processing
        
        Parameters:
        -----------
        optimization_kwargs : Dict[str, Any]
            Optimization parameters for lmfit
        verbose : bool, default=True
            Whether to print progress information
        parallel : bool, default=False
            Whether to use parallel processing for Monte Carlo realizations
        n_workers : int, default=4
            Number of workers for parallel processing
        """
        # Store optimization kwargs
        self.optimization_kwargs = optimization_kwargs
        
        if verbose:
            mode_str = "parallel" if parallel else "serial"
            print(f"Starting inversion optimization in {mode_str} mode...")
        
        # Choose optimization method
        if parallel:
            # Apply parallel _optim to each row
            results = self.measurements.apply(
                lambda row: self._optim_parallel(row, verbose=verbose, n_workers=n_workers),
                axis=1
            )
        else:
            # Apply standard _optim to each row
            results = self.measurements.apply(
                lambda row: self._optim(row, verbose=verbose),
                axis=1
            )
        
        if verbose:
            print("Inversion completed!")
        
        return results

def uniformly_bent_cylinder_optimized(
    n_segments_array: awk.Array,
    radius_of_curvature_ratio: float,
    taper_order: float,
) -> Tuple[awk.Array, awk.Array, awk.Array, awk.Array, awk.Array]:
    """
    Optimized version of uniformly_bent_cylinder with better performance for large arrays
    """
    # Pre-compute constants
    gamma = 0.5 / radius_of_curvature_ratio
    norm_ratio = radius_of_curvature_ratio * 2
    eps = np.finfo(float).eps
    
    # Process each segment count efficiently
    results = []
    for n_segments in n_segments_array:
        # Create z coordinates
        z = np.linspace(-1.0, 1.0, n_segments)
        
        # Compute taper
        taper_i = np.sqrt(1 - z**taper_order)
        
        # Bend the cylinder
        z_curved = np.sin(gamma) * z
        x_curved = 1 - np.sqrt(1 - z_curved**2)
        
        # Normalize curvature
        z_norm = z_curved * norm_ratio
        x_norm = x_curved * norm_ratio
        
        # Calculate slope
        gamma_tilt_i = np.arctan2(z_norm, x_norm)
        
        # Calculate differences efficiently
        dz = np.diff(z_norm)
        dx = np.diff(x_norm) + eps
        
        # Calculate tilt angles
        alpha_tilt = np.arctan(dz / dx)
        alpha_tilt_final = np.arctan(dz[-1] / dx[-1])
        alpha_tilts = np.concatenate([alpha_tilt, [alpha_tilt_final]])
        
        # Beta tilt angles
        beta_tilt_i = np.where(alpha_tilts >= 0.0, alpha_tilts - np.pi / 2, alpha_tilts + np.pi / 2)
        
        # Position vectors
        r_pos_i = np.sqrt(x_norm**2 + z_norm**2)
        
        # First derivatives
        dr_first = np.sqrt(dx[0]**2 + dz[0]**2)
        dr_rest = np.sqrt(dx**2 + dz**2)
        dr_pos_i = np.concatenate([[dr_first], dr_rest])
        
        results.append((taper_i, gamma_tilt_i, beta_tilt_i, r_pos_i, dr_pos_i))
    
    # Convert to awkward arrays
    taper = awk.Array([result[0] for result in results])
    gamma_tilt = awk.Array([result[1] for result in results])
    beta_tilt = awk.Array([result[2] for result in results])
    r_pos = awk.Array([result[3] for result in results])
    dr_pos = awk.Array([result[4] for result in results])
    
    return taper, gamma_tilt, beta_tilt, r_pos, dr_pos

def pcdwba_ultra_optimized(
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
    """
    Ultra-optimized PCDWBA with minimal awkward array operations and maximum numpy vectorization
    """
    # Pre-compute constants once
    sqrt_2pi = np.sqrt(2 * np.pi)
    eps = np.finfo(float).eps
    pi_over_180 = np.pi / 180.0
    
    n_freqs = len(center_frequencies)
    n_theta = orientation_distribution["bins"]
    n_length = length_distribution["bins"]
    
    # Pre-compute orientation arrays - these are used heavily
    theta_values = np.linspace(theta_mean - 3.1 * theta_sd, theta_mean + 3.1 * theta_sd, n_theta)
    cos_theta = np.cos(theta_values * pi_over_180)
    sin_theta = np.sin(theta_values * pi_over_180)
    orientation_interval = np.diff(theta_values).mean()
    
    # Pre-compute reflection coefficient
    Gamma = reflection_coefficient(g, h)
    
    # Generate frequencies and wavenumbers with minimal awkward operations
    frequencies = generate_frequency_interval(center_frequencies, length_sd_norm, frequency_interval)
    k_c = wavenumber(center_frequencies, sound_speed_sw)
    ka_c = k_c * length_mean / length_radius_ratio
    
    # Process each frequency independently to avoid large awkward arrays
    Sv_results = []
    
    for freq_idx in range(n_freqs):
        freq_array = frequencies[freq_idx]
        k_f_i = wavenumber(freq_array, sound_speed_sw)
        ka_f_i = k_f_i * length_mean / length_radius_ratio
        
        n_k_i = len(ka_f_i)
        kL_max_i = np.nanmax(k_f_i) * length_mean * (1 + 3.1 * length_sd_norm)
        n_int_i = int(max(n_integration, np.ceil(kL_max_i * n_wavelength / (2 * np.pi))))
        
        # Generate shape for this frequency only
        taper_i, gamma_tilt_i, beta_tilt_i, r_pos_i, dr_pos_i = uniformly_bent_cylinder(
            awk.Array([n_int_i]), radius_of_curvature_ratio, taper_order
        )
        
        # Convert to numpy arrays immediately
        taper_np = np.array(taper_i[0])
        gamma_tilt_np = np.array(gamma_tilt_i[0])
        beta_tilt_np = np.array(beta_tilt_i[0])
        r_pos_np = np.array(r_pos_i[0])
        dr_pos_np = np.array(dr_pos_i[0])
        
        # Vectorized calculations with pure numpy
        ka_f_tapered = ka_f_i[:, np.newaxis] * taper_np[np.newaxis, :]
        ka_f_norm = np.linspace(0, np.nanmax(ka_f_tapered), len(ka_f_i))
        
        # Batch compute Bessel functions once
        J1_vals = j1(ka_f_norm)
        
        # Vectorized delta calculations
        delta_gamma_cos = gamma_tilt_np[:, np.newaxis] * cos_theta[np.newaxis, :]
        delta_theta_cos = (beta_tilt_np[:, np.newaxis] * cos_theta[np.newaxis, :] +
                          r_pos_np[:, np.newaxis] * sin_theta[np.newaxis, :])
        
        # Phase calculation
        M1 = 2 * ka_f_i[:, np.newaxis] * delta_gamma_cos
        phase = np.exp(1j * M1)
        M2 = Gamma * dr_pos_np
        
        # Optimized Bessel interpolation for single frequency
        ARG = 2 * ka_f_tapered[:, :, np.newaxis] * delta_theta_cos[np.newaxis, :, :] + eps
        ARG_flat = ARG.ravel(order="F").reshape((n_k_i * len(taper_np), n_theta), order="F")
        
        # Fast interpolation with pre-allocated result
        J1_interp = np.zeros_like(ARG_flat)
        for th in range(n_theta):
            J1_interp[:, th] = np.interp(ARG_flat[:, th], ka_f_norm, J1_vals)
        
        # Normalize and reshape
        J1_norm = (J1_interp / ARG_flat).reshape((n_k_i, len(taper_np), n_theta), order="F")
        
        # Form function calculation
        f_j = ka_f_tapered[:, :, np.newaxis]**2 * J1_norm * phase[:, :, np.newaxis]
        
        # Einstein summation
        f_bs = np.einsum('ijk,j->ik', f_j, M2) + eps
        
        # Orientation averaging
        f_bs_orientation = orientation_interval * np.sum(np.abs(f_bs)**2, axis=1) / sqrt_2pi
        
        Sv_results.append(f_bs_orientation)
    
    # Length distribution (vectorized across all frequencies)
    length_values = np.linspace(
        length_mean * (1 - 3.1 * length_sd_norm),
        length_mean * (1 + 3.1 * length_sd_norm), 
        n_length
    )
    
    length_norm = length_values / length_mean
    length_interval = np.diff(length_norm).mean()
    
    # Pre-compute PDF
    PDF_length = (
        length_interval * 
        np.exp(-0.5 * (length_norm - 1)**2 / length_sd_norm**2) / 
        (sqrt_2pi * length_sd_norm)
    )
    
    # Final integration for all frequencies
    sigma_bs = np.zeros(n_freqs)
    for i in range(n_freqs):
        sigma_bs[i] = np.sum(
            length_norm**2 * PDF_length * 
            np.interp(length_norm * ka_c[i], np.arange(len(Sv_results[i])), Sv_results[i]**2)
        ) * length_mean**2
    
    Sv_prediction = 10 * np.log10(number_density * sigma_bs)
    ed = time.time()
    print(ed-st)
    return Sv_prediction
