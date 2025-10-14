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

import matplotlib.pyplot as plt
import numpy as np

from lmfit import Parameters, minimize, report_fit


def gauss(x, amp, cen, sigma):
    """Gaussian lineshape."""
    return amp * np.exp(-(x-cen)**2 / (2.*sigma**2))


def gauss_dataset(params, i, x):
    """Calculate Gaussian lineshape from parameters for data set."""
    amp = params[f'amp_{i+1}']
    cen = params[f'cen_{i+1}']
    sig = params[f'sig_{i+1}']
    return gauss(x, amp, cen, sig)


def objective(params, x, data):
    """Calculate total residual for fits of Gaussians to several data sets."""
    ndata, _ = data.shape
    resid = 0.0*data[:]

    # make residual per data set
    for i in range(ndata):
        resid[i, :] = data[i, :] - gauss_dataset(params, i, x)

    # now flatten this to a 1D array, as minimize() needs
    return resid.flatten()

np.random.seed(2021)
x = np.linspace(-1, 2, 151)
data = []
for _ in np.arange(5):
    amp = 0.60 + 9.50*np.random.rand()
    cen = -0.20 + 1.20*np.random.rand()
    sig = 0.25 + 0.03*np.random.rand()
    dat = gauss(x, amp, cen, sig) + np.random.normal(size=x.size, scale=0.1)
    data.append(dat)
data = np.array(data)
fit_params = Parameters()
for iy, y in enumerate(data):
    fit_params.add(f'amp_{iy+1}', value=0.5, min=0.0, max=200)
    fit_params.add(f'cen_{iy+1}', value=0.4, min=-2.0, max=2.0)
    fit_params.add(f'sig_{iy+1}', value=0.3, min=0.01, max=3.0)


for iy in (2, 3, 4, 5):
    fit_params[f'sig_{iy}'].expr = 'sig_1'

out = minimize(objective, fit_params, method="least_squares", args=(x, data), verbose=2)
report_fit(out.params)


import line_profiler
import cProfile
import pstats
import io
import importlib
import echopop.inversion.pcdwba
importlib.reload(echopop.inversion.pcdwba)
# from echopop.inversion.scattering_models import pcdwba, pcdwba_fbs

def profile_pcdwba(*args, **kwargs):
    pr = cProfile.Profile()
    pr.enable()
    result = pcdwba(*args, **kwargs)
    pr.disable()
    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats("cumulative")
    ps.print_stats(40)  # Show top 40 lines
    print(s.getvalue())
    return result

profile = line_profiler.LineProfiler()
profile.add_function(pcdwba)
profile.enable_by_count()

res = pcdwba(
    center_frequencies=np.array([18e3, 38e3, 70e3, 120e3, 200e3]),
    **scattering_parameters.values,
    **self.model_settings,
    **self.model_settings["environment"]
)
profile.print_stats()


profile = line_profiler.LineProfiler()
profile.add_function(pcdwba_fbs)
profile.enable_by_count()

f_bs = pcdwba_fbs(
    taper_order,
    length_sd_norm,
    length_mean,
    length_radius_ratio,
    radius_of_curvature_ratio,
    theta_radians,
    k_f,
    ka_f,
    g,
    h,
    n_integration,
    n_wavelength,
)
profile.print_stats()