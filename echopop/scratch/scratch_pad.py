import pandas as pd
import numpy as np
from typing import Any, Dict, List, Union, Optional, Literal
from pydantic import ValidationError
from pathlib import Path
from typing import Any, Dict
import numpy as np
import os
import numpy.typing as npt
import pandas as pd
from pandarallel import pandarallel
from tqdm import tqdm
from lmfit import Parameters, Minimizer
from echopop.inversion.inversion_base import InversionBase
from echopop.inversion.operations import impute_missing_sigma_bs
from echopop import validators as val
from echopop.nwfsc_feat import ingest_nasc, utils
from echopop import acoustics
from dask.diagnostics import ProgressBar
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from tqdm import tqdm
import dask.dataframe as dd
from echopop.core.echoview import (
    ECHOVIEW_DATABASE_EXPORT_FILESET,
    ECHOVIEW_EXPORT_ROW_SORT,
    ECHOVIEW_TO_ECHOPOP,
)
from echopop.validators.inversion import ValidateInversionMatrix, ValidateBuildModelArgs, SCATTERING_MODEL_PARAMETERS
from echopop.typing import InvParameters, MCInvParameters
from echopop.nwfsc_feat import ingest_sv
from scipy.special import j1
import hashlib
from echopop.typing import InvParameters
from echopop.typing.inversion import ModelInputParameters
from echopop.validators import ValidateInversionMatrix
from echopop.inversion.inversion_matrix_krill import prepare_minimizer, optim, monte_carlo_initialize
import warnings
import time
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

sv_data = ingest_sv.ingest_echoview_sv(sv_path=sv_path, 
                                       center_frequencies=center_frequencies, 
                                       transect_pattern=transect_pattern, 
                                       aggregate_method=method,
                                       impute_coordinates=True)

# !!! EQUIVALENT TO `reduce_dataset` from inversion extension PR
utils.apply_filters(sv_data, include_filter={"transect_num": np.linspace(1, 30, 30)})

####################################################################################################
def prepare_minimizer_safe(*args, **kwargs):
    import numpy as np
    from echopop.typing import MCInvParameters
    from echopop.inversion.inversion_matrix_krill import prepare_minimizer
    try:
        return prepare_minimizer(*args, **kwargs)
    except Exception as e:
        return str(e)  # Return the error message for inspection

class InversionMatrix(InversionBase):
    """
    !!! DOCSTRING
    """
    def __new__(
        cls,      
        data: pd.DataFrame,   
        simulation_settings: Dict[str, Any],
        verbose: bool = True,
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
        self.verbose = verbose
        
        # Generate
        return self

    def __init__(
        self, 
        data: pd.DataFrame,   
        simulation_settings: Dict[str, Any],
        verbose: bool = True,
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
            self.simulation_settings["rng"] = np.random.default_rng(mc_seed)            

    def _set_minimizers(self):
        
        # Definine a new column for the `lmfit.Minimizer` class
        self.measurements["minimizer"] = np.array(np.nan).astype(object)

        # Create list of `Minimizer` objects depending on number of defined realizations
        # # ---- Check for parallelization
        # if self.simulation_settings["parallel"]:
        #     if self.verbose:
        #         "Initializing optimizers using Monte Carlo methods."
        #     ddf = dd.from_pandas(self.measurements, npartitions=min(os.cpu_count() or 4, len(self.measurements)))
        #     with ProgressBar():
        #         results = ddf.map_partitions(lambda df: df.apply(
        #             lambda row: prepare_minimizer(row["sv_mean"], self.model_params, self.model_settings, self.simulation_settings, self.verbose),
        #             axis=1
        #         )).compute()
        #     # Prepare argument tuples for each row
        #     arg_tuples = [
        #         (row["sv_mean"], self.model_params, self.model_settings, self.simulation_settings, self.verbose)
        #         for _, row in self.measurements.iterrows()
        #     ]
        #     results = []

        #     ddf = dd.from_pandas(self.measurements, npartitions=min(os.cpu_count() or 4, len(self.measurements)))
        #     results = ddf.map_partitions(lambda df: df.apply(
        #         lambda row: prepare_minimizer(row["sv_mean"], self.model_params, self.model_settings, self.simulation_settings, self.verbose),
        #         axis=1
        #     )).compute()

        #     def run_prepare_minimizer(args):
        #         return prepare_minimizer(*args)

        #     with ProcessPoolExecutor(max_workers=16) as executor:
        #         for result in tqdm(executor.map(run_prepare_minimizer, arg_tuples), total=len(arg_tuples)):
        #             results.append(result)
        #     # with ThreadPoolExecutor() as executor:
        #     with ThreadPoolExecutor(max_workers=16) as executor:
        #         for result in tqdm(executor.map(lambda args: prepare_minimizer(*args), arg_tuples), 
        #                            total=len(arg_tuples)):
        #             results.append(result)
        #     self.measurements["minimizer"] = results
        # ---- Sequential
        # else:
        self.measurements["minimizer"] = self.measurements["sv_mean"].apply(
            prepare_minimizer,
            axis=1,
            args = (
                self.model_params,
                self.model_settings,
                self.simulation_settings,  
                self.verbose,
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
        model_parameters: InvParameters,
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
        self.model_params = valid_args["model_parameters"]
        self.model_settings.update(valid_args["model_settings"])

        # Update the model settings to include the model function
        self.model_settings["model_function"] = SCATTERING_MODEL_PARAMETERS.get(
            self.model_settings["type"]
        )["function"]

        # Scale the parameters, if defined 
        if self.simulation_settings["scale_parameters"] and not self.model_params.is_scaled:
            self.model_params.scale()

        # Add the `lmfit.Minimizer` objects to the dataset
        self._set_minimizers()
            
    def invert(self, optimization_kwargs):
        
        #
        inversion_results = self.measurements[["sv_mean", "minimizer", "label"]].apply(
            optim,
            axis=1,
            args=(
                self.model_params,
                self.simulation_settings,
                optimization_kwargs,
                self.verbose,
            ),
            result_type="expand",
        )
        return inversion_results


############################################################################################
# # Prepare the optimizer parameterization
# measurements_df = prepare_minimizer(
#     data_df=measurements_df,
#     center_frequencies=center_frequencies,
#     aggregate=aggregate,
#     sv_threshold=sv_threshold,
#     scattering_parameters=self.inversion_config["scattering_parameters"],
#     simulation_parameters=self.inversion_config["processing_parameters"]["simulation"],
#     processing_parameters=self.inversion_config["processing_parameters"]["acoustic_models"],
# )        
Sv_measured = self.measurements["sv_mean"].loc[1]
import time
from echopop.inversion.operations import reflection_coefficient
from lmfit import minimize
import warnings
self = INVERSION

sv = self.measurements["sv_mean"].iloc[0]
res = prepare_minimizer(sv, self.model_params, self.model_settings, self.simulation_settings, self.verbose)
print("direct call ->", type(res), repr(res))

res_apply = self.measurements["sv_mean"].apply(
    lambda s: prepare_minimizer(s, self.model_params, self.model_settings, self.simulation_settings, self.verbose)
)
print("apply with lambda ->", res_apply)

min_list = [
    prepare_minimizer(s, self.model_params, self.model_settings, self.simulation_settings, self.verbose)
    for s in self.measurements["sv_mean"]
]
self.measurements["minimizer"] = pd.Series(min_list, index=self.measurements.index)
print("assigned, first:", type(self.measurements['minimizer'].iloc[0]))

Sv_measured = self.measurements[["sv_mean", "minimizer", "label"]].loc[1]
obj = Sv_measured["minimizer"].iloc[0][0]
print("container type:", type(obj), "is_list:", isinstance(obj, (list, tuple)))
print("repr[:500]:", repr(obj)[:500])
print("len(container):", len(obj))
for i, m in enumerate(obj):
    print("--- item", i, "type:", type(m))
    # lmfit.Minimizer exposes .params
    if hasattr(m, "params"):
        print("  params type:", type(m.params), "n_params:", len(m.params))
        print("  first params keys:", list(m.params.keys())[:8])
        print("  sample param values:", {k: m.params[k].value for k in list(m.params.keys())[:5]})
    else:
        print("  item has no .params, repr:", repr(m)[:300])
scattering_parameters = self.model_params
simulation_settings = self.simulation_settings
OPTIMIZATION_KWARGS = {
    "max_nfev": 200,
    "method": "nelder",
    "loss": "huber",
    "xtol": 1e-8,
    "diff_step": 1e-9,
    "verbose": 1,
    "burnin": {
        "method": "nelder",
        "max_nfev": 200,
        "tol": 1e-3,
    }
}
optimization_kwargs = OPTIMIZATION_KWARGS
optimization_kwargs["diff_step"]: 1e-9
parameters_meta: InvParameters = self.model_params
verbose = True

# Catch start time in case `verbose=True`
start_time = time.time()

# Find which values are below the defined threshold
valid_idx = np.argwhere(Sv_measured["sv_mean"] > -999.).flatten()
center_frequencies = np.array(Sv_measured["sv_mean"].index.values[valid_idx], dtype=float)

# Only run if the correct number of frequencies are valid
if len(valid_idx) >= simulation_settings["minimum_frequency_count"]:

    # Assign message string
    frequency_msg = ""

    # Convert to a numpy array
    Sv = Sv_measured["sv_mean"].to_numpy()[valid_idx]

    # Optimize over all realizations [initialize]
    # ---- Fit errors
    fit_errors = np.ones(simulation_settings["mc_realizations"]) * np.nan
    # ---- Parameter sets (column names/structure)
    parameter_fits = (
        pd.DataFrame(
            scattering_parameters.values,
            index=range(simulation_settings["mc_realizations"]),
        )
        * np.nan
    )

    # Iterate through the realizations
    for realization in range(simulation_settings["mc_realizations"]):
        with warnings.catch_warnings():
            warnings.filterwarnings(action="ignore", category=RuntimeWarning)
            # ---- Optimize
            minimizer = Sv_measured["minimizer"].iloc[0][realization]
            parameters_optimized = minimizer.minimize(
                method="least_squares",
                **optimization_kwargs,
            )        
            params_opt_dict = parameters_optimized.params.valuesdict()
            if parameters_meta.is_scaled:
                params_opt_dict = parameters_meta.unscale_dict(params_opt_dict)
            parameter_fits.loc[realization] = pd.Series(params_opt_dict)
            fit_errors[realization] = parameters_optimized.residual





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
    "mc_realizations": 10,
    "scale_parameters": True, 
    "minimum_frequency_count": 2,
    "reference_frequency": 120e3,
}

OPTIMIZATION_KWARGS = {
    "max_nfev": 500,
    "method": "lbfgsb",
    "gtol": 1e-7,
    "ftol": 1e-7,
    "xtol": 1e-7,
    "diff_step": 1e-6,
    "burnin": {
        "method": "nelder",
        "max_nfev": 200,
        "tol": 1e-3,
    }
}
####################################################################################################
parameters_inv = InvParameters(MODEL_PARAMETERS)
sv_data_sub = utils.apply_filters(sv_data, include_filter={"transect_num": [1, 2, 3, 4]})
MOD = InversionMatrix(sv_data_sub, SIMULATION_SETTINGS)
MOD.build_scattering_model(parameters_inv, MODEL_SETTINGS)
MOD.invert(OPTIMIZATION_KWARGS)


# microbenchmark (run in your REPL)
import time, statistics
def bench(fn, n=5):
    times = []
    for _ in range(n):
        t0 = time.perf_counter()
        fn()
        times.append(time.perf_counter() - t0)
    return statistics.mean(times), statistics.stdev(times)

_model_params = self.model_params
_model_settings = self.model_settings
_simulation_settings = self.simulation_settings
_verbose = self.verbose
_fn = prepare_minimizer

# option A: pandas.apply (what you measured)
def run_applyA():
    self.measurements["sv_mean"].apply(
        _fn,
        axis=1,
        args=(_model_params, _model_settings, _simulation_settings, _verbose),
    )

def run_applyB():
    self.measurements["sv_mean"].apply(
        _fn,
        axis=1,
        args=(self.model_params, self.model_settings, self.simulation_settings, self.verbose),
    )


# option B: iterrows (your earlier listcomp)
def run_iterrows():
    _ = [
        _fn(row["sv_mean"], _model_params, _model_settings, _simulation_settings, _verbose)
        for _, row in self.measurements.iterrows()
    ]

# option C: itertuples (usually fastest for row access)
def run_iloc():
    _ = [
        prepare_minimizer(self.measurements["sv_mean"].iloc[i], _model_params, _model_settings, _simulation_settings, _verbose)
        for i in range(len(self.measurements))
    ]

# option D: iterate numpy column (fastest if only one column needed)
col = self.measurements["sv_mean"].to_numpy().ravel()
def run_numpycol():
    _ = [_fn(sv["sv_mean"], _model_params, _model_settings, _simulation_settings, _verbose) for sv in col]

out = []
for row in self.measurements["sv_mean"].itertuples(index=False):
    print(row)
    out.append(row["sv_mean"])

print("apply:", bench(run_applyA, n=10))
print("apply:", bench(run_applyB, n=10))
print("iterrows:", bench(run_iterrows, n=10))
print("numpycol:", bench(run_iloc, n=10))