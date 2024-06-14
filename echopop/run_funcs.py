import numpy as np
import pandas as pd
import copy
from scipy.stats import norm
from echopop.survey import Survey
from echopop.spatial.mesh import griddify_lag_distances
import math

survey = Survey( init_config_path = "/Users/blucca/github/echopop/config_files/initialization_config.yml" ,
                 survey_year_config_path = "/Users/blucca/github/echopop/config_files/survey_year_2019_config.yml" )
self = survey
survey.transect_analysis()
survey.stratified_analysis(bootstrap_ci = 0.95, bootstrap_ci_method = "empirical")
survey.kriging_analysis()
survey.stratified_analysis("kriging")
self = survey

transect_data = self.analysis["kriging"]["transect_df"]
estimates = transect_data["biomass_density"].to_numpy()
settings_dict = self.analysis["settings"]["kriging"]
variogram_parameters = {
    "max_range": 0.06,
    "lag_resolution": 0.002,
    "n_lags": 30,   
    "azimuth_range": 360.0,  
    "force_lag_zero": True
}

lags, gamma_h, lag_counts, lag_covariance = empirical_variogram(transect_data, variogram_parameters, 
                                                                settings_dict)

variogram_parameters.update({"nugget": 0.0, "decay_power": 1.5})
init_parameters = {"distance_lags": lags, 
                   "sill": {
                       "init": 0.91, 
                       "bounds": (0.0, np.inf),
                       "vary": True
                    },
                    "hole_effect_range": {
                        "init": 0.0,
                        "bounds": (0.0, np.inf),
                        "vary": True
                    }, 
                    "correlation_range": {
                        "init": 0.004,
                        "bounds": (0.0, np.inf),
                        "vary": True
                    }
}
bounds = {
    "lower": {"sill": 0.0, "hole_effect_range": 0.0, }
}

variogram_parameters["model"] = ["bessel", "exponential"]
fit_parameters = ["nugget", "sill", "decay_power", "correlation_range", "hole_effect_range"]
optimization_settings = {
    "max_fun_evaluations": 500,
    "cost_fun_tolerance": 1e-6,
    "solution_tolerance": 1e-4,
    "gradient_tolerance": 1e-4,
    "finite_step_size": 1e-8,
    "trust_region_solver": "exact",
    "x_scale": "jacobian",
    "jacobian_approx": "forward",
}

OPTIMIZATION_OPTIONS = {
    "max_fun_evaluations": "max_nfev",
    "cost_fun_tolerance": "ftol",
    "solution_tolerance": "xtol",
    "gradient_tolerance": "gtol",
    "finite_step_size": "diff_step",
    "trust_region_solver": "exact",
    "x_scale": "x_scale",
    "jacobian_approx": "jac"
}


options = {
    'max_nfev': 2000,                  # Maximum number of function evaluations
    'ftol': 1e-6,                      # Tolerance on cost function
    'xtol': 1e-4,                      # Tolerance on solution
    'gtol': 1e-4,                      # Tolerance on gradient norm
    'verbose': 2,                      # Display convergence messages
    'diff_step': 1e-8,                 # Step size for finite difference approximation
    'tr_solver': 'exact',              # Solver for trust-region subproblems
    'x_scale': 'jac',                  # Jacobian scaling
    'jac': '2-point',                  # Finite difference approximation for Jacobian (equivalent to 'forward')
}

variogram_parameters = {
    "max_range": 0.06,
    "lag_resolution": 0.002,
    "n_lags": 30,   
    "azimuth_range": 360.0,  
    "force_lag_zero": True
}

# ARGS (temp)
fit_variogram: bool = True # - Fit variogram, yay or nay
model: Union[list[str], str] = ["bessel", "exponential"]
fit_parameters: Optional[Union[list[str], str]] = ["nugget", "sill", "decay_power", "correlation_range", "hole_effect_range"]
verbose: bool = True
force_lag_zero: bool = True
azimuth_range: float = 360.0
n_lags: Union[int, np.array] = 30
lag_resolution: Optional[float] = None
hole_effect_range: Optional[float] = None
correlation_range: Optional[float] = None
decay_power: Optional[float] = None
nugget: Optional[float] = None
sill: Optional[float] = None
enhance_semivariance: Optional[bool] = None
max_range: float = 0.06
init_parameters: dict = None
optimization_settings: dict = None
variable = "biomass_density"

self.input["statistics"]["variogram"]["model_config"]
################
self.analysis["settings"].update(
    {
        "variogram": {
            "fit_variogram": fit_variogram,
            "variable": variable,
            "verbose": verbose
        }
    }
)

# Pull input variogram settings
default_values = self.input["statistics"]["variogram"]["model_config"]

variogram_parameters = {
    "azimuth_range": azimuth_range,
    "correlation_range": correlation_range if not None else default_values["correlation_range"],
    "decay_power": decay_power if not None else default_values["decay_power"],
    "fit_variogram": fit_variogram,
    "force_lag_zero": force_lag_zero,
    "hole_effect_range": hole_effect_range if not None else default["hole_effect_range"],
    "lag_reolution": lag_resolution if not None else default_values["lag_resolution"],
    "max_range": max_range,
    "model": model,
    "n_lags": n_lags,
    "nugget": nugget if not None else default_values["nugget"],
    "sill": sill if not None else default_values["sill"],
}
          


optimization_settings
lower_bounds = [0.0, 0.0, 0.0, 0.0 , 0.0]
upper_bounds = [np.inf, np.inf, np.inf, np.inf, np.inf]


# Check against 

fit_parameters in list(function_arguments.keys())


# Validate that the model exists within the available API
# Parse the variogram parameters dictionary for the model name
model_name = variogram_parameters["model"]
# ---- Search for correct model
if len(model_name) > 1:
    args = inspect.signature(VARIOGRAM_MODELS["composite"][tuple(model_name)])
else:
    args = inspect.signature(VARIOGRAM_MODELS["single"][model_name])


inspect.getargs(VARIOGRAM_MODELS["composite"][tuple(model_name)])

# Get argument values/names
arg_names = args.parameters.values()
arg_names.astype(str)
list(args.parameters.keys())



sill_init, nugget_init, length_scale_init = initialize_variogram_parameters(gamma_h, lags)
variogram_parameters
# LEAST-SQUARES FITTING
def fit_variogram(lags: np.ndarray,
                  gamma_h: np.ndarray,
                  lag_counts: np.nddaray,
                  variogram_parameters: dict,
                  parameter_bounds: bool = True,
                  fit_parameters: Optional[list[str]] = None):
    
    # Construct 
    
    # Extract model function name 
    model_id = variogram_parameters["model"]


    #
    if fit_parameters is not None:
        # ---- Check against signature arguments for `model_id`
        if len(model_id) > 1:
            model_function = VARIOGRAM_MODELS["composite"][tuple(model_id)]
        else: 
            model_function = VARIOGRAM_MODELS["single"][model_id]
        
    # Initialize parameters when user-input is absent
    sill_init, nugget_init, length_scale_init = initialize_variogram_parameters(gamma_h, lags)
    tuple(variogram_parameters["model"])
    if len(variogram_parameters["model"]) > 1:
        inspect.signature(VARIOGRAM_MODELS["composite"][tuple(variogram_parameters["model"])])

    fun_name = variogram_parameters["model"]

    for models in VARIOGRAM_MODELS.items():
        if fun_name in models:
            models[fun_name], inspect.signature(models[fun_name])
        

    # ---- Add to variogram parameters, if necessary
    if "initial_parameters" not in variogram_parameters.keys():
        variogram_parameters.update({
            "initial_parameters": {

            }
        })
    variogram_parameters.update({
        ""
    })
    
    # Calculate the variogram weights (for fitting)
    variogram_weights = lag_counts / lag_counts.sum()

# Vertically stack the lags, semivariance, and weights
data_stack = np.vstack((lags_zero, semivariance_zero, variogram_weights))

# Index lag distances that are within the parameterized range
within_range = np.where(lags_zero <= max_range)[0]

# Truncate the data stack
truncated_stack = data_stack[:, within_range]

#### ARGUMENTS
max_iterations = 2e3

def bessel_exponential(distance_lags: np.ndarray, nugget, sill, lcsl, decay_power, hole):
    """
    Calculates the composite J-Bessel and exponential semivariogram model at defined lagged
    distances

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    variogram_parameters: dict
        A dictionary containing required parameters for calculating the semivariogram
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = sill - nugget

    # Calculate the spatial decay term
    decay = 1.0 - np.exp(
        -(
            (distance_lags / lcsl)
            ** decay_power
        )
    )

    # Calculate the hole effect
    hole_effect = special.j0(hole * distance_lags)

    # Compute the composite J-Bessel and exponential semivariogram
    return partial_sill * (decay * hole_effect) + nugget

from lmfit import Model, minimize, create_params, fit_report
from scipy import special



vmodel = Model(bessel_exponential)
pars = create_params(nugget=dict(value=0.0, 
                                 min=0), 
                     sill=dict(value=0.7 * semivariance.max(),
                               min=0.0), 
                     lcsl=dict(value=length_scale, 
                               min=0.0), 
                     decay_power=dict(value=1.5, 
                                      min=0.0), 
                     hole=dict(value=0.000, 
                               min=0.0))

def cost_function(pars, data_in):
    vals = pars.valuesdict()
    nugget = vals["nugget"]
    sill = vals["sill"]
    lcsl = vals["lcsl"]
    decay_power = vals["decay_power"]
    hole = vals["hole"]

    x = data_in[0, :]
    w = data_in[2, :]

    yr = bessel_exponential(x, nugget=nugget, sill=sill, lcsl=lcsl, decay_power=decay_power, hole=hole)
    
    # return np.mean((((yr - y)**2) * w))
    # return np.sum(((np.abs((yr - y)) * w)))
    return (yr - y) * w

def cost_function(x, y, w, nugget, sill, lcsl, decay_power, hole):
    # vals = pars.valuesdict()
    # nugget = vals["nugget"]
    # sill = vals["sill"]
    # lcsl = vals["lcsl"]
    # decay_power = vals["decay_power"]
    # hole = vals["hole"]

    # x = data_in[0, :]
    # w = data_in[2, :]

    yr = bessel_exponential(x, nugget, sill, lcsl, decay_power, hole)
    
    # return np.mean((((yr - y)**2) * w))
    # return np.sum(((np.abs((yr - y)) * w)))
    return (yr - y) * w

from scipy.optimize import least_squares
from scipy.optimize import minimize
def cost_function(params, x, y, w):
    nugget, sill, lcsl, decay_power, hole = params
    yr = bessel_exponential(x, nugget, sill, lcsl, decay_power, hole)
    return (yr - y) * w
options = {
    'maxiter': 2000,              # Equivalent to MaxIter
    # 'max_nfev': 500,
    'disp': True,                 # Display convergence messages
    'xtol': 1e-4,                 # Tolerance on solution
    'finite_diff_rel_step': 1e-8, # Step size for finite difference approximation (forward difference)
    'gtol': 1e-4,                 # Gradient tolerance
    'initial_tr_radius': 0.01,    # Initial trust region radius
}
initial_guess = [0.0, 0.69, 0.004, 1.5, 0.0]
lower_bounds = [0.0, 0.0, 0.0, 0.0, 0.0]
upper_bounds = [np.inf, np.inf, np.inf, np.inf, np.inf]
bounds = list(zip(lower_bounds, upper_bounds))
options = {
    'max_nfev': 500,                   # Equivalent to MaxFunEvals
    # 'max_iter': 2000,                  # Equivalent to MaxIter
    'ftol': 1e-4,                      # Equivalent to TolFun, not specified
    'xtol': 1e-4,                      # Equivalent to TolX
    'gtol': None,
    'verbose': 2,                      # No direct equivalent to Display 'notify'; set to 0 for no output
    'diff_step': 1e-8,                 # Not specified in MATLAB settings; default value used
    'tr_solver': 'exact',              # Not specified in MATLAB settings; default value used
    'jac':"2-point",                     # Equivalent to Jacobian, not specified
    'x_scale': 'jac',
    'method': 'trf',                   # Not specified in MATLAB settings; default value used
}
options = {
    'max_nfev': 2000,                  # Maximum number of function evaluations
    'ftol': 1e-6,                      # Tolerance on cost function
    'xtol': 1e-4,                      # Tolerance on solution
    'gtol': 1e-4,                      # Tolerance on gradient norm
    'verbose': 2,                      # Display convergence messages
    'diff_step': 1e-8,                 # Step size for finite difference approximation
    'tr_solver': 'exact',              # Solver for trust-region subproblems
    'x_scale': 'jac',                  # Jacobian scaling
    'jac': '2-point',                  # Finite difference approximation for Jacobian (equivalent to 'forward')
}


# result = least_squares(cost_function, initial_guess, args=(truncated_stack[0], truncated_stack[1], truncated_stack[2]), bounds=(lower_bounds, upper_bounds), **options)
# result = minimize(cost_function, initial_guess, args=(truncated_stack[0], truncated_stack[1], truncated_stack[2]), method='trust-constr', bounds=bounds, options=options)

result = least_squares(cost_function, initial_guess, args=(truncated_stack[0], truncated_stack[1], truncated_stack[2]), bounds=(lower_bounds, upper_bounds), **options)

nugget_opt, sill_opt, lcsl_opt, decay_power_opt, hole_opt = result.x
print(f"Optimized parameters: nugget = {nugget_opt}, sill = {sill_opt}, lcsl = {lcsl_opt}, decay_power = {decay_power_opt}, hole = {hole_opt}")

new_out = bessel_exponential(truncated_stack[0], nugget_opt, sill_opt, lcsl_opt, decay_power_opt, hole_opt)
x = truncated_stack[0]
y = truncated_stack[1]

import matplotlib.pyplot as plt
plt.scatter(x,y)
plt.plot(x, new_out, label="Echopop -- constrainted Levenberg-Marquardt: direct fit")
plt.plot(x, out1, color="orange", label=r"EchoPro -- constrainted Levenberg-Marquardt w/ cost-function: $(\gamma_{fit} - \gamma)w$")
plt.ylabel(r'Semivariance [$\gamma$]')
plt.xlabel("Distance lags (nmi)")
plt.legend(loc="lower right")
plt.show()


from lmfit import Minimizer, Parameters, report_fit

# Initial guess for the parameters
params = Parameters()
params.add('nugget', value=0.0, min=0.0, max=nugget_max)
params.add('sill', value=0.69, min=sill_min, max=sill_max)
params.add('lcsl', value=0.004, min=0.0)
params.add('decay_power', value=1.5, min=1.2, max=4.0)
params.add('hole', value=0, min=0.0, max=8.0)

# Configure optimization options
max_nfev = 500  # Maximum number of function evaluations
max_iter = 2000  # Maximum number of iterations
ftol = 1e-2      # Function tolerance
xtol = 1e-2      # Step tolerance
gtol = 1e-4    # Gradient tolerance
def residual(params, x, y, w):
    nugget = params['nugget']
    sill = params['sill']
    lcsl = params['lcsl']
    decay_power = params['decay_power']
    hole = params['hole']
    
    yr = bessel_exponential(x, nugget, sill, lcsl, decay_power, hole)
    return (yr - y) * w
# Perform the least squares optimization
minimizer = Minimizer(residual, params, fcn_args=(truncated_stack[0], truncated_stack[1], truncated_stack[2]))
result = minimizer.minimize(method='tcf', options={'max_nfev': 500, 'ftol': ftol, 'xtol': xtol, 'gtol': gtol})
report_fit(result)

new_out = bessel_exponential(truncated_stack[0], nugget_opt, sill_opt, lcsl_opt, decay_power_opt, hole_opt)
x = truncated_stack[0]
y = truncated_stack[1]

import matplotlib.pyplot as plt
plt.scatter(x,y)
plt.plot(x, new_out, label=r"Echopop -- constrainted Levenberg-Marquardt w/ cost-function: $(\gamma_{fit} - \gamma)w$")
plt.plot(x, out1, color="orange", label=r"EchoPro -- constrainted Levenberg-Marquardt w/ cost-function: $(\gamma_{fit} - \gamma)w$")
plt.ylabel(r'Semivariance [$\gamma$]')
plt.xlabel("Distance lags (nmi)")
plt.legend(loc="lower right")
plt.show()



vmodel = Model(bessel_exponential)
params = vmodel.make_params(nugget=dict(value=0.0, 
                                 min=0), 
                            sill=dict(value=0.7 * semivariance.max(),
                                    min=0.0), 
                            lcsl=dict(value=length_scale, 
                                    min=0.0), 
                            decay_power=dict(value=1.5, 
                                            min=0.0), 
                            hole=dict(value=0.000, 
                                    min=0.0))
fit = vmodel.fit(truncated_stack[1], params, distance_lags = truncated_stack[0])
fit.fit_report()
pars = create_params(nugget=dict(value=0.0, 
                                 min=0), 
                     sill=dict(value=0.7 * semivariance.max(),
                               min=0.0), 
                     lcsl=dict(value=length_scale, 
                               min=0.0), 
                     decay_power=dict(value=1.5, 
                                      min=0.0), 
                     hole=dict(value=0.000, 
                               min=0.0))


def cost_function(x, y, w, nugget, sill, lcsl, decay_power, hole):
    vals = pars.valuesdict()
    nugget = vals["nugget"]
    sill = vals["sill"]
    lcsl = vals["lcsl"]
    decay_power = vals["decay_power"]
    hole = vals["hole"]

    x = data_in[0, :]
    w = data_in[2, :]

    yr = bessel_exponential(x, nugget=nugget, sill=sill, lcsl=lcsl, decay_power=decay_power, hole=hole)
    
    # return np.mean((((yr - y)**2) * w))
    # return np.sum(((np.abs((yr - y)) * w)))
    return (yr - y) * w

def cost_function(x, w, nugget, sill, lcsl, decay_power, hole):
    yr = bessel_exponential(x, nugget=nugget, sill=sill, lcsl=lcsl, decay_power=decay_power, hole=hole)
    return yr * w

# out = minimize(cost_function, pars, kws={"data_in": truncated_data}, method="least_squares", xtol=1e-4, max_nfev=2e3)
# fit_report(out)
import numdifftools as nd
from scipy.optimize import curve_fit
x = truncated_stack[0]
y = truncated_stack[1]
result = vmodel.fit(y, pars, distance_lags=x, method="lm", max_nfev=2000, xtol=1e-6, ftol=1e-6)
fit_report(result)
fprime = lambda x: scipy.optimize.approx_fprime(x, f, 0.01)
Hfun = nd.Hessdiag(cost_function)
result1 = minimize(cost_function, pars, kws={"data_in": truncated_stack}, method="trust-ncg", jac=fprime, hess=Hfun, tol=1e-4, options={"damp": 0.1})
dely = result.eval_uncertainty(sigma=3)
om = result1.params.valuesdict()
out = bessel_exponential(x, om["nugget"], om["sill"], om["lcsl"], om["decay_power"], om["hole"])
fit_report(result1)
out1 = bessel_exponential(x, 0.00015684, 0.94284, 0.0079618, 1.4986, 2.2204e-14)
params, _ = curve_fit(bessel_exponential, x, y, p0=[0.0, 0.69, 0.004, 1.5, 0.00], bounds=(0, np.inf))
bessel_exponential(x, params)
import matplotlib.pyplot as plt
# matplotlib.rc('text', usetex = True)

plt.scatter(x,y)
# plt.fill_between(x, result.best_fit-dely, result.best_fit+dely, color="#ABABAB",
#                  label=r'3-$\sigma$ uncertainty band')
plt.plot(x, result.best_fit, label="Echopop -- constrainted Levenberg-Marquardt: direct fit")
plt.plot(x, out1, color="orange", label=r"EchoPro -- constrainted Levenberg-Marquardt w/ cost-function: $(\gamma_{fit} - \gamma)w$")
# plt.plot(x, out, color="red", label = r"Echopop -- constrainted Levenberg-Marquardt w/ cost-function: $\frac{\Sigma|(\gamma_{fit} - \gamma)w|}{n_{lags}}$")
plt.plot(x, out, color="red", label = r"Echopop -- constrainted Levenberg-Marquardt w/ cost-function: $(\gamma_{fit} - \gamma)w$")
plt.ylabel(r'Semivariance [$\gamma$]')
plt.xlabel("Distance lags (nmi)")
plt.legend(loc="lower right")
plt.show()

np.mean(np.abs(result.best_fit - out1))
np.mean(np.abs(out - out1))

plt.plot(x, result.best_fit - out)
plt.show()

fit_report(result)
plt.plot(result.eval_components(x=x))

out = minimize(cost_function, )

np.abs(np.real(result[0]))
np.abs(np.real(result[1]))
np.abs(np.real(result[2]))
np.abs(np.real(result[3]))
np.abs(np.real(result[4]))

indx = np.where(cnt >= 1)[0]
cnt_nz = cnt[indx]

mask = (angularity >= lower_bound[:, np.newaxis]) & (angularity < upper_bound[:, np.newaxis])

# Use np.where to find the indices where the mask is True for each row
ang_indx = np.array([np.where(row)[0] for row in mask])

ang_indx = np.where((angularity >= angularity - 0.5 * azimuth_range) & (angularity < angularity + 0.5 * azimuth_range))

bio = transect_data[variable].values
bio_reshaped = bio[:, np.newaxis]
field_matrix = bio_reshaped - bio
np.fill_diagonal(field_matrix, np.nan)



np.subtract.outer(transect_data[variable], transect_data[variable])
from scipy.spatial.distance import pdist, squareform

coords = transect_data[["x", "y"]].values
# values = transect_data[variable].values
values = bio
dists = squareform(pdist(coords, "euclidean"))
diffs = squareform(pdist(values[:, np.newaxis], "euclidean")) ** 2

bin_indices = np.digitize(dists, lags)
bin_sums = np.bincount(bin_indices.ravel(), weights=diffs.ravel(), minlength=30)
bin_counts = np.bincount(bin_indices.ravel(), minlength=len(lags))

with np.errstate(invalid="ignore"):
    semivariances = bin_sums / bin_counts
semivariances = semivariances[:len(lags) - 1]
h = (lags[:-1] + lag_resolution / 2)


def exponential_model(h, nugget, sill, range_):
    return nugget + sill * (1 - np.exp(-h / range_))

from scipy.optimize import curve_fit
initial = [0, np.var(semivariances), np.max(h) / 2]
params, _ = curve_fit(exponential_model, h, semivariances, p0=initial, bounds=(0, np.inf))

def bessel_exponential(distance_lags: np.ndarray, variogram_parameters: dict):
    """
    Calculates the composite J-Bessel and exponential semivariogram model at defined lagged
    distances

    Parameters
    ----------
    distance_lags: np.ndarray
        An array of lag distances
    variogram_parameters: dict
        A dictionary containing required parameters for calculating the semivariogram
    """

    # Calculate the partial sill (or the sill minus the nugget)
    partial_sill = variogram_parameters["sill"] - variogram_parameters["nugget"]

    # Calculate the spatial decay term
    decay = 1.0 - np.exp(
        -(
            (distance_lags / variogram_parameters["correlation_range"])
            ** 1.5
        )
    )

    # Calculate the hole effect
    hole_effect = special.j0(variogram_parameters["hole_effect_range"] * distance_lags)

    # Compute the composite J-Bessel and exponential semivariogram
    return partial_sill * (decay * hole_effect) + variogram_parameters["nugget"]


def params_list_to_dict(params_list):
    return {
        "nugget": params_list[0],
        "sill": params_list[1],
        "correlation_range": params_list[2],
        "decay_power": params_list[3],
        "hole_effect_range": params_list[4],
    }

def params_dict_to_list(params_dict):
    return [
        params_dict["nugget"],
        params_dict["sill"],
        params_dict["decay_power"], 
        params_dict["correlation_range"],
        params_dict["hole_effect_range"],
    ]

initial_params = list(map(variogram_parameters.get, ["nugget", "sill", "correlation_range", "decay_power", "hole_effect_range"]))
from scipy import special

def model_function(h, nugget, sill, correlation_range, decay_power, hole_effect_range):
    params_dict = {
        "nugget": nugget,
        "sill": sill,
        "correlation_range": correlation_range,
        "decay_power": decay_power,
        "hole_effect_range": hole_effect_range,
    }
    return bessel_exponential(h, params_dict)

params, _ = curve_fit(model_function, h, semivariances, p0=initial_params)


variogram_parameters["nugget", "sill"]

variogram_parameters.get(["nugget", "sill", "correlation_range", "hole_effect_range"])
initial_guess = params_dict_to_list(variogram_parameters)

import matplotlib.pyplot as plt

plt.plot(h, semivariances)
plt.show()

semivariances = np.zeros(len(lags)-1)

for i in range(len(lags) -1):
    mask = (dists >= lags[i]) & (dists < lags[i+1])
    if np.any(mask):
        semivariances[i] = np.mean(diffs[mask]) / 2

from echopop.analysis import process_transect_data, stratified_summary
from echopop.spatial.transect import transect_array
from echopop.biology import (
    age1_metric_proportions,
    distribute_length_age,
    filter_species,
    fit_length_weight_relationship,
    fit_length_weights,
    number_proportions,
    partition_transect_age,
    quantize_number_counts,
    quantize_weights,
    weight_proportions,
)
from echopop.spatial.krige import kriging
from echopop.spatial.mesh import crop_mesh, mesh_to_transects, stratify_mesh
from echopop.spatial.projection import transform_geometry
from echopop.spatial.transect import (
    correct_transect_intervals,
    edit_transect_columns,
    save_transect_coordinates,
    summarize_transect_strata,
    transect_distance,
)
from echopop.statistics import stratified_transect_statistic, bootstrap_confidence_intervals

dataframe_list = [
            input_dict["biology"]["length_df"],
            input_dict["biology"]["specimen_df"],
            input_dict["biology"]["catch_df"],
        ]
species_id = [ 22500 , 23010 , 23202 ]
input_dict = self.input
analysis_dict = self.analysis["transect"]
settings_dict = self.analysis["settings"]
configuration_dict = self.config
spatial_dict = input_dict["spatial"]
dataframe_list[0]["species_id"].isin(species_id)

analysis_dict = self.analysis
results_dict = self.results
spatial_dict = self.input["spatial"]
settings_dict = self.analysis["settings"]["stratified"]

import matplotlib.pyplot as plt
dat = self.analysis["stratified"]["transect"]["stratified_replicates_df"]
import scipy as sp

sp.stats.shapiro(dat[ 'survey_cv' ])

plt.hist( dat[ 'survey_cv' ] , bins = 40)
plt.show()

analysis_dict = self.analysis
kriged_mesh = self.results["kriging"]["mesh_results_df"]
settings_dict = self.analysis["settings"]["kriging"]

aged_age_length_table = aged_pivot
unaged_length_table = unaged_pivot
aged_length_totals = aged_length_totals
unaged_apportioned_table = unaged_apportioned_values


def impute_kriged_values(aged_age_length_table: pd.DataFrame,
                         unaged_length_table: pd.DataFrame,
                         aged_length_totals: pd.DataFrame,
                         unaged_apportioned_table: pd.DataFrame):
    
    # Extract the biological variable name (independent of area)
    biology_col = settings_dict["variable"].replace("_density", "")

    # Imputation is required when unaged values are present but aged values are absent at shared
    # length bins! This requires an augmented implementation to address this accordingly
    # ---- Sum across all age bins (of the aged fish data) to generate totals for each row (i.e.
    # ---- length bin)
    summed_aged_length_totals = aged_age_length_table.T.sum()
    # ---- Extract the indices of the summed totals that equal 0.0 for male and female fish
    # -------- Male
    male_zero_aged = np.where(summed_aged_length_totals.loc["male"] == 0.0)[0]
    # -------- Female
    female_zero_aged = np.where(summed_aged_length_totals.loc["female"] == 0.0)[0]
    # ---- Extract the inverse where biological totals are present
    # -------- Male
    male_nonzero_aged = np.where(summed_aged_length_totals.loc["male"] != 0.0)[0]
    # -------- Female
    female_nonzero_aged = np.where(summed_aged_length_totals.loc["female"] != 0.0)[0]
    # ---- Pivot the unaged data and find male and female values that are non-zero
    # -------- Male
    male_nonzero_unaged = unaged_length_table["male"].iloc[male_zero_aged] != 0.0
    # -------- Convert to index
    male_nonzero_unaged_idx = male_zero_aged[male_nonzero_unaged]
    # -------- Female
    female_nonzero_unaged = unaged_length_table["female"].iloc[female_zero_aged] != 0.0
    # -------- Convert to index
    female_nonzero_unaged_idx = female_zero_aged[female_nonzero_unaged]
    # ---- Re-pivot the unaged apportioned values (if necessary)
    if (len(male_nonzero_unaged) > 0) | (len(female_nonzero_unaged)) > 0:
        unaged_values_pvt = (
            unaged_apportioned_table.copy()
            .unstack()
            .reset_index(name="values")
            .pivot_table(
                index=["length_bin"], columns=["sex", "age_bin"], values="values", observed=False
            )
        )
        # ---- Find the closest indices that can be used for nearest-neighbors imputation
        if len(male_nonzero_unaged) > 0:
            # -------- Male
            imputed_male = male_nonzero_aged[
                np.argmin(
                    np.abs(male_zero_aged[male_nonzero_unaged][:, np.newaxis] - male_nonzero_aged),
                    axis=1,
                )
            ]
            # ---- Update the values
            unaged_values_pvt.iloc[
                male_nonzero_unaged_idx, unaged_values_pvt.columns.get_loc("male")
            ] = (
                unaged_length_table["male"].iloc[male_nonzero_unaged_idx].to_numpy()
                * aged_age_length_table.loc["male"].iloc[imputed_male].T
                / aged_length_totals["male"].iloc[imputed_male]
            ).T
        if len(female_nonzero_unaged) > 0:
            # -------- Female
            imputed_female = female_nonzero_aged[
                np.argmin(
                    np.abs(
                        female_zero_aged[female_nonzero_unaged][:, np.newaxis] - female_nonzero_aged
                    ),
                    axis=1,
                )
            ]
            # ---- Update the values
            unaged_values_pvt.iloc[
                female_nonzero_unaged_idx, unaged_values_pvt.columns.get_loc("female")
            ] = (
                unaged_length_table["female"].iloc[female_nonzero_unaged_idx].to_numpy()
                * aged_age_length_table.loc["female"].iloc[imputed_female].T
                / aged_length_totals["female"].iloc[imputed_female]
            ).T
        # ---- Update the original unaged apportioned table
        unaged_apportioned_table = (
            unaged_values_pvt.unstack()
            .reset_index(name="values")
            .pivot_table(
                index=["length_bin"], columns=["age_bin", "sex"], values="values", observed=False
            )
        )
        # ---- Alert message (if verbose = T)
        if settings_dict["verbose"]:
            # ---- Male:
            if len(male_nonzero_unaged) > 0:
                # ---- Get interval values
                intervals_list = [
                    str(interval)
                    for interval in male_nonzero_unaged.index[male_nonzero_unaged].values
                ]
                # ---- Print
                print(
                    f"""Imputed apportioned unaged male {biology_col} at length bins:\n"""
                    f"""{', '.join(intervals_list)}"""
                )
            # ---- Female:
            if len(female_nonzero_unaged) > 0:
                # ---- Get interval values
                intervals_list = [
                    str(interval)
                    for interval in female_nonzero_unaged.index[female_nonzero_unaged].values
                ]
                # ---- Print
                print(
                    f"""Imputed apportioned unaged female {biology_col} at length bins:\n"""
                    f"""{', '.join(intervals_list)}"""
                )
    # ---- Sum the aged and unaged estimates together
    return (
        (unaged_apportioned_table + aged_age_length_table.unstack("sex"))
        .unstack()
        .reset_index(name=f"{biology_col}_apportioned")
    )