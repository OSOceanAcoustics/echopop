import numpy as np
import pandas as pd
import copy
from scipy.stats import norm
from echopop.survey import Survey
from echopop.spatial.mesh import griddify_lag_distances
import math

survey = Survey( init_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config.yml" ,
                 survey_year_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_config.yml" )

survey.transect_analysis()
survey.stratified_analysis(bootstrap_ci = 0.95, bootstrap_ci_method = "empirical")
survey.kriging_analysis()
survey.stratified_analysis("kriging")
self = survey
analysis_dict = self.analysis
results_dict = self.results
spatial_dict = self.input["spatial"]
settings_dict = self.analysis["settings"]["stratified"]
variogram_parameters = self.analysis["settings"]["kriging"]["variogram_parameters"]
n_lags = 30
# max_range = variogram_parameters["range"]
max_range = 0.06
lag_resolution = variogram_parameters["lag_resolution"]
variable = settings_dict["variable"]
transect_data = self.analysis["kriging"]["transect_df"]
estimates = transect_data["biomass_density"].to_numpy()

# Define the range of possible azimuth angles
azimuth_range = 360

# Compute the range and lags 
# ---- Estimate range based on the lag resolution and number of lags
max_range = lag_resolution * n_lags
# ---- Compute the lags
lags = np.arange(1, n_lags) * lag_resolution

# Calculate the lag distance matrix among transect data
transect_distance_matrix = griddify_lag_distances(transect_data, transect_data)

# Compute angles among coordinates
# ---- Extract the x- and y-coordinates
# -------- x
x_coord = transect_data["x"].to_numpy()
# -------- y
y_coord = transect_data["y"].to_numpy()
# ---- Copmute the differences
# -------- x
x_distance = np.subtract.outer(x_coord, x_coord)
# -------- y
y_distance = np.subtract.outer(y_coord, y_coord)
# ---- Fill the diagonals 
# -------- x
np.fill_diagonal(x_distance, np.nan)
# -------- y
np.fill_diagonal(y_distance, np.nan)
# ---- Calculate angle
angularity = np.arctan(y_distance / x_distance) * 180.0 / np.pi + 180 % 180

# Pre-allocate lag counts
lag_counts = np.zeros(n_lags - 1)

# Pre-allocate summed deviation per lag
lag_deviations = np.zeros_like(lag_counts)

# Pre-allocate summed estimates
lag_estimates = np.zeros_like(lag_counts)

# Pre-allocate summed estimats squared
lag_estimates_squared = np.zeros_like(lag_counts)

# Tally head- and tail-indices
# ---- Head
head_index = np.zeros((len(x_coord), n_lags - 1))
# ---- Tail
tail_index = np.zeros((len(x_coord), n_lags - 1))

# Iterate through the coordinates to calculate the angles
for i in range(len(x_coord)):
    # ---- Extract the iterated distances
    distances = transect_distance_matrix[i][i + 1:]
    # ---- Extract the iterated angles
    angles = angularity[i][i + 1:]
    # ---- Find valid angles
    angles_index = np.where(
        (angles >= -180.0) & (angles < 180.0)
    )[0]
    # ---- Compute the semivariogram
    if len(angles_index) > 0:
        # ---- Filter the values
        # -------- x        
        delta_x = x_coord[i] - x_coord[angles_index]
        # -------- y
        delta_y = y_coord[i] - y_coord[angles_index]
        # -------- estimates
        estimates_filter = estimates[angles_index]
        # ---- Calculate the new distance
        delta_xy = np.sqrt(delta_x ** 2 + delta_y ** 2)
        # ---- Sort the distances from closest to furthest
        distance_index = np.argsort(delta_xy)
        # ---- Apply sorted index to biological variables
        estimates_sorted = estimates_filter[distance_index]
        # ---- Apply sorted index to distances
        distances_sorted = delta_xy[distance_index]
        # ---- Quantize the distances into the equivalent lags
        equivalent_lags = np.round(distances_sorted / lag_resolution).astype(int) + 1
        # ---- Iterate through the k-lags
        # -------- Initialize "k_start"
        k_start = 0
        for k in np.arange(1, n_lags):
            # ---- Find the quantized lags 
            lag_k_index = np.where(equivalent_lags[k_start:len(delta_x)] == k)[0] + k_start
            # ---- Get the length of the indexed lags
            lag_k_count = len(lag_k_index)
            # ---- If values are present, advance
            if lag_k_count > 0:
                # ---- Add to running tallied lag-count
                lag_counts[k - 1] += lag_k_count
                # ---- Advance "k_start"
                k_start += lag_k_count
                # ---- Calculate differences among indexed estimates
                delta_estimates = estimates[i] - estimates_sorted[lag_k_index]
                # ---- Sum estimate for each lag bin
                lag_estimates[k - 1] += estimates_sorted[lag_k_index].sum()
                # ---- Squared sum estimate for each lag bin
                lag_estimates_squared[k - 1] += (estimates_sorted[lag_k_index] ** 2).sum()
                # ---- Sum the squared estimate differences
                lag_deviations[k - 1] += (delta_estimates ** 2).sum()
                # ---- Update head index
                head_index[i, k - 1] = lag_k_count
                # ---- Update lag index
                lag_index = distance_index[lag_k_index] + i
                # ---- Update tail index
                tail_index[lag_index, k - 1] += 1

# Pre-allocate the sigma and mean vectors
sigma_head = np.zeros_like(lag_counts)
mean_head = np.zeros_like(lag_counts)
sigma_tail = np.zeros_like(lag_counts)
mean_tail = np.zeros_like(lag_counts)

# Iterate through the computed values necessary for estimating the mean and variance for each lag
for k in np.arange(1, n_lags):
    # ---- Iterate through the head vector
    head_k_index = np.where(head_index[:, k - 1] >= 1)[0]
    head_k_count = len(head_k_index)
    # ---- Advance if values are present
    if head_k_count > 0:
        # ---- Compute the PDF
        pdf = head_index[head_k_index, k - 1] / lag_counts[k - 1]
        # ---- Calculate the head mean
        mean_head[k - 1] = (estimates[head_k_index] * pdf).sum()
        # ---- Calculate the head sigma
        sigma_head[k - 1] = np.sqrt(
            (((estimates[head_k_index] - mean_head[k - 1]) ** 2) * pdf).sum()
        )
    # ---- Iterate through the tail vector
    tail_k_index = np.where(tail_index[:, k - 1] >= 1)[0]
    tail_k_count = len(tail_k_index)
    # ---- Advance if values are present
    if tail_k_count > 0:
        # ---- Compute the PDF
        pdf = tail_index[tail_k_index, k - 1] / lag_counts[k - 1]
        # ---- Calculate the head mean
        mean_tail[k - 1] = (estimates[tail_k_index] * pdf).sum()
        # ---- Calculate the head sigma
        sigma_tail[k - 1] = np.sqrt(
            (((estimates[tail_k_index] - mean_tail[k - 1]) ** 2) * pdf).sum()
        )
        
# Initialize mean and sigma tail vectors that represent updated estimates
mean_tail_up = np.zeros_like(mean_tail)
sigma_tail_up = np.zeros_like(sigma_tail)

# Initialize the partial sill and semivariance
partial_sill = np.zeros_like(mean_tail_up)
semivariance = np.zeros_like(mean_tail_up)

# Compute the semivariance        
lags_non_zero_counts = np.where(lag_counts >= 1)[0]
if len(lags_non_zero_counts) > 0:
    counts_non_zero = lag_counts[lags_non_zero_counts]
    estimate_mean = lag_estimates[lags_non_zero_counts] / counts_non_zero
    estimate_variance = lag_estimates_squared[lags_non_zero_counts] / counts_non_zero - estimate_mean**2
    mean_tail_up[lags_non_zero_counts] = estimate_mean
    sigma_tail[lags_non_zero_counts] = np.sqrt(np.abs(estimate_variance))
    partial_sill[lags_non_zero_counts] = sigma_tail[lags_non_zero_counts] * sigma_head[lags_non_zero_counts]
    semivariance[lags_non_zero_counts] = 0.5 * lag_deviations[lags_non_zero_counts] / (counts_non_zero * partial_sill[lags_non_zero_counts] + np.finfo(float).eps)
else:
    semivariance = np.full(n_lags - 1, np.nan)
        
# Compute the cross variogram sill
# ---- Find non-zero index
non_zero_index = np.where((sigma_head > 0) & (sigma_tail > 0))[0]
cross_sill = (sigma_head[non_zero_index] * sigma_tail[non_zero_index]).mean()

# Compute the variogram sill
sill = np.std(estimates) ** 2

# Compute empirical boundary estimates for the sill
# ---- Lower
sill_min = 0.7 * semivariance.min()
# ---- Upper
sill_max = 1.2 * semivariance.max()

# Compute empirical boundary for the nugget
# ---- Upper
nugget_max = semivariance.max()

# Estimate the length scale
# ---- Find the maximum value within the first half of the variogram
semivariogram_max = np.argmax(semivariance[:int(np.round(0.5 * len(semivariance)))])
# ---- Find the lag distance corresponding to the 6 dB attenuation from the maximum
lag_6db = np.argmin(np.abs(semivariance[:semivariogram_max] - 0.3))
# ---- Assign the length scale
if lags[lag_6db] <= 0.01 * lags.max():
    # ---- Minimum value
    length_scale = 0.1 * lags.max()
else:
    length_scale = lags[lag_6db]

# Prepend 0's (or other values) where necessary
# ---- Semivariance
semivariance_zero = np.concatenate([[0], semivariance])
# ---- Lags
lags_zero = np.concatenate([[0], lags])
# ---- Lag counts
lag_counts_zero = np.concatenate([[len(estimates) - 1], lag_counts])

# LEAST-SQUARES FITTING

# Calculate the variogram weights (for fitting)
variogram_weights = lag_counts_zero / lag_counts_zero.sum()

semivariance_zero = (
    np.concatenate([[0],
                    pd.read_csv("C:/Users/Brandyn/Documents/gamma_h.csv").values.flatten()])
)

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


# result = least_squares(cost_function, initial_guess, args=(truncated_stack[0], truncated_stack[1], truncated_stack[2]), bounds=(lower_bounds, upper_bounds), **options)
# result = minimize(cost_function, initial_guess, args=(truncated_stack[0], truncated_stack[1], truncated_stack[2]), method='trust-constr', bounds=bounds, options=options)

result = least_squares(cost_function, initial_guess, args=(truncated_stack[0], truncated_stack[1], truncated_stack[2]), bounds=(lower_bounds, upper_bounds), **options)

nugget_opt, sill_opt, lcsl_opt, decay_power_opt, hole_opt = result.x
print(f"Optimized parameters: nugget = {nugget_opt}, sill = {sill_opt}, lcsl = {lcsl_opt}, decay_power = {decay_power_opt}, hole = {hole_opt}")

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

plt.scatter(x,y)
plt.plot(x, new_out, label="Echopop -- constrainted Levenberg-Marquardt: direct fit")
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