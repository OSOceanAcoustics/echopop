import numpy as np
import pandas as pd
import copy
from scipy.stats import norm
from echopop.survey import Survey

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
from echopop.spatial.mesh import griddify_lag_distances
import math

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

# Prepare the data for determining the best-fit/optimized parameters
# ---- 

griddify_lag_distances(transect_data.iloc[i], transect_data.iloc[i + 1:])
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


# ---- Create empty "z" coordinate 
z = np.zeros_like(transect_distance_matrix)
# ---- Calculate  angularity
x_coord = transect_data["x"].to_numpy()
y_coord = transect_data["y"].to_numpy()
x_distance = np.subtract.outer(x_coord, x_coord)
# ---- Differences across y-coordinates
y_distance = np.subtract.outer(y_coord, y_coord)
np.fill_diagonal(x_distance, np.nan)
np.fill_diagonal(y_distance, np.nan)
# ---- Calculate angle
angularity = np.arctan(y_distance / x_distance) * 180.0 / np.pi

# Filter the angles across the entire `angularity` array
# ---- Define possible range of azimuth angles
azimuth_range = 360
lower_bound = angularity - 0.5 * azimuth_range
upper_bound = angularity + 0.5 * azimuth_range
ang_indx = []
ang_indx[1000].size
# Iterate over each row of the angularity array
for row in range(angularity.shape[0]):
    ang_z = angularity[row]
    row_indices = np.where((ang_z >= ang_z - 0.5 * azimuth_range) & (ang_z < ang_z + 0.5 * azimuth_range))[0]
    ang_indx.append(row_indices)

import numpy as np

variable = "biomass_density"
bio = transect_data[variable].values
azm = 0
dip = 0
angz = azm
angd = dip
nlag = n_lags
res = lag_resolution
nd = 9278
x = x_coord
y = y_coord
z = np.zeros_like(x)
vr = bio
azm_deg = azm
dip_deg = dip
azm_res = 360
dip_res = 180
ndir = 1
nl_cnt = int(np.max([1, np.round(0.01 * nd * ndir)]))

# Initialize arrays
cnt = np.zeros(nlag - 1)
gamma_sum = np.zeros(nlag - 1)
var_sum = np.zeros(nlag - 1)
var_mean = np.zeros(nlag - 1)
var2_sum = np.zeros(nlag - 1)
gammah = np.full(nlag - 1, np.nan)
c0 = np.zeros(nlag - 1)
head_indx = np.zeros((nd, nlag - 1))
tail_indx = np.zeros((nd, nlag - 1))

dx = np.subtract.outer(x, x)
dy = np.subtract.outer(y, y)
dz = np.subtract.outer(z, z)
dr = np.sqrt(dx**2 + dy**2)
ang_z = np.degrees(np.arctan2(dy, dx))
ang_d = np.degrees(np.arctan2(dz, dr))

neg_indx = ang_z < azm_deg - 0.5 * azm_res
ang_z[neg_indx] += 180

mask = (ang_z >= azm_deg - 0.5 * azm_res) & (ang_z < azm_deg + 0.5 * azm_res) & \
       (ang_d >= dip_deg - 0.5 * dip_res) & (ang_d < dip_deg + 0.5 * dip_res)

for i in range(nd):
    if (i * nd + 1) % nl_cnt == 0:
        percent = f"{(i * nd + 1) / (nl_cnt * nd) * 100:.2f} percent done"
        print(percent)

    valid_mask = mask[i, i+1:]
    ni = np.sum(valid_mask)
    if ni > 0:
        xi, yi, zi, var_i = x[i+1:][valid_mask], y[i+1:][valid_mask], z[i+1:][valid_mask], vr[i+1:][valid_mask]
        dxi, dyi, dzi = x[i] - xi, y[i] - yi, z[i] - zi
        d = np.sqrt(dxi**2 + dyi**2 + dzi**2)
        indx_sort = np.argsort(d)
        var_sort = var_i[indx_sort]
        dsort = d[indx_sort]
        dindx = np.round(dsort / res).astype(int) + 1
        k_acc = 0
        for k in np.arange(1, n_lags):
            indxk = np.where(dindx[k_acc:ni] == k)[0] + k_acc
            nk = len(indxk)
            if nk > 0:
                cnt[k - 1] += nk
                k_acc += nk
                var_dif = vr[i] - var_sort[indxk]
                var_sum[k - 1] += np.sum(var_sort[indxk])
                var2_sum[k - 1] += np.sum(var_sort[indxk]**2)
                gamma_sum[k - 1] += np.sum(var_dif**2)
                head_indx[i, k - 1] = nk
                jindx = indx_sort[indxk] + i + 1
                tail_indx[jindx, k - 1] += 1

sigma_head = np.zeros(nlag - 1)
mean_head = np.zeros(nlag - 1)
sigma_tail = np.zeros(nlag - 1)
mean_tail = np.zeros(nlag - 1)

for k in np.arange(1, n_lags):
    indx1 = np.where(head_indx[:, k - 1] >= 1)[0]
    if len(indx1) > 0:
        npk = np.nansum(head_indx[:, k - 1])
        pdf = head_indx[indx1, k - 1] / npk
        mean_head[k - 1] = np.nansum(vr[indx1] * pdf)
        sigma_head[k - 1] = np.sqrt(np.nansum(((vr[indx1] - mean_head[k - 1])**2) * pdf))

    indx1 = np.where(tail_indx[:, k - 1] >= 1)[0]
    if len(indx1) > 0:
        npk = np.nansum(tail_indx[:, k - 1])
        pdf = tail_indx[indx1, k - 1] / npk
        mean_tail[k - 1] = np.nansum(vr[indx1] * pdf)
        sigma_tail[k - 1] = np.sqrt(np.nansum(((vr[indx1] - mean_tail[k - 1])**2) * pdf))

indx = np.where(cnt >= 1)[0]
if len(indx) > 0:
    cnt_nz = cnt[indx]
    var_mean[indx] = var_sum[indx] / cnt_nz
    var_std2 = var2_sum[indx] / cnt_nz - var_mean[indx]**2
    mean_tail[indx] = var_mean[indx]
    sigma_tail[indx] = np.sqrt(np.abs(var_std2))
    c0[indx] = sigma_tail[indx] * sigma_head[indx]
    gammah[indx] = 0.5 * gamma_sum[indx] / (cnt_nz * np.mean(c0[indx]) + np.finfo(float).eps)
    gammah[indx] = 0.5 * gamma_sum[indx] / (cnt_nz * c0[indx] + np.finfo(float).eps)
else:
    gammah = np.full(nlag - 1, np.nan)

max_value = 1.1889
indx = np.where(gammah > max(2.5, max_value))[0]
gammah[indx] = np.nan

gam = np.zeros((nlag, 1))
dis = np.zeros((nlag, 1))
np_arr = np.zeros((nlag, 1))
hm = np.zeros((nlag, 1))
tm = np.zeros((nlag, 1))
hv = np.zeros((nlag, 1))
tv = np.zeros((nlag, 1))
xx = np.zeros((nlag, 1))
yy = np.zeros((nlag, 1))
zz = np.zeros((nlag, 1))

indx_ij = np.arange(0, nlag)
gam[indx_ij] = np.concatenate([np.array([0]), gammah])[:, np.newaxis]
dis[indx_ij] = lags[:-1][:, np.newaxis]
np_arr[indx_ij] = np.concatenate([np.array([nd - 1]), cnt])[:, np.newaxis]
hm[indx_ij] = np.concatenate([np.array([0]), mean_head])[:, np.newaxis]
tm[indx_ij] = np.concatenate([np.array([0]), mean_tail])[:, np.newaxis]
hv[indx_ij] = np.concatenate([np.array([0]), sigma_head])[:, np.newaxis]
tv[indx_ij] = np.concatenate([np.array([0]), sigma_tail])[:, np.newaxis]
xx[indx_ij] = (lags[:-1] * np.cos(angz * np.pi / 180) * np.cos(angd * np.pi / 180))[:, np.newaxis]
yy[indx_ij] = (lags[:-1] * np.sin(angz * np.pi / 180) * np.cos(angd * np.pi / 180))[:, np.newaxis]
zz[indx_ij] = (lags[:-1] * np.sin(angd * np.pi / 180))[:, np.newaxis]


variable = "biomass_density"
bio = transect_data[variable].values
azm = 0
dip = 0
angz = azm
angd = dip
nlag = n_lags
res = lag_resolution
nd = 9278
x = x_coord
y = y_coord
z = np.zeros_like(x)
vr = bio
azm_deg = azm
dip_deg = dip
azm_res = 360
dip_res = 180
ndir = 1
nl_cnt = int(np.max([1, np.round(0.01*nd*ndir)]))

# Initialize arrays
cnt = np.zeros(nlag - 1)
gamma_sum = np.zeros(nlag - 1)
var_sum = np.zeros(nlag - 1)
var_mean = np.zeros(nlag - 1)
var2_sum = np.zeros(nlag - 1)
gammah = np.full(nlag - 1, np.nan)
c0 = np.zeros(nlag - 1)
head_indx = np.zeros((nd, nlag - 1))
tail_indx = np.zeros((nd, nlag - 1))
ss = []

for i in range(nd):
    if (i * nd + 1) % nl_cnt == 0:
        percent = f"{(i * nd + 1) / (nl_cnt * nd) * 100:.2f} percent done"
        print(percent)

    dx = x[i + 1:] - x[i]
    dy = y[i + 1:] - y[i]
    dz = z[i + 1:] - z[i]
    dr = np.sqrt(dx**2 + dy**2)
    
    ang_z = np.arctan(dy / dx) * 180 / np.pi
    ang_d = np.arctan(dz / dr) * 180 / np.pi

    neg_indx = np.where(ang_z < azm_deg - 0.5 * azm_res)[0]
    ang_z[neg_indx] += 180

    ang_indx = np.where((ang_z >= azm_deg - 0.5 * azm_res) & (ang_z < azm_deg + 0.5 * azm_res) &
                        (ang_d >= dip_deg - 0.5 * dip_res) & (ang_d < dip_deg + 0.5 * dip_res))[0]

    if len(ang_indx) > 0:
        xi, yi, zi, var_i = x[ang_indx], y[ang_indx], z[ang_indx], vr[ang_indx]
        ni = len(xi)
        dxi, dyi, dzi = x[i] - xi, y[i] - yi, z[i] - zi
        d = np.sqrt(dxi**2 + dyi**2 + dzi**2)
        indx_sort = np.argsort(d)
        var_sort = var_i[indx_sort]
        dsort = d[indx_sort]
        dindx = np.round(dsort / res).astype(int) + 1
        k_acc = 0
        for k in np.arange(1, n_lags):
            indxk = np.where(dindx[k_acc:ni] == k)[0] + k_acc
            nk = len(indxk)
            if nk > 0:
                cnt[k - 1] += nk
                k_acc += nk
                var_dif = vr[i] - var_sort[indxk]
                var_sum[k - 1] += np.sum(var_sort[indxk])
                var2_sum[k - 1] += np.sum(var_sort[indxk]**2)
                gamma_sum[k - 1] += np.sum(var_dif**2)
                head_indx[i, k - 1] = nk
                jindx = indx_sort[indxk] + i
                tail_indx[jindx, k - 1] += 1
    
sigma_head = np.zeros(nlag - 1)
mean_head = np.zeros(nlag - 1)
sigma_tail = np.zeros(nlag - 1)
mean_tail = np.zeros(nlag - 1)

for k in np.arange(1, n_lags):
    indx1 = np.where(head_indx[:, k - 1] >= 1)[0]
    hl = len(indx1)
    if hl > 0:
        npk = np.nansum(head_indx[:, k - 1])
        pdf = head_indx[indx1, k-1] / npk
        mean_head[k - 1] = np.nansum(vr[indx1] * pdf)
        sigma_head[k - 1] = np.sqrt(np.nansum(((vr[indx1] - mean_head[k - 1])**2) * pdf))

    indx1 = np.where(tail_indx[:, k - 1] >= 1)[0]
    hl = len(indx1)
    if hl > 0:
        npk = np.nansum(tail_indx[:, k - 1])
        pdf = tail_indx[indx1, k - 1] / npk
        mean_tail[k - 1] = np.nansum(vr[indx1] * pdf)
        sigma_tail[k - 1] = np.sqrt(np.nansum(((vr[indx1] - mean_tail[k - 1])**2) * pdf))

indx = np.where(cnt >= 1)[0]
if len(indx) > 0:
    cnt_nz = cnt[indx]
    var_mean[indx] = var_sum[indx] / cnt_nz
    var_std2 = var2_sum[indx] / cnt_nz - var_mean[indx]**2
    mean_tail[indx] = var_mean[indx]
    sigma_tail[indx] = np.sqrt(np.abs(var_std2))
    c0[indx] = sigma_tail[indx] * sigma_head[indx]
    gammah[indx] = 0.5 * gamma_sum[indx] / (cnt_nz * np.mean(c0[indx]) + np.finfo(float).eps)
    gammah[indx] = 0.5 * gamma_sum[indx] / (cnt_nz * c0[indx] + np.finfo(float).eps)
else:
    gammah = np.full(nlag - 1, np.nan)

max_value = 1.1889
indx = np.where(gammah > max(2.5, max_value))[0]
gammah[indx] = np.nan

gam = np.zeros((nlag, 1))
dis = np.zeros((nlag, 1))
np_arr = np.zeros((nlag, 1))
hm = np.zeros((nlag, 1))
tm = np.zeros((nlag, 1))
hv = np.zeros((nlag, 1))
tv = np.zeros((nlag, 1))
xx = np.zeros((nlag, 1))
yy = np.zeros((nlag, 1))
zz = np.zeros((nlag, 1))

indx_ij = np.arange(0, nlag)
gam[indx_ij] = np.concatenate([np.array([0]), gammah])[:, np.newaxis]
dis[indx_ij] = lags[:-1][:, np.newaxis]
np_arr[indx_ij] = np.concatenate([np.array([nd - 1]), cnt])[:, np.newaxis]
hm[indx_ij] = np.concatenate([np.array([0]), mean_head])[:, np.newaxis]
tm[indx_ij] = np.concatenate([np.array([0]), mean_tail])[:, np.newaxis]
hv[indx_ij] = np.concatenate([np.array([0]), sigma_head])[:, np.newaxis]
tv[indx_ij] = np.concatenate([np.array([0]), sigma_tail])[:, np.newaxis]
xx[indx_ij] = (lags[:-1] * np.cos(angz * np.pi / 180) * np.cos(angd * np.pi / 180))[:, np.newaxis]
yy[indx_ij] = (lags[:-1] * np.sin(angz * np.pi / 180) * np.cos(angd * np.pi / 180))[:, np.newaxis]
zz[indx_ij] = (lags[:-1] * np.sin(angd * np.pi / 180))[:, np.newaxis]

max_sill = 1.2 * gammah.max()
min_sill = 0.7 * gammah.max()
max_nugt = gammah.max()
min_nugt = gammah.min()
indx1 = np.where(( hv.flatten() > 0) & (tv.flatten() > 0))[0]
indx1 = np.argmin(np.abs(gam[indx1-1] - 0.3))
lcsl = lags[indx1]

cnt = np_arr
wgt = cnt / cnt.sum()
# model_para = [min_nugt, min_sill, lcsl, variogram_parameters["decay_power"], 
#               variogram_parameters["hole_effect_range"]]
model_para = [0.0002, 0.6935, 0.0040, 1.4986, 0.000]
R = 0.06
lags1 = np.concatenate([lags, [0.06]])
indx = np.where(lags1 <= R)[0]

data_in = np.vstack((lags1.flatten(), gam.flatten(), wgt.flatten()))
truncated_data = data_in[:, indx]
optimset = "fminbnd"
MaxIter = 2e3
# LB = [min_nugt*0.99, min_sill*0.99, 0, 1, 0]
LB = np.zeros(5)
UB = [max_nugt*1.01, max_sill*1.01, np.inf, 2, np.inf]
LevenbergMarquardt=1



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

def cost_function(model_para, data_in):
    x = data_in[0, :]
    y = data_in[1, :]
    yr = bessel_exponential(x, model_para[0], model_para[1], model_para[2], model_para[3], model_para[4])
    w = data_in[2, :]
    return (yr - y) * w

from scipy.optimize import curve_fit
from scipy.optimize import least_squares
from scipy.optimize import minimize
from scipy import special
bounds = [(low, np.inf) for low in LB]
# result = least_squares(cost_function, model_para, args=(truncated_data,), bounds = (LB, UB))
opts = {
    "max_nfev": 2000
}
x = truncated_data[0, :]
y = truncated_data[1, :]
result, _ = curve_fit(bessel_exponential, x, y, bounds=(LB, UB), max_nfev=2000)
result, _ = curve_fit(bessel_exponential, x, y, method="lm")
model_para = [0.02, 1.0, 0.07, 1.5, 0.000]
result = least_squares(cost_function, 
                       model_para, 
                       args=(truncated_data,), 
                    #    bounds=(LB, np.inf), 
                       max_nfev=2000, method="lm")
result = minimize(cost_function, model_para, args=(data_in,))
# popt, pcov = curve_fit(bessel_exponential, distance_lags, semivariances, p0=model_para, bounds=(LB, UB))
result
result.x
result
out = bessel_exponential(x, 2.83188e-5, 9.35407160e-01, 7.66394436e-03, 1.67616771e+00, 9.76531531e-08)
out1 = bessel_exponential(x, 2.52638102e-02, 0.94209, 0.0007, 1.4958, 4.4408e-14)
out2 = bessel_exponential(x, 0.00015684, 0.94284, 0.0079618, 1.4986, 2.2204e-14)
C = 0.94209 - 2.52638102e-02
C * ( 1 - np.exp( - (x/0.0007) ** 1.4958)) * special.j0(4.4408e-14 * x) + 2.52638102e-02
bessel_exponential(x, 2.52638102e-02, 0.94209, 0.0007, 1.4958, 4.4408e-14)
bessel_exponential(0.002, 0.00015684, 0.94284, 0.0079618, 1.4986, 2.2204e-14)
C = 0.94284 - 0.00015684
C * ( 1 - np.exp( - (0.058/0.0079618) ** 1.4958)) * special.j0(2.2204e-14 * x) + 0.00015684

bessel_exponential(x, 0.00015684, 0.94284, 0.0079618, 10, 2.2204e-14)

import matplotlib.pyplot as plt
plt.scatter(x, y)
plt.plot(x, out)
plt.plot(x, out2, color='orange')
plt.show()

from lmfit import Model, minimize, create_params, fit_report

vmodel = Model(bessel_exponential)
# params = vmodel.make_params(nugget=dict(value=0.002, min=0.0), sill=dict(value=0.90, min=0.0), 
#                             lcsl=dict(value=0.007, min=0.0), decay_power=dict(value=1.5, min=1.0), 
#                             hole=dict(value=0.0, min=0.0))

pars = create_params(nugget=dict(value=0.0005, 
                                 min=0, 
                                 max=max_nugt), 
                     sill=dict(value=0.9428, 
                               min=min_sill, 
                               max=max_sill), 
                     lcsl=dict(value=0.008, min=0.0, max=np.inf, vary=False), 
                     decay_power=dict(value=1.4958, min=0, max=np.inf), 
                     hole=dict(value=0.000, min=0.0, max=np.inf))

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
    return np.mean(((np.abs((yr - y)) * w)))
    # return (yr - y) * w

# out = minimize(cost_function, pars, kws={"data_in": truncated_data}, method="least_squares", xtol=1e-4, max_nfev=2e3)
# fit_report(out)

result = vmodel.fit(y, pars, distance_lags=x - lag_resolution, method="trf")
fit_report(result)
result1 = minimize(cost_function, pars, kws={"data_in": truncated_data}, method="trf")
dely = result.eval_uncertainty(sigma=3)
om = result1.params.valuesdict()
out = bessel_exponential(x - lag_resolution, om["nugget"], om["sill"], om["lcsl"], om["decay_power"], om["hole"])
fit_report(result1)
out1 = bessel_exponential(x - lag_resolution, 0.00015684, 0.94284, 0.0079618, 1.4986, 2.2204e-14)

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('text', usetex = True)

plt.scatter(x - lag_resolution,y)
# plt.fill_between(x, result.best_fit-dely, result.best_fit+dely, color="#ABABAB",
#                  label=r'3-$\sigma$ uncertainty band')
plt.plot(x - lag_resolution, result.best_fit, label="Echopop -- constrainted Levenberg-Marquardt: direct fit")
plt.plot(x - lag_resolution, out1, color="orange", label="EchoPro -- constrainted Levenberg-Marquardt w/ cost-function")
# plt.plot(x, out, color="red", label = r"Echopop -- constrainted Levenberg-Marquardt w/ cost-function: $\frac{\Sigma|(\gamma_{fit} - \gamma)w|}{n_{lags}}$")
plt.plot(x - lag_resolution, out, color="red", label = "Echopop -- constrainted Levenberg-Marquardt w/ cost-function")
plt.ylabel('Semivariance')
plt.xlabel("Distance lags (nmi)")
plt.legend(loc="lower right")
plt.show()

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