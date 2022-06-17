import numba as nb
import numpy as np
import math


@nb.njit(nb.float64[:, :](nb.float64[:], nb.float64[:]), parallel=True)
def nb_subtract_outer(a, b):
    res = np.empty((a.shape[0], b.shape[0]), dtype=np.float64)
    for i in nb.prange(a.shape[0]):
        for j in range(b.shape[0]):
            res[i, j] = a[i] - b[j]

    return res


@nb.njit(nb.float64[:](nb.float64[:], nb.float64[:]), parallel=True, fastmath=True)
def nb_dis_vec(a, b):
    res = np.empty(a.shape, dtype=np.float64)
    for i in nb.prange(a.shape[0]):
        res[i] = math.sqrt(a[i] * a[i] + b[i] * b[i])
    return res


@nb.njit(nb.float64[:, :](nb.float64[:, :], nb.float64[:, :]), parallel=True, fastmath=True)
def nb_dis_mat(a, b):
    res = np.empty(a.shape)
    for i in nb.prange(a.shape[0]):
        for j in range(a.shape[1]):
            res[i, j] = math.sqrt(a[i, j] * a[i, j] + b[i, j] * b[i, j])
    return res


@nb.njit(nb.float64[:](nb.float64[:], nb.float64[:]), parallel=True, fastmath=True)
def nb_diff_sqrd(a, b):
    res = np.empty(a.shape, dtype=np.float64)
    for i in nb.prange(a.shape[0]):
        res[i] = (a[i] - b[i]) ** 2
    return res


## Functions for CV analysis

@nb.njit
def func(distance, field, num_ind):
    # transect-length weighting factor of the transects in the stratum
    wgt = distance / np.mean(distance)

    # normalized biomass of the transects in the stratum
    rhom_trans_stratum = field / distance

    # transect-length-normalized mean density of the stratum
    rhom = np.nansum(field * distance) / np.nansum(distance)

    # variance of the transect-length weighted biomass within the stratum
    if num_ind != 1:
        var_rhom = np.nansum(wgt ** 2 * (rhom_trans_stratum - rhom) ** 2) / (num_ind * (num_ind - 1))
    else:
        var_rhom = np.nansum(wgt ** 2 * (rhom_trans_stratum - rhom) ** 2) / (num_ind * num_ind)

    #     stratum_field = np.nansum(field)

    return rhom, var_rhom  # , stratum_field


@nb.njit
def func2(JH_fac, num_transects, s_e_ind, distance, field, total_transect_area):
    rhom = np.empty(s_e_ind.shape[0], dtype=np.float64)
    var_rhom = np.empty(s_e_ind.shape[0], dtype=np.float64)

    for i in range(s_e_ind.shape[0]):
        num_ind = round(JH_fac * num_transects[i])

        inds = np.arange(num_transects[i])
        sel_ind = np.random.choice(inds, num_ind, replace=False)

        start = s_e_ind[i][0]
        end = s_e_ind[i][1]

        rhom[i], var_rhom[i] = func(distance[start:end][sel_ind],
                                    field[start:end][sel_ind], num_ind)

    # area weighted variance of the "transect-length weighted biomass"
    CV = np.sqrt(np.nansum(var_rhom * total_transect_area ** 2)) / np.nansum(
        total_transect_area * rhom)

    return CV  # , rhom, var_rhom, stratum_field


@nb.njit(parallel=True)
def func3(nr, JH_fac, num_transects,
          s_e_ind, distance, field, total_transect_area):
    CV_JH_vals = np.empty(nr, dtype=np.float64)

    for i in nb.prange(nr):
        CV = func2(JH_fac, num_transects, s_e_ind, distance, field, total_transect_area)

        CV_JH_vals[i] = CV

    return np.nanmean(CV_JH_vals)