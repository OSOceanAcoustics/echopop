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


##############################################
# The below functions are for CV analysis    #
##############################################

@nb.njit
def compute_mean_var_density(distance, field, num_ind):
    """
    Computes the transect-length-normalized mean density
    of the stratum and its associated variance for the
    provided field value.

    Parameters
    ----------
    distance : Numpy array
        1D array of distances between (mean latitude, min longitude)
        and (mean latitude, max longitude).
    field : Numpy array
        1D array of field values.
    num_ind : int
        The number of samples selected within the transect

    Returns
    -------
    rhom : Numpy array
        1D numpy array of the transect-length-normalized mean
        density of the stratum
    var_rhom : Numpy array
        1D numpy array of the variance of the transect-length
        weighted field within the stratum

    """

    # transect-length weighting factor of the transects in the stratum
    wgt = distance / np.mean(distance)

    # normalized field of the transects in the stratum
    rhom_trans_stratum = field / distance

    # transect-length-normalized mean density of the stratum
    rhom = np.nansum(field * distance) / np.nansum(distance)

    # variance of the transect-length weighted field within the stratum
    if num_ind != 1:
        var_rhom = np.nansum(wgt ** 2 * (rhom_trans_stratum - rhom) ** 2) / (num_ind * (num_ind - 1))
    else:
        var_rhom = np.nansum(wgt ** 2 * (rhom_trans_stratum - rhom) ** 2) / (num_ind * num_ind)

    return rhom, var_rhom


@nb.njit
def seed(val):
    """ seeds the random number generator """
    np.random.seed(val)


@nb.njit
def compute_cv_value(JH_fac, num_transects, s_e_ind, distance, field, total_transect_area):
    """
    Computes the CV value for the strata i.e. the area weighted variance of the
    transect-length weighted field.

    Parameters
    ----------
    JH_fac : float
        Portion of points to select within each stratum
    num_transects : Numpy array
        1D array specifying the number of transects
        within each stratum.
    s_e_ind : Numpy array
        2D array specifying the indices of the distance and
        field arrays that correspond to each stratum.
    distance : Numpy array
        1D array of distances between (mean latitude, min longitude)
        and (mean latitude, max longitude).
    field : Numpy array
        1D array of field values.
    total_transect_area : Numpy array
        1D array specifying the total area covered by the stratum

    Returns
    -------
    CV : float
        CV value for the provided distance and field values

    """

    rhom = np.empty(s_e_ind.shape[0], dtype=np.float64)
    var_rhom = np.empty(s_e_ind.shape[0], dtype=np.float64)

    for i in range(s_e_ind.shape[0]):

        # randomly select samples within the stratum
        num_ind = round(JH_fac * num_transects[i])
        inds = np.arange(num_transects[i])
        sel_ind = np.random.choice(inds, num_ind, replace=False)

        # start and end indices of the stratum for the distance and field arrays
        start = s_e_ind[i][0]
        end = s_e_ind[i][1]
        rhom[i], var_rhom[i] = compute_mean_var_density(distance[start:end][sel_ind],
                                                        field[start:end][sel_ind], num_ind)

    # area weighted variance of the "transect-length weighted field"
    CV = np.sqrt(np.nansum(var_rhom * total_transect_area ** 2)) / np.nansum(
        total_transect_area * rhom)

    return CV


@nb.njit
def compute_jolly_hampton(nr, JH_fac, num_transects,
          s_e_ind, distance, field, total_transect_area, seed_val):
    """
    Computes the Jolly-Hampton CV value using
    nr iterations.

    Parameters
    ----------
    nr : int
        The number of iterations to run the algorithm for
    JH_fac : float
        Portion of points to select within each stratum
    num_transects : Numpy array
        1D array specifying the number of transects
        within each stratum.
    s_e_ind : Numpy array
        2D array specifying the indices of the distance and
        field arrays that correspond to each stratum.
    distance : Numpy array
        1D array of distances between (mean latitude, min longitude)
        and (mean latitude, max longitude).
    field : Numpy array
        1D array of field values.
    total_transect_area : Numpy array
        1D array specifying the total area covered by the stratum
    seed_val : int
        Seed value for the random number generator

    Returns
    -------
    The NaN mean of the nr computed CV values

    """

    CV_JH_vals = np.empty(nr, dtype=np.float64)

    if seed_val is not None:
        seed(seed_val)

    for i in range(nr):
        CV_JH_vals[i] = compute_cv_value(JH_fac, num_transects,
                                         s_e_ind, distance, field,
                                         total_transect_area)

    return np.nanmean(CV_JH_vals)
