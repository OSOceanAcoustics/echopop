import numba as nb
import numpy as np
import math
from typing import Tuple


@nb.njit(nb.float64[:, :](nb.float64[:], nb.float64[:]), parallel=True)
def nb_subtract_outer(a: nb.float64[:], b: nb.float64[:]) -> nb.float64[:, :]:
    """
    Performs an outer subtraction between the inputs
    ``a`` and ``b``.

    Parameters
    ----------
    a : nb.float64[:]
        The first array
    b : nb.float64[:]
        The second array

    Returns
    -------
    res : nb.float64[:, :]
        An outer subtraction between ``a`` and ``b`` e.g
        ``res[i_0, ... i_M-1, j_0, ... j_M-1] = a[i] - b[j]``
    """
    res = np.empty((a.shape[0], b.shape[0]), dtype=np.float64)
    for i in nb.prange(a.shape[0]):
        for j in range(b.shape[0]):
            res[i, j] = a[i] - b[j]

    return res


@nb.njit(nb.float64[:](nb.float64[:], nb.float64[:]), fastmath=True, parallel=True)
def nb_dis_vec(a: nb.float64[:], b: nb.float64[:]) -> nb.float64[:]:
    """
    Calculates the distance between the elements
    of the input vectors.

    Parameters
    ----------
    a : nb.float64[:]
        The first vector
    b : nb.float64[:]
        The second vector

    Returns
    -------
    res : nb.float64[:]
        A vector representing the distance
        between the elements of ``a`` and
        ``b`` e.g. ``res[i] = sqrt(a[i]^2 + b[i]^2)``
    """
    res = np.empty(a.shape, dtype=np.float64)
    for i in nb.prange(a.shape[0]):
        res[i] = math.sqrt(a[i] * a[i] + b[i] * b[i])
    return res


@nb.njit(nb.float64[:, :](nb.float64[:, :], nb.float64[:, :]), fastmath=True, parallel=True)
def nb_dis_mat(a: nb.float64[:, :], b: nb.float64[:, :]) -> nb.float64[:, :]:
    """
    Calculates the distance between the elements
    of the input matrices.

    Parameters
    ----------
    a : nb.float64[:, :]
        The first matrix
    b : nb.float64[:, :]
        The second matrix

    Returns
    -------
    res : nb.float64[:, :]
        A matrix representing the distance between
        the elements of ``a`` and ``b`` e.g.
        ``res[i,j] = sqrt(a[i,j]^2 + b[i,j]^2)``
    """
    res = np.empty(a.shape)
    for i in nb.prange(a.shape[0]):
        for j in range(a.shape[1]):
            res[i, j] = math.sqrt(a[i, j] * a[i, j] + b[i, j] * b[i, j])
    return res


@nb.njit(nb.float64[:](nb.float64[:], nb.float64[:]), fastmath=True, parallel=True)
def nb_diff_sqrd(a: nb.float64[:], b: nb.float64[:]) -> nb.float64[:]:
    """
    Calculates the squared difference between
    the elements of the inputs.

    Parameters
    ----------
    a : nb.float64[:]
        The first vector
    b : nb.float64[:]
        The second vector

    Returns
    -------
    res : nb.float64[:]
        A vector representing the squared
        difference between ``a`` and ``b``
        e.g. ``res[i] = (a[i] - b[i])^2``
    """
    res = np.empty(a.shape, dtype=np.float64)
    for i in nb.prange(a.shape[0]):
        res[i] = (a[i] - b[i]) ** 2
    return res


##############################################
# The below functions are for CV analysis    #
##############################################

@nb.njit
def compute_mean_var_density(distance: np.ndarray, field: np.ndarray,
                             num_ind: int) -> Tuple[float, float]:
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
    rhom : float
        The transect-length-normalized mean density of the stratum
    var_rhom : float
        The transect-length weighted field within the stratum
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
def seed(val: int):
    """
    Seeds the random number generator.

    Parameters
    ----------
    val : int
        Seed value for the generator
    """
    np.random.seed(val)


@nb.njit
def compute_cv_value(jh_fac: float, num_transects: np.ndarray,
                     s_e_ind: np.ndarray, distance: np.ndarray,
                     field: np.ndarray, total_transect_area: np.ndarray) -> float:
    """
    Computes the CV value for the strata i.e. the area weighted
    variance of the transect-length weighted field.

    Parameters
    ----------
    jh_fac : float
        Portion of points to select within each stratum
    num_transects : np.ndarray
        1D array specifying the number of transects
        within each stratum.
    s_e_ind : np.ndarray
        2D array specifying the indices of the distance and
        field arrays that correspond to each stratum.
    distance : np.ndarray
        1D array of distances between (mean latitude, min longitude)
        and (mean latitude, max longitude).
    field : np.ndarray
        1D array of field values.
    total_transect_area : np.ndarray
        1D array specifying the total area covered by the stratum

    Returns
    -------
    CV : float
        CV value for the provided distance and field values
    """

    rhom = np.empty(s_e_ind.shape[0], dtype=np.float64)
    var_rhom = np.empty(s_e_ind.shape[0], dtype=np.float64)

    for i in range(s_e_ind.shape[0]):

        # for subsets of data, it is possible to get 0 transects in a region.
        # If this is encountered, we set the values to NaN
        if int(num_transects[i]) == int(0):

            rhom[i] = np.nan
            var_rhom[i] = np.nan

        else:

            # randomly select samples within the stratum
            num_ind = round(jh_fac * num_transects[i])
            inds = np.arange(num_transects[i])
            sel_ind = np.random.choice(inds, num_ind, replace=False)

            # start and end indices of the stratum for the distance and field arrays
            start = s_e_ind[i][0]
            end = s_e_ind[i][1]
            rhom[i], var_rhom[i] = compute_mean_var_density(distance[start:end][sel_ind],
                                                            field[start:end][sel_ind], num_ind)

    # area weighted variance of the "transect-length weighted field"
    cv = np.sqrt(np.nansum(var_rhom * total_transect_area ** 2)) / np.nansum(
        total_transect_area * rhom)

    return cv


@nb.njit
def compute_jolly_hampton(nr: int, jh_fac: float, num_transects: np.ndarray,
                          s_e_ind: np.ndarray, distance: np.ndarray,
                          field: np.ndarray, total_transect_area: np.ndarray,
                          seed_val: int):
    """
    Computes the Jolly-Hampton CV value using
    nr iterations.

    Parameters
    ----------
    nr : int
        The number of iterations to run the algorithm for
    jh_fac : float
        Portion of points to select within each stratum
    num_transects : np.ndarray
        1D array specifying the number of transects
        within each stratum.
    s_e_ind : np.ndarray
        2D array specifying the indices of the distance and
        field arrays that correspond to each stratum.
    distance : np.ndarray
        1D array of distances between (mean latitude, min longitude)
        and (mean latitude, max longitude).
    field : np.ndarray
        1D array of field values.
    total_transect_area : np.ndarray
        1D array specifying the total area covered by the stratum
    seed_val : int
        Seed value for the random number generator

    Returns
    -------
    The NaN mean of the nr computed CV values
    """

    cv_jh_vals = np.empty(nr, dtype=np.float64)

    if seed_val is not None:
        seed(seed_val)

    for i in range(nr):
        cv_jh_vals[i] = compute_cv_value(jh_fac, num_transects,
                                         s_e_ind, distance, field,
                                         total_transect_area)

    return np.nanmean(cv_jh_vals)
