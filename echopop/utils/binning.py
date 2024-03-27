"""
Functions used across modules to bin data
"""
from typing import List

import numpy as np

# TODO: should we rename this script to something more general?


def get_bin_ind(input_data: np.ndarray, centered_bins: np.ndarray) -> List[np.ndarray]:
    """
    This function manually computes bin counts given ``input_data``. This
    function is computing the histogram of ``input_data`` using
    bins that are centered, rather than bins that are on the edge.
    The first value is between negative infinity and the first bin
    center plus the bin width divided by two. The last value is
    between the second to last bin center plus the bin width
    divided by two to infinity.


    Parameters
    ----------
    input_data: np.ndarray
        The data to create a histogram of.
    centered_bins: np.ndarray
        An array that specifies the bin centers.

    Returns
    -------
    hist_ind: list
        The index values of input_data corresponding to the histogram

    """

    # get the distance between bin centers
    bin_diff = np.diff(centered_bins) / 2.0

    # fill the first bin
    hist_ind = [np.argwhere(input_data <= centered_bins[0] + bin_diff[0]).flatten()]

    for i in range(len(centered_bins) - 2):
        # get values greater than lower bound
        g_lb = centered_bins[i] + bin_diff[i] < input_data

        # get values less than or equal to the upper bound
        le_ub = input_data <= centered_bins[i + 1] + bin_diff[i + 1]

        # fill bin
        hist_ind.append(np.argwhere(g_lb & le_ub).flatten())

    # fill in the last bin
    hist_ind.append(
        np.argwhere(input_data > centered_bins[-2] + bin_diff[-1]).flatten()
    )

    return hist_ind
