import numpy as np
import warnings
import sys


class SemiVariogram:
    """
    This class contains a routine that calculates a
    standardized semi-variogram and routines for
    obtaining the best semi-variogram model for the
    estimated semi-variogram.

    Parameters
    ----------
    x : Numpy array
        A 1D array representing the x coordinate for
        the semi-variogram calculation.
    y : Numpy array
        A 1D array representing the y coordinate for
        the semi-variogram calculation.
    field : Numpy array
        A 1D array representing the field for the
        semi-variogram calculation (e.g. biomass density).
    """

    def __init__(self, x, y, field):

        self.x = x
        self.y = y
        self.field = field

    def calculate_semi_variogram(self, bins):
        """
        Calculates the semi-variogram standardized by the
        standard deviation of the head multiplied by
        the standard deviation of the tail for each
        lag. This calculation assumes that the mesh
        points (x,y) are isotropic and the search
        area is omnidirectional.

        Parameters
        ----------
        bins

        Returns
        -------

        """

        x_diff = np.subtract.outer(self.x, self.x)
        y_diff = np.subtract.outer(self.y, self.y)

        field_head, field_tail = np.meshgrid(self.field, self.field, indexing='ij')
        field_diff_sqrd = np.power(field_head - field_tail, 2)

        # find the distance between points
        dis = np.sqrt(np.power(x_diff, 2) + np.power(y_diff, 2))

        # obtain the upper triangular portion of dis
        dis_2 = np.triu(dis, k=1)  # TODO: we can put this in a sparse form

        gamma_standardized = []
        for i in range(len(bins) - 1):
            # get indices of distances that are between bins[i] and bins[i+1] (i.e. in the lag)
            ind_in_lag = np.argwhere((bins[i] <= dis_2) & (dis_2 < bins[i + 1]) & (dis_2 != 0))

            # indices in the lag
            x_ind = ind_in_lag[:, 0]
            y_ind = ind_in_lag[:, 1]

            # calculate the semi-variogram value
            gamma = 0.5 * np.mean(field_diff_sqrd[x_ind, y_ind])

            # standardize gamma by the standard deviation of the head
            # multiplied by the standard deviation of the tail
            std_head = np.std(field_head[x_ind, y_ind])
            std_tail = np.std(field_tail[x_ind, y_ind])
            gamma_standardized.append(gamma / (std_head * std_tail))

        return np.array(gamma_standardized)




