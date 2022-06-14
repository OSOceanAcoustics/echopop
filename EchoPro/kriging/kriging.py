import numpy as np
import pandas as pd
import warnings
import sys
import folium
from matplotlib.colors import to_hex
from ..semivariogram import SemiVariogram as SV


class Kriging:
    """
    This class constructs all data necessary
    for kriging and performs kriging

    Parameters
    ----------
    EPro : EchoPro object
        An initialized EchoPro object. Note that any change to
        self.EPro will also change this object.
    """

    def __init__(self, EPro = None):

        self.EPro = EPro

    def run_kriging(self, x_mesh, x_data, y_mesh, y_data):
        kmax = 10
        k_min = 3
        R = 0.0226287
        ratio = 0.001
        nugget = 0.0
        sill = 0.95279
        ls = 0.0075429
        exp_pow = 1.5
        ls_hole_eff = 0.0

        # compute the distance between the mesh points and transect points
        x_diff = np.subtract.outer(x_mesh, x_data)
        y_diff = np.subtract.outer(y_mesh, y_data)
        dis = np.sqrt(x_diff * x_diff + y_diff * y_diff)

        # sort dis up to the kmax smallest elements in each row
        dis_sort_ind = np.argpartition(dis, kmax, axis=1)

        # select only the kmax smallest elements in each row
        dis_kmax_ind = dis_sort_ind[:, :kmax]

        ep = []

        # does Ordinary Kriging, follow Journel and Huijbregts, p. 307
        for row in range(dis_kmax_ind.shape[0]):

            sel_ind = dis_kmax_ind[row, :]
            R_ind = np.argwhere(dis[row, sel_ind] <= R).flatten()

            if len(R_ind) < k_min:
                R_ind = np.argsort(dis[row, sel_ind]).flatten()[:k_min]

            # TODO: put in the M2_unity statement

            # TODO: replace with automatic call to correct model
            M20 = SV.generalized_exp_bessel(dis[row, sel_ind[R_ind]], sill, ls, exp_pow, ls_hole_eff, nugget)

            M2 = np.concatenate([M20, np.array([1.0])])  # for Ordinary Kriging

            #     print(f"M2 = {M2}")

            x1 = x_data[sel_ind[R_ind]]
            y1 = y_data[sel_ind[R_ind]]

            # compute the distance between the points
            x1_diff = np.subtract.outer(x1, x1)
            y1_diff = np.subtract.outer(y1, y1)
            dis1 = np.sqrt(x1_diff * x1_diff + y1_diff * y1_diff)

            # TODO: replace with automatic call to correct model
            K0 = SV.generalized_exp_bessel(dis1, sill, ls, exp_pow, ls_hole_eff, nugget)

            # Add column and row of ones for Ordinary Kriging
            K = np.concatenate([K0, np.ones((len(x1), 1))], axis=1)
            K = np.concatenate([K, np.ones((1, len(x1) + 1))], axis=0)

            # do an inplace fill of diagonal
            np.fill_diagonal(K, 0.0)

            # compute SVD
            u, s, vh = np.linalg.svd(K, full_matrices=True)

            kindx = np.argwhere(np.abs(s / s[0]) > ratio).flatten()

            s_inv = 1.0 / s[kindx]

            k_inv = np.matmul(vh.T[:, kindx], np.diag(s_inv))
            k_inv = np.matmul(k_inv, u[:, kindx].T)

            lamb = np.dot(k_inv, M2)

            #     Vp(i)=sum_nan(lambda(1:nk).*var1)*M2_unity;

            ep.append(np.nansum(lamb * M2))


