import numpy as np
import folium
import branca.colormap as cm
from ..numba_modules import nb_subtract_outer, nb_dis_mat


class Kriging:
    """
    This class constructs all data necessary
    for kriging and performs kriging

    Parameters
    ----------
    survey : Survey
        An initialized Survey object. Note that any change to
        self.survey will also change this object.
    """

    def __init__(self, survey=None, params={}):

        self.survey = survey

        # Kriging parameters
        self.k_max = params['k_max']
        self.k_min = params['k_min']
        self.R = params['R']
        self.ratio = params['ratio']

        # # parameters for semi-variogram model
        self.s_v_params = params['s_v_params']

        # grab appropriate semi-variogram model
        self.s_v_model = params['s_v_model']

    def __compute_k_smallest_distances(self, x_mesh, x_data,
                                       y_mesh, y_data, k_max):
        """
        Computes the distance between the data
        and the mesh points. Using the distance it then
        selects the kmax smallest values amongst the
        columns.

        Parameters
        ----------
        x_mesh : Numpy array
            1D array denoting the x coordinates of the mesh
        x_data : Numpy array
            1D array denoting the x coordinates of the data
        y_mesh : Numpy array
            1D array denoting the y coordinates of the mesh
        y_data : Numpy array
            1D array denoting the y coordinates of the data
        k_max : int
            Denotes the maximum number of data points within the
            search radius.

        Returns
        -------
        dis : Numpy array
            2D numpy array representing the distance between each
            mesh points and each transect point
        dis_kmax_ind : Numpy array
            A 2D Numpy array index array representing the kmax
            closest transect points to each mesh point.
        """

        # compute the distance between the mesh points and transect points
        x_diff = nb_subtract_outer(x_mesh, x_data)
        y_diff = nb_subtract_outer(y_mesh, y_data)
        dis = nb_dis_mat(x_diff, y_diff)

        # sort dis up to the kmax smallest elements in each row
        dis_sort_ind = np.argpartition(dis, k_max, axis=1)

        # select only the kmax smallest elements in each row
        dis_kmax_ind = dis_sort_ind[:, :k_max]

        return dis, dis_kmax_ind

    @staticmethod
    def __get_indices_and_weight(dis_kmax_ind, row, dis, R, k_min):
        """
        Obtains the indices of dis that are in R, outside R, and
        the k_max smallest values. Additionally, obtains the weight
        for points outside the transect region.


        Parameters
        ----------
        dis_kmax_ind : Numpy array
            The indices of dis that represent the k_max smallest
            values
        row : int
            Row index of dis_kmax_ind being considered
        dis : Numpy array
            2D numpy array representing the distance between each
            mesh points and each transect point
        R : float
            Search radius for kriging
        k_min : int
            Denotes the minimum number of data points within the
            search radius.

        Returns
        -------
        R_ind : Numpy array
            Indices of dis that are in R
        R_ind_not : Numpy array
            Indices of dis that are outside R
        sel_ind : Numpy array
            Indices of dis representing the k_max smallest values
        M2_weight : float
            The weight for points outside the transect region.
        """

        # weight for points outside of transect region
        M2_weight = 1.0

        # get indices of dis
        sel_ind = dis_kmax_ind[row, :]

        # get all indices within R
        R_ind = np.argwhere(dis[row, sel_ind] <= R).flatten()

        if len(R_ind) < k_min:
            # get the k_min smallest distance values' indices
            R_ind = np.argsort(dis[row, sel_ind]).flatten()[:k_min]

            # get the indices of sel_ind[R_ind] that are outside of R
            R_ind_not = np.argwhere(dis[row, sel_ind[R_ind]] > R).flatten()

            # TODO: should we change this to how Chu does it?
            # tapered function to handle extrapolation
            M2_weight = np.exp(-np.nanmean(dis[row, sel_ind[R_ind]]) / R)
        else:
            R_ind_not = []

        return R_ind, R_ind_not, sel_ind, M2_weight

    @staticmethod
    def __get_M2_K(x_data, y_data, dis, row, dis_sel_ind,
                   s_v_params, s_v_model):
        """

        Parameters
        ----------
        x_data : Numpy array
            1D array denoting the x coordinates of the data
        y_data : Numpy array
            1D array denoting the y coordinates of the data
        dis : Numpy array
            2D numpy array representing the distance between each
            mesh points and each transect point
        row : int
            Row index of dis_kmax_ind being considered
        dis_sel_ind : Numpy array
            Indices of dis within the search radius
        s_v_params : dict
            Dictionary specifying the parameter values for the
            semi-variogram model.
        s_v_model : function
            A Semi-variogram model from the SemiVariogram class

        Returns
        -------
        M2 : Numpy array
            1D array representing the vector in Kriging
        K : Numpy array
            2D array representing the matrix in Kriging
        """

        # calculate semi-variogram value
        M20 = s_v_model(dis[row, dis_sel_ind], **s_v_params)

        # TODO: put in statements for Objective mapping and
        #  Universal Kriging w/ Linear drift
        M2 = np.concatenate([M20, np.array([1.0])])  # for Ordinary Kriging

        # select those data within the search radius
        x1 = x_data[dis_sel_ind]
        y1 = y_data[dis_sel_ind]

        # compute the distance between the points
        x1_diff = np.subtract.outer(x1, x1)
        y1_diff = np.subtract.outer(y1, y1)
        dis1 = np.sqrt(x1_diff * x1_diff + y1_diff * y1_diff)

        # calculate semi-variogram value
        K0 = s_v_model(dis1, **s_v_params)

        # TODO: put in statements for Objective mapping and
        #  Universal Kriging w/ Linear drift
        # Add column and row of ones for Ordinary Kriging
        K = np.concatenate([K0, np.ones((len(x1), 1))], axis=1)
        K = np.concatenate([K, np.ones((1, len(x1) + 1))], axis=0)

        # do an inplace fill of diagonal
        np.fill_diagonal(K, 0.0)

        return M2, K

    @staticmethod
    def __compute_lambda_weights(M2, K, ratio):
        """
        Computes the lambda weights of Kriging.

        Parameters
        ----------
        M2 : Numpy array
            1D array representing the vector in Kriging
        K : Numpy array
            2D array representing the matrix in Kriging
        ratio : float
            Acceptable ratio for the singular values divided
            by the largest singular value.

        Returns
        -------
        lamb : Numpy array
            Lambda weights of Kriging
        """

        # compute SVD
        u, s, vh = np.linalg.svd(K, full_matrices=True)

        kindx = np.argwhere(np.abs(s / s[0]) > ratio).flatten()

        s_inv = 1.0 / s[kindx]

        k_inv = np.matmul(vh.T[:, kindx], np.diag(s_inv))
        k_inv = np.matmul(k_inv, u[:, kindx].T)

        lamb = np.dot(k_inv, M2)

        return lamb

    @staticmethod
    def __compute_kriging_vals(field_data, M2, lamb, M2_weight,
                               R_ind, R_ind_not, dis_sel_ind):
        """
        Computes the Kriged values, Kriging variance, and
        Kriging sample variance.

        Parameters
        ----------
        field_data : Numpy array
            1D array denoting the field values at the (x ,y)
            coordinates of the data (e.g. biomass density).
        M2 : Numpy array
            1D array representing the vector in Kriging
        lamb : Numpy array
            Lambda weights of Kriging
        M2_weight : float
            The weight for points outside the transect region.
        R_ind : Numpy array
            Indices of dis that are in R
        R_ind_not : Numpy array
            Indices of dis that are outside R
        dis_sel_ind : Numpy array
            Indices of dis within the search radius

        Returns
        -------
        ep_val : float
            Kriging variance
        eps_val : float
            Kriging sample variance
        vp_val : float
            Kriged value
        """

        # obtain field values for indices within R
        if len(R_ind_not) > 0:
            field_vals = field_data[dis_sel_ind]
            field_vals[R_ind_not] = 0.0  # accounts for less than k_min points
        else:
            field_vals = field_data[dis_sel_ind]

        vp_val = np.nansum(lamb[:len(R_ind)] * field_vals) * M2_weight
        ep_val = np.nansum(lamb * M2)

        if abs(vp_val) < np.finfo(float).eps:
            eps_val = np.nan
        else:
            field_var = np.nanvar(field_vals, ddof=1)
            eps_val = np.sqrt(ep_val * field_var) / abs(vp_val)

        # TODO: Do we count the anomalies like Chu does?

        return ep_val, eps_val, vp_val

    def run_kriging(self, x_mesh, x_data, y_mesh, y_data, field_data,
                    k_max, k_min, R, ratio, s_v_params, s_v_model):
        """
        Runs Kriging using the provided mesh and data. Currently,
        only Ordinary Kriging has been implemented.

        Parameters
        ----------
        x_mesh : Numpy array
            1D array denoting the x coordinates of the mesh
        x_data : Numpy array
            1D array denoting the x coordinates of the data
        y_mesh : Numpy array
            1D array denoting the y coordinates of the mesh
        y_data : Numpy array
            1D array denoting the y coordinates of the data
        field_data : Numpy array
            1D array denoting the field values at the (x ,y)
            coordinates of the data (e.g. biomass density).
        k_max : int
            Denotes the maximum number of data points within the
            search radius.
        k_min : int
            Denotes the minimum number of data points within the
            search radius.
        R : float
            Search radius for kriging
        ratio : float
            Acceptable ratio for the singular values divided
            by the largest singular value.
        s_v_params : dict
            Dictionary specifying the parameter values for the
            semi-variogram model.
        s_v_model : function
            A Semi-variogram model from the SemiVariogram class

        Returns
        -------
        ep_arr : 1D Numpy array
            Kriging variance for each mesh coordinate
        eps_arr : 1D Numpy array
            Kriging sample variance for each mesh coordinate
        vp_arr : 1D Numpy array
            Kriged value for each mesh coordinate
        """

        # TODO: think about making kriging mesh class an input

        dis, dis_kmax_ind = self.__compute_k_smallest_distances(x_mesh, x_data,
                                                                y_mesh, y_data, k_max)

        ep_arr = np.empty(dis_kmax_ind.shape[0])
        eps_arr = np.empty(dis_kmax_ind.shape[0])
        vp_arr = np.empty(dis_kmax_ind.shape[0])

        # TODO: look into parallelizing this for loop
        # does Ordinary Kriging, follow Journel and Huijbregts, p. 307
        for row in range(dis_kmax_ind.shape[0]):

            R_ind, R_ind_not, sel_ind, M2_weight = \
                self.__get_indices_and_weight(dis_kmax_ind, row, dis, R, k_min)

            # indices of dis within the search radius
            dis_sel_ind = sel_ind[R_ind]

            M2, K = self.__get_M2_K(x_data, y_data, dis, row, dis_sel_ind,
                                    s_v_params, s_v_model)

            lamb = self.__compute_lambda_weights(M2, K, ratio)

            ep_val, eps_val, vp_val = self.__compute_kriging_vals(field_data, M2, lamb,
                                                                  M2_weight, R_ind,
                                                                  R_ind_not, dis_sel_ind)

            ep_arr[row] = ep_val
            eps_arr[row] = eps_val
            vp_arr[row] = vp_val

        # zero-out all vp values that are nan or negative # TODO: Is this necessary?
        neg_nan_ind = np.argwhere((vp_arr < 0) | np.isnan(vp_arr)).flatten()
        vp_arr[neg_nan_ind] = 0.0

        return ep_arr, eps_arr, vp_arr

    def plot_kriging_results(self, x_mesh, y_mesh, krig_val):

        # TODO: formalize this function more (add kwargs, doc string, ...)

        fmap = folium.Map(location=[44.61, -125.66], zoom_start=4)

        data = np.hstack([x_mesh[:, None], y_mesh[:, None]])

        # colormap = cm.LinearColormap(colors=['#000080', '#3385ff'],
        #                              vmin=0.0, vmax=np.max(vp_arr))

        colormap = cm.LinearColormap(colors=['#3385ff', '#FF0000'],
                                     vmin=0.0, vmax=np.max(krig_val))

        for i in range(len(krig_val)):
            coordinates = (data[i, 0], data[i, 1])
            color_val = colormap(krig_val[i])

            # Place the markers with a color value
            fmap.add_child(folium.CircleMarker(location=coordinates,
                                               radius=1,
                                               color=color_val))

        fmap.add_child(colormap)

        return fmap


