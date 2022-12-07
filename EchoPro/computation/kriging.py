import numpy as np
import folium
import branca.colormap as cm
from ..numba_modules import nb_subtract_outer, nb_dis_mat
from typing import Callable, Tuple
from ..data_loader import KrigingMesh
import geopandas as gpd


class Kriging:
    """
    This class constructs all data necessary
    for Kriging and performs Kriging

    Parameters
    ----------
    survey : Survey
        An initialized Survey object. Note that any change to
        ``self.survey`` will also change this object.
    k_max: int
        The maximum number of data points within the search radius
    k_min: int
        The minimum number of data points within the search radius
    R: float
        Search radius for Kriging
    ratio: float
        Acceptable ratio for the singular values divided by the largest
        singular value.
    s_v_params: dict
        Dictionary specifying the parameter values for the semi-variogram model.
    s_v_model: Callable
        a Semi-variogram model from the ``SemiVariogram`` class
    """

    def __init__(self, survey=None, k_max: int = None, k_min: int = None,
                 R: float = None, ratio: float = None, s_v_params: dict = None,
                 s_v_model: Callable = None):

        self.survey = survey

        # Kriging parameters
        self.k_max = k_max
        self.k_min = k_min
        self.R = R
        self.ratio = ratio

        # parameters for semi-variogram model
        self.s_v_params = s_v_params

        # grab appropriate semi-variogram model
        self.s_v_model = s_v_model

    def _compute_k_smallest_distances(self, x_mesh: np.ndarray,
                                      x_data: np.ndarray, y_mesh: np.ndarray,
                                      y_data: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Computes the distance between the data
        and the mesh points. Using the distance, it then
        selects the kmax smallest values amongst the
        columns.

        Parameters
        ----------
        x_mesh : np.ndarray
            1D array denoting the x coordinates of the mesh
        x_data : np.ndarray
            1D array denoting the x coordinates of the data
        y_mesh : np.ndarray
            1D array denoting the y coordinates of the mesh
        y_data : np.ndarray
            1D array denoting the y coordinates of the data

        Returns
        -------
        dis : np.ndarray
            2D array representing the distance between each
            mesh point and each transect point
        dis_kmax_ind : np.ndarray
            A 2D array index array representing the kmax
            closest transect points to each mesh point.
        """

        # compute the distance between the mesh points and transect points
        x_diff = nb_subtract_outer(x_mesh, x_data)
        y_diff = nb_subtract_outer(y_mesh, y_data)
        dis = nb_dis_mat(x_diff, y_diff)

        # sort dis up to the kmax smallest elements in each row
        dis_sort_ind = np.argpartition(dis, self.k_max, axis=1)

        # select only the kmax smallest elements in each row
        dis_kmax_ind = dis_sort_ind[:, :self.k_max]

        return dis, dis_kmax_ind

    def _get_indices_and_weight(self, dis_kmax_ind: np.ndarray, row: int,
                                dis: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float]:
        """
        Obtains the indices of ``dis`` that are in ``R``, outside ``R``, and
        the ``k_max`` smallest values. Additionally, obtains the weight
        for points outside the transect region.


        Parameters
        ----------
        dis_kmax_ind : np.ndarray
            The indices of ``dis`` that represent the ``k_max`` smallest
            values
        row : int
            Row index of ``dis_kmax_ind`` being considered
        dis : np.ndarray
            2D numpy array representing the distance between each
            mesh point and each transect point

        Returns
        -------
        R_ind : np.ndarray
            Indices of ``dis`` that are in ``R``
        R_ind_not : np.ndarray
            Indices of ``dis`` that are outside ``R``
        sel_ind : np.ndarray
            Indices of ``dis`` representing the ``k_max`` smallest values
        M2_weight : float
            The weight for points outside the transect region.
        """

        # weight for points outside of transect region
        M2_weight = 1.0

        # get indices of dis
        sel_ind = dis_kmax_ind[row, :]

        # get all indices within R
        R_ind = np.argwhere(dis[row, sel_ind] <= self.R).flatten()

        if len(R_ind) < self.k_min:

            # get the k_min smallest distance values' indices
            R_ind = np.argsort(dis[row, sel_ind]).flatten()[:self.k_min]

            # get the indices of sel_ind[R_ind] that are outside of R
            R_ind_not = np.argwhere(dis[row, sel_ind[R_ind]] > self.R).flatten()

            # TODO: should we change this to how Chu does it?
            # tapered function to handle extrapolation
            M2_weight = np.exp(-np.nanmean(dis[row, sel_ind[R_ind]]) / self.R)
        else:
            R_ind_not = []

        return R_ind, R_ind_not, sel_ind, M2_weight

    def _get_M2_K(self, x_data: np.ndarray, y_data: np.ndarray, dis: np.ndarray,
                  row: int, dis_sel_ind: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """

        Parameters
        ----------
        x_data : np.ndarray
            1D array denoting the x coordinates of the data
        y_data : np.ndarray
            1D array denoting the y coordinates of the data
        dis : np.ndarray
            2D numpy array representing the distance between each
            mesh point and each transect point
        row : int
            Row index of ``dis_kmax_ind`` being considered
        dis_sel_ind : np.ndarray
            Indices of ``dis`` within the search radius

        Returns
        -------
        M2 : np.ndarray
            1D array representing the vector in Kriging
        K : np.ndarray
            2D array representing the matrix in Kriging
        """

        # calculate semi-variogram value
        M20 = self.s_v_model(dis[row, dis_sel_ind], **self.s_v_params)

        # TODO: Should we put in statements for Objective mapping and Universal Kriging w/ Linear drift?
        M2 = np.concatenate([M20, np.array([1.0])])  # for Ordinary Kriging

        # select those data within the search radius
        x1 = x_data[dis_sel_ind]
        y1 = y_data[dis_sel_ind]

        # compute the distance between the points
        x1_diff = np.subtract.outer(x1, x1)
        y1_diff = np.subtract.outer(y1, y1)
        dis1 = np.sqrt(x1_diff * x1_diff + y1_diff * y1_diff)

        # calculate semi-variogram value
        K0 = self.s_v_model(dis1, **self.s_v_params)

        # TODO: Should we put in statements for Objective mapping and Universal Kriging w/ Linear drift?
        # Add column and row of ones for Ordinary Kriging
        K = np.concatenate([K0, np.ones((len(x1), 1))], axis=1)
        K = np.concatenate([K, np.ones((1, len(x1) + 1))], axis=0)

        # do an inplace fill of diagonal
        np.fill_diagonal(K, 0.0)

        return M2, K

    def _compute_lambda_weights(self, M2: np.ndarray, K: np.ndarray) -> np.ndarray:
        """
        Computes the lambda weights of Kriging.

        Parameters
        ----------
        M2 : np.ndarray
            1D array representing the vector in Kriging
        K : np.ndarray
            2D array representing the matrix in Kriging

        Returns
        -------
        lamb : np.ndarray
            Lambda weights computed from Kriging
        """

        # compute SVD
        u, s, vh = np.linalg.svd(K, full_matrices=True)

        kindx = np.argwhere(np.abs(s / s[0]) > self.ratio).flatten()

        s_inv = 1.0 / s[kindx]

        k_inv = np.matmul(vh.T[:, kindx], np.diag(s_inv))
        k_inv = np.matmul(k_inv, u[:, kindx].T)

        lamb = np.dot(k_inv, M2)

        return lamb

    @staticmethod
    def _compute_kriging_vals(field_data: np.ndarray, M2: np.ndarray, lamb: np.ndarray,
                              M2_weight: float, R_ind: np.ndarray, R_ind_not: np.ndarray,
                              dis_sel_ind: np.ndarray) -> Tuple[float, float, float]:
        """
        Computes the Kriged values, Kriging variance, and
        Kriging sample variance.

        Parameters
        ----------
        field_data : np.ndarray
            1D array denoting the field values at the (x ,y)
            coordinates of the data (e.g. biomass density).
        M2 : np.ndarray
            1D array representing the vector in Kriging
        lamb : np.ndarray
            Lambda weights of Kriging
        M2_weight : float
            The weight for points outside the transect region.
        R_ind : np.ndarray
            Indices of dis that are in R
        R_ind_not : np.ndarray
            Indices of dis that are outside R
        dis_sel_ind : np.ndarray
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

        # calculate Kriging value and variance
        vp_val = np.nansum(lamb[:len(R_ind)] * field_vals) * M2_weight
        ep_val = np.nansum(lamb * M2)

        # calculate Kriging sample variance
        if abs(vp_val) < np.finfo(float).eps:
            eps_val = np.nan
        else:
            field_var = np.nanvar(field_vals, ddof=1)
            eps_val = np.sqrt(ep_val * field_var) / abs(vp_val)

        # TODO: Do we count the anomalies like Chu does?

        return ep_val, eps_val, vp_val

    def run_kriging(self, x_mesh: np.ndarray, x_data: np.ndarray,
                    y_mesh: np.ndarray, y_data: np.ndarray,
                    field_data: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        A low-level interface that runs Kriging using the provided
        mesh and data.

        Parameters
        ----------
        x_mesh : np.ndarray
            1D array denoting the x coordinates of the mesh
        x_data : np.ndarray
            1D array denoting the x coordinates of the data
        y_mesh : np.ndarray
            1D array denoting the y coordinates of the mesh
        y_data : np.ndarray
            1D array denoting the y coordinates of the data
        field_data : np.ndarray
            1D array denoting the field values at the (x ,y)
            coordinates of the data (e.g. biomass density).

        Returns
        -------
        ep_arr : np.ndarray
            1D array representing the Kriging variance for each mesh coordinate
        eps_arr : np.ndarray
            1D array representing the Kriging sample variance for each mesh coordinate
        vp_arr : np.ndarray
            1D array representing the Kriged value for each mesh coordinate

        Notes
        -----
        Currently, this routine only runs Ordinary Kriging.
        """

        dis, dis_kmax_ind = self._compute_k_smallest_distances(x_mesh, x_data,
                                                               y_mesh, y_data)

        # initialize arrays that store calculated Kriging values
        ep_arr = np.empty(dis_kmax_ind.shape[0])
        eps_arr = np.empty(dis_kmax_ind.shape[0])
        vp_arr = np.empty(dis_kmax_ind.shape[0])

        # TODO: This loop can be parallelized, if necessary
        # does Ordinary Kriging, follow Journel and Huijbregts, p. 307
        for row in range(dis_kmax_ind.shape[0]):

            R_ind, R_ind_not, sel_ind, M2_weight = self._get_indices_and_weight(dis_kmax_ind, row, dis)

            # indices of dis within the search radius
            dis_sel_ind = sel_ind[R_ind]

            M2, K = self._get_M2_K(x_data, y_data, dis, row, dis_sel_ind)

            lamb = self._compute_lambda_weights(M2, K)

            ep_val, eps_val, vp_val = self._compute_kriging_vals(field_data, M2, lamb,
                                                                 M2_weight, R_ind,
                                                                 R_ind_not, dis_sel_ind)

            # store important calculated values
            ep_arr[row] = ep_val
            eps_arr[row] = eps_val
            vp_arr[row] = vp_val

        # zero-out all vp values that are nan or negative # TODO: Is this necessary?
        neg_nan_ind = np.argwhere((vp_arr < 0) | np.isnan(vp_arr)).flatten()
        vp_arr[neg_nan_ind] = 0.0

        return ep_arr, eps_arr, vp_arr

    def run_biomass_kriging(self, krig_mesh: KrigingMesh) -> None:
        """
        A high-level interface that sets up and runs
        Kriging using the areal biomass density.
        The results are then stored in the ``Survey``
        object as ``krig_results_gdf``.


        Parameters
        ----------
        krig_mesh : KrigingMesh
            Object representing the Kriging mesh

        Notes
        -----
        To run this routine, one must first compute the areal biomass
        density using ``compute_biomass_density``.
        """

        if not isinstance(krig_mesh, KrigingMesh):
            raise ValueError("You must provide a KrigingMesh object!")

        if (not isinstance(self.survey.bio_calc.final_biomass_table, gpd.GeoDataFrame)) \
                and ('areal_biomass_density_adult' not in self.survey.bio_calc.final_biomass_table):
            raise ValueError("The areal biomass density must be calculated before running this routine!")

        ep_arr, eps_arr, vp_arr = self.run_kriging(
            krig_mesh.transformed_mesh_df['x_mesh'].values,
            krig_mesh.transformed_transect_df['x_transect'].values,
            krig_mesh.transformed_mesh_df['y_mesh'].values,
            krig_mesh.transformed_transect_df['y_transect'].values,
            self.survey.bio_calc.final_biomass_table['areal_biomass_density_adult'].values.flatten())

        # collect all important Kriging results
        results_gdf = krig_mesh.mesh_gdf.copy()
        results_gdf['krig_biomass_vp'] = vp_arr
        results_gdf['krig_biomass_ep'] = ep_arr
        results_gdf['krig_biomass_eps'] = eps_arr
        results_gdf["area_calc"] = self.survey.params['kriging_A0'] * results_gdf['fraction_cell_in_polygon']
        results_gdf["krig_biomass_vals"] = results_gdf['krig_biomass_vp'] * results_gdf["area_calc"]

        self.survey.bio_calc.krig_results_gdf = results_gdf

    # @staticmethod
    # def plot_kriging_results(krig_results_gdf: gpd.GeoDataFrame,
    #                          krig_field_name: str) -> folium.Map:
    #     """
    #     Constructs a Folium plot depicting the ``krig_field_name``
    #     values at each mesh point.
    #
    #     Parameters
    #     ----------
    #     krig_results_gdf: gpd.GeoDataFrame
    #         Dataframe containing a geometry column that holds the
    #         mesh coordinates and the column ``krig_field_name``
    #     krig_field_name: str
    #         The name of the column in ``krig_results_gdf`` containing
    #         the Kriging values to plot at each mesh point
    #
    #     Returns
    #     -------
    #     fmap : folium.Map
    #         A Folium plot with the ``krig_field_name`` values at each
    #         mesh point.
    #     """
    #
    #     # create folium map
    #     fmap = folium.Map(location=[44.61, -125.66], zoom_start=4)
    #
    #     # collect the appropriate data from the input Dataframe
    #     x_mesh = krig_results_gdf.geometry.x.values
    #     y_mesh = krig_results_gdf.geometry.y.values
    #     krig_val = krig_results_gdf[krig_field_name].values
    #
    #     # create a colormap for the values
    #     colormap = cm.LinearColormap(colors=['#3385ff', '#FF0000'],
    #                                  vmin=0.0, vmax=np.max(krig_val))
    #
    #     # plot each mesh point and the corresponding krig_field_name value
    #     for i in range(len(krig_val)):
    #
    #         # get the color of the marker we want to plot
    #         color_val = colormap(krig_val[i])
    #
    #         # get coordinate of marker we are plotting
    #         coordinates = (y_mesh[i], x_mesh[i])
    #
    #         # Place the markers with a color value
    #         fmap.add_child(folium.CircleMarker(location=coordinates,
    #                                            radius=1,
    #                                            color=color_val))
    #
    #     fmap.add_child(colormap)  # adds color bar to map
    #
    #     return fmap

    # TODO: Move plotting to visualization sub-package


