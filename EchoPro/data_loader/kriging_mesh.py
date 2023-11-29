from typing import Tuple, Union

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy import interpolate
from shapely.geometry import Polygon
from shapely.ops import unary_union

from ..utils.input_checks_read import check_and_read


class KrigingMesh:
    """
    This class loads the mesh data, provides
    functions to manipulate the mesh, and
    provides functions to plot the mesh.

    Parameters
    ----------
    survey : Survey
        An initialized ``Survey`` object.

    Notes
    -----
    Any change to ``self.survey`` will also change
    the input survey object.
    """

    def __init__(self, survey=None):

        self.survey = survey

        # expected columns for the mesh Dataframe
        self.mesh_cols_types = {
            "centroid_latitude": float,
            "centroid_longitude": float,
            "fraction_cell_in_polygon": np.float64,
        }

        # expected columns for the smoothed contour Dataframe
        self.contour_cols_types = {"latitude": float, "longitude": float}

        # initialize mesh parameters
        self.transformed_transect_df = None
        self.transect_d_x = None
        self.transect_d_y = None
        self.transformed_mesh_df = None

        self._load_mesh()
        self._load_smoothed_contour()

    def _load_mesh(self) -> None:
        """
        Loads the full mesh of the region being considered.
        Action is completed by reading in an Excel file provided
        by the user defined parameter ``'mesh_filename'``.
        Finally, constructs the GeoPandas Dataframe representing
        the full mesh and assigns it as the class variable ``mesh_gdf``.
        """

        df = check_and_read("kriging/mesh", self.mesh_cols_types, self.survey.params)

        # construct geopandas DataFrame to simplify downstream processes
        gdf = gpd.GeoDataFrame(
            df,
            geometry=gpd.points_from_xy(
                df["centroid_longitude"], df["centroid_latitude"]
            ),
        )

        # assign class variable
        self.mesh_gdf = gdf

    def _load_smoothed_contour(self) -> None:
        """
        Loads the smoothed contour of the region being considered.
        Action is completed by reading in an Excel file determined
        by the user defined parameter ``'smoothed_contour_filename'``.
        Finally, constructs the GeoPandas Dataframe representing
        the smoothed contour and assigns it as the class variable
        ``smoothed_contour_gdf``.
        """

        df = check_and_read(
            "kriging/smoothed_contour",
            self.contour_cols_types,
            self.survey.params
        )

        # construct geopandas DataFrame to simplify downstream processes
        df = gpd.GeoDataFrame(
            df, geometry=gpd.points_from_xy(df.longitude, df.latitude)
        )

        # assign class variable
        self.smoothed_contour_gdf = df

    @staticmethod
    def _get_coordinate_mean(
        df: Union[gpd.GeoDataFrame, pd.DataFrame]
    ) -> gpd.GeoDataFrame:
        """
        Creates a GeoPandas Dataframe representing
        the coordinate (latitude and longitude) mean
        of the provided Dataframe.

        Parameters
        ----------
        df : Union[gpd.GeoDataFrame, pd.DataFrame]
            All transect points with index transect and
            columns: Latitude, Longitude, and geometry.

        Returns
        -------
        gpd.GeoDataFrame
            Dataframe containing the mean coordinate
            point (latitude, longitude) of each transect with
            index representing the transect and geometry column
            of shapely points.
        """

        # get the mean latitude and longitude based on the transects
        df_tran_mean = df[["latitude", "longitude"]].groupby(level=0).mean()

        return gpd.GeoDataFrame(
            df_tran_mean,
            geometry=gpd.points_from_xy(df_tran_mean.longitude, df_tran_mean.latitude),
        )

    def get_polygon_of_transects(
        self, gdf: gpd.GeoDataFrame, n_close: int, nm_to_buffer: float = 1.25
    ) -> Polygon:
        """
        This function constructs a polygon that contains
        all transects.

        Parameters
        ----------
        gdf : gpd.GeoDataFrame
            All transect points with index transect and
            columns: Latitude, Longitude, and geometry
        n_close : int
            The number of closest transects to include in the Polygon.
            This value includes the transect under consideration. Thus,
            n_close should be greater than or equal to 2.
        nm_to_buffer : float
            The number of nautical miles to buffer the constructed Polygon

        Returns
        -------
        Polygon
            Constructed Polygon that contains all transect data

        Notes
        -----
        How it works: Using the mean (latitude, longitude)
        point of each transect, we find the ``n_close`` closest
        transects to each transect and construct a Polygon
        of these transects. Then, we compute the convex
        hull of this Polygon, which creates the smallest
        convex Polygon containing all the points in the object.
        Lastly, we take the unary union of all constructed
        convex Polygons.

        The connectivity of this polygon is determined by ``n_close``.
        """

        gdf_tran_mean = self._get_coordinate_mean(gdf)

        # for each transect construct the smallest convex
        # Polygon containing all the points in the transect
        transect_polygons = []
        for transect in gdf_tran_mean.index:

            # obtain n_close closest transects
            closest_trans = gdf_tran_mean.geometry.distance(
                gdf_tran_mean.loc[transect, "geometry"]
            ).nsmallest(n_close)

            # create polygon encasing closest_trans
            full_pol = Polygon(list(gdf.loc[closest_trans.index, "geometry"]))
            transect_polygons.append(full_pol.convex_hull)

        # obtain Polygon surrounding all transects
        pol = unary_union(transect_polygons)

        # 1 degree of latitude equals 60nm
        buf_val = (1.0 / 60.0) * nm_to_buffer

        # buffer Polygon by a number of nm
        return pol.buffer(buf_val)

    def reduce_grid_points(self, transect_polygon: Polygon) -> gpd.GeoDataFrame:
        """
        Reduces the full mesh points provided to the ``KrigingMesh``
        class by selecting those points that are within the
        transect polygon.

        Parameters
        ----------
        transect_polygon : Polygon
            A Polygon that contains all transect data.

        Returns
        -------
        gpd.GeoDataFrame
            Dataframe representing the reduced mesh points
        """

        # get bool mask of points that are within the polygon
        in_poly = self.mesh_gdf["geometry"].within(transect_polygon)

        # select gdf rows based on bool mask
        return self.mesh_gdf.loc[in_poly].copy()

    def align_longitude(
        self, gdf: gpd.GeoDataFrame, lon_ref: float = -124.78338
    ) -> gpd.GeoDataFrame:
        """
        This function applies a transformation to the
        longitude column of the provided Dataframe so that
        the anisotropic signature of the animal biomass
        distribution can be approximately characterized by two
        perpendicular (principal) correlation scales: one is
        along the isobaths and the other is across the isobaths.

        Parameters
        ----------
        gdf : gpd.GeoDataFrame
            Dataframe with columns specifying Latitude, Longitude,
            and geometry.
        lon_ref : float
            An arbitrary scalar, or a reference longitude
            (e.g., the mean longitude of the 200m isobath)

        Returns
        -------
        transformed_gdf: gpd.GeoDataFrame
            A copy of ``gdf`` with the Longitude column transformed
            according to the interpolated (Latitude, Longitude) values
            pulled from the file specified by the parameter
            ``smoothed_contour_filename``.

        Notes
        -----
        Extreme caution should be used here. This transformation
        was specifically designed for a NWFSC application!
        """

        # construct an interpolation between points
        f = interpolate.interp1d(
            self.smoothed_contour_gdf["latitude"],
            self.smoothed_contour_gdf["longitude"],
            kind="linear",
            bounds_error=False,
        )

        # TODO: do we need to drop NaNs after interpolating?
        #  Investigate this further.

        # apply longitude transformation and store values
        transformed_gdf = gdf.copy()
        transformed_gdf["longitude_transformed"] = (
            gdf.geometry.x - f(gdf.geometry.y) + lon_ref
        )
        transformed_gdf["geometry"] = gpd.points_from_xy(
            transformed_gdf["longitude_transformed"], gdf.geometry.y
        )

        return transformed_gdf

    @staticmethod
    def apply_distance_transformation(
        gdf: gpd.GeoDataFrame,
        d_x: float,
        d_y: float,
        x_offset: float = -124.78338,
        y_offset: float = 45.0,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Transforms the input coordinates from degrees
        to distance i.e. puts them into a distance
        coordinate system.

        Parameters
        ----------
        gdf : gpd.GeoDataFrame
            Dataframe with geometry column.
        d_x : float
            the distance between the maximum longitude
            value and the minimum longitude value.
        d_y : float
            the distance between the maximum latitude
            value and the minimum latitude value
        x_offset : float
            An arbitrary scalar, or a reference longitude
            (e.g., the mean longitude of the 200m isobath)
        y_offset : float
            An arbitrary scalar, or a reference latitude
            (e.g., the mean latitude of the 200m isobath)

        Returns
        -------
        x : np.ndarray
            1D array representing the x coordinate of the distance
            coordinate system
        y : np.ndarray
            1D array representing the y coordinate of the distance
            coordinate system

        Notes
        -----
        Extreme caution should be used here. This transformation
        was specifically designed for a NWFSC application!
        """

        # constant that converts degrees to radians
        DEG2RAD = np.pi / 180.0

        # initial coordinates without distance transformation
        x = gdf.geometry.x
        y = gdf.geometry.y

        # transform the coordinates so they are in terms of distance
        x = np.cos(DEG2RAD * y) * (x - x_offset) / d_x
        y = (y - y_offset) / d_y

        return x.values.flatten(), y.values.flatten()

    def _transform_transect_data(
        self,
        lon_ref: float = -124.78338,
        x_offset: float = -124.78338,
        y_offset: float = 45.0,
    ) -> None:
        """
        Applies a coordinate transformation to ``survey.bio_calc.transect_results_gdf``
        by first aligning the longitude along the smoothed contour data
        specified by the configuration parameter ``'smoothed_contour_filename'``
        and then transforming the input coordinates from degrees to distance.

        Parameters
        ----------
        lon_ref : float
            An arbitrary scalar, or a reference longitude
            (e.g., the mean longitude of the 200m isobath)
            used in aligning the longitude
        x_offset : float
            An arbitrary scalar, or a reference longitude
            (e.g., the mean longitude of the 200m isobath)
            used in transforming from degrees to distance
        y_offset : float
            An arbitrary scalar, or a reference latitude
            (e.g., the mean latitude of the 200m isobath)
            used in transforming from degrees to distance

        Notes
        -----
        This function constructs the following class variables:
        - ``transformed_transect_df`` pd.Dataframe representing the
        transformed transect data
        - ``transect_d_x`` the distance between the maximum longitude
        value and the minimum longitude value (after aligning the longitude)
        - ``transect_d_y`` the distance between the maximum latitude
        value and the minimum latitude value (after aligning the longitude)
        """

        if isinstance(self.survey.bio_calc.transect_results_gdf, gpd.GeoDataFrame):
            # apply transformations to transect points
            transect_df = self.align_longitude(
                self.survey.bio_calc.transect_results_gdf, lon_ref
            )

            # compute distances for each transect
            d_x = transect_df.geometry.x.max() - transect_df.geometry.x.min()
            d_y = transect_df.geometry.y.max() - transect_df.geometry.y.min()

            x_transect, y_transect = self.apply_distance_transformation(
                transect_df, d_x, d_y, x_offset, y_offset
            )

            # store transformed points
            transect_df["x_transect"] = x_transect
            transect_df["y_transect"] = y_transect
            self.transformed_transect_df = transect_df

            # store distance information
            self.transect_d_x = d_x
            self.transect_d_y = d_y
        else:
            raise RuntimeError(
                "survey.bio_calc.transect_results_gdf has not been constructed yet. One "
                "must compute the biomass density before running this function!"
            )

    def apply_coordinate_transformation(
        self,
        coord_type: str = "transect",
        lon_ref: float = -124.78338,
        x_offset: float = -124.78338,
        y_offset: float = 45.0,
    ) -> None:
        """
        Applies a coordinate transformation to either ``survey.bio_calc.transect_results_gdf``
        or ``self.mesh_gdf`` by first aligning the longitude along the
        smoothed contour data specified by the configuration
        parameter ``'smoothed_contour_filename'`` and then
        transforming the input coordinates from degrees to distance.

        Parameters
        ----------
        coord_type : str
            The type of coordinate points to transform.
            Possible options:

            - ``'transect'`` specifies that one should
              copy and transform ``survey.bio_calc.transect_results_gdf``
            - ``'mesh'`` specifies that one should copy '
              and transform `self.mesh_gdf``

        lon_ref : float
            An arbitrary scalar, or a reference longitude
            (e.g., the mean longitude of the 200m isobath)
            used in aligning the longitude
        x_offset : float
            An arbitrary scalar, or a reference longitude
            (e.g., the mean longitude of the 200m isobath)
            used in transforming from degrees to distance
        y_offset : float
            An arbitrary scalar, or a reference latitude
            (e.g., the mean latitude of the 200m isobath)
            used in transforming from degrees to distance

        Notes
        -----
        Extreme caution should be used here. This transformation
        was specifically designed for a NWFSC application!

        This function constructs the following class variables
        for each ``coord_type``:

        - If ``coord_type='transect'``
            - ``transformed_transect_df`` DataFrame representing the
              transformed transect data
            - ``transect_d_x`` the distance between the maximum longitude
              value and the minimum longitude value (after aligning the
              longitude) for transect data
            - ``transect_d_y`` the distance between the maximum latitude
              value and the minimum latitude value (after aligning the
              longitude) for the transect data

        - If ``coord_type='mesh'``
            - ``transformed_mesh_df`` DataFrame representing the
              transformed mesh data
            - Additionally, if ``coord_type='transect'`` has not been
              run, this function will create all class variables associated
              with this input.

        """

        if coord_type == "transect":

            self._transform_transect_data(lon_ref, x_offset, y_offset)

        elif coord_type == "mesh":

            if not self.transect_d_x:
                self._transform_transect_data(lon_ref, x_offset, y_offset)

            # apply transformations to mesh points
            mesh_df = self.align_longitude(self.mesh_gdf)
            x_mesh, y_mesh = self.apply_distance_transformation(
                mesh_df, self.transect_d_x, self.transect_d_y, x_offset, y_offset
            )

            # store transformed mesh for downstream processes
            mesh_df["x_mesh"] = x_mesh
            mesh_df["y_mesh"] = y_mesh
            self.transformed_mesh_df = mesh_df
        else:
            raise ValueError("Unrecognized coordinate type.")
