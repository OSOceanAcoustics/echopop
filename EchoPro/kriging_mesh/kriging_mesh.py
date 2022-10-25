import numpy as np
import pandas as pd
from matplotlib.colors import to_hex
from matplotlib import cm
import folium
import geopandas as gpd
from shapely.geometry import Polygon
from shapely.ops import unary_union
from scipy import interpolate
from typing import Union, Tuple


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

        # Default parameters for the folium map
        self.folium_map_kwargs = {
            'location': [44.61, -125.66],
            'zoom_start': 4,
            'tiles': 'CartoDB positron'
        }

        # expected columns for the mesh Dataframe
        self.mesh_cols = {'Latitude of centroid', 'Longitude of centroid', 'Area (km^2)', 'Cell portion'}

        # expected columns for the smoothed contour Dataframe
        self.contour_cols = {'Latitude', 'Longitude'}

        # initialize mesh parameters
        self.transformed_transect_df = None
        self.transect_d_x = None
        self.transect_d_y = None
        self.transformed_mesh_df = None

        self._load_mesh()
        self._load_smoothed_contour()

    def _check_mesh_df(self, mesh_df: pd.DataFrame) -> None:
        """
        Ensures that the appropriate columns are
        contained in the mesh Dataframe.

        TODO: should we add more in-depth checks here?
        """

        if len(set(mesh_df.columns).intersection(self.mesh_cols)) != len(self.mesh_cols):
            raise NameError("Mesh dataframe does not contain all expected columns!")

    def _check_smoothed_contour_df(self, contour_df: pd.DataFrame) -> None:
        """
        Ensures that the appropriate columns are
        contained in the smoothed contour Dataframe.

        TODO: should we add more in-depth checks here?
        """

        if len(set(contour_df.columns).intersection(self.contour_cols)) != len(self.contour_cols):
            raise NameError("Smoothed contour dataframe does not contain all expected columns!")

    def _load_mesh(self) -> None:
        """
        Loads the full mesh of the region being considered.
        Action is completed by reading in an Excel file provided
        by the user defined parameter ``'mesh_filename'``.
        Finally, constructs the GeoPandas Dataframe representing
        the full mesh and assigns it as the class variable ``mesh_gdf``.
        """

        df = pd.read_excel(self.survey.params['data_root_dir'] + self.survey.params['mesh_filename'],
                           sheet_name=self.survey.params['mesh_sheetname'])
        self._check_mesh_df(df)

        # obtaining those columns that are required
        df = df[['Latitude of centroid', 'Longitude of centroid', 'Area (km^2)', 'Cell portion']].copy()

        # set data types of dataframe
        df = df.astype({'Latitude of centroid': float,
                        'Longitude of centroid': float,
                        'Area (km^2)': float,
                        'Cell portion': np.float64})

        # construct geopandas DataFrame to simplify downstream processes
        gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['Longitude of centroid'],
                                                               df['Latitude of centroid']))

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

        df = pd.read_excel(self.survey.params['data_root_dir'] + self.survey.params['smoothed_contour_filename'],
                           sheet_name=self.survey.params['smoothed_contour_sheetname'])
        self._check_smoothed_contour_df(df)

        # obtaining those columns that are required
        df = df[['Latitude', 'Longitude']].copy()

        # set data types of dataframe
        df = df.astype({'Latitude': float,
                        'Longitude': float})

        # construct geopandas DataFrame to simplify downstream processes
        df = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.Longitude, df.Latitude))

        # assign class variable
        self.smoothed_contour_gdf = df

    @staticmethod
    def _get_coordinate_mean(df: Union[gpd.GeoDataFrame, pd.DataFrame]) -> gpd.GeoDataFrame:
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
        df_tran_mean = df[["Latitude", "Longitude"]].groupby(level=0).mean()

        return gpd.GeoDataFrame(df_tran_mean,
                                geometry=gpd.points_from_xy(df_tran_mean.Longitude,
                                                            df_tran_mean.Latitude))

    def get_polygon_of_transects(self, gdf: gpd.GeoDataFrame,
                                 n_close: int, nm_to_buffer: float = 1.25) -> Polygon:
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
                gdf_tran_mean.loc[transect, 'geometry']).nsmallest(n_close)

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
        in_poly = self.mesh_gdf['geometry'].within(transect_polygon)

        # select gdf rows based on bool mask
        return self.mesh_gdf.loc[in_poly].copy()

    def align_longitude(self, gdf: gpd.GeoDataFrame,
                        lon_ref: float = -124.78338) -> gpd.GeoDataFrame:
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
        f = interpolate.interp1d(self.smoothed_contour_gdf['Latitude'],
                                 self.smoothed_contour_gdf['Longitude'],
                                 kind='linear', bounds_error=False)

        # TODO: do we need to drop NaNs after interpolating?
        #  Investigate this further.

        # apply longitude transformation and store values
        transformed_gdf = gdf.copy()
        transformed_gdf['longitude_transformed'] = gdf.geometry.x - f(gdf.geometry.y) + lon_ref
        transformed_gdf['geometry'] = gpd.points_from_xy(
            transformed_gdf['longitude_transformed'],
            gdf.geometry.y
        )

        return transformed_gdf

    @staticmethod
    def apply_distance_transformation(gdf: gpd.GeoDataFrame, d_x: float,
                                      d_y: float, x_offset: float = -124.78338,
                                      y_offset: float = 45.0) -> Tuple[np.ndarray, np.ndarray]:
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

    def _transform_transect_data(self, lon_ref: float = -124.78338,
                                 x_offset: float = -124.78338,
                                 y_offset: float = 45.0) -> None:
        """
        Applies a coordinate transformation to ``survey.bio_calc.final_biomass_table``
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

        if isinstance(self.survey.bio_calc.final_biomass_table, gpd.GeoDataFrame):
            # apply transformations to transect points
            transect_df = self.align_longitude(self.survey.bio_calc.final_biomass_table, lon_ref)

            # compute distances for each transect
            d_x = transect_df.geometry.x.max() - transect_df.geometry.x.min()
            d_y = transect_df.geometry.y.max() - transect_df.geometry.y.min()

            x_transect, y_transect = self.apply_distance_transformation(transect_df,
                                                                        d_x, d_y,
                                                                        x_offset, y_offset)

            # store transformed points
            transect_df['x_transect'] = x_transect
            transect_df['y_transect'] = y_transect
            self.transformed_transect_df = transect_df

            # store distance information
            self.transect_d_x = d_x
            self.transect_d_y = d_y
        else:
            raise RuntimeError("survey.bio_calc.final_biomass_table has not been constructed yet. One "
                               "must compute the biomass density before running this function!")

    def apply_coordinate_transformation(self, coord_type: str = 'transect',
                                        lon_ref: float = -124.78338,
                                        x_offset: float = -124.78338,
                                        y_offset: float = 45.0) -> None:
        """
        Applies a coordinate transformation to either ``survey.bio_calc.final_biomass_table``
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
            copy and transform ``survey.bio_calc.final_biomass_table``
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

        if coord_type == 'transect':

            self._transform_transect_data(lon_ref, x_offset, y_offset)

        elif coord_type == 'mesh':

            if not self.transect_d_x:
                self._transform_transect_data(lon_ref, x_offset, y_offset)

            # apply transformations to mesh points
            mesh_df = self.align_longitude(self.mesh_gdf)
            x_mesh, y_mesh = self.apply_distance_transformation(mesh_df,
                                                                self.transect_d_x,
                                                                self.transect_d_y,
                                                                x_offset,
                                                                y_offset)

            # store transformed mesh for downstream processes
            mesh_df['x_mesh'] = x_mesh
            mesh_df['y_mesh'] = y_mesh
            self.transformed_mesh_df = mesh_df
        else:
            raise ValueError("Unrecognized coordinate type.")

    def get_folium_map(self, map_kwargs: dict = None) -> folium.Map:
        """
        Grabs the folium map.

        Parameters
        ----------
        map_kwargs : dict
            Dictionary of kwargs for the folium.Map function.

        Returns
        -------
        folium.Map
            A folium map object with provided specifications.

        Notes
        -----
        If ``map_kwargs`` is empty, then the default dictionary
        ``map_kwargs={'location': [44.61, -125.66], 'zoom_start': 4,
        'tiles': 'CartoDB positron'}`` will be used.
        """

        if map_kwargs is None:
            map_kwargs = self.folium_map_kwargs

        return folium.Map(**map_kwargs)

    def plot_points(self, gdf: gpd.GeoDataFrame,
                    fobj: Union[folium.Map, folium.map.FeatureGroup] = None,
                    color: str = 'hex', cmap_column: str = None,
                    marker_kwargs: dict = {}) -> Union[folium.Map, folium.map.FeatureGroup]:
        """
        Allows for a simple way to plot and
        visualize mesh points on a folium map.

        Parameters
        ----------
        gdf : gpd.GeoDataFrame
            Contains a geometry column representing the points to plot
        fobj : Union[folium.Map, folium.map.FeatureGroup]
            Folium object to plot the points on
        color : str
            The color of the markers representing the points. If
            color='hex', then a matplotlib color map will be created.
        cmap_column : str
            Column of gdf that the colormap should correspond
            to. This is only used if color='hex'.
        marker_kwargs : dict
            Dictionary of kwargs that should be provided to
            folium.CircleMarker

        Returns
        -------
        Union[folium.Map, folium.map.FeatureGroup]
            Folium object with points attached to it

        Notes
        -----
        If ``fobj`` is not provided, then a default folium map
        will be created using the class function ``get_folium_map``.
        """

        if fobj is None:
            fobj = self.get_folium_map()

        if color == 'hex':

            if cmap_column in gdf.index.names:
                gdf = gdf.reset_index()
            elif cmap_column not in gdf:
                raise RuntimeError(f"gdf does not contain an index or column with name {cmap_column}!")

            # create a color map using hex values based off a value in gdf
            uniq_vals = gdf[cmap_column].unique()
            cmap = cm.get_cmap('viridis', len(uniq_vals))
            hex_color_options = {rgb[0]: to_hex(rgb[1])
                                 for rgb in zip(uniq_vals, cmap(uniq_vals))}

        # plot each point in gdf on a folium object
        # TODO: Can we do this more efficiently?
        for _, row in gdf.iterrows():

            if color == 'hex':
                color_val = hex_color_options[row[cmap_column]]
            else:
                color_val = color

            # Place the markers with specific color
            fobj.add_child(folium.CircleMarker(location=(row.geometry.y, row.geometry.x),
                                               radius=1,
                                               color=color_val,
                                               **marker_kwargs))

        return fobj

    def plot_layered_points(self) -> folium.Map:
        """
        This function constructs a layered Folium plot.
        The layers correspond to the full set of mesh
        points, the ``final_biomass_table`` points with
        color corresponding to the transect number, and
        the smoothed contour points (e.g. 200m isobath).

        Returns
        -------
        fmap : folium.Map
            A Map with layers corresponding to the
            various points plotted.

        Notes
        -----
        This function is considered a high-level plotting
        tool. Thus, fine control of color and the type of
        points to plot has not been implemented. If one
        wants low-level control for plotting, the class
        function ``plot_points`` should be used.
        """

        fmap = self.get_folium_map()

        # plot mesh points and add them to fmap
        folium_layer = folium.FeatureGroup(name='mesh')
        folium_layer = self.plot_points(self.mesh_gdf, folium_layer, color='gray')
        folium_layer.add_to(fmap)

        # plot the transect points and add them to fmap
        folium_layer = folium.FeatureGroup(name='transects')
        folium_layer = self.plot_points(self.survey.bio_calc.final_biomass_table, folium_layer,
                                        cmap_column='Transect', color='hex')
        folium_layer.add_to(fmap)

        # plot smoothed contour points and add them to fmap
        folium_layer = folium.FeatureGroup(name='smoothed contour')
        folium_layer = self.plot_points(self.smoothed_contour_gdf, folium_layer, color='blue')
        folium_layer.add_to(fmap)

        # add layer control to fmap
        folium.LayerControl().add_to(fmap)

        return fmap
