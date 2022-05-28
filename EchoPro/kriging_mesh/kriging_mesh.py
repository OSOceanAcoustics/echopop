import numpy as np
import pandas as pd
import warnings
import sys
import folium
from matplotlib.colors import to_hex
from matplotlib import cm
import folium
import geopandas
from shapely.geometry import Polygon, Point, MultiPoint, MultiPolygon
from shapely.geometry import LineString
import shapely
from shapely.ops import unary_union


class KrigingMesh:
    """
    This class loads the mesh data, provides
    functions to manipulate the mesh, and
    provides functions to plot the mesh.

    Parameters
    ----------
    EPro : EchoPro object
        An initialized EchoPro object. Note that any change to
        self.EPro will also change this object.
    """

    def __init__(self, EPro = None):

        self.EPro = EPro

        self.__load_mesh()
        self.__load_smoothed_contour()

    def __load_mesh(self):
        """
        The returned Dataframe is a
        GeoPandas Dataframe
        """

        df = pd.read_excel(
            self.EPro.params['data_root_dir'] + self.EPro.params['filename_grid_cell'],
            sheet_name='krigedgrid2_5nm_forChu')  # TODO: make the sheet name an input

        # obtaining those columns that are required
        df = df[['Latitude of centroid', 'Longitude of centroid', 'Area (km^2)', 'Cell portion']].copy()

        # set data types of dataframe
        df = df.astype({'Latitude of centroid': float,
                        'Longitude of centroid': float,
                        'Area (km^2)': float,
                        'Cell portion': np.float64})

        df = geopandas.GeoDataFrame(df,
                                    geometry=geopandas.points_from_xy(
                                        df['Longitude of centroid'],
                                        df['Latitude of centroid']
                                    )
                                    )

        self.mesh_gdf = df

    def __load_smoothed_contour(self):
        """
        The returned Dataframe is a
        GeoPandas Dataframe
        """

        df = pd.read_excel(
            self.EPro.params['data_root_dir'] + self.EPro.params['filename_smoothed_contour'],
            sheet_name='Smoothing_EasyKrig')  # TODO: make the sheet name an input

        # obtaining those columns that are required
        df = df[['Latitude', 'Longitude']].copy()

        # set data types of dataframe
        df = df.astype({'Latitude': float,
                        'Longitude': float})

        df = geopandas.GeoDataFrame(df,
                                    geometry=geopandas.points_from_xy(df.Longitude, df.Latitude))

        self.smoothed_contour_gdf = df

    def get_polygon_of_transects(self, gdf, gdf_tran_mean, n_close, nm_to_buffer=1.25):
        """
        This function constructs a polygon that contains
        all transects.

        How it works: Using the mean latitude, longitude
        point of each transect, we find the `nsmall` closest
        transects to each transect and construct a Polygon
        of these transects. Then, we compute the convex
        hull of this Polygon, which creates the smallest
        convex Polygon containing all the points in the object.
        Lastly, we take the unary union of all constructed
        convex Polygons.

        Parameters
        ----------
        gdf : GeoPandas Dataframe
            All transect points with index transect and
            geometry column of shapely points.
        gdf_tran_mean : GeoPandas Dataframe
            Contains the mean coordinate point (latitude, longitude)
            of each transect with index representing the transect
            and geometry column of shapely points.
        n_close : int
            The number of closest transects to include in the Polygon.
            This value includes the transect under consideration. Thus,
            n_close should be greater than or equal to 2.
        nm_to_buffer : float
            The number of nautical miles to buffer the constructed Polygon

        Returns
        -------
        Polygon that contains all transect data

        Note: The connectivity of this polygon
        is determined by n_close.
        """

        transect_polygons = []
        for transect in gdf_tran_mean.index:
            closest_trans = gdf_tran_mean.geometry.distance(
                gdf_tran_mean.loc[transect,
                                  'geometry']).nsmallest(n_close)

            full_pol = Polygon(list(gdf.loc[closest_trans.index, "geometry"]))
            transect_polygons.append(full_pol.convex_hull)

        pol = unary_union(transect_polygons)

        # 1 degree of latitude equals 60nm
        # buffer Polygon by a number of nm
        buf_val = (1.0 / 60.0) * nm_to_buffer

        return pol.buffer(buf_val)

    def reduce_grid_points(self,
                           bio_df, bio_df_lat_name, bio_df_lon_name,
                           gdf, gdf_lat_name, gdf_lon_name,
                           n_close):
        """

        Parameters
        ----------
        bio_df : Pandas Dataframe
            Contains the position of the trawl points and has columns
            `bio_data_lat_name`, `bio_data_lon_name`, and index
            `Transects`.
        bio_df_lat_name : str
            String specifying the name of the latitude column for bio_df
        bio_df_lon_name : str
            String specifying the name of the longitude column for bio_df
        gdf : GeoPandas Dataframe
            Contains the
        gdf_lat_name : str
            String specifying the name of the latitude column for gdf
        gdf_lon_name : str
            String specifying the name of the longitude column for gdf
        n_close :
            The number of closest transects to include in the Polygon
            construction. This value includes the transect under
            consideration. Thus, it should be greater than or equal to 2.

        Returns
        -------

        Note: The connectivity of the final polygon
        is determined by n_close.
        """

        bio_gdf = geopandas.GeoDataFrame(bio_df,
                                         geometry=geopandas.points_from_xy(
                                             bio_df[bio_df_lon_name],
                                             bio_df[bio_df_lat_name]))

        # get the mean latitude and longitude values of bio_df along the transects
        bio_df_tran_mean = bio_gdf[[bio_df_lon_name, bio_df_lat_name]].groupby(level=0).mean()

        # get GeoPandas form of bio_df_tran_mean
        bio_gdf_tran_mean = geopandas.GeoDataFrame(bio_df_tran_mean,
                                                   geometry=geopandas.points_from_xy(
                                                       bio_df[bio_df_lon_name],
                                                       bio_df[bio_df_lat_name]))

        self.get_polygon_of_transects(bio_gdf, bio_gdf_tran_mean, n_close)


    def get_folium_map(self,
                       map_kwargs={
                           'location': [44.61, -125.66],
                           'zoom_start': 4, 'tiles':
                               'CartoDB positron'}):
        """
        grabs the folium map with specifications
        Returns
        -------

        """

        fmap = folium.Map( **map_kwargs)

        return fmap

    def plot_points(self, fmap, geo_df, cmap_column=None, color='hex'):
        """
        Parameters
        ----------
        fmap
            Folium map to plot the points on
        geo_df : GeoPandas Dataframe
            Contains the data to be plotted with a geometry column
        color : str
            The color of the markers representing the points. If
            color='hex', then a matplotlib color map will be created.
        cmap_column : str
            Column of geo_df that the colormap should correspond
            to. This is only used if color='hex'.
        # TODO: make kwargs for folium.Map and folium.CircleMarker

        Returns
        -------
        folium.Map
        """

        if color == 'hex':
            uniq_vals = geo_df[cmap_column].unique()
            cmap = cm.get_cmap('viridis', len(uniq_vals))
            hex_color_options = {rgb[0]: to_hex(rgb[1])
                                 for rgb in zip(uniq_vals, cmap(uniq_vals))}

        for index, row in geo_df.iterrows():

            coordinates = (row.geometry.xy[1][0], row.geometry.xy[0][0])

            if color == 'hex':
                color_val = hex_color_options[row[cmap_column]]
            else:
                color_val = color

            # Place the markers with the popup labels and data
            fmap.add_child(folium.CircleMarker(location=coordinates,
                                               radius=1,
                                               color=color_val))

        return fmap
