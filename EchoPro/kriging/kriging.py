import numpy as np
import pandas as pd
import warnings
import sys
import folium
from matplotlib.colors import to_hex


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

        self.__load_mesh()
        self.__load_smoothed_contour()

    def __load_mesh(self):

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

        self.mesh_df = df

    def __load_smoothed_contour(self):

        df = pd.read_excel(
            self.EPro.params['data_root_dir'] + self.EPro.params['filename_smoothed_contour'],
            sheet_name='Smoothing_EasyKrig')  # TODO: make the sheet name an input

        # obtaining those columns that are required
        df = df[['Latitude', 'Longitude']].copy()

        # set data types of dataframe
        df = df.astype({'Latitude': float,
                        'Longitude': float})

        self.smoothed_contour_df = df

    def plot_mesh(self, geo_df, location, cmap_column, cmap):
        """

        Parameters
        ----------
        geo_df : GeoPandas Dataframe
            Contains the data to be plotted with a geometry column
        location : list
            A length two list representing the location of folium.Map
        cmap_column : str
            Column of geo_df that the colormap should correspond to
        cmap : matplotlib.cm
            Matplotlib colormap
        # TODO: make kwargs for folium.Map and folium.CircleMarker

        Returns
        -------
        folium.Map
        """

        uniq_vals = geo_df[cmap_column].unique()
        hex_color_options = {rgb[0]: to_hex(rgb[1]) for rgb in zip(uniq_vals, cmap(uniq_vals))}

        fmap = folium.Map(location=location, zoom_start=4, tiles='CartoDB positron')

        for index, row in geo_df.iterrows():

            coordinates = (row.geometry.xy[1][0], row.geometry.xy[0][0])

            # Place the markers with the popup labels and data
            fmap.add_child(folium.CircleMarker(location=coordinates,
                                               radius=1,
                                               color=hex_color_options[row[cmap_column]]))

        return fmap
