import numpy as np
import pandas as pd
import warnings
import sys


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

    def __load_mesh(self):

        mesh_df = pd.read_excel(self.EPro.params['data_root_dir'] + self.EPro.params['filename_grid_cell'],
                                sheet_name='krigedgrid2_5nm_forChu')  # TODO: make the sheet name an input

        self.mesh_df = mesh_df

    # import geopandas
    # import matplotlib.pyplot as plt
    #
    # gdf = geopandas.GeoDataFrame(df,
    #                              geometry=geopandas.points_from_xy(df.Longitude, df.Latitude))
    #
    # gdf
    #
    # import folium
    # map = folium.Map(location=[44.61, -125.66], zoom_start=4, tiles='CartoDB positron')
    #
    # # Create a geometry list from the GeoDataFrame
    # geo_df_list = [[point.xy[1][0], point.xy[0][0]] for point in gdf.geometry]
    #
    # for coordinates in geo_df_list:
    #     # Place the markers with the popup labels and data
    #     map.add_child(folium.CircleMarker(location=coordinates, radius=1, color='blue'))
    #
    # map





