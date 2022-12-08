from matplotlib.colors import to_hex
from matplotlib import cm
import folium
import branca.colormap as bcm
import numpy as np
import geopandas as gpd
from typing import Union, Optional
from ..data_loader.kriging_mesh import KrigingMesh

# Default parameters for the folium map
folium_map_kwargs = {
    'location': [44.61, -125.66],
    'zoom_start': 4,
    'tiles': 'CartoDB positron'
}


def get_folium_map(map_kwargs: Optional[dict] = None) -> folium.Map:
    """
    Grabs the folium map.

    Parameters
    ----------
    map_kwargs : dict, optional
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
        map_kwargs = folium_map_kwargs

    return folium.Map(**map_kwargs)


def plot_points(gdf: gpd.GeoDataFrame,
                fobj: Union[folium.Map, folium.map.FeatureGroup] = None,
                color: str = 'hex', cmap_column: str = None,
                marker_kwargs: dict = {},
                map_kwargs: Optional[dict] = None) -> Union[folium.Map, folium.map.FeatureGroup]:
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
    map_kwargs : dict, optional
        Dictionary of kwargs for the folium.Map function.

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
        fobj = get_folium_map(map_kwargs)

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


def plot_layered_points(krig_mesh_obj: KrigingMesh) -> folium.Map:
    """
    This function constructs a layered Folium plot.
    The layers correspond to the full set of mesh
    points, the ``final_biomass_table`` points with
    color corresponding to the transect number, and
    the smoothed contour points (e.g. 200m isobath).

    Parameters
    ----------
    krig_mesh_obj : KrigingMesh
        An object specifying the Kriging mesh

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

    fmap = get_folium_map()

    # plot mesh points and add them to fmap
    folium_layer = folium.FeatureGroup(name='mesh')
    folium_layer = plot_points(krig_mesh_obj.mesh_gdf, folium_layer, color='gray')
    folium_layer.add_to(fmap)

    # plot the transect points and add them to fmap
    folium_layer = folium.FeatureGroup(name='transects')
    folium_layer = plot_points(krig_mesh_obj.survey.bio_calc.final_biomass_table, folium_layer,
                               cmap_column='transect_num', color='hex')
    folium_layer.add_to(fmap)

    # plot smoothed contour points and add them to fmap
    folium_layer = folium.FeatureGroup(name='smoothed contour')
    folium_layer = plot_points(krig_mesh_obj.smoothed_contour_gdf, folium_layer, color='blue')
    folium_layer.add_to(fmap)

    # add layer control to fmap
    folium.LayerControl().add_to(fmap)

    return fmap


# Visualization function for Kriging
def plot_kriging_results(krig_results_gdf: gpd.GeoDataFrame,
                         krig_field_name: str) -> folium.Map:
    """
    Constructs a Folium plot depicting the ``krig_field_name``
    values at each mesh point.

    Parameters
    ----------
    krig_results_gdf: gpd.GeoDataFrame
        Dataframe containing a geometry column that holds the
        mesh coordinates and the column ``krig_field_name``
    krig_field_name: str
        The name of the column in ``krig_results_gdf`` containing
        the Kriging values to plot at each mesh point

    Returns
    -------
    fmap : folium.Map
        A Folium plot with the ``krig_field_name`` values at each
        mesh point.
    """

    # create folium map
    fmap = folium.Map(location=[44.61, -125.66], zoom_start=4)

    # collect the appropriate data from the input Dataframe
    x_mesh = krig_results_gdf.geometry.x.values
    y_mesh = krig_results_gdf.geometry.y.values
    krig_val = krig_results_gdf[krig_field_name].values

    # create a colormap for the values
    colormap = bcm.LinearColormap(colors=['#3385ff', '#FF0000'],
                                  vmin=0.0, vmax=np.max(krig_val))

    # plot each mesh point and the corresponding krig_field_name value
    for i in range(len(krig_val)):

        # get the color of the marker we want to plot
        color_val = colormap(krig_val[i])

        # get coordinate of marker we are plotting
        coordinates = (y_mesh[i], x_mesh[i])

        # Place the markers with a color value
        fmap.add_child(folium.CircleMarker(location=coordinates,
                                           radius=1,
                                           color=color_val))

    fmap.add_child(colormap)  # adds color bar to map

    return fmap