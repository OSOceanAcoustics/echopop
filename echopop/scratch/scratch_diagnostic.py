# ==================================================================================================
# MESH REGION PLOTTING
# --------------------
from echopop.nwfsc_feat import FEAT
import geoviews as gv
import holoviews as hv
import pandas as pd
import numpy as np
import bokeh
import panel as pn
from IPython.display import clear_output, display, display_png

transect_data=df_nasc_no_age1_prt.copy()
mesh_data=df_mesh.copy()
latitude_resolution=1.25/60.
crop_function=FEAT.fun.transect_ends_crop
transect_mesh_region_function=FEAT.parameters.transect_mesh_region_2019

renderer = hv.plotting.mpl.MPLRenderer.instance(dpi=120)
hv.renderer("bokeh")

# Get the mesh region assignment [Tuple[mesh, transect]]
me, tr = crop_function(transect_data, mesh_data, latitude_resolution, transect_mesh_region_function)

pn.extension()
hv.extension("bokeh")
gv.extension("bokeh")


PLOT = gv.Points(tr, kdims=["longitude", "latitude"], vdims=["transect_num", "mesh_region"]).opts(color="mesh_region", show_legend=True)
renderer.get_plot(PLOT)
gv.render(PLOT)