from pathlib import Path
from echopop.workflows.nwfsc_feat import functions as feat
from echopop.graphics.utils import dataframe_to_geodataframe
from echopop.graphics.transect_map import get_transect_lines
import echopop.ingest as ingestion
import geopandas as gpd
import geoviews as gv
import holoviews as hv
import hvplot.pandas  # noqa: F401; import is required even though it is not directly called
import pandas as pd
import shapely as sg
from bokeh.models import HoverTool
from bokeh.themes import Theme

DATA_ROOT = Path("C:/Data/EchopopData/echopop_2011")
# NASC EXPORTS FILE(S)
NASC_EXPORTS_FILES = DATA_ROOT / "raw_nasc/"
# NASC EXPORTS SHEET
NASC_EXPORTS_SHEET = "Sheet1"

df_intervals, df_exports = ingestion.nasc.merge_echoview_nasc(
    nasc_path = NASC_EXPORTS_FILES,
    filename_transect_pattern = "T(\\d+(?:\\.\\d+)?)(?=-)",
    default_transect_spacing = 10.0,
    default_latitude_threshold = 60.0,
)

# TRANSECT REGION HAUL KEY NAME MAPPING
TRANSECT_REGION_FILE_RENAME = {
    "transect": "transect_num",
    "region id": "region_id",
    "assigned haul": "haul_num",
}

TRANSECT_REGION_HAUL_FILE = (
    DATA_ROOT / "Stratification/US&CAN_T_reg_haul_final.csv"
)

# LOAD
df_transect_region_haul_key = ingestion.nasc.read_transect_region_haul_key(
    filename=TRANSECT_REGION_HAUL_FILE,
    sheetname=None,
    rename_dict=TRANSECT_REGION_FILE_RENAME
)

df_nasc_old = ingestion.nasc.consolidate_echvoiew_nasc(
    df_merged=df_exports,
    interval_df=df_intervals,
    region_class_names=["Age-1 Hake", "Age-1 Hake Mix", "Hake", "Hake Mix"],
    impute_region_ids=True,
    transect_region_haul_key_df=df_transect_region_haul_key
)

# TRANSECT BOUNDARY FILE
TRANSECT_BOUNDARY_FILE = DATA_ROOT / "Kriging_files/Kriging_grid_files/Transect Bounds to 2011.xlsx"
# TRANSECT BOUNDARY SHEET
TRANSECT_BOUNDARY_SHEET = "1995-2011"
# SURVEY FILTER
SURVEY_FILTER = "survey == 201103"

df_nasc_new = feat.filter_transect_intervals(
    nasc_df=df_nasc_old, 
    transect_filter_df=TRANSECT_BOUNDARY_FILE,
    transect_filter_sheet=TRANSECT_BOUNDARY_SHEET,
    subset_filter=SURVEY_FILTER
)

diff = df_nasc_old.merge(df_nasc_new, how='outer', indicator=True)
rows_not_in_new = diff[diff['_merge'] == 'left_only'].drop(columns=['_merge'])
print(rows_not_in_new)

# Class-level flag to track renderer initialization
_is_bokeh_initialized = False

# Class-level flag to track renderer initialization
_is_theme_initialized = False

# Internal default JSON theme
_json_theme = {
    "attrs": {
        "Title": {
            "align": "left",
            "text_font_size": "20px",
            "text_color": "black",
        },
        "Axis": {
            "axis_label_text_font_style": "bold",
            "axis_label_text_color": "black",
            "axis_label_text_font_size": "18px",
            "major_label_text_font_size": "15px",
            "major_label_text_color": "black",
        },
        "ColorBar": {
            "title_text_font_style": "normal",
            "title_text_font_size": "18px",
            "title_text_color": "black",
            "major_label_text_font_size": "16px",
            "major_label_text_color": "black",
        },
        "Legend": {
            "title_text_font_style": "bold",
            "title_text_font_size": "16px",
            "title_text_color": "black",
        },
    },
}

hv.extension("bokeh")

hv.renderer("bokeh").theme = Theme(json=_json_theme)

# Format filtered data
projection: str = "epsg:4326"
filtered_gdf = dataframe_to_geodataframe(df_nasc_new, projection, ("longitude", "latitude"))

# Format values removed
removed_gdf = dataframe_to_geodataframe(rows_not_in_new, projection, ("longitude", "latitude"))

# Convert to LINESTRINGs
filtered_lines = get_transect_lines(filtered_gdf)

# LINE layer
filtered_lines_layer = filtered_lines.hvplot.paths(
    geo=True, line_width=2, label="On-effort", line_color="black"
).opts()

# POINT layer
off_effort = gv.Points(
    removed_gdf,
    label="Off-effort intervals"
).opts(
    color="red",
    line_color="black",
    size=4,
    line_width=0.1
)

# FULL overlay
plot = (gv.tile_sources.CartoLight * filtered_lines_layer * off_effort).opts(
    width=600, 
    height=400, 
    legend_position="bottom_left", 
    show_legend=True,
    xlabel="Longitude (\u00b0E)",
    ylabel="Latitude (\u00b0N)",
    title="Filtering out off-effort transect intervals",
    active_tools=["pan", "wheel_zoom"],
)

# SAVE
hv.save(plot, "../_static/offeffort_transect_filtering.html")