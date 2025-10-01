# noqa: D100
"""
Echoview export data processing utilities.

This module contains constants and mappings for processing Echoview export data.
"""

# Column name translation map for Echoview export data.
#
# This dictionary maps raw input column names from external Echoview (*.csv) files to standardized
# internal names used consistently throughout the Echopop codebase.
#
# Key Mappings
# ------------
# - 'date_s' → 'ping_date' : date of first ping in the interval/domain (YYYYMMDD format)
# - 'exclude_below_line_depth_mean' → 'max_depth' : mean maximum depth of interval/domain
# - 'lat_s' → 'latitude' : latitude of first ping in the interval/domains (DD.ddddd format)
# - 'lon_s' → 'longitude' : longitude of first ping in the interval/domains (DD.ddddd format)
# - 'prc_nasc' → 'nasc' : region-integrated NASC (m^2 nmi^-2)
# - 'time_s' → 'ping_time' : time of first ping in the interval/domain (HH:mm:ss.SSSS format)
# - 'vl_end' → 'distance_e' : vessel log distance of the last ping in the interval/domain (nmi)
# - 'vl_start' → 'disance_s' : vessel log distance of the first ping in the interval/domain (nmi)
ECHOVIEW_TO_ECHOPOP = {
    "date_s": "ping_date",
    "exclude_below_line_depth_mean": "max_depth",
    "exclude_below_depth_mean": "max_depth",
    "lat_s": "latitude",
    "lon_s": "longitude",
    "prc_nasc": "nasc",
    "time_s": "ping_time",
    "vl_end": "distance_e",
    "vl_start": "distance_s",
}

# Valid Echoview export sorting columns
#
# This dictionary is used for sorting and reindexing the Echoview exports when ingested
#
# Key Mappings
# ------------
# - 'interval' : Interval number
# - 'layer' : Layer number
# - 'transect_num' : Transect number
ECHOVIEW_EXPORT_ROW_SORT = {
    "interval",
    "layer",
    "transect_num",
}


# Valid Echoview export filetypes
#
# Key Mappings
# ------------
# - 'analysis' : File containing metadata for the overall fileset when exported in the database
# format. This file may be deprecated to an optional file since it is otherwise not required in
# Echopop's current version
# - 'cells' : File that contains the region classifications, interval-layer mapped 'Sv_mean' and
# 'PRC_NASC' values
# - 'intervals' : File that contains georeferenced intervals, datetimes, and distance
# - 'layers' : File that contains the layer depths and dimensions
ECHOVIEW_DATABASE_EXPORT_FILESET = {"analysis", "cells", "intervals", "layers"}
