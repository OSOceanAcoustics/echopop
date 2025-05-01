"""
Column name translation map for Echoview export data.

This dictionary maps raw input column names from external Echoview (*.csv) files to standardized
internal names used consistently throughout the Echopop codebase.    

Key Mappings
------------
- 'date_s' → 'ping_date' : date of first ping in the interval/domain (YYYYMMDD format)
- 'exclude_below_line_depth_mean' → 'max_depth' : mean maximum depth of interval/domain
- 'lat_s' → 'latitude' : latitude of first ping in the interval/domains (DD.ddddd format)
- 'lon_s' → 'longitude' : longitude of first ping in the interval/domains (DD.ddddd format)
- 'prc_nasc' → 'nasc' : region-integrated NASC (m^2 nmi^-2)
- 'time_s' → 'ping_time' : time of first ping in the interval/domain (HH:mm:ss.SSSS format)
- 'vl_end' → 'vessel_log_end' : vessel log distance of the last ping in the interval/domain 
(nmi)
- 'vl_start' → 'date_of_birth' : vessel log distance of the first ping in the interval/domain 
(nmi)
"""
ECHOVIEW_TO_ECHOPOP = {
    "date_s": "ping_date",
    "exclude_below_line_depth_mean": "max_depth",
    "lat_s": "latitude",
    "lon_s": "longitude",
    "prc_nasc": "nasc",
    "time_s": "ping_time",
    "vl_end": "vessel_log_end",
    "vl_start": "vessel_log_start",
}
    