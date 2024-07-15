import zarr
import xarray as xr
import shutil 
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

specimen_df = pd.DataFrame(
    {
        "haul_num": np.repeat([1,2,3], 4),
        "station": "specimen",
        "sex": np.tile(["male", "female"], 6),
        "length": np.array([11, 11, 11, 18, 21, 23, 13, 11, 19, 25, 18, 9]), 
        "weight": np.array([11, 14, 16, 18, 21, 23, 13, 11, 19, 25, 18, 9]) / 3.5,
    },
)

length_df = pd.DataFrame(
    {
        "haul_num": np.repeat([1,2,3], 4),
        "station": "length",
        "sex": np.tile(["male", "female"], 6),
        "length": np.array([16, 15, 19, 14, 9, 10, 18, 15, 16, 22, 17, 11]), 
        "length_count": np.array([103, 123, 257, 106, 52, 329, 131, 72, 101, 212, 93, 81]),
    },
)

catch_df = pd.DataFrame(
    {
        "haul_num": np.array([1, 2, 3]),
        "weight": np.array([503.12, 684.32, 978.54])
    }
)

TS_SLOPE = 20.0
TS_INTERCEPT = -68.0

####
# CONCATENATE FILE SOURCES
specimen_reframed = specimen_df.groupby(["haul_num", "station", "sex", "length"])["length"].value_counts().to_frame("length_count").reset_index()
specimen_reframed
# MELD
all_lengths = pd.concat([length_df, specimen_reframed])
# COMBINE 
comb_lengths = all_lengths.groupby(["haul_num", "sex", "length"])["length_count"].sum().to_frame("length_count").reset_index()


# CONVERT TO TS
comb_lengths["ts"] = TS_SLOPE * np.log10(comb_lengths["length"]) + TS_INTERCEPT
# TO SIGMA_BS
comb_lengths["sigma_bs"] = 10 ** (comb_lengths["ts"] / 10)
# WEIGHTED MEAN SIGMA_BS
sigma_mean = np.average(comb_lengths["sigma_bs"], weights=comb_lengths["length_count"])

### 
# INTEGRATE NASC
path2file = "C:/Users/15052/Downloads/win_1720457505_1720460000_NASC.zarr"

Path(path2file).exists()
xds = xr.open_dataset(path2file, engine="zarr")
xds
xdf = xds.to_dataframe().reset_index()
xdf["NASC"] = xdf["NASC"].fillna(0.0)
# convert frequency
xdf["frequency_nominal"] = (xdf["frequency_nominal"] * 1e-3).astype(int)
# filter
xdf_38 = xdf[xdf["frequency_nominal"] == nasc_frequency]

xdf_38.plot.scatter(x="distance", y="depth", c="NASC")
plt.show()

xdf_int = xdf_38.groupby(["distance", "longitude", "latitude"])["NASC"].sum().reset_index()

plt.scatter(xdf_int["longitude"], xdf_int["latitude"], c=xdf_int["NASC"])
plt.plot(xdf_int["longitude"], xdf_int["latitude"])
plt.show()

# CONVERT TO NUMBER DENSITY
xdf_int["number_density"] = xdf_int["NASC"] / (4.0 * np.pi * sigma_mean)


###################
from typing import Union
from pathlib import Path
import copy
import yaml

# from echopop.acoustics import ts_length_regression, to_dB, to_linear
# from echopop.live.core import DATA_STRUCTURE


### INIT CONFIG
initialization_config = "C:/Users/15052/Documents/GitHub/echopop/config_files/live_initialization_config.yml"

# Initialize `meta` attribute
meta = copy.deepcopy(LIVE_DATA_STRUCTURE["meta"])

# Loading the configuration settings and definitions that are used to
# initialize the Survey class object
config = yaml.safe_load(Path(initialization_config).read_text())

nasc_frequency = config["acoustics"]["nasc_frequency"]