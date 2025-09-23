<<<<<<< HEAD
VARIOGRAM_MODEL_PARAMETER_MAP = {
    "Bessel-Exponential": ["bessel", "exponential"],
    "Bessel-Gaussian": ["bessel", "gaussian"],
    "Cosine-Exponential": ["cosine", "exponential"],
    "Cosine-Gaussian": ["cosine", "gaussian"],
    "Cubic": "cubic",
    "Exponential": "exponential",
    "Exponential-Linear": ["exponential", "linear"],
    "Gaussian": "gaussian",
    "Gaussian-Linear": ["gaussian", "linear"],
    "J-Bessel": "jbessel",
    "K-Bessel": "kbessel",
    "Linear": "linear",
    "Matérn": "matern",
    "Nugget": "nugget",
    "Pentaspherical": "pentaspherical",
    "Power law": "power",
    "Rational quadratic": "quadratic",
    "Cardinal sine": "sinc",
    "Spherical": "spherical",
=======
import abc
import numpy as np
import pandas as pd
from typing import Callable, Union, Dict, List, Literal, Optional, Any
from functools import reduce

# Import the existing acoustics functions
from ..acoustics import ts_length_regression, to_linear, to_dB, impute_missing_sigma_bs
from echopop.nwfsc_feat import utils
from typing import Optional, Tuple

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy import interpolate

#######

########
from echopop.survey import Survey
from echopop.biology import age1_metric_proportions, impute_kriged_values, reallocate_kriged_age1
from echopop.spatial.transect import correct_transect_intervals
from echopop.analysis import process_transect_data

survey = Survey(init_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config_2019.yml",
                survey_year_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_single_biodata_config.yml")
# survey = Survey(init_config_path = "C:/Users/Brandyn Lucca/Documents/GitHub/echopop/config_files/initialization_config_2019.yml",
#                 survey_year_config_path = "C:/Users/Brandyn Lucca/Documents/GitHub/echopop/config_files/survey_year_2019_single_biodata_config.yml")
survey.load_acoustic_data(ingest_exports="echoview")
survey.load_survey_data()
survey.transect_analysis()
survey.stratified_analysis()
survey.kriging_analysis()

import inspect
import warnings
from typing import Any, Callable, Dict, Optional, Tuple, Union

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj
import verde as vd
from cartopy.mpl.geoaxes import GeoAxes
from matplotlib.axes import Axes
from matplotlib.colors import LogNorm, Normalize, SymLogNorm
from matplotlib.ticker import FixedLocator
import cartopy.feature as cfeature
from cartopy import crs as ccrs
import shapely
import cartopy.feature as cfeature
from geopy.distance import distance
import pandas as pd
import numpy as np
import geopandas as gpd
import shapely
import shapely.geometry as sg
import shapely.geometry as sg
from bokeh.themes.theme import Theme
import geoviews.tile_sources as gvts
import inspect
from echopop.nwfsc_feat import utils
from matplotlib.ticker import FixedLocator
import matplotlib as mpl
from pathlib import Path
plt.close('all')
mpl.use("tkagg")

save_directory = DATA_ROOT / "reports_test"
filename = "EchoPro_kriged_output_0.xlsx"
save_directory
kriged_data = df_kriged_results.copy()
sigma_bs_data = invert_hake.sigma_bs_strata.copy()
kriged_stratum = "geostratum_ks"
sigma_bs_stratum = "stratum_ks"
bio_data = dict_df_bio["specimen"].dropna(subset=["age", "length", "weight"]).copy()
filename = "aged_len_haul_counts_table.xlsx"
sheetnames = {"male": "Sheet1", "female": "Sheet2", "all": "Sheet3"}

"title": "Aged Length-Haul Counts ({SEX})",
"filename": "aged_len_haul_counts_table.xlsx",

Path(save_directory)
self.verbose = True
self = FEATReports(save_directory)
class FEATReports:

    def __init__(self, save_directory: Union[str, Path], verbose: bool = True):

        # Ensure Path-typing
        save_directory = Path(save_directory)

        # Validate existence -- create path if missing
        try:
            if not save_directory.exists():
                save_directory.mkdir(parents=True, exist_ok=True)
            # ---- Store attribute
            self.save_directory = save_directory
        except Exception as e:
            raise FileNotFoundError(
                f"The save directory '{save_directory.as_posix()}' could not be accessed or "
                f"created due to the following error(s):\n{e}"
            )

    def aged_length_haul_counts_report(
        bio_data: pd.DataFrame,
        filename: str,
        sheetnames: Dict[str, str],
    ) -> None:

        # Create copy
        bio_data = bio_data.copy()

        # Filter out unsexed -- i.e. keep only sexed fish
        bio_data_cleaned = bio_data.loc[bio_data["sex"].isin(["male", "female"])]

        # Create a pivot table
        bio_pvt = bio_data_cleaned.pivot_table(
            columns=["sex", "length_bin"],
            index=["haul_num"],
            values="length",
            aggfunc="count",
            observed=False
        )

        # Add a category for "all"
        bio_pvt_all = bio_pvt["male"] + bio_pvt["female"]
        # ---- Convert into a MultiIndex
        bio_pvt_all.columns = pd.MultiIndex.from_product([["all"], bio_pvt_all.columns])
        # ---- Concatenate
        bio_pvt_full = pd.concat([bio_pvt, bio_pvt_all], axis=1)        

    def kriged_mesh_results_report(
        self, 
        filename: str,
        sheetname: str,
        kriged_data: pd.DataFrame,
        kriged_stratum: str,
        kriged_variable: str,
        sigma_bs_data: pd.DataFrame,
        sigma_bs_stratum: str,
    ) -> None:

        # Create DataFrame copies
        kriged_data = kriged_data.copy()
        sigma_bs_data = sigma_bs_data.copy()

        # Update the filepath
        filepath = self.save_directory / filename

        # Check for kriged variable
        if kriged_variable not in kriged_data.columns:
            raise KeyError(
                f"Column '{kriged_variable}' not found in the kriged mesh DataFrame."
            )
        
        # Validate kriged mesh stratum name
        if (
            kriged_stratum not in kriged_data.columns and 
            kriged_stratum not in kriged_data.index.names
        ):
            raise KeyError(
                f"Column '{kriged_stratum}' not found in the kriged mesh DataFrame."
            )
        else:
            # ---- Reset index
            kriged_data.reset_index(inplace=True)
            # ---- Rename to generic
            kriged_data.rename(columns={kriged_stratum: "stratum"}, inplace=True)
            # ---- Set new index
            kriged_data.set_index("stratum", inplace=True)

        # Validate sigma_bs stratum name
        if (
            sigma_bs_stratum not in sigma_bs_data.columns and 
            sigma_bs_stratum not in sigma_bs_data.index.names
        ):
            raise KeyError(
                f"Column '{sigma_bs_stratum}' not found in the sigma_bs DataFrame."
            )
        else:
            # ---- Reset index
            sigma_bs_data.reset_index(inplace=True)
            # ---- Rename to generic
            sigma_bs_data.rename(columns={sigma_bs_stratum: "stratum"}, inplace=True)
            # ---- Set new index
            sigma_bs_data.set_index("stratum", inplace=True)
            # ---- Reindex to align with kriged mesh DataFrame
            sigma_bs_data = sigma_bs_data.reindex(kriged_data.index)

        # Add column to kriged data
        kriged_data["sig_b"] = 4. * np.pi * sigma_bs_data["sigma_bs"]
        # ---- Reset its index
        kriged_data.reset_index(inplace=True)

        # Add `kriged_sd` to the kriged data
        kriged_data["krig_SD"] = kriged_data[kriged_variable] * kriged_data["cell_cv"]

        # Subset the columns that are required for the report
        output_df = kriged_data.filter([
            "longitude", "latitude", "stratum", "nasc", "abundance_male", "abundance_female", 
            "abundance", "biomass_male", "biomass_female", "biomass", "sig_b", "cell_cv", 
            "krig_SD"
        ])

        # Rename the columns to match required output
        output_df.rename(
            columns={
                "latitude": "Lat",
                "longitude": "Lon",
                "nasc": "NASC",
                "abundance_male": "ntk_male",
                "abundance_female": "ntk_female",
                "abundance": "ntk_total",
                "biomass_male": "wgt_male",
                "biomass_female": "wgt_female",
                "biomass": "wgt_total",
                "cell_cv": "krig_CV",
            },
            inplace=True
        )

        # Save the *.xlsx sheet
        output_df.to_excel(filepath, sheet_name=sheetname, index=None)

        # Verbosity
        if self.verbose:
            print(
                f"Kriged mesh report saved to '{filepath.as_posix()}'."
            )

bio_pvt_full.T.loc["male", :]
# Subset the data into the specific sex
# haul_sex = dataframe.loc[sex, :]
haul_sex = bio_pvt_full.loc[:, "male"]
haul_sex.columns = pd.Series(haul_sex.columns).apply(lambda x: x.mid).values.to_numpy() 



# Stack the table
haul_stk = haul_sex.stack(future_stack=True).reset_index(name="count")

# Convert the length bins into number lengths
haul_stk.loc[:, "length"] = haul_stk["length_bin"].apply(lambda x: x.mid).astype(float)

# Repivot the table with margin subtotals
haul_agg = haul_stk.pivot_table(
    index=["length"],
    columns=["haul_num"],
    values="count",
    margins=True,
    margins_name="Subtotal",
    aggfunc="sum",
)

(save_directory / "AHHH" / "teehee.xlsx").resolve()

##########################
mesh_data = df_kriged_results.copy()
geostratum_df = df_dict_geostrata["inpfc"].copy()
geostrata_df=geostratum_df.copy()
stratify_by = ["geostratum_inpfc"]
variable = "biomass"
mesh_transects_per_latitude = 5
num_replicates = 10000
strata_transect_proportion = 0.75
##########################

## ~~~
transects_per_latitude = mesh_transects_per_latitude
## ~~~

model_params = {
    "transects_per_latitude": 5,
    "strata_transect_proportion": 0.75,
    "num_replicates": 100
>>>>>>> afc648a (Initial FEAT report things)
}

DEFAULT_VARIOGRAM_PARAMETERS = {
    "correlation_range": {
        "name": r"Correlation range ($a$)" ,
        "min": 1e-10, "value": 1.0, "max": 99999, "vary": False, "step": 1.0
    },
    "decay_power": {
        "name": r"Decay power ($\omega$)",
        "min": 1e-10, "value": 1.0, "max": 2.0, "vary": False, "step": 0.1
    },
    "enhance_semivariance": {
        "name": "Enhance semivariance", 
        "value": True
    },
    "hole_effect_range": {
        "name": r"Hole effect range ($a_h$)",
        "min": 0.0, "value": 0.0, "max": 99999, "vary": False, "step": 1.0
    },
    "sill": {
        "name": r"Sill ($C$)",
        "min": 1e-10, "value": 1.0, "max": 99999, "vary": False, "step": 1.0
    },
    "nugget": {
        "name": r"Nugget ($C_0$)",
        "min": 0.0, "value": 0.0, "max": 99999, "vary": False, "step": 1.0
    },
    "smoothness_parameter": {
        "name": r"Matérn shape parameter ($\nu$)",
        "min": 0.0, "value": 0.5, "max": 10.0, "vary": False, "step": 0.1
    },
    "shape_parameter": {
        "name": r"Scale ($\beta$)",
        "min": 1e-10, "value": 1.0, "max": 100.0, "vary": False, "step": 1.0
    },
    "power_exponent": {
        "name": r"Power ($\omega$)",
        "min": 1e-10, "value": 1.0, "max": 2.0 - 1e-10, "vary": False, "step": 0.1
    }
}

VARIABLE_MAP = {
    "Number density (animals nmi^-2)": "number_density", 
    "Biomass density (kg nmi^-2)": "biomass_density", 
    "NASC (m^2 nmi^-2)": "nasc",
}


class CompleteVariogramGUI(param.Parameterized):
    """Complete interactive GUI for variogram analysis with base parameter controls"""
    
    def __init__(self, data=None, **params):
        super().__init__(**params)
        self.data = data
        self.vgm = None  # Will be created based on GUI parameters
        self.empirical_results = {}
        self.model_results = {}
        self.optimization_results = {}

import holoviews as hv
hv.extension('bokeh')
from bokeh.models import HoverTool
from IPython.display import display_html
self = CompleteVariogramGUI(data)

self.vgm = vgm


####
lags, gamma, lag_counts = (self.vgm.lags, self.vgm.gamma, self.vgm.lag_counts)

data = pd.DataFrame({
    "distance": lags, 
    "semivariance": gamma, 
    "lag_count": lag_counts, 
    "lag_number": range(1, len(lags) + 1)
}) 

# Scale point sizes based on lag counts
if len(lag_counts) > 1 and lag_counts.max() > lag_counts.min():
    size_scale = lambda x: ((x - x.min()) / (x.max() - x.min()) + 1) * 15
    data['point_size'] = size_scale(lag_counts)
else:
    data['point_size'] = 15

scatter = hv.Scatter(
    data, 
    kdims=['distance'], 
    vdims=['semivariance', 'lag_count', 'lag_number', 'point_size']
).opts(
    size='point_size',
    color='blue',
    line_color='black',
    line_width=1,
    width=700,
    height=400,
    xlabel='Lag distance [h]',
    ylabel='Semivariance [γ]',
    title='Empirical Variogram',
    show_grid=True,
    # Explicitly specify tools to prevent duplication and include reset
    tools=['hover', 'pan', 'wheel_zoom', 'box_zoom', 'reset'],
    active_tools=['pan']
)

# Custom hover tool configuration
hover = HoverTool(
    tooltips=[
        ('Lag #', '@lag_number'),
        ('Distance', '@distance{0.000}'),
        ('Semivariance', '@semivariance{0.000}'),
        ('Count', '@lag_count')
    ]
)

hv.output(scatter)
