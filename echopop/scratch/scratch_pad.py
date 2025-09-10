import abc
import numpy as np
import pandas as pd
from typing import Callable, Union, Dict, List, Optional, Any
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

self = survey
stratum = "ks"
exclude_age1 = True
species_id = 22500
input_dict, analysis_dict, configuration_dict, settings_dict = self.input, self.analysis["transect"], self.config, self.analysis["settings"]
analysis_dict, results_dict, spatial_dict, settings_dict = (
    self.analysis,
            self.results,
            self.input["spatial"],
            self.analysis["settings"]["stratified"]
)
transect_data, transect_summary, strata_summary, settings_dict = (
    transect_data, transect_summary, strata_summary, settings_dict
)
population_dict = survey_dict
bootstrap_samples, population_statistic, ci_percentile, boot_ci_method, boot_ci_method_alt, adjust_bias, estimator_name = (
        bootstrap_dict["density"]["stratum"],
        population_dict["density"]["stratum"],
        settings_dict["bootstrap_ci"],
        settings_dict["bootstrap_ci_method"],
        settings_dict["bootstrap_ci_method_alt"],
        settings_dict["bootstrap_adjust_bias"],
        f"STRATUM {var_name.upper()} DENSITY",
)

rhom_all = (strata_transect_areas.to_numpy() * mean_arr_np).sum(axis=1) / strata_transect_areas.to_numpy().sum()





            / strata_transect_areas.sum()).mean()

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
}

mp = model_params

self = JollyHampton(model_params)
data_df = self.create_virtual_transects(mesh_data, geostratum_df, stratify_by=["geostratum_inpfc"], variable="biomass")

# Generate bins
latitude_bins = np.concatenate([[-90.0], geostratum_df["northlimit_latitude"], [90.0]])

data_df = mesh_data.copy()

# Partition the mesh values into virtual transects
mesh_data.loc[:, "latitude"] = (
    np.round(mesh_data.loc[:, "latitude"] * transects_per_latitude + 0.5)
    / transects_per_latitude
)

# Get the unique latitude values
unique_latitude_transect_key = pd.DataFrame(
    {
        "latitude": np.unique(mesh_data["latitude"]),
        "transect_num": np.arange(0, len(np.unique(mesh_data["latitude"])), 1) + 1,
    }
).set_index("latitude")

# Temporarily set `mesh_data` index
mesh_data.set_index("latitude", inplace=True)

# Append the transect numbers
mesh_data["transect_num"] = unique_latitude_transect_key

# Reset the index
unique_latitude_transect_key.reset_index(inplace=True)

# Reset
mesh_data.reset_index(inplace=True)

# Create equivalent transect dataframe needed for the stratified summary analysis
virtual_transect_data = mesh_data[
    stratify_by + ["transect_num", "longitude", "latitude", "area", variable]
].rename(columns={"area": "area_interval"})

# Cut the virtual latitudes into their unique strata
mesh_data.loc[:, stratify_by] = pd.cut(
    mesh_data["latitude"],
    latitude_bins,
    labels=list(geostratum_df["stratum_num"]) + [1],
    ordered=False,
)


####
i = 1
import geopy.distance

df = virtual_transect_data.sort_values(["transect_num"]).set_index(["transect_num"])#.loc[i]
df_red = df.groupby(level=0).apply(
    lambda x: geopy.distance.distance(
        (x.latitude.min(), x.longitude.min()), (x.latitude.max(), x.longitude.max())
    ).nm
).to_frame("distance")

df_red.loc[:, "area"] = np.where(
    df_red.index.isin([virtual_transect_data.transect_num.min(), 
                       virtual_transect_data.transect_num.max()]),
    df_red["distance"] * np.diff(unique_latitude_transect_key["latitude"]).mean() * 60,
    df_red["distance"] * np.diff(unique_latitude_transect_key["latitude"]).mean() * 60 / 2
)

# Flip the latitude key
unique_latitude_transect_key.set_index("transect_num", inplace=True)

# Set the latitude
df_red.loc[:, "latitude"] = unique_latitude_transect_key["latitude"]

# 
df_red = load_data.join_geostrata_by_latitude(df_red, 
                                              geostratum_df, 
                                              stratum_name=stratify_by[0])

# Sum variable
df_red[variable] = virtual_transect_data.groupby(["transect_num"])[variable].sum()

############

self = JollyHampton(mp)

# Enumerate the number of transects per stratum
strata_transect_counts = df_red.groupby(stratify_by, observed=False)["distance"].count()

# Calculate the total transect area per stratum
strata_transect_areas = df_red.groupby(stratify_by, observed=False)["area"].sum()

# Index by stratum
stratum_index = df_red[stratify_by].reset_index().set_index(stratify_by)

# Calculate the number of transects per stratum
num_transects_to_sample = np.round(
    strata_transect_counts * strata_transect_proportion
).astype(int)

# Offset term used for later variance calculation
sample_offset = np.where(num_transects_to_sample == 1, 0, 1)

# Calculate effective sample size/degrees of freedom for variance calculation
sample_dof = num_transects_to_sample * (num_transects_to_sample - sample_offset)

import awkward as awk

rng = np.random.default_rng()  # Use the new Generator API
np.random.default_rng(None).choice([1, 2, 3, 4], 1)
n_replicates = 10000

# For each stratum, generate replicate samples of transect indices
transect_samples = [
    np.sort(
        np.array([
            rng.choice(stratum_index.loc[j, "transect_num"].values, 
                       size=num_transects_to_sample.loc[j], replace=False)
            for _ in range(num_replicates)
        ])
    )
    for j in num_transects_to_sample.index
]

# Map integer positions to actual index values for each stratum
# transect_samples = [
#     np.sort(df_red.index.values[sample])  # sample is an array of integer positions
#     for sample in transect_samples_idx
# ]

samples_ak = awk.Array(transect_samples)  # shape: (n_strata, num_replicates, variable n_transects)

# Use advanced indexing to get the sampled values
# ---- Distances
sampled_distances = awk.Array([
    [df_red.loc[replicate, "distance"].to_numpy() for replicate in stratum]
    for stratum in transect_samples
])
# ---- Areas
sampled_areas = awk.Array([
    [df_red.loc[replicate, "area"].to_numpy() for replicate in stratum]
    for stratum in transect_samples
])
# ---- Variable
sampled_values = awk.Array([
    [df_red.loc[replicate, variable].to_numpy() for replicate in stratum]
    for stratum in transect_samples
])

# Now, all calculations are vectorized over awkward arrays:
length_arr   = awk.sum(sampled_distances, axis=-1)  # shape: (n_strata, n_replicates)
area_arr     = awk.sum(sampled_areas, axis=-1)
total_arr    = awk.sum(sampled_values, axis=-1)

#
stratified_weights = sampled_distances / awk.mean(sampled_distances, axis=-1, keepdims=True)

#
biology_adjusted = sampled_values / sampled_distances

#
mean_arr     = awk.sum(sampled_values * sampled_distances, axis=-1) / length_arr
biology_adjusted = sampled_values / sampled_distances

stratified_weights = sampled_distances / awk.mean(sampled_distances, axis=-1, keepdims=True)

squared_deviation = (biology_adjusted - mean_arr[..., None]) ** 2
squared_deviation_wgt = awk.sum(stratified_weights**2 * squared_deviation, axis=-1)

# sample_dof: a 1D array or Series, one per stratum
variance_arr = squared_deviation_wgt / awk.Array(sample_dof.values)[:, None]

# Convert and transpose all arrays in one step
length_arr_np, mean_arr_np, area_arr_np, total_arr_np, variance_arr_np = [
    awk.to_numpy(arr).T
    for arr in [length_arr, mean_arr, area_arr, total_arr, variance_arr]
]

mean_arr_np.sum(axis=1)
DEG = (mean_arr_np * strata_transect_areas.to_numpy())
DEG * 1e-9

transect_samples[0][0] = np.array([1, 3, 5, 6, 7])
transect_samples[1][0] = np.array([8,    10,    12,    13,    14,    15,    16,    17,    20,    21,    22,    23,    24,    25,    26,    27,    29]) 
transect_samples[2][0] = np.array([31, 32, 33, 35, 36, 37, 38, 39, 41])                                                                               
transect_samples[3][0] = np.array([43, 45, 46, 47, 48, 49, 50, 52, 53, 54, 55]) 
transect_samples[4][0] = np.array([57, 58, 59, 61, 62, 63, 64, 66, 67, 69, 70])
transect_samples[5][0] = np.array([71, 72, 73, 74, 75, 78, 79, 80, 81, 83, 84, 86, 87, 88, 90, 91, 93, 95, 96, 97, 98, 99, 101, 102])

mean_arr_np.min(axis=0) * 1e-7; mean_arr_np.mean(axis=0) * 1e-7; mean_arr_np.max(axis=0) * 1e-7
total_arr_np.min(axis=0) * 1e-8; total_arr_np.mean(axis=0) * 1e-8; total_arr_np.max(axis=0) * 1e-8
variance_arr_np.min(axis=0) * 1e-14; variance_arr_np.mean(axis=0) * 1e-14; variance_arr_np.max(axis=0) * 1e-14

# ---- By stratum (density)
unweighted_stratum_density = mean_arr_np / length_arr_np

# ---- By stratum (total)
unweighted_stratum_total = unweighted_stratum_density * strata_transect_areas.to_numpy()

# ---- By survey (total)
unweighted_survey_total = unweighted_stratum_total.sum(axis=1)

# ---- By survey (density)
unweighted_survey_density = unweighted_survey_total / strata_transect_areas.sum()

# ---- Proportional stratum distributions
unweighted_stratum_proportions = total_arr_np / total_arr_np.sum(axis=1, keepdims=True)

# ---- Transect-length weighted coefficient of variation (CV)
weighted_variance = (variance_arr_np * strata_transect_areas.to_numpy()**2).sum(axis=1)
weighted_stdev = np.sqrt(weighted_variance)
weighted_mean = (mean_arr_np * strata_transect_areas.to_numpy()).sum(axis=1)
bootstrap_cv = weighted_stdev / weighted_mean

weighted_mean.mean() * 1e-11

(unweighted_stratum_density.T * length_arr_np.sum(axis=1)).sum(axis=1) * 1e-9

    # Create bootstrapped results dictionary
    bootstrap_dict = {
        "density": {"stratum": unweighted_stratum_density, "survey": unweighted_survey_density},
        "total": {"stratum": unweighted_stratum_total, "survey": unweighted_survey_total},
        "proportions": unweighted_stratum_proportions,
        "cv": bootstrap_cv,
    }

    # Estimate the confidence intervals (CIs) and biases for the survey data using the bootstrapped
    # results
    bootstrapped_cis = bootstrap_confidence_intervals(bootstrap_dict, survey_dict, settings_dict)

    # Output the related summary statistics
    # ---- Save the output resampled distributions
    resampled_distributions = pd.DataFrame(
        {
            "realization": np.arange(1, transect_replicates + 1),
            "unweighted_survey_density": unweighted_survey_density,
            "unweighted_survey_total": unweighted_survey_total,
            "weighted_survey_total": weighted_mean,
            "weighted_survey_variance": weighted_variance,
            "survey_cv": bootstrap_cv,
        }
    )
    # ---- Save the stratified results
    stratified_results = {
        "variable": settings_dict["variable"],
        "ci_percentile": 0.95,
        "num_transects": strata_copy["transect_count"].sum(),
        "stratum_area": area_array,
        "total_area": total_area,
        "estimate": {
            "strata": {
                "density": stratum_density_means,
                "total": stratum_total,
                "proportion": stratum_proportions,
            },
            "survey": {
                "density": survey_density_mean,
                "total": survey_total,
                "cv": bootstrap_cv.mean(),
            },
        },
        "ci": bootstrapped_cis["ci"],
        "bias": bootstrapped_cis["bias"],
    }
    # ---- Return outputs
    return resampled_distributions, stratified_results
    

        # # Create copy
        # data_df = data_df.copy()

        # # Partition the dataset into virtual transects based on latitude
        # data_df.loc[:, "latitude"] = (
        #     np.round(data_df.loc[:, "latitude"] * mp["transects_per_latitude"] + 0.5) /
        #     mp["transects_per_latitude"]
        # )

        # # Create unique key pairs for latitude and transect
        # unique_latitude_transect_key = pd.DataFrame(
        #     {
        #         "latitude": np.unique(data_df["latitude"]),
        #         "transect_num": np.arange(0, len(np.unique(data_df["latitude"])), 1) + 1,
        #     }
        # ).set_index("latitude")

        # # Temporarily set `mesh_data` index
        # data_df.set_index("latitude", inplace=True)

        # # Append the transect numbers
        # data_df["transect_num"] = unique_latitude_transect_key

        # # Reset the key index
        # unique_latitude_transect_key.reset_index(inplace=True)

        # # Reset the dataset index
        # data_df.reset_index(inplace=True)

        # # Create equivalent transect dataframe needed for the stratified summary analysis
        # virtual_df = data_df[
        #     stratify_by + ["transect_num", "longitude", "latitude", "area", variable]
        # ].rename(
        #     columns={"area": "area_interval"}
        # ).sort_values(["transect_num"]).set_index(["transect_num"])

        # # Compute the transect distances
        # virtual_transect_data = virtual_df.groupby(level=0).apply(
        #     lambda x: geopy.distance.distance(
        #         (x.latitude.min(), x.longitude.min()), (x.latitude.max(), x.longitude.max())
        #     ).nm
        # ).to_frame("distance")

        # # Set the areas
        # virtual_transect_data.loc[:, "area"] = np.where(
        #     virtual_transect_data.index.isin([
        #         data_df.transect_num.min(), data_df.transect_num.max()
        #     ]),
        #     virtual_transect_data["distance"] * np.diff(
        #         unique_latitude_transect_key["latitude"]
        #     ).mean() * 60,
        #     virtual_transect_data["distance"] * np.diff(
        #         unique_latitude_transect_key["latitude"]
        #     ).mean() * 60 / 2
        # )

        # # Flip the latitude key
        # unique_latitude_transect_key.set_index("transect_num", inplace=True)

        # # Set the latitude
        # virtual_transect_data.loc[:, "latitude"] = unique_latitude_transect_key["latitude"]

        # # Stratify the virtual transects
        # virtual_transect_data = load_data.join_geostrata_by_latitude(
        #     virtual_transect_data, 
        #     geostrata_df, 
        #     stratum_name=stratify_by[0]
        # )

        # # Sum the biological variable
        # virtual_transect_data[variable] = virtual_df.groupby(["transect_num"])[variable].sum()

        # # Return the DataFrame
        # return virtual_transect_data

    def _prepare_bootstrap_arrays(
        self,
        data_df: pd.DataFrame,
        stratum_index: pd.DataFrame,
        num_transects_to_sample: pd.DataFrame,
        variable: str,
    ) -> Tuple[awk.Array, awk.Array, awk.Array]:

        # For each stratum, generate replicate samples of transect indices
        transect_samples = [
            np.sort(
                np.array([
                    self.rng.choice(stratum_index.loc[j, "transect_num"].values, 
                            size=num_transects_to_sample.loc[j], replace=False)
                    for _ in range(self.model_params["num_replicates"])
                ])
            )
            for j in num_transects_to_sample.index
        ]

        # Use advanced indexing to get the sampled values
        # ---- Distances
        sampled_distances = awk.Array([
            [data_df.loc[replicate, "distance"].to_numpy() for replicate in stratum]
            for stratum in transect_samples
        ])
        # ---- Areas
        sampled_areas = awk.Array([
            [data_df.loc[replicate, "area"].to_numpy() for replicate in stratum]
            for stratum in transect_samples
        ])
        # ---- Variable
        sampled_values = awk.Array([
            [data_df.loc[replicate, variable].to_numpy() for replicate in stratum]
            for stratum in transect_samples
        ])

        # Return a tuple of the uneven arrays
        return sampled_distances, sampled_areas, sampled_values     

    @staticmethod
    def _get_variance(
        values: awk.Array,
        value_mean: awk.array,
        distance: awk.Array,
        weights: awk.Array,
        dof: np.array,
    ) -> awk.Array:

        # Adjust the values based on distance
        values_adjusted = values / distance

        # Calculate the squared deviation
        squared_deviation = (values_adjusted - value_mean[..., None]) ** 2

        # Sum the weighted squared deviations
        squared_deviation_wgt = awk.sum(weights**2 * squared_deviation, axis=-1)

        # Return the variance
        return squared_deviation_wgt / awk.Array(dof.values)[:, None]


    def stratified_bootstrap(
        self,
        data_df: pd.DataFrame,
        stratify_by: List[str],
        variable: str,
    ):

        # Get the model parameters
        mp = self.model_params

        # Enumerate the number of transects per stratum
        strata_transect_counts = data_df.groupby(stratify_by, observed=False)["distance"].count()

        # Calculate the total transect area per stratum
        strata_transect_areas = data_df.groupby(stratify_by, observed=False)["area"].sum()

        # Index by stratum
        stratum_index = data_df[stratify_by].reset_index().set_index(stratify_by)

        # Calculate the number of transects per stratum
        num_transects_to_sample = np.round(
            strata_transect_counts * mp["strata_transect_proportion"]
        ).astype(int)

        # Offset term used for later variance calculation
        sample_offset = np.where(num_transects_to_sample == 1, 0, 1)

        # Calculate effective sample size/degrees of freedom for variance calculation
        sample_dof = num_transects_to_sample * (num_transects_to_sample - sample_offset)

        # Compute the resampling arrays
        sampled_distances, sampled_areas, sampled_values = self._prepare_bootstrap_arrays(
            data_df, stratum_index, num_transects_to_sample, variable
        )

        # Compute the stratified weights
        stratified_weights = sampled_distances / awk.mean(sampled_distances, axis=-1, keepdims=True)

        # Compute the value mean
        self.value_mean_array = awk.sum(
            sampled_values * sampled_distances, axis=-1
        ) / awk.sum(sampled_distances, axis=-1)

        # Calculate the variance
        variance = self._get_variance(
            sampled_values, self.value_mean_array, sampled_distances, stratified_weights, sample_dof
        )

        # Convert and transpose all arrays in one step
        length_arr_np, mean_arr_np, area_arr_np, total_arr_np, variance_arr_np = [
            awk.to_numpy(arr).T
            for arr in [length_arr, mean_arr, area_arr, total_arr, variance_arr]
        ]

        # ---- By stratum (density)
        unweighted_stratum_density = mean_arr_np / length_arr_np

        # ---- By stratum (total)
        unweighted_stratum_total = unweighted_stratum_density * strata_transect_areas.to_numpy()

        # ---- By survey (total)
        unweighted_survey_total = unweighted_stratum_total.sum(axis=1)

        # ---- By survey (density)
        unweighted_survey_density = unweighted_survey_total / strata_transect_areas.sum()

        # ---- Proportional stratum distributions
        unweighted_stratum_proportions = total_arr_np / total_arr_np.sum(axis=1, keepdims=True)

        # ---- Transect-length weighted coefficient of variation (CV)
        weighted_variance = (variance_arr_np * strata_transect_areas.to_numpy()**2).sum(axis=1)


        # Sum the quantities over each stratum
        # ---- Transect lengths
        distance_totals = awk.sum(sampled_distances, axis=-1)
        # ---- Areas
        area_totals = awk.sum(sampled_areas, axis=-1)
        # ---- Biological variable
        variable_totals = awk.sum(sampled_values, axis=-1)


        

        # Output the related summary statistics
        # ---- Save the output resampled distributions
        resampled_distributions = pd.DataFrame(
            {
                "realization": np.arange(1, transect_replicates + 1),
                "unweighted_survey_density": unweighted_survey_density,
                "unweighted_survey_total": unweighted_survey_total,
                "weighted_survey_total": weighted_mean,
                "weighted_survey_variance": weighted_variance,
                "survey_cv": bootstrap_cv,
            }
        )

        # Sum the quantities over each stratum
        # ---- Transect lengths
        distance_totals = awk.sum(sampled_distances, axis=-1)
        # ---- Areas
        area_totals = awk.sum(sampled_areas, axis=-1)
        # ---- Biological variable
        variable_totals = awk.sum(sampled_values, axis=-1)

        
