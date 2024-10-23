"""
General analysis orchestration functions that bundle related functions and procedures
"""

import copy
import warnings

import numpy as np
import pandas as pd

from .acoustics import aggregate_sigma_bs, nasc_to_biomass
from .biology import (
    distribute_length_age,
    filter_species,
    fit_length_weight_relationship,
    fit_length_weights,
    impute_kriged_values,
    number_proportions,
    partition_transect_age,
    quantize_number_counts,
    quantize_weights,
    reallocate_kriged_age1,
    weight_proportions,
)
from .spatial.krige import kriging
from .spatial.mesh import crop_mesh, mesh_to_transects, stratify_mesh
from .spatial.projection import transform_geometry
from .spatial.transect import (
    edit_transect_columns,
    save_transect_coordinates,
    summarize_transect_strata,
    transect_spatial_features,
)
from .spatial.variogram import (
    empirical_variogram,
    initialize_initial_optimization_values,
    initialize_optimization_config,
    initialize_variogram_parameters,
    optimize_variogram,
)
from .statistics import stratified_transect_statistic
from .utils.validate_dict import (
    KrigingAnalysis,
    KrigingParameterInputs,
    MeshCrop,
    VariogramBase,
    VariogramEmpirical,
)


def process_transect_data(
    input_dict: dict, analysis_dict: dict, settings_dict: dict, configuration_dict: dict
) -> dict:
    """
    Process acoustic and biological data collected along each transect across an entire survey

    Parameters
    ----------
    input_dict: dict
        A dictionary containing the loaded survey data.
    analysis_dict: dict
        A dictionary containing processed biological and transect data.
    configuration_dict: dict
        Dictionary that contains all of the `Survey`-object configurations found within
        the `config` attribute.
    settings_dict: dict
        Dictionary that contains all of the analysis settings that detail specific algorithm
        arguments and user-defined inputs.
    """

    # Filter out non-target species
    length_data, specimen_data, catch_data = filter_species(
        [
            input_dict["biology"]["length_df"],
            input_dict["biology"]["specimen_df"],
            input_dict["biology"]["catch_df"],
        ],
        settings_dict["transect"]["species_id"],
    )
    # ---- For cases where all samples were aged (i.e. in `specimen_data` and absent from
    # ---- `length_data`), these hauls are removed from `catch_data`
    catch_data = catch_data[catch_data.haul_num.isin(length_data.haul_num)]

    # Save the transect coordinate information
    analysis_dict.update(
        {
            "coordinates": save_transect_coordinates(
                input_dict["acoustics"]["nasc_df"], settings_dict["transect"]
            )
        }
    )

    # Calculate mean sigma_bs per individual haul, KS stratum, and INPFC stratum
    analysis_dict["acoustics"]["sigma_bs"].update(
        aggregate_sigma_bs(length_data, specimen_data, configuration_dict, settings_dict)
    )

    # Fit length-weight regression required for biomass calculation
    analysis_dict["biology"]["weight"].update(
        fit_length_weight_relationship(
            specimen_data, input_dict["biology"]["distributions"]["length_bins_df"]
        )
    )

    # Count the number of specimens across age and length bins
    analysis_dict["biology"]["distributions"].update(
        quantize_number_counts(
            specimen_data, length_data, settings_dict["transect"]["stratum_name"]
        )
    )

    # Calculate the number proportions
    analysis_dict["biology"]["proportions"].update(
        {"number": number_proportions(analysis_dict["biology"]["distributions"])}
    )

    # Sum the weights of both aged and unaged specimens across length and weight bins
    # ---- Extract the length-weight fit
    length_weight_df = analysis_dict["biology"]["weight"]["length_weight_regression"][
        "weight_fitted_df"
    ]
    # ---- Quantize the weights
    analysis_dict["biology"]["distributions"].update(
        {
            "weight": quantize_weights(
                specimen_data,
                length_data,
                length_weight_df,
                settings_dict["transect"]["stratum_name"],
            )
        }
    )

    # Calculate the average weights among male, female, and all fish across strata
    analysis_dict["biology"]["weight"].update(
        {
            "weight_stratum_df": fit_length_weights(
                analysis_dict["biology"]["proportions"]["number"],
                analysis_dict["biology"]["weight"],
            )
        }
    )

    # Calculate the weight proportions
    analysis_dict["biology"]["proportions"].update(
        {
            "weight": weight_proportions(
                catch_data,
                analysis_dict["biology"]["proportions"]["number"],
                length_weight_df,
                analysis_dict["biology"]["distributions"]["weight"],
                settings_dict["transect"]["stratum_name"],
            )
        }
    )

    # Return the analysis dictionary
    return analysis_dict


def acoustics_to_biology(
    input_dict: dict, analysis_dict: dict, configuration_dict: dict, settings_dict: dict
) -> tuple[pd.DataFrame, dict]:
    """
    Convert acoustic transect data into population-level metrics such as abundance and biomass.

    Parameters
    ----------
    input_dict: dict
        A dictionary containing the loaded survey data.
    analysis_dict: dict
        A dictionary containing processed biological and transect data.
    configuration_dict: dict
        Dictionary that contains all of the `Survey`-object configurations found within
        the `config` attribute.
    settings_dict: dict
        Dictionary that contains all of the analysis settings that detail specific algorithm
        arguments and user-defined inputs.
    """

    # Convert NASC into number density (animals/nmi^2), biomass density (kg/nmi^2), abundance
    # (# animals), and biomass (kg) for all fish, sexed (male/female) fish, and unsexed fish
    strata_adult_proportions, nasc_to_biology = nasc_to_biomass(
        input_dict, analysis_dict, configuration_dict, settings_dict
    )

    # Distribute abundance and biomass estimates over length and age bins
    analysis_dict["biology"]["population"].update(
        {
            "tables": distribute_length_age(
                nasc_to_biology, analysis_dict["biology"]["proportions"], settings_dict
            )
        }
    )

    # Reapportion transect results to separate age-1 and age-2+ fish, generate age-1
    # abundance distributions for unaged fish, and generate biomass summary
    adult_transect, biomass_summary, unaged_age1_abundance = partition_transect_age(
        nasc_to_biology,
        analysis_dict["biology"]["weight"],
        settings_dict,
        analysis_dict["biology"]["population"],
        strata_adult_proportions,
    )

    # Append outputs to their appropriate positions within the analysis attribute
    # ---- Adult transect data
    analysis_dict["acoustics"].update({"adult_transect_df": adult_transect})
    # ---- Age-1 abundances (unaged fish)
    analysis_dict["biology"]["population"]["tables"]["abundance"].update(
        {"unaged_age1_abundance_df": unaged_age1_abundance}
    )

    # Return the updated analysis attribute and `biomass_summary`
    return biomass_summary, analysis_dict


def stratified_summary(
    analysis_dict: dict, results_dict: dict, spatial_dict: dict, settings_dict: dict
) -> tuple[pd.DataFrame, dict]:
    """
    Calculate stratified statistics (with and without resampling) for the entire survey.

    Parameters
    ----------
    analysis_dict: dict
        A dictionary containing processed biological and transect data.
    results_dict: dict
        A dictionary containing georeferenced results (i.e. kriging).
    spatial_dict: dict
        A dictionary containing dataframes that define the KS and INPFC stratum limits.
    settings_dict: dict
        Dictionary that contains all of the analysis settings that detail specific algorithm
        arguments and user-defined inputs.
    """

    # Toggle the correct transect data
    # ---- The 'original' acoustic transect data
    if settings_dict["dataset"] == "transect":
        # ---- Define the and prepare the processed and georeferenced transect data
        transect_data = edit_transect_columns(analysis_dict["transect"], settings_dict)
        # ---- Summarize transect spatial information
        # -------- Transect distances
        transect_summary = transect_spatial_features(transect_data)
        # -------- Counts and area coverage per stratum
        strata_summary = summarize_transect_strata(transect_summary)
    # ---- Kriged data
    elif settings_dict["dataset"] == "kriging":
        # ---- Convert the kriged mesh results into virtual transects
        transect_data, transect_summary, strata_summary = mesh_to_transects(
            results_dict["kriging"], spatial_dict, settings_dict
        )

    # Compute the stratified mean, variance, and coefficient of variation (CV)
    # ---- This includes the statistical (Gaussian) estimates (mean) and 95% confidence interval
    # ---- for each statistic
    replicates, stratified_results = stratified_transect_statistic(
        transect_data, transect_summary, strata_summary, settings_dict
    )

    # Update the analysis attribute with the resampled/bootstrapped replicates
    analysis_dict["stratified"].update(
        {f"{settings_dict['dataset']}": {"stratified_replicates_df": replicates}}
    )

    # Return the outputs
    return stratified_results, analysis_dict


def variogram_analysis(
    variogram_parameters: dict,
    default_variogram_parameters: dict,
    optimization_parameters: dict,
    initialize_variogram: dict,
    transect_dict: dict,
    settings_dict: dict,
    isobath_df: pd.DataFrame,
):

    # Validate the relevant empirical variogram parameters
    empirical_variogram_params = VariogramEmpirical.create(**settings_dict)

    # Initialize and validate the variogram model parameters
    valid_variogram_params = initialize_variogram_parameters(
        variogram_parameters, default_variogram_parameters
    )

    # Initialize and validate the optimization parameters
    valid_optimization_params = initialize_optimization_config(optimization_parameters)

    # Initialize and validate the initial values/boundary inputs
    valid_initial_values = initialize_initial_optimization_values(
        initialize_variogram, valid_variogram_params
    )

    # Prepare the transect data
    # ---- Create a copy of the transect dictionary
    transect_input = copy.deepcopy(transect_dict)
    # ---- Edit the transect data
    transect_data = edit_transect_columns(transect_input, settings_dict)

    # Standardize the transect coordinates, if necessary
    if settings_dict["standardize_coordinates"]:
        # ---- Transform geometry
        transect_data, _, _ = transform_geometry(transect_data, isobath_df, settings_dict)
        # ---- Print message if verbose
        if settings_dict["verbose"]:
            # ---- Print alert
            print(
                "Longitude and latitude coordinates (WGS84) converted to standardized "
                "coordinates (x and y)."
            )
    else:
        # ---- x
        transect_data["x"] = "longitude"
        # ---- y
        transect_data["y"] = "latitude"

    # Compute the empirical variogram
    lags, gamma_h, lag_counts, _ = empirical_variogram(
        transect_data, {**valid_variogram_params, **empirical_variogram_params}, settings_dict
    )

    # Least-squares fitting
    # ---- Consolidate the optimization dictionaries into a single one
    optimization_settings = {
        "parameters": valid_initial_values,
        "config": valid_optimization_params,
    }
    # ---- Optimize parameters
    best_fit_variogram, initial_fit, optimized_fit = optimize_variogram(
        lag_counts, lags, gamma_h, optimization_settings, **valid_variogram_params
    )

    # Return a dictionary of results
    return {
        "best_fit_parameters": best_fit_variogram,
        "initial_fit": {
            "parameters": dict(zip(initial_fit[0], initial_fit[1])),
            "MAD": initial_fit[2],
        },
        "optimized_fit": {
            "parameters": dict(zip(optimized_fit[0], optimized_fit[1])),
            "MAD": optimized_fit[2],
        },
    }


def krige(input_dict: dict, analysis_dict: dict, settings_dict: dict) -> tuple[pd.DataFrame, dict]:
    """
    Interpolate spatial data using ordinary kriging.

    Parameters
    ----------
    input_dict: dict
        A dictionary containing the loaded survey data.
    analysis_dict: dict
        A dictionary containing processed biological and transect data.
    settings_dict: dict
        Dictionary that contains all of the analysis settings that detail specific algorithm
        arguments and user-defined inputs.
    """

    # Validate cropping method parameters
    validated_cropping_methods = MeshCrop.create(**settings_dict["cropping_parameters"])
    # ---- Update the dictionary
    settings_dict["cropping_parameters"].update({**validated_cropping_methods})

    # Validate the variogram parameters
    valid_variogram_parameters = VariogramBase.create(**settings_dict["variogram_parameters"])
    # ---- Update the dictionary
    settings_dict["variogram_parameters"].update({**valid_variogram_parameters})

    # Validate kriging parameters
    valid_kriging_parameters = KrigingParameterInputs.create(
        **{**settings_dict["kriging_parameters"], **settings_dict["variogram_parameters"]}
    )
    # ---- Update the dictionary
    settings_dict["kriging_parameters"].update({**valid_kriging_parameters})

    # Validate the additional kriging arguments
    _ = KrigingAnalysis(**settings_dict["variogram_parameters"])

    # Extract kriging mesh data
    mesh_data = input_dict["statistics"]["kriging"]["mesh_df"]

    # Extract the reference grid (200 m isobath)
    isobath_data = input_dict["statistics"]["kriging"]["isobath_200m_df"]

    # Define the and prepare the processed and georeferenced transect data
    transect_data = edit_transect_columns(analysis_dict["transect"], settings_dict)

    # Add kriging parameters to the settings config
    settings_dict.update(
        {
            "kriging_parameters": {
                **input_dict["statistics"]["kriging"]["model_config"],
                **valid_kriging_parameters,
            },
            "variogram_parameters": {
                **settings_dict["variogram_parameters"],
                **valid_variogram_parameters,
            },
        },
    )

    # Crop the mesh grid if the kriged data will not be extrapolated
    if settings_dict["extrapolate"]:
        # ---- Else, extract original mesh dataframe
        mesh_df = mesh_data.copy()
        # ---- Extract longitude column name
        mesh_longitude = [col for col in mesh_df.columns if "lon" in col.lower()][0]
        # ---- Latitude
        mesh_latitude = [col for col in mesh_df.columns if "lat" in col.lower()][0]
        # ---- Rename the dataframe
        mesh_full = mesh_df.copy().rename(
            columns={f"{mesh_longitude}": "longitude", f"{mesh_latitude}": "latitude"}
        )
    else:
        # ---- Compute the cropped mesh
        mesh_full = crop_mesh(transect_data, mesh_data, validated_cropping_methods)
        if (settings_dict["verbose"]) and (
            validated_cropping_methods["crop_method"] == "convex_hull"
        ):
            # ---- Print alert
            print(
                f"Kriging mesh cropped to prevent extrapolation beyond the defined "
                f"`mesh_buffer_distance` value ({settings_dict['mesh_buffer_distance']} nmi)."
            )

    # Standardize the x- and y-coordinates, if necessary
    if settings_dict["standardize_coordinates"]:
        # ---- Transform transect data geometry (generate standardized x- and y-coordinates)
        transect_data, d_x, d_y = transform_geometry(transect_data, isobath_data, settings_dict)
        # ---- Transform mesh grid geometry (generate standardized x- and y-coordinates)
        mesh_full, _, _ = transform_geometry(mesh_full, isobath_data, settings_dict, d_x, d_y)
        if settings_dict["verbose"]:
            # ---- Print alert
            print(
                "Longitude and latitude coordinates (WGS84) converted to standardized "
                "coordinates (x and y)."
            )
    else:
        # ---- Else, duplicate the transect longitude and latitude coordinates as 'x' and 'y'
        # -------- x
        transect_data["x"] = transect_data["longitude"]
        # -------- y
        transect_data["y"] = transect_data["latitude"]
        # ---- Duplicate the mesh grid longitude and latitude coordinates as 'x' and 'y'
        # -------- x
        mesh_full["x"] = mesh_full["longitude"]
        # -------- y
        mesh_full["y"] = mesh_full["latitude"]
    # --- Append to the analysis attribute
    analysis_dict.update({"kriging": {"mesh_df": mesh_full, "transect_df": transect_data}})

    # Kriged results
    kriged_results = kriging(
        analysis_dict["kriging"]["transect_df"], analysis_dict["kriging"]["mesh_df"], settings_dict
    )

    # Stratified the kriging mesh
    kriged_results["mesh_results_df"] = stratify_mesh(
        input_dict, kriged_results["mesh_results_df"], settings_dict
    )

    # Return kriged (interpolated) results
    return kriged_results, analysis_dict


def apportion_kriged_values(
    analysis_dict: dict, kriged_mesh: pd.DataFrame, settings_dict: dict
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Apportion spatially distributed kriged biological estimates over length and age bins.

    Parameters
    ----------
    analysis_dict: dict
        A dictionary containing processed biological and transect data.
    kriged_mesh: pd.DataFrame
        A dataframe containing the spatially distributed kriged estimates.
    settings_dict: dict
        Dictionary that contains all of the analysis settings that detail specific algorithm
        arguments and user-defined inputs.
    """

    # Sum the kriged weights for each stratum
    # ---- Extract stratum column name
    stratum_col = settings_dict["stratum_name"]
    # ---- Extract the biological variable (independent of area)
    biology_col = settings_dict["variable"].replace("_density", "")
    # ---- Sum the kriged values for each stratum
    kriged_strata = kriged_mesh.groupby([stratum_col], observed=False)[biology_col].sum()

    # Extract the weight proportions from the analysis object
    proportions_dict = analysis_dict["transect"]["biology"]["proportions"]
    # ---- Aged
    aged_proportions = proportions_dict["weight"]["aged_weight_proportions_df"].copy()
    # ---- Unaged
    unaged_proportions = proportions_dict["weight"]["unaged_weight_proportions_df"].copy()
    # ---- Aged-unaged sexed weight proportions
    aged_unaged_sex_proportions = proportions_dict["weight"][
        "aged_unaged_sex_weight_proportions_df"
    ].copy()[[stratum_col, "sex", "weight_proportion_overall_unaged"]]

    # Compute the apportioned unaged kriged biological values per stratum
    # ---- Merge the unaged proportions
    unaged_sexed_apportioned = unaged_proportions.merge(aged_unaged_sex_proportions)
    # ---- Set index to stratum column
    unaged_sexed_apportioned.set_index([stratum_col], inplace=True)
    # ---- Append the stratum-aggregated values
    unaged_sexed_apportioned[f"{biology_col}_apportioned_unaged"] = (
        unaged_sexed_apportioned["weight_proportion"]
        * unaged_sexed_apportioned["weight_proportion_overall_unaged"]
        * kriged_strata
    )

    # Distribute biological values over the overall proportions (i.e. relative to aged and unaged
    # fish) for aged fish
    # ---- Set index to stratum column
    aged_proportions.set_index([stratum_col], inplace=True)
    # ---- Compute the distributed values
    aged_proportions[f"{biology_col}_apportioned"] = (
        aged_proportions["weight_proportion_overall"] * kriged_strata
    ).fillna(0.0)

    # Distribute the aged biological distributions over unaged length distributions to estimate
    # aged distributions
    # ---- Pivot aged data
    aged_pivot = aged_proportions.reset_index().pivot_table(
        index=["sex", "length_bin"],
        columns=["age_bin"],
        values=f"{biology_col}_apportioned",
        aggfunc="sum",
        observed=False,
    )
    # ---- Calculate the total biomass values for each sex per length bin
    aged_length_totals = aged_pivot.sum(axis=1).unstack("sex")
    # ---- Pivot unaged data
    unaged_pivot = unaged_sexed_apportioned.reset_index().pivot_table(
        index=["length_bin"],
        columns=["sex"],
        values=f"{biology_col}_apportioned_unaged",
        aggfunc="sum",
        observed=False,
    )
    # ---- Calculate the new unaged biological values distributed over age
    unaged_apportioned_values = (
        unaged_pivot * aged_pivot.unstack("sex") / aged_length_totals
    ).fillna(0)

    # Imputation is required when unaged values are present but aged values are absent at shared
    # length bins! This requires an augmented implementation to address this accordingly
    kriged_table = impute_kriged_values(
        aged_pivot, unaged_pivot, aged_length_totals, unaged_apportioned_values, settings_dict
    )
    # ---- Duplicate so there is an 'all' category
    kriged_full_table = pd.concat([kriged_table, kriged_table.copy().assign(sex="all")])
    # ---- Consolidate
    kriging_full_table = (
        kriged_full_table.groupby(["length_bin", "age_bin", "sex"], observed=False)[
            f"{biology_col}_apportioned"
        ]
        .sum()
        .reset_index(name=f"{biology_col}_apportioned")
    )

    # Additional reapportionment if age-1 fish are excluded
    if settings_dict["exclude_age1"]:
        kriging_full_table = reallocate_kriged_age1(kriging_full_table, settings_dict)
        # ---- Validate that apportioning age-1 values over all adult values did not 'leak'
        # -------- Previous apportioned totals by sex
        previous_totals = kriged_full_table.groupby(["sex"])[f"{biology_col}_apportioned"].sum()
        # -------- New apportioned totals by sex
        new_totals = kriging_full_table.groupby(["sex"])[f"{biology_col}_apportioned"].sum()
        # -------- Check (1 kg tolerance)
        if np.any((previous_totals - new_totals) > 1e-6):
            warnings.warn(
                f"Apportioned kriged {biology_col} for age-1 not fully distributed over all age-2+ "
                f"age bins."
            )

    # Check equality between original kriged estimates and (imputed) apportioned estimates
    if (kriged_table[f"{biology_col}_apportioned"].sum() - kriged_strata.sum()) > 1e-6:
        # ---- If not equal, generate warning
        warnings.warn(
            f"Apportioned kriged {biology_col} does not equal the total kriged mesh "
            f"{biology_col}! Check for cases where kriged values may only be present in aged "
            f"(`self.results['kriging']['tables']['aged_tbl']`) or unaged ("
            f"(`self.results['kriging']['tables']['unaged_tbl']`) distributions for each sex."
        )

    # Return output
    return aged_pivot, unaged_pivot, kriging_full_table
