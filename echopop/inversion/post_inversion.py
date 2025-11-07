from typing import Literal

import numpy as np
import pandas as pd

from .utils import _extract_parameters_optimized


def estimate_population(
    inverted_data: pd.DataFrame,
    nasc_data: pd.DataFrame,
    density_sw: float,
    reference_frequency: float,
    aggregate_method: Literal["cells", "interval", "transect"],
):
    """
    Estimate population characteristics from inverted acoustic data.

    This function computes population-level estimates including body size,
    density, and biomass from inverted acoustic scattering parameters.
    It combines acoustic data with biological parameter estimates to
    provide comprehensive population assessments.

    Parameters
    ----------
    inverted_data : pd.DataFrame
        DataFrame containing inverted scattering parameters with columns
        including 'parameters' containing biological model parameters
    nasc_data : pd.DataFrame
        DataFrame containing Nautical Area Scattering Coefficient (NASC)
        values with frequency as column level
    density_sw : float
        Seawater density in kg/m³, typically ~1026 kg/m³
    reference_frequency : float
        Reference acoustic frequency in Hz for population estimation
    aggregate_method : Literal["cells", "interval", "transect"]
        Method for spatial aggregation of population estimates:
        - "cells": Individual acoustic cells
        - "interval": Depth/range intervals
        - "transect": Survey transect lines

    Returns
    -------
    pd.DataFrame
        Combined DataFrame containing:
        - nasc: NASC values at reference frequency
        - number_density: Estimated areal number density (individuals/m²)
        - biomass_density: Estimated areal biomass density (kg/m²)
        - Biological parameters: length_mean, radius_mean, body_volume,
          body_density, body_weight, etc.

    Notes
    -----
    The function assumes organisms can be modeled as uniformly bent cylinders
    for volume estimation. Body weight calculation uses:

    .. math::
        W = \\rho_{sw} \\cdot g \\cdot \\pi r^2 L

    where W is body weight, ρ_sw is seawater density, g is the relative
    density factor, r is mean radius, and L is mean length.
    """

    # Unpackage the parameters using optimized extraction
    parameters = _extract_parameters_optimized(inverted_data)

    # Compute the average body radius for a single animal
    parameters["radius_mean"] = parameters["length_mean"] / parameters["length_radius_ratio"]

    # Estimate the average volume assuming an uniformly bent cylinder [single animal]
    parameters["body_volume"] = np.pi * parameters["radius_mean"] ** 2 * parameters["length_mean"]

    # Compute the animal body density [single animal]
    parameters["body_density"] = density_sw * parameters["g"]

    # Compute average body weight [single animal]
    parameters["body_weight"] = parameters["body_volume"] * parameters["body_density"]

    # Subset the dataset further to the single frequency
    acoustic_data = (
        inverted_data.loc[
            :, inverted_data.columns.get_level_values("frequency") == reference_frequency
        ]
    ).droplevel("frequency", axis=1)

    # Repeat the same for the coordinates
    reference_nasc = (
        nasc_data.loc[:, nasc_data.columns.get_level_values("frequency") == reference_frequency]
    )[reference_frequency].to_frame("nasc")

    # Calculate the areal number density
    if aggregate_method == "transect":
        reference_nasc.loc[:, "number_density"] = calculate_transect_number_density(
            acoustic_data, reference_nasc, parameters
        ).to_numpy()
    else:
        reference_nasc["number_density"] = calculate_intervals_number_density(
            acoustic_data, parameters
        )

    # Compute areal biomass density
    reference_nasc["biomass_density"] = (
        reference_nasc["number_density"] * parameters["body_weight"]
    ).fillna(0.0)

    # Find shared indices
    idx = list(set(reference_nasc.index.names).intersection(parameters.index.names))

    # Reset indices
    reference_nasc = reference_nasc.reset_index().set_index(idx)

    # Concatenate the parameters
    return pd.concat(
        [
            reference_nasc,
            parameters.reset_index()
            .drop("number_density", axis=1)
            .set_index(idx)
            .reindex(reference_nasc.index),
        ],
        axis=1,
    )


def calculate_transect_number_density(
    acoustic_data: pd.DataFrame,
    reference_nasc: pd.DataFrame,
    parameters: pd.DataFrame,
):
    """
    Compute areal number density aggregated by survey transect.

    This function estimates the areal number density of organisms by
    weighting volumetric densities by NASC contributions and aggregating
    across survey transects. The method accounts for varying layer
    thickness and spatial sampling intensity.

    Parameters
    ----------
    acoustic_data : pd.DataFrame
        DataFrame containing acoustic measurements including:
        - nasc: Nautical Area Scattering Coefficient values
        - thickness_mean: Mean layer thickness in meters
    reference_nasc : pd.DataFrame
        DataFrame with NASC values at reference frequency for weighting
    parameters : pd.DataFrame
        DataFrame containing biological parameters including:
        - number_density: Volumetric number density (individuals/m³)

    Returns
    -------
    pd.Series
        Areal number density in individuals/m² aggregated by transect,
        converted from nautical miles to meters (factor of 1852²)

    Notes
    -----
    The areal number density calculation follows:

    .. math::
        N_A = N_V \\cdot w \\cdot h \\cdot n \\cdot 1852^2

    where N_A is areal density, N_V is volumetric density, w is NASC weight,
    h is layer thickness, n is interval count, and 1852² converts from
    nautical miles² to meters².

    NASC weighting ensures proper spatial representation:

    .. math::
        w = \\frac{\\text{NASC}_{ref}}{\\text{NASC}_{acoustic}}
    """

    # Create copies
    acoustic_df = acoustic_data.copy()
    reference_df = reference_nasc.copy()

    # Find shared indices
    idx = list(set(reference_nasc.index.names).intersection(acoustic_data.index.names))

    # Reset the indices for each
    acoustic_df = acoustic_df.reset_index().set_index(idx)
    reference_df = reference_df.reset_index().set_index(idx)

    # Compute the NASC weights
    reference_df["nasc_weight"] = reference_df["nasc"] / acoustic_df["nasc"]

    # Get the counts of each transect
    interval_counts = reference_df.groupby(level="transect_num")["nasc"].count()

    # Compute the areal number density
    areal_number_density = (
        parameters["number_density"]
        * reference_df["nasc_weight"]
        * acoustic_df["thickness_mean"]
        * interval_counts
    ).fillna(0.0) * 1852**2

    # Return
    return areal_number_density


def calculate_intervals_number_density(
    acoustic_data: pd.DataFrame,
    parameters: pd.DataFrame,
):
    """
    Compute areal number density for individual depth/range intervals.

    This function calculates areal number density by converting volumetric
    density estimates to areal density using layer thickness information.
    Unlike transect aggregation, this preserves interval-level resolution.

    Parameters
    ----------
    acoustic_data : pd.DataFrame
        DataFrame containing acoustic measurements with:
        - thickness_mean: Mean layer thickness in meters for each interval
    parameters : pd.DataFrame
        DataFrame containing biological parameters with:
        - number_density: Volumetric number density (individuals/m³)

    Returns
    -------
    pd.Series
        Areal number density in individuals/m² for each interval,
        converted from nautical miles to meters (factor of 1852²)

    Notes
    -----
    The conversion from volumetric to areal density is:

    .. math::
        N_A = N_V \\cdot h \\cdot 1852^2

    where N_A is areal density (individuals/m²), N_V is volumetric density
    (individuals/m³), h is layer thickness (m), and 1852² converts from
    nautical miles² to meters².

    This method maintains spatial resolution at the interval level, making
    it suitable for detailed vertical distribution analysis.
    """

    # Calculate the layer density
    layer_density = parameters["number_density"] * acoustic_data["thickness_mean"]

    # Return
    return layer_density * 1852**2
