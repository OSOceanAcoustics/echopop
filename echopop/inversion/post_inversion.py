from typing import Literal

import numpy as np
import xarray as xr

from .utils import _extract_parameters_optimized


def estimate_population(
    inverted_data: xr.Dataset,
    nasc_data: xr.DataArray,
    density_sw: float,
    reference_frequency: float,
    aggregate_method: Literal["cells", "interval", "transect"],
) -> xr.Dataset:
    """
    Estimate population characteristics from inverted acoustic data.

    This function computes population-level estimates including body size, density, and biomass
    from inverted acoustic scattering parameters. It combines acoustic data with biological
    parameter estimates to provide comprehensive population assessments.

    Parameters
    ----------
    inverted_data : xr.Dataset
        Dataset containing inverted scattering parameters with 'parameters' coordinate containing
        biological model parameters.
    nasc_data : xr.DataArray
        DataArray containing Nautical Area Scattering Coefficient (NASC) values with frequency
        dimension.
    density_sw : float
        Seawater density in kg/m³, typically ~1026 kg/m³
    reference_frequency : float
        Reference acoustic frequency in Hz for population estimation.
    aggregate_method : Literal["cells", "interval", "transect"]
        Method for spatial aggregation of population estimates:

        - "cells": Individual acoustic cells

        - "interval": Depth/range intervals

        - "transect": Survey transect lines

    Returns
    -------
    xr.Dataset
        Dataset containing:

        - nasc: NASC values at reference frequency

        - number_density: Estimated areal number density (individuals/m²)

        - biomass_density: Estimated areal biomass density (kg/m²)

        - Biological parameters: length_mean, radius_mean, body_volume, body_density,
          body_weight, etc.

    Notes
    -----
    The function assumes organisms can be modeled as uniformly bent cylinders for volume
    estimation. Body weight calculation uses:

    .. math::
        W = \\rho_{sw} \\cdot g \\cdot \\pi r^2 L

    where W is body weight, ρ_sw is seawater density, g is the relative density factor, r is mean
    radius, and L is mean length.
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

    # Subset the datasets further to the single frequency
    acoustic_data = inverted_data.sel(frequency=reference_frequency)
    reference_nasc = nasc_data.sel(frequency=reference_frequency)

    # Calculate the areal number density
    if aggregate_method == "transect":
        number_density = calculate_transect_number_density(
            acoustic_data, reference_nasc, parameters
        )
    else:
        number_density = calculate_intervals_number_density(acoustic_data, parameters)

    # Convert NASC to drop 'point' coordinate
    nasc_cnv = reference_nasc.to_series().to_frame("nasc").reset_index()

    # Convert the number densities and parameter sets and drop the point coordinate
    number_density_cnv = number_density.to_series().to_frame("number_density").reset_index()

    # Get the shared column names (excluding 'point')
    shared_cols = [
        col for col in number_density_cnv.columns if col in nasc_cnv.columns and col != "point"
    ]

    # Do a left-join to broadcast density values
    merged_dens = nasc_cnv.merge(number_density_cnv, on=shared_cols, how="left").fillna(0.0)

    # Reformat to Dataset structure
    # ---- Get overlapping coordinates
    overlapping_coords = [c for c in reference_nasc.coords if c in merged_dens.columns]
    # ---- Find valid ordering
    coord_order = [
        c
        for c in ["transect_num", "longitude", "latitude", "interval", "layer"]
        if c in overlapping_coords
    ]
    # ---- Set the index
    output_indexed = merged_dens.set_index(coord_order)
    parameters_indexed = parameters.reset_index().set_index(
        [c for c in parameters.index.names if c in coord_order]
    )

    # Calculate biomass
    output_indexed["biomass_density"] = (
        output_indexed["number_density"] * parameters_indexed["body_weight"]
    ).fillna(0.0)

    # Create Dataset by converting each column separately
    data_vars = {}
    for col in ["nasc", "number_density", "biomass_density", "parameters"]:
        if col in output_indexed.columns:
            data_vars[col] = (["point"], output_indexed[col].values)
    # ---- Create coordinates from the MultiIndex using from_pandas_multiindex
    coords = xr.Coordinates.from_pandas_multiindex(output_indexed.index, "point")

    return xr.Dataset(data_vars, coords=coords)


def calculate_transect_number_density(
    acoustic_data: xr.Dataset,
    reference_nasc: xr.DataArray,
    parameters: xr.Dataset,
):
    """
    Compute areal number density aggregated by survey transect.

    This function estimates the areal number density of organisms by weighting volumetric densities
    by NASC contributions and aggregating across survey transects. The method accounts for varying
    layer thickness and spatial sampling intensity.

    Parameters
    ----------
    acoustic_data : xr.Dataset
        Dataset containing acoustic measurements including:

        - nasc: Nautical Area Scattering Coefficient values

        - thickness_mean: Mean layer thickness in meters

    reference_nasc : xr.DataArray
        DataArray with NASC values at reference frequency for weighting.

    parameters : xr.Dataset
        Dataset containing biological parameters including:

        - number_density: Volumetric number density (individuals/m³)

    Returns
    -------
    xr.DataArray
        Areal number density in individuals/m² aggregated by transect, converted from nautical
        miles to meters (factor of 1852²).

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

    # Get the transects that exist in acoustic_data
    valid_transects = acoustic_data["transect_num"].values

    # Filter reference_nasc to only include valid transects
    reference_filtered = reference_nasc.where(
        reference_nasc["transect_num"].isin(valid_transects), drop=True
    )

    # Reindex acoustic data to match reference_nasc structure
    acoustic_aligned = acoustic_data.sel(transect_num=reference_filtered["transect_num"])
    parameters_aligned = parameters.to_xarray().sel(transect_num=reference_filtered["transect_num"])

    # Compute the NASC weights
    nasc_weight = reference_filtered.values / acoustic_aligned["nasc"].values

    # Count intervals per transect
    interval_counts = reference_filtered.groupby("transect_num").count()

    # Broadcast interval counts back to reference dimensions
    interval_counts_broadcast = xr.DataArray(
        [
            interval_counts.sel(transect_num=t).values
            for t in reference_filtered["transect_num"].values
        ],
        coords=reference_filtered.coords,
        dims=reference_filtered.dims,
    )

    # Compute the areal number density using values
    areal_number_density_vals = (
        parameters_aligned["number_density"].values
        * nasc_weight
        * acoustic_aligned["thickness_mean"].values
        * interval_counts_broadcast.values
    ) * 1852**2

    # Create output DataArray with reference_filtered structure
    areal_number_density = xr.DataArray(
        areal_number_density_vals, coords=reference_filtered.coords, dims=reference_filtered.dims
    )

    # Return
    return areal_number_density


def calculate_intervals_number_density(
    acoustic_data: xr.Dataset,
    parameters: xr.Dataset,
) -> xr.DataArray:
    """
    Compute areal number density for individual depth/range intervals.

    This function calculates areal number density by converting volumetric density estimates to
    areal density using layer thickness information. Unlike transect aggregation, this preserves
    interval-level resolution.

    Parameters
    ----------
    acoustic_data : xr.Dataset
        Dataset containing acoustic measurements with:

        - thickness_mean: Mean layer thickness in meters for each interval

    parameters : xr.Dataset
        Dataset containing biological parameters with:

        - number_density: Volumetric number density (individuals/m³)

    Returns
    -------
    xr.DataArray
        Areal number density in individuals/m² for each interval, converted from nautical miles to
        meters (factor of 1852²).

    Notes
    -----
    The conversion from volumetric to areal density is:

    .. math::
        N_A = N_V \\cdot h \\cdot 1852^2

    where N_A is areal density (individuals/m²), N_V is volumetric density (individuals/m³), h is
    layer thickness (m), and 1852² converts from nautical miles² to meters².

    This method maintains spatial resolution at the interval level, making it suitable for
    detailed vertical distribution analysis.
    """

    # Calculate the layer density using values to avoid index conflict
    layer_density_vals = (
        parameters["number_density"].values * acoustic_data["thickness_mean"].values
    ) * 1852**2

    # Create output DataArray with acoustic_data structure
    layer_density = xr.DataArray(
        layer_density_vals, coords=acoustic_data.coords, dims=acoustic_data.dims
    )

    # Return
    return layer_density
