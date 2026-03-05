"""
This module corrects biological distributions based on length-based net selectivity.
"""

from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import xarray as xr

from ..validators import ValidateSelectivityParams


def logistic_selectivity(
    length_bins: xr.DataArray,
    slope: Optional[float] = None,
    intercept: Optional[float] = None,
    l50: Optional[float] = None,
    sr: Optional[float] = None,
    minimum_selectivity: float = 1e-6,
) -> xr.DataArray:
    """
    Compute logistic net selectivity by length bin.

    Selectivity can be parameterized using either:

    1) slope/intercept:

       S(l) = 1 / (1 + exp(-(intercept + slope*l)))

    2) L50%/SR:

       S(l) = (1 + exp(k*(L50 - l)/SR))^-1, where k = 2*log(3)

    Parameters
    ----------
    length_bins : xr.DataArray
        Coordinate-like DataArray containing interval-based length bins.
    slope : float, optional
        Logistic slope coefficient.
    intercept : float, optional
        Logistic intercept coefficient.
    l50 : float, optional
        Length-at-50%-retention parameter.
    sr : float, optional
        Selection range (length interval between 25% and 75% retention).
    minimum_selectivity : float, default 1e-6
        Lower bound applied to selectivity values to avoid division by zero.

    Returns
    -------
    xr.DataArray
        Selectivity values indexed by ``length_bin``.

    Raises
    ------
    ValueError
        If parameterization inputs are invalid.

    Notes
    -----

    - This function supports both common parameterizations used in selectivity studies.

    - Exactly one parameter pair must be provided: ``(slope, intercept)`` or ``(l50, sr)``.

    References
    ----------
    Wileman, D. A., Ferro, R. S. T., Fonteyne, R., & Millar, R. B. (Eds.). (1996). *Manual of
    methods of measuring the selectivity of towed fishing gears*. ICES Cooperative Research Report
    No. 215.
    """
    # Convert interval bins to midpoint lengths
    length_mid = xr.DataArray(
        [interval.mid for interval in length_bins.values],
        dims=["length_bin"],
        coords={"length_bin": length_bins.values},
    )

    # Compute selectivity
    uses_slope_intercept = (slope is not None) and (intercept is not None)
    if uses_slope_intercept:
        selectivity = 1.0 / (1.0 + np.exp(-(intercept + slope * length_mid)))
    else:
        if sr is None or sr <= 0:
            raise ValueError("'sr' must be > 0.")
        k = 2.0 * np.log(3.0)
        selectivity = (1.0 + np.exp(k * (l50 - length_mid) / sr)) ** -1

    # ---- Apply lower bound
    return xr.where(selectivity < minimum_selectivity, minimum_selectivity, selectivity)


def _to_l50_sr(
    raw: Dict[str, float],
    k: float,
) -> Tuple[float, float]:
    """
    Validate and convert a selectivity parameter dict to an ``(l50, sr)`` pair.

    Parameters
    ----------
    raw : dict
        Must contain exactly one of ``{'l50', 'sr'}`` or ``{'slope', 'intercept'}``.
    k : float
        Logistic steepness constant, typically ``2 * ln(3)``.

    Returns
    -------
    tuple of (float, float)
        ``(l50, sr)`` ready for use in the logistic selectivity formula.
    """
    validated = ValidateSelectivityParams.create(**raw)
    if validated.get("l50") is not None:
        return float(validated["l50"]), float(validated["sr"])
    slope = float(validated["slope"])
    intercept = float(validated["intercept"])
    return -intercept / slope, k / slope


def assign_selectivity_expansion(
    data: pd.DataFrame,
    gear_selectivity_params: Union[Dict[str, float], Dict[str, Dict[str, float]]],
    minimum_selectivity: float = 1e-12,
    *,
    gear_col: str = "gear",
    length_col: str = "length",
) -> pd.DataFrame:
    """
    Compute per-fish logistic selectivity expansion factors and add them as a new column.

    The caller is responsible for ensuring ``data`` already contains both a length column and a
    gear column (e.g., via a prior merge of the individual-level and haul-level frames). This
    function only performs the selectivity computation.

    Parameters
    ----------
    data : pd.DataFrame
        Individual-level DataFrame containing at least ``length_col`` and ``gear_col``.
    gear_selectivity_params : dict
        Accepts three shapes:

        1. **Flat** — a single parameter set applied uniformly to every row:

           ``{'l50': float, 'sr': float}``  or  ``{'slope': float, 'intercept': float}``

        2. **Single-gear nested** — one gear name mapped to its parameter dict:

           ``{'AWT': {'l50': float, 'sr': float}}``

        3. **Multi-gear nested** — multiple gear names each mapped to their own parameter dict:

           ``{'AWT': {'l50': float, 'sr': float}, 'MFT': {'slope': float, 'intercept': float}}``

        In the nested forms (2 and 3) each value must contain **exactly one** complete set.
        Mixed or incomplete sets raise ``ValueError`` via ``ValidateSelectivityParams``.
        When a row's gear is not found in the nested mapping, its
        ``selectivity_expansion`` will be ``NaN``.
    minimum_selectivity : float, default 1e-12
        Lower bound applied to ``S(L)`` before inversion to prevent division by zero.
    gear_col : str, default ``'gear'``
        Column identifying the gear type for each individual.
    length_col : str, default ``'length'``
        Column containing individual fish lengths.

    Returns
    -------
    pd.DataFrame
        Copy of ``data`` with one new column added:

        - ``selectivity_expansion`` : ``1 / S(L)`` — the per-fish selectivity expansion factor.

    Raises
    ------
    ValueError
        If any entry in ``gear_selectivity_params`` fails parameter validation.

    Notes
    -----
    Slope/intercept parameterizations are converted to ``(l50, sr)`` once before any per-row
    work:

    .. math::

        L_{50} = -\\beta_0 / \\beta_1, \\quad SR = 2\\ln(3) / \\beta_1

    The per-fish selectivity expansion factor is then:

    .. math::

        \\frac{1}{\\max\\left(\\left[1 + e^{\\frac{2\\ln(3)(L_{50} - L)}{SR}}\\right]^{-1},\\ S_{\\min}\\right)}

    Rows whose gear is not found in ``gear_selectivity_params`` will have ``NaN`` for
    ``selectivity_expansion``.
    """
    _k = 2.0 * np.log(3.0)

    # Detect whether a single flat dict or a gear-keyed nested dict was supplied
    is_nested = any(isinstance(v, dict) for v in gear_selectivity_params.values())

    if not is_nested:
        # Single parameter set — apply uniformly to every row
        l50, sr = _to_l50_sr(gear_selectivity_params, _k)  # type: ignore[arg-type]
        L = data[length_col].to_numpy(dtype=float)
        l50_arr = np.full(len(L), l50)
        sr_arr = np.full(len(L), sr)
    else:
        # Validate and normalise all gear params to (l50, sr)
        gear_l50_sr: Dict[str, Tuple[float, float]] = {}
        for gear_name, raw_params in gear_selectivity_params.items():
            gear_l50_sr[gear_name] = _to_l50_sr(raw_params, _k)  # type: ignore[arg-type]

        # Map each row's gear to its (l50, sr) pair, then compute S(L) vectorised
        l50_arr = data[gear_col].map(
            lambda g: gear_l50_sr[g][0] if g in gear_l50_sr else np.nan
        ).to_numpy(dtype=float)
        sr_arr = data[gear_col].map(
            lambda g: gear_l50_sr[g][1] if g in gear_l50_sr else np.nan
        ).to_numpy(dtype=float)
        L = data[length_col].to_numpy(dtype=float)

    S = np.maximum((1.0 + np.exp(_k * (l50_arr - L) / sr_arr)) ** -1, minimum_selectivity)

    out = data.copy()
    out["selectivity_expansion"] = 1.0 / S
    return out


def _extract_proportion_arrays(
    number_data: Union[xr.Dataset, Dict[str, Union[xr.Dataset, xr.DataArray]]],
) -> List[xr.DataArray]:
    """
    Internal utility function for validating and extracting the input proportions.
    """
    # Normalize to a list of candidate DataArrays
    arrays: List[xr.DataArray] = []

    # Dataset
    if isinstance(number_data, xr.Dataset):
        if "proportion_overall" in number_data:
            arrays = [number_data["proportion_overall"]]
        else:
            arrays = [array for array in number_data.data_vars.values()]
    # Dict[Dataset]
    elif isinstance(number_data, dict):
        for value in number_data.values():
            if "proportion_overall" in value:
                arrays.append(value["proportion_overall"])
            else:
                arrays.extend([array for array in value.data_vars.values()])
    else:
        raise TypeError("'number_data' must be an xarray Dataset, or dictionary of Datasets.")

    # Keep only arrays that include length_bin
    arrays = [array for array in arrays if "length_bin" in array.dims]
    if len(arrays) == 0:
        raise ValueError("No proportion arrays with 'length_bin' were found in 'number_data'.")

    # ---- Coerce to numeric arrays for stable arithmetic (handles object arrays)
    arrays = [xr.where(array.isnull(), 0.0, array).astype(float) for array in arrays]

    return arrays


def _sum_aligned(arrays: List[xr.DataArray]) -> xr.DataArray:
    """
    Internal utility function for summing aligned DataArrays.
    """
    # Return if already a single array
    if len(arrays) == 1:
        return arrays[0].fillna(0.0)

    # Align arrays
    aligned_arrays = xr.align(*arrays, join="outer", fill_value=0.0)
    # ---- Replace in-array NaN values before addition to avoid NaN propagation
    aligned_arrays = tuple(array.fillna(0.0) for array in aligned_arrays)

    # Iterate and sum
    total = aligned_arrays[0]
    for array in aligned_arrays[1:]:
        total = total + array

    # Return with NaN filled with 0.0
    return total.fillna(0.0)


def _normalize_by_dims(da: xr.DataArray, dims: List[str]) -> xr.DataArray:
    """
    Internal utility function for normalizing the proportions over the defined dimensions such that
    they always sum to 1.0. This is a required calculation for the corrected proportions since the
    inflation of smaller bins (e.g., tails) results in the sum exceeding 1.0.
    """
    # Normalize over dimensions not explicitly retained
    sum_dims = [dim for dim in da.dims if dim not in dims]
    if len(sum_dims) == 0:
        return da.fillna(0.0)

    # Normalize over the defined dimensions
    return (da / da.sum(dim=sum_dims)).fillna(0.0)


def _construct_corrected_distribution(
    combined_marginal: xr.DataArray,
    normalization_dims: List[str],
    minimum_selectivity: Optional[float] = None,
    slope: Optional[float] = None,
    intercept: Optional[float] = None,
    l50: Optional[float] = None,
    sr: Optional[float] = None,
) -> xr.DataArray:
    """
    Internal utility function for constructing the selectivity-corrected length distribution.
    """
    # Oserved length distribution
    obs = _normalize_by_dims(combined_marginal, normalization_dims)

    # Selectivity-corrected length distribution
    selectivity = logistic_selectivity(
        length_bins=obs["length_bin"],
        slope=slope,
        intercept=intercept,
        l50=l50,
        sr=sr,
        minimum_selectivity=minimum_selectivity,
    )

    # Pure weighting correction (zeros remain zeros)
    return _normalize_by_dims(obs * (1.0 / selectivity), normalization_dims)


def _compute_length_age_marginal(
    proportion_arrays: List[xr.DataArray], keep_dims: List[str]
) -> Optional[xr.DataArray]:
    """
    Internal utility function for calculate the length-age marginals.
    """
    # Keep only arrays that include age information
    age_arrays = [array for array in proportion_arrays if "age_bin" in array.dims]
    if len(age_arrays) == 0:
        return None

    # Marginalize each array to length-age proportions
    length_age_arrays = [
        array.sum(
            dim=[dim for dim in array.dims if dim not in [*keep_dims, "length_bin", "age_bin"]]
        )
        for array in age_arrays
    ]

    # Align and combine all length-age marginals
    return _sum_aligned(length_age_arrays)


def _build_adjusted_outputs(
    observed_len_proportion: xr.DataArray,
    corrected_len_proportion: xr.DataArray,
    normalization_dims: List[str],
    length_age_marginal: Optional[xr.DataArray] = None,
    corrected_length_age_proportion: Optional[xr.DataArray] = None,
) -> xr.Dataset:
    """
    Internal utility function for formatting the output Dataset structure of the length-based
    selectivity-corrected proportions.
    """

    # Build corrected length proportions
    output = xr.Dataset(
        {
            "length_proportions_observed": observed_len_proportion,
            "length_proportions_corrected": corrected_len_proportion,
        }
    )

    # Optionally build corrected length-age distributions
    if length_age_marginal is not None and corrected_length_age_proportion is not None:
        output["length_age_proportions_observed"] = _normalize_by_dims(
            length_age_marginal,
            normalization_dims,
        )
        output["length_age_proportions_corrected"] = corrected_length_age_proportion

    # Canonical downstream variable names
    corrected_overall = (
        output["length_age_proportions_corrected"]
        if "length_age_proportions_corrected" in output
        else output["length_proportions_corrected"]
    )
    output["proportion_overall"] = corrected_overall
    output["proportion"] = corrected_overall

    return output


def correct_number_proportions(
    number_data: Union[xr.Dataset, Dict[str, Union[xr.Dataset]]],
    selectivity_params: Dict[str, float],
    stratum_dim: Optional[Union[str, List[str]]] = None,
) -> xr.Dataset:
    """
    Adjust number proportions using a logistic net selectivity correction.

    This function applies selectivity to a single combined length distribution
    (aged + unaged), then uses aged age-composition information to reconstruct
    adjusted aged and unaged outputs.

    Parameters
    ----------
    number_data : xr.Dataset, Dict[str, xr.Dataset]
        Observed length-binned counts in the form of proportions (i.e., normalized to sum to 1).
        This can either be a single Dataset, or a dictionary of Datasets. Each Dataset must, at a
        minimum, include the coordinate ``length_bin`` and the variable ``proportion_overall``.
    selectivity_params : Dict[str, float]
        A dictionary of selectivity parameters that must contain either:

            - 'intercept' (float): Logistic regression intercept coefficient.

            - 'slope' (float): Logistic regression slope coefficient.

        or:

            - 'l50' (float): Length-at-50%-retention parameter (the length at which selectivity is
               0.5).

            - 'sr' (float): Selection range (the length interval between 25% and 75% retention).

    minimum_selectivity : float, default 1e-12
        Lower bound applied to selectivity values to avoid division-by-zero errors.
    stratum_dim : str or list of str, optional
        Dimensions that define the normalization groups. When provided, both observed and adjusted
        length proportions are normalized to sum to 1 within each stratum. Other inferred
        coordinates are still retained in the output arrays.

    Notes
    -----

    - Selectivity is applied strictly as a weighting factor.

    - Zero-probability bins remain zero (no gap-filling behavior).

    Returns
    -------
    xr.Dataset
        Dataset containing adjusted number distributions with key DataArrays:

        - ``length_proportion_adjusted`` over ``length_bin``

        - ``length_age_proportion_adjusted`` over ``length_bin`` and ``age_bin`` (when available)

        Additional DataArrays include adjusted counts and overall proportions.
    """

    # Parse and validate proportion arrays
    proportion_arrays = _extract_proportion_arrays(number_data)

    # Resolve retained dimensions and build combined length marginal
    all_dims = sorted(set().union(*[set(array.dims) for array in proportion_arrays]))
    keep_dims = [dim for dim in all_dims if dim not in ["length_bin", "age_bin"]]

    # Resolve normalization dimensions
    # ---- Stratified: each stratum sums to 1
    if stratum_dim:
        # ---- Convert to list if needed
        normalization_dims = [stratum_dim] if isinstance(stratum_dim, str) else list(stratum_dim)
    # ---- Grouped: each combination of non-length, non-age coordinates sum to 1
    else:
        normalization_dims = keep_dims

    # Marginalize each array to their respective length-bin proportions
    length_arrays = [
        array.sum(dim=[dim for dim in array.dims if dim not in [*keep_dims, "length_bin"]])
        for array in proportion_arrays
    ]
    # ---- Align and combine all length marginals
    combined_length_marginal = _sum_aligned(length_arrays)

    # Construct the originally observed and corrected length distributions
    # ---- Observed
    observed_length = _normalize_by_dims(combined_length_marginal, normalization_dims)
    # ---- Validate selectivity parameters
    valid_params = ValidateSelectivityParams.create(**selectivity_params)
    # ---- Correct
    corrected_length = _construct_corrected_distribution(
        combined_length_marginal, normalization_dims=normalization_dims, **valid_params
    )

    # When age is present, optionally build the adjusted 2D length-age distribution
    # ---- Marginalize
    combined_length_age_marginal = _compute_length_age_marginal(proportion_arrays, keep_dims)
    # ---- Apply correction, if applicable
    corrected_length_age = None
    if combined_length_age_marginal is not None:
        corrected_length_age = _construct_corrected_distribution(
            combined_length_age_marginal, normalization_dims=normalization_dims, **valid_params
        )

    # Assemble adjusted outputs
    return _build_adjusted_outputs(
        observed_len_proportion=observed_length,
        corrected_len_proportion=corrected_length,
        normalization_dims=normalization_dims,
        length_age_marginal=combined_length_age_marginal,
        corrected_length_age_proportion=corrected_length_age,
    )


def _build_aged_adjusted_outputs(
    observed_len_proportion: xr.DataArray,
    corrected_len_proportion: xr.DataArray,
    normalization_dims: List[str],
) -> xr.Dataset:
    """
    Internal utility function for formatting the output Dataset structure of the length-age-based
    selectivity-corrected proportions.
    """

    # Normalize the observed and corrected distributions
    proportions_weight_observed = _normalize_by_dims(observed_len_proportion, normalization_dims)
    proportions_weight_corrected = _normalize_by_dims(corrected_len_proportion, normalization_dims)

    # Handle the top-level length-based proportions
    if "age_bin" in proportions_weight_corrected.dims:
        length_props = {
            "length_proportions_observed": proportions_weight_observed.sum(
                dim="age_bin", skipna=True
            ),
            "length_proportions_corrected": proportions_weight_corrected.sum(
                dim="age_bin", skipna=True
            ),
            "length_age_proportions_observed": proportions_weight_observed,
            "length_age_proportions_corrected": proportions_weight_corrected,
        }
    else:
        length_props = {
            "length_proportions_observed": proportions_weight_observed,
            "length_proportions_corrected": proportions_weight_corrected,
        }
    # ---- Add downstream variable names for compatibility
    length_props["proportion"] = proportions_weight_corrected
    length_props["proportion_overall"] = proportions_weight_corrected

    # Return the Dataset
    return xr.Dataset(length_props)


def correct_weight_proportions(
    corrected_number_proportions: xr.Dataset,
    mean_length_binned_weights: xr.DataArray,
    stratum_dim: Optional[Union[str, List[str]]] = None,
) -> xr.Dataset:
    """
    Compute weight proportions from sensitivity-corrected number proportions and mean length-binned
    weights.

    Parameters
    ----------
    corrected_number_proportions : xr.Dataset
        Selectivity-corrected number proportions containing that must include the variable
        ``length_proportions_corrected`` and optionally ``length_age_proportions_corrected``.
    mean_length_binned_weights : xr.DataArray
        Mean fitted weights by ``length_bin`` based on the fitted length-weight relationship.
    stratum_dim : str or list of str, optional
        Dimensions that define the normalization groups. When provided, both observed and adjusted
        length proportions are normalized to sum to 1 within each stratum. Other inferred
        coordinates are still retained in the output arrays.

    Returns
    -------
    xr.Dataset

    """

    # Validate that 'length_proportions_corrected`, at a minimum, exists
    if not {"length_proportions_corrected", "length_proportions_observed"} <= set(
        corrected_number_proportions.data_vars
    ):
        raise KeyError(
            "Variables 'length_proportions_corrected' and 'length_proportions_observed' must be "
            "present as variables in 'corrected_number_proportions'."
        )

    # Validate that full set for length-age are together, if present
    # ---- Also select the representative proportions
    if any(
        [
            v
            for v in ["length_age_proportions_corrected", "length_age_proportions_observed"]
            if v in list(corrected_number_proportions.data_vars)
        ]
    ):
        if not {"length_age_proportions_corrected", "length_age_proportions_observed"} <= set(
            corrected_number_proportions.data_vars
        ):
            raise KeyError(
                "When provided, variables 'length_age_proportions_corrected' and "
                "'length_age_proportions_observed' must both be provided in "
                "'corrected_number_proportions'."
            )
        else:
            proportions_corrected = corrected_number_proportions["length_age_proportions_corrected"]
            proportions_observed = corrected_number_proportions["length_age_proportions_observed"]
    else:
        proportions_corrected = corrected_number_proportions["length_proportions_corrected"]
        proportions_observed = corrected_number_proportions["length_proportions_observed"]

    # Collapse the non-length dimensions and find mismatches
    if not set(mean_length_binned_weights.coords) <= set(proportions_corrected.coords):
        raise KeyError(
            "Coordinates in 'mean_length_binned_weights' must either match or be a subset of those "
            "in 'corrected_number_proportions'."
        )

    # Get the weight proportions from the number proportions
    weight_observed = (
        xr.where(proportions_observed.isnull(), 0.0, proportions_observed).astype(float)
        * mean_length_binned_weights
    )
    weight_corrected = (
        xr.where(proportions_corrected.isnull(), 0.0, proportions_corrected).astype(float)
        * mean_length_binned_weights
    )

    # Resolve retained dimensions and build combined length marginal
    keep_dims = [dim for dim in weight_corrected.dims if dim not in ["length_bin", "age_bin"]]

    # Resolve normalization dimensions
    # ---- Stratified: each stratum sums to 1
    if stratum_dim:
        # ---- Convert to list if needed
        normalization_dims = [stratum_dim] if isinstance(stratum_dim, str) else list(stratum_dim)
    # ---- Grouped: each combination of non-length, non-age coordinates sum to 1
    else:
        normalization_dims = keep_dims

    # Format and return the output Datasetthe output
    return _build_aged_adjusted_outputs(weight_observed, weight_corrected, normalization_dims)
