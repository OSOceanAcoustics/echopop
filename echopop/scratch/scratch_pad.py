import abc
import numpy as np
import pandas as pd
from typing import Union, Dict, List, Optional, Any

# Import the existing acoustics functions
from ..acoustics import ts_length_regression, to_linear, to_dB, impute_missing_sigma_bs


# ==============================================================================
# MY IMPLEMENTATIONS OF INVERSION CLASSES
# ==============================================================================
model_parameters = {
    "ts_length_regression": {
        "slope": 20.,
        "intercept": -68.
    },
    "stratify_by": "stratum_ks",
    "strata": df_dict_strata["ks"].stratum_num.unique(),
    "impute_missing_strata": True,
}

strata_options = np.array([1, 2, 3, 4, 5])
sigma_bs_stratum = pd.DataFrame(
    {
        "stratum_num": [1, 2, 3],
        "species_id": np.repeat(94832, 3),
        "sigma_bs_mean": [1.0, 2.0, 3.0],
    }
)


# ==============================================================================
# TRANSECT INTERVAL CORRECTION FUNCTIONS
# ==============================================================================

def correct_transect_intervals_refactored(
    nasc_data: pd.DataFrame, 
    interval_threshold: float = 0.05,
    return_columns: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    Calculate along-transect intervals and impute erroneous values for NASC data
    
    Refactored version specifically designed for DataFrames with columns:
    ['stratum_ks', 'transect_num', 'region_id', 'distance_s', 'distance_e',
     'latitude', 'longitude', 'transect_spacing', 'layer_mean_depth',
     'layer_height', 'bottom_depth', 'nasc', 'haul_num', 'group',
     'stratum_inpfc', 'stratum name', 'year', 'nasc_proportion',
     'mean length', 'geostratum_inpfc', 'geostratum_ks', 'number_density']

    Parameters
    ----------
    nasc_data : pd.DataFrame
        DataFrame containing NASC data with distance and spacing information.
        Must contain columns: 'distance_s', 'distance_e', 'transect_spacing'
    interval_threshold : float, default 0.05
        Along-transect interval threshold for detecting erroneous values.
        Values that deviate from the median interval by more than this threshold
        will be corrected using distance_e - distance_s calculation.
    return_columns : Optional[List[str]], default None
        List of specific columns to return. If None, returns all columns plus
        the calculated 'interval' and 'interval_area' columns.

    Returns
    -------
    pd.DataFrame
        DataFrame with added 'interval' and 'interval_area' columns, optionally
        filtered to specific columns if return_columns is specified.

    Examples
    --------
    >>> # Basic usage with all columns
    >>> result = correct_transect_intervals_refactored(df_nasc_all_ages)
    >>> print(result.columns)
    Index(['stratum_ks', 'transect_num', ..., 'interval', 'interval_area'])
    
    >>> # Return only specific columns
    >>> columns_to_keep = ['transect_num', 'latitude', 'longitude', 'stratum_ks', 
    ...                    'haul_num', 'nasc', 'interval_area']
    >>> result = correct_transect_intervals_refactored(
    ...     df_nasc_all_ages, return_columns=columns_to_keep
    ... )
    >>> print(result.columns)
    Index(['transect_num', 'latitude', 'longitude', 'stratum_ks', 'haul_num', 
           'nasc', 'interval_area'])

    Notes
    -----
    This function calculates the along-track transect interval length and areas.
    It identifies and corrects potentially erroneous values at transect endpoints
    by comparing intervals to the median and replacing outliers with direct
    distance calculations (distance_e - distance_s).
    
    The interval calculation uses diff(periods=-1) to compute forward differences,
    making each interval represent the distance to the next measurement point.
    
    Interval area is calculated as: interval * transect_spacing
    """
    
    # Validate required columns
    required_cols = ['distance_s', 'distance_e', 'transect_spacing']
    missing_cols = [col for col in required_cols if col not in nasc_data.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    # Create dataframe copy to avoid modifying original
    data_copy = nasc_data.copy()

    # Calculate the along-transect interval distance
    # Use forward difference to get distance to next point
    data_copy["interval"] = data_copy["distance_s"].diff(periods=-1).abs()
    
    # Handle the final interval (no next point available)
    # Use distance_e - distance_s for the last measurement
    final_interval = data_copy["distance_e"].iloc[-1] - data_copy["distance_s"].iloc[-1]
    data_copy["interval"] = data_copy["interval"].fillna(final_interval)

    # Identify and correct erroneous interval lengths
    # Calculate median interval for comparison
    median_interval = np.median(data_copy["interval"])
    
    # Find intervals that deviate significantly from median
    deviation_mask = np.abs(data_copy["interval"] - median_interval) > interval_threshold
    
    # Replace erroneous intervals with direct distance calculation
    data_copy.loc[deviation_mask, "interval"] = (
        data_copy.loc[deviation_mask, "distance_e"] - 
        data_copy.loc[deviation_mask, "distance_s"]
    )

    # Calculate the interval area (interval distance * transect width)
    data_copy["interval_area"] = data_copy["interval"] * data_copy["transect_spacing"]

    # Return specified columns or all columns
    if return_columns is not None:
        # Ensure interval and interval_area are included if not specified
        columns_to_return = list(return_columns)
        if "interval" not in columns_to_return:
            columns_to_return.append("interval")
        if "interval_area" not in columns_to_return:
            columns_to_return.append("interval_area")
        
        # Filter to only columns that exist in the DataFrame
        available_columns = [col for col in columns_to_return if col in data_copy.columns]
        return data_copy[available_columns]
    else:
        return data_copy


def get_transect_summary_stats(nasc_data: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate summary statistics for transect intervals
    
    Parameters
    ----------
    nasc_data : pd.DataFrame
        DataFrame with 'interval' and 'interval_area' columns from 
        correct_transect_intervals_refactored()
        
    Returns
    -------
    pd.DataFrame
        Summary statistics including mean, median, std, min, max for
        interval distances and areas
        
    Examples
    --------
    >>> corrected_data = correct_transect_intervals_refactored(df_nasc_all_ages)
    >>> stats = get_transect_summary_stats(corrected_data)
    >>> print(stats)
                  mean    median       std       min       max
    interval      2.45      2.50      0.15      2.10      3.20
    interval_area 4.90      5.00      0.30      4.20      6.40
    """
    
    if 'interval' not in nasc_data.columns or 'interval_area' not in nasc_data.columns:
        raise ValueError("DataFrame must contain 'interval' and 'interval_area' columns. "
                        "Run correct_transect_intervals_refactored() first.")
    
    # Calculate summary statistics
    stats = nasc_data[['interval', 'interval_area']].agg([
        'mean', 'median', 'std', 'min', 'max', 'count'
    ]).round(4)
    
    return stats.T  # Transpose for better readability