from typing import Union
import pandas as pd
import numpy as np

def impute_missing_sigma_bs(
    unique_strata: Union[list, np.ndarray], sigma_bs_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Imputes sigma_bs for strata without measurements or values

    Parameters
    ----------
    unique_strata : Union[list, np.ndarray]
        An array comprising all expected stratum numbers.
    sigma_bs_df : pd.DataFrame
        DataFrame containing the mean `sigma_bs` calculated for each stratum.
        Must have stratum as index and 'sigma_bs' column.

    Returns
    -------
    pd.DataFrame
        DataFrame with imputed sigma_bs values for missing strata

    Examples
    --------
    >>> # Create sigma_bs data for some strata
    >>> sigma_bs_df = pd.DataFrame({
    ...     'sigma_bs': [0.001, 0.003, 0.002]
    ... }, index=pd.Index([1, 3, 5], name='stratum'))
    >>>
    >>> # Define all expected strata
    >>> all_strata = [1, 2, 3, 4, 5]
    >>>
    >>> # Impute missing strata (2 and 4)
    >>> result = impute_missing_sigma_bs(all_strata, sigma_bs_df)
    >>> print(result)
           sigma_bs
    stratum
    1       0.001000
    2       0.002000  # interpolated between 1 and 3
    3       0.003000
    4       0.002500  # interpolated between 3 and 5
    5       0.002000

    Notes
    -----
    This function iterates through all stratum layers to impute either the
    nearest neighbor interpolation or mean sigma_bs for strata that are missing values.

    For missing strata, the function:
    1. Finds the nearest stratum below and above the missing stratum
    2. Interpolates the sigma_bs value as the mean of these neighbors
    3. If no neighbors exist on one side, uses the available neighbor value
    """

    # Extract the stratum index name
    stratum_name = sigma_bs_df.index.name

    # Collect present strata
    present_strata = np.unique(sigma_bs_df.index)

    # Invert to retrieve missing strata
    missing_strata = set(unique_strata).difference(set(present_strata))

    # Impute values for missing strata
    if len(missing_strata) > 0:

        # Concatenate the existing data with the missing strata
        sigma_bs_stratum_impute = pd.concat(
            [
                sigma_bs_df,
                pd.DataFrame({stratum_name: list(missing_strata), "sigma_bs": np.nan}).set_index(
                    stratum_name
                ),
            ],
        ).sort_values(stratum_name)

        # Find strata intervals to impute over
        for i in missing_strata:

            # Get the stratum indices below and above the missing stratum
            strata_floor = present_strata[present_strata < i]
            strata_ceil = present_strata[present_strata > i]

            new_stratum_below = (
                np.max(strata_floor) if strata_floor.size > 0 else np.min(strata_ceil)
            )
            new_stratum_above = (
                np.min(strata_ceil) if strata_ceil.size > 0 else np.max(strata_floor)
            )

            # Get the indexed values
            sigma_bs_indexed = sigma_bs_stratum_impute.loc[[new_stratum_below, new_stratum_above]]

            # Impute the mean to the missing value
            sigma_bs_stratum_impute.loc[i] = sigma_bs_indexed.mean()

        # Return the imputed values
        return sigma_bs_stratum_impute
    else:
        return sigma_bs_df