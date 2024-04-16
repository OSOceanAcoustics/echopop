import numpy as np
import scipy.stats as st


def stratified_transect_statistic(
    transects, strata, sample_fraction, num_replicates, parameter: str
):
    """
    Calculates stratified mean statistics for a set of transects


    Parameters
    ----------
    transects: pd.DataFrame
        DataFrame comprising a variety of spatial metrics for transect data
    strata: pd.DataFrame
        DataFrame comprising summary features of latitude (INPFC) delimited strata
    sample_fraction: np.float64
        Value representing the proportion of transects that are resampled from the
        overall dataset within each strata
    num_replicates: int
        The number of iterations/realizations used for bootstrapping
    parameter: str
        Parameter-of-interest that will be summarized (e.g. biomass)

    Notes
    -----
    This function calculates the stratified summary statistics for biomass within
    `echopop.survey.stratified_summary()`.
    """

    # Convert specific DataFrame columns to arrays for speed
    distance = transects["transect_distance"].values
    value = transects[parameter].values
    num_transects = strata["num_transects"].values
    total_transect_area = strata.set_index("stratum_inpfc")["total_transect_area"]

    # Calculate the number of transects within each stratum based on the
    # sampling faction defined from the configuration file
    # ---- Number of transects
    num_transects_to_sample = np.round(sample_fraction * num_transects).astype(int)

    # ---- Offset term used for later variance calculation
    sample_offset = np.where(num_transects_to_sample == 1, 0, 1)

    # ---- Calculate effective sample size/degrees of freedom for variance calculation
    sample_dof = num_transects_to_sample * (num_transects_to_sample - sample_offset)

    # Pre-allocate and pre-compute the cumulative sum of numbered transects per strata
    # ---- Transect indices
    cum_num_transects = np.concatenate(([0], np.cumsum(num_transects)))

    # ---- Stratified statistics
    mean_arr = np.empty(num_replicates)
    variance_arr = np.empty(num_replicates)

    # Iterate across all replicate iterations/realizations
    for i in range(num_replicates):

        # Pre-allocate the stratum-specific means and variances
        rho_j = np.empty_like(total_transect_area)  # mean
        var_j = np.empty_like(total_transect_area)  # variance

        # Iterate across all strata
        for j in strata.stratum_inpfc - 1:

            # Resample (without replacement) based on binned indices
            # ---- Define start and end transects within each stratum
            start, end = cum_num_transects[j], cum_num_transects[j + 1]

            # ---- Resample without replacement
            sel_inds = np.random.choice(
                np.arange(start, end), num_transects_to_sample[j], replace=False
            )

            # Define stratified weights
            stratified_weights = distance[sel_inds] / np.mean(distance[sel_inds])

            # Weighted value (e.g. biomass)
            value_distance_density = value[sel_inds] / distance[sel_inds]

            # Compute mean and variance
            rho_j[j] = np.nansum(value[sel_inds] * stratified_weights) / np.nansum(
                stratified_weights
            )
            var_j[j] = (
                np.nansum((stratified_weights**2 * (value_distance_density - rho_j[j]) ** 2))
                / sample_dof[j]
            )

        # Calculate the overall weighted means and variances for later calculations
        # ---- Mean
        mean_arr[i] = np.nansum(strata.total_transect_area * rho_j)

        # ---- Variance
        variance_arr[i] = np.sqrt(np.nansum(var_j * strata.total_transect_area**2))

    # Calculate the summary statistics
    stratified_results = {
        "biomass": {
            "mean": {
                "estimate": np.mean(mean_arr),
                "confidence_interval": confidence_interval(mean_arr),
            },
            "variance": {
                "estimate": np.mean(variance_arr),
                "confidence_interval": confidence_interval(variance_arr),
            },
            "CV": {
                "estimate": np.mean(variance_arr / mean_arr),
                "confidence_interval": confidence_interval(variance_arr / mean_arr),
            },
        }
    }

    # Return output
    return stratified_results


def confidence_interval(values):
    """
    Calculates the 95% confidence interval (Normal) for a given array

    Parameters
    ----------
    values: np.array
        An array of values

    Notes
    -----
    This function calculates the 95% confidence interval (assuming a Normal) distribution
    for the bootstrapped stratified samples. This is done as opposed to using the percentile
    method for estimate the intervals.
    """
    return np.mean(values) + np.array([-1, 1]) * st.norm.ppf(0.975) * np.std(values)
