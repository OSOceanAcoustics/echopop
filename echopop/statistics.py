import warnings
from typing import Literal, Optional, Union

import numpy as np
import pandas as pd
import scipy.stats as st

from .spatial.transect import transect_array


def stratified_transect_statistic(
    transect_data: pd.DataFrame,
    transect_summary: pd.DataFrame,
    strata_summary: pd.DataFrame,
    settings_dict: dict,
) -> tuple[pd.DataFrame, dict]:
    """
    Calculates stratified mean statistics for a set of transects


    Parameters
    ----------
    transect_data: pd.DataFrame
        Dataframe comprising georeferenced biological data collected from survey transects
    transect_summary: pd.DataFrame
        DataFrame comprising a variety of spatial metrics for transect data
    strata_summary: pd.DataFrame
        DataFrame comprising summary features of latitude (INPFC) delimited strata
    settings_dict: dict
        Dictionary containing algorithm arguments that define the proportion of transects resampled
        from the overall dataset within each strata (`transect_sample`) and the number of
        iterations/realizations used for bootstrapping/resampling (`transect_replicates`).

    Notes
    -----
    This function calculates the stratified summary statistics for biomass within
    `echopop.survey.stratified_summary()`.
    """

    # Extract algorithm arguments
    # ---- Number of replicates
    transect_replicates = settings_dict["transect_replicates"]
    # ---- Transect sampling fraction
    transect_sample = settings_dict["transect_sample"]
    # ---- Get stratum column name
    stratum_col = settings_dict["stratum_name"]
    # ---- Get the variable name
    var_name = settings_dict["variable"]

    # Get indexed transect distance
    transect_distances = transect_summary.set_index(["transect_num"])["transect_distance"]
    # ---- Drop any transects where distance is 0.0 (i.e. from a single mesh node)
    if np.any(transect_distances == 0.0):
        # ---- Pick out transects where distance = 0.0 nmi
        zero_distances = transect_distances[transect_distances == 0.0].index.to_numpy()
        # ---- Update `transect_distances`
        transect_distances = transect_distances[transect_distances > 0.0]
        # ---- Update `transect_data`
        transect_data = transect_data[~transect_data["transect_num"].isin(zero_distances)]
        # ---- Get the 'poor' transect strata
        zero_distances_strata = (
            transect_summary.loc[zero_distances].groupby([stratum_col], observed=False).size()
        )
        # ---- Update `transect_summary`
        transect_summary = transect_summary[~transect_summary["transect_num"].isin(zero_distances)]
        # ---- Update `strata_summary`
        # -------- Set index
        strata_summary.set_index([stratum_col], inplace=True)
        # -------- Subtract the 'poor' transects from the total transect counts
        strata_summary["transect_count"] = strata_summary["transect_count"] - zero_distances_strata
        # -------- Reset index
        strata_summary.reset_index(inplace=True)

        if settings_dict["verbose"]:
            if settings_dict["dataset"] == "kriging":
                # ---- Generate warning
                warnings.warn(
                    "Invalid virtual transects detected with distances of 0.0 nmi:"
                    f"\n[virtual transect numbers] -> {zero_distances}."
                    "\nThese will be removed from the stratified analysis.",
                    stacklevel=1,
                )
            else:
                # ---- Generate warning
                warnings.warn(
                    "Invalid survey transects detected with distances of 0.0 nmi:"
                    f"\n[virtual transect numbers] -> {zero_distances}."
                    "\nThese will be removed from the stratified analysis.",
                    stacklevel=1,
                )

    # Calculate the number of transects per stratum
    num_transects_to_sample = np.round(
        strata_summary.set_index(stratum_col)["transect_count"] * transect_sample
    ).astype(int)

    # Offset term used for later variance calculation
    sample_offset = np.where(num_transects_to_sample == 1, 0, 1)

    # Calculate effective sample size/degrees of freedom for variance calculation
    sample_dof = num_transects_to_sample * (num_transects_to_sample - sample_offset)

    # Transect areas
    transect_areas = transect_summary.groupby(["transect_num"])["transect_area"].sum()

    # Get indexed total transect area
    total_transect_area = strata_summary.set_index(stratum_col)["transect_area_total"]

    # Get indexed biological value
    biological_values = transect_data.groupby(["transect_num"])[var_name].sum()

    # Get indexed transect numbers
    transect_numbers = transect_summary.set_index(stratum_col)["transect_num"]

    # Calculate the mean density
    transect_summary["density"] = transect_data.groupby(["transect_num"])[
        settings_dict["variable"]
    ].sum()

    # Pre-allocate the stratum-specific metrics
    # ---- Mean
    mean_arr = np.zeros([transect_replicates, len(total_transect_area.index)])
    # ---- Variance
    variance_arr = np.zeros_like(mean_arr)
    # ---- Transect length
    length_arr = np.zeros_like(mean_arr)
    # ---- Transect area
    area_arr = np.zeros_like(mean_arr)
    # ---- Sum/integrated total across full stratified area/region
    total_arr = np.zeros_like(mean_arr)

    # Iterate across all iterations/realizations
    for j in total_transect_area.index:

        # Create an index array/matrix containing resampled (without replacement) transect numbers
        transect_numbers_arr = np.array(
            [
                np.random.choice(
                    transect_numbers[j].values, num_transects_to_sample[j], replace=False
                )
                for i in range(transect_replicates)
            ]
        )

        # Assign the indexed transect distance and biological variables to each transect
        # ---- Transect lengths
        distance_replicates = np.apply_along_axis(
            transect_array, 0, transect_numbers_arr, transect_distances
        )
        # -------- Calculate the summed transect length
        length_arr[:, j - 1] = distance_replicates.sum(axis=1)
        # ---- Transect areas
        area_replicates = np.apply_along_axis(
            transect_array, 0, transect_numbers_arr, transect_areas
        )
        # -------- Calculate the summed transect length
        area_arr[:, j - 1] = area_replicates.sum(axis=1)
        # ---- Biological variable
        biology_replicates = np.apply_along_axis(
            transect_array, 0, transect_numbers_arr, biological_values
        )

        # Calculate the stratified weights for along-transect values (within transect weights)
        stratified_weights = (
            distance_replicates.transpose() / distance_replicates.mean(axis=1)
        ).transpose()

        # Standardize the biological values by their respective distances
        biology_adjusted = biology_replicates / distance_replicates

        # Calculate the mean transect-length-weighted biological values
        mean_arr[:, j - 1] = (biology_replicates * distance_replicates).sum(axis=1) / length_arr[
            :, j - 1
        ]

        # Sum the total of the biology variable
        total_arr[:, j - 1] = biology_replicates.sum(axis=1)

        # Calculate the variance of the transect-length-weighted biological values
        # ---- Calculate the sqauared deviation of the mean
        squared_deviation = (biology_adjusted.transpose() - mean_arr[:, j - 1]).transpose() ** 2
        # ---- Sum of all weighted squared deviations
        squared_deviation_wgt = (stratified_weights**2 * squared_deviation).sum(axis=1)
        # ---- Compute the variance by incorporating the degrees of freedom
        variance_arr[:, j - 1] = squared_deviation_wgt / sample_dof[j]

    # Compute summary statistics first
    # ---- Convert transect area to an array
    area_array = total_transect_area.to_numpy()
    # ---- Sum the total area
    total_area = area_array.sum()

    # Compute the "population" (i.e. original data) statistics
    # This is necessary for constructing the bootstrapped confidence intervals
    # ---- Mean density
    if settings_dict["variable"] == "nasc":
        # ---- Compute sum per transect line first
        line_density = transect_data.groupby(["transect_num"])[var_name].sum().to_frame()
        # ---- Create copy of `transect_summary` and set index
        line_length = transect_summary.copy().set_index("transect_num")
        # ---- Add stratum
        line_density[stratum_col] = line_length[stratum_col]
        # ---- Convert to the density
        line_density["density"] = line_density[var_name] / line_length["transect_distance"]
        # ---- Reset index
        line_density.reset_index(inplace=True)
        # ---- Calculate mean per stratum
        stratum_density_means = (
            line_density.groupby([stratum_col])["density"].mean().to_numpy().flatten()
        )
        # ---- Calculate mean per survey
        survey_density_mean = stratum_density_means.mean()
    else:
        # ---- Get density column name
        density_name = [col for col in transect_data.columns if "_density" in col]
        # ---- Calculate mean per stratum
        stratum_density_means = (
            transect_data.groupby([stratum_col], observed=False)[density_name]
            .mean()
            .to_numpy()
            .flatten()
        )
        # ---- Calculate mean per survey
        survey_density_mean = stratum_density_means.mean()
    # ---- Total
    # -------- By stratum
    stratum_total = transect_data.groupby([stratum_col], observed=False)[var_name].sum().to_numpy()
    # -------- By survey
    survey_total = stratum_total.sum()
    # ---- Compute the stratum total proportions relative to survey sum
    stratum_proportions = stratum_total / survey_total

    # Create population/survey dictionary
    survey_dict = {
        "density": {"stratum": stratum_density_means, "survey": survey_density_mean},
        "total": {"stratum": stratum_total, "survey": survey_total},
        "proportions": stratum_proportions,
    }

    # Reverse-engineer the density and total stratified estimates from bootstrapped results
    # ---- By stratum (density)
    unweighted_stratum_density = mean_arr / length_arr
    # ---- By stratum (total)
    unweighted_stratum_total = unweighted_stratum_density * area_array
    # ---- By survey (total)
    unweighted_survey_total = unweighted_stratum_total.sum(axis=1)
    # ---- By survey (density)
    unweighted_survey_density = unweighted_survey_total / total_area

    # Reverse-engineer the proportional stratum distributions
    unweighted_stratum_proportions = total_arr / total_arr.sum(axis=1).reshape(-1, 1)

    # Compute the transect-length weighted coefficient of variation (CV)
    # ---- Compute the variance
    weighted_variance = (variance_arr * area_array**2).sum(axis=1)
    # ---- Convert to the standard deviation
    weighted_stdev = np.sqrt(weighted_variance)
    # ---- Compute the mean
    weighted_mean = (mean_arr * area_array).sum(axis=1)
    # ---- Compute the CV
    bootstrap_cv = weighted_stdev / weighted_mean

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
        "num_transects": strata_summary["transect_count"].sum(),
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


def empirical_ci(bootstrap_samples: np.ndarray, ci_args: dict):
    """
    Empirical bootstrap interval.
    """
    # Extract the interval
    interval = ci_args["interval"]

    # Compute the empirical/basic CI
    # ---- Compute the delta deviation term
    delta = bootstrap_samples - ci_args["population_statistic"]
    # ---- Compute the percentiles
    percentiles = np.quantile(delta, interval)
    # ---- Return the output
    return bootstrap_samples.mean() + percentiles


def percentile_ci(bootstrap_samples: np.ndarray, ci_args: dict):
    """
    Percentile bootstrap interval.
    """
    # Extract the interval
    interval = ci_args["interval"]

    # Return the exact percentiles of the bootstrapped distribution
    return np.quantile(bootstrap_samples, interval)


def standard_ci(bootstrap_samples: np.ndarray, ci_args: dict):
    """
    Normal/standard bootstrap interval.
    """
    # Extract the interval
    interval = ci_args["interval"]

    # Return the parametric confidence interval (assuming a Normal distribution)
    return bootstrap_samples.mean() + st.norm.ppf(interval) * bootstrap_samples.std(ddof=1)


def student_ci(bootstrap_samples: np.ndarray, ci_args: dict):
    """
    Studentized bootstrap interval assuming either a t-distribution or using a jackknife approach.
    """
    # Extract the interval
    interval = ci_args["interval"]

    # Compute the studentized confidence interval
    if ci_args["method"] == "standard":
        # ---- Return the parametric confidence interval (assuming a t-distribution)
        return bootstrap_samples.mean() + st.t.ppf(
            interval, len(bootstrap_samples) - 1
        ) * bootstrap_samples.std(ddof=1)
    else:
        # ---- Use the jackknife approach to confer additional robustness
        # -------- Pre-initialize the t-statistic list
        t_statistic = []
        # -------- Iterate through the bootstrapped sample using a jackknife approach
        for i in range(len(bootstrap_samples)):
            # -------- Create a pseudo-sample via a leave-one-out (jackknife) approach
            pseudo_sample = np.delete(bootstrap_samples, i)
            # -------- Calculate the pseudo-statistic
            pseudo_statistic = pseudo_sample.mean()
            # -------- Bootstrapped standard error from the pseudo-sample
            pseudo_se = pseudo_sample.std(ddof=1)
            # -------- Compute the t-statistic (using the left-out value)
            t_statistic.append((bootstrap_samples[i] - pseudo_statistic) / pseudo_se)
        # ---- Convert the list to an array
        t_array = np.array(t_statistic)
        # ---- Extract the t-distribution intervals
        t_interval = np.quantile(t_array, interval)
        # ---- Return the jackknifed (non-parameteric) studentized confidence interval
        return bootstrap_samples.mean() + t_interval * bootstrap_samples.std(ddof=1)


def bc_ci(bootstrap_samples: np.ndarray, ci_args: dict):
    """
    Bias-corrected (BC) bootstrap interval.
    """
    # Extract the interval
    interval = ci_args["interval"]

    # Calculate the initial Z-statistic observed in the population data
    # ---- Compute the ECDF and proportion of values below population statistic
    ecdf = np.mean(bootstrap_samples < ci_args["population_statistic"])
    # -------- Navigate edge cases but otherwise proceed to computing the Z-statistic
    if ecdf in [0.0, 1.0]:
        if "alternative_approach" in ci_args.keys():
            # ---- Pass new dictionary instructions
            temp_ci_args = ci_args.copy()
            # ---- Add new method argument, if necessary
            if "alternative_method" in ci_args.keys():
                temp_ci_args.update({"method": ci_args["alternative_method"]})
                # ---- Run the alternative model
                results = ci_args["alternative_approach"](bootstrap_samples, temp_ci_args)
        else:
            # -------- Generate Error string that will be iteratively constructed
            # -------- If samples extracted from 2D array
            if "group_name" in ci_args.keys():
                z0_error = (
                    f"Significant skewness detected among bootstrapped {ci_args['estimator_name']} "
                    f"in group index {ci_args['group_name']} when calculating the bootstrapped "
                    f"confidence intervals via {ci_args['boot_ci_method']}! "
                )
            # -------- If samples extracted from 1D array
            else:
                z0_error = (
                    f"Significant skewness detected among bootstrapped {ci_args['estimator_name']} "
                    f"when calculating the bootstrapped confidence intervals via "
                    f"{ci_args['boot_method']}!"
                )
            # -------- Finish Error message
            z0_error += (
                f"The estimated Z-statistic (`z0`) was either degenerate or too skewed to produce "
                f"reasonable interval estimates. Consider using a different interval approximation "
                f"method (`boot_ci_method={ci_args['boot_ci_method']}) or provide an alternative "
                f"(`boot_ci_method_alt`)."
            )
            # -------- Print
            raise ValueError(z0_error)
    else:
        # ---- Compute Z-statistic
        z0 = st.norm.ppf(ecdf)

        # Estimate the bias-corrected percentiles from the Normal cumulative density function (CDF)
        percentiles_cdf = st.norm.cdf(2 * z0 + st.norm.ppf(interval))
        # ---- Compute the adjusted confidence interval
        results = np.quantile(bootstrap_samples, percentiles_cdf)

    # Return results
    return results


def bca_ci(bootstrap_samples: np.ndarray, ci_args: dict):
    """
    Bias-Corrected and Accelerated (BCa) bootstrap interval.

    Notes
    ----------
    The acceleration constant (a) is estimated using jackknife resampling.
    """
    # Extract the interval
    interval = ci_args["interval"]

    # Calculate the initial Z-statistic observed in the population data
    # ---- Compute the ECDF and proportion of values below population statistic
    ecdf = np.mean(bootstrap_samples < ci_args["population_statistic"])
    # -------- Navigate edge cases but otherwise proceed to computing the Z-statistic
    if ecdf in [0.0, 1.0]:
        if "alternative_approach" in ci_args.keys():
            # ---- Pass new dictionary instructions
            temp_ci_args = ci_args.copy()
            # ---- Add new method argument, if necessary
            if "alternative_method" in ci_args.keys():
                temp_ci_args.update({"method": ci_args["alternative_method"]})
                # ---- Run the alternative model
                results = ci_args["alternative_approach"](bootstrap_samples, temp_ci_args)
        else:
            # -------- Generate Error string that will be iteratively constructed
            # -------- If samples extracted from 2D array
            if "group_name" in ci_args.keys():
                z0_error = (
                    f"Significant skewness detected among bootstrapped {ci_args['estimator_name']} "
                    f"in group index {ci_args['group_name']} when calculating the bootstrapped "
                    f"confidence intervals via {ci_args['boot_ci_method']}! "
                )
            # -------- If samples extracted from 1D array
            else:
                z0_error = (
                    f"Significant skewness detected among bootstrapped {ci_args['estimator_name']} "
                    f"when calculating the bootstrapped confidence intervals via "
                    f"{ci_args['boot_method']}!"
                )
            # -------- Finish Error message
            z0_error += (
                f"The estimated Z-statistic (`z0`) was either degenerate or too skewed to produce "
                f"reasonable interval estimates. Consider using a different interval approximation "
                f"method (`boot_ci_method={ci_args['boot_ci_method']}) or provide an alternative "
                f"(`boot_ci_method_alt`)."
            )
            # -------- Print
            raise ValueError(z0_error)
    else:
        # ---- Compute Z-statistic
        z0 = st.norm.ppf(ecdf)

        # Approximate the acceleration constant using a jackknife approach
        # ---- Define jackknife samples
        jackknife_samples = np.array(
            [np.delete(bootstrap_samples, i) for i in range(len(bootstrap_samples))]
        )
        # ---- Calculate the mean of each sample
        jackknife_means = jackknife_samples.mean(axis=1)
        # ---- Calculate sample-mean deviations
        deviations = jackknife_means.mean() - jackknife_means
        # ---- Estimate the acceleration constant (`a_hat``)
        a_hat = (deviations**3).sum() / (6.0 * (deviations**2).sum() ** 1.5)

        # Estimate the bias-corrected percentiles from the Normal cumulative density function (CDF)
        percentiles_cdf = st.norm.cdf(
            z0 + (z0 + st.norm.ppf(interval)) / (1 - a_hat * (z0 + st.norm.ppf(interval)))
        )

        # Compute the adjusted confidence interval
        results = np.quantile(bootstrap_samples, percentiles_cdf)

    # Return results
    return results


def bootstrap_confidence_intervals(
    bootstrap_dict: dict, population_dict: dict, settings_dict: dict
) -> dict:
    """
    Compute the bootstrap confidence interval.

    Notes
    ----------
    Percentile method[1]_
    t-jackknife[4]_
    t-standard[4]_
    Empirical[1]_
    Bias-Corrected and Accelerated (BCa)[2]_
    Bias-Corrected (BC)[1]_
    Standard[3]_

    .. [1] Efron, B. (1981). Nonparametric standard errors and confidence intervals. *Canadian
       Journal of Statistics*, *9*(2), 139-158. https://doi.org/10.2307/3314608
    .. [2] Efron, B., and Tibshirani, R.J. (1993). *An introduction to the Bootstrap*. Springer US.
        https://doi.org/10.1007/978-1-4899-4541-9
    .. [3] Efron, B., and Tibshirani, R.J. (1986). Bootstrap methods for standard errors, confidence
       intervals, and other measures of statistical accuracy. *Statistical Science*, *1*(1).
       https://doi.org/10.1214/ss/1177013815
    .. [4] DiCiccio, T.J., and Efron, B. (1996). Bootstrap confidence intervals. *Statistical
       Science*, *11*(3). https://doi.org/10.1214/ss/1032280214
    """
    # Extract variable name
    var_name = settings_dict["variable"]
    # Compute the confidence interval (CI) and bias estimates for density measurements
    # ---- Stratum
    stratum_density_cis, stratum_density_bias = confidence_interval(
        bootstrap_dict["density"]["stratum"],
        population_dict["density"]["stratum"],
        settings_dict["bootstrap_ci"],
        settings_dict["bootstrap_ci_method"],
        settings_dict["bootstrap_ci_method_alt"],
        settings_dict["bootstrap_adjust_bias"],
        f"STRATUM {var_name.upper()} DENSITY",
    )
    # ---- Survey
    survey_density_cis, survey_density_bias = confidence_interval(
        bootstrap_dict["density"]["survey"],
        population_dict["density"]["survey"],
        settings_dict["bootstrap_ci"],
        settings_dict["bootstrap_ci_method"],
        settings_dict["bootstrap_ci_method_alt"],
        settings_dict["bootstrap_adjust_bias"],
        f"SURVEY {var_name.upper()} DENSITY",
    )

    # Compute the confidence interval (CI) and bias estimates for total measurements
    # ---- Stratum
    stratum_total_cis, stratum_total_bias = confidence_interval(
        bootstrap_dict["total"]["stratum"],
        population_dict["total"]["stratum"],
        settings_dict["bootstrap_ci"],
        settings_dict["bootstrap_ci_method"],
        settings_dict["bootstrap_ci_method_alt"],
        settings_dict["bootstrap_adjust_bias"],
        f"STRATUM {var_name.upper()} TOTAL",
    )
    # ---- Survey
    survey_total_cis, survey_total_bias = confidence_interval(
        bootstrap_dict["total"]["survey"],
        population_dict["total"]["survey"],
        settings_dict["bootstrap_ci"],
        settings_dict["bootstrap_ci_method"],
        settings_dict["bootstrap_ci_method_alt"],
        settings_dict["bootstrap_adjust_bias"],
        f"SURVEY {var_name.upper()} TOTAL",
    )

    # Compute the confidence interval (CI) and bias estimates for proportional estimates
    # ---- Proportions (stratum)
    stratum_proportion_cis, stratum_proportion_bias = confidence_interval(
        bootstrap_dict["proportions"],
        population_dict["proportions"],
        settings_dict["bootstrap_ci"],
        settings_dict["bootstrap_ci_method"],
        settings_dict["bootstrap_ci_method_alt"],
        settings_dict["bootstrap_adjust_bias"],
        f"STRATUM {var_name.upper()} PROPORTIONS",
    )
    # ---- CV
    survey_cv_cis, survey_cv_bias = confidence_interval(
        bootstrap_dict["cv"],
        bootstrap_dict["cv"].mean(),
        settings_dict["bootstrap_ci"],
        settings_dict["bootstrap_ci_method"],
        settings_dict["bootstrap_ci_method_alt"],
        settings_dict["bootstrap_adjust_bias"],
        "SURVEY CV",
    )

    # Organize the output results
    ci_results = {
        "ci": {
            "strata": {
                "density": stratum_density_cis,
                "total": stratum_total_cis,
                "proportion": stratum_proportion_cis,
            },
            "survey": {
                "density": survey_density_cis,
                "total": survey_total_cis,
                "cv": survey_cv_cis,
            },
        },
        "bias": {
            "strata": {
                "density": stratum_density_bias,
                "total": stratum_total_bias,
                "proportion": stratum_proportion_bias,
            },
            "survey": {
                "density": survey_density_bias,
                "total": survey_total_bias,
                "cv": survey_cv_bias,
            },
        },
    }

    # Return the results
    return ci_results


# Bootstrap confidence interval methods API
BOOTSTRAP_CI_METHODS = {
    "single_method": {
        "bc": bc_ci,
        "bca": bca_ci,
        "empirical": empirical_ci,
        "percentile": percentile_ci,
        "standard": standard_ci,
    },
    "multimethod": {"t": student_ci},
}


def confidence_interval(
    bootstrap_samples: np.ndarray[float],
    population_statistic: np.ndarray[float],
    ci_percentile: float = 0.95,
    boot_ci_method: Literal[
        "BC", "BCa", "empirical", "percentile", "standard", "t-jackknife", "t-standard"
    ] = "BCa",
    boot_ci_method_alt: Optional[
        Literal["empirical", "percentile", "standard", "t-jackknife", "t-standard"]
    ] = None,
    adjust_bias: bool = True,
    estimator_name: Optional[str] = None,
) -> tuple[Union[np.ndarray, list[np.ndarray]], Union[float, np.ndarray[float]]]:
    """
    Calculates the confidence interval for a given array using bootstrapped methods.

    Parameters
    ----------
    bootstrap_samples: np.ndarray
        Bootstrap replicates resampled from stratified data.
    population_statistic: np.ndarray
        An array of values representing the full dataset assuming no stratification.
    ci_percentile: float
        Confidence interval percentile.
    boot_ci_method: Literal["BC", "BCa", "empirical", "percentile", "standard", "t-jackknife",
                            "t-standard"]
        Method for computing the bootstrap interval.
    boot_ci_method_alt: Literal["empirical", "percentile", "standard", "t-jackknife", "t-standard"]
        An optional argument that states the "back-up" method to use for computing the bootstrap
        intervals in case of distributions not being compatible with specific algorithms (e.g. BCa).
    adjust_bias: bool
        Flags whether to de-bias the CI estimates.
    estimator_name: str
        An optional string that is for error-tracing purposes.

    Notes
    ----------
    Percentile method[1]_
    t-jackknife[4]_
    t-standard[4]_
    Empirical[1]_
    Bias-Corrected and Accelerated (BCa)[2]_
    Bias-Corrected (BC)[1]_
    Standard[3]_

    .. [1] Efron, B. (1981). Nonparametric standard errors and confidence intervals. *Canadian
       Journal of Statistics*, *9*(2), 139-158. https://doi.org/10.2307/3314608
    .. [2] Efron, B., and Tibshirani, R.J. (1993). *An introduction to the Bootstrap*. Springer US.
        https://doi.org/10.1007/978-1-4899-4541-9
    .. [3] Efron, B., and Tibshirani, R.J. (1986). Bootstrap methods for standard errors, confidence
       intervals, and other measures of statistical accuracy. *Statistical Science*, *1*(1).
       https://doi.org/10.1214/ss/1177013815
    .. [4] DiCiccio, T.J., and Efron, B. (1996). Bootstrap confidence intervals. *Statistical
       Science*, *11*(3). https://doi.org/10.1214/ss/1032280214
    """

    # Compute the interval percentiles
    interval = np.array([(1 - ci_percentile) / 2, (1 + ci_percentile) / 2])

    # Bundle the CI arguments
    ci_args = {
        "population_statistic": population_statistic,
        "interval": interval,
        "boot_ci_method": boot_ci_method,
        "estimator_name": (estimator_name if estimator_name is not None else "Unnamed variable"),
        "adjust_bias": adjust_bias,
    }

    # Convert to lowercase and parse reference API for the correct method
    # ---- Split the `boot_ci_method`
    ci_method = boot_ci_method.split("-")
    # ---- Parse the reference dictionary for the correct model
    # -------- Approximations with single methods
    if len(ci_method) == 1:
        # -------- Amend the ci_method string
        ci_method_name = "".join(ci_method).lower()
        # -------- Search for matching function
        ci_method_function = BOOTSTRAP_CI_METHODS["single_method"][ci_method_name]
    # -------- Approximations with multiple methods
    else:
        # -------- Amend the ci_method string
        ci_method_name = ci_method[0].lower()
        # -------- Population the ci_args list
        ci_args.update({"method": ci_method[1]})
        # -------- Search for matching function
        ci_method_function = BOOTSTRAP_CI_METHODS["multimethod"][ci_method_name]

    # Define the alternative method if parameterized
    if boot_ci_method_alt is not None:
        # ---- Split the alternative `boot_ci_method` if present
        ci_method_alt = boot_ci_method_alt.split("-")
        # ---- Cannot be the same as the primary model!
        if np.all(ci_method_alt == ci_method):
            warnings.warn(
                "Identical primary and alternative bootstrap interval approximation methods "
                "selected. Alternative method will be removed to avoid recursion."
            )
        else:
            if len(ci_method_alt) == 1:
                # ---- Amend the ci_method_alt string
                ci_method_alt_name = "".join(ci_method_alt).lower()
                # ---- Search for matching function
                ci_method_alt_function = BOOTSTRAP_CI_METHODS["single_method"][ci_method_alt_name]
            else:
                # -------- Amend the ci_method string
                ci_method_alt_name = ci_method_alt[0].lower()
                # -------- Population the ci_args list
                ci_args.update({"alternative_method": ci_method_alt[1]})
                # -------- Search for matching function
                ci_method_alt_function = BOOTSTRAP_CI_METHODS["multimethod"][ci_method_alt_name]
            # ---- Add alternative model to ci_args
            ci_args.update({"alternative_approach": ci_method_alt_function})

    # Check the dimensions of `bootstrap_samples` to determine a 1D- or 2D-array approach
    # ---- 1D
    if bootstrap_samples.ndim == 1:
        # ---- Bias
        bias = bootstrap_samples.mean() - population_statistic
        # ---- Subtract bias, if necessary
        if adjust_bias:
            bootstrap_samples_adj = bootstrap_samples - bias
            # ---- Confidence interval
            ci_estimate = ci_method_function(bootstrap_samples_adj, ci_args)
        else:
            # ---- Confidence interval
            ci_estimate = ci_method_function(bootstrap_samples, ci_args)
    # ---- 2D (this assumes that the shape has replicates distributed across rows)
    elif bootstrap_samples.ndim == 2:
        # ---- Bias
        bias = bootstrap_samples.mean(axis=0) - population_statistic
        # ---- Confidence interval
        # -------- Pre-allocate a list
        ci_estimate = []
        for i in range(bootstrap_samples.shape[1]):
            # -------- Push an indexed population statistic
            ci_args["population_statistic"] = population_statistic[i]
            # -------- Push an indexed group value
            ci_args.update({"group_name": i})
            # ---- Subtract bias, if necessary
            if adjust_bias:
                bootstrap_samples_adj = bootstrap_samples[:, i] - bias[i]
                # -------- Compute the CI
                ci_estimate.append(ci_method_function(bootstrap_samples_adj, ci_args))
            else:
                # -------- Compute the CI
                ci_estimate.append(ci_method_function(bootstrap_samples[:, i], ci_args))
    # ---- Invalid
    else:
        raise ValueError("Bootstrapped data must be a 1D array or 2D matrix!")

    # Return the (corrected/adjusted) CI and bias estimates
    return ci_estimate, bias
