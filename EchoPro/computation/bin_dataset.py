"""
Constructs an xarray Dataset filled with binned parameters
from survey data which are useful for downstream processes,
such as calculating variables over the Kriging mesh points,
length-age defined variables, and creating variables for reports.
"""

from typing import List, Tuple

import numpy as np
import pandas as pd
import xarray as xr


def _get_bin_ind(input_data: np.ndarray, bin_edges: np.ndarray) -> List[np.ndarray]:
    """
    This function manually finds those indices in ``input_data``
    that are in each bin.

    Parameters
    ----------
    input_data: np.ndarray
        The data to bin
    bin_edges: np.ndarray
        An array that specifies the bin edges.

    Returns
    -------
    hist_ind: list
        The index values of ``input_data`` corresponding to the histogram

    Notes
    -----
    The construction of the bin counts differs from the method
    produced in `transect_results.py`. This is because the Matlab
    version of EchoPro is inconsistent in how it is binning.
    """

    # initialize list that will hold the indices
    hist_ind = []

    for i in range(len(bin_edges) - 1):

        # get values greater or equal than lower bound
        g_lb = bin_edges[i] <= input_data

        # get values less than the upper bound
        le_ub = input_data < bin_edges[i + 1]

        # fill bin
        hist_ind.append(np.argwhere(g_lb & le_ub).flatten())

    # fill in the last bin
    hist_ind.append(np.argwhere(input_data >= bin_edges[-1]).flatten())

    return hist_ind


def _get_bin_ind_age(input_data: np.ndarray, age_bins: np.ndarray) -> List[np.ndarray]:
    """
    This function manually finds the indices of ``input_data``
    for each age bin. An age bin here is described as all values
    equal to the provided age

    Parameters
    ----------
    input_data: np.ndarray
        The data to bin
    age_bins: np.ndarray
        An array that specifies the age for each bin

    Returns
    -------
    hist_ind: list
        The index values of ``input_data`` corresponding to the histogram

    Notes
    -----
    The construction of the bin counts differs from the method
    produced in `transect_results.py`. This is because the Matlab
    version of EchoPro is inconsistent in how it is binning ages.
    """

    # initialize list that will hold indices
    hist_ind = []

    for age in age_bins:

        # determine indices with a value equal to age
        hist_ind.append(np.argwhere(input_data == age).flatten())

    return hist_ind


def _initialize_ds(
    stratum_ind: np.ndarray, len_bin: np.ndarray, age_bin: np.ndarray
) -> xr.Dataset:
    """
    Initializes the parameter Dataset.

    Parameters
    ----------
    stratum_ind: np.ndarray
        A one dimensional array specifying all strata values
    len_bin: np.ndarray
        A one dimensional array specifying the length bin values
    age_bin: np.ndarray
        A one dimensional array specifying all age bin values

    Returns
    -------
    ds: xr.Dataset
        The parameter Dataset with all variables initialized
    """

    # initialize variable that will hold all Dataset initialized variables
    data_vars_dict = {
        "total_weight": ("stratum_num", np.zeros(len(stratum_ind))),
        "aged_proportion": ("stratum_num", np.zeros(len(stratum_ind))),
        "unaged_proportion": ("stratum_num", np.zeros(len(stratum_ind))),
        "station_1_N": ("stratum_num", np.zeros(len(stratum_ind))),
        "weight_len_all_normalized": (
            ["stratum_num", "len_bin"],
            np.zeros((len(stratum_ind), len(len_bin))),
        ),
    }

    # add all variables that have only male and female versions
    for sex in ["M", "F"]:
        data_vars_dict[f"num_{sex}"] = ("stratum_num", np.zeros(len(stratum_ind)))
        data_vars_dict[f"station_1_N_{sex}"] = (
            "stratum_num",
            np.zeros(len(stratum_ind)),
        )
        data_vars_dict[f"unaged_{sex}_wgt_proportion"] = (
            "stratum_num",
            np.zeros(len(stratum_ind)),
        )
        data_vars_dict[f"len_dist_station1_normalized_{sex}"] = (
            ["stratum_num", "len_bin"],
            np.zeros((len(stratum_ind), len(len_bin))),
        )

    # add all variables that have male, female, and all versions
    for sex in ["M", "F", "all"]:
        data_vars_dict[f"len_age_weight_prop_{sex}"] = (
            "stratum_num",
            np.zeros(len(stratum_ind)),
        )
        data_vars_dict[f"len_age_dist_{sex}"] = (
            ["stratum_num", "len_bin", "age_bin"],
            np.zeros((len(stratum_ind), len(len_bin), len(age_bin))),
        )
        data_vars_dict[f"len_age_weight_dist_{sex}"] = (
            ["stratum_num", "len_bin", "age_bin"],
            np.zeros((len(stratum_ind), len(len_bin), len(age_bin))),
        )

        data_vars_dict[f"len_age_weight_dist_{sex}_normalized"] = (
            ["stratum_num", "len_bin", "age_bin"],
            np.zeros((len(stratum_ind), len(len_bin), len(age_bin))),
        )

    # initialize Dataset that will hold length age distributions
    ds = xr.Dataset(
        data_vars=data_vars_dict,
        coords={
            "stratum_num": ("stratum_num", stratum_ind),
            "len_bin": ("len_bin", len_bin),
            "age_bin": ("age_bin", age_bin),
        },
    )

    return ds


def _get_len_wgt_distributions(survey) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Computes the length-weight distributions for each gender
    and all animals using the specimen data

    Parameters
    ----------
    survey : Survey
        An initialized Survey object

    Returns
    -------
    len_wgt_M: np.ndarray
        The length-weight distribution for males
    len_wgt_F: np.ndarray
        The length-weight distribution for females
    len_wgt_all: np.ndarray
        The length-weight distribution for all genders

    Notes
    -----
    These calculations use all specimen data, rather than the
    specimen data with NA values dropped. This is necessary to
    match the Matlab output.
    """

    len_wgt_M = survey.bio_calc._generate_length_val_conversion(
        len_name="length",
        val_name="weight",
        df=survey.bio_calc.specimen_df[survey.bio_calc.specimen_df["sex"] == 1],
    )

    len_wgt_F = survey.bio_calc._generate_length_val_conversion(
        len_name="length",
        val_name="weight",
        df=survey.bio_calc.specimen_df[survey.bio_calc.specimen_df["sex"] == 2],
    )

    len_wgt_all = survey.bio_calc._generate_length_val_conversion(
        len_name="length",
        val_name="weight",
        df=survey.bio_calc.specimen_df,
    )

    return len_wgt_M, len_wgt_F, len_wgt_all


def _set_age_distribution_data(
    survey, stratum: int, ds: xr.Dataset, df_M: pd.DataFrame, df_F: pd.DataFrame
) -> None:
    """
    Computes distributions for each age using the input DataFrames.
    Additionally, assigns these computed quantities to the Dataset ``ds``.

    Parameters
    ----------
    survey : Survey
        An initialized Survey object
    stratum: int
        The stratum corresponding to the input data (used when
        assigning values to ``ds``)
    ds: xr.Dataset
        The Dataset where computed quantities should be assigned to
    df_M: pd.DataFrame
        A DataFrame specifying length, age, and weight measurements
        of a male animal within a stratum
    df_F: pd.DataFrame
        A DataFrame specifying length, age, and weight measurements
        of a female animal within a stratum
    """

    # account for the case when df is a Series
    if isinstance(df_M, pd.Series):
        # get numpy arrays of length, age, and weight
        input_arr_len_M = np.array([df_M.length])
        input_arr_age_M = np.array([df_M.age])
        input_arr_wgt_M = np.array([df_M.weight])

        input_arr_len_F = np.array([df_F.length])
        input_arr_age_F = np.array([df_F.age])
        input_arr_wgt_F = np.array([df_F.weight])

    else:
        # get numpy arrays of length, age, and weight
        input_arr_len_M = df_M.length.values
        input_arr_age_M = df_M.age.values
        input_arr_wgt_M = df_M.weight.values

        input_arr_len_F = df_F.length.values
        input_arr_age_F = df_F.age.values
        input_arr_wgt_F = df_F.weight.values

    # bin the ages
    # TODO: binning is occurring differently than in transect_results.py! It may be
    #  better and more consistent to use the function self.krig.survey.bio_calc._get_bin_ind
    age_bins_ind_M = _get_bin_ind_age(input_arr_age_M, survey.bio_calc.bio_hake_age_bin)
    age_bins_ind_F = _get_bin_ind_age(input_arr_age_F, survey.bio_calc.bio_hake_age_bin)

    # round input lengths
    # TODO: this is necessary to match the Matlab output and may not be necessary!
    input_arr_len_M = np.round(input_arr_len_M)
    input_arr_len_F = np.round(input_arr_len_F)

    # compute distributions for each age bin
    for age_bin in range(len(age_bins_ind_M)):

        # bin those lengths that correspond to the lengths in the given age bin
        # TODO: binning is occurring differently than in transect_results.py!

        for sex in ["M", "F"]:

            if sex == "M":
                input_arr_len = input_arr_len_M
                age_bins_ind = age_bins_ind_M
                input_arr_wgt = input_arr_wgt_M
            else:
                input_arr_len = input_arr_len_F
                age_bins_ind = age_bins_ind_F
                input_arr_wgt = input_arr_wgt_F

            len_bin_ind = _get_bin_ind(
                input_arr_len[age_bins_ind[age_bin]],
                survey.bio_calc.bio_hake_len_bin,
            )

            # get the distribution of weight for a particular age bin
            ds[f"len_age_weight_dist_{sex}"].sel(stratum_num=stratum)[
                :, age_bin
            ] = np.array(
                [np.sum(input_arr_wgt[age_bins_ind[age_bin]][i]) for i in len_bin_ind]
            )

            # get the distribution of lengths for a particular age bin
            ds[f"len_age_dist_{sex}"].sel(stratum_num=stratum)[:, age_bin] = np.array(
                [len(i) for i in len_bin_ind]
            )

        # get the distribution of weight for a particular age bin for all genders
        ds.sel(stratum_num=stratum).len_age_weight_dist_all[:, age_bin] = (
            ds.sel(stratum_num=stratum).len_age_weight_dist_M[:, age_bin]
            + ds.sel(stratum_num=stratum).len_age_weight_dist_F[:, age_bin]
        )

        # get the distribution of lengths for a particular age bin for all genders
        ds.sel(stratum_num=stratum).len_age_dist_all[:, age_bin] = (
            ds.sel(stratum_num=stratum).len_age_dist_M[:, age_bin]
            + ds.sel(stratum_num=stratum).len_age_dist_F[:, age_bin]
        )

    # obtain normalized distributions
    for sex in ["all", "M", "F"]:
        ds[f"len_age_weight_dist_{sex}_normalized"].sel(stratum_num=stratum)[:, :] = ds[
            f"len_age_weight_dist_{sex}"
        ].sel(stratum_num=stratum) / np.nansum(
            ds[f"len_age_weight_dist_{sex}"].sel(stratum_num=stratum)
        )


def _set_total_weight(
    survey, haul_nums: np.ndarray, stratum: int, ds: xr.Dataset
) -> float:
    """
    Computes the total weight within the stratum using data from
    both stations. Additionally, assigns the variable ``total_weight``
    to ``ds``.

    Parameters
    ----------
    survey : Survey
        An initialized Survey object
    haul_nums: np.ndarray
        A one dimensional array of haul numbers within the stratum
        under consideration
    stratum: int
        The stratum corresponding to the input data (used when
        assigning values to ``ds``)
    ds: xr.Dataset
        The Dataset where computed quantities should be assigned to

    Returns
    -------
    wgt_station_1: float
        The total weight at station 1 for the stratum
    """

    # calculate the total weight within the stratum for station 1
    wgt_station_1 = np.array(
        [
            survey.catch_df.loc[j][
                "haul_weight"
            ].sum()  # TODO: make catch_df a bio_calc variable?
            for j in haul_nums
            if (j in survey.catch_df.index) and (j in survey.length_df.index)
        ]
    ).sum()

    # calculate the total weight within the stratum for station 2
    wgt_station_2 = survey.bio_calc.specimen_df.loc[stratum]["weight"].sum()

    # the total weight within the stratum
    ds.total_weight.loc[stratum] = wgt_station_1 + wgt_station_2

    return wgt_station_1


def _get_length_based_wgt_interp(
    haul_nums: np.ndarray,
    df: pd.DataFrame,
    len_bin: np.ndarray,
    len_wgt: np.ndarray,
) -> float:
    """
    Obtains the weight of animals within a stratum using ``length_df``
    data and interpolation.

    Parameters
    ----------
    haul_nums: np.ndarray
        A one dimensional array of haul numbers within the stratum
        under consideration
    df: pd.DataFrame
        A DataFrame containing data from ``length_df``
    len_bin: np.ndarray
        Length bin values
    len_wgt: np.ndarray
        The length-weight distribution corresponding to ``df``

    Returns
    -------
    final_wgt: float
        The total weight of animals within the stratum
    """

    # initialize weight value
    final_wgt = 0.0

    # obtain the hauls that are within df
    len_haul = [round(df.loc[j]["length"]) for j in haul_nums if j in df.index]

    if len_haul:

        # get the number of lengths for each haul
        len_haul_counts = np.concatenate(
            [df.loc[j]["length_count"].values for j in haul_nums if j in df.index]
        )

        # calculate the weight of animals within the stratum
        len_haul = np.concatenate(len_haul)
        final_wgt = (np.interp(len_haul, len_bin, len_wgt) * len_haul_counts).sum()

    return final_wgt


def _get_length_df_based_wgt(
    survey,
    stratum: int,
    ds: xr.Dataset,
    haul_nums: np.ndarray,
    length_df_M: pd.DataFrame,
    length_df_F: pd.DataFrame,
    len_bin: np.ndarray,
    len_wgt_M: np.ndarray,
    len_wgt_F: np.ndarray,
    len_wgt_all: np.ndarray,
) -> Tuple[float, float]:
    """
    Computes the weight of a particular stratum based off of the
    length data. Additionally, assigns values to the ``ds`` variable
    ``weight_len_all_normalized``.

    Parameters
    ----------
    survey : Survey
        An initialized Survey object
    stratum: int
        The stratum corresponding to the input data (used when
        assigning values to ``ds``)
    ds: xr.Dataset
        The Dataset where computed quantities should be assigned to
    haul_nums: np.ndarray
        A one dimensional array of haul numbers within the stratum
        under consideration
    length_df_M: pd.DataFrame
        Male data obtained from ``length_df``
    length_df_F: pd.DataFrame
        Female data obtained from ``length_df``
    len_bin: np.ndarray
        Length bin values
    len_wgt_M: np.ndarray
        The length-weight distribution for males
    len_wgt_F: np.ndarray
        The length-weight distribution for females
    len_wgt_all: np.ndarray
        The length-weight distribution for all genders

    Returns
    -------
    male_wgt: float
        The total weight of males within the stratum
    female_wgt: float
        The total weight of females within the stratum
    """

    # calculate the weight of male and female animals within the stratum
    male_wgt = _get_length_based_wgt_interp(haul_nums, length_df_M, len_bin, len_wgt_M)
    female_wgt = _get_length_based_wgt_interp(
        haul_nums, length_df_F, len_bin, len_wgt_F
    )

    # obtain all haul numbers within the stratum and in length_df
    hauls_in_all = [j for j in haul_nums if j in survey.length_df.index]

    if hauls_in_all:

        for sex in ["M", "F"]:

            if sex == "M":
                df = length_df_M
            else:
                df = length_df_F

            # get normalized length distribution of length_df gender data
            ds[f"len_dist_station1_normalized_{sex}"].sel(stratum_num=stratum)[
                :
            ] = survey.bio_calc._get_distribution_lengths_station_1(
                df.loc[hauls_in_all]
            )

            # store the number of animals in station 1 for the stratum
            ds[f"station_1_N_{sex}"].loc[stratum] = df.loc[hauls_in_all][
                "length_count"
            ].sum()

        # get normalized length distribution of length_df data
        len_dist_station1_normalized = (
            survey.bio_calc._get_distribution_lengths_station_1(
                survey.length_df.loc[hauls_in_all]
            )
        )

        # store the number of animals in station 1 for the stratum
        ds.station_1_N.loc[stratum] = survey.length_df.loc[hauls_in_all][
            "length_count"
        ].sum()

        # obtain the weight per unit length distribution
        weight_len_all = len_wgt_all * len_dist_station1_normalized

        # normalized weight per unit length distribution
        ds.sel(stratum_num=stratum).weight_len_all_normalized[:] = (
            weight_len_all / weight_len_all.sum()
        )

    return female_wgt, male_wgt


def _set_proportion_parameters(
    stratum: int,
    ds: xr.Dataset,
    male_wgt: float,
    female_wgt: float,
    wgt_station_1: float,
) -> None:
    """
    Calculates and assigns proportion parameters to ``ds`` for a stratum.

    Parameters
    ----------
    stratum: int
        The stratum corresponding to the input data (used when
        assigning values to ``ds``)
    ds: xr.Dataset
        The Dataset where computed quantities should be assigned to
    male_wgt: float
        The total weight of males within the stratum
    female_wgt: float
        The total weight of females within the stratum
    wgt_station_1: float
        The total weight at station 1 for the stratum
    """

    # calculate and assign the len_age_weight proportions
    for sex in ["all", "M", "F"]:
        ds[f"len_age_weight_prop_{sex}"].loc[stratum] = np.nansum(
            ds[f"len_age_weight_dist_{sex}"].sel(stratum_num=stratum)
        ) / ds.total_weight.sel(stratum_num=stratum)

    # calculate and assign aged proportions
    ds.aged_proportion.loc[stratum] = (
        ds.len_age_weight_prop_M.loc[stratum] + ds.len_age_weight_prop_F.loc[stratum]
    )

    # calculate and assign the unaged proportions
    ds.unaged_proportion.loc[stratum] = 1.0 - ds.aged_proportion.loc[stratum]

    # obtain the normalized weight of station 1 for males and females
    if (male_wgt != 0.0) and (female_wgt != 0.0):
        nM_wgt1 = wgt_station_1 * male_wgt / (male_wgt + female_wgt)
        nF_wgt1 = wgt_station_1 * female_wgt / (male_wgt + female_wgt)
    else:
        nM_wgt1 = 0.0
        nF_wgt1 = 0.0

    # obtain length and gender based weight proportion
    Len_M_wgt_proportion = nM_wgt1 / ds.total_weight.sel(stratum_num=stratum).values
    Len_F_wgt_proportion = nF_wgt1 / ds.total_weight.sel(stratum_num=stratum).values

    # obtain the proportion of males and females
    if (Len_M_wgt_proportion == 0.0) and (Len_F_wgt_proportion == 0.0):
        M_proportion = 0.5
        F_proportion = 0.5
    else:
        M_proportion = Len_M_wgt_proportion / (
            Len_M_wgt_proportion + Len_F_wgt_proportion
        )
        F_proportion = Len_F_wgt_proportion / (
            Len_M_wgt_proportion + Len_F_wgt_proportion
        )

    # calculate and assign the unaged weight proportion of males and females
    ds.unaged_M_wgt_proportion.loc[stratum] = (
        ds.unaged_proportion.loc[stratum].values * M_proportion
    )
    ds.unaged_F_wgt_proportion.loc[stratum] = (
        ds.unaged_proportion.loc[stratum].values * F_proportion
    )


def generate_bin_ds(survey) -> xr.Dataset:
    """
    Creates a Dataset containing parameters that are
    necessary to compute Kriging result variables at
    each mesh point. The parameters created are associated
    with the binning of data over length, age, and weight.

    Parameters
    ----------
    survey : Survey
        An initialized Survey object

    Returns
    -------
    ds: xr.Dataset
        A Dataset containing useful parameters
    """

    # obtain specimen DataFrames without NA values
    # TODO: This is necessary to match the Matlab output
    #  in theory this should be done when we load the df,
    #  however, this changes the results slightly.
    spec_drop = survey.bio_calc.specimen_df.dropna(how="any")
    spec_drop_M = spec_drop[spec_drop["sex"] == 1]
    spec_drop_F = spec_drop[spec_drop["sex"] == 2]

    # obtain gender based length DataFrames
    # TODO: does not use bio_calc.length_df value (it would give incorrect
    #  answers as that df drops values)
    length_df_M = survey.length_df[survey.length_df["sex"] == 1]
    length_df_F = survey.length_df[survey.length_df["sex"] == 2]

    # get all unique stratum values
    stratum_ind = spec_drop.index.unique()

    # get age and length bins (created to reduce clutter)
    len_bin = survey.params["bio_hake_len_bin"]
    age_bin = survey.params["bio_hake_age_bin"]

    # initialize Dataset that will hold all parameters
    ds = _initialize_ds(stratum_ind, len_bin, age_bin)

    # obtain a mapping of hauls to strata
    haul_vs_stratum = survey.bio_calc.strata_df.reset_index()[
        ["haul_num", "stratum_num"]
    ]

    # get length-weight distributions (includes all ages in quantity)
    len_wgt_M, len_wgt_F, len_wgt_all = _get_len_wgt_distributions(survey)

    for i in stratum_ind:

        # obtain haul numbers that are in the stratum i
        haul_nums = haul_vs_stratum[haul_vs_stratum["stratum_num"] == i][
            "haul_num"
        ].values

        # calculate and set the number of animals in the stratum based on gender
        ds.num_M.loc[i] = len(spec_drop_M.loc[i])
        ds.num_F.loc[i] = len(spec_drop_F.loc[i])

        # set age distribution related data
        _set_age_distribution_data(
            survey=survey,
            stratum=i,
            ds=ds,
            df_M=spec_drop_M.loc[i],
            df_F=spec_drop_F.loc[i],
        )

        # assign the total weight in both station to ds and obtain the weight in station 1
        wgt_station_1 = _set_total_weight(survey, haul_nums, i, ds)

        # get length_df based weight for males and females
        female_wgt, male_wgt = _get_length_df_based_wgt(
            survey,
            i,
            ds,
            haul_nums,
            length_df_M,
            length_df_F,
            len_bin,
            len_wgt_M,
            len_wgt_F,
            len_wgt_all,
        )

        # calculate and assign proportion parameters to ds
        _set_proportion_parameters(i, ds, male_wgt, female_wgt, wgt_station_1)

    return ds
