import numpy as np
import pandas as pd

from echopop.computation.acoustics import to_dB, to_linear, ts_length_regression


def test_strata_mean_sigma_bs(mock_survey):

    # Re-parameterize `specimen_df` with dummy data
    mock_survey.biology["specimen_df"] = pd.DataFrame(
        {
            "stratum_num": np.repeat([0, 1, 2, 4, 5], 4),
            "haul_num": np.repeat([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 2),
            "species_id": np.append(np.repeat([19350], 19), 43130),
            "length": np.linspace(10, 100, 20),
            "weight": np.linspace(1, 5, 20),
        }
    )

    # Re-parameterize `length_df` with dummy data
    mock_survey.biology["length_df"] = pd.DataFrame(
        {
            "stratum_num": np.repeat([0, 1, 2, 4, 5], 4),
            "haul_num": np.repeat([1, 2, 3, 4, 5], 4),
            "species_id": np.append(np.repeat([19350], 19), 43130),
            "length": np.linspace(10, 100, 20),
            "length_count": np.linspace(10, 100, 20),
        }
    )

    # Re-parameterize `strata_df` with dummy data
    mock_survey.spatial["strata_df"] = pd.DataFrame({"stratum_num": [0, 1, 2, 3, 4, 5, 6]})

    # Dummy parameters
    mock_survey.config["TS_length_regression_parameters"]["pacific_hake"] = {
        "species_code": 19350,
        "TS_L_slope": 20.0,
        "TS_L_intercept": -84.1,
        "length_units": "cm",
    }

    # Define dummy `species_id` code
    species_id = 19350

    # Evaluate object for later comparison
    mock_survey.strata_mean_sigma_bs(species_id)

    # ----------------------------------
    # Run tests: `strata_mean_sigma_bs`
    # ----------------------------------
    # Evaluate whether non-specified `species_id` code are removed
    # ---- `specimen_df`: 20 rows -> 19 rows
    specimen_df_copy = mock_survey.biology["specimen_df"].copy()
    specimen_df_copy = specimen_df_copy[specimen_df_copy.species_id == species_id]
    assert specimen_df_copy.shape[0] == 19

    # ---- `length_df`: 20 rows -> 19 rows
    length_df_copy = mock_survey.biology["length_df"].copy()
    length_df_copy = length_df_copy[length_df_copy.species_id == species_id]
    assert length_df_copy.shape[0] == 19

    # ----------------------------------
    # Next step: Concatenate the two dataframes
    # ---- Re-bin `specimen_df_copy`
    spec_df_reframed = (
        specimen_df_copy.groupby(["haul_num", "stratum_num", "species_id", "length"])
        .apply(lambda x: len(x["length"]))
        .reset_index(name="length_count")
    )

    # ---- Concatenate
    all_length_df = pd.concat([spec_df_reframed, length_df_copy], join="inner")
    assert all_length_df.shape[0] == 38

    # ----------------------------------
    # TS-length parameterization & modeling
    # ---- Fit parameters
    ts_length_parameters = mock_survey.config["TS_length_regression_parameters"]["pacific_hake"]
    slope = ts_length_parameters["TS_L_slope"]
    intercept = ts_length_parameters["TS_L_intercept"]

    # ---- Predict `TS`
    all_length_df["TS"] = ts_length_regression(all_length_df["length"], slope, intercept)
    assert np.isclose(np.median(all_length_df["TS"]), -49.675)

    # ---- Linearize to `sigma_bs`
    all_length_df["sigma_bs"] = to_linear(all_length_df["TS"])
    assert np.isclose(all_length_df["sigma_bs"].mean(), 1.3395e-5)
    assert mock_survey.acoustics["sigma_bs"]["length_binned"].equals(all_length_df)

    # ----------------------------------
    # Calculate mean `sigma_bs` per `haul_num` and then `stratum_num`
    # ---- `haul_num`
    mean_haul_sigma_bs = (
        all_length_df.groupby(["haul_num", "stratum_num", "species_id"])[
            ["sigma_bs", "length_count"]
        ]
        .apply(lambda x: np.average(x["sigma_bs"], weights=x["length_count"]))
        .to_frame("sigma_bs_mean")
        .reset_index()
    )
    assert mock_survey.acoustics["sigma_bs"]["haul_mean"].equals(mean_haul_sigma_bs)
    assert mean_haul_sigma_bs.shape[0] == 14
    assert np.isclose(mean_haul_sigma_bs.sigma_bs_mean.mean(), 1.549e-5)

    # ---- `stratum_num`
    mean_strata_sigma_bs = (
        mean_haul_sigma_bs.groupby(["stratum_num", "species_id"])["sigma_bs_mean"]
        .mean()
        .reset_index()
    )
    assert mean_strata_sigma_bs.shape[0] == 5
    assert np.allclose(
        mean_strata_sigma_bs.sigma_bs_mean,
        np.array([1.659e-6, 5.238e-6, 1.195e-5, 2.145e-5, 3.254e-5]),
    )
    assert any(~mean_strata_sigma_bs.stratum_num.isin([3, 6]))

    # ----------------------------------
    # Add to object as dictionary
    strata_mean_dictionary = {
        "length_binned": all_length_df,
        "haul_mean": mean_haul_sigma_bs,
        "strata_mean": mean_strata_sigma_bs,
    }

    # ----------------------------------
    # `impute_missing_sigma_bs`
    # Run function -- missing strata (3, 6)
    # ---- Collect strata numbers
    strata_options = np.unique(mock_survey.spatial["strata_df"].copy().stratum_num)

    # ---- Pull `strata_mean` dataframe from dictionary
    strata_mean = strata_mean_dictionary["strata_mean"]

    # Evaluate imputed values
    # ---- Check mismatch between present/absent values
    present_strata = np.unique(strata_mean["stratum_num"]).astype(int)
    missing_strata = strata_options[~(np.isin(strata_options, present_strata))]
    assert all(np.array([3, 6]) == missing_strata)

    # Iterate through for imputation
    if len(missing_strata) > 0:

        # Fill missing values with `np.nan`
        sigma_bs_impute = pd.concat(
            [
                strata_mean,
                pd.DataFrame(
                    {
                        "stratum_num": missing_strata,
                        "species_id": np.repeat(
                            np.unique(strata_mean.species_id), len(missing_strata)
                        ),
                        "sigma_bs_mean": np.repeat(np.nan, len(missing_strata)),
                    }
                ),
            ]
        ).sort_values("stratum_num")

        # Loop over `np.nan` values
        for i in missing_strata:
            strata_floor = present_strata[present_strata < i]
            strata_ceil = present_strata[present_strata > i]

            new_stratum_below = np.max(strata_floor) if strata_floor.size > 0 else None
            new_stratum_above = np.min(strata_ceil) if strata_ceil.size > 0 else None

            sigma_bs_indexed = sigma_bs_impute[
                sigma_bs_impute["stratum_num"].isin([new_stratum_below, new_stratum_above])
            ]

            sigma_bs_impute.loc[sigma_bs_impute.stratum_num == i, "sigma_bs_mean"] = (
                sigma_bs_indexed["sigma_bs_mean"].mean()
            )

    # Test against `mock_survey`
    assert mock_survey.acoustics["sigma_bs"]["strata_mean"].equals(sigma_bs_impute)


def test_ts_linear_regression():

    # Dummy variables/inputs
    length_arr = np.array([21.1, 41.2, 81.5, 19.4, 2.3, 6.8, 16.7, 101.4])
    slope = 20.0
    intercept = -54.3

    # Generate TS array from function
    ts_arr = ts_length_regression(length_arr, slope, intercept)

    # Produce TS array manually
    ts_test_arr = slope * np.log10(length_arr) + intercept

    # Check equality
    assert all(ts_arr == ts_test_arr)


def test_acoustic_unit_conversion():

    # Dummy variables/inputs
    ts_values = np.array([-101.4, -33.2, -49.8, -50.2, -81.7, -20.1, -39.6])

    # Generate value arrays from functions -- sigma_bs
    sigma_bs_arr = to_linear(ts_values)
    sigma_bs_test_arr = 10 ** (ts_values / 10.0)
    assert all(sigma_bs_arr == sigma_bs_test_arr)

    # Generate value arrays from functions -- TS
    ts_arr = to_dB(sigma_bs_arr)
    assert all(ts_arr == ts_values)
