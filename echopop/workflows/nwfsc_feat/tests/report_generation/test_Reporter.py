import numpy as np
import pandas as pd

from echopop.workflows.nwfsc_feat import Reporter


def test_Reporter(tmp_path):
    """
    Test Reporter class initialization
    """

    # Initialize
    reports = Reporter(tmp_path, verbose=False)

    # Run tests
    assert reports.save_directory == tmp_path
    assert reports.verbose is False

    # Define sheetnames
    SHEETNAME = "Test"
    SHEETNAMES = {"male": "Male", "female": "Female", "all": "All"}

    # Run age-length haul counts report
    specimen_df = pd.DataFrame(
        {
            "length": np.random.random(12),
            "sex": np.tile(["male", "female", "all"], 4),
            "haul_num": np.repeat([1, 2, 3, 4], 3),
            "length_bin": np.tile(pd.IntervalIndex.from_tuples([(2, 4), (4, 6), (6, 8)]), 4),
        }
    )
    specimen_df["length_bin"] = pd.CategoricalIndex(specimen_df["length_bin"])

    reports.aged_length_haul_counts_report("haul.xlsx", SHEETNAMES, specimen_df)
    assert (tmp_path / "haul.xlsx").exists()

    # Run total age-length haul counts report
    length_df = pd.DataFrame(
        {
            "length": np.random.random(12),
            "length_count": np.random.randint(1, 20, size=12),
            "sex": np.tile(["male", "female", "all"], 4),
            "haul_num": np.repeat([1, 2, 3, 4], 3),
            "length_bin": np.tile(pd.IntervalIndex.from_tuples([(2, 4), (4, 6), (6, 8)]), 4),
        }
    )
    length_df["length_bin"] = pd.CategoricalIndex(length_df["length_bin"])

    reports.total_length_haul_counts_report(
        "total_haul.xlsx", SHEETNAMES, {"specimen": specimen_df, "length": length_df}
    )
    assert (tmp_path / "total_haul.xlsx").exists()

    # Basic transect population results report
    transect_df = pd.DataFrame(
        {
            "transect_num": [1, 2, 3],
            "region_id": [999] * 3,
            "distance_s": [1, 2, 3],
            "distance_e": [2, 3, 4],
            "longitude": np.random.random(3),
            "latitude": np.random.random(3),
            "stratum_arbitrary": [10, 11, 12],
            "bottom_depth": np.random.random(3),
            "nasc": np.random.random(3),
            "abundance_male": np.random.random(3),
            "abundance_female": np.random.random(3),
            "abundance": np.random.random(3),
            "number_density": np.random.random(3),
            "number_density_male": np.random.random(3),
            "number_density_female": np.random.random(3),
            "biomass": np.random.random(3),
            "biomass_male": np.random.random(3),
            "biomass_female": np.random.random(3),
            "biomass_density_male": np.random.random(3),
            "biomass_density_female": np.random.random(3),
            "biomass_density": np.random.random(3),
            "layer_mean_depth": np.random.random(3),
            "layer_height": np.random.random(3),
            "transect_spacing": np.random.random(3),
            "distance_interval": np.random.random(3),
            "nasc_proportion": np.random.random(3),
        }
    )
    sigma_bs_stratum_df = pd.DataFrame(
        {
            "sigma_bs": np.random.random(3),
        },
        index=pd.Index([10, 11, 12], name="stratum_arbitrary"),
    )
    weight_stratum_df = pd.DataFrame(
        {
            "all": np.random.random(3),
            "female": np.random.random(3),
            "male": np.random.random(3),
        },
        index=pd.Index([10, 11, 12], name="stratum_arbitrary"),
    )
    weight_stratum_df.columns.name = "sex"

    reports.transect_population_results_report(
        filename="transect_results.xlsx",
        sheetname=SHEETNAME,
        transect_data=transect_df,
        weight_strata_data=weight_stratum_df,
        sigma_bs_stratum=sigma_bs_stratum_df,
        stratum_name="stratum_arbitrary",
    )
    assert (tmp_path / "transect_results.xlsx").exists()

    # Aged kriged biomass report
    kriged_df = pd.DataFrame(
        {
            "longitude": np.random.random(3),
            "latitude": np.random.random(3),
            "biomass": np.random.random(3),
            "biomass_male": np.random.random(3),
            "biomass_female": np.random.random(3),
            "geostratum_arbitrary": [10, 11, 12],
        }
    )

    weight_df = pd.DataFrame(
        [np.random.random(18)] * 3,
        columns=pd.MultiIndex.from_product(
            [
                pd.Index(["female", "male"], name="sex"),
                pd.Index(
                    pd.Categorical(pd.IntervalIndex.from_tuples([(0, 1), (1, 2), (2, 3)])),
                    name="age_bin",
                ),
                pd.Index([10, 11, 12], name="stratum_arbitrary"),
            ]
        ),
        index=pd.Index(
            pd.Categorical(pd.IntervalIndex.from_tuples([(2, 4), (4, 6), (6, 8)])),
            name="length_bin",
        ),
    )

    # ---- No exclude filter
    reports.kriged_aged_biomass_mesh_report(
        filename="aged_kriged_biomass_nofilter.xlsx",
        sheetnames=SHEETNAMES,
        kriged_data=kriged_df,
        weight_data=weight_df,
        kriged_stratum_link={"geostratum_arbitrary": "stratum_arbitrary"},
        exclude_filter={},
    )
    assert (tmp_path / "aged_kriged_biomass_nofilter.xlsx").exists()
    # ---- With exclude filter
    reports.kriged_aged_biomass_mesh_report(
        filename="aged_kriged_biomass_yesfilter.xlsx",
        sheetnames=SHEETNAMES,
        kriged_data=kriged_df,
        weight_data=weight_df,
        kriged_stratum_link={"geostratum_arbitrary": "stratum_arbitrary"},
        exclude_filter={"age_bin": 0.5},
    )
    assert (tmp_path / "aged_kriged_biomass_yesfilter.xlsx").exists()

    # Run kriged length-age abundance table report
    aged_df = pd.DataFrame(
        [np.random.random(3)] * 18,
        index=pd.MultiIndex.from_product(
            [
                pd.Index(
                    pd.Categorical(pd.IntervalIndex.from_tuples([(2, 4), (4, 6), (6, 8)])),
                    name="length_bin",
                ),
                pd.Index(
                    pd.Categorical(pd.IntervalIndex.from_tuples([(0, 1), (1, 2), (2, 3)])),
                    name="age_bin",
                ),
                pd.Index(["male", "female"], name="sex"),
            ]
        ),
        columns=pd.Index([10, 11, 12], name="stratum_arbitrary"),
    )
    unaged_df = pd.DataFrame(
        [np.random.random(3)] * 6,
        index=pd.MultiIndex.from_product(
            [
                pd.Index(
                    pd.Categorical(pd.IntervalIndex.from_tuples([(2, 4), (4, 6), (6, 8)])),
                    name="length_bin",
                ),
                pd.Index(["male", "female"], name="sex"),
            ]
        ),
        columns=pd.Index([10, 11, 12], name="stratum_arbitrary"),
    )

    # ---- No exclude filter
    reports.kriged_length_age_abundance_report(
        filename="kriged_length_age_abundance_nofilter.xlsx",
        sheetnames=SHEETNAMES,
        datatables={"aged": aged_df, "unaged": unaged_df},
        exclude_filter={},
    )
    assert (tmp_path / "kriged_length_age_abundance_nofilter.xlsx").exists()
    # ---- With exclude filter
    reports.kriged_length_age_abundance_report(
        filename="kriged_length_age_abundance_yesfilter.xlsx",
        sheetnames=SHEETNAMES,
        datatables={"aged": aged_df, "unaged": unaged_df},
        exclude_filter={"age_bin": 0.5},
    )
    assert (tmp_path / "kriged_length_age_abundance_yesfilter.xlsx").exists()

    # Run kriged length-aged biomass table report
    # ---- No exclude filter
    reports.kriged_length_age_biomass_report(
        filename="kriged_length_age_biomass_nofilter.xlsx",
        sheetnames=SHEETNAMES,
        datatables={"aged": aged_df, "unaged": unaged_df},
        exclude_filter={},
    )
    assert (tmp_path / "kriged_length_age_biomass_nofilter.xlsx").exists()
    # ---- With exclude filter
    reports.kriged_length_age_biomass_report(
        filename="kriged_length_age_biomass_yesfilter.xlsx",
        sheetnames=SHEETNAMES,
        datatables={"aged": aged_df, "unaged": unaged_df},
        exclude_filter={"age_bin": 0.5},
    )
    assert (tmp_path / "kriged_length_age_biomass_yesfilter.xlsx").exists()

    # Run kriging input report
    reports.kriging_input_report(
        "kriging_input.xlsx",
        SHEETNAME,
        transect_df,
    )
    assert (tmp_path / "kriging_input.xlsx").exists()

    # Run transect aged biomass report
    # ---- No exclude filter
    reports.transect_aged_biomass_report(
        filename="aged_transect_biomass_nofilter.xlsx",
        sheetnames=SHEETNAMES,
        transect_data=transect_df,
        weight_data=weight_df,
        exclude_filter={},
    )
    assert (tmp_path / "aged_transect_biomass_nofilter.xlsx").exists()
    # ---- With exclude filter
    reports.transect_aged_biomass_report(
        filename="aged_transect_biomass_yesfilter.xlsx",
        sheetnames=SHEETNAMES,
        transect_data=transect_df,
        weight_data=weight_df,
        exclude_filter={"age_bin": 0.5},
    )
    assert (tmp_path / "aged_transect_biomass_yesfilter.xlsx").exists()

    # Run transect length-age abundance table report
    reports.transect_length_age_abundance_report(
        filename="transect_length_age_abundance.xlsx",
        sheetnames=SHEETNAMES,
        datatables={"aged": aged_df, "unaged": unaged_df},
    )
    assert (tmp_path / "transect_length_age_abundance.xlsx").exists()

    # Run biomass length-age abundance table report
    reports.transect_length_age_biomass_report(
        filename="transect_length_age_biomass.xlsx", sheetnames=SHEETNAMES, datatable=aged_df
    )
    assert (tmp_path / "transect_length_age_biomass.xlsx").exists()
