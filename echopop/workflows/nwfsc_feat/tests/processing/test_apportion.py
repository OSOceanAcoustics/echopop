import numpy as np
import pandas as pd
import pytest

from echopop.workflows.nwfsc_feat import apportionment as apportion


def test_mesh_biomass_to_nasc(
    apportion_mesh,
    apportion_weight_proportions,
    apportion_stratum_weights,
    apportion_stratum_sigma_bs,
):
    """
    Test biomass to NASC conversion
    """

    # Create copy of apportion mesh DataFrame
    mesh_data_df = apportion_mesh.copy()

    # Add biomass
    mesh_data_df["biomass"] = mesh_data_df["biomass_density"] * mesh_data_df["area"]

    # This should not raise an error (operation done in-place)
    apportion.mesh_biomass_to_nasc(
        mesh_data_df=mesh_data_df,
        biodata=apportion_weight_proportions,
        mesh_biodata_link={"mesh_stratum": "bio_stratum"},
        stratum_weights_df=apportion_stratum_weights,
        stratum_sigma_bs_df=apportion_stratum_sigma_bs,
        group_by=["contrast"],
    )

    # Verify in-place operation
    test_result = apportion.mesh_biomass_to_nasc(
        mesh_data_df=mesh_data_df,
        biodata=apportion_weight_proportions,
        mesh_biodata_link={"mesh_stratum": "bio_stratum"},
        stratum_weights_df=apportion_stratum_weights,
        stratum_sigma_bs_df=apportion_stratum_sigma_bs,
        group_by=["contrast"],
    )

    # Check
    assert test_result is None

    # Check shapes
    assert mesh_data_df.shape == (3, 10)

    # Check columns
    assert sorted(mesh_data_df.columns) == sorted(
        [
            "mesh_stratum",
            "biomass_density",
            "area",
            "biomass",
            "biomass_A",
            "biomass_B",
            "abundance",
            "abundance_A",
            "abundance_B",
            "nasc",
        ]
    )

    # Check values
    assert np.allclose(mesh_data_df["biomass"], [1e7, 4e7, 9e7])
    assert np.allclose(mesh_data_df["biomass_A"], [4.0e6, 2.0e7, 6.3e7])
    assert np.allclose(mesh_data_df["biomass_B"], [6.0e6, 2.0e7, 2.7e7])
    assert np.allclose(mesh_data_df["abundance"], [1.0e7, 2.0e7, 3.0e7])
    assert np.allclose(mesh_data_df["abundance_A"], [4.0e6, 1.0e7, 2.1e7])
    assert np.allclose(mesh_data_df["abundance_B"], [6.0e6, 1.0e7, 9.0e6])
    assert np.allclose(np.round(mesh_data_df["nasc"], 4), [125663.7061, 25132.7412, 3769.9112])


def test_distribute_population_estimates_kriged(
    apportion_mesh_with_nasc, apportion_weight_proportions, apportion_number_proportions
):
    """
    Test distribution of kriged estimates
    """

    # Create copy of apportion mesh DataFrame
    mesh_data_df = apportion_mesh_with_nasc.copy()

    # Add biomass
    mesh_data_df["biomass"] = mesh_data_df["biomass_density"] * mesh_data_df["area"]

    # Test-distribute biomass
    biomass_tables = apportion.distribute_population_estimates(
        data=mesh_data_df,
        proportions=apportion_weight_proportions,
        variable="biomass",
        group_by=["contrast", "index_bin", "extra_bin"],
        stratify_by=["bio_stratum"],
        data_proportions_link={"mesh_stratum": "bio_stratum"},
    )

    # Check shapes
    assert len(biomass_tables) == 2
    assert [len(df) for df in biomass_tables.values()] == [12, 6]

    # Check columns
    assert ["subgroup1", "subgroup2"] == list(biomass_tables.keys())

    # Check values
    assert all(biomass_tables["subgroup1"].sum() == np.array([5.0e6, 3.2e7, 4.5e7]))
    assert all(biomass_tables["subgroup2"].sum() == np.array([5.0e6, 8.0e6, 4.5e7]))

    # Test the same for abundances using a different format for DataFrame
    abundance_tables = apportion.distribute_population_estimates(
        data=mesh_data_df,
        proportions=apportion_number_proportions,
        variable="abundance",
        group_by=["contrast", "index_bin", "extra_bin"],
        stratify_by=["bio_stratum"],
        data_proportions_link={"mesh_stratum": "bio_stratum"},
    )

    # Check shapes
    assert len(abundance_tables) == 2
    assert [len(df) for df in abundance_tables.values()] == [12, 6]

    # Check columns
    assert ["subgroup1", "subgroup2"] == list(abundance_tables.keys())

    # Check values
    assert all(abundance_tables["subgroup1"].sum() == np.array([5e6, 16e6, 15e6]))
    assert all(abundance_tables["subgroup2"].sum() == np.array([5e6, 4e6, 15e6]))


def test_distribute_unaged_from_aged(apportion_biomass_table):
    """
    Test group-standardization of kriged estimates
    """

    # Standardize 'subgroup2' based on 'subgroup1' -- No imputation
    std_tbl_no_imp = apportion.distribute_unaged_from_aged(
        population_table=apportion_biomass_table["subgroup2"],
        reference_table=apportion_biomass_table["subgroup1"],
        group_by=["contrast"],
        impute=False,
    )

    # Check shape
    assert std_tbl_no_imp.shape == (3, 4)

    # Check indices
    assert list(std_tbl_no_imp.index.names) == ["index_bin"]
    assert all(std_tbl_no_imp.index == pd.Index([5, 10, 15]))
    assert list(std_tbl_no_imp.columns.names) == ["extra_bin", "contrast"]
    assert all(
        std_tbl_no_imp.columns
        == pd.MultiIndex.from_tuples([(1, "A"), (2, "A"), (1, "B"), (2, "B")])
    )

    # Check values
    # --- Expected values
    expected_no_imp = pd.DataFrame(
        {
            "extra_bin": [1, 2, 1, 2],
            "contrast": ["A", "A", "B", "B"],
            "value": [3.94e7, 1.56e7, 0.00, 3.00e6],
        }
    ).set_index(["extra_bin", "contrast"])["value"]

    assert all(std_tbl_no_imp.sum() == expected_no_imp)

    # Standardize 'subgroup2' based on 'subgroup1' -- Imputation
    # ---- Mock a missing value
    apportion_biomass_table_imp = apportion_biomass_table.copy()
    apportion_biomass_table_imp["subgroup1"].loc[[("A", 5, 1), ("B", 15, 2)]] = 0.0
    # ---- Run standardization
    std_tbl_imp = apportion.distribute_unaged_from_aged(
        population_table=apportion_biomass_table_imp["subgroup2"],
        reference_table=apportion_biomass_table_imp["subgroup1"],
        group_by=["contrast"],
        impute=True,
        impute_variable=["extra_bin"],
    )

    # Check inequality
    assert not std_tbl_imp.equals(std_tbl_no_imp)

    # Check shape
    std_tbl_imp.shape == (3, 4)

    # Check indices
    assert list(std_tbl_imp.index.names) == ["index_bin"]
    assert all(std_tbl_imp.index == pd.Index([5, 10, 15]))
    assert list(std_tbl_imp.columns.names) == ["extra_bin", "contrast"]
    assert all(
        std_tbl_imp.columns == pd.MultiIndex.from_tuples([(1, "A"), (2, "A"), (1, "B"), (2, "B")])
    )

    # Check values
    # --- Expected values
    expected_imp = pd.DataFrame(
        {
            "extra_bin": [1, 2, 1, 2],
            "contrast": ["A", "A", "B", "B"],
            "value": [2.44e7, 3.06e7, 0.00, 3.00e6],
        }
    ).set_index(["extra_bin", "contrast"])["value"]

    assert all(std_tbl_imp.sum() == expected_imp)


def test_sum_population_tables(apportion_biomass_table_with_standardized):
    """
    Test functionality for combining the various apportionment tables
    """

    # Define arguments that will be used for manual verification in tests
    TABLE_INDEX = ["index_bin"]
    TABLE_COLUMNS = ["extra_bin", "contrast"]

    # Combine tables
    df_biomass_table = apportion.sum_population_tables(
        population_table=apportion_biomass_table_with_standardized,
        table_names=["subgroup1", "standardized_subgroup2"],
        table_index=TABLE_INDEX,
        table_columns=TABLE_COLUMNS,
    )

    # Check shape
    assert df_biomass_table.shape == (3, 4)

    # Check indices
    assert list(df_biomass_table.index.names) == ["index_bin"]
    assert all(df_biomass_table.index == pd.Index([5, 10, 15]))
    assert list(df_biomass_table.columns.names) == ["extra_bin", "contrast"]
    assert all(
        df_biomass_table.columns
        == pd.MultiIndex.from_tuples([(1, "A"), (1, "B"), (2, "A"), (2, "B")])
    )

    # Check values
    # ---- Break out each defined group for value verification
    DF_A = apportion_biomass_table_with_standardized["subgroup1"].copy()
    DF_B = apportion_biomass_table_with_standardized["standardized_subgroup2"].copy()
    DF_SUM = DF_A.sum(axis=1).unstack(TABLE_COLUMNS) + DF_B

    assert np.allclose(df_biomass_table - DF_SUM, 0.0)


def test_reallocate_excluded_estimates(apportion_combined_biomass_table):
    """
    Test functionality for redistributing apportioned values within a table
    """

    # Create copy for multiple tests
    TEST_DF = apportion_combined_biomass_table.copy()

    # Empty exclusion filter
    df_no_filter = apportion.reallocate_excluded_estimates(
        population_table=TEST_DF,
        exclusion_filter={},
        group_by=["contrast"],
    )

    # Check equality with `TEST_DF`
    assert df_no_filter.equals(TEST_DF)
    assert all(df_no_filter.index == TEST_DF.index)
    assert sorted(df_no_filter.columns) == sorted(TEST_DF.columns)

    # No grouping -- distribute across all contrasts not excluded
    df_no_grouping = apportion.reallocate_excluded_estimates(
        population_table=TEST_DF,
        exclusion_filter={"contrast": "A"},
        group_by=[],
    )

    # Check shape and indices with `TEST_DF`
    assert df_no_grouping.shape == TEST_DF.shape
    assert all(df_no_grouping.index == TEST_DF.index)
    assert sorted(df_no_grouping.columns) == sorted(TEST_DF.columns)

    # Check total sum
    assert np.isclose(df_no_grouping.sum().sum() - TEST_DF.sum().sum(), 0.0)

    # Check sliced sums
    assert all(df_no_grouping.sum(axis=1) != TEST_DF.sum(axis=1))
    assert all(df_no_grouping[[(1, "A"), (2, "A")]] == 0.0)

    # Grouping -- defined by column in exclusion filter
    df_grouping_column = apportion.reallocate_excluded_estimates(
        population_table=TEST_DF,
        exclusion_filter={"contrast": "A"},
        group_by=["contrast"],
    )

    # Check shape and indices with `TEST_DF`
    assert df_grouping_column.shape == TEST_DF.shape
    assert all(df_grouping_column.index == TEST_DF.index)
    assert sorted(df_grouping_column.columns) == sorted(TEST_DF.columns)

    # Check total sum
    assert np.isclose(df_grouping_column.sum().sum() - TEST_DF.sum().sum(), 0.0)

    # Check sliced sums
    assert all(df_grouping_column.sum(axis=1) != TEST_DF.sum(axis=1))
    assert all(df_grouping_column[[(1, "A"), (2, "A")]] == 0.0)

    # Grouping -- defined by extra bin
    df_grouping_col2 = apportion.reallocate_excluded_estimates(
        population_table=TEST_DF,
        exclusion_filter={"extra_bin": 1},
        group_by=["contrast"],
    )

    # Check shape and indices with `TEST_DF`
    assert df_grouping_col2.shape == TEST_DF.shape
    assert all(df_grouping_col2.index == TEST_DF.index)
    assert sorted(df_grouping_col2.columns) == sorted(TEST_DF.columns)

    # Check total sum
    assert np.isclose(df_grouping_col2.sum().sum() - TEST_DF.sum().sum(), 0.0)

    # Check sliced sums
    assert all(df_grouping_col2.sum(axis=1) != TEST_DF.sum(axis=1))
    assert all(df_grouping_col2[1] == 0.0)

    # Grouping -- defined by index
    df_grouping_index = apportion.reallocate_excluded_estimates(
        population_table=TEST_DF,
        exclusion_filter={"index_bin": 5},
        group_by=["contrast"],
    )

    # Check shape and indices with `TEST_DF`
    assert df_grouping_index.shape == TEST_DF.shape
    assert all(df_grouping_index.index == TEST_DF.index)
    assert sorted(df_grouping_index.columns) == sorted(TEST_DF.columns)

    # Check total sum
    assert np.isclose(df_grouping_index.sum().sum() - TEST_DF.sum().sum(), 0.0)

    # Check sliced sums
    assert all(df_grouping_index.sum(axis=1) != TEST_DF.sum(axis=1))
    assert np.isclose(df_grouping_index.loc[5].sum(), 0.0)

    # Full filtering -- raises warning

    with pytest.warns(UserWarning) as record:
        apportion.reallocate_excluded_estimates(
            population_table=TEST_DF,
            exclusion_filter={"contrasts": ["A", "B"]},
            group_by=["contrast"],
        )
    assert len(record) == 1
