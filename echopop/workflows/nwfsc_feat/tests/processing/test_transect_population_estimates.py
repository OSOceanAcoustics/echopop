import numpy as np
import pandas as pd
import pytest

from echopop.workflows.nwfsc_feat import apportionment as apportion, biology


def test_remove_group_from_estimates_nasc_only(sample_transect_dataset, age1_nasc_proportions):
    """Test partitioning with only NASC proportions."""

    # NASC-only
    result = apportion.remove_group_from_estimates(
        transect_data=sample_transect_dataset, group_proportions={"nasc": age1_nasc_proportions}
    )

    # Check shape/columns
    assert set(result.columns) == set(["stratum_ks", "transect_num", "nasc", "area_interval"])

    # Check that NASC was partitioned correctly for stratum 1 (proportion 0.1)
    stratum_1_rows = result[result["stratum_ks"] == 1]
    expected_nasc_values = [90.0, 135.0]
    assert np.allclose(stratum_1_rows["nasc"].values, expected_nasc_values)

    # Check that NASC was partitioned correctly for stratum 2 (proportion 0.15)
    stratum_2_rows = result[result["stratum_ks"] == 2]
    expected_nasc_values = [170.0, 102.0]
    assert np.allclose(stratum_2_rows["nasc"].values, expected_nasc_values)

    # Check that NASC was partitioned correctly for stratum 2 (proportion 0.15)
    stratum_3_rows = result[result["stratum_ks"] == 3]
    expected_nasc_values = [144.0]
    assert np.allclose(stratum_3_rows["nasc"].values, expected_nasc_values)


def test_remove_group_from_estimates_abundance_only(
    sample_transect_dataset, age1_abundance_proportions
):
    """Test partitioning with only abundance proportions."""

    # Abundance-only
    result = apportion.remove_group_from_estimates(
        transect_data=sample_transect_dataset,
        group_proportions={"abundance": age1_abundance_proportions},
    )

    # Check shape/columns
    assert set(result.columns) == set(
        [
            "stratum_ks",
            "transect_num",
            "number_density",
            "abundance",
            "abundance_male",
            "abundance_female",
            "area_interval",
        ]
    )

    # Check abundance
    assert np.allclose(
        result.filter(["abundance", "abundance_female", "abundance_male"]),
        pd.DataFrame(
            {
                "abundance": [440.0, 660.0, 820.0, 492.0, 702.0],
                "abundance_female": [220.0, 330.0, 410.0, 246.0, 351.0],
                "abundance_male": [220.0, 330.0, 410.0, 246.0, 351.0],
            }
        ),
    )


def test_remove_group_from_estimates_biomass_only(
    sample_transect_dataset, age1_biomass_proportions
):
    """Test partitioning with only biomass proportions."""

    # Biomass-only
    result = apportion.remove_group_from_estimates(
        transect_data=sample_transect_dataset,
        group_proportions={"biomass": age1_biomass_proportions},
    )

    # Check shape/columns
    assert set(result.columns) == set(
        [
            "stratum_ks",
            "transect_num",
            "biomass_density",
            "biomass",
            "biomass_male",
            "biomass_female",
            "area_interval",
        ]
    )

    # Check biomass
    assert np.allclose(
        result.filter(["biomass", "biomass_female", "biomass_male"]),
        pd.DataFrame(
            {
                "biomass": [230.0, 345.0, 430.0, 258.0, 369.0],
                "biomass_female": [115.0, 172.5, 215.0, 129.0, 184.5],
                "biomass_male": [115.0, 172.5, 215.0, 129.0, 184.5],
            }
        ),
    )


def test_remove_group_from_estimates_all_variables(
    sample_transect_dataset,
    age1_nasc_proportions,
    age1_abundance_proportions,
    age1_biomass_proportions,
):
    """Test partitioning with all three variable types."""

    # All variables
    result = apportion.remove_group_from_estimates(
        transect_data=sample_transect_dataset,
        group_proportions={
            "nasc": age1_nasc_proportions,
            "abundance": age1_abundance_proportions,
            "biomass": age1_biomass_proportions,
        },
    )

    # Verify all variables were partitioned
    assert all(
        [
            col in result.columns
            for col in [
                "nasc",
                "abundance",
                "abundance_female",
                "abundance_male",
                "number_density",
                "biomass",
                "biomass_male",
                "biomass_female",
                "biomass_density",
            ]
        ]
    )
    assert len(result) == len(sample_transect_dataset)


def test_remove_group_from_estimates_empty_proportions(
    sample_transect_dataset,
):
    """Test partitioning with all three variable types."""

    # Input dictionary CANNOT be empty
    with pytest.raises(ValueError):
        assert apportion.remove_group_from_estimates(
            transect_data=sample_transect_dataset,
            group_proportions={},
        )


def test_compute_abundance_basic(biology_test_dataset):
    """Test basic abundance calculation."""
    test_data = biology_test_dataset.copy()

    biology.compute_abundance(test_data)

    # Verify abundance = area_interval * number_density
    expected_abundance = test_data["area_interval"] * test_data["number_density"]
    assert np.allclose(test_data["abundance"], expected_abundance)


def test_compute_abundance_with_proportions(
    biology_test_dataset, number_proportions_partition_dict
):
    """Test abundance calculation with number proportions."""
    test_data = biology_test_dataset.copy()

    biology.compute_abundance(
        test_data,
        stratify_by=["stratum_ks"],
        group_by=["sex"],
        number_proportions=number_proportions_partition_dict,
    )

    # Check that grouped abundance columns were created
    assert "abundance_female" in test_data.columns
    assert "abundance_male" in test_data.columns
    assert "abundance" in test_data.columns


def test_compute_biomass_with_weights(biology_test_dataset, average_weight_df):
    """Test biomass calculation with average weights."""
    test_data = biology_test_dataset.copy()

    biology.compute_biomass(
        dataset=test_data,
        stratify_by=["stratum_ks"],
        group_by=["sex"],
        df_average_weight=average_weight_df,
    )

    # Check that biomass columns exist
    assert all(
        [
            col in test_data.columns
            for col in [
                "biomass",
                "biomass_density",
                "biomass_density_female",
                "biomass_density_male",
            ]
        ]
    )


def test_matrix_multiply_grouped_table(biology_test_dataset, average_weight_df):
    """Test matrix multiplication with grouped table."""
    test_data = biology_test_dataset.copy()
    test_data["abundance_female"] = [100, 120, 150, 180, 110, 140]
    test_data["abundance_male"] = [100, 120, 150, 180, 110, 140]

    # Set up the table with proper index
    weight_table = average_weight_df.copy()
    weight_table["all"] = weight_table.mean(axis=1)

    biology.matrix_multiply_grouped_table(
        dataset=test_data,
        table=weight_table,
        variable="abundance",
        output_variable="biomass",
        group="sex",
    )

    # Check that biomass columns were created
    assert "biomass" in test_data.columns
    assert "biomass_female" in test_data.columns
    assert "biomass_male" in test_data.columns
