import numpy as np
import pandas as pd
import pytest

from echopop.workflows.nwfsc_feat import apportionment as apportion


@pytest.fixture
def apportion_mesh():
    """Mesh dataset for apportionment tests"""

    return pd.DataFrame(
        {"mesh_stratum": [1, 2, 3], "biomass_density": [1e6, 2e6, 3e6], "area": [10, 20, 30]}
    )


@pytest.fixture
def apportion_number_proportions():
    """Number proportions for apportionment tests"""

    subgroup1 = pd.DataFrame(
        {
            "bio_stratum": np.tile([1, 2, 3], 6),
            "contrast": np.repeat(["A", "B"], 9),
            "index_bin": np.concatenate([np.repeat([5, 10, 15], 3), np.repeat([5, 10, 15], 3)]),
            "extra_bin": np.tile([1, 2], 9),
            "proportion_overall": [
                0.10,
                0.20,
                0.10,
                0.10,
                0.10,
                0.00,
                0.00,
                0.00,
                0.10,
                0.20,
                0.00,
                0.00,
                0.00,
                0.20,
                0.10,
                0.10,
                0.30,
                0.20,
            ],
        }
    )
    subgroup2 = pd.DataFrame(
        {
            "bio_stratum": np.tile([1, 2, 3], 6),
            "contrast": np.repeat(["A", "B"], 9),
            "index_bin": np.concatenate([np.repeat([5, 10, 15], 3), np.repeat([5, 10, 15], 3)]),
            "proportion_overall": [
                0.10,
                0.20,
                0.20,
                0.00,
                0.00,
                0.20,
                0.10,
                0.00,
                0.10,
                0.30,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
            ],
        }
    )

    return {"subgroup1": subgroup1, "subgroup2": subgroup2}


@pytest.fixture
def apportion_weight_proportions():
    """Weight proportions for apportionment tests"""

    subgroup1 = (
        pd.DataFrame(
            {
                "bio_stratum": np.tile([1, 2, 3], 6),
                "contrast": np.repeat(["A", "B"], 9),
                "index_bin": np.concatenate([np.repeat([5, 10, 15], 3), np.repeat([5, 10, 15], 3)]),
                "extra_bin": np.tile([1, 2], 9),
                "proportions": [
                    0.10,
                    0.20,
                    0.10,
                    0.10,
                    0.10,
                    0.00,
                    0.00,
                    0.00,
                    0.10,
                    0.20,
                    0.00,
                    0.00,
                    0.00,
                    0.20,
                    0.10,
                    0.10,
                    0.30,
                    0.20,
                ],
            }
        )
        .pivot_table(
            index=["contrast", "index_bin", "extra_bin"],
            columns=["bio_stratum"],
            values="proportions",
        )
        .fillna(0.0)
    )
    subgroup2 = pd.DataFrame(
        {
            "bio_stratum": np.tile([1, 2, 3], 6),
            "contrast": np.repeat(["A", "B"], 9),
            "index_bin": np.concatenate([np.repeat([5, 10, 15], 3), np.repeat([5, 10, 15], 3)]),
            "proportions": [
                0.10,
                0.20,
                0.20,
                0.00,
                0.00,
                0.20,
                0.10,
                0.00,
                0.10,
                0.30,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
            ],
        }
    ).pivot_table(index=["contrast", "index_bin"], columns=["bio_stratum"], values="proportions")

    return {"subgroup1": subgroup1, "subgroup2": subgroup2}


@pytest.fixture
def apportion_stratum_weights():
    """Average weights per bin for each stratum for apportionment tests"""

    return pd.DataFrame(
        {
            "bio_stratum": [1, 2, 3],
            "global": [1.0, 2.0, 3.0],
        }
    ).set_index(
        "bio_stratum"
    )["global"]


@pytest.fixture
def apportion_stratum_sigma_bs():
    """Average sigma_bs for each stratum for apportionment tests"""

    return pd.DataFrame(
        {
            "bio_stratum": [1, 2, 3],
            "sigma_bs": [1e-3, 1e-4, 1e-5],
        }
    ).set_index(["bio_stratum"])


@pytest.fixture
def apportion_mesh_with_nasc(
    apportion_mesh,
    apportion_weight_proportions,
    apportion_stratum_weights,
    apportion_stratum_sigma_bs,
):
    """Pre-calculated mesh NASC for apportionment tests"""

    # Create copy of apportion mesh DataFrame
    mesh_data_df = apportion_mesh.copy()

    # Compute biomass
    mesh_data_df["biomass"] = mesh_data_df["biomass_density"] * mesh_data_df["area"]

    # Convert biomass into the various derived population estimates
    apportion.mesh_biomass_to_nasc(
        mesh_data_df=mesh_data_df,
        biodata=apportion_weight_proportions,
        mesh_biodata_link={"mesh_stratum": "bio_stratum"},
        stratum_weights_df=apportion_stratum_weights,
        stratum_sigma_bs_df=apportion_stratum_sigma_bs,
        group_by=["contrast"],
    )

    # Return
    return mesh_data_df


@pytest.fixture
def apportion_biomass_table(apportion_mesh_with_nasc, apportion_weight_proportions):
    """Pre-calculated apportioned biomass tables for apportionment tests"""

    # Create copy of apportion mesh DataFrame
    mesh_data_df = apportion_mesh_with_nasc.copy()

    # Add biomass
    mesh_data_df["biomass"] = mesh_data_df["biomass_density"] * mesh_data_df["area"]

    biomass_tables = apportion.distribute_population_estimates(
        data=mesh_data_df,
        proportions=apportion_weight_proportions,
        variable="biomass",
        group_by=["contrast", "index_bin", "extra_bin"],
        stratify_by=["bio_stratum"],
        data_proportions_link={"mesh_stratum": "bio_stratum"},
    )

    return biomass_tables


@pytest.fixture
def apportion_biomass_table_with_standardized(apportion_biomass_table):
    """Pre-calculate apportioned biomass tables with standardized group for apportionment tests"""

    # Generate copy
    table_copy = apportion_biomass_table.copy()

    # Add standardized table
    table_copy["standardized_subgroup2"] = apportion.distribute_unaged_from_aged(
        population_table=apportion_biomass_table["subgroup2"],
        reference_table=apportion_biomass_table["subgroup1"],
        group_by=["contrast"],
        impute=True,
        impute_variable=["extra_bin"],
    )

    return table_copy


@pytest.fixture
def apportion_combined_biomass_table(apportion_biomass_table_with_standardized):
    """Pre-calculate combined apportioned biomass table for apportionment tests"""

    # Combine tables
    return apportion.sum_population_tables(
        population_table=apportion_biomass_table_with_standardized,
        table_names=["subgroup1", "standardized_subgroup2"],
        table_index=["index_bin"],
        table_columns=["extra_bin", "contrast"],
    )
