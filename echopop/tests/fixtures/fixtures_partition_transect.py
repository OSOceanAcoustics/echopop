import pandas as pd
import pytest


@pytest.fixture
def sample_transect_dataset():
    """Create a sample transect dataset for testing."""
    return pd.DataFrame(
        {
            "stratum_ks": [1, 1, 2, 2, 3],
            "transect_num": [1, 1, 2, 2, 3],
            "nasc": [100.0, 150.0, 200.0, 120.0, 180.0],
            "number_density": [50.0, 75.0, 100.0, 60.0, 90.0],
            "abundance": [500.0, 750.0, 1000.0, 600.0, 900.0],
            "abundance_female": [250.0, 375.0, 500.0, 300.0, 450.0],
            "abundance_male": [250.0, 375.0, 500.0, 300.0, 450.0],
            "biomass_density": [25.0, 37.5, 50.0, 30.0, 45.0],
            "biomass": [250.0, 375.0, 500.0, 300.0, 450.0],
            "biomass_female": [125.0, 187.5, 250.0, 150.0, 225.0],
            "biomass_male": [125.0, 187.5, 250.0, 150.0, 225.0],
            "area_interval": [10.0, 10.0, 10.0, 10.0, 10.0],
        }
    )


@pytest.fixture
def age1_nasc_proportions():
    """Create age-1 NASC proportions for testing."""
    return pd.Series([0.1, 0.15, 0.2], index=[1, 2, 3], name="proportion").rename_axis("stratum_ks")


@pytest.fixture
def age1_abundance_proportions():
    """Create age-1 abundance proportions for testing."""
    return pd.Series([0.12, 0.18, 0.22], index=[1, 2, 3], name="proportion").rename_axis(
        "stratum_ks"
    )


@pytest.fixture
def age1_biomass_proportions():
    """Create age-1 biomass proportions for testing."""
    return pd.Series([0.08, 0.14, 0.18], index=[1, 2, 3], name="proportion").rename_axis(
        "stratum_ks"
    )


@pytest.fixture
def number_proportions_partition_dict():
    """Create number proportions dictionary for biology functions."""
    aged_data = pd.DataFrame(
        {
            "stratum_ks": [1, 1, 2, 2, 3, 3],
            "sex": ["female", "male", "female", "male", "female", "male"],
            "proportion_overall": [0.3, 0.2, 0.35, 0.25, 0.28, 0.32],
        }
    ).set_index(["stratum_ks", "sex"])

    return {"aged": aged_data}


@pytest.fixture
def average_weight_df():
    """Create average weight DataFrame for biology functions."""
    df = pd.DataFrame(
        {"female": [0.5, 0.6, 0.55], "male": [0.45, 0.52, 0.48], "all": [0.475, 0.56, 0.515]},
        index=pd.Index([1, 2, 3], name="stratum_ks"),
    )
    df.columns.names = ["sex"]
    return df


@pytest.fixture
def biology_test_dataset():
    """Create a dataset for testing biology functions."""
    return pd.DataFrame(
        {
            "stratum_ks": [1, 1, 2, 2, 3, 3],
            "number_density": [100.0, 120.0, 150.0, 180.0, 110.0, 140.0],
            "number_density_female": [50.0, 60.0, 75.0, 90.0, 55.0, 70.0],
            "number_density_male": [50.0, 60.0, 75.0, 90.0, 55.0, 70.0],
            "area_interval": [5.0, 8.0, 6.0, 7.0, 5.5, 6.5],
            "abundance": [500.0, 960.0, 900.0, 1260.0, 605.0, 910.0],
        }
    )
