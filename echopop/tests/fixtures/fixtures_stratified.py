import numpy as np
import pandas as pd
import pytest

from echopop.survey import stratified


# Test data fixtures
@pytest.fixture
def sample_ci_grid_data():
    """Create sample gridded data for testing."""
    np.random.seed(1234)
    n_points = 1000

    data = pd.DataFrame(
        {
            "latitude": np.random.uniform(45, 50, n_points),
            "longitude": np.random.uniform(-130, -120, n_points),
            "area": np.random.uniform(1, 5, n_points),
            "biomass": np.random.exponential(100, n_points),
        }
    )

    return data


@pytest.fixture
def jolly_hampton():
    """Create JollyHampton instance for testing."""
    model_params = {
        "transects_per_latitude": 5,
        "strata_transect_proportion": 0.75,
        "num_replicates": 10,
    }
    return stratified.JollyHampton(model_params, resample_seed=1234)


@pytest.fixture
def sample_survey_ci_data():
    """Create sample survey data for testing."""
    np.random.seed(1234)

    # Create sample data with multiple strata and transects
    n_points = 200
    data = pd.DataFrame(
        {
            "geostratum_inpfc": np.random.choice(["stratum1", "stratum2"], n_points),
            "transect_num": np.random.randint(1, 21, n_points),
            "longitude": np.random.uniform(-130, -120, n_points),
            "latitude": np.random.uniform(45, 50, n_points),
            "transect_distance": np.random.uniform(5, 15, n_points),
            "transect_area": np.random.uniform(25, 75, n_points),
            "biomass": np.random.exponential(100, n_points),
        }
    )

    return data


@pytest.fixture
def prepared_jolly_hampton():
    """Create JollyHampton with mock data for bootstrap testing."""
    model_params = {
        "transects_per_latitude": 5,
        "strata_transect_proportion": 0.8,
        "num_replicates": 5,  # Small for testing
    }

    jh = stratified.JollyHampton(model_params, resample_seed=1234)
    jh.variable = "biomass"

    # Create mock transect summary
    strata_idx = pd.Index(
        ["stratum1", "stratum1", "stratum1", "stratum2", "stratum2", "stratum2"],
        name="geostratum_inpfc",
    )
    transect_idx = pd.Index([1, 2, 3, 4, 5, 6], name="transect_num")

    multi_idx = pd.MultiIndex.from_arrays([strata_idx, transect_idx])
    transect_summary = pd.DataFrame(
        {
            "distance": [10.0, 12.0, 8.0, 15.0, 9.0, 11.0],
            "area": [50.0, 60.0, 40.0, 75.0, 45.0, 55.0],
            "biomass": [1000, 1200, 800, 1500, 900, 1100],
        },
        index=multi_idx,
    )

    jh.transect_summary = transect_summary

    # Create mock strata summary
    strata_summary = pd.DataFrame(
        {
            "transect_counts": [3, 3],
            "num_transects_to_sample": [2, 2],
            "distance": [30.0, 35.0],
            "area": [150.0, 175.0],
            "biomass": [3000, 3500],
        },
        index=pd.Index(["stratum1", "stratum2"], name="geostratum_inpfc"),
    )

    jh.strata_summary = strata_summary

    return jh
