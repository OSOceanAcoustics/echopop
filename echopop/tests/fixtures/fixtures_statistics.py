import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def sample_ci_data():
    """Generate sample bootstrap data for testing."""
    np.random.seed(999)
    samples = np.random.normal(100, 10, 1000)
    population_stat = 95.0
    interval = np.array([0.025, 0.975])
    return samples, population_stat, interval


@pytest.fixture
def skewed_data():
    """Generate highly skewed data that should cause BC/BCa to fail."""
    samples = np.array([100] * 1000)  # All same value
    population_stat = 95.0
    interval = np.array([0.025, 0.975])
    return samples, population_stat, interval


@pytest.fixture
def bootstrap_dataframe():
    """Create sample bootstrap DataFrame."""
    np.random.seed(999)
    n_replicates = 100

    data = {
        "biomass": np.random.normal(1000, 100, n_replicates),
        "density": np.random.normal(50, 5, n_replicates),
        "cv": np.random.normal(0.15, 0.02, n_replicates),
    }

    return pd.DataFrame(data)


@pytest.fixture
def population_dataframe():
    """Create sample population DataFrame."""
    data = {"biomass": [950], "density": [48], "cv": [0.14]}
    return pd.DataFrame(data)
