import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def plotting_mesh_data():
    np.random.seed(999)

    return pd.DataFrame(
        {
            "longitude": np.linspace(-125, -124, 10),
            "latitude": np.linspace(45, 46, 10),
            "biomass": np.random.rand(10),
        }
    )


@pytest.fixture
def plotting_transect_data():
    np.random.seed(999)
    df = pd.DataFrame(
        {
            "longitude": np.linspace(-125, -124, 10),
            "latitude": np.linspace(45, 46, 10),
            "transect_num": [1] * 5 + [2] * 5,
            "biomass": np.random.rand(10),
        }
    )
    return df


@pytest.fixture
def plotting_heatmap_data():
    np.random.seed(999)
    cols = pd.IntervalIndex.from_breaks([0, 1, 2, 3])
    idx = pd.IntervalIndex.from_breaks([10, 20, 30, 40])
    df = pd.DataFrame(np.random.rand(3, 3), index=idx, columns=cols)
    df.index.name = "length_bin"
    df.columns.name = "age_bin"
    return df
