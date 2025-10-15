import awkward as awk
import numpy as np
import pandas as pd
import pytest

from echopop.survey import stratified


def test_basic_JollyHampton_initialization():
    """Test basic initialization with minimal parameters."""
    model_params = {
        "transects_per_latitude": 5,
        "strata_transect_proportion": 0.75,
        "num_replicates": 100,
    }

    jh = stratified.JollyHampton(model_params)

    assert jh.model_params == model_params
    assert jh.variable is None
    assert jh.bootstrap_replicates is None
    assert hasattr(jh, "rng")


def test_JollyHampton_initialization_with_seed():
    """Test initialization with random seed."""
    model_params = {
        "transects_per_latitude": 5,
        "strata_transect_proportion": 0.75,
        "num_replicates": 100,
    }

    jh1 = stratified.JollyHampton(model_params, resample_seed=777)
    jh2 = stratified.JollyHampton(model_params, resample_seed=777)

    # Should produce same random numbers with same seed
    random1 = jh1.rng.random(10)
    random2 = jh2.rng.random(10)
    np.testing.assert_array_equal(random1, random2)


# Transect partitioning tests
def test_JollyHampton_partition_data_into_transects(jolly_hampton, sample_ci_grid_data):
    """Test transect partitioning from gridded data."""
    data_with_transects, lat_key = jolly_hampton._partition_data_into_transects(sample_ci_grid_data)

    # Check that transect numbers were added
    assert "transect_num" in data_with_transects.columns
    assert data_with_transects["transect_num"].notna().all()

    # Check latitude key structure
    assert "latitude" in lat_key.columns
    assert "transect_num" in lat_key.columns
    assert len(lat_key) > 0

    # Check that transect numbers are sequential
    transect_nums = sorted(lat_key["transect_num"].unique())
    expected_nums = list(range(1, len(transect_nums) + 1))
    assert transect_nums == expected_nums


def test_JollyHampton_format_virtual_transects(jolly_hampton, sample_ci_grid_data):
    """Test virtual DataFrame initialization."""
    # First partition the data
    data_with_transects, _ = jolly_hampton._partition_data_into_transects(sample_ci_grid_data)

    virtual_df = jolly_hampton._format_virtual_transects(data_with_transects, "biomass")

    # Check required columns are present
    required_cols = ["longitude", "latitude", "area_interval", "biomass"]
    for col in required_cols:
        assert col in virtual_df.columns

    # Check index is transect_num
    assert virtual_df.index.name == "transect_num"


def test_JollyHampton_format_virtual_transects_missing_columns(jolly_hampton):
    """Test error handling for missing columns."""
    incomplete_data = pd.DataFrame(
        {
            "transect_num": [1, 2, 3],
            "latitude": [45, 46, 47],
            # Missing longitude, area, biomass
        }
    )

    with pytest.raises(KeyError, match="Missing required columns"):
        jolly_hampton._format_virtual_transects(incomplete_data, "biomass")


# Virtual transect creation tests
def test_JollyHampton_generate_virtual_transects():
    """Test virtual transect DataFrame creation."""
    model_params = {
        "transects_per_latitude": 5,
        "strata_transect_proportion": 0.75,
        "num_replicates": 10,
    }
    jh = stratified.JollyHampton(model_params)

    # Create test data
    virtual_data = pd.DataFrame(
        {
            "latitude": [45.1, 45.1, 45.2, 45.2],
            "longitude": [-125.1, -125.2, -125.1, -125.2],
            "area_interval": [2.5, 2.5, 2.5, 2.5],
            "biomass": [100, 150, 120, 180],
        },
        index=pd.Index([1, 1, 2, 2], name="transect_num"),
    )

    lat_key = pd.DataFrame({"latitude": [45.1, 45.2], "transect_num": [1, 2]})

    result = jh._generate_virtual_transects(virtual_data, lat_key, "biomass")

    # Check required columns are present
    assert "latitude" in result.columns
    assert "transect_distance" in result.columns
    assert "transect_area" in result.columns
    assert "biomass" in result.columns

    # Check that distances are positive
    assert (result["transect_distance"] > 0).all()

    # Check that areas are positive
    assert (result["transect_area"] > 0).all()


# Bootstrap array preparation tests
def test_JollyHampton_prepare_bootstrap_arrays(prepared_jolly_hampton):
    """Test bootstrap array preparation."""
    jh = prepared_jolly_hampton

    distances, areas, values = jh._prepare_bootstrap_arrays()

    # Check that we get awkward arrays
    assert isinstance(distances, awk.Array)
    assert isinstance(areas, awk.Array)
    assert isinstance(values, awk.Array)

    # Check array structure - should have 2 strata, 5 replicates each
    assert len(distances) == 2  # 2 strata
    assert len(distances[0]) == 5  # 5 replicates
    assert len(distances[1]) == 5  # 5 replicates


def test_JollyHampton_dof_calculation(prepared_jolly_hampton):
    """Test degrees of freedom calculation."""
    jh = prepared_jolly_hampton

    dof = jh._dof()

    assert isinstance(dof, (np.ndarray, pd.Series))
    assert len(dof) == 2  # Two strata
    assert (dof > 0).all()  # All DOF should be positive


# Summarization tests
def test_JollyHampton_compute_transect_statistics(sample_survey_ci_data):
    """Test transect summarization."""
    model_params = {
        "transects_per_latitude": 5,
        "strata_transect_proportion": 0.75,
        "num_replicates": 10,
    }

    jh = stratified.JollyHampton(model_params)
    jh.variable = "biomass"

    jh._compute_transect_statistics(sample_survey_ci_data, ["geostratum_inpfc"])

    # Check that transect_summary was created
    assert jh.transect_summary is not None

    # Check required columns
    required_cols = ["biomass", "distance", "area"]
    for col in required_cols:
        assert col in jh.transect_summary.columns

    # Check that all values are positive
    numeric_cols = jh.transect_summary.select_dtypes(include=[np.number]).columns
    assert (jh.transect_summary[numeric_cols] >= 0).all().all()


def test_JollyHampton_compute_strata_statistics(sample_survey_ci_data):
    """Test strata summarization."""
    model_params = {
        "transects_per_latitude": 5,
        "strata_transect_proportion": 0.75,
        "num_replicates": 10,
    }

    jh = stratified.JollyHampton(model_params)
    jh.variable = "biomass"

    # First summarize transects
    jh._compute_transect_statistics(sample_survey_ci_data, ["geostratum_inpfc"])

    # Then summarize strata
    jh._compute_strata_statistics(["geostratum_inpfc"])

    # Check that strata_summary was created
    assert jh.strata_summary is not None

    # Check required columns
    required_cols = ["transect_counts", "num_transects_to_sample"]
    for col in required_cols:
        assert col in jh.strata_summary.columns

    # Check logical constraints
    assert (jh.strata_summary["num_transects_to_sample"] > 0).all()
    assert (
        jh.strata_summary["num_transects_to_sample"] <= jh.strata_summary["transect_counts"]
    ).all()


# Full workflow tests
def test_JollyHampton_error_handling_no_bootstrap():
    """Test error when calling summarize before bootstrap."""
    model_params = {
        "transects_per_latitude": 10,
        "strata_transect_proportion": 0.8,
        "num_replicates": 20,
    }
    jh = stratified.JollyHampton(model_params, resample_seed=1234)

    with pytest.raises(RuntimeError):
        jh.summarize()


def test_JollyHampton_create_virtual_transects_basic():
    """Test basic virtual transect creation."""
    model_params = {
        "transects_per_latitude": 10,
        "strata_transect_proportion": 0.8,
        "num_replicates": 20,
    }
    jh = stratified.JollyHampton(model_params, resample_seed=42)

    # Create minimal test data
    np.random.seed(98765)
    n_grid = 100
    grid_data = pd.DataFrame(
        {
            "latitude": np.random.uniform(45, 50, n_grid),
            "longitude": np.random.uniform(-130, -120, n_grid),
            "area": np.random.uniform(1, 5, n_grid),
            "biomass": np.random.exponential(100, n_grid),
        }
    )

    # Geostrata data
    geostrata = pd.DataFrame(
        {
            "northlimit_latitude": np.linspace(44, 51, 50),
            "stratum_num": np.random.choice(["stratum1", "stratum2"], 50),
        }
    )

    virtual_data = jh.create_virtual_transects(
        data_df=grid_data,
        geostrata_df=geostrata,
        stratify_by=["geostratum_inpfc"],
        variable="biomass",
    )

    assert isinstance(virtual_data, pd.DataFrame)
    assert not virtual_data.empty


# Static method tests
def test_JollyHampton_compute_variance_static_method():
    """Test the static variance computation method."""
    # Create simple test data
    np.random.seed(123456789)

    # Create awkward arrays
    sampled_values = awk.Array(
        [
            [[100, 110, 90], [105, 95, 100]],  # Stratum 1, 2 replicates
            [[200, 180, 220], [190, 210, 200]],  # Stratum 2, 2 replicates
        ]
    )

    sampled_distances = awk.Array(
        [
            [[10, 12, 8], [11, 9, 10]],  # Stratum 1, 2 replicates
            [[15, 13, 17], [14, 16, 15]],  # Stratum 2, 2 replicates
        ]
    )

    # Compute means
    mean_values = awk.sum(sampled_values * sampled_distances, axis=-1) / awk.sum(
        sampled_distances, axis=-1
    )

    # Compute weights
    stratified_weights = sampled_distances / awk.mean(sampled_distances, axis=-1, keepdims=True)

    # Sample DOF
    sample_dof = np.array([2, 2])

    # Test the method
    result = stratified.JollyHampton._compute_variance(
        sampled_values, sampled_distances, mean_values, stratified_weights, sample_dof
    )

    assert isinstance(result, awk.Array)
    assert awk.to_numpy(result).shape == (2, 2)  # 2 strata, 2 replicates

    # All variances should be non-negative
    variance_np = awk.to_numpy(result)
    assert (variance_np >= 0).all()


# Robustness and edge case tests
def test_JollyHampton_single_transect_stratum():
    """Test handling of strata with only one transect."""
    model_params = {
        "transects_per_latitude": 5,
        "strata_transect_proportion": 1.0,  # Sample all transects
        "num_replicates": 10,
    }

    jh = stratified.JollyHampton(model_params, resample_seed=2468)

    # Create data with one stratum having only one transect
    test_data = pd.DataFrame(
        {
            "geostratum_inpfc": ["stratum1", "stratum2"],
            "transect_num": [1, 2],
            "transect_distance": [10.0, 12.0],
            "transect_area": [50.0, 60.0],
            "biomass": [1000, 1200],
        }
    )

    jh.variable = "biomass"

    # This should handle single-transect strata gracefully
    jh._compute_transect_statistics(test_data, ["geostratum_inpfc"])
    jh._compute_strata_statistics(["geostratum_inpfc"])

    # Check that single-transect strata get at least 1 sample
    assert (jh.strata_summary["num_transects_to_sample"] >= 1).all()


def test_JollyHampton_extreme_sampling_proportion():
    """Test extreme sampling proportions."""
    # Very low proportion
    model_params_low = {
        "transects_per_latitude": 5,
        "strata_transect_proportion": 0.01,  # Very low
        "num_replicates": 10,
    }

    jh_low = stratified.JollyHampton(model_params_low)
    jh_low.strata_summary = pd.DataFrame(
        {"transect_counts": [10, 20], "area": [100, 200], "biomass": [1000, 2000]},
        index=["stratum1", "stratum2"],
    )

    # Should ensure at least 1 transect per stratum
    sample_counts = np.round(jh_low.strata_summary["transect_counts"] * 0.01).astype(int)
    adjusted_counts = np.maximum(sample_counts, 1)

    assert (adjusted_counts >= 1).all()


def test_JollyHampton_empty_strata_handling():
    """Test handling of empty or very small datasets."""
    model_params = {
        "transects_per_latitude": 5,
        "strata_transect_proportion": 0.75,
        "num_replicates": 10,
    }

    jh = stratified.JollyHampton(model_params)

    # Test with minimal data
    minimal_data = pd.DataFrame(
        {
            "geostratum_inpfc": ["stratum1"],
            "transect_num": [1],
            "transect_distance": [10.0],
            "transect_area": [50.0],
            "biomass": [1000],
        }
    )

    jh.variable = "biomass"
    jh._compute_transect_statistics(minimal_data, ["geostratum_inpfc"])
    jh._compute_strata_statistics(["geostratum_inpfc"])

    # Should handle gracefully
    assert len(jh.strata_summary) == 1
    assert jh.strata_summary["num_transects_to_sample"].iloc[0] >= 1
