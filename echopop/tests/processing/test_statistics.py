import numpy as np
import pandas as pd
import pytest

from echopop.survey import statistics as stats


def test_percentile_method(sample_ci_data):
    """Test percentile bootstrap method."""
    samples, population_stat, interval = sample_ci_data
    result = stats.percentile(samples, interval)

    # Assert shape and typing
    assert len(result) == 2
    assert result[0] < result[1]  # Lower bound < upper bound
    assert np.all(np.isfinite(result))  # No NaN/inf values

    # Assert values
    assert np.allclose(result, np.array([79.6233578, 120.4478217]))


def test_normal_method(sample_ci_data):
    """Test normal approximation method."""
    samples, population_stat, interval = sample_ci_data
    result = stats.normal(samples, interval)

    # Assert shape and typing
    assert len(result) == 2
    assert result[0] < result[1]
    assert np.all(np.isfinite(result))

    # Assert values
    assert np.allclose(result, np.array([79.55786428, 120.06676714]))


def test_empirical_method(sample_ci_data):
    """Test empirical bootstrap method."""
    samples, population_stat, interval = sample_ci_data
    result = stats.empirical(samples, interval, np.array([population_stat]))

    # Assert shape and typing
    assert len(result) == 2
    assert result[0] < result[1]
    assert np.all(np.isfinite(result))

    # Assert values
    assert np.allclose(result, np.array([84.43567348, 125.26013739]))


def test_student_standard_method(sample_ci_data):
    """Test studentized t-distribution method."""
    samples, population_stat, interval = sample_ci_data
    result = stats.student_standard(samples, interval)

    # Assert shape and typing
    assert len(result) == 2
    assert result[0] < result[1]
    assert np.all(np.isfinite(result))

    # Assert values
    assert np.allclose(result, np.array([79.53329522, 120.0913362]))


def test_student_jackknife_method(sample_ci_data):
    """Test jackknife studentized method."""
    samples, population_stat, interval = sample_ci_data
    result = stats.student_jackknife(samples, interval)

    # Assert shape and typing
    assert len(result) == 2
    assert result[0] < result[1]
    assert np.all(np.isfinite(result))

    # Assert values
    assert np.allclose(result, np.array([79.57453103, 120.49950457]))


def test_bc_method_normal_case(sample_ci_data):
    """Test bias-corrected method with normal data."""
    samples, population_stat, interval = sample_ci_data
    result = stats.bc(samples, interval, population_stat, "test_estimator")

    # Assert shape and typing
    assert len(result) == 2
    assert result[0] < result[1]
    assert np.all(np.isfinite(result))

    # Assert values
    assert np.allclose(result, np.array([67.86816843, 109.43886671]))


def test_bc_method_skewed_failure(skewed_data):
    """Test that BC method raises error for highly skewed data."""
    samples, population_stat, interval = skewed_data

    with pytest.raises(ValueError, match="Significant skewness detected"):
        stats.bc(samples, interval, population_stat, "test_estimator")


def test_bca_method_normal_case(sample_ci_data):
    """Test bias-corrected and accelerated method with normal data."""
    samples, population_stat, interval = sample_ci_data
    result = stats.bca(samples, interval, population_stat, "test_estimator")

    # Assert shape and typing
    assert len(result) == 2
    assert result[0] < result[1]
    assert np.all(np.isfinite(result))

    # Assert values
    assert np.allclose(result, np.array([68.32433457, 109.85050844]))


def test_bca_method_skewed_failure(skewed_data):
    """Test that BCa method raises error for highly skewed data."""
    samples, population_stat, interval = skewed_data

    with pytest.raises(ValueError, match="Significant skewness detected"):
        stats.bca(samples, interval, population_stat, "test_estimator")


def test_small_sample_jackknife():
    """Test jackknife method with small sample size."""
    np.random.seed(500)
    samples = np.random.normal(100, 10, 10)  # Small sample
    interval = np.array([0.025, 0.975])

    result = stats.student_jackknife(samples, interval)

    # Assert shape and  typing
    assert len(result) == 2
    assert np.all(np.isfinite(result))

    # Assert values
    assert np.allclose(result, np.array([76.98222532, 124.62199098]))


def test_zero_variance_handling():
    """Test handling of zero variance in jackknife method."""
    samples = np.array([100.0] * 20 + [100.1])  # Nearly constant
    interval = np.array([0.025, 0.975])

    # Should not raise error, should handle gracefully
    result = stats.student_jackknife(samples, interval)
    assert len(result) == 2
    assert np.all(np.isfinite(result))


# Confidence interval function tests
def test_basic_confidence_interval(bootstrap_dataframe, population_dataframe):
    """Test basic confidence interval calculation."""
    result = stats.confidence_interval(
        bootstrap_samples=bootstrap_dataframe,
        population_values=population_dataframe,
        ci_method="percentile",
        ci_percentile=0.95,
    )

    # Check structure
    assert isinstance(result, pd.DataFrame)
    assert result.shape == (1, 12)
    assert "biomass" in result.columns
    assert "density" in result.columns
    assert "cv" in result.columns


def test_all_ci_methods(bootstrap_dataframe, population_dataframe):
    """Test that all CI methods work without errors."""
    methods = ["percentile", "normal", "empirical", "t", "t-jackknife"]

    for method in methods:
        result = stats.confidence_interval(
            bootstrap_samples=bootstrap_dataframe,
            population_values=population_dataframe,
            ci_method=method,
            ci_percentile=0.95,
        )

        assert isinstance(result, pd.DataFrame)
        assert not result.isnull().any().any()  # No NaN values


def test_invalid_ci_method(bootstrap_dataframe, population_dataframe):
    """Test error handling for invalid CI method."""
    with pytest.raises(ValueError, match="ci_method 'invalid' not supported"):
        stats.confidence_interval(
            bootstrap_samples=bootstrap_dataframe,
            population_values=population_dataframe,
            ci_method="invalid",
        )


def test_invalid_population_values(bootstrap_dataframe):
    """Test error handling for invalid population values type."""
    with pytest.raises(TypeError, match="must be a pandas.Series or pandas.DataFrame"):
        stats.confidence_interval(
            bootstrap_samples=bootstrap_dataframe,
            population_values=np.array([1, 2, 3]),  # Wrong type
            ci_method="percentile",
        )


def test_mismatched_columns(bootstrap_dataframe):
    """Test handling of mismatched column names."""
    # Should still work but may produce warnings or handle gracefully
    # This tests the robustness of the column matching logic
    try:
        result = stats.confidence_interval(
            bootstrap_samples=bootstrap_dataframe[["biomass"]],  # Only one column
            population_values=pd.DataFrame({"biomass": [950]}),
            ci_method="percentile",
        )
        assert isinstance(result, pd.DataFrame)
    except (KeyError, IndexError):
        # This is acceptable behavior for mismatched columns
        pass


def test_different_ci_percentiles(bootstrap_dataframe, population_dataframe):
    """Test different confidence interval percentiles."""
    percentiles = [0.90, 0.95, 0.99]

    for pct in percentiles:
        result = stats.confidence_interval(
            bootstrap_samples=bootstrap_dataframe,
            population_values=population_dataframe,
            ci_method="percentile",
            ci_percentile=pct,
        )

        assert isinstance(result, pd.DataFrame)
        assert not result.isnull().any().any()


def test_single_column_input():
    """Test confidence interval with single column input."""
    np.random.seed(12345)
    bootstrap_data = pd.DataFrame({"biomass": np.random.normal(1000, 100, 100)})
    population_data = pd.DataFrame({"biomass": [950]})

    result = stats.confidence_interval(
        bootstrap_samples=bootstrap_data,
        population_values=population_data,
        ci_method="percentile",
        ci_percentile=0.95,
    )

    assert isinstance(result, pd.DataFrame)
    assert "biomass" in result.columns


def test_bootstrap_methods_dictionary():
    """Test that all methods in dictionary are callable."""
    for method_name, method_func in stats.BOOTSTRAP_CI_METHODS.items():
        assert callable(method_func), f"Method {method_name} is not callable"

    # Test expected methods are present
    expected_methods = ["bc", "bca", "empirical", "normal", "percentile", "t", "t-jackknife"]
    for method in expected_methods:
        assert method in stats.BOOTSTRAP_CI_METHODS, f"Method {method} missing from dictionary"


# Edge case tests
def test_empty_samples():
    """Test handling of empty sample arrays."""
    with pytest.raises((ValueError, IndexError)):
        stats.percentile(np.array([]), np.array([0.025, 0.975]))


def test_single_value_samples():
    """Test handling of single-value samples."""
    samples = np.array([100.0])
    interval = np.array([0.025, 0.975])

    # Most methods should handle this, though results may be degenerate
    result = stats.percentile(samples, interval)
    assert len(result) == 2

    # Assert values
    assert all(result == np.array([100.0]))


def test_extreme_intervals():
    """Test extreme confidence interval bounds."""
    np.random.seed(666)
    samples = np.random.normal(100, 10, 1000)

    # Very wide interval
    wide_interval = np.array([0.001, 0.999])
    result = stats.percentile(samples, wide_interval)

    # Assert shape
    assert result[0] < result[1]

    # Assert value
    assert np.allclose(result, np.array([68.95992362, 133.42305677]))

    # Very narrow interval
    narrow_interval = np.array([0.49, 0.51])
    result = stats.percentile(samples, narrow_interval)

    # Assert shape
    assert result[0] <= result[1]  # May be equal for very narrow intervals

    # Assert value
    assert np.allclose(result, np.array([99.41458156, 99.66322173]))


def test_large_sample_size():
    """Test performance with large sample sizes."""
    np.random.seed(666)
    samples = np.random.normal(100, 10, 10000)  # Large sample
    interval = np.array([0.025, 0.975])

    result = stats.percentile(samples, interval)

    # Assert shape and typing
    assert len(result) == 2
    assert np.all(np.isfinite(result))

    # Assert values
    assert np.allclose(result, np.array([80.20937026, 119.3864021]))
