import numpy as np
import pandas as pd
import pytest
from lmfit import Parameters

from echopop.nwfsc_feat.geostatistics import Geostats
from echopop.nwfsc_feat.projection import reproject_dataset


# ==================================================================================================
# Test Geostat class initialization
# ---------------------------------
def test_geostats_init_basic(sample_coordinates_df_with_biomass, sample_mesh_data_with_lonlat):
    """Test basic Geostat initialization."""

    # Create sample parameter dictionaries
    variogram_params = {
        "model": ["bessel", "exponential"],
        "n_lags": 30,
        "lag_resolution": 0.5,
    }

    kriging_params = {
        "kmin": 3,
        "kmax": 8,
        "search_radius": 5.0,
    }

    geostat = Geostats(
        data_df=sample_coordinates_df_with_biomass,
        mesh_df=sample_mesh_data_with_lonlat,
        variogram_params=variogram_params,
        kriging_params=kriging_params,
    )

    # Check that attributes are set correctly
    assert geostat.data_df.equals(sample_coordinates_df_with_biomass)
    assert geostat.mesh_df.equals(sample_mesh_data_with_lonlat)
    assert geostat.variogram_params == variogram_params
    assert geostat.kriging_params == kriging_params
    assert geostat.projection == "epsg:4326"
    assert geostat.coordinates == ("longitude", "latitude")
    assert geostat.projection_coordinates == ("longitude", "latitude")


def test_geostats_init_custom_params(sample_coordinates_df_xy, sample_mesh_data_xy):
    """Test Geostats initialization with custom parameters."""

    variogram_params = {"model": ["exponential"], "n_lags": 20}
    kriging_params = {"kmin": 5, "kmax": 10}

    geostat = Geostats(
        data_df=sample_coordinates_df_xy,
        mesh_df=sample_mesh_data_xy,
        variogram_params=variogram_params,
        kriging_params=kriging_params,
        projection="epsg:32610",
        coordinate_names=("x", "y"),
    )

    assert geostat.projection == "epsg:32610"
    assert geostat.coordinates == ("x", "y")
    assert geostat.projection_coordinates == ("x", "y")


# ==================================================================================================
# Test coordinate projection methods
# ---------------------------------
def test_geostats_project_coordinates_normalize(
    sample_coordinates_df_with_biomass, sample_mesh_data_with_lonlat
):
    """Test coordinate projection with normalization."""

    geostat = Geostats(
        data_df=sample_coordinates_df_with_biomass,
        mesh_df=sample_mesh_data_with_lonlat,
        variogram_params={},
        kriging_params={},
    )

    # Project coordinates with normalization
    geostat.project_coordinates(x_offset=-124.0, y_offset=47.0, normalize=True)

    # Check that coordinates were projected
    assert "x" in geostat.data_df.columns
    assert "y" in geostat.data_df.columns
    assert "x" in geostat.mesh_df.columns
    assert "y" in geostat.mesh_df.columns
    assert geostat.projection_coordinates == ("x", "y")


def test_geostats_project_coordinates_crs_error(
    sample_coordinates_df_with_biomass, sample_mesh_data_with_lonlat
):
    """Test that CRS projection raises error when crs_out not specified."""

    geostat = Geostats(
        data_df=sample_coordinates_df_with_biomass,
        mesh_df=sample_mesh_data_with_lonlat,
        variogram_params={},
        kriging_params={},
    )

    with pytest.raises(ValueError, match="crs_out must be specified"):
        geostat.project_coordinates(normalize=False)


# ==================================================================================================
# Test mesh cropping
# ---------------------------------
def test_geostats_crop_mesh(sample_coordinates_df_with_biomass, sample_mesh_data_with_lonlat):
    """Test mesh cropping functionality."""

    geostat = Geostats(
        data_df=sample_coordinates_df_with_biomass,
        mesh_df=sample_mesh_data_with_lonlat,
        variogram_params={},
        kriging_params={},
    )

    # First project coordinates
    geostat.project_coordinates(normalize=True)

    # Define a simple crop function that returns half the mesh
    def simple_crop_function(data_df, mesh_df, **kwargs):
        return mesh_df.iloc[: len(mesh_df) // 2]

    original_mesh_size = len(geostat.mesh_df)
    geostat.crop_mesh(simple_crop_function)

    # Check that mesh was cropped
    assert len(geostat.mesh_df) == original_mesh_size // 2


def test_geostats_crop_mesh_with_tuple_return(
    sample_coordinates_df_with_biomass, sample_mesh_data_with_lonlat
):
    """Test mesh cropping with function that returns tuple."""

    geostat = Geostats(
        data_df=sample_coordinates_df_with_biomass,
        mesh_df=sample_mesh_data_with_lonlat,
        variogram_params={},
        kriging_params={},
    )

    geostat.project_coordinates(normalize=True)

    # Define crop function that returns tuple
    def crop_function_with_tuple(data_df, mesh_df, **kwargs):
        cropped_mesh = mesh_df.iloc[: len(mesh_df) // 2]
        return cropped_mesh, "additional_info"

    original_mesh_size = len(geostat.mesh_df)
    geostat.crop_mesh(crop_function_with_tuple)

    # Check that mesh was cropped (first element of tuple)
    assert len(geostat.mesh_df) == original_mesh_size // 2


# ==================================================================================================
# Test empirical variogram calculation
# ---------------------------------
def test_geostats_calculate_empirical_variogram(
    sample_coordinates_df_with_biomass, sample_mesh_data_with_lonlat
):
    """Test empirical variogram calculation."""

    geostat = Geostats(
        data_df=sample_coordinates_df_with_biomass,
        mesh_df=sample_mesh_data_with_lonlat,
        variogram_params={"n_lags": 10, "lag_resolution": 0.02},
        kriging_params={},
    )

    # Project coordinates first
    geostat.project_coordinates(normalize=True)

    # Calculate empirical variogram
    geostat.calculate_empirical_variogram(
        variable="biomass_density", azimuth_filter=True, force_lag_zero=True
    )

    # Check that variogram attributes are set
    assert hasattr(geostat, "lags")
    assert hasattr(geostat, "gamma")
    assert hasattr(geostat, "lag_counts")
    assert hasattr(geostat, "lag_covariance")
    assert hasattr(geostat, "variable")
    assert geostat.variable == "biomass_density"

    # Check that arrays have expected properties
    assert isinstance(geostat.lags, np.ndarray)
    assert isinstance(geostat.gamma, np.ndarray)
    assert len(geostat.lags) == len(geostat.gamma)


# ==================================================================================================
# Test variogram model fitting
# ---------------------------------
def test_geostats_fit_variogram_model(sample_coordinates_df_latlon):
    """Test variogram model fitting."""

    # Create a simple mesh for testing
    mesh_df = sample_coordinates_df_latlon.copy()[["longitude", "latitude"]].iloc[:25]

    geostat = Geostats(
        data_df=sample_coordinates_df_latlon,
        mesh_df=mesh_df,
        variogram_params={"model": "exponential", "n_lags": 100, "lag_resolution": 0.005},
        kriging_params={},
    )

    # Calculate empirical variogram (no projection needed since coordinates are already x,y)
    geostat.calculate_empirical_variogram(variable="biomass_density")

    # Create parameters for fitting
    params = Parameters()
    params.add("nugget", value=0.1, min=0)
    params.add("sill", value=1.0, min=0)
    params.add("correlation_range", value=0.5, min=0)

    # Fit variogram model
    geostat.fit_variogram_model(params)

    # Check that fitting attributes are set
    assert hasattr(geostat, "best_fit_variogram_params")
    assert hasattr(geostat, "variogram_fit_initial")
    assert hasattr(geostat, "variogram_fit_optimized")
    assert isinstance(geostat.best_fit_variogram_params, dict)

    # Compare the values
    assert np.isclose(geostat.variogram_fit_initial, 0.004606002049)
    assert np.isclose(geostat.variogram_fit_optimized, 0.001522485501)


# ==================================================================================================
# Test kriging
# ---------------------------------
def test_geostats_krige(sample_coordinates_df_latlon, sample_mesh_data_latlon):
    """Test kriging interpolation."""

    # Create a simple mesh for testing
    geostat = Geostats(
        data_df=sample_coordinates_df_latlon,
        mesh_df=sample_mesh_data_latlon,
        variogram_params={"model": "exponential", "n_lags": 5, "lag_resolution": 0.5},
        kriging_params={"k_min": 3, "k_max": 8, "search_radius": 1.5, "anisotropy": 0.001},
    )

    # Standardize coordinates
    geostat.project_coordinates()

    # Complete workflow (no projection needed since coordinates are already x,y)
    geostat.calculate_empirical_variogram(variable="biomass_density")

    # Fit variogram
    params = Parameters()
    params.add("nugget", value=0.1, min=0)
    params.add("sill", value=1.0, min=0)
    params.add("correlation_range", value=0.5, min=0)
    geostat.fit_variogram_model(params)

    # Perform kriging
    results = geostat.krige(default_mesh_cell_area=1.0)

    # Check results
    assert isinstance(results, pd.DataFrame)
    assert "biomass_density" in results.columns
    assert len(results) > 0
    assert hasattr(geostat, "survey_cv")

    # Check values
    assert np.isclose(geostat.survey_cv, 0.009740188342)


# ==================================================================================================
# Test full workflow integration
# ---------------------------------
def test_geostats_full_workflow(sample_coordinates_df_latlon, sample_mesh_data_latlon):
    """Test complete geostatistics workflow."""

    # Create a simple mesh for testing
    geostat = Geostats(
        data_df=sample_coordinates_df_latlon,
        mesh_df=sample_mesh_data_latlon,
        variogram_params={"model": "exponential", "n_lags": 5, "lag_resolution": 0.5},
        kriging_params={"k_min": 3, "k_max": 8, "search_radius": 1.5, "anisotropy": 0.2},
    )

    # Step 1: Crop mesh (simple crop)
    def simple_crop(data_df, mesh_df, **kwargs):
        return mesh_df.iloc[: len(mesh_df) // 2]

    geostat.crop_mesh(simple_crop)

    # Step 2: Project new coordinates
    geostat.project_coordinates()

    # Step 3: Calculate empirical variogram
    geostat.calculate_empirical_variogram(variable="biomass_density")

    # Step 4: Fit variogram model
    params = Parameters()
    params.add("nugget", value=0.1, min=0)
    params.add("sill", value=1.0, min=0)
    params.add("correlation_range", value=0.5, min=0)
    geostat.fit_variogram_model(params)

    # Step 5: Perform kriging
    results = geostat.krige(default_mesh_cell_area=1.0)

    # Verify complete workflow
    assert isinstance(results, pd.DataFrame)
    assert "biomass_density" in results.columns
    assert geostat.variable == "biomass_density"
    assert geostat.projection_coordinates == ("x", "y")
    assert len(geostat.mesh_df) == len(sample_mesh_data_latlon) // 2  # cropped

    # Check values
    assert np.isclose(geostat.survey_cv, 0.014007013796)


# ==================================================================================================
# Test reproject_dataset function
# ---------------------------------
def test_reproject_dataset_basic():
    """Test basic reproject_dataset functionality."""

    # Create sample DataFrame
    df = pd.DataFrame(
        {
            "longitude": [-124.0, -123.8, -123.6],
            "latitude": [47.0, 47.2, 47.4],
            "value": [10, 20, 30],
        }
    )

    # Test with CRS projection (this is what the function actually does)
    result = reproject_dataset(
        data_df=df,
        crs_out="epsg:32610",  # UTM Zone 10N
        coordinate_names=("longitude", "latitude"),
        projection="epsg:4326",
    )

    # Check that x, y columns were added
    assert "x" in result.columns
    assert "y" in result.columns
    assert "value" in result.columns
    assert len(result) == len(df)

    # Check that x, y values are different from original (projected)
    assert not result["x"].equals(df["longitude"])
    assert not result["y"].equals(df["latitude"])


def test_reproject_dataset_custom_coordinates():
    """Test reproject_dataset with custom coordinate names."""

    df = pd.DataFrame({"x_coord": [-124.0, -123.8], "y_coord": [47.0, 47.2], "value": [10, 20]})

    result = reproject_dataset(
        data_df=df,
        crs_out="epsg:32610",
        coordinate_names=("x_coord", "y_coord"),
        projection="epsg:4326",
    )

    assert "x" in result.columns
    assert "y" in result.columns
    assert "value" in result.columns


def test_reproject_dataset_different_input_projection():
    """Test reproject_dataset with different input projection."""

    # Create sample DataFrame with UTM coordinates
    df = pd.DataFrame(
        {"x": [500000, 505000, 510000], "y": [5200000, 5205000, 5210000], "value": [10, 20, 30]}
    )

    # Project from UTM back to WGS84
    result = reproject_dataset(
        data_df=df,
        crs_out="epsg:4326",
        coordinate_names=("x", "y"),
        projection="epsg:32610",  # UTM Zone 10N input
    )

    # Check that x, y columns were added (should be lat/lon values)
    assert "x" in result.columns
    assert "y" in result.columns
    assert "value" in result.columns
    assert len(result) == len(df)


def test_reproject_dataset_preserves_other_columns():
    """Test that reproject_dataset preserves all other columns."""

    df = pd.DataFrame(
        {
            "longitude": [-124.0, -123.8],
            "latitude": [47.0, 47.2],
            "biomass_density": [10.5, 12.3],
            "transect_num": [1, 2],
            "depth": [100, 150],
        }
    )

    result = reproject_dataset(
        data_df=df, crs_out="epsg:32610", coordinate_names=("longitude", "latitude")
    )

    # Check that all original columns are preserved
    assert "biomass_density" in result.columns
    assert "transect_num" in result.columns
    assert "depth" in result.columns
    assert "longitude" in result.columns
    assert "latitude" in result.columns

    # Check that new x, y columns were added
    assert "x" in result.columns
    assert "y" in result.columns
