import numpy as np
import pandas as pd

from echopop.geostatistics import transform_coordinates


# ==================================================================================================
# Test transform_coordinates
# --------------------------
def test_transform_coordinates_basic(sample_coordinates_df):
    """Test basic coordinate transformation."""
    result_df, delta_x, delta_y = transform_coordinates(
        data=sample_coordinates_df,
        x_offset=0.0,
        y_offset=0.0,
    )

    # Check that x and y columns are created
    assert "x" in result_df.columns
    assert "y" in result_df.columns

    # Check that delta values are returned
    assert delta_x is not None
    assert delta_y is not None
    assert isinstance(delta_x, (int, float))
    assert isinstance(delta_y, (int, float))

    # Check that values are numeric
    assert result_df["x"].dtype in [np.float64, np.float32]
    assert result_df["y"].dtype in [np.float64, np.float32]

    # Check that transformation maintains relative ordering
    assert len(result_df) == len(sample_coordinates_df)


def test_transform_coordinates_with_offsets(sample_coordinates_df):
    """Test coordinate transformation with offsets."""
    x_offset = -124.0
    y_offset = 46.0

    result_df, delta_x, delta_y = transform_coordinates(
        data=sample_coordinates_df,
        x_offset=x_offset,
        y_offset=y_offset,
    )

    # Check that x and y columns are created
    assert "x" in result_df.columns
    assert "y" in result_df.columns

    # Check that delta values are positive
    assert delta_x > 0
    assert delta_y > 0

    # Check that transformed coordinates are reasonable (not strict [-1, 1] range)
    # The function applies cosine transformation and scaling, not normalization
    assert np.isfinite(result_df["x"]).all()
    assert np.isfinite(result_df["y"]).all()
    assert result_df["x"].max() > result_df["x"].min()
    assert result_df["y"].max() > result_df["y"].min()


def test_transform_coordinates_with_reference(sample_coordinates_df, sample_reference_df):
    """Test coordinate transformation with reference DataFrame."""
    result_df, delta_x, delta_y = transform_coordinates(
        data=sample_coordinates_df,
        reference=sample_reference_df,
        x_offset=-124.0,
        y_offset=46.0,
    )

    # Check that x and y columns are created
    assert "x" in result_df.columns
    assert "y" in result_df.columns

    # Check that delta values are returned
    assert delta_x is not None
    assert delta_y is not None


def test_transform_coordinates_with_deltas(sample_coordinates_df):
    """Test coordinate transformation with provided deltas."""
    delta_x_input = 2.0
    delta_y_input = 1.5

    result_df, delta_x, delta_y = transform_coordinates(
        data=sample_coordinates_df,
        x_offset=0.0,
        y_offset=0.0,
        delta_x=delta_x_input,
        delta_y=delta_y_input,
    )

    # Check that provided deltas are returned
    assert delta_x == delta_x_input
    assert delta_y == delta_y_input


def test_transform_coordinates_custom_column_names():
    """Test coordinate transformation with custom coordinate names."""
    data = pd.DataFrame(
        {
            "lon": [-124.5, -124.3, -124.1],
            "lat": [46.2, 46.4, 46.6],
            "value": [10.5, 15.2, 12.8],
        }
    )

    result_df, delta_x, delta_y = transform_coordinates(
        data=data,
        x_offset=0.0,
        y_offset=0.0,
        coordinate_names=("lon", "lat"),
    )

    # Check that x and y columns are created
    assert "x" in result_df.columns
    assert "y" in result_df.columns

    # Check original columns are preserved
    assert "lon" in result_df.columns
    assert "lat" in result_df.columns
    assert "value" in result_df.columns
