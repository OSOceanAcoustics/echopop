import numpy as np
import pandas as pd

from echopop.workflows.nwfsc_feat import functions as ingest_nasc


# ==================================================================================================
# Test ASFC NASC to FEAT NASC DataFrame conversion
# ------------------------------------------------
def test_afsc_nasc_to_feat_conversion():
    """
    Test the conversion of ASFC NASC data to FEAT NASC DataFrame format.
    """

    # Mock DataFrame
    mock_df = pd.DataFrame(
        {
            "transect_num": [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11],
            "distance": np.repeat(np.linspace(1, 11, 11), 2) + np.tile(np.array([0.00, 0.25]), 11),
            "transect_spacing": np.concatenate([np.repeat(1.0, 11), np.repeat(np.nan, 11)]),
        }
    )

    # Test conversion function with no filters
    feat_nasc_df = ingest_nasc.convert_afsc_nasc_to_feat(
        mock_df, default_interval_distance=0.25, default_transect_spacing=2.00
    )

    # Check shape and columns
    assert feat_nasc_df.shape == (mock_df.shape[0], mock_df.shape[1] + 1)
    assert list(feat_nasc_df.columns) == [
        "transect_num",
        "distance_s",
        "transect_spacing",
        "distance_e",
    ]

    # Check values
    assert all(feat_nasc_df.loc[np.isnan(mock_df["transect_spacing"]), "transect_spacing"] == 2.00)
    assert all(feat_nasc_df["distance_e"] - feat_nasc_df["distance_s"] == 0.25)

    # Test conversion function with filters
    feat_nasc_filtered_df = ingest_nasc.convert_afsc_nasc_to_feat(
        mock_df,
        default_interval_distance=0.25,
        default_transect_spacing=2.00,
        inclusion_filter={"transect_num": np.arange(1, 8)},
        exclusion_filter={"transect_spacing": 1.0},
    )

    # Check shape
    assert feat_nasc_filtered_df.shape == (3, 4)

    # Check values
    assert all(
        feat_nasc_filtered_df.loc[np.isnan(mock_df["transect_spacing"]), "transect_spacing"] == 2.00
    )
    assert all(feat_nasc_filtered_df["distance_e"] - feat_nasc_filtered_df["distance_s"] == 0.25)
