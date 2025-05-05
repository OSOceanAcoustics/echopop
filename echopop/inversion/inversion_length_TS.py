from typing import Dict, Literal

import pandas as pd

from .inversion_base import InversionBase


class InversionLengthTS(InversionBase):
    """
    Class to perform inversion using length-TS regression.
    """

    def __init__(self, df_model_params: pd.DataFrame):
        super().__init__(df_model_params)

        # Set inversion method
        self.inversion_method = "length_TS_regression"

        # Check df_model_params
        # -- check if df_model_params contain all required parameters
        # -- for length-TS regression these are slope and intercept

    # same as the current aggregate_sigma_bs but with explicit inputs
    def get_stratified_sigma_bs(self, df_length: pd.DataFrame) -> pd.DataFrame:
        # Add a stratum column to df_length if it doesn't exist
        if "stratum" not in df_length.columns:
            df_length["stratum"] = 0

        # Calculate predicted sigma_bs based on length and length_count
        # -- convert TS to the linear domain ('sigma_bs')

        # Calculate mean sigma_bs for each stratum
        # -- impute sigma_bs values if necessary
        df_sigma_bs_stratum: pd.DataFrame

        return df_sigma_bs_stratum

    def invert(self, df_nasc: pd.DataFrame, df_length: pd.DataFrame) -> pd.DataFrame:
        """
        Parameters
        ----------
        df_nasc : pd.DataFrame
            Dataframe with NASC values.
            At the minimum must have the following columns:

            - ("lat", "lon") or ("latitude", "longitude")
            - NASC

            Can optionally has a column "stratum".
            "stratum" must simultanenously exist in `df_nasc` and `df_length`.

        df_length : pd.DataFrame
            Dataframe with scatterer length distribution.
            At the minimum must have the following columns:

            - "length"
            - "length_count"

            Can optionally has a column "stratum".
            "stratum" must simultanenously exist in `df_nasc` and `df_length`.

        Returns
        -------
        df_nasc : pd.DataFrame
            Dataframe with an added number density column
        """

        # Compute mean sigma_bs for each stratum
        df_sigma_bs_stratum = self.get_stratified_sigma_bs(df_length)

        # Add a stratum column to df_nasc if it doesn't exist
        if "stratum" not in df_nasc.columns:
            df_nasc["stratum"] = 0

        # Make sure both df_nasc and df_length have the same stratum values
        if not df_nasc["stratum"].equals(df_length["stratum"]):
            raise ValueError("Stratum values in df_nasc and df_length do not match.")

        # Perform inversion (compute number density based on stratified sigma_bs)

        pass
