import abc
import pandas as pd


class InversionBase(abc.ABC):
    """
    Class to handle inversion.

    Allow stratum-specific model parameters.
    """

    def __init__(self, df_model_params: pd.DataFrame):
        
        # Ingest model parameters
        # -- this is a dataframe for convenience to manage strata
        self.df_model_params = df_model_params  

        # Add a "stratum" column if not already exists
        if "stratum" not in df_model_params:
            df_model_params["stratum"] = 0  # enable groupby/index-based ops in invert()

        self.inversion_method = ""  # inversion method set in child class

    @abc.abstractmethod
    def invert(self, df_nasc: pd.DataFrame) ->  pd.DataFrame:
        """
        Perform inversion on each row of the input dataframe.

        To perform inversion on non-stratum groupings, pre-process the dataframe
        so that each row contains the minimum unit inversion will be performed on.

        Parameters
        ----------
        df_nasc : pd.DataFrame
            Input dataframe with NASC to perform inversion on

        Returns
        -------
        Input dataframe with the inverted number density as an added column
        """
        pass