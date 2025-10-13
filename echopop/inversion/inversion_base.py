import abc

import pandas as pd


class InversionBase(abc.ABC):
    """
    Abstract base class for handling acoustic inversion methods.
    
    This class provides a framework for different types of acoustic inversions
    by establishing common interfaces and shared functionality for parameter
    management and stratification handling.
    
    Parameters
    ----------
    model_parameters : Dict[str, Any]
        Dictionary containing model configuration parameters. Common keys include:
        - 'stratify_by': str or List[str] - columns to stratify by
        - 'strata': array-like - specific strata to process
        - 'impute_missing_strata': bool - whether to impute missing strata
        
    Attributes
    ----------
    model_params : Dict[str, Any]
        Processed model parameters
    inversion_method : str
        String identifier for the specific inversion method (set by subclasses)
        
    Examples
    --------
    >>> # Example parameters for length-TS inversion
    >>> params = {
    ...     "ts_length_regression": {"slope": 20.0, "intercept": -68.0},
    ...     "stratify_by": "stratum_ks",
    ...     "strata": [1, 2, 3, 4, 5],
    ...     "impute_missing_strata": True
    ... }
    >>> 
    >>> # Create concrete inversion class (subclass)
    >>> inverter = InversionLengthTS(params)
    >>> print(inverter.inversion_method)
    length_TS_regression
        
    Notes
    -----
    This is an abstract base class and cannot be instantiated directly.
    Subclasses must implement the abstract `invert` method.
    
    The class automatically converts single-string 'stratify_by' parameters
    to lists for consistent handling across different inversion methods.
    """

    def __init__(self, model_parameters):

        # Ingest model parameters
        self.model_params = model_parameters

        # Modify "stratify_by" if needed
        if "stratify_by" in self.model_params:
            if isinstance(self.model_params["stratify_by"], str):
                self.model_params["stratify_by"] = [self.model_params["stratify_by"]]

        # Initialize method
        self.inversion_method = ""    

    @abc.abstractmethod        
    def invert(self, df_nasc: pd.DataFrame) -> pd.DataFrame:
        """
        Perform inversion on input dataframe to convert NASC to number density.

        This abstract method must be implemented by subclasses to define the
        specific inversion algorithm (e.g., length-TS, age-TS, species-specific).

        Parameters
        ----------
        df_nasc : pd.DataFrame
            Input dataframe with NASC (Nautical Area Scattering Coefficient) values
            to perform inversion on. Must contain 'nasc' column and any required
            stratification columns.

        Returns
        -------
        pd.DataFrame
            Input dataframe with added 'number_density' column containing the
            inverted number density estimates.
            
        Notes
        -----
        To perform inversion on non-stratum groupings, pre-process the dataframe
        so that each row contains the minimum unit inversion will be performed on.
        
        Subclasses should implement the specific inversion algorithm appropriate
        for their method (length-based, age-based, etc.).
        """
        pass