import pandas as pd
import numpy as np
from typing import Any, Dict, List, Optional, Union

from .inversion_base import InversionBase
from ..nwfsc_feat import utils
from .. import acoustics

class InversionLengthTS(InversionBase):
    """
    Class to perform acoustic inversion using length-TS (Target Strength) regression.
    
    This class implements acoustic inversion by relating fish length to acoustic
    backscatter through empirical TS-length relationships. It calculates stratified
    mean backscattering cross-sections and converts NASC to number density.
    
    Parameters
    ----------
    model_parameters : Dict[str, Any]
        Dictionary containing model configuration. Required keys:
        - 'ts_length_regression': Dict with 'slope' and 'intercept' for TS equation
        - 'stratify_by': str or List[str] - stratification columns (e.g., 'stratum_ks')
        
        Optional keys:
        - 'strata': array-like - specific strata to process
        - 'impute_missing_strata': bool - whether to impute missing strata values
        
    Attributes
    ----------
    sigma_bs_haul : pd.Series or None
        Cached haul-level backscattering cross-sections after set_haul_sigma_bs()
        
    Examples
    --------
    >>> # Define model parameters
    >>> params = {
    ...     "ts_length_regression": {"slope": 20.0, "intercept": -68.0},
    ...     "stratify_by": "stratum_ks",
    ...     "strata": [1, 2, 3, 4, 5],
    ...     "impute_missing_strata": True
    ... }
    >>> 
    >>> # Initialize inverter
    >>> inverter = InversionLengthTS(params)
    >>> 
    >>> # Set haul-level sigma_bs from biological data
    >>> inverter.set_haul_sigma_bs([specimen_df, length_df])
    >>> 
    >>> # Perform inversion on NASC data
    >>> result = inverter.invert(nasc_df)
    >>> print(result['number_density'].sum())
    1234567.89
        
    Notes
    -----
    The inversion process follows these steps:
    1. Quantize length data and compute TS using the regression equation
    2. Convert TS to linear backscattering cross-section (sigma_bs)
    3. Calculate stratified mean sigma_bs values
    4. Impute missing strata if specified
    5. Convert NASC to number density: number_density = NASC / (4π * sigma_bs)
    
    The TS-length relationship follows: TS = slope * log10(length) + intercept
    """
    
    def __init__(self, model_parameters: Dict[str, Any]):
        super().__init__(model_parameters)
        
        # Set inversion method
        self.inversion_method = "length_TS_regression"
        
        # Initialize internal values     
        self.sigma_bs_haul = None
        self.sigma_bs_strata = None
        
    def set_haul_sigma_bs(self, 
                          df_length: Optional[Union[pd.DataFrame, 
                                                    List[pd.DataFrame], 
                                                    Dict[str, pd.DataFrame]]] = None):
        """
        Compute and cache the mean linear scattering coefficient (sigma_bs) for each haul.
        
        This method processes length data from biological samples to calculate haul-specific
        backscattering cross-sections using the TS-length regression relationship. The results
        are cached internally for later use in stratified sigma_bs calculations.
        
        Parameters
        ----------
        df_length : pd.DataFrame, list of pd.DataFrame, or dict of pd.DataFrame
            Length data containing fish measurements. Can be:
            - Single DataFrame with length measurements per haul
            - List of DataFrames to be concatenated
            - Dictionary of DataFrames (values will be used)
            
            Required columns:
            - All columns specified in model_params["stratify_by"] 
            - "haul_num": Haul identifier
            - "length": Fish length measurements
            - "length_count" (optional): Pre-aggregated counts per length
            
        Notes
        -----
        The method performs the following steps:
        1. Quantizes length data by stratification variables and haul
        2. Applies TS-length regression to convert lengths to target strength
        3. Converts TS to linear backscattering cross-section (sigma_bs)
        4. Calculates length-weighted average sigma_bs per haul
        5. Caches results in self.sigma_bs_haul for later use
        
        The cached results can be accessed later via get_stratified_sigma_bs()
        without providing df_length parameter.
        
        Examples
        --------
        >>> # Single DataFrame
        >>> inverter.set_haul_sigma_bs(specimen_df)
        >>> 
        >>> # Multiple DataFrames
        >>> inverter.set_haul_sigma_bs([specimen_df, length_df])
        >>> 
        >>> # Dictionary of DataFrames
        >>> data_dict = {"bio": specimen_df, "length": length_df}
        >>> inverter.set_haul_sigma_bs(data_dict)
        """
        
        # Prepare the calculation depending on if `df` is a single DataFrame or Dictionary
        if isinstance(df_length, pd.DataFrame):
            df_length = [df_length]
        elif isinstance(df_length, dict):
            df_length = [d for d in df_length.values()]
            
        # Quantize the lengths
        haul_counts = pd.concat(
            [utils.quantize_length_data(d, self.model_params["stratify_by"] + ["haul_num"]) 
             for d in df_length],
            axis=1
        )
        
        # Fill any mismatched indices
        haul_counts = haul_counts.fillna(0.).sum(axis=1).reset_index(name="length_count")
        
        # Compute the average TS
        haul_counts["TS"] = acoustics.ts_length_regression(
            haul_counts["length"], 
            self.model_params["ts_length_regression"]["slope"],
            self.model_params["ts_length_regression"]["intercept"]
        )
        
        # Linearize
        haul_counts["sigma_bs"] = 10. ** (haul_counts["TS"] / 10.)
        
        # Weighted average across hauls
        sigma_bs_haul = (
            haul_counts
            .groupby(self.model_params["stratify_by"] + ["haul_num"])[["sigma_bs", "length_count"]]
            .apply(lambda x: np.average(x.sigma_bs, weights=x.length_count))        )
        
        # Store
        self.sigma_bs_haul = sigma_bs_haul
        
    def get_stratified_sigma_bs(
        self, 
        df_length: Optional[Union[pd.DataFrame, 
                                  List[pd.DataFrame], 
                                  Dict[str, pd.DataFrame]]] = None
    ) -> pd.DataFrame:
        """
        Calculate stratified mean backscattering cross-sections (sigma_bs) by stratum.
        
        This method computes stratum-level sigma_bs values either from cached haul data
        (if set_haul_sigma_bs was called previously) or by processing new length data.
        The resulting values are used in acoustic inversion to convert NASC to number density.
        
        Parameters
        ----------
        df_length : pd.DataFrame, list of pd.DataFrame, dict of pd.DataFrame, or None
            Length data for computing sigma_bs. If None, uses cached haul-level data
            from previous set_haul_sigma_bs() call. If provided, processes the data
            to compute fresh sigma_bs values.
            
            Required columns (when provided):
            - All columns specified in model_params["stratify_by"]
            - "length": Fish length measurements
            - "length_count" (optional): Pre-aggregated counts per length
            
        Returns
        -------
        pd.DataFrame
            DataFrame with stratified sigma_bs values. Index contains stratum identifiers,
            and 'sigma_bs' column contains the mean backscattering cross-section values
            for each stratum.
            
        Notes
        -----
        This method has two modes of operation:
        
        1. **Cached mode** (df_length=None):
           - Uses haul-level sigma_bs from set_haul_sigma_bs()
           - Averages across hauls within each stratum
           - Faster for repeated calculations
           
        2. **Direct calculation mode** (df_length provided):
           - Processes length data directly
           - Quantizes by stratification variables only (no haul grouping)
           - Computes length-weighted average sigma_bs per stratum
           
        The sigma_bs calculation follows:
        1. Convert length to TS using regression: TS = slope * log10(length) + intercept
        2. Convert TS to linear scale: sigma_bs = 10^(TS/10)
        3. Calculate weighted average by length frequency within each stratum
        
        Examples
        --------
        >>> # Using cached haul data
        >>> inverter.set_haul_sigma_bs(specimen_df)
        >>> sigma_bs = inverter.get_stratified_sigma_bs()
        >>> print(sigma_bs)
                   sigma_bs
        stratum_ks         
        1           0.00012
        2           0.00015
        3           0.00011
        
        >>> # Direct calculation from new data
        >>> sigma_bs = inverter.get_stratified_sigma_bs(specimen_df)
        
        Raises
        ------
        AttributeError
            If df_length is None but no cached haul data exists (set_haul_sigma_bs 
            was not called previously).
        """
        
        # Grab mean hauls if `df_length` is not specified
        if df_length is None and self.sigma_bs_haul is not None:
            sigma_bs_strata = self.sigma_bs_haul.unstack(
                self.model_params["stratify_by"]
            ).mean(axis=0)
        # ---- Otherwise, apply the appropriate groupby operation
        else:
            # ---- Prepare the calculation depending on if `df` is a single DataFrame or Dictionary
            if isinstance(df_length, pd.DataFrame):
                df_length = [df_length]
            elif isinstance(df_length, dict):
                df_length = [d for d in df_length.values()]
                
            # Quantize the length counts across all datasets per length value
            df_length_counts = (
                pd.concat([utils.quantize_length_data(d, self.model_params["stratify_by"]) 
                        for d in df_length], axis=1)
                .fillna(0.).sum(axis=1).reset_index(name="length_count")
            )
            
            # Compute the average TS
            df_length_counts["TS"] = acoustics.ts_length_regression(
                df_length_counts["length"], 
                self.model_params["ts_length_regression"]["slope"],
                self.model_params["ts_length_regression"]["intercept"]
            )
            
            # Linearize
            df_length_counts["sigma_bs"] = 10. ** (df_length_counts["TS"] / 10.)
            
            # Aggregate by stratum
            sigma_bs_strata = (
                df_length_counts
                .groupby(["stratum_ks"])[["length_count", "sigma_bs"]].apply(
                    lambda x: np.average(x.sigma_bs, weights=x.length_count)
                )
            )

        # Store
        self.sigma_bs_strata = sigma_bs_strata.to_frame("sigma_bs")
            
        # Return the average sigma_bs per stratum
        return sigma_bs_strata.to_frame("sigma_bs")
        
    def invert(self, 
               df_nasc: pd.DataFrame, 
               df_length: Optional[Union[pd.DataFrame, 
                                  List[pd.DataFrame], 
                                  Dict[str, pd.DataFrame]]] = None):
        
        # Create copy
        df_nasc = df_nasc.copy()
        
        # When stratified               
        if "stratify_by" in self.model_params:
            
            # Get the unique strata
            if "strata" in self.model_params:
                unique_strata = np.unique(self.model_params["strata"])
            else:
                unique_strata = np.unique(df_nasc[self.model_params["stratify_by"]])
                
            # Get the mean linear scattering coefficnet for each stratum
            sigma_bs_strata = self.get_stratified_sigma_bs(df_length) 
            
            # Impute, if defined
            if (
                "impute_missing_strata" in self.model_params and 
                self.model_params["impute_missing_strata"]
            ):
                sigma_bs_strata = utils.impute_missing_sigma_bs(unique_strata, sigma_bs_strata)
            
            # Index `df_nasc` to align with `sigma_bs_strata` 
            df_nasc.set_index(self.model_params["stratify_by"], inplace=True)
            
            # Compute number density
            df_nasc["number_density"] = (
                df_nasc["nasc_proportion"].fillna(0.) * df_nasc["nasc"] / 
                (4 * np.pi * sigma_bs_strata["sigma_bs"].reindex_like(df_nasc))
            ).astype(float).fillna(0.)
            
            # Return the NASC DataFrame with the inverted number densities
            return df_nasc.reset_index()