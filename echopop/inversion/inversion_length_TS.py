import pandas as pd
import numpy as np
from typing import Any, Dict, List, Optional, Union
from pydantic import field_validator, Field

from .inversion_base import InversionBase
from .operations import impute_missing_sigma_bs
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
        - 'expected_strata': array-like - specific strata to process
        - 'impute_missing_strata': bool - whether to impute missing strata values
        - 'haul_replicates': bool - whether to use haul numbers/ids as replicates instead of 
        individuals
        
    Attributes
    ----------
    sigma_bs_haul : pd.DataFrame or None
        Haul-level backscattering cross-sections
    sigma_bs_strata : pd.DataFrame or None
        Stratum-level backscattering cross-sections
        
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
        # Validate and parse the model parameters
        validated_params = ValidateLengthTS.create(**model_parameters)
        super().__init__(validated_params)
        
        # Set inversion method
        self.inversion_method = "length_TS_regression"
        
        # Initialize attributes
        self.sigma_bs_haul = None
        self.sigma_bs_strata = None
        
    def set_haul_sigma_bs(
        self, 
        df_length: Union[pd.DataFrame, List[pd.DataFrame]]
    ) -> pd.DataFrame:
        """
        Compute the mean linear scattering coefficient (sigma_bs) for each haul.
        
        This method processes length data from biological samples to calculate haul-specific
        backscattering cross-sections using the TS-length regression relationship. 
        
        Parameters
        ----------
        df_length : pd.DataFrame or list of pd.DataFrame
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
        
        Examples
        --------
        >>> # Single DataFrame
        >>> inverter.set_haul_sigma_bs(specimen_df)
        >>> 
        >>> # Multiple DataFrames
        >>> inverter.set_haul_sigma_bs([specimen_df, length_df])
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
            .apply(lambda x: np.average(x.sigma_bs, weights=x.length_count))        
        )
        
        # Return the output
        return sigma_bs_haul.to_frame("sigma_bs")
        
    def get_stratified_sigma_bs(
        self, 
        df_length: Union[pd.DataFrame, List[pd.DataFrame]],
        haul_replicates: bool = True,
    ) -> None:
        """
        Calculate stratified mean backscattering cross-sections (sigma_bs) by stratum.
        
        This method computes stratum-level sigma_bs values from either all individuals directly, or 
        by hauls to account for pseudoreplication. The latter would treat each haul as an 
        individual replicate (i.e. the unit of replication is the haul)[1]_. The resulting values 
        are used in acoustic inversion to convert NASC to number density (animals nmi^-2).
        
        Parameters
        ----------
        df_length : pd.DataFrame or list of pd.DataFrame
            Length data for computing sigma_bs.
            
            Required columns (when provided):
            - All columns specified in model_params["stratify_by"]
            - "length": Fish length measurements
            - "length_count" (optional): Pre-aggregated counts per length
            - "haul_num": When `haul_replicates=True`, otherwise it is optional.

        haul_replicates : bool
            When True, individual hauls are used as replicates instead of individual lengths. This 
            involves averaging sigma_bs for each haul.
            
        Returns
        -------
        pd.DataFrame
            DataFrame with stratified sigma_bs values. Index contains stratum identifiers,
            and 'sigma_bs' column contains the mean backscattering cross-section values
            for each stratum.
            
        Notes
        -----
        Using hauls as the unit of replication is recommended to avoid pseudoreplication when
        individuals within hauls are not independent[1]_.
           
        The sigma_bs calculation follows:
        1. Convert length to TS using regression: TS = slope * log10(length) + intercept
        2. Convert TS to linear scale: sigma_bs = 10^(TS/10)
        3. Calculate weighted average by length frequency within each stratum
        
        Examples
        --------
        >>> # Direct calculation from new data
        >>> sigma_bs = inverter.get_stratified_sigma_bs(specimen_df)
        
        References
        ----------
        .. [1] Hurlbert, S.H. (1984). Pseudoreplication and the Design of Ecological Field 
        Experiments. *Ecological Monographs*, 54(2), 187-211. https://doi.org/10.2307/1942661
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

        # Use hauls as replicates
        if haul_replicates:
            # ---- Store the values
            self.sigma_bs_haul = self.set_haul_sigma_bs(df_length)
            # ---- Set the mean sigma_bs for each stratum
            self.sigma_bs_strata = self.sigma_bs_haul["sigma_bs"].unstack(
                self.model_params["stratify_by"]
            ).mean(axis=0).to_frame("sigma_bs")
        # Use individuals as replicates
        else:
            # ---- Prepare the calculation depending on if `df` is a single DataFrame
            if isinstance(df_length, pd.DataFrame):
                df_length = [df_length]                
            # ---- Quantize the length counts across all datasets per length value
            df_length_counts = (
                pd.concat([utils.quantize_length_data(d, self.model_params["stratify_by"]) 
                        for d in df_length], axis=1)
                .fillna(0.).sum(axis=1).reset_index(name="length_count")
            )            
            # ---- Compute the average TS
            df_length_counts["TS"] = acoustics.ts_length_regression(
                df_length_counts["length"], 
                self.model_params["ts_length_regression"]["slope"],
                self.model_params["ts_length_regression"]["intercept"]
            )            
            # ---- Linearize
            df_length_counts["sigma_bs"] = 10. ** (df_length_counts["TS"] / 10.)
            # ---- Aggregate by stratum
            self.sigma_bs_strata = (
                df_length_counts
                .groupby(self.model_params["stratify_by"])[["length_count", "sigma_bs"]].apply(
                    lambda x: np.average(x.sigma_bs, weights=x.length_count)
                )
            ).to_frame("sigma_bs")
        
    def invert(self, 
               df_nasc: pd.DataFrame, 
               df_length: Union[pd.DataFrame, List[pd.DataFrame]]):
        
        # Create copy
        df_nasc = df_nasc.copy()
        
        # When stratified               
        if "stratify_by" in self.model_params:
            
            # Get the unique strata
            if self.model_params["expected_strata"] is not None:
                unique_strata = np.unique(self.model_params["expected_strata"])
            else:
                unique_strata = np.unique(df_nasc[self.model_params["stratify_by"]])
                
            # Get the mean linear scattering coefficient for each stratum
            self.get_stratified_sigma_bs(df_length,
                                         self.model_params["haul_replicates"]) 
            
            # Impute, if defined
            if (
                "impute_missing_strata" in self.model_params and 
                self.model_params["impute_missing_strata"]
            ):
                self.sigma_bs_strata = impute_missing_sigma_bs(unique_strata,
                                                               self.sigma_bs_strata)
            
            # Index `df_nasc` to align with `sigma_bs_strata` 
            df_nasc.set_index(self.model_params["stratify_by"], inplace=True)
            
            # Compute number density
            df_nasc["number_density"] = (
                df_nasc["nasc_proportion"].fillna(0.) * df_nasc["nasc"] / 
                (4 * np.pi * self.sigma_bs_strata["sigma_bs"].reindex_like(df_nasc))
            ).astype(float).fillna(0.)
            
            # Return the NASC DataFrame with the inverted number densities
            return df_nasc.reset_index()

####################################################################################################
# Validation
# ----------

class ValidateLengthTS(
    utils.InputModel, 
    arbitrary_types_allowed=True, 
    title="TS-length inversion model parameters"
):
    """
    Validation model for TS-length inversion parameters used by InversionLengthTS.

    This Pydantic model validates and documents the configuration parameters required by the 
    InversionLengthTS class for acoustic inversion using length-TS regression.

    Parameters
    ----------
    ts_length_regression : utils.TSLRegressionParameters
        Regression parameters for converting fish length to target strength (TS).
    stratify_by : List[str]
        List of column names used for data stratification (e.g., 'stratum_ks').
    expected_strata : np.ndarray, optional
        Array of expected strata identifiers to process.
    impute_missing_strata : bool, default=True
        Whether to impute missing strata values during inversion.
    haul_replicates : bool, default=True
        Whether to use hauls as the statistical unit of replication (recommended to avoid
        pseudoreplication[1]_).

    Notes
    -----
    This model is intended for use with the InversionLengthTS class, which performs
    acoustic inversion by relating fish length to acoustic backscatter using empirical
    TS-length relationships.

    Using hauls as the unit of replication is recommended to avoid pseudoreplication when
    individuals within hauls are not independent[1]_.

    References
    ----------
    .. [1] Hurlbert, S.H. (1984). Pseudoreplication and the Design of Ecological Field Experiments.
       *Ecological Monographs*, 54(2), 187-211. https://doi.org/10.2307/1942661
    """

    ts_length_regression: utils.TSLRegressionParameters
    stratify_by: List[str]
    expected_strata: Optional[np.ndarray[np.number]] = Field(default=None)
    impute_missing_strata: bool = Field(default=True)
    haul_replicates: bool = Field(default=True)

    @field_validator("stratify_by", mode="before")
    def validate_stratify_by(cls, v):
        if isinstance(v, str):
            v = [v]
        return v

    @field_validator("expected_strata", mode="before")
    def validate_expected_strata(cls, v):
        if v is None:
            return v

        if isinstance(v, list):
            return np.array(v)

        return v
