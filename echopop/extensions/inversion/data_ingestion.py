from pydantic import BaseModel
from pandera import DataFrameModel
from lmfit import Parameters
from typing import Any, Dict, List, Union

import pandas as pd

####################################################################################################
# Validation / preparation
####################################################################################################

def inversion_configuration_validator() -> Dict[str, Any]:
    """
    Pydantic model for validating configuration parameters
    """
    pass

def dataset_validator() -> pd.DataFrame:
    """
    Pandera model for validating dataset values
    """
    pass

def prepare_scattering_model_inputs() -> Dict[str, Any]: 
    """
    Prepare scattering model parameter inputs
    """
    # == functions/set_para.m
    # == functions/inversion_para.m
    
    # PHASE 1) INGEST VALUES FROM CONFIGURATION AND DATASET FILE
    # PHASE 2) COMPUTE INTERMEDIATE VARIABLES (e.g. acoustic wavenumber, position matrix)
    # PHASE 3) PASS TO SCATTERER CLASS
    #          --> EXTERNAL TO THIS FUNCTION
    pass    

####################################################################################################
# Data ingestion
####################################################################################################

def yaml_configuration_reader() -> Dict[str, Union[float, int, Parameters, pd.DataFrame, str]]:
    """
    Read and validate the input parameterization YAML configuration
    """
    # == functions/load_para_data.m
    # == functions/load_geo_phy_para.m
    # == functions/get_simu_para.m

    # PHASE 1) READ CONFIGURATION FILE
    # PHASE 2) VALIDATE CONFIGURATION SETTINGS
    # PHASE 3) PASS CONFIGURATION SETTINGS TO SCATTERER CLASS AS APPROPRIATE `lmfit::PARAMETERS`
    #          --> EXTERNAL TO THIS FUNCTION
    pass

def dataset_reader() -> pd.DataFrame:
    """
    Read aggregate acoustic backscatter measurements
    """
    # == functions/get_acoustic_data.m
    # == functions/load_MOCNESS_data.m
    # == functions/load_BIOMAPPER_data.m
    
    # PHASE 1) READ IN FILES
    # PHASE 2) VALIDATE DATASET FILE
    # PHASE 3) PASS DATASET TO SCATTERER CLASS
    #          --> EXTERNAL TO THIS FUNCTION
    pass