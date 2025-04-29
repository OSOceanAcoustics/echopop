import pandas as pd
from .inversion_base import InversionBase


class InversionMatrix(InversionBase):
    """
    Class to perform matrix inversion for krill.

    Reference:
    Chu, D., Lawson, G. L., and Wiebe, P. H. (2016).
    “Estimation of biological parameters of marine organisms using linear and nonlinear 
    acoustic scattering model-based inversion methods,” 
    The Journal of the Acoustical Society of America, 139, 2885–2895. doi:10.1121/1.4948759
    """

    def __init__(self, df_model_params: pd.DataFrame):
        super().__init__(df_model_params)
        
        # Set inversion method
        self.inversion_method = "krill_matrix_inversion"

        # Check df_model_params
        # -- check if df_model_params contain all required parameters
        # -- for matrix inversion these are the krill shape parameters
        
    def invert(self, df_nasc: pd.DataFrame) -> pd.DataFrame:

        # Krill inversion ops
        # -- note any non-stratum grouping should be performed OUTSIDE of this class

        pass
