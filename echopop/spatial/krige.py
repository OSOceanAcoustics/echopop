import numpy as np
import pandas as pd

def kriging_lambda( anisotropy: float , 
                    lagged_semivariogram: np.ndarray , 
                    kriging_matrix: np.ndarray , ):
    """
    Apply singular value decomposition (SVD) to compute kriging (lambda) weights    

    Parameters
    ----------
    anisotropy: np.float64
        Anisotropy ratio.
    lagged_semivariogram: np.array
        Lagged semivariogram
    kriging_matrix: np.array
        Kriging matrix.
    """
    # Singular value decomposition (SVD)
    # ---- U: left singular vectors (directions of maximum variance)
    # ---- Sigma: singular values (amount of variance captured by each singular vector, U)
    # ---- VH: conjugate transpose of the right singular vectors
    U , Sigma , VH = np.linalg.svd( kriging_matrix , full_matrices = True )
    
    # Create Sigma mask informed by the ratio-threshold
    # ---- The ratio between all singular values and their respective
    # ---- maximum is used to apply a mask that informs which values
    # ---- are used to tabulate the kriging weights (aka lambda)
    Sigma_mask = np.abs( Sigma / Sigma[ 0 ] ) > anisotropy

    # Inverse masked semivariogram (K)
    K_inv = np.matmul(
        np.matmul( VH.T[:, Sigma_mask] , np.diag( 1.0 / Sigma[ Sigma_mask ] ) ) ,                    
        U[ : , Sigma_mask ].T
    )

    # Calculate kriging weights (lambda)
    return np.dot( K_inv , lagged_semivariogram )
