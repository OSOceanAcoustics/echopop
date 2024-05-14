import numpy as np
import pandas as pd

from ..spatial.mesh import griddify_lag_distances
from ..spatial.transect import define_western_extent

def kriging( transect_data: pd.DataFrame ,
             mesh_data: pd.DataFrame ,
             settings_dict: dict ):
    """
    Use kriging to interoplate data
    
    Parameters
    ----------
    transect_data: pd.DataFrame
        Dataframe including georeferenced data
    mesh_data: pd.DataFrame
        Grid data that has been transformed
    settings_dict: dict
        Kriging and variogram model parameters
    """

    # Generate the distance matrix for each mesh point relative to all transect coordinates
    distance_matrix = griddify_lag_distances( mesh_data , transect_data )

    # Calculate the western extent of the transect data
    western_extent = define_western_extent( transect_data )



def search_radius_mask( distance_matrix: np.ndarray , search_radius: float ) :

    # Create copy of matrix
    matrix_copy = distance_matrix.copy( )
    
    # Find values beyond the maximum search radius and mask (assign NaN) to those beyond
    # ---- Generate mask
    mask = matrix_copy > search_radius 
    # ---- Set values outside the radius to NaN
    matrix_copy[ mask ] = np.nan

    # Return output
    return matrix_copy

def count_within_radius( distance_matrirx_masked: np.ndarray ):

    # Create copy of matrix
    matrix_copy = distance_matrirx_masked.copy( )

    # Create boolean matrix of values outside search radius
    nan_mask = np.isnan( matrix_copy )

    # Count the number of values within the search radius and return output
    return np.sum( ~ nan_mask , axis = 1 )

def adaptive_search_radius( distance_matrix: np.ndarray ,
                            mesh_data: pd.DataFrame ,
                            western_extent: pd.DataFrame ,
                            kriging_parameters: dict ) :
    """
    Find the indices of the k-th nearest points (relative to a reference coordinate) required
    for computing the lagged semivariogram   

    Parameters
    ----------
    distance_matrix: np.ndarray
        An array/matrix that includes the distances of each mesh points from every 
        georeferenced along-transect interval
    tranect_extent: pd.DataFrame
        A DataFrame that includes the x- and y-coordinates of the western-most coordinates for 
        each survey transect
    k_min: int 
        The minimum number of nearest neighbors. Mesh points with fewer than `k_min` valid distances
        are supplemented with extrapolated ranges numbering up to `k_min`.
    k_max: int
        The maximum number of nearest neighbors.
    R: float 
        The search radius used to identify between `k_min` and `k_max` nearest neighbors.
    """     

    # Extract key search radius parameters
    # ---- k_min
    k_min = kriging_parameters[ 'kmin' ]
    # ---- k_max
    k_max = kriging_parameters[ 'kmax' ]
    # ---- Search radius (distance)
    search_radius = kriging_parameters[ 'search_radius' ]

    # Generate the search radius mask
    distance_matrix_masked = (
        search_radius_mask( distance_matrix ,
                            search_radius )
    )
    
    # Identify mesh points that require extrapolation
    # ---- Count the number of values within the search radius
    valid_distances = count_within_radius( distance_matrix_masked )
    # ---- Identify rows where the number of valid points are less than `k_min`
    sparse_radii = np.hstack( np.where( valid_distances < k_min ) )

    # Calculate the closest grid points and their indices
    # ---- Closest indices
    local_indices = distance_matrix.argsort( axis = 1 )
    # ---- Map the distance matrix to the local indices
    local_points = np.take_along_axis( distance_matrix , local_indices , axis = 1 )

    # Initialize matrices
    # ---- Within-radius (WR) samples
    wr_indices = local_indices[ : , : k_max ].astype( float )                             
    # ---- Out-of-sample (OOS) indices
    oos_indices = np.full( ( len( valid_distances ) , k_min ) , np.nan )
    # ---- OOS weights
    oos_weights = np.ones( len( valid_distances ) )

    # For points where there are fewer than `k_min` points within the search radius, extrapolate 
    if np.size( sparse_radii ) > 0 :
        # Extract the indices requiring extrapolation
        nearby_indices = local_indices[ sparse_radii ][ : , : k_min ]
        # Index the mesh grid coordinates for bounding the search radius expansion/extrapolation
        # ---- y-coordinates with array transformation to access matrix operations
        mesh_y = mesh_data[ 'y' ][ sparse_radii ].to_numpy( ).reshape(-1, 1)
        # ---- x-coordinates
        mesh_x = mesh_data[ 'x' ][ sparse_radii ].to_numpy( )

        # Update local points
        # ---- Fill NaN values
        local_points[ sparse_radii , k_min : ] = np.nan 

        # Calculate the mesh distance from the western boundary of the survey transects
        # ---- Find closest point
        mesh_western_distance = np.abs( mesh_y - western_extent['y'].to_numpy( ) ).argmin( axis = 1)
        # ---- Calculate the western limits (x-axis)
        western_limit = western_extent.iloc[ np.ravel( mesh_western_distance ) ][ 'x' ]
        # ---- Compute bounding threshold (for tapered extrapolation function)
        western_threshold = western_limit - search_radius

        # Taper function for extrapolating values outside the search radius
        if np.any( mesh_x < western_threshold ) :
            # ---- Create a thresholded mask for lazy operations
            western_limit_mask = mesh_x < western_threshold
            # ---- Index these values
            extrapolation_index = sparse_radii[ western_limit_mask ]
            # ---- Compute the OOS kriging weights
            oos_mean = np.apply_along_axis( np.nanmean , 
                                            1 , 
                                            local_points[ extrapolation_index , : k_min ] )
            # ---- Exponentiate the OOS mean
            oos_exp = np.exp( - oos_mean / search_radius )
            # ---- Update the OOS weights
            oos_weights[ extrapolation_index ] = oos_exp
            # ---- Get the outside indices that correspond to this tapered extrapolation
            sparse_extrapolation_index = nearby_indices[ western_limit_mask ].astype( float )
            # ---- Apply indices as a mask to the NaN-masked distance matrix
            extrapolated_distance = (
                np.take_along_axis( distance_matrix_masked[ sparse_radii ][ western_limit_mask ] ,
                                    sparse_extrapolation_index.astype( int ) , axis = 1 )
            )
            # ---- Create NaN mask
            extrapolated_nan_mask = ~ np.isnan( extrapolated_distance )
            # -------- Apply mask to indices
            sparse_extrapolation_index_nan = sparse_extrapolation_index.copy( )
            sparse_extrapolation_index_nan[ extrapolated_nan_mask ] = np.nan
            # -------- Update `out_of_sample_indices` matrix
            oos_indices[ extrapolation_index ] = np.sort( sparse_extrapolation_index_nan )
            # ---- Get inside indices that apply to these points
            # -------- Create NaN mask for within-sample values
            interpolated_nan_mask = np.isnan( extrapolated_distance )
            # -------- Apply mask to indices
            sparse_interpolation_index_nan = sparse_extrapolation_index.copy( )
            sparse_interpolation_index_nan[ interpolated_nan_mask ] = np.nan
            # -------- Pad NaN to match `within_sample_indices` matrix
            sparse_interpolation_pad = np.pad( sparse_interpolation_index_nan ,
                                               [ ( 0 , 0 ) , ( 0 , k_max - k_min ) ] ,
                                               mode = 'constant' ,
                                               constant_values = np.nan )
            # -------- Updated `within_sample_indices` matrix
            wr_indices[ extrapolation_index ] = np.sort( sparse_interpolation_pad )

            print( f"""Extrapolation applied to kriging mesh points ({len( sparse_radii )} of {wr_indices.shape[ 0 ]}):
                    * {len( valid_distances[ valid_distances == 0 ] )} points had 0 valid range estimates without extrapolation
                    * {len( valid_distances[ ( valid_distances != 0 ) & ( valid_distances< k_min ) ] )} points had at least 1 valid point but fewer than {k_min} valid neighbors""")

    # Return output
    return local_points[ : , : k_max ] , wr_indices , oos_indices , oos_weights

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
