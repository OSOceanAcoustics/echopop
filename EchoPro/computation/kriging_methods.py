import numpy as np
import pandas as pd
from ..computation.spatial import lag_distance_griddify , local_search_index
from ..computation.variogram_models import variogram

def compute_kriging_weights( ratio , M2 , K ):

    ### Singular value decomposition (SVD)
    # U: left singular vectors (directions of maximum variance)
    # Sigma: singular values (amount of variance captured by each singular vector, U)
    # VH: conjugate transpose of the right singular vectors
    U , Sigma , VH = np.linalg.svd( K , full_matrices = True )
    
    ### Create Sigma mask informed by the ratio-threshold
    # The ratio between all singular values and their respective
    # maximum is used to apply a mask that informs which values
    # are used to tabulate the kriging weights (aka lambda)
    Sigma_mask = np.abs( Sigma / Sigma[ 0 ] ) > ratio

    ### Inverse masked semivariogram (K)
    K_inv = np.matmul(
        np.matmul( VH.T[:, Sigma_mask] , np.diag( 1.0 / Sigma[ Sigma_mask ] ) ) ,                    
        U[ : , Sigma_mask ].T
    )

    ### Calculate kriging weights (lambda)
    return np.dot( K_inv , M2 )

def compute_kriging_matrix( x_coordinates ,
                            y_coordinates ,
                            variogram_parameters ):
    
    ### Calculate local distance matrix of within-range samples
    local_distance_matrix = lag_distance_griddify( x_coordinates , y_coordinates )
    
    ### Calculate the covariance/kriging matrix (without the constant term)
    kriging_matrix_initial = variogram( distance_lags = local_distance_matrix , 
                                        variogram_parameters = variogram_parameters )

    ### Expand the covariance/kriging matrix with a constant
    # In Ordinary Kriging, this should be '1'
    # Columns
    kriging_matrix = np.concatenate( [ kriging_matrix_initial , 
                                       np.ones( ( x_coordinates.size , 1 ) ) ] , 
                                    axis = 1 )
    
    # Rows
    kriging_matrix = np.concatenate( [ kriging_matrix , 
                                       np.ones( ( 1 , x_coordinates.size + 1 ) ) ] , 
                                    axis = 0 )
    
    ### Diagonal fill (0.0)
    # TODO: Should we put in statements for Objective mapping and Universal Kriging w/ Linear drift?  # noqa
    # Add column and row of ones for Ordinary Kriging
    np.fill_diagonal( kriging_matrix , 0.0 )

    return kriging_matrix

def range_index_threshold( local_point_grid ,
                           distance_matrix ,
                           R ,
                           k_min ):
    
    ### Calculate the within-sample grid indices
    inside_indices = np.where( distance_matrix[ local_point_grid ] <= R )[ 0 ]

    ### Expand search radius if the number of local points are insufficient
    if len( inside_indices ) < k_min:

        ### Sort the closest within-sample points
        inside_indices = np.argsort( distance_matrix[ local_point_grid ] )[ : k_min ]

        ### Look beyond the sample to collect indices beyond the range threshold
        outside_indices = np.where( distance_matrix[ local_point_grid[ inside_indices ] ] > R )[ 0 ]

        ### Calculate the appropraite weights for these values beyond the threshold
        # TODO: should we change this to how Chu does it?
        # tapered function to handle extrapolation
        out_of_sample_weights = np.exp( - np.nanmean( distance_matrix[ local_point_grid[ inside_indices ] ] ) / R )

    else:
        ### If the sample size is sufficient
        outside_indices = [ ]
        out_of_sample_weights = 1.0

    ### Carriage return
    return inside_indices , outside_indices , out_of_sample_weights


def compute_kriging_statistics( point_values ,
                                lagged_semivariogram ,
                                kriging_weights , 
                                inside_indices ,
                                outside_indices ,
                                out_of_sample_weights ):
    
    ### Remove any extrapolation
    if len( outside_indices ) > 0:
        point_values[ outside_indices ] = 0.0
    
    ### Calculate locally weighted kriging mean
    local_mean = np.nansum( kriging_weights[: len( inside_indices ) ] * point_values ) * out_of_sample_weights

    ### Calculate locally weighted kriging prediction variance
    local_prediction_variance = np.nansum( kriging_weights * lagged_semivariogram )

    ### Calculate locally weighted point sample variance
    if abs( local_mean ) < np.finfo( float ).eps:
        local_sample_variance = np.nan
    else:
        local_arithmetic_variance = np.nanvar( point_values , ddof = 1 ) # Non-spatial, arithmetic sample variance
        local_sample_variance = np.sqrt( local_prediction_variance * local_arithmetic_variance ) / abs( local_mean )

    ### Carriage return
    return local_mean , local_prediction_variance , local_sample_variance

def ordinary_kriging( dataframe ,
                      transformed_mesh ,
                      variogram_parameters ,
                      kriging_parameters ,                      
                      variable = 'B_a_adult' ):
    
    ### Calculate the kriging distance matrix and corresponding indices
    distance_matrix , local_point_grid  = local_search_index( transformed_mesh , 
                                                              dataframe , 
                                                              kriging_parameters[ 'kmax' ] )
    
    ### Initialize kriging grid, weights, and results prior to loop    
    kriging_prediction_variance = np.empty( local_point_grid.shape[ 0 ] ) # Prediction variance    
    kriging_sample_variance = np.empty( local_point_grid.shape[ 0 ] ) # Sample variance     
    kriging_mean = np.empty( local_point_grid.shape[ 0 ] ) # Mean

    ### Pull target data variable
    variable_data = dataframe.loc[ : , variable ].values

    ### Iterate through the local point grid to evaluate the kriged/interpolated results
    for row in range( local_point_grid.shape[ 0 ] ):

        ### Parse the k-closest points based on the defined range-threshold
        inside_indices , outside_indices , out_of_sample_weights = (
            range_index_threshold( local_point_grid[ row , : ] ,
                                   distance_matrix[ row , : ] ,
                                   variogram_parameters[ 'range' ] , 
                                   kriging_parameters[ 'kmin' ] )
            )
        
        ### Index the within-sample distance indices
        distance_within_indices = local_point_grid[ row , : ][ inside_indices ]

        ### Calculate the theoretical (semi)variogram at defined lag distances
        modeled_semivariogram = variogram( distance_lags = distance_matrix[ row , distance_within_indices ] , 
                                           variogram_parameters = variogram_parameters )

        ### For Ordinary Kriging, we shift the lags by 1 to include a constant term (1.0)
        # TODO: Should we put in statements for Objective mapping and Universal Kriging w/ Linear drift?  # noqa
        lagged_semivariogram = np.concatenate( [ modeled_semivariogram , [ 1.0 ] ] )

        ### Calculate the covariance/kriging matrix
        kriging_matrix = compute_kriging_matrix( dataframe.x_transformed[ distance_within_indices ].values ,
                                                 dataframe.y_transformed[ distance_within_indices ].values , 
                                                 variogram_parameters = variogram_parameters )

        ### Use singular value decomposition (SVD) to solve for the 
        ### kriging weights (lambda)
        kriging_weights = compute_kriging_weights( kriging_parameters[ 'anisotropy' ] , lagged_semivariogram , kriging_matrix )

        ### Compute Kriging value and variance
        # Index the appropriate data values
        point_values = variable_data[ distance_within_indices ]
        
        ### Calculate the kriged mean, predication variance, and sample variance
        kriged_stats = compute_kriging_statistics( point_values , 
                                                   lagged_semivariogram , 
                                                   kriging_weights , 
                                                   inside_indices , 
                                                   outside_indices ,
                                                   out_of_sample_weights )
        
        ### Fill in the iterable values from the output (tuple) of 
        ### `compute_kriging_statistics()`
        kriging_mean[ row ] = kriged_stats[ 0 ]
        kriging_prediction_variance[ row ] = kriged_stats[ 1 ]
        kriging_sample_variance[ row ] = kriged_stats[ 2 ]

    ### Remove nonsense values (NaN, negative)
    kriging_mean = np.where( ( kriging_mean < 0 ) | np.isnan( kriging_mean ) , 0.0 , kriging_mean )

    ### Carriage return
    return kriging_mean , kriging_prediction_variance , kriging_sample_variance

def kriging_interpolation( dataframe ,
                           transformed_mesh ,
                           dataframe_mesh ,
                           dataframe_geostrata ,
                           variogram_parameters ,
                           kriging_parameters ,                      
                           variable = 'B_a_adult' ):
    
    ### Discretize latitudinal bins
    latitude_bins = np.concatenate( [ [ -90.0 ] , dataframe_geostrata.northlimit_latitude , [ 90.0 ] ] )

    ### Discretize mesh data into the same strata
    dataframe_mesh[ 'stratum_num' ] = pd.cut( dataframe_mesh.centroid_latitude ,
                                              latitude_bins ,
                                              labels = list( dataframe_geostrata.stratum_num ) + [ 1 ] ,
                                              ordered = False )
    
    ### Run kriging
    kriged_results = ordinary_kriging( dataframe , transformed_mesh , variogram_parameters , kriging_parameters )

    ### Append results to the mesh dataframe
    dataframe_mesh[ 'B_a_adult_mean' ] = kriged_results[ 0 ]
    dataframe_mesh[ 'B_a_adult_prediction_variance' ] = kriged_results[ 1 ]
    dataframe_mesh[ 'B_a_adult_sample_variance' ] = kriged_results[ 2 ]

    ### Calculate cell area
    kriged_result_df = (
        dataframe_mesh
        .assign( cell_area_nmi2 = lambda x: kriging_parameters[ 'A0' ] * x.fraction_cell_in_polygon ,
                 B_adult_kriged = lambda x: x.B_a_adult_mean * x.cell_area_nmi2 )
    )

    ### Carriage return
    return kriged_result_df