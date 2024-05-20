import numpy as np
import pandas as pd
import scipy.stats as st

from .spatial.transect import transect_array

def stratified_transect_statistic( transect_data: pd.DataFrame , 
                                   transect_summary: pd.DataFrame , 
                                   strata_summary: pd.DataFrame , 
                                   settings_dict: dict ) -> tuple[ pd.DataFrame , dict ]:    
    """
    Calculates stratified mean statistics for a set of transects
    

    Parameters
    ----------
    transect_data: pd.DataFrame
        Dataframe comprising georeferenced biological data collected from survey transects
    transect_summary: pd.DataFrame
        DataFrame comprising a variety of spatial metrics for transect data
    strata_summary: pd.DataFrame
        DataFrame comprising summary features of latitude (INPFC) delimited strata
    settings_dict: dict
        Dictionary containing algorithm arguments that define the proportion of transects resampled
        from the overall dataset within each strata (`transect_sample`) and the number of 
        iterations/realizations used for bootstrapping/resampling (`transect_replicates`).

    Notes
    -----
    This function calculates the stratified summary statistics for biomass within 
    `echopop.survey.stratified_summary()`.
    """ 
    
    # Extract algorithm arguments
    # ---- Number of replicates
    transect_replicates = settings_dict[ 'transect_replicates' ]
    # ---- Transect sampling fraction
    transect_sample = settings_dict[ 'transect_sample' ]
    # ---- Get stratum column name
    stratum_col = settings_dict[ 'stratum_name' ]

    # Calculate the number of transects per stratum
    num_transects_to_sample = (
        np.round(
            strata_summary
            .set_index( stratum_col )[ 'transect_count' ]
            * transect_sample
        ).astype( int )
    )

    # Offset term used for later variance calculation
    sample_offset = np.where( num_transects_to_sample == 1 , 0 , 1 )

    # Calculate effective sample size/degrees of freedom for variance calculation
    sample_dof =  num_transects_to_sample * ( num_transects_to_sample - sample_offset )

    # Get indexed total transect area
    total_transect_area = (
        strata_summary
        .set_index( stratum_col )[ 'transect_area_total' ]
    ) 

    # Get indexed transect distance
    transect_distances = transect_summary.set_index( [ 'transect_num' ] )[ 'transect_distance' ]

    # Get indexed biological value
    biological_values = (
        transect_data
        .groupby( [ 'transect_num' ] )[ settings_dict[ 'variable' ] ]
        .sum( )
    )

    # Get indexed transect numbers
    transect_numbers = (
        transect_summary
        .set_index( stratum_col )[ 'transect_num' ]
    )

    # Pre-allocate the stratum-specific metrics
    # ---- Mean
    mean_arr = np.zeros( [ transect_replicates ,
                          len( total_transect_area.index ) ] )
    # ---- Variance
    variance_arr = np.zeros_like( mean_arr ) 
    # ---- Transect length
    length_arr = np.zeros_like( mean_arr ) 
    # ---- Sum/integrated total across full stratified area/region
    total_arr = np.zeros_like( mean_arr )

    # Iterate across all iterations/realizations
    for j in total_transect_area.index: 

        # Create an index array/matrix containing resampled (without replacement) transect numbers
        transect_numbers_arr = np.array( [ np.random.choice( transect_numbers[j].values , 
                                                             num_transects_to_sample[ j ] , 
                                                             replace = False ) 
                                          for i in range( transect_replicates ) ] )
        
        # Assign the indexed transect distance and biological variables to each transect
        # ---- Transect lengths
        distance_replicates = np.apply_along_axis( transect_array , 
                                                   0 , 
                                                   transect_numbers_arr , 
                                                   transect_distances )
        # -------- Calculate the summed transect length
        length_arr[ : , j - 1 ] = distance_replicates.sum( axis = 1 )
        # ---- Biological variable
        biology_replicates = np.apply_along_axis( transect_array , 
                                                  0 , 
                                                  transect_numbers_arr , 
                                                  biological_values )
        
        # Calculate the stratified weights for along-transect values (within transect weights)
        stratified_weights = array_math( distance_replicates , 
                                         distance_replicates.mean( axis = 1 ) , 
                                         operation = "/" )
        
        # Standardize the biological values by their respective distances
        biology_adjusted = biology_replicates / distance_replicates

        # Calculate the mean transect-length-weighted biological values
        mean_arr[ : , j - 1 ] = (
            biology_replicates * distance_replicates
        ).sum( axis = 1 ) / length_arr[ : , j - 1 ]

        # Sum the total of the biology variable
        total_arr[ : , j - 1 ] = biology_replicates.sum( axis = 1 )
        
        # Calculate the variance of the transect-length-weighted biological values
        # ---- Calculate the sqauared deviation of the mean
        squared_deviation = array_math( biology_adjusted , mean_arr[ : , j - 1 ] , "-" ) ** 2
        # ---- Sum of all weighted squared deviations
        squared_deviation_wgt = ( stratified_weights ** 2 * squared_deviation ).sum( axis = 1 )
        # ---- Compute the variance by incorporating the degrees of freedom
        variance_arr[ : , j - 1 ] = squared_deviation_wgt / sample_dof[ j ]

    # Calculate the area-weight variance of the resulting transect-length-weighted biological 
    # variable    
    # ---- Variance
    replicate_variance = ( total_transect_area.to_numpy( ) ** 2 * variance_arr ).sum( axis = 1 ) 
    # ---- Convert to standard deviation
    replicate_stdev = np.sqrt( replicate_variance )
    # ---- Calculate the replicate mean
    replicate_mean = ( total_transect_area.to_numpy( ) * mean_arr ).sum( axis = 1 )
    # ---- Convert to CV 
    replicate_cv = replicate_stdev / replicate_mean
    
    # Calculate the replicate grand total of the (weighted) biological variable
    # ---- Overall mean of transect integrated sums
    replicate_sum_mean = (
        ( total_transect_area.to_numpy( ) * total_arr ).sum( axis = 1 ) / total_transect_area.sum( )
    )
    # ---- Overall total using the overall mean
    replicate_total = total_arr.sum( axis = 1 )

    # Output the related summary statistics
    # ---- Save the output resampled distributions    
    resampled_distributions = pd.DataFrame( { 'realization': np.arange( 1 , transect_replicates + 1 ) ,
                                               'mean_standardized': replicate_mean ,
                                               'standard_deviation': replicate_stdev ,
                                               'variance': replicate_variance ,
                                               'cv': replicate_cv ,
                                               'mean': replicate_sum_mean ,
                                               'total': replicate_total } )
    # ---- Save the stratified results
    stratified_results = {
        'variable': settings_dict[ 'variable' ] ,
        'ci_percentile': 0.95 ,
        'mean': {
            'weighted_estimate': np.mean( replicate_mean ) ,
            'unweighted_estimate': np.mean( replicate_sum_mean ) ,
            'weighted_confidence_interval': confidence_interval( replicate_mean ) ,
            'unweighted_confidence_interval': confidence_interval( replicate_sum_mean )
        } ,
        'variance': {
            'estimate': np.mean( replicate_variance ) ,
            'confidence_interval': confidence_interval( replicate_variance )
        } ,
        'cv': {
            'estimate': np.mean( replicate_cv ) ,
            'confidence_interval': confidence_interval( replicate_cv )
        } ,
        'total': {
            'estimate': np.mean( replicate_total ) ,
            'confidence_interval': confidence_interval( replicate_total )
        }
    }
    # ---- Return outputs
    return resampled_distributions , stratified_results

def array_math( left_array: np.ndarray ,
                right_array: np.ndarray ,
                operation: str ) -> np.ndarray :
    """
    Helper function that applies arithmetric operations across the dimensions of two arrays.
    
    Parameters
    ----------
    left_array: np.ndarray
        The first array containing data.
    right_array: np.ndarray
        The second array containing data.
    operation: str
        A string value representing a basic arithmetic operation ('+', '+', '-', '/')
    """     
    # Define operations dictionary key
    operations = {
        "+": np.add,
        "-": np.subtract,
        "*": np.multiply,
        "/": np.divide
    }

    # Apply the operation and return the output
    return operations[ operation ]( left_array.transpose( ) , right_array ).transpose( )

def confidence_interval( values: np.ndarray ) -> np.ndarray:
    """
    Calculates the 95% confidence interval (Normal) for a given array    

    Parameters
    ----------
    values: np.ndarray
        An array of values

    Notes
    -----
    This function calculates the 95% confidence interval (assuming a Normal) distribution
    for the bootstrapped stratified samples. This is done as opposed to using the percentile
    method for estimate the intervals.
    """ 
    return np.mean( values ) + np.array( [ -1 , 1 ] ) * st.norm.ppf( 0.975 ) * np.std( values )