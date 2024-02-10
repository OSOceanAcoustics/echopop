import pandas as pd
import geopy.distance

def correct_transect_intervals( dataframe: pd.DataFrame ,
                                threshold: np.float64 = 0.05 ):
    """
    Calculate along-transect intervals and impute erroneous values

    Parameters
    ----------
    dataframe: pd.DataFrame
        DataFrame

    Notes
    -----
    This function calculates the along-track transect interval length and areas.
    It then 'searches' for possible erroneous values at the end of each line 
    and replaces/imputes with alternative lengths/areas.
    """
    return (
            dataframe
                        # Calculate along-transect interval distances
                .pipe( lambda df: df.assign( interval = df[ 'vessel_log_start' ].diff( periods = -1 ).abs( ) ) 
                                    .replace( np.nan , df[ 'vessel_log_end' ].iloc[ -1 ] - df[ 'vessel_log_start' ].iloc[ -1 ] ) )
                # Replace likely erroneous interval lengths associated with the edges of each transect
                .pipe
                ( lambda df: df.assign( median_interval = np.median( df[ 'interval' ] ) )
                                    .assign( interval = lambda x: np.where( np.abs( x[ 'interval' ] - x[ 'median_interval' ] > threshold ) ,
                                                                            x.vessel_log_end - x.vessel_log_start ,
                                                                            x.interval ) ) )
                # Calculate interval area
                .pipe( lambda df: df.assign( interval_area = df[ 'interval' ] * df[ 'transect_spacing' ] ) )                            
                # Keep dataframe tidy by only retaining the necessary columns/variables
                .loc[ : , [ 'latitude' , 'longitude' , 'transect_num' , 'stratum_num' , 'haul_num' , 
                            'interval' , 'interval_area' , 'NASC_all_ages' , 'NASC_no_age1' ] ]
            )


def calculate_start_end_coordinates( group ,
                                     contrast ):
    """
    Calculates latitude/longitude boundary box    

    Parameters
    ----------
    group: pd.DataFrameGroupBy
        Grouped DataFrame
    contrast: List
        Target contrast for grouping

    Notes
    -----
    This function calculates the bounding rectangle surrounding the latitude/longitude values
    for each grouped value (e.g. transect)
    """ 

    return (
        group
        .groupby( contrast )
        .apply( lambda x: pd.Series( { 'minimum_longitude': x['longitude'].min() , 
                                       'maximum_longitude': x['longitude'].max() ,
                                       'center_latitude': x[ 'latitude' ].mean() } ) )
        .reset_index( )
    )

def calculate_transect_distance( dataframe ,
                                 contrast = 'transect_num' ):
    """
    Calculates spatial features of each transect    

    Parameters
    ----------
    dataframe: pd.DataFrame
        DataFrame
    contrast: List
        Target contrast for grouping

    Notes
    -----
    This function calculates the bounding rectangle surrounding the latitude/longitude values
    for each transect and stratum, the average spacing between transects, approximate areas 
    relative to each transect, and the distance for each transect
    """ 
    
    ### Calculate mean transect spacinge
    transect_spacing = dataframe.groupby( contrast )[ 'transect_spacing' ].mean().reset_index()
    
    ### Calculate minimum/maximum longitude, mean latitude, spacing, and area
    return (
            dataframe
            .pipe( lambda df: calculate_start_end_coordinates( df , [ contrast ] ) )
            .assign(transect_distance=lambda x: x.apply( lambda row: geopy.distance.distance( 
                    ( row[ 'center_latitude' ] , row[ 'minimum_longitude' ] ) , 
                    ( row[ 'center_latitude' ] , row[ 'maximum_longitude' ] ) ).nm , 
                    axis=1 ) )
            .merge( transect_spacing , on = [ contrast ] )
            .assign( transect_area = lambda x: x.transect_distance * x.transect_spacing )
    )
