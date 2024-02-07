import pandas as pd
import geopy.distance

def calculate_bounds( group ,
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
    
    transect_spacing_area = dataframe.groupby( contrast )[ 'transect_spacing' ].mean().reset_index()
    
    return (
            dataframe
            .pipe( lambda df: calculate_bounds( df , [ contrast ] ) )
            .assign(transect_distance=lambda x: x.apply( lambda row: geopy.distance.distance( 
                    ( row[ 'center_latitude' ] , row[ 'minimum_longitude' ] ) , 
                    ( row[ 'center_latitude' ] , row[ 'maximum_longitude' ] ) ).nm , 
                    axis=1 ) )
            .merge( transect_spacing_area , on = [ contrast ] )
            .assign( transect_area = lambda x: x.transect_distance * x.transect_spacing )
    )
