import numpy as np
import pandas as pd
from typing import Union, List
from ..computation.spatial import correct_transect_intervals
from ..computation.operations import group_merge

def filter_species( dataframe_list: Union[List[pd.DataFrame], pd.DataFrame] , 
                    species_id: np.float64 ):
    """
    Filter species in one or multiple biological datasets

    Parameters
    ----------
    dataframe_list: Union[List[pd.DataFrame], pd.DataFrame]
        A list of dataframes or a single dataframe containing biological data/measurements
    species_id: np.float64
        Numeric code representing a particular species of interest
    """   
    
    ### If a single dataframe, convert to a list 
    if isinstance( dataframe_list , pd.DataFrame ):
        dataframe_list = [ dataframe_list ]

    ### Filter out the species-of-interest
    filtered_dataframes = tuple( df[ df.species_id == species_id ] for df in dataframe_list )

    ### Return output
    return filtered_dataframes

def index_sex_weight_proportions( biology_dict: dict ):
    """
    Generate dataframe containing sex-stratified weight proportions    

    Parameters
    ----------
    biology_dict: dict
        Biology data attribute dictionary 
    """     
    
    # Age-stratified weight proportions
    age_stratified_proportions = biology_dict[ 'weight' ][ 'proportions' ][ 'age_weight_proportions_df' ]

    # Age-stratified & sex-indexed weight proportions
    age_sex_stratified_proportions = biology_dict[ 'weight' ][ 'proportions' ][ 'sex_age_weight_proportions_df' ]

    # Concatenate the two to add a 'total' category
    return (
        pd.concat(
            [ ( age_sex_stratified_proportions
                .rename( columns = { 'weight_sex_stratum_proportion': 'weight_proportion' } ) ) ,
              ( age_stratified_proportions.assign( sex = 'total' )
                .rename( columns = { 'weight_stratum_proportion': 'weight_proportion' } ) ) ]
        )
    )

def index_transect_age_sex_proportions( acoustics_dict: dict ,
                                        biology_dict: dict ,
                                        info_strata: pd.DataFrame ):
    """
    Prepares the age- and sex-stratified dataframe for biomass calculation    

    Parameters
    ----------
    acoustics_dict: dict
        Acoustic data attribute dictionary 
    biology_dict: dict
        Biology data attribute dictionary 
    infra_strata: pd.DataFrame
        Dataframe containing strata definitions
    """     
    
    ### Prepare initial dataframes used for calculation population statistics
    # Construct georeferenced dataframe containing NASC data
    nasc_interval_df = correct_transect_intervals( acoustics_dict[ 'nasc' ][ 'nasc_df' ] )

    ### Call additional dataframes needed to merge with the NASC data and subsequently calculate
    ### population-level metrics (and later statistics)
    # Sex-stratum-indexed proportions and average weight
    weight_sex_strata = biology_dict[ 'weight' ][ 'weight_strata_df' ]

    # Stratum-averaged sigma_bs
    sigma_bs_strata = acoustics_dict[ 'sigma_bs' ][ 'strata_mean' ]

    # Adult NASC proportions for each stratum (number)
    # !!! TODO: Currently only uses 'age_1_excluded' -- this should become an argument that toggles
    ## This is not a major issue since both 'NASC_*' columns can be pivoted to create a single 
    ## NASC column so the column name does not have to be hard-coded. This could then correspond to
    ## the configuration settings in some way, or this may be where the argument comes into play where
    ## the dataframe can be simply filtered based on the input/selection.
    # between excluding and including age-1 fish
    nasc_number_proportions = (
        biology_dict[ 'weight' ][ 'proportions' ][ 'age_proportions_df' ]
    )

    # Adult NASC proportions for each stratum (weight)
    nasc_weight_proportions = (
        biology_dict[ 'weight' ][ 'proportions' ][ 'age_weight_proportions_df' ]
    )

    ### Consolidate dataframes that will be added into a list
    dataframes_to_add = [ nasc_interval_df , sigma_bs_strata , weight_sex_strata , nasc_number_proportions , 
                          nasc_weight_proportions ]
    

    ## Merge the relevant dataframes
    return (
        nasc_interval_df
        # Merge stratum information ( join = 'outer' since missing values will be filled later on)
        .merge( info_strata , on = [ 'stratum_num' , 'haul_num' ] , how = 'outer' )
        # Drop unused hauls
        .dropna( subset = 'transect_num' )
        # Fill NaN w/ 0's for 'fraction_hake'
        .assign( fraction_hake = lambda x: x[ 'fraction_hake' ].fillna( 0 ) )
        # Group merge
        .group_merge( dataframes_to_add = dataframes_to_add , inner_on = 'age' , outer_on = 'stratum_num' )
        )

def sum_strata_weight( haul_data: pd.DataFrame ,
                       specimen_data: pd.DataFrame ):
    """
    Sum haul weights from different datasets    

    Parameters
    ----------
    haul_data: pd.DataFrame
        Dataframe containing haul catch weight data
    specimen_data: pd.DataFrame
        Dataframe containing specimen weights
    """ 

    ### Calculate station weights -- stratified by strata
    weight_strata_station = pd.concat(
         [
              # Haul data 
              haul_data
              .groupby( [ 'stratum_num' ] )
              .apply( lambda df: pd.Series( { 'stratum_weight': df.haul_weight.sum( ) } ) )
              .reset_index( )
              .assign( station = 1 ) ,
              # Specimen data
              specimen_data 
              .groupby( [ 'stratum_num' ] )
              .apply( lambda df: pd.Series( { 'stratum_weight': df.weight.sum( ) } ) )
              .reset_index()
              .assign( station = 2 )
         ]
    )

    ### Sum weights for each stratum
    weight_strata = (
            weight_strata_station
            .groupby( [ 'stratum_num' ] )
            .apply( lambda x: pd.Series( { 'weight_stratum_total': x.stratum_weight.sum( ) } ) )            
            .reset_index( )
        )
    
    ### Carriage return
    return weight_strata_station , weight_strata

def sum_strata_length_age_sex_weight( haul_data: pd.DataFrame ,
                                      specimen_data: pd.DataFrame ,
                                      length_intervals: np.ndarray ,
                                      age_intervals: np.ndarray ):
    
    ### Sum weights from each station and strata 
    weight_strata_station , weight_strata = sum_strata_weight( haul_data ,
                                                               specimen_data )
    
    ### Calculate summed length-age-sex weight bins across strata
    weight_length_age_sex_stratum = (
        specimen_data
        .dropna( how = 'any' )
        .loc[ lambda x: x.sex != 'unsexed' ]
        .pipe( lambda df: pd.concat( [ df , 
                                       df.copy( ).assign( sex = 'all' ) ] ) )
        .bin_variable( age_intervals , 'age' )
        .bin_variable( length_intervals , 'length' )      
        .groupby( [ 'stratum_num' , 'age_bin' , 'length_bin' , 'sex' ] )
        .apply( lambda df: pd.Series( { 'summed_weight_all': df.weight.sum( ) ,
                                        'summed_weight_adult': df.loc[ df.age > 1 ].weight.sum( ) } ) ) 
        .reset_index( )
        .replace( np.nan , 0 )
        .assign( total_weight_sex_all = lambda df: df.groupby( [ 'stratum_num' , 'sex' ] )[ 'summed_weight_all' ].transform( sum ) ,
                 total_weight_sex_adult = lambda df: df.groupby( [ 'stratum_num' , 'sex' ] )[ 'summed_weight_adult' ].transform( sum ) ,
                 proportion_weight_all = lambda df: df.summed_weight_all / df.total_weight_sex_all ,
                 proportion_weight_adult = lambda df: df.summed_weight_adult / df.total_weight_sex_adult )            
    )

    ### Carriage return
    return weight_strata , weight_strata_station , weight_length_age_sex_stratum

def normalize_length_age_sex_weight_proportions( weight_length_age_sex_stratum: pd.DataFrame ,
                                                 weight_strata: pd.DataFrame ):
    
    ### Normalize the age-length-sex indexed proportions
    dist_weight_sum = ( 
       weight_length_age_sex_stratum 
        .merge( weight_strata , on = [ 'stratum_num' ] ) 
        .groupby( [ 'stratum_num' , 'sex' ] ) 
        .apply( lambda df: pd.Series( {  
            'proportion_normalized': ( df.proportion_weight_all * ( df.summed_weight_all / df.weight_stratum_total ).sum( ) ).sum( ) 
        } ) ) 
        .reset_index( ) 
    ) 

    ### Carriage return
    return dist_weight_sum

def calculate_aged_proportions( weight_length_age_sex_stratum: pd.DataFrame ,
                                weight_strata: pd.DataFrame ):
    
    ### Calculate aged / unaged proportions 
    aged_proportions = ( 
        weight_length_age_sex_stratum
        .loc[ lambda df: df.sex != 'all' ] 
        .groupby( [ 'stratum_num' ] ) 
        .apply( lambda df: pd.Series( { 'weight_aged_total': df.summed_weight_all.sum( ) } ) )  
        .reset_index( ) 
        .merge( weight_strata , on = [ 'stratum_num' ] )  
        .assign( proportion_aged_total = lambda x: x.weight_aged_total / x.weight_stratum_total , 
                 proportion_unaged_total = lambda x: np.round( 1.0 - x.proportion_aged_total , decimals = 10 ) ) 
    ) 

    ### Carriage return
    return aged_proportions

def normalize_haul_sex_weights( length_weight_df: pd.DataFrame ,
                                length_df: pd.DataFrame ,
                                weight_strata: pd.DataFrame ,
                                weight_strata_station: pd.DataFrame ,
                                aged_proportions: pd.DataFrame ,
                                length_intervals: np.ndarray ):
    
    ### Calculate interpolated weights based on length bins for each sex per haul 
    length_weight_fit = (
        length_weight_df
        .copy( )
        .assign( length_bin_value  = lambda x: x[ 'length_bin' ].apply( lambda y: y.mid ) )
    )

     # Sum haul weights per sex per stratum 
    haul_weights = ( 
        length_df 
        .bin_variable( length_intervals , 'length' ) 
        .loc[ lambda x: x.group == 'sexed' ] 
        .pipe( lambda df: pd.concat( [ df , df.assign( sex = 'all' ) ] ) )   
        .merge( length_weight_fit , on = [ 'sex' , 'length_bin' ] ) 
        .groupby( [ 'stratum_num' , 'haul_num' , 'sex' ] ) 
        .apply( lambda df: pd.Series( { 'weight_interp': ( np.interp( df.length,  
                                                                    length_weight_fit.loc[ lambda x: x.sex.isin( df.sex ) ][ 'length_bin_value' ] ,  
                                                                    length_weight_fit.loc[ lambda x: x.sex.isin( df.sex ) ][ 'weight_modeled' ] ) * 
                                                            df.length_count ).sum( ) } ) ) 
        .groupby( [ 'stratum_num' , 'sex' ] ) 
        .apply( lambda df: pd.Series( { 'summed_haul_weights': df.weight_interp.sum( ) } ) ) 
        .reset_index( ) 
    ) 
  
    ### Normalize haul weights (Station 1) 
    haul_sex_weights_normalized = ( 
        haul_weights 
        .merge( weight_strata_station.loc[ lambda x: x.station == 1 ] ,  
                on = 'stratum_num' ) 
        .assign( weight_normalized_station_1 = lambda df: ( 
            df 
            .groupby( [ 'stratum_num' ] ) 
            .apply( lambda strata: strata[ 'stratum_weight' ] * strata[ 'summed_haul_weights' ] / 
                                    ( df[ ( df.stratum_num == strata.name ) & ( df.sex == 'male' ) ][ 'summed_haul_weights' ].iloc[ 0 ] + 
                                      df[ ( df.stratum_num == strata.name ) & ( df.sex == 'female' ) ][ 'summed_haul_weights' ].iloc[ 0 ] ) ) 
            .reset_index( drop = True ) ) ) 
        .merge( weight_strata , on = [ 'stratum_num' ] ) 
        .assign( proportion_normalized_station_1 = lambda df: ( 
            df.weight_normalized_station_1 / df.weight_stratum_total 
        ) ) 
        .assign( sex_proportion = lambda df: ( 
            df 
            .groupby( [ 'stratum_num' ] ) 
            .apply( lambda strata: strata[ 'proportion_normalized_station_1' ] / 
                                    ( df[ ( df.stratum_num == strata.name ) & ( df.sex == 'male' ) ][ 'proportion_normalized_station_1' ].iloc[ 0 ] + 
                                      df[ ( df.stratum_num == strata.name ) & ( df.sex == 'female' ) ][ 'proportion_normalized_station_1' ].iloc[ 0 ] ) ) 
            .reset_index( drop = True ) ) ) 
        .merge( aged_proportions.loc[ : , [ 'stratum_num' , 'proportion_unaged_total' ] ] ,  
                on = [ 'stratum_num' ] )         
        .assign( proportion_unaged_sex = lambda x: x.sex_proportion * x.proportion_unaged_total ) 
    ) 

    ### Carriage return
    return haul_sex_weights_normalized