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

    ### Process biological datasets
    # ---- Process haul data ( Station 2 - unsexed - unaged )
    haul_strata_weight = (
        haul_data
        .groupby( 'stratum_num' )[ 'haul_weight' ]
        .sum( )
        .to_frame( 'stratum_weight' )
        .reset_index( )
        .assign( group = 'unaged' )
    ) 

    # ---- Process specimen data ( Station 1 - sexed - aged )
    specimen_strata_weight = (
        specimen_data
        .groupby( 'stratum_num' )[ 'weight' ]
        .sum( )
        .to_frame( 'stratum_weight' )
        .reset_index( )
        .assign( group = 'aged' )
    )

    ### Merge the two dataframes
    # Yields the summed weights per stratum and station (or aged/unaged designation)
    weight_strata_aged_unaged = pd.concat( [ haul_strata_weight , specimen_strata_weight ] )
 
    ### Sum weights for each stratum
    weight_strata = (
        weight_strata_aged_unaged
        .groupby( "stratum_num" )[ "stratum_weight" ]
        .sum( )
        .to_frame( 'weight_stratum_all' )
        .reset_index( )
    )
    
    ### Carriage return
    return weight_strata_aged_unaged , weight_strata

def compute_index_aged_weight_proportions( specimen_data: pd.DataFrame ,
                                           length_intervals: np.ndarray ,
                                           age_intervals: np.ndarray ):
    """
    Calculate length-age binned weight proportions    

    Parameters
    ----------
    specimen_data: pd.DataFrame
        Dataframe containing specimen weights
    length_intervals: np.ndarray
        Array containing length bins/intervals
    age_intervals: np.ndarray
        Array containing age bins/intervals
    """ 

    ### Process the specimen data 
    # ---- Drop unaged and unsexed fish
    # ==== !!! TODO: pending what FEAT says, weights associated with 
    # ==== missing ages should be added into the 'unaged' category. 
    # ==== This would further mean apportioning these into `weight_strata_aged_uanged`
    specimen_data_filtered = (
        specimen_data[ specimen_data.sex != 'unsexed' ]
        .dropna( how = 'any' , subset = 'age' )
    ) 
    
    # ---- Bin length and age measurements
    specimen_data_filtered = (
        specimen_data_filtered
        # ---- Bin length
        .bin_variable( length_intervals , 'length' )
        # ---- Age bin
        .bin_variable( age_intervals , 'age' )
    )

    ### Sum weights within each length and age bin for each sex within each stratum
    specimen_binned_weight = (
        specimen_data_filtered
        # ---- Group weight summations across stratum/species/sex/length/age
        .groupby( [ 'stratum_num' , 'species_id' , 'sex' , 'length_bin' , 'age_bin' ] )
        # ---- Sum the weights 
        .apply( lambda df: pd.Series( { 'weight_all': df.weight.sum( ) ,
                                        'weight_adult': df.loc[ df.age > 1 ].weight.sum( ) } ) ) 
        # ---- Fill empty/non-existent values with 0's
        .fillna( 0 )
        .reset_index( )
    )

    ### Calculate the relative weight proportions of each length-age bin for each sex within each stratum
    proportions_weight_length_age_sex = (
        specimen_binned_weight
        # ---- Calculate total sex-specific weights for each stratum
        .assign( total_weight_sex_all = lambda df: df.groupby( [ 'stratum_num' , 'species_id' , 'sex' ] )[ 'weight_all' ].transform( sum ) ,
                 total_weight_sex_adult = lambda df: df.groupby( [ 'stratum_num' , 'species_id' , 'sex' ] )[ 'weight_adult' ].transform( sum ) )
        # ---- Calculate the weight proportions within each sex: from Matlab --> Len_age_key_wgt_*n
        .assign( proportion_weight_sex_all = lambda df: df.weight_all / df.total_weight_sex_all ,
                 proportion_weight_sex_adult = lambda df: df.weight_adult / df.total_weight_sex_adult )
    )

    ### Fill empty/non-existent values with 0's
    proportions_weight_length_age_sex[ 'proportion_weight_sex_all' ] = (
        proportions_weight_length_age_sex[ 'proportion_weight_sex_all' ].fillna( 0 )
    )

    proportions_weight_length_age_sex[ 'proportion_weight_sex_adult' ] = (
        proportions_weight_length_age_sex[ 'proportion_weight_sex_adult' ].fillna( 0 )
    )

    ### Return output
    return proportions_weight_length_age_sex

def compute_summed_aged_proportions( proportions_weight_length_age_sex: pd.DataFrame ,
                                     weight_strata: pd.DataFrame ):
    """
   Compute the aged proportions across all fish and specific sexes for each
   length and age bin  

    Parameters
    ----------
    proportions_weight_length_age_sex: pd.DataFrame
        Dataframe containing sexed weight proportions distributed across
        each age and length bin
    weight_strata: pd.DataFrame
        Dataframe contained summed weights of both aged and unaged fish
    """   

    ### Caclulate the summed weights for all and adult fish 
    # ---- Group by `stratum_num` and `sex` and sum
    # ---- Also initializes the correct shape for the output dataframe
    aged_sex_weights = (
        proportions_weight_length_age_sex
        .groupby( [ 'stratum_num' , 'sex' ] )
        .agg(
            weight_aged_sex_all = ( 'weight_all' , 'sum' ) ,
            weight_aged_sex_adult = ( 'weight_adult' , 'sum' )
        )
        .reset_index( )
    )

    ### Calculate the weight proportions
    # ---- Merge `aged_sex_weights` with `weight_strata` to get strata weights
    aged_sex_proportions = pd.merge( aged_sex_weights , 
                                     weight_strata , 
                                     on = [ 'stratum_num' ] , 
                                     how = 'left' )

    # ---- Proportions for all fish
    aged_sex_proportions[ 'proportion_weight_all' ] = (
        aged_sex_proportions.weight_aged_sex_all / aged_sex_proportions.weight_stratum_all
    )

    # ---- Proportions for adult fish
    aged_sex_proportions[ 'proportion_weight_adult' ] = (
        aged_sex_proportions.weight_aged_sex_adult / aged_sex_proportions.weight_stratum_all
    )

    # ---- Fill empty/NaN values with 0's and drop unused columns
    aged_sex_proportions = aged_sex_proportions.fillna( 0 ).filter( regex = '^(?!weight_).*' )

    ### Calculate aged proportions
    aged_proportions = (
        aged_sex_proportions
        .groupby( [ 'stratum_num' ] )
        # ---- Sum proportions from each sex for each stratum
        .agg( {
            'proportion_weight_all': 'sum' ,
            'proportion_weight_adult': 'sum'
        } )
        .reset_index( )
    )

    ### Return output
    return aged_sex_proportions , aged_proportions

def distribute_aged_weight_proportions( proportions_weight_length_age_sex: pd.DataFrame ,
                                        aged_sex_proportions: pd.DataFrame ):
    """
    Distribute overall weight proportions across each sex and age/length bins   

    Parameters
    ----------
    proportions_weight_length_age_sex: pd.DataFrame
        Dataframe containing sexed weight proportions distributed across
        each age and length bin
    aged_sex_proportions: pd.DataFrame
        Dataframe contained weight proportions of sexed fish for aged fish
    """     

    ### Calculate the normalized age-length-sex aged fish weight proportions
    # ---- Merge `proportions_weight_length_age_sex` with `aged_sex_proportions`
    distributed_aged_weight_proportions = pd.merge( proportions_weight_length_age_sex , 
                                                    aged_sex_proportions , 
                                                    on = [ 'stratum_num' , 'sex' ] , 
                                                    how = 'left' )
    
    # ---- Normalized weight proportions for all fish
    distributed_aged_weight_proportions[ 'normalized_proportion_weight_all' ] = (
        distributed_aged_weight_proportions.proportion_weight_sex_all * distributed_aged_weight_proportions.proportion_weight_all
    )

    # ---- Normalized weight proportions for jus adult fish
    distributed_aged_weight_proportions[ 'normalized_proportion_weight_adult' ] = (
        distributed_aged_weight_proportions.proportion_weight_sex_all * distributed_aged_weight_proportions.proportion_weight_adult
    )

    # ---- Remove unnecessary columns and return the dataframe
    return distributed_aged_weight_proportions.filter( regex = '^(?!weight_|total_weight_|proportion_).*' )

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

# def normalize_haul_sex_weights( length_weight_df: pd.DataFrame ,
#                                 length_df: pd.DataFrame ,
#                                 weight_strata: pd.DataFrame ,
#                                 weight_strata_station: pd.DataFrame ,
#                                 aged_proportions: pd.DataFrame ,
#                                 length_intervals: np.ndarray ):
    
#     ### Calculate interpolated weights based on length bins for each sex per haul 
#     length_weight_fit = (
#         length_weight_df
#         .copy( )
#         .assign( length_bin_value  = lambda x: x[ 'length_bin' ].apply( lambda y: y.mid ) )
#     )

#      # Sum haul weights per sex per stratum 
#     haul_weights = ( 
#         length_df 
#         .bin_variable( length_intervals , 'length' ) 
#         .loc[ lambda x: x.group == 'sexed' ] 
#         .pipe( lambda df: pd.concat( [ df , df.assign( sex = 'all' ) ] ) )   
#         .merge( length_weight_fit , on = [ 'sex' , 'length_bin' ] ) 
#         .groupby( [ 'stratum_num' , 'haul_num' , 'sex' ] ) 
#         .apply( lambda df: pd.Series( { 'weight_interp': ( np.interp( df.length,  
#                                                                     length_weight_fit.loc[ lambda x: x.sex.isin( df.sex ) ][ 'length_bin_value' ] ,  
#                                                                     length_weight_fit.loc[ lambda x: x.sex.isin( df.sex ) ][ 'weight_modeled' ] ) * 
#                                                             df.length_count ).sum( ) } ) ) 
#         .groupby( [ 'stratum_num' , 'sex' ] ) 
#         .apply( lambda df: pd.Series( { 'summed_haul_weights': df.weight_interp.sum( ) } ) ) 
#         .reset_index( ) 
#     ) 
  
#     ### Normalize haul weights (Station 1) 
#     haul_sex_weights_normalized = ( 
#         haul_weights 
#         .merge( weight_strata_station.loc[ lambda x: x.station == 1 ] ,  
#                 on = 'stratum_num' ) 
#         .assign( weight_normalized_station_1 = lambda df: ( 
#             df 
#             .groupby( [ 'stratum_num' ] ) 
#             .apply( lambda strata: strata[ 'stratum_weight' ] * strata[ 'summed_haul_weights' ] / 
#                                     ( df[ ( df.stratum_num == strata.name ) & ( df.sex == 'male' ) ][ 'summed_haul_weights' ].iloc[ 0 ] + 
#                                       df[ ( df.stratum_num == strata.name ) & ( df.sex == 'female' ) ][ 'summed_haul_weights' ].iloc[ 0 ] ) ) 
#             .reset_index( drop = True ) ) ) 
#         .merge( weight_strata , on = [ 'stratum_num' ] ) 
#         .assign( proportion_normalized_station_1 = lambda df: ( 
#             df.weight_normalized_station_1 / df.weight_stratum_total 
#         ) ) 
#         .assign( sex_proportion = lambda df: ( 
#             df 
#             .groupby( [ 'stratum_num' ] ) 
#             .apply( lambda strata: strata[ 'proportion_normalized_station_1' ] / 
#                                     ( df[ ( df.stratum_num == strata.name ) & ( df.sex == 'male' ) ][ 'proportion_normalized_station_1' ].iloc[ 0 ] + 
#                                       df[ ( df.stratum_num == strata.name ) & ( df.sex == 'female' ) ][ 'proportion_normalized_station_1' ].iloc[ 0 ] ) ) 
#             .reset_index( drop = True ) ) ) 
#         .merge( aged_proportions.loc[ : , [ 'stratum_num' , 'proportion_unaged_total' ] ] ,  
#                 on = [ 'stratum_num' ] )         
#         .assign( proportion_unaged_sex = lambda x: x.sex_proportion * x.proportion_unaged_total ) 
#     ) 

#     ### Carriage return
#     return haul_sex_weights_normalize

def compute_index_unaged_number_proportions( length_data: pd.DataFrame ,
                                             length_intervals: np.ndarray ):
    """
    Calculate length binned weight proportions among unaged fish  

    Parameters
    ----------
    length_data: pd.DataFrame
        Dataframe containing length data measured from unaged fish
    length_intervals: np.ndarray
        Array containing length bins/intervals
    """ 

    ### Calculate number proportion
    # ---- Remove unsexed 
    length_data_filtered = length_data[ length_data.sex != 'unsexed' ].copy( )

    # ---- Bin length measurements
    length_data_filtered = (
        length_data_filtered
        # ---- Bin length
        .bin_variable( length_intervals , 'length' ) 
    )

    # ---- Sum the number of individuals within each bin 
    proportions_unaged_length = (
        length_data_filtered
        # ---- Group number summations across stratum/species/sex/length
        .groupby( [ 'stratum_num' , 'species_id' , 'length_bin' ] )
        # ---- Sum count
        .agg( number_all = ( 'length_count' , 'sum' ) )
        # ---- Fill empty/non-existent values with 0's
        .fillna( 0 )
        .reset_index( )
    )

    ### Calculate the number proportions
    # --- Stratum total counts
    proportions_unaged_length[ 'stratum_number_all' ] = (
        proportions_unaged_length.groupby( [ 'stratum_num' ] )[ 'number_all' ].transform( sum )
    )

    # ---- Proportions of each sex-length bin pair relative to `stratum_number_all`
    proportions_unaged_length[ 'proportion_number_all' ] = (
        proportions_unaged_length.number_all / proportions_unaged_length.stratum_number_all
    ).fillna( 0 )

    ### Filter out unnecessary columns and return output
    return proportions_unaged_length.filter( regex = '^(?!number_|stratum_number_).*' )

def compute_summed_unaged_weight_proportions( proportions_unaged_length: pd.DataFrame ,
                                              length_weight_df: pd.DataFrame ):
    """
   Compute the unaged weight proportions across all fish and specific sexes for each
   length bin  

    Parameters
    ----------
    proportions_unaged_length: pd.DataFrame
        Dataframe containing the total number proportions distributed across
        each length bin
    length_weight_df: pd.DataFrame
        Dataframe containing the modeled weights for each sex and length bin
        fit from the length-weight regression (fit for all animals)
    """   

    ### Extract length-weight regression results calculated for all animals
    length_weight_all = length_weight_df[ length_weight_df.sex == 'all' ][ [ 'length_bin' , 'weight_modeled' ] ]

    ### Calculate the weight proportion of each length bin that is 'weighted' by the number proportions
    # ---- Merge `proportions_unaged_length_sex` and `length_weight_all`
    proportions_unaged_weight_length = pd.merge( proportions_unaged_length ,
                                                 length_weight_all ,
                                                 on = [ 'length_bin' ] ,
                                                 how = 'left' )
    
    # ---- Calculate the weight proportion (`w_ln_all_array` in the original Matlab code)
    proportions_unaged_weight_length[ 'weight_proportion_all' ] = (
        proportions_unaged_weight_length.weight_modeled * proportions_unaged_weight_length.proportion_number_all
    )

    ### Sum weight proportions for each strata to get the weight per length bin
    # ---- Calculate the summed weights per stratum (`w_ln_array_sum` in the original Matlab code)
    proportions_unaged_weight_length[ 'weight_stratum_proportion' ] = (
       proportions_unaged_weight_length.groupby( [ 'stratum_num' ] )[ 'weight_proportion_all' ].transform( sum )
    )
    
    # ---- Calculate the weight per length distribution bin (`w_ln_all_N` in the original Matlab code)
    proportions_unaged_weight_length[ 'proportion_weight_length' ] = (
        proportions_unaged_weight_length.weight_proportion_all / proportions_unaged_weight_length.weight_stratum_proportion
    )

    ### Drop unnecessary columns and return output
    return proportions_unaged_weight_length.filter( regex = '^(?!weight_|proportion_number_).*' )

def compute_unaged_sex_proportions( length_data: pd.DataFrame ,
                                    length_intervals: pd.DataFrame ,
                                    length_weight_df: pd.DataFrame ,
                                    weight_strata: pd.DataFrame ,
                                    weight_strata_aged_unaged: pd.DataFrame ):
    """
    Calculate the normalized weight proportions for unaged sexed fish

    Parameters
    ----------
    length_data: pd.DataFrame
        Dataframe containing length data measured from unaged fish
    length_intervals: np.ndarray
        Array containing length bins/intervals
    length_weight_df: pd.DataFrame
        Dataframe containing the modeled weights for each sex and length bin
        fit from the length-weight regression (fit for all animals)
    weight_strata: pd.DataFrame
        Summed haul weights for each stratum
    weight_strata_aged_uanged: pd.DataFrame
        Summed haul weights for specifically aged and unaged fish across strata
    """ 

    ### Calculate sexed weights within each stratum
    # ---- Remove unsexed 
    length_data_filtered = length_data[ length_data.sex != 'unsexed' ].copy( )

    # ---- Extract sexed length-weight regression fits
    length_weight_sex = length_weight_df[ length_weight_df.sex != 'all' ][ [ 'length_bin' , 'sex' , 'weight_modeled' ] ]

    # ---- Bin length measurements
    length_data_filtered = (
        length_data_filtered
        # ---- Bin length
        .bin_variable( length_intervals , 'length' ) 
    )

    # ---- Extract bin length values (i.e. mid)
    length_data_filtered[ 'length_bin_value' ] = (
        length_data_filtered[ 'length_bin' ].apply( lambda x: x.mid )
    )

    # ---- Merge `length_data_filtered ` and `length_weight_sex`
    sexed_unaged_model_weights = pd.merge( length_data_filtered ,
                                           length_weight_sex ,
                                           on = [ 'length_bin' , 'sex' ] ,
                                           how = 'left' )
    
    # ---- Interpolate weights over length values for each sex (weighted by `length_count`)
    interpolated_weights = (
        sexed_unaged_model_weights
        .groupby( [ 'stratum_num' , 'sex' ] ) 
        .apply( lambda df: pd.Series( {
                'weight_interp': ( np.interp( df.length ,
                                              df.length_bin_value ,
                                              df.weight_modeled ) * df.length_count ).sum( )
            } ) )
        .reset_index( )
    )

    ### Calculate the unaged weight proportion relative to the haul catches
    # ---- Parse out unaged haul totals
    weight_strata_unaged = weight_strata_aged_unaged[ weight_strata_aged_unaged.group == 'unaged' ]

    # ---- Merge the unaged weight strata with the interpolated weights
    proportions_unaged_weight_sex = pd.merge( interpolated_weights ,
                                              weight_strata_unaged ,
                                              on = [ 'stratum_num' ] ,
                                              how = 'left' )
    
    # ---- Sum `weight_interp` over each stratum
    proportions_unaged_weight_sex[ 'stratum_sex_weight_interp' ] = (
        proportions_unaged_weight_sex.groupby( ['stratum_num' ] )[ 'weight_interp' ].transform( sum ) 
    )

    # ---- Calculate the normalized weights for each sex
    proportions_unaged_weight_sex[ 'stratum_sex_weight_normalized' ] = (
        proportions_unaged_weight_sex.stratum_weight *
        proportions_unaged_weight_sex.weight_interp / proportions_unaged_weight_sex.stratum_sex_weight_interp
    )

    # ---- Merge `proportions_unaged_weight_sex` with `weight strata`
    proportions_unaged_weight_sex_normalized = pd.merge( proportions_unaged_weight_sex ,
                                                         weight_strata ,
                                                         on = [ 'stratum_num' ] ,
                                                         how = 'left' )

    # ---- Sum the normalized stratum sex weight across each stratum
    proportions_unaged_weight_sex[ 'proportions_weight_sex' ] = (
        proportions_unaged_weight_sex_normalized.stratum_sex_weight_normalized /
        proportions_unaged_weight_sex_normalized.weight_stratum_all
    )

    # ---- Sum 'proportions_weight_sex'
    proportions_unaged_weight_sex[ 'proportions_weight_sex_total' ] = (
        proportions_unaged_weight_sex.groupby( [ 'stratum_num' ] )[ 'proportions_weight_sex' ].transform( sum )
    )

    # ---- Calculate the final sex proportion
    proportions_unaged_weight_sex[ 'proportion_weight_sex' ] = (
        proportions_unaged_weight_sex.proportions_weight_sex / proportions_unaged_weight_sex.proportions_weight_sex_total
    )

    ### Remove unnecessary columns and return output
    return proportions_unaged_weight_sex[ [ 'stratum_num' , 'sex' , 'proportion_weight_sex' ] ]