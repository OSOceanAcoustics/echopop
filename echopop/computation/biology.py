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

def transect_number_proportions( aged_data: pd.DataFrame ,
                                 unaged_data: pd.DataFrame , ) :
    """_summary_

    Parameters
    ----------
    aged_data: pd.DataFrame
        DataFrame containing aged specimen data
    unaged_data: pd.DataFrame
        Dataframe containing unaged specimen data
    """
    # Bin counts by sex
    # ---- Aged
    aged_number_distribution = (
        pd.concat( [ aged_data , aged_data.assign( sex = 'all' ) ] )
        .count_variable( contrasts = [ 'stratum_num' , 'sex' , 'length_bin' , 'age_bin' ] , 
                         variable = 'length' , 
                         fun = 'size' )
    )

    # ---- Filter out unsexed data for parallel number counts and drop any NA's
    aged_number_distribution_filtered = (
        pd.concat( [ aged_data[ aged_data.sex != 'unsexed' ] ,
                     aged_data[ aged_data.sex != 'unsexed' ].assign( sex = 'all' ) ] )
        .dropna( subset = [ 'length' , 'age' , 'weight' ] )
        .count_variable( contrasts = [ 'stratum_num' , 'sex' , 'length_bin' , 'age_bin' ] , 
                         variable = 'length' , 
                         fun = 'size' )        
    ) 
    # ---- Unaged
    unaged_number_distribution = (
        pd.concat( [ unaged_data ,
                     unaged_data.assign( sex = 'all' ) ] )
        .count_variable( contrasts = [ 'stratum_num' , 'sex' , 'length_bin' ] , 
                         variable = 'length_count' , 
                         fun = 'sum' )
    )

    # Calculate total numbers among aged samples
    # ---- Collect unique strata
    strata_unique = pd.DataFrame(
        {
            'stratum_num': np.unique( np.concatenate( [ aged_number_distribution.stratum_num ,
                                                        unaged_number_distribution.stratum_num ] ) ) ,
        } ,
    )  
    # ---- Collect unique sex
    sex_unique = pd.DataFrame(
        {
            'sex': np.unique( np.concatenate( [ aged_number_distribution.sex ,
                                                unaged_number_distribution.sex ] ) ) ,
        } ,
    )
    # ---- Initialize dataframe
    count_total = (
        strata_unique
        .merge( sex_unique , how = 'cross' )
        .set_index( [ 'stratum_num' , 'sex' ] )
    )
    # ---- Aged (overall)
    count_total[ 'total_overall_aged' ] = (
        aged_number_distribution
        .groupby( [ 'stratum_num' , 'sex' ] )[ 'count' ].sum( )
    )
    # ---- Aged (filtered)
    count_total[ 'total_filtered_aged' ] = (
        aged_number_distribution_filtered
        .groupby( [ 'stratum_num' , 'sex' ] )[ 'count' ].sum( )
    )
    # ---- Unaged
    count_total[ 'total_overall_unaged' ] = (
        unaged_number_distribution
        .groupby( [ 'stratum_num' , 'sex' ] )[ 'count' ].sum( )
    )
    # ---- Fill NaN
    count_total.fillna( 0 , inplace = True )
    count_total = count_total.reset_index( ).set_index( 'stratum_num' )
    # ---- Grand totals
    count_total[ 'total_overall' ] = (
        count_total[ count_total.sex == 'all' ].total_filtered_aged
         + count_total[ count_total.sex == 'all' ].total_overall_unaged
    )
    count_total = count_total.reset_index( )
    # Calculate number proportions
    # ---- Aged (number distributed over age and length)
    aged_number_proportion = (
        aged_number_distribution_filtered[ aged_number_distribution_filtered.sex.isin( [ 'male' , 'female' , 'all' ] ) ]        
        .merge( count_total[ [ 'stratum_num' , 'sex' , 'total_filtered_aged' , 'total_overall' ] ] , 
                on = [ 'stratum_num' , 'sex' ] )
    )
    # -------- Aged-specific proportion (sex)
    # aged_number_proportion[ 'proportion_number_sex_aged' ] = (
    #     aged_number_proportion[ 'count' ] / aged_number_proportion[ 'total_filtered_aged' ]
    # )
    # -------- Aged-specific proportion
    aged_number_proportion[ 'proportion_number_aged' ] = (
        aged_number_proportion[ 'count' ] / aged_number_proportion[ 'total_filtered_aged' ]
    )

    # -------- Overall proportion
    aged_number_proportion[ 'proportion_number_overall_aged' ] = (
        aged_number_proportion[ 'count' ] / aged_number_proportion[ 'total_overall' ]
    )

    # -------- Gather aged (sexed) number proportions
    sex_number_proportions = (
        aged_number_proportion
        .groupby( [ 'stratum_num' , 'sex' ] , observed = False )[ 'proportion_number_overall_aged' ].sum( )
        .reset_index( )
    )
    # ---- Unaged (number distributed over length)
    unaged_number_proportion = (
        unaged_number_distribution
        .merge( count_total[ [ 'stratum_num' , 'sex' , 'total_overall_unaged' , 'total_overall' ] ] , 
                on = [ 'stratum_num' , 'sex' ] )
    )   
    
    # -------- Unaged-specific proportion
    unaged_number_proportion[ 'proportion_number_unaged' ] = (
        unaged_number_proportion[ 'count' ] / unaged_number_proportion[ 'total_overall_unaged' ]
    )

    # -------- Overall proportion
    unaged_number_proportion[ 'proportion_number_overall_unaged' ] = (
        unaged_number_proportion[ 'count' ] / unaged_number_proportion[ 'total_overall' ]
    )
    # -------- Gather unaged (sexed) number proportions
    # ------------ Merge
    sex_number_proportions = (
        sex_number_proportions
        .merge( unaged_number_proportion
                .groupby( [ 'stratum_num' , 'sex' ] )[ 'proportion_number_overall_unaged' ]
                .sum( ).reset_index( ) , how = 'outer' ).fillna( 0.0 )
    )
    # ------------ Sum overall total across aged and unaged samples
    sex_number_proportions[ 'proportion_number_overall' ] = (
        sex_number_proportions.proportion_number_overall_aged
        + sex_number_proportions.proportion_number_overall_unaged
    )

    # Calculate overall number proportions across age
    age_number_proportions = (
        aged_number_proportion
        .groupby( [ 'stratum_num' , 'age_bin' ] , observed = False )[ 'proportion_number_aged' ].sum( )
        .reset_index( name = 'proportion_number' )
    )

    # Output
    return aged_number_proportion , unaged_number_proportion , sex_number_proportions , age_number_proportions

def transect_stratum_weights( aged_proportions_df ,
                              unaged_proportions_df ,
                              sex_proportions_df ,
                              fitted_weight_df ):

    # Sum number proportions of aged specimens per stratum 
    aged_unaged_proportions = (
        aged_proportions_df[ aged_proportions_df.sex == 'all' ]
        .groupby( [ 'stratum_num' ] )[ 'proportion_number_overall_aged' ].sum( )
        .reset_index( name = 'number_proportion_aged' )
    )

    # Calculate unaged proportions per stratum
    aged_unaged_proportions[ 'number_proportion_unaged' ] = (
        1.0 - aged_unaged_proportions[ 'number_proportion_aged' ]
    )

    # Calculate the mixed aged and unaged number proportions
    # ---- Merge aged and unaged number proportions
    stratum_proportions_sexed = (
        aged_unaged_proportions
        .merge( sex_proportions_df , on = [ 'stratum_num' ] )
    ) 
    # ---- Calculate unaged number proportions per sex per stratum
    stratum_proportions_sexed[ 'proportion_unaged' ] = (
        stratum_proportions_sexed.number_proportion_unaged / (
            stratum_proportions_sexed.number_proportion_unaged
            + stratum_proportions_sexed.proportion_number_overall_aged
        )
    )
    # ---- Calculate aged number proportions per sex per stratum
    stratum_proportions_sexed[ 'proportion_aged' ] = (
        stratum_proportions_sexed[ 'proportion_number_overall_aged' ] / (
            stratum_proportions_sexed[ 'proportion_number_overall_aged' ]
            + stratum_proportions_sexed[ 'proportion_unaged' ]
        )
    )
    # ---- Reduce columns and ensure sex is in 'male/female/all'
    stratum_proportions_sexed = (
        stratum_proportions_sexed[ stratum_proportions_sexed.sex != 'unsexed' ]
        [ [ 'stratum_num' , 'sex' , 'proportion_aged' , 'proportion_unaged' ] ]
    )

    # Combine the aged-unaged (or station-specific) proportions for calculations
    # ---- Wide-to-long DataFrame
    station_proportions = pd.wide_to_long( stratum_proportions_sexed , 
                                           stubnames = "proportion" , 
                                           i = [ 'stratum_num' , 'sex' ] , 
                                           j = 'group' ,
                                           sep = "_" ,
                                           suffix = "\\w+" ).reset_index( )
    # ---- Convert to Table (to replicate indexed matrix operations)
    station_proportions_table = station_proportions.pivot_table( index = [ 'group' , 'sex' ]  , 
                                                                 columns = [ 'stratum_num' ] , 
                                                                 values = 'proportion' ).fillna( 0.0 )
    
    # Calculate the number length proportions that will be later converted into weight
    # ---- Aged length bins
    aged_length_distribution = (
        aged_proportions_df
        .groupby( [ 'stratum_num' , 'sex' , 'length_bin' ] , observed = False )
        [ 'proportion_number_aged' ].sum( )
        .reset_index( name = 'number_proportion' )
    )
    # ---- Unaged length bins
    unaged_length_distribution = (
        unaged_proportions_df[ unaged_proportions_df.sex != 'unsexed' ]
        [ [ 'stratum_num' , 'sex' , 'length_bin' , 'proportion_number_unaged' ] ]
        .rename( columns = { 'proportion_number_unaged': 'number_proportion' } )
    )
    # ---- Concatenate the two datasets
    length_number_proportions = pd.concat( [ aged_length_distribution.assign( group = 'aged' ) ,
                                             unaged_length_distribution.assign( group = 'unaged' ) ] )    
    # ---- Convert to Table (to replicate indexed matrix operations)
    length_proportions_table = (
        length_number_proportions
        .pivot_table( index = [ 'group' , 'sex' , 'length_bin' ] ,
                     columns = [ 'stratum_num' ] ,
                     values = 'number_proportion' ,
                     observed = False )
        .fillna( 0.0 )
    )
    
    # Convert the fitteed weights into a Table (to replicate index matrix operations)
    fitted_weight_table = fitted_weight_df.pivot_table( index = [ 'sex' , 'length_bin' ] ,
                                                        values = 'weight_modeled' ,
                                                        observed = False )
    
    # Calculate the average weights for male, female, and all fish within each stratum
    # ---- All
    weight_all = (
        fitted_weight_table.loc[ 'all' ][ 'weight_modeled' ].values
        .dot(
            length_proportions_table.loc[ 'aged' , 'all' ] 
            * station_proportions_table.loc[ 'aged' , 'all' ]
            + length_proportions_table.loc[ 'unaged' , 'all' ] 
            * station_proportions_table.loc[ 'unaged' , 'all' ]
        )
    ) 
    # ---- Male (Note: the modeled weight calculated for all fish is used instead of the 
    # ---- male-specific values)
    weight_male = (
        fitted_weight_table.loc[ 'all' ][ 'weight_modeled' ].values
        .dot(
            length_proportions_table.loc[ 'aged' , 'male' ] 
            * station_proportions_table.loc[ 'aged' , 'male' ]
            + length_proportions_table.loc[ 'unaged' , 'male' ] 
            * station_proportions_table.loc[ 'unaged' , 'male' ]
        )
    ) 
    # ---- Female (Note: the modeled weight calculated for all fish is used instead of the 
    # ---- female-specific values)
    weight_female = (
        fitted_weight_table.loc[ 'all' ][ 'weight_modeled' ].values
        .dot(
            length_proportions_table.loc[ 'aged' , 'female' ] 
            * station_proportions_table.loc[ 'aged' , 'female' ]
            + length_proportions_table.loc[ 'unaged' , 'female' ] 
            * station_proportions_table.loc[ 'unaged' , 'female' ]
        )
    ) 
    # ---- Combine the stratum-averaged weights for each sex and all fish
    stratum_weight_df = pd.DataFrame(
        {
            'stratum_num': np.tile( np.unique( station_proportions.stratum_num ) , 
                                    len( np.unique( station_proportions.sex ) ) ) ,
            'sex': np.repeat( [ 'all' , 'male' , 'female' ] , 
                              len( np.unique( station_proportions.stratum_num ) ) ) ,
            'average_weight': np.concatenate( [ weight_all , weight_male , weight_female ] ) ,
        } ,
    )

    # Return output
    return stratum_weight_df

def transect_weight_proportions( aged_weights_df , unaged_weights_df , unaged_number_proportion_df , length_weight_df ) :

    # Sum the weights for each stratum from the aged and unaged biological datasets
    # ---- Aged specimens
    aged_stratum_weights = (
        aged_weights_df
        .count_variable( contrasts = [ 'stratum_num' ] ,
                         variable = 'weight' ,
                         fun = 'sum' )
        .rename( columns = { 'count': 'stratum_weight' } )        
    )
    # ---- Unaged catches
    unaged_stratum_weights = (
        unaged_weights_df
        .count_variable( contrasts = [ 'stratum_num' ] ,
                         variable = 'haul_weight' ,
                         fun = 'sum' )
        .rename( columns = { 'count': 'stratum_weight' } )
    )
    # ---- Sum together to get grand total weight per stratum
    total_stratum_weights = (
        aged_stratum_weights.set_index( [ 'stratum_num' ] )
        + unaged_stratum_weights.set_index( [ 'stratum_num' ] )
    ).reset_index( ).rename( columns = { 'stratum_weight': 'overall_stratum_weight' } )

    # Remove specimen data with key missing information
    aged_weights_filtered = (
        aged_weights_df[ aged_weights_df.sex != 'unsexed' ]
        .dropna( subset = [ 'length' , 'weight' , 'age' ] )
    ) 

    # Sum specimen weights across multiple indices/scales
    # ---- Weights over each sex, stratum, length, and age
    aged_weight_bins = (
        aged_weights_filtered 
        .count_variable( contrasts = [ 'stratum_num' , 'sex' , 'length_bin' , 'age_bin' ] , 
                         variable = 'weight' , 
                         fun = 'sum' )
        .rename( columns = { 'count': 'weight' } )
    )
    # ---- Weights of specifically sexed aged fish per sex and stratum
    aged_stratum_sex_weights = (
        aged_weight_bins[ aged_weight_bins.sex.isin( [ 'male' , 'female' ] ) ]
        .count_variable( contrasts = [ 'stratum_num' , 'sex' ] ,
                         variable = 'weight' ,
                         fun = 'sum' )
        .rename( columns = { 'count': 'sex_stratum_weight' } )
    )

    # Calculate weight proportions (aged and unaged samples)
    # ---- Merge binned weights with the summed weights for each sex among aged fish for each 
    # ---- stratum
    aged_weight_proportions = (
        aged_weight_bins
        .merge( aged_stratum_sex_weights , on = [ 'stratum_num' , 'sex' ] )
    )
    # ---- Calculate within-sex weight proportion across length and age bins for each stratum
    aged_weight_proportions[ 'weight_aged_sex_proportion' ] = (
        aged_weight_proportions[ 'weight' ]
        / aged_weight_proportions[ 'sex_stratum_weight' ]
    )
    # ---- Additional merge with summed weights among aged fish for each stratum
    aged_weight_proportions = aged_weight_proportions.merge( aged_stratum_weights ,
                                                             on = [ 'stratum_num' ] )
    # ---- Calculate the aged-specific weight proportions for each stratum 
    aged_weight_proportions[ 'weight_aged_proportion' ] = (
        aged_weight_proportions[ 'weight' ]
        / aged_weight_proportions[ 'stratum_weight' ]
    )
    # ---- Additional merge with summed weights among all fish collected during the survey for 
    # ---- each stratum
    aged_weight_proportions = aged_weight_proportions.merge( total_stratum_weights ,
                                                             on = [ 'stratum_num' ] )
    # ---- Calculate the aged-specific weight proportions relative to aged and unaged fish weights
    aged_weight_proportions[ 'weight_aged_overall_proportion' ] = (
        aged_weight_proportions[ 'weight' ]
        / aged_weight_proportions[ 'overall_stratum_weight' ]
    )
    # ---- Calculate the weight proportions for all aged sexed fish
    aged_weight_proportions_all = (
        aged_weight_proportions
        [ [ 'stratum_num' , 'length_bin' , 'age_bin' ,
            'weight_aged_proportion' , 'weight_aged_overall_proportion' ] ]
        .groupby( [ 'stratum_num' , 'length_bin' , 'age_bin' ] , observed = False )
        .sum( )
        .reset_index( )
        .rename( columns = { 'weight_aged_proportion': 'weight_aged_sex_proportion' } )
        .assign( sex = 'all' )
    )
    # ---- Concatenate the aged sex-specific and overall proportions into a single DataFrame
    aged_weight_proportions_overall = (
        pd.concat( [ aged_weight_proportions.filter( items = aged_weight_proportions_all.columns ) ,
                     aged_weight_proportions_all ] )
        [ [ 'stratum_num' , 'sex' , 'length_bin' , 'age_bin' , 'weight_aged_sex_proportion' ,
            'weight_aged_overall_proportion' ] ]
    )
    
    # Estimate the weight distributions of the unaged fish
    # ---- Merge the length-weight regression coefficients with the unaged number proportions
    unaged_weight_bins = (
        unaged_number_proportion_df[ unaged_number_proportion_df.sex.isin( [ 'male' , 'female' ] ) ]
        .merge( length_weight_df , 
                on = [ 'sex' , 'length_bin' ] , 
                how = 'outer' )
    )
    # ---- Sum the weight for each bin based on the count
    unaged_weight_bins[ 'weight' ] = (
       unaged_weight_bins[ 'weight_modeled' ] * unaged_weight_bins[ 'count' ] 
    ) 
    # ---- Append the summed weight for each length bin
    unaged_weight_bins[ 'weight_binned_all' ] = (
        unaged_weight_bins
        .groupby( [ 'stratum_num' , 'length_bin' ] , observed = False )[ 'weight' ]
        .transform( 'sum' )
    )
    # ---- Sum fitted weights for each sex
    unaged_sex_weights = (
        unaged_weight_bins.groupby( [ 'stratum_num' , 'sex' ] , observed = False )
        .apply( lambda x: x[ 'weight' ].sum( ) / x[ 'weight_binned_all' ].sum( ) , 
                include_groups = False )
        .reset_index( name = 'proportion_unaged_sex_weight' )
    )
    # ---- Merge with unaged stratum weights estimated from survey-wide catches
    unaged_sex_weights = (
        unaged_sex_weights
        .merge( unaged_stratum_weights , on = [ 'stratum_num' ] , how = 'outer' )
        .dropna( subset = [ 'sex' , 'proportion_unaged_sex_weight' ] )
    )
    # ---- Re-weight the summed weights for each sex
    unaged_sex_weights[ 'weight_unaged_sex' ] = (
        unaged_sex_weights[ 'proportion_unaged_sex_weight' ]
        * unaged_sex_weights[ 'stratum_weight' ]
    )
    # ---- Calculate the unaged sex weight proportions relative to all aged and unaged fish
    # -------- Merge with total stratum weights
    unaged_sex_weights = (
        unaged_sex_weights.merge( total_stratum_weights , on = [ 'stratum_num' ] )
    )
    # -------- Compute proportion
    unaged_sex_weights[ 'weight_unaged_overall_proportion' ] = (
        unaged_sex_weights[ 'weight_unaged_sex' ]
        / unaged_sex_weights[ 'overall_stratum_weight' ]
    )

    return aged_weight_proportions_overall , unaged_sex_weights

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
    return weight_strata

def calculate_aged_unaged_proportions( specimen_data: pd.DataFrame ,
                                       weight_strata: pd.DataFrame ):
    """
   Compute the overall weight proportions for aged and uanged fish within each stratum 

    Parameters
    ----------
    specimen_data: pd.DataFrame
        Dataframe containing specimen weights
    weight_strata: pd.DataFrame
        Dataframe contained summed weights of both aged and unaged fish
    """   

    ### Calculate adult:all age ratio for each stratum
    # ---- Drop unaged fish from Station 2 data where fish are expected to be aged
    specimen_data_filtered = specimen_data[ specimen_data.sex != 'unsexed' ].dropna( how = 'any' , subset = 'age' ) 

    # ---- Sum the weight of adult fish
    specimen_adult = specimen_data_filtered[ specimen_data_filtered.age > 1 ].groupby( [ 'stratum_num' ] )[ 'weight' ].sum( )

    # ---- Sum the weight of all fish
    specimen_all = specimen_data_filtered.groupby( [ 'stratum_num' ] )[ 'weight' ].sum( ).to_frame( )

    # ---- Calculate the ratio between the two groups
    specimen_all[ 'weight_ratio' ] = ( specimen_adult / specimen_all.weight ).fillna( 0 )

    ### Calculate the aged proportions of all fish
    # ---- Merge with `weight_strata`
    aged_proportions = specimen_all.reset_index( ).merge( weight_strata , on = [ 'stratum_num' ] )

    # ---- Calculate the aged proportion of all fish
    aged_proportions[ 'proportion_aged_weight_all' ] = (
        aged_proportions.weight / aged_proportions.weight_stratum_all
    )

    # ---- Calculate the aged proportion of adult fish
    aged_proportions[ 'proportion_aged_weight_adult' ] = (
        aged_proportions.proportion_aged_weight_all * aged_proportions.weight_ratio
    )

    ### Calculate the unaged weight proportions
    # ---- All fish
    aged_proportions[ 'proportion_unaged_weight_all' ] = (
        1.0 - aged_proportions.proportion_aged_weight_all
    )

    # ---- Adult fish
    aged_proportions[ 'proportion_unaged_weight_adult' ] = (
        1.0 - aged_proportions.proportion_aged_weight_adult
    )

    return aged_proportions.drop( [ 'weight_ratio' ] , axis = 1 )

def aged_weight_proportions( specimen_data: pd.DataFrame ,
                             length_intervals: np.ndarray ,
                             age_intervals: np.ndarray ):
    """
    Calculate length-age binned weight proportions based on `specimen_df` (Station 2)   

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
    # ---- Drop unaged fish
    # ==== !!! TODO: pending what FEAT says, weights associated with 
    # ==== missing ages should be added into the 'unaged' category. 
    # ==== This would further mean apportioning these into `weight_strata_aged_uanged`
    specimen_data_filtered = specimen_data[ specimen_data.sex != 'unsexed' ].dropna( how = 'any' , subset = 'age' )
    
    # ---- Bin length and age measurements
    specimen_data_filtered = (
        specimen_data_filtered
        # ---- Bin length
        .bin_variable( length_intervals , 'length' )
        # ---- Age bin
        .bin_variable( age_intervals , 'age' )
    )

    ### Sum weights within each length and age bin for each sex within each stratum
    # ---- Create separate 'weight_adult' column
    specimen_data_filtered[ 'weight_adult' ] = (
        np.where( specimen_data_filtered.age > 1 , specimen_data_filtered.weight , 0.0 )
    )
    
    # ---- Calculate aggregate sums of fish weights
    specimen_binned_weight = (
        specimen_data_filtered
        # ---- Group weight summations across stratum/species/sex/length/age
        .groupby( [ 'stratum_num' , 'species_id' , 'sex' , 'length_bin' , 'age_bin' ] ,
                  observed = False  )
        # ---- Sum the weights 
        .agg( weight_all =( 'weight' , 'sum' )  , weight_adult =( 'weight_adult' , 'sum' ) ) 
        # ---- Fill empty/non-existent values with 0's
        .fillna( 0 )
        .reset_index( )
    )
    
    ### Calculate the relative weight proportions of each length-age bin for each sex within each stratum
    # ---- Sum the total weights for age-1+ and age-2+ 
    proportions_weight_length_age_sex = (
        specimen_binned_weight
        # ---- Calculate total sex-specific weights for each stratum
        .assign( total_weight_sex_all = lambda df: df.groupby( [ 'stratum_num' , 'species_id' , 'sex' ] )[ 'weight_all' ].transform( 'sum' ) ,
                 total_weight_sex_adult = lambda df: df.groupby( [ 'stratum_num' , 'species_id' , 'sex' ] )[ 'weight_adult' ].transform( 'sum' ) )
    )

    # ---- Sum the weight proportions for age-1+
    proportions_weight_length_age_sex[ 'proportion_weight_sex_all' ] = (
        proportions_weight_length_age_sex.weight_all / proportions_weight_length_age_sex.total_weight_sex_all
    ).fillna( 0 )

     # ---- Sum the weight proportions for age-2+
    proportions_weight_length_age_sex[ 'proportion_weight_sex_adult' ] = (
        proportions_weight_length_age_sex.weight_adult / proportions_weight_length_age_sex.total_weight_sex_adult
    ).fillna( 0 )

    ### Return output
    return proportions_weight_length_age_sex

def aged_sex_weight_proportions( proportions_weight_length_age_sex: pd.DataFrame ,
                                 aged_proportions: pd.DataFrame ):
    """
    Compute the weight proportion for aged/unaged fish of either sex (2x2 choice) across all fish.
    Also compute the weight proportion of only adult aged/unaged fish of either sex across all fish.

    Parameters
    ----------
    proportions_weight_length_age_sex: pd.DataFrame
        Dataframe containing sexed weight proportions distributed across
        each age and length bin
    aged_proportions: pd.DataFrame
        Dataframe containing the weight proportions of aged fish within each stratum
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
    aged_sex_proportions = aged_sex_weights.merge( aged_proportions ,
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
    aged_sex_proportions = aged_sex_proportions.fillna( 0 ).filter( regex = "^(?=proportion_aged_|proportion_weight|stratum_num|sex).*" )

    ### Return output 
    return aged_sex_proportions.reset_index( drop = True )

def distribute_aged_weight_proportions( proportions_weight_length_age_sex: pd.DataFrame ,
                                        aged_sex_proportions: pd.DataFrame ):
    """
    Distribute overall weight proportions across each sex and age/length bins   

    Parameters
    ----------
    proportions_weight_length_age_sex: pd.DataFrame
        Dataframe containing sexed weight proportions distributed across
        each age and length bin for age-1+ (`proportion_weight_sex_all`) and 
        age-2+ (`proportion_weight_sex_adult`) fish
    aged_sex_proportions: pd.DataFrame
        Dataframe contained weight proportions of sexed fish for aged fish
    """     

    ### Calculate the normalized age-length-sex aged fish weight proportions
    # ---- Merge `proportions_weight_length_age_sex` with `aged_sex_proportions`
    distributed_aged_weight_proportions = proportions_weight_length_age_sex.merge( aged_sex_proportions ,
                                                                                   on = [ 'stratum_num' , 'sex' ] , 
                                                                                   how = 'left' )
    
    # ---- Normalized weight proportions for age-1+ fish
    distributed_aged_weight_proportions[ 'normalized_proportion_weight_all' ] = (
        distributed_aged_weight_proportions.proportion_weight_sex_all * distributed_aged_weight_proportions.proportion_weight_all
    )

    # ---- Normalized weight proportions for age-2+ fish
    distributed_aged_weight_proportions[ 'normalized_proportion_weight_adult' ] = (
        distributed_aged_weight_proportions.proportion_weight_sex_adult * distributed_aged_weight_proportions.proportion_weight_adult
    )   

    # ---- Remove unnecessary columns and return the dataframe
    return distributed_aged_weight_proportions.filter( regex = '^(?!weight_|total_weight_|proportion_).*' )

def calculate_aged_biomass( kriging_biomass_df: pd.DataFrame ,
                            specimen_data: pd.DataFrame ,
                            length_distribution: pd.DataFrame ,
                            age_distribution: pd.DataFrame ,
                            aged_proportions: pd.DataFrame ):
    """
    Calculate the kriged biomass distributed over length and age for each sex and among
    all fish (sexed)

    Parameters
    ----------
    kriging_biomass_df: pd.DataFrame
        Dataframe containing kriged biomass estimates
    specimen_data: pd.DataFrame
        Dataframe containing specimen weights
    length_distribution: np.ndarray
        Array containing length bins/intervals
    age_distribution: np.ndarray
        Array containing age bins/intervals
    aged_proportions: pd.DataFrame
        Dataframe containing the weight proportions of aged fish within each stratum
    """

    ### Sum aged fish weights across age and length bins for each sex, then calculate weight proportions within each sex 
    proportions_weight_length_age_sex = aged_weight_proportions( specimen_data ,
                                                                 length_distribution ,
                                                                 age_distribution )
    
    ### Calculate the weight proportion of aged/unaged fish of either sex (2x2 choice) across all fish 
    # ---- belonging to each sex
    aged_sex_proportions = aged_sex_weight_proportions( proportions_weight_length_age_sex ,
                                                        aged_proportions )
    
    ### Calculate weight proportions of aged fish distributed over age-length bins
    # ---- for each sex within each stratum relative to the summed weights of each
    # ---- stratum (i.e. aged + unaged weights)    
    distributed_sex_length_age_proportions = distribute_aged_weight_proportions( proportions_weight_length_age_sex ,
                                                                                 aged_sex_proportions )
    
    ### Sum 'kriging_biomass_df' across each stratum for appropriate stratum-specific apportionment
    kriged_stratum_biomass = kriging_biomass_df.groupby( [ 'stratum_num' ] , observed = False )[ 'B_adult_kriged' ].sum( )

    ### Apportion sexed biomass across age-length
    # ---- Initialize dataframe
    stratum_sexed_kriged_biomass = distributed_sex_length_age_proportions.set_index( 'stratum_num' )
    stratum_sexed_kriged_biomass[ 'biomass_kriged' ] = kriged_stratum_biomass

    # ---- Apportioned biomass for age-1+ fish to different sex-length-age bins
    stratum_sexed_kriged_biomass[ 'biomass_sexed_aged_all' ] = (  
        stratum_sexed_kriged_biomass.biomass_kriged
        * stratum_sexed_kriged_biomass.normalized_proportion_weight_all
    )

    # ---- Apportioned biomass for age-2+ fish to different sex-length-age bins
    stratum_sexed_kriged_biomass[ 'biomass_sexed_aged_adult' ] = (  
        stratum_sexed_kriged_biomass.biomass_kriged
        * stratum_sexed_kriged_biomass.normalized_proportion_weight_adult
    )

    # ---- Sum across strata to produce 'grand totals' for each sex
    apportioned_sexed_kriged_biomass = (
        stratum_sexed_kriged_biomass
        .groupby( [ 'species_id' , 'sex' , 'length_bin' , 'age_bin' ] ,
                  observed = False )
        .agg( { 'biomass_sexed_aged_all': 'sum' ,
                'biomass_sexed_aged_adult': 'sum' } )
        .reset_index( )
    )

    # ---- Sum across strata to produce 'grand totals' across all fish
    apportioned_total_kriged_biomass = (
        apportioned_sexed_kriged_biomass
        .groupby( [ 'species_id' , 'length_bin' , 'age_bin' ] ,
                  observed = False )
        .agg( biomass_aged_all = ( 'biomass_sexed_aged_all' , 'sum' ) ,
              biomass_aged_adult = ( 'biomass_sexed_aged_adult' , 'sum' ) )
        .reset_index( )
    )

    ### Return output (tuple)
    return apportioned_sexed_kriged_biomass , apportioned_total_kriged_biomass

def unaged_number_proportions( length_data: pd.DataFrame ,
                               length_intervals: np.ndarray ):
    """
    Calculate length binned number proportions among unaged fish  

    Parameters
    ----------
    length_data: pd.DataFrame
        Dataframe containing length data measured from unaged fish
    length_intervals: np.ndarray
        Array containing length bins/intervals
    """ 

    ### Calculate number proportion
    # ---- Bin length measurements
    length_data_binned = (
        length_data
        # ---- Bin length
        .bin_variable( length_intervals , 'length' ) 
    )

    # ---- Sum the number of individuals within each bin 
    proportions_unaged_length = (
        length_data_binned
        # ---- Group number summations across stratum/species/sex/length
        .groupby( [ 'stratum_num' , 'species_id' , 'length_bin' ] ,
                  observed = False )
        # ---- Sum count
        .agg( number_all = ( 'length_count' , 'sum' ) )
        # ---- Fill empty/non-existent values with 0's
        .fillna( 0 )
        .reset_index( )
    )

    ### Calculate the number proportions
    # --- Stratum total counts
    proportions_unaged_length[ 'stratum_number_all' ] = (
        proportions_unaged_length.groupby( [ 'stratum_num' ] )[ 'number_all' ].transform( 'sum' )
    ).fillna( 0 )

    # ---- Proportions of each sex-length bin pair relative to `stratum_number_all`
    proportions_unaged_length[ 'proportion_number_all' ] = (
        proportions_unaged_length.number_all / proportions_unaged_length.stratum_number_all
    ).fillna( 0 )

    ### Filter out unnecessary columns and return output
    return proportions_unaged_length.filter( regex = '^(?!number_|stratum_number_).*' )

def unaged_weight_proportions( proportions_unaged_length: pd.DataFrame ,
                               length_weight_df: pd.DataFrame ):
    """
    Compute the weight proportion of each length bin within each stratum for all unaged fish samples. 

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
    proportions_unaged_weight_length = proportions_unaged_length.merge( length_weight_all ,
                                                                        on = [ 'length_bin' ] ,
                                                                        how = 'left' )
    
    # ---- Calculate the estimated weight for each length bin within each stratum (`w_ln_all_array` in the original Matlab code)  
    proportions_unaged_weight_length[ 'weight_proportion_all' ] = (  
        proportions_unaged_weight_length.weight_modeled * proportions_unaged_weight_length.proportion_number_all  
    )  

    ### Calculate weight proportions for each length bin for all strata  
    # ---- Calculate the summed estimated weights across all length bins within each stratum (`w_ln_array_sum` in the original Matlab code)  
    proportions_unaged_weight_length[ 'weight_stratum_proportion' ] = (  
       proportions_unaged_weight_length.groupby( [ 'stratum_num' ] )[ 'weight_proportion_all' ].transform( 'sum' )  
    )  
    
    # ---- Calculate the weight proportion of each length bin within each stratum (`w_ln_all_N` in the original Matlab code)  
    proportions_unaged_weight_length[ 'proportion_weight_length' ] = (  
        proportions_unaged_weight_length.weight_proportion_all / proportions_unaged_weight_length.weight_stratum_proportion  
    )  

    ### Drop unnecessary columns and return output
    return proportions_unaged_weight_length.filter( regex = '^(?!weight_|proportion_number_).*' )

def unaged_sex_weight_proportions( length_data: pd.DataFrame ,
                                   length_intervals: pd.DataFrame ,
                                   regression_parameters: pd.DataFrame ,
                                   aged_proportions: pd.DataFrame ):
    """
    Calculate the weight proportion of different sexes within each stratum for all unaged fish samples.

    Parameters
    ----------
    length_data: pd.DataFrame
        Dataframe containing length data measured from unaged fish
    length_intervals: np.ndarray
        Array containing length bins/intervals
    regression_parameters: pd.DataFrame
        Dataframe containing the best-fit coefficients for the length-weight log-linear relationship for
        each sex
    aged_proportions: pd.DataFrame
        Dataframe containing the weight proportions of unaged fish within each stratum
    """ 

    ### Prepare data
    # ---- Bin length measurements
    length_data_binned = (
        length_data
        # ---- Bin length
        .bin_variable( length_intervals , 'length' ) 
    )

    # ---- Extract bin length values (i.e. mid)
    length_data_binned[ 'length_bin_value' ] = (
        length_data_binned[ 'length_bin' ].apply( lambda x: x.mid )
    )

    # ---- Merge `length_data_filtered ` and `length_weight_sex`
    sexed_unaged_model_weights = length_data_binned.merge( regression_parameters ,
                                                           on = [ 'sex' ] ,
                                                           how = 'left' )
    
    ### Interpolate weights over length values for each sex (weighted by `length_count`)
    # ---- Convert length to weight using the fit regression parameters
    sexed_unaged_model_weights[ 'weight_modeled' ] = (
        10 ** sexed_unaged_model_weights.initial
        * sexed_unaged_model_weights.length ** sexed_unaged_model_weights.rate
    )

    # ---- Group by each stratum and sex to calculate the average weight per length bin
    fit_strata_weights = (
        sexed_unaged_model_weights
        .groupby( [ 'stratum_num' , 'sex' ] ,
                  observed = False )
        .apply( lambda df: ( df.weight_modeled * df.length_count ).sum( ) ,
                include_groups = False )
        .reset_index( name = 'weight_modeled' )
    )

    ### Calculate the unaged weight proportion relative to the haul catches
    # ---- Merge the unaged weight strata with the interpolated weights
    proportions_unaged_weight_sex = fit_strata_weights.merge( aged_proportions , 
                                                              on = [ 'stratum_num' ] ,
                                                              how = 'left' )
    
    # ---- Sum `weight_interp` over each stratum
    proportions_unaged_weight_sex[ 'stratum_sex_weight_modeled' ] = (
        proportions_unaged_weight_sex.groupby( ['stratum_num' ] )[ 'weight_modeled' ].transform( 'sum' ) 
    )

    # ---- Calculate the normalized weights for each sex
    proportions_unaged_weight_sex[ 'stratum_sex_weight_normalized' ] = (
        proportions_unaged_weight_sex[ 'weight' ] *
        proportions_unaged_weight_sex.weight_modeled / proportions_unaged_weight_sex.stratum_sex_weight_modeled
    )

    # ---- Sum the normalized stratum sex weight across each stratum
    proportions_unaged_weight_sex[ 'proportions_weight_sex' ] = (
        proportions_unaged_weight_sex.stratum_sex_weight_normalized /
        proportions_unaged_weight_sex.weight_stratum_all
    )

    # ---- Sum 'proportions_weight_sex'
    proportions_unaged_weight_sex[ 'proportions_weight_sex_total' ] = (
        proportions_unaged_weight_sex.groupby( [ 'stratum_num' ] )[ 'proportions_weight_sex' ].transform( 'sum' )
    )

    # ---- Calculate the final sex proportion
    proportions_unaged_weight_sex[ 'proportion_weight_sex' ] = (
        proportions_unaged_weight_sex.proportions_weight_sex / proportions_unaged_weight_sex.proportions_weight_sex_total
    )

    ### Remove unnecessary columns and return output
    return proportions_unaged_weight_sex[ [ 'stratum_num' , 'sex' , 'proportion_weight_sex' ] ]

def calculate_unaged_biomass( kriging_biomass_df: pd.DataFrame ,
                              length_data: pd.DataFrame ,
                              length_distribution: np.ndarray ,
                              length_weight_df: pd.DataFrame ,
                              regression_parameters: pd.DataFrame ,
                              aged_proportions: pd.DataFrame ):
    """
    Calculate the kriged biomass distributed over length and age for each sex and among
    all fish (sexed)

    Parameters
    ----------
    kriging_biomass_df: pd.DataFrame
        Dataframe containing kriged biomass estimates
    length_data: np.ndarray
        Dataframe containing length data measured from unaged fish
    length_distribution: pd.DataFrame
        Dataframe containing the weight proportions for each length and age bin computed for all fish
    length_weight_df: pd.DataFrame
        Dataframe containing the modeled weights for each sex and length bin
        fit from the length-weight regression (fit for all animals)
    regression_parameters: pd.DataFrame
        Dataframe containing the best-fit coefficients for the length-weight log-linear relationship for
        each sex
    aged_proportions: pd.DataFrame
        Dataframe containing the weight proportions of unaged fish within each stratum
    """

    ### Calculate number proportion
    proportions_unaged_length = unaged_number_proportions( length_data , length_distribution )

    ### Calculate the weight proportion of unaged fish of each length bin (W_Ln_ALL in the original Matlab code)  
    proportions_unaged_weight_length = unaged_weight_proportions( proportions_unaged_length ,
                                                                  length_weight_df )

    ### Calculate sex-specific weight proportions for unaged fish (within unaged fish)
    proportions_unaged_weight_sex = unaged_sex_weight_proportions( length_data ,
                                                                   length_distribution ,
                                                                   regression_parameters ,
                                                                   aged_proportions )

    ### Sum 'kriging_biomass_df' across each stratum for appropriate stratum-specific apportionment
    kriged_stratum_biomass = kriging_biomass_df.groupby( [ 'stratum_num' ] , observed = False )[ 'B_adult_kriged' ].sum( )

    ### Calculate sex-specific biomass apportionment
    # ---- Index `unaged_proportions`
    indexed_unaged_proportions = aged_proportions.filter( regex = "^(?=.*stratum_num|.*unaged).*(?!.*aged)" ).set_index( 'stratum_num' )

    # ---- Index `proportions_unaged_weight_length`
    indexed_proportions_unaged_weight_length = proportions_unaged_weight_length.set_index( 'stratum_num' )

    # ---- Initialize dataframe
    stratum_sexed_kriged_biomass = proportions_unaged_weight_sex.merge( indexed_proportions_unaged_weight_length ,
                                                                        on = [ 'stratum_num' ] ,
                                                                        how = 'left' ).set_index( 'stratum_num' )
    stratum_sexed_kriged_biomass[ 'proportion_unaged_weight_all' ] = indexed_unaged_proportions.proportion_unaged_weight_all
    stratum_sexed_kriged_biomass[ 'proportion_unaged_weight_adult' ] = indexed_unaged_proportions.proportion_unaged_weight_adult
    stratum_sexed_kriged_biomass[ 'biomass_kriged' ] = kriged_stratum_biomass

    # ---- Calculate the overall sexed proportions taking into account aged and unaged fish (all age class)
    stratum_sexed_kriged_biomass[ 'proportion_sexed_weight_all' ] = (
        stratum_sexed_kriged_biomass.proportion_weight_sex *
        stratum_sexed_kriged_biomass.proportion_unaged_weight_all
    )

    # ---- Calculate the overall sexed proportions taking into account aged and unaged fish (adults)
    stratum_sexed_kriged_biomass[ 'proportion_sexed_weight_adult' ] = (
        stratum_sexed_kriged_biomass.proportion_weight_sex *
        stratum_sexed_kriged_biomass.proportion_unaged_weight_adult
    )

    # ---- Apportioned biomass for age-1+ fish
    stratum_sexed_kriged_biomass[ 'biomass_sexed_unaged_all' ] = (
        stratum_sexed_kriged_biomass.biomass_kriged
        * stratum_sexed_kriged_biomass.proportion_sexed_weight_all
        * stratum_sexed_kriged_biomass.proportion_weight_length
    )

    # ---- Apportioned biomass for age-2+ fish
    stratum_sexed_kriged_biomass[ 'biomass_sexed_unaged_adult' ] = (
        stratum_sexed_kriged_biomass.biomass_kriged
        * stratum_sexed_kriged_biomass.proportion_sexed_weight_adult
        * stratum_sexed_kriged_biomass.proportion_weight_length
    )

    # ---- 'Grand total' for age-1+ and age-2+ sexed unaged fish
    apportioned_sexed_kriged_biomass = (
            stratum_sexed_kriged_biomass
            .groupby( [ 'species_id' , 'sex' , 'length_bin'  ] ,
                      observed = False )
            .agg( { 'biomass_sexed_unaged_all': 'sum' ,
                    'biomass_sexed_unaged_adult': 'sum' } )
            .reset_index( )
        )
    
    ### Calculate biomass proportions across length
    # ---- Initialize dataframe
    stratum_kriged_biomass = indexed_proportions_unaged_weight_length
    stratum_kriged_biomass[ 'proportion_unaged_weight_all' ] = indexed_unaged_proportions.proportion_unaged_weight_all
    stratum_kriged_biomass[ 'proportion_unaged_weight_adult' ] = indexed_unaged_proportions.proportion_unaged_weight_adult
    stratum_kriged_biomass[ 'biomass_kriged' ] = kriged_stratum_biomass

    # ---- Apportioned biomass for all fish
    stratum_kriged_biomass[ 'biomass_unaged_all' ] = (
        stratum_kriged_biomass.biomass_kriged
        * stratum_kriged_biomass. proportion_unaged_weight_all
        * stratum_kriged_biomass.proportion_weight_length
    )

    # ---- Apportioned biomass for adult fish
    stratum_kriged_biomass[ 'biomass_unaged_adult' ] = (
        stratum_kriged_biomass.biomass_kriged
        * stratum_kriged_biomass.proportion_unaged_weight_adult
        * stratum_kriged_biomass.proportion_weight_length
    )

    ### Return output (tuple)
    return apportioned_sexed_kriged_biomass
    
def apply_age_bins( aged_sex_biomass: pd.DataFrame ,
                    unaged_sex_biomass: pd.DataFrame ):
    """
    Redistribute unaged biomass over the defined age distribution

    Parameters
    ----------
    aged_sex_biomass: pd.DataFrame
        Dataframe containing biomass data of sexed aged fish
    unaged_sex_biomass: pd.DataFrame
        Dataframe containing biomass data of sexed unaged fish
    """ 

    ### Calculate summed length-bin biomass (i.e. sum across age bins)
    # ---- All fish
    aged_sex_biomass[ 'summed_aged_biomass_all' ] = (
        aged_sex_biomass.groupby( [ 'sex' , 'length_bin' ] , observed = False )[ 'biomass_sexed_aged_all' ].transform( 'sum' )
    )

    # ---- Adult fish
    aged_sex_biomass[ 'summed_aged_biomass_adult' ] = (
        aged_sex_biomass.groupby( [ 'sex' , 'length_bin' ] , observed = False )[ 'biomass_sexed_aged_adult' ].transform( 'sum' )
    )
    
    ### Merge unaged biomass length bins with aged bioamss length-age bins
    aged_unaged_sexed_biomass = aged_sex_biomass.merge( unaged_sex_biomass ,
                                                        on = [ 'length_bin' , 'sex' , 'species_id' ] ,
                                                        how = 'left' )

    ### Redistribute unaged biomass over the length-age bins for each sex
    # ---- All fish
    aged_unaged_sexed_biomass[ 'biomass_sexed_unaged_all' ] = (
        aged_unaged_sexed_biomass.biomass_sexed_unaged_all * 
        aged_unaged_sexed_biomass.biomass_sexed_aged_all /
        aged_unaged_sexed_biomass.summed_aged_biomass_all
    ).fillna( 0 )

    # ---- Adult fish
    aged_unaged_sexed_biomass[ 'biomass_sexed_unaged_adult' ] = (
        aged_unaged_sexed_biomass.biomass_sexed_unaged_adult * 
        aged_unaged_sexed_biomass.biomass_sexed_aged_adult /
        aged_unaged_sexed_biomass.summed_aged_biomass_adult
    ).fillna( 0 )

    ### Remove unnecessary columns 
    aged_unaged_sexed_biomass.drop( [ 'biomass_sexed_aged_all' , 
                                      'biomass_sexed_aged_adult' ,
                                      'summed_aged_biomass_all' ,
                                      'summed_aged_biomass_adult' ] ,
                                    axis = 1 ,
                                    inplace = True )
    
    ### Sum sexes together to retrieve the total biomass 
    aged_unaged_biomass = (
        aged_unaged_sexed_biomass
        .groupby( [ 'length_bin' , 'age_bin' , 'species_id' ] ,
                  observed = False )
        .agg( biomass_unaged_all = ( 'biomass_sexed_unaged_all' , 'sum' ) ,
              biomass_unaged_adult = ( 'biomass_sexed_unaged_adult' , 'sum' ) )
        .reset_index( )
    )    

    ### Return output (tuple)
    return aged_unaged_sexed_biomass , aged_unaged_biomass