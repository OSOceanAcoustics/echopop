import numpy as np
import pandas as pd
from typing import Union, List
from .spatial.transect import correct_transect_intervals
from .utils.operations import group_interpolator_creator

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

def fit_length_weight_relationship( specimen_data: pd.DataFrame , length_distribution: pd.DataFrame ) :
    """
    Fit a length-weight relationship across discrete bins

    Parameters
    ----------
    species_id : np.float64
        Numeric code representing a particular species of interest

    Notes
    -----
    This function first fits a length-weight regression based on measured 
    values and then produces an array of fitted weight values based on 
    binned length values.  
    The length-weight relationship produced here are used for later 
    biomass calculations and apportionment.
    """    

    # Gather specimen measurements to represent 'all' fish
    specimen_data_all = specimen_data.assign( sex = 'all' )

    # Combine sexed and 'all' specimens
    # ---- Vertical concatenation
    specimen_data_all = pd.concat( [ specimen_data[ specimen_data.group_sex == 'sexed' ] ,
                                     specimen_data_all ] )
    # ---- Remove bad values
    specimen_data_all.dropna( subset = [ 'length' , 'weight' ] , inplace = True )

    # Fit length-weight linear regression by male, female, and all fish
    length_weight_regression_df = (
        specimen_data_all
        .groupby( [ 'species_id' , 'sex' ] )
        .apply( lambda df: pd.Series( np.polyfit( np.log10( df[ 'length' ] ) ,
                                                  np.log10( df[ 'weight' ] ) , 1 ) ,
                                     index = [ 'rate' , 'initial' ] ) ,
                include_groups = False )
        .reset_index( )
    )

    # Predict weights for binned lengths
    # ---- Initialize dataframe
    weight_fitted_df = length_distribution.copy( )
    # ---- Expand/merge with length-weight regression coefficients
    weight_fitted_df = weight_fitted_df.merge( length_weight_regression_df , how = 'cross' )
    # ---- Predict weight per bin
    weight_fitted_df[ 'weight_modeled' ] = (
        10.0 ** weight_fitted_df[ 'initial' ] 
        * weight_fitted_df[ 'length_bins' ] ** weight_fitted_df[ 'rate' ]
    )
    # ---- Drop unused columns
    weight_fitted_df = weight_fitted_df.filter( [ 'length_intervals' , 'species_id' , 'sex' , 
                                                  'weight_modeled' ] )
    
    # Adjust for cases where there are too few (< 5) specimens within a given length bin
    # ---- Count number of specimens across length bins
    weight_fitted_distribution_df = (
        specimen_data_all
        .count_variable( contrasts = [ 'species_id' , 'sex' , 'length_bin' ] , 
                         variable = 'length' , fun = 'size' )
        .set_index( [ 'species_id' , 'sex' , 'length_bin' ] )
    )
    # ---- Get mean weight per bin as well
    weight_fitted_distribution_df[ 'weight_mean' ] = (
        specimen_data_all
        .groupby( [ 'species_id' , 'sex' , 'length_bin' ] , observed = False )[ 'weight' ]
        .mean( ).fillna( 0.0 )
    )
    # ---- Merge with the fitted weights
    weight_fitted_distribution_df = (
        weight_fitted_distribution_df.merge( weight_fitted_df ,
                                             left_on = [ 'species_id' , 'sex' , 'length_bin' ] ,
                                             right_on = [ 'species_id' , 'sex' , 'length_intervals' ] )
    )
    # ---- Find fitted weights accounting for low sample sizes
    weight_fitted_distribution_df[ 'weight_fitted' ] = np.where( weight_fitted_distribution_df[ 'count' ] < 5 ,
                                                                 weight_fitted_distribution_df[ 'weight_modeled' ] ,
                                                                 weight_fitted_distribution_df[ 'weight_mean' ] )
    # ---- Pull out unused columns
    weight_fitted_distribution_df = weight_fitted_distribution_df.filter( [ 'species_id' , 'sex' ,
                                                                            'length_intervals' ,
                                                                            'weight_fitted' ] )
    # Return output
    return (
        {
            'length_weight_regression': {
                'parameters_df': length_weight_regression_df ,
                'weight_fitted_df': weight_fitted_distribution_df ,
            } 
        } 
    )
    
def quantize_number_counts( specimen_data , length_data , stratum ) :

    # Sub-select the correct stratum definition
    if stratum.lower( ) == 'ks': 
        stratum_col = 'stratum_num'
    else:
        stratum_col = 'stratum_inpfc'

    # Bin counts by sex
    # ---- Aged
    aged_number_distribution = (
        pd.concat( [ specimen_data , specimen_data.assign( sex = 'all' ) ] )
        .count_variable( contrasts = [ stratum_col , 'species_id' , 'sex' , 'length_bin' , 'age_bin' ] , 
                         variable = 'length' , 
                         fun = 'size' )
    )
    # ---- Filter out unsexed data for parallel number counts and drop any NA's
    aged_number_distribution_filtered = (
        pd.concat( [ specimen_data[ specimen_data.sex != 'unsexed' ] ,
                     specimen_data[ specimen_data.sex != 'unsexed' ].assign( sex = 'all' ) ] )
        .dropna( subset = [ 'length' , 'age' , 'weight' ] )
        .count_variable( contrasts = [ stratum_col , 'species_id' , 'sex' , 'length_bin' , 'age_bin' ] , 
                         variable = 'length' , 
                         fun = 'size' )        
    ) 
    # ---- Unaged
    unaged_number_distribution = (
        pd.concat( [ length_data ,
                     length_data.assign( sex = 'all' ) ] )
        .count_variable( contrasts = [ stratum_col , 'species_id' , 'sex' , 'length_bin' ] , 
                         variable = 'length_count' , 
                         fun = 'sum' )
    )

    # Return output
    return (
        {
            'binned_aged_counts_df': aged_number_distribution ,
            'binned_aged_counts_filtered_df': aged_number_distribution_filtered ,
            'binned_unaged_counts_df': unaged_number_distribution
        }
    )

def number_proportions( count_dict ) :

    # Get the name of the stratum column
    stratum_col = [ col for col in count_dict[ 'binned_aged_counts_df' ].columns if 'stratum' in col.lower( ) ][ 0 ]

    # Calculate total numbers among aged samples
    # ---- Collect unique strata
    strata_unique = pd.DataFrame(
        {
            'stratum_num': np.unique( np.concatenate( [ count_dict[ 'binned_aged_counts_df' ][ stratum_col ] ,
                                                        count_dict[ 'binned_unaged_counts_df' ][ stratum_col ] ] ) ) ,
        } ,
    )  
    # ---- Collect unique sex
    sex_unique = pd.DataFrame(
        {
            'sex': np.unique( np.concatenate( [ count_dict[ 'binned_aged_counts_df' ].sex ,
                                                count_dict[ 'binned_unaged_counts_df' ].sex ] ) ) ,
        } ,
    )
    # ---- Collect unique species id
    species_unique = pd.DataFrame(
        {
            'species_id': np.unique( np.concatenate( [ count_dict[ 'binned_aged_counts_df' ].species_id ,
                                                       count_dict[ 'binned_unaged_counts_df' ].species_id ] ) ) ,
        } ,
    )
    # ---- Initialize dataframe
    count_total = (
        strata_unique
        .merge( sex_unique , how = 'cross' )
        .merge( species_unique , how = 'cross' )
        .set_index( [ 'stratum_num' , 'species_id' , 'sex' ] )
    )
    # ---- Aged (overall)
    count_total[ 'total_overall_aged' ] = (
        count_dict[ 'binned_aged_counts_df' ]
        .groupby( [ stratum_col , 'species_id' , 'sex' ] )[ 'count' ].sum( )
    )
    # ---- Aged (filtered)
    count_total[ 'total_filtered_aged' ] = (
        count_dict[ 'binned_aged_counts_filtered_df' ]
        .groupby( [ stratum_col , 'species_id' , 'sex' ] )[ 'count' ].sum( )
    )
    # ---- Unaged
    count_total[ 'total_overall_unaged' ] = (
        count_dict[ 'binned_unaged_counts_df' ]
        .groupby( [ stratum_col , 'species_id' , 'sex' ] )[ 'count' ].sum( )
    )
    # ---- Fill NaN
    count_total.fillna( 0 , inplace = True )
    count_total = count_total.reset_index( ).set_index( stratum_col )
    # ---- Grand totals
    count_total[ 'total_overall' ] = (
        count_total[ count_total.sex == 'all' ].total_filtered_aged
         + count_total[ count_total.sex == 'all' ].total_overall_unaged
    )
    count_total = count_total.reset_index( )

    # Calculate number proportions
    # ---- Aged (number distributed over age and length)
    aged_number_proportion = (
        count_dict[ 'binned_aged_counts_filtered_df' ]
        [ count_dict[ 'binned_aged_counts_filtered_df' ].sex.isin( [ 'male' , 'female' , 'all' ] ) ]        
        .merge( count_total[ [ stratum_col , 'species_id' , 'sex' , 'total_filtered_aged' , 'total_overall' ] ] , 
                on = [ stratum_col , 'species_id' , 'sex' ] )
    )
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
        .groupby( [ stratum_col , 'species_id' , 'sex' ] , observed = False )[ 'proportion_number_overall_aged' ].sum( )
        .reset_index( )
    )
    # ---- Unaged (number distributed over length)
    unaged_number_proportion = (
        count_dict[ 'binned_unaged_counts_df' ]
        .merge( count_total[ [ stratum_col , 'species_id' , 'sex' , 'total_overall_unaged' , 'total_overall' ] ] , 
                on = [ stratum_col , 'species_id' , 'sex' ] )
    )      
    # -------- Unaged-specific proportion
    unaged_number_proportion[ 'proportion_number_unaged' ] = (
        unaged_number_proportion[ 'count' ] / unaged_number_proportion[ 'total_overall_unaged' ]
    )
    # -------- Overall proportion
    unaged_number_proportion[ 'proportion_number_overall_unaged' ] = (
        unaged_number_proportion[ 'count' ] / unaged_number_proportion[ 'total_overall' ]
    )
    # ---- Gather unaged (sexed) number proportions
    # -------- Merge
    sex_number_proportions = (
        sex_number_proportions
        .merge( unaged_number_proportion
                .groupby( [ stratum_col , 'species_id' , 'sex' ] )[ 'proportion_number_overall_unaged' ]
                .sum( ).reset_index( ) , how = 'outer' ).fillna( 0.0 )
    )
    # -------- Sum overall total across aged and unaged samples
    sex_number_proportions[ 'proportion_number_overall' ] = (
        sex_number_proportions.proportion_number_overall_aged
        + sex_number_proportions.proportion_number_overall_unaged
    )

    # Calculate overall number proportions across age
    age_number_proportions = (
        aged_number_proportion
        .groupby( [ stratum_col , 'species_id' , 'age_bin' ] , observed = False )[ 'proportion_number_aged' ].sum( )
        .reset_index( name = 'proportion_number' )
    )

    # Return output
    return (
        {
            'age_proportions_df': age_number_proportions ,
            'sex_proportions_df': sex_number_proportions ,
            'aged_length_proportions_df': aged_number_proportion ,
            'unaged_length_proportions_df': unaged_number_proportion ,
        }
    )

def fit_length_weights( proportions_dict: dict , length_weight_dict: dict ) :


    aged_proportions = proportions_dict[ 'aged_length_proportions_df' ]
    unaged_proportions = proportions_dict[ 'unaged_length_proportions_df' ]
    sex_proportions = proportions_dict[ 'sex_proportions_df' ]

    fitted_weight = length_weight_dict[ 'length_weight_regression' ][ 'weight_fitted_df' ]

    # Get the name of the stratum column
    stratum_col = [ col for col in aged_proportions.columns if 'stratum' in col.lower( ) ][ 0 ]

    # Sum number proportions of aged specimens per stratum 
    aged_unaged_proportions = (
        aged_proportions[ aged_proportions.sex == 'all' ]
        .groupby( [ stratum_col ] )[ 'proportion_number_overall_aged' ].sum( )
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
        .merge( sex_proportions , on = [ stratum_col ] )
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
        [ [ stratum_col , 'sex' , 'proportion_aged' , 'proportion_unaged' ] ]
    )

    # Combine the aged-unaged (or station-specific) proportions for calculations
    # ---- Wide-to-long DataFrame
    station_proportions = pd.wide_to_long( stratum_proportions_sexed , 
                                           stubnames = "proportion" , 
                                           i = [ stratum_col , 'sex' ] , 
                                           j = 'group' ,
                                           sep = "_" ,
                                           suffix = "\\w+" ).reset_index( )
    # ---- Convert to Table (to replicate indexed matrix operations)
    station_proportions_table = station_proportions.pivot_table( index = [ 'group' , 'sex' ]  , 
                                                                 columns = [ stratum_col ] , 
                                                                 values = 'proportion' ).fillna( 0.0 )
    
    # Calculate the number length proportions that will be later converted into weight
    # ---- Aged length bins
    aged_length_distribution = (
        aged_proportions
        .groupby( [ stratum_col , 'sex' , 'length_bin' ] , observed = False )
        [ 'proportion_number_aged' ].sum( )
        .reset_index( name = 'number_proportion' )
    )
    # ---- Unaged length bins
    unaged_length_distribution = (
        unaged_proportions[ unaged_proportions.sex != 'unsexed' ]
        [ [ stratum_col , 'sex' , 'length_bin' , 'proportion_number_unaged' ] ]
        .rename( columns = { 'proportion_number_unaged': 'number_proportion' } )
    )
    # ---- Concatenate the two datasets
    length_number_proportions = pd.concat( [ aged_length_distribution.assign( group = 'aged' ) ,
                                             unaged_length_distribution.assign( group = 'unaged' ) ] )    
    # ---- Convert to Table (to replicate indexed matrix operations)
    length_proportions_table = (
        length_number_proportions
        .pivot_table( index = [ 'group' , 'sex' , 'length_bin' ] ,
                     columns = [ stratum_col ] ,
                     values = 'number_proportion' ,
                     observed = False )
        .fillna( 0.0 )
    )
    
    # Convert the fitteed weights into a Table (to replicate index matrix operations)
    fitted_weight_table = fitted_weight.pivot_table( index = [ 'sex' , 'length_intervals' ] ,
                                                     values = 'weight_fitted' ,
                                                     observed = False )
    
    # Calculate the average weights for male, female, and all fish within each stratum
    # ---- All
    weight_all = (
        fitted_weight_table.loc[ 'all' ][ 'weight_fitted' ].values
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
        fitted_weight_table.loc[ 'all' ][ 'weight_fitted' ].values
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
        fitted_weight_table.loc[ 'all' ][ 'weight_fitted' ].values
        .dot(
            length_proportions_table.loc[ 'aged' , 'female' ] 
            * station_proportions_table.loc[ 'aged' , 'female' ]
            + length_proportions_table.loc[ 'unaged' , 'female' ] 
            * station_proportions_table.loc[ 'unaged' , 'female' ]
        )
    ) 
    # ---- Combine the stratum-averaged weights for each sex and all fish
    fitted_weight_df = pd.DataFrame(
        {
            f"{stratum_col}": np.tile( np.unique( station_proportions[ stratum_col ] ) , 
                                    len( np.unique( station_proportions.sex ) ) ) ,
            'sex': np.repeat( [ 'all' , 'male' , 'female' ] , 
                              len( np.unique( station_proportions[ stratum_col ] ) ) ) ,
            'average_weight': np.concatenate( [ weight_all , weight_male , weight_female ] ) ,
        }
    )

    # Return output
    return fitted_weight_df    

def quantize_weights( specimen_data: pd.DataFrame ,
                      length_data: pd.DataFrame ,
                      length_weight_df: pd.DataFrame ) :
    
    # Prepare data
    # ---- Get the name of the stratum column
    stratum_col = [ col for col in specimen_data.columns if 'stratum' in col.lower( ) ][ 0 ]

    # Generate sex-specific interpolators for fitted length-weight values for unaged fish (station 1)
    # ---- Parse the male- and female-specific fitted weight values
    length_weight_sex = (
        length_weight_df.copy( )[ length_weight_df[ 'sex' ].isin( [ 'male' , 'female' ] ) ]
    )
    # ---- Extract 'length' from the interval categories
    length_weight_sex.loc[ : , 'length' ] = (
        length_weight_sex.loc[ : , 'length_intervals' ].apply( lambda x: x.mid )
    )
    # ---- Create interpolator functions
    interpolators = group_interpolator_creator( grouped_data = length_weight_sex ,
                                                independent_var = 'length' ,
                                                dependent_var = 'weight_fitted' ,
                                                contrast = 'sex' )
    # ---- Create helper/lambda function
    def weight_interpolator( dataframe_row ):
        sex = dataframe_row[ 'sex' ]
        length = dataframe_row[ 'length' ]
        if sex in interpolators:
            return interpolators[ sex ]( length )
        else: 
            return None
    # ---- Extract only sexed fish from the unaged (station 1) length dataset
    length_data_sexed = length_data[ length_data[ 'sex' ].isin( [ 'male' , 'female' ] ) ].copy( )
    # ---- Add interpolated weights to the general length dataset
    length_data_sexed.loc[ : , 'weight_interp' ] = (
        length_data_sexed.apply( weight_interpolator , axis = 1 ) 
        * length_data_sexed[ 'length_count' ]
    )
    # ---- Convert interpolated weights (summed across length counts) into a table
    length_table_sexed = (
        length_data_sexed.pivot_table( columns = [ stratum_col , 'sex' ] , 
                                    index = [ 'length_bin' ] , 
                                    values = 'weight_interp' , 
                                    aggfunc = 'sum' , 
                                    observed = False  )
    )

    # Remove specimen data with missing data required for this analysis
    # ---- Drop unsexed fish
    specimen_data_filtered = specimen_data[ specimen_data.group_sex != 'unsexed' ]
    # ---- Remove NaN
    specimen_data_filtered = (
        specimen_data_filtered.dropna( subset = [ 'length' , 'weight' , 'age' ] )
    )
    # ---- Convert to a table
    specimen_table_sexed = (
        specimen_data_filtered
        .pivot_table( columns = [ stratum_col , 'sex' , 'age_bin' ] , 
                    index = [ 'length_bin' ] , 
                    values = 'weight' , 
                    aggfunc = 'sum' , 
                    observed = False )
    )

    # Return a dictionary of the quantized weights
    return (
        {
            'aged_length_weight_tbl': specimen_table_sexed ,
            'unaged_length_weight_tbl': length_table_sexed
        }
    )

def weight_proportions( specimen_data: pd.DataFrame , 
                        catch_data: pd.DataFrame , 
                        proportions_dict: dict ,
                        length_weight_df: pd.DataFrame ,
                        distributions_dict: dict ):
    """
    Calculate the total and sex-specific mean weight for each stratum

    Parameters
    ----------
    species_id : np.float64
        Numeric code representing a particular species of interest

    Notes
    -----
    This function produces the proportion of male and female, 
    and the average weight of male, female, and total (male, female, and unsexed fish).  
    The average weight is estimated using the length-weight relationship 
    fitted in ``fit_binned_length_weight_relationship``.  
    """ 

    # ---- Get the name of the stratum column
    stratum_col = [ col for col in specimen_data.columns if 'stratum' in col.lower( ) ][ 0 ]
    
    # Calculate the sexed and total stratum weights for each sex among aged fish
    # ---- Extract the aged/specimen quantized weights
    aged_weights_binned = distributions_dict[ 'aged_length_weight_tbl' ].copy( )
    # ---- Sum these weights for each sex (male/female)    
    aged_weights_sex = (
        aged_weights_binned.sum( ).unstack( [ 'age_bin' ] ).sum( axis = 1 ).unstack( 0 )
    )
    # ---- Calculate the stratum totals
    aged_strata_weights = (
        aged_weights_sex.sum( ).reset_index( name = 'stratum_weight' )
    )

    # Calculate the sexed and total stratum weights for each sex among unaged fish
    # ---- Sum the net haul weights from station 1/unaged fish
    catch_strata_weights = catch_data.count_variable( contrasts = [ stratum_col ] , 
                                                      variable = 'haul_weight' , 
                                                      fun = 'sum' ) 
    # ---- Rename resulting columns for both
    catch_strata_weights.rename( columns = { 'count': 'stratum_weight' } , inplace = True )

    # Sum the sexed and total weights from the weight-fitted unaged data
    # ---- Extract the unaged/length quantized weights
    unaged_weights_binned =  distributions_dict[ 'unaged_length_weight_tbl' ].copy( )
    # ---- Calculate the total weight per stratum per sex
    unaged_weights_sex = unaged_weights_binned.sum( )
    # ---- Calculate the stratum totals
    unaged_strata_weights = unaged_weights_sex.unstack( 0 ).sum( axis = 0 )
    # ---- Standardize the unaged sexed weights
    unaged_weights_sex_standardized = (
        ( unaged_weights_sex / unaged_strata_weights ).unstack( 0 ) 
        * catch_strata_weights[ 'stratum_weight' ].to_numpy( )
    )

    # Calculate the specimen (aged) weight proportions 
    # ---- Re-pivot the aged weight bins table
    aged_weights_binned_pvt = (
        aged_weights_binned
        .unstack()
        .reset_index( name = 'weight' )
        .pivot_table( columns = [ stratum_col ] , 
                      index = [ 'age_bin' , 'length_bin' , 'sex' ] , 
                      values = 'weight' , 
                      observed = False )
    )
    # ---- Divide by the aged stratum weights (relative to only aged fish)
    aged_weight_proportions_pvt = (
        aged_weights_binned_pvt / aged_strata_weights[ 'stratum_weight' ].to_numpy( )
    )
    # ---- Pivot back to the desired format
    aged_weight_proportions = (
        aged_weight_proportions_pvt.stack()
        .reset_index( name = 'weight_proportion' )
        .pivot_table( columns = [ stratum_col , 'sex' , 'age_bin' ] ,
                      index = 'length_bin' ,
                      values = 'weight_proportion' ,
                      observed = False )
    )
    # ---- Calculate the internal (i.e. only aged fish) for each sex
    within_aged_sex_proportions = aged_weight_proportions.sum( ).unstack( 'age_bin' ).sum( axis = 1 )
    
    # Calculate the total strata weights
    total_strata_weights = (
        pd.concat( [ aged_strata_weights , 
                      catch_strata_weights  ] )
        .pivot_table( columns = [ stratum_col ] , aggfunc = 'sum' , values = 'stratum_weight' ,
                      observed = False )
    )
    
    # Calculate the weight proportions relative to the overall stratum weights
    # ---- Aged
    # -------- Reformat into dataframe and merge with total stratum weights
    aged_weights_binned_df = (
        aged_weights_binned_pvt.stack()
        .to_frame( 'specimen_weight' )
        .reset_index( ).merge( total_strata_weights.T.reset_index( ) , on = stratum_col )
    )
    # -------- Calculate proportions
    aged_weights_binned_df[ 'weight_proportion_overall' ] = (
        aged_weights_binned_df[ 'specimen_weight' ] /
        aged_weights_binned_df[ 'stratum_weight' ]
    )
    # -------- Consolidate to calculate the sexed proportions per stratum
    aged_weight_sex_proportions = (
        aged_weights_binned_df
        .groupby( [ stratum_col , 'sex' ] )[ 'weight_proportion_overall' ].sum( )
    )    
    # ---- Unaged
    # -------- Reformat into dataframe and merge with total stratum weights 
    unaged_weights_sex_standardized_df = (
        unaged_weights_sex_standardized.stack()
        .to_frame('catch_weight')
        .reset_index( ).merge( total_strata_weights.T.reset_index( ) , on = stratum_col )
    )
    # -------- Calculate proportions
    unaged_weights_sex_standardized_df[ 'weight_proportion_overall' ] = (
        unaged_weights_sex_standardized_df[ 'catch_weight' ] /
        unaged_weights_sex_standardized_df[ 'stratum_weight' ]
    )
    # -------- Back-calculate the sexed weight proportions relative to just unaged fish
    # ------------ Aggregate proportions
    unaged_total_sex_proportions = (
        unaged_weights_sex_standardized_df
        .pivot_table( columns = [ 'sex' ] , 
                     index = [ stratum_col ] , 
                     values = 'weight_proportion_overall'  ).sum( axis = 1 )
    )
    # ------------ Re-compute the proportions
    unaged_weight_sex_proportions = (
        unaged_weights_sex_standardized_df
        .pivot_table( index = [ 'sex' ] , 
                      columns = [ stratum_col ] , 
                      values = 'weight_proportion_overall'  )
        / unaged_total_sex_proportions.to_numpy( )
    )

    # Compute the overall length-binned weight distributions among unaged fish
    # ---- Extract the fitted weight values calculated for all fish
    length_weight_all = length_weight_df[ length_weight_df[ 'sex' ] == 'all' ]
    # ---- Generate the fitted weight array
    fitted_weights = length_weight_all[ 'weight_fitted' ].to_numpy( )
    # ---- Extract the number proportions computed for unaged fish
    unaged_number_proportions = proportions_dict[ 'unaged_length_proportions_df' ]
    # ---- Filter out values besides those computed for 'all' fish
    unaged_number_proportions = (
        unaged_number_proportions[ unaged_number_proportions[ 'sex' ] == 'all' ]
    )
    # ---- Convert to a table
    unaged_number_proportions_tbl = (
        unaged_number_proportions.pivot_table( columns = [ stratum_col ] ,
                                               index = [ 'length_bin' ] ,
                                               values = 'proportion_number_unaged' ,
                                               aggfunc = 'sum' ,
                                               observed = False )
    )
    # ---- Apportion the averaged weights
    unaged_apportioned_weights = unaged_number_proportions_tbl.T * fitted_weights
    # ---- Compute the average weight proportions per length bin per stratum
    unaged_length_weights = (
        unaged_apportioned_weights.T / unaged_apportioned_weights.sum( axis = 1 )
    )
    # ---- Convert back to a DataFrame
    unaged_weight_proportions_df = (
        unaged_length_weights.unstack( ).reset_index( name = 'weight_proportion' )
    )
    
    # Calculate the aged and unaged weight proportions
    # ---- Aged
    aged_proportions = aged_weight_sex_proportions.unstack( 'sex' ).sum( axis = 1 )
    # ---- Unaged
    unaged_proportions = 1 - aged_proportions
    # -------- Re-weight the unaged sexed proportions
    unaged_weight_sex_proportions_overall = ( 
        ( unaged_weight_sex_proportions * unaged_proportions ).astype( float ).fillna( 0.0 )
    )
    
    # Format the outputs
    # ---- Aged: stratum-sex-age-length relative to aged and total weights
    aged_overall_df = (
        aged_weight_proportions.unstack( )
        .reset_index( name = 'weight_proportions' )
        .merge( 
            aged_weights_binned_df[ [ stratum_col , 'age_bin' , 'length_bin' , 
                                      'sex' , 'weight_proportion_overall' ] ] )
    )
    # ---- Aged: stratum-sex relative to total weights
    aged_sex_df = (
        within_aged_sex_proportions.reset_index( name = 'weight_proportion_aged' )
        .set_index( [ stratum_col , 'sex' ] )
    )
    # ---- Add the aged sex proportiosn relative to the overall survey
    aged_sex_df[ 'weight_proportion_overall_aged'] = (
        aged_weight_sex_proportions
    )
    # ---- Consolidate the aged and unaged sexed dataframes
    # -------- Initialize the dataframe
    aged_unaged_sex_proportions = aged_sex_df
    # --------- Add the within-unaged weight proportions
    aged_unaged_sex_proportions[ 'weight_proportion_unaged' ] = (
        unaged_weight_sex_proportions.unstack( )
    )
    # --------- Add the overall-unaged weight proportions
    aged_unaged_sex_proportions[ 'weight_proportion_overall_unaged' ] = (
        unaged_weight_sex_proportions_overall.unstack( )
    )
    # ---- Overall aged and unaged proportions
    aged_unaged_proportions = aged_proportions.reset_index( name = 'aged_proportions' )
    # -------- Add unaged proportions
    aged_unaged_proportions[ 'unaged_proportions' ] = unaged_proportions
    
    # Return output
    return (
        {
            'aged_weight_proportions_df': aged_overall_df ,
            'unaged_weight_proportions_df': unaged_weight_proportions_df.reset_index( ) ,
            'aged_unaged_sex_weight_proportions_df': (
                aged_unaged_sex_proportions.astype(float).reset_index().fillna( 0.0 )
            ) ,
            'aged_unaged_weight_proportions_df': aged_unaged_proportions
        }
    )

def age1_metric_proportions( distributions_dict: dict ,
                             proportions_dict: dict ,
                             TS_L_parameters: dict ,
                             settings_dict: dict ):

    # Get stratum column name
    stratum_col = settings_dict[ 'transect' ][ 'stratum_name' ]

    # Calculate the age1 number proportion
    # ---- Initialize the dataframe
    age_proportions = proportions_dict[ 'number' ][ 'aged_length_proportions_df' ].copy( )
    # ---- Only consider 'all' fish
    age_proportions = age_proportions[ ( age_proportions[ 'sex' ] == 'all' ) & ( age_proportions[ stratum_col ] != 0 ) ].reset_index( )
    # ---- Convert to center of length bin
    age_proportions.loc[ : , 'length_mid' ] = age_proportions.loc[ : , 'length_bin' ].apply( lambda df: df.mid ).astype( float )
    # ---- Convert into a pivot table
    age_proportions_table = age_proportions.pivot_table( index = [ 'length_bin' ] , 
                                                         columns = [ stratum_col , 'age_bin' ] , 
                                                         values = 'proportion_number_aged' , 
                                                         aggfunc = 'sum' , 
                                                         observed = False )
    # ---- Sum the number proportions for each age bin within each stratum
    age1_proportions = age_proportions_table.sum( )[ : , 1 ].to_numpy( )

    # Calculate the new length-averaged sigma_bs for each stratum
    # ---- Extract the length-binned values
    length_bins = distributions_dict[ 'length_bins_df' ].copy( )
    # ---- Square the length values
    length_bins[ 'length_sq' ] = length_bins[ 'length_bins' ] ** 2.0
    # ---- Multiply by the TS-length regression coefficient (in the linear domain)
    length_bins[ 'length_sq' ] = (
        length_bins[ 'length_sq' ] * 10 ** ( TS_L_parameters[ 'TS_L_intercept' ] / 10.0 )
    )
    # ---- Repivot the number proportion data
    age_proportions_alt_table = age_proportions.pivot_table( index = [ 'length_bin' ] , 
                                                             columns = [ stratum_col ] , 
                                                             values = 'proportion_number_aged' , 
                                                             aggfunc = 'sum' , 
                                                             observed = False )
    # ---- Dot product to calculate the new average sigma_bs for all ages
    updated_sigma_bs = length_bins[ 'length_sq' ].values.dot( age_proportions_alt_table )

    # 
    # ---- Filter out adult (age-2+ fish)
    age_proportions_age1_table = age_proportions[ age_proportions[ 'age_bin' ] == pd.Interval( left = 0.5 , right = 1.5 ) ]
    # ---- Repivot for just age-1 fish
    age_proportions_age1_table = age_proportions_age1_table.pivot_table( index = [ 'length_bin' ] , 
                                                                        columns = [ stratum_col ] , 
                                                                        values = 'proportion_number_aged' , 
                                                                        aggfunc = 'sum' , 
                                                                        observed = False )
    # ---- Dot product to calculate the average sigma_bs for age-1 fish
    age1_sigma_bs = length_bins[ 'length_sq' ].values.dot( age_proportions_age1_table )
    # ---- Calculate age-1 NASC proportioon per stratum
    age1_nasc_proportions = age1_sigma_bs / updated_sigma_bs

    unage_proportions = proportions_dict[ 'number' ][ 'unaged_length_proportions_df' ]
    # ---- Only consider 'all' fish
    unage_proportions = unage_proportions[ unage_proportions[ 'sex' ] == 'all' ].reset_index( )
    # ---- Convert into a pivot table
    unage_proportions_table = unage_proportions.pivot_table( columns = [ stratum_col ] , 
                                                            index = [ 'length_bin' ] , 
                                                            values = 'proportion_number_unaged' , 
                                                            aggfunc = 'sum' , 
                                                            observed = False )

    min_index = np.where( length_bins[ 'length_bins' ] == 10.0 )[ 0 ]
    max_index = length_bins[ 'length_bins' ].size

    # Calculate thresholds derived from the summed length distributions of age-1 fish
    # ---- General length distribution
    age1_length_distribution_threshold = (
        unage_proportions_table.iloc[ np.arange( min_index[ 0 ] , max_index ) , : ] 
        * age_proportions_age1_table.iloc[ np.arange( min_index[ 0 ] , max_index ) , : ]         
    ).sum( )

    # ---- Just aged length distribution (age-1)
    age1_specific_length_distribution_threshold = age_proportions_age1_table.sum( )

    age_weight_proportions = proportions_dict[ 'weight' ][ 'aged_weight_proportions_df' ]
    age_weight_proportions = age_weight_proportions[ age_weight_proportions[ stratum_col ] != 0 ]
    
    age_weight_proportions_table = age_weight_proportions.pivot_table( index = [ 'length_bin' ] ,
                                                                    columns = [ 'age_bin' , stratum_col ] ,
                                                                    values = 'weight_proportions' ,
                                                                    aggfunc = 'sum' ,
                                                                    observed = False )

    age_weight_proportions_repivot = age_weight_proportions.pivot_table( index = [ 'age_bin' ] ,
                                                                        columns = [ stratum_col ] ,
                                                                        values = 'weight_proportions' ,
                                                                        aggfunc = 'sum' ,
                                                                        observed = False )

    age1_weight_proportions = np.where( ( age1_length_distribution_threshold <= 1e-10 ) & 
                                        ( age1_specific_length_distribution_threshold <= 1e-10 ) ,
                                        0.0 ,
                                        age_weight_proportions_table[ 1 ].sum( )
                                        / age_weight_proportions_repivot.sum( ) )
    
    # Return output
    # ---- Create DataFrame
    apportioned_age1 = pd.DataFrame( { f"{stratum_col}": np.unique( age_proportions[ stratum_col ] ) } )
    # ---- Number proportions
    apportioned_age1[ 'number_proportion' ] = age1_proportions
    # ---- Weight proportions
    apportioned_age1[ 'weight_proportion' ] = age1_weight_proportions
    # ---- NASC
    apportioned_age1[ 'nasc_proportion' ] = age1_nasc_proportions
    # ---- Return
    return apportioned_age1

def distribute_length_age( nasc_biology_df: pd.DataFrame ,
                           proportions_dict: dict ,
                           settings_dict: dict ) :

    # Get the name of the stratum column
    stratum_col = settings_dict[ 'transect' ][ 'stratum_name' ]

    # Extract the correct number proportions
    # ---- Unaged, length (Station 1)
    unaged_number_proportions = proportions_dict[ 'number' ][ 'unaged_length_proportions_df' ]
    # -------- Drop unusexed sex categories
    unaged_number_proportions = (
        unaged_number_proportions[ unaged_number_proportions.sex != 'unsexed' ]
    ) 
    # ---- Aged, length (Station 2)
    aged_number_proportions = proportions_dict[ 'number' ][ 'aged_length_proportions_df' ]
    # -------- Drop unusexed sex categories
    aged_number_proportions = (
        aged_number_proportions[ aged_number_proportions.sex != 'unsexed' ]
    )     

    # Extract the correct number proportions
    # ---- Aged, length (Station 2)
    aged_weight_proportions = proportions_dict[ 'weight' ][ 'aged_weight_proportions_df' ]
    # -------- Sum to create a total/complete key
    aged_weight_proportions_all = (
        aged_weight_proportions
        .groupby( [ stratum_col , 'length_bin' , 'age_bin' ] ,
                  observed = False )[ 'weight_proportions' ]
        .sum( )
        .reset_index( )
    )
    # ---- Concatenate with the full dataset
    aged_weight_proportions = (
        pd.concat( [ aged_weight_proportions.drop( 'weight_proportion_overall' ,
                                                    axis = 1 ) ,
                    aged_weight_proportions_all ] ) 
    )

    # Apportion biomass by sex, length, and age bins
    nasc_abundance = nasc_biology_df[ [ stratum_col , 'transect_num' , 'longitude' , 'latitude' ,
                                    'abundance' ] ]
    # ---- Sum the abundance for each stratum
    abundance_strata = nasc_abundance.groupby( [ stratum_col ] )[ 'abundance' ].sum( )
    # ---- Initialize apportioned unaged abundance
    unaged_apportioned_abundance = unaged_number_proportions.set_index( [ stratum_col ] )
    # ---- Merge with the grouped sex proportions for each stratum
    unaged_apportioned_abundance[ 'abundance_total' ] = abundance_strata
    # ---- Sexed abundance for unaged fish
    unaged_apportioned_abundance[ 'abundance_unaged' ] = (
        unaged_apportioned_abundance[ 'abundance_total' ] 
        * unaged_apportioned_abundance[ 'proportion_number_overall_unaged' ]
    ).fillna( 0.0 )
    # ---- Drop unused columns
    unaged_apportioned_abundance = (
        unaged_apportioned_abundance[ [  'sex' , 'length_bin' , 'abundance_unaged' ] ]
        .reset_index( )
    ) 
    # ---- Convert to pivot table
    unaged_apportioned_abundance_tbl = (
        unaged_apportioned_abundance.pivot_table( index = [ 'sex' , 'length_bin' ] ,
                                                  columns = [ stratum_col ] ,
                                                  values = 'abundance_unaged' ,
                                                  aggfunc = 'sum' ,
                                                  observed = False )
    ) 
    # ---- Initialize apportioned aged abundance
    aged_apportioned_abundance = aged_number_proportions.set_index( [ stratum_col ] )
    # ---- Merge with the grouped sex proportions for each stratum
    aged_apportioned_abundance[ 'abundance_total' ] = abundance_strata
    # ---- Sexed abundance for aged fish
    aged_apportioned_abundance[ 'abundance_aged' ] = (
        aged_apportioned_abundance[ 'abundance_total' ] 
        * aged_apportioned_abundance[ 'proportion_number_overall_aged' ]
    ).fillna( 0.0 )    
    # ---- Drop unused columns
    aged_apportioned_abundance = (
        aged_apportioned_abundance[ [ 'sex' , 'length_bin' , 'age_bin' , 'abundance_aged' ] ]
        .reset_index( )        
    ) 
    # ---- Convert to pivot table
    aged_apportioned_abundance_tbl = (
        aged_apportioned_abundance.pivot_table( index = [ 'sex' , 'length_bin' ] ,
                                                columns = [ stratum_col , 'age_bin' ] ,
                                                values = 'abundance_aged' ,
                                                aggfunc = 'sum' ,
                                                observed = False )
    ) 

    # Apportion abundance by sex, length, and age bins
    # ---- Extract abundance-specific variable
    nasc_biomass = nasc_biology_df[ [ stratum_col , 'transect_num' , 'longitude' , 'latitude' ,
                                     'biomass' , 'biomass_female' , 'biomass_male' ] ].copy( )
    # -------- Adjust column name
    nasc_biomass.rename( columns = { 'biomass': 'biomass_all' } , inplace = True )
    # -------- Pivot from wide to long
    nasc_biomass = pd.wide_to_long( nasc_biomass , 
                                    stubnames = 'biomass' , 
                                    i = [ stratum_col , 'transect_num' , 'longitude' , 'latitude' ] , 
                                    j = 'sex' , sep = "_" , 
                                    suffix=r"\w+" ).reset_index( )    
    # ---- Sum the abundance for each stratum
    biomass_strata = nasc_biomass.groupby( [ stratum_col , 'sex' ] )[ 'biomass' ].sum( )
    # -------- Reset the index
    biomass_strata = biomass_strata.reset_index( )
    # ---- Initialize apportioned aged biomass
    aged_apportioned_biomass = aged_weight_proportions_all.merge( biomass_strata , 
                                                                  on = [ stratum_col ] )
    # ---- Sexed abundance for unaged fish
    aged_apportioned_biomass[ 'biomass_aged' ] = (
        aged_apportioned_biomass[ 'biomass' ] 
        * aged_apportioned_biomass[ 'weight_proportions' ]
    ).fillna( 0.0 )
    # ---- Drop unused columns
    aged_apportioned_biomass = (
        aged_apportioned_biomass[ [ stratum_col , 'sex' , 'length_bin' , 'age_bin' , 
                                    'biomass_aged' ] ]
    ) 
    # ---- Convert to pivot table
    aged_apportioned_biomass_tbl = (
        aged_apportioned_biomass.pivot_table( index = [ 'sex' , 'length_bin' ] ,
                                              columns = [ stratum_col , 'age_bin' ] ,
                                              values = 'biomass_aged' ,
                                              aggfunc = 'sum' ,
                                              observed = False )
    ) 

    # Return outputs
    return (
        {
            'abundance': {
                'aged_abundance_df': aged_apportioned_abundance_tbl ,
                'unaged_abundance_df': unaged_apportioned_abundance_tbl
            } ,
            'biomass': {
                'aged_biomass_df': aged_apportioned_biomass_tbl
            }
        }
    )

def partition_transect_age( nasc_biology_df: pd.DataFrame ,
                            fitted_weight_df: pd.DataFrame ,
                            settings_dict: dict ,
                            population_dict: dict ,
                            strata_adult_proportions_df: pd.DataFrame ) :
    
    # Get the name of the stratum column
    stratum_col = settings_dict[ 'transect' ][ 'stratum_name' ]

    # Procedure where age-1 fish were excluded from the acoustic survey biological distributions
    # ---- Merge adult proportions with acoustically derived georeferenced biological data
    adult_data = nasc_biology_df.merge( strata_adult_proportions_df , on = [ stratum_col ] )
    # ---- Adjust the full number density, biomass density, and NASC estimates
    # -------- NASC
    adult_data[ 'nasc' ] = (
        adult_data[ 'nasc' ] * adult_data[ 'nasc_proportion' ]
    )   
    # -------- Number density
    adult_data[ 'number_density' ] = (
        adult_data[ 'number_density' ] * adult_data[ 'number_proportion' ]
    )
    adult_data[ 'number_density_female' ] = (
        adult_data[ 'number_density_female' ] * adult_data[ 'number_proportion' ]
    )
    adult_data[ 'number_density_male' ] = (
        adult_data[ 'number_density_male' ] * adult_data[ 'number_proportion' ]
    )
    # -------- Abundance
    adult_data[ 'abundance' ] = (
        adult_data[ 'abundance' ] * adult_data[ 'number_proportion' ]
    )
    adult_data[ 'abundance_female' ] = (
        adult_data[ 'abundance_female' ] * adult_data[ 'number_proportion' ]
    )
    adult_data[ 'abundance_male' ] = (
        adult_data[ 'abundance_male' ] * adult_data[ 'number_proportion' ]
    )
    # -------- Biomass density
    adult_data[ 'biomass_density' ] = (
        adult_data[ 'biomass_density' ] * adult_data[ 'weight_proportion' ]
    )  
    adult_data[ 'biomass_density_female' ] = (
        adult_data[ 'biomass_density_female' ] * adult_data[ 'weight_proportion' ]
    )  
    adult_data[ 'biomass_density_male' ] = (
        adult_data[ 'biomass_density_male' ] * adult_data[ 'weight_proportion' ]
    )  
    # -------- Biomass
    adult_data[ 'biomass' ] = (
        adult_data[ 'biomass' ] * adult_data[ 'weight_proportion' ]
    )  
    adult_data[ 'biomass_female' ] = (
        adult_data[ 'biomass_female' ] * adult_data[ 'weight_proportion' ]
    )  
    adult_data[ 'biomass_male' ] = (
        adult_data[ 'biomass_male' ] * adult_data[ 'weight_proportion' ]
    )  
    # ---- Drop unused column names
    adult_data = adult_data.filter( regex = "^(?!.*(_proportion|_unsexed))")

    # Adjust the population abundance tables to include only age-1 fish
    # ---- Extract abundances distributed for unaged lengths
    abundance_unaged_length = population_dict[ 'tables' ][ 'abundance' ][ 'unaged_abundance_df' ]
    # ---- Convert `strata_adult_proportions_df` into a similarly indexed table
    strata_adult_proportions_table = (
        strata_adult_proportions_df.pivot_table( index = stratum_col ).T
    )
    # ---- Convert the table to represent age-1 proportions
    strata_age1_proportions_table = 1.0 - strata_adult_proportions_table
    # ---- Distribute the age-1 proportions across the unaged abundances
    abundance_unaged_age1_tbl = (
        strata_age1_proportions_table.loc[ 'number_proportion' ] * abundance_unaged_length
    )

    fitted_weight = fitted_weight_df[ 'length_weight_regression' ][ 'weight_fitted_df' ]

    # ---- Extract abundances distributed for unaged lengths
    abundance_aged_length = population_dict[ 'tables' ][ 'abundance' ][ 'aged_abundance_df' ]
    # ---- Reindex
    abundance_age1_length = (
        abundance_aged_length.T
        .reset_index( )
        .set_index( [ 'age_bin' , stratum_col ] ).loc[ 1 ].T
    )
    # ---- Extract fitted weights for all fish
    fitted_weight_all = fitted_weight[ fitted_weight.sex == 'all' ].drop( 'sex' , axis = 1 )
    # ---- Calculate the summed biomass across all strata for age-1 fish
    # -------- All fish
    fitted_biomass_all = (
        fitted_weight[ fitted_weight.sex == 'all' ]
        .drop( 'sex' , axis = 1 )[ 'weight_fitted' ]
        .values.dot( abundance_age1_length.loc[ 'all' , : ]).sum( )
    )
    # -------- Female fish
    fitted_biomass_female = (
        fitted_weight[ fitted_weight.sex == 'all' ]
        .drop( 'sex' , axis = 1 )[ 'weight_fitted' ]
        .values.dot( abundance_age1_length.loc[ 'female' , : ]).sum( )
    )
    # -------- Male fish
    fitted_biomass_male = (
        fitted_weight[ fitted_weight.sex == 'all' ]
        .drop( 'sex' , axis = 1 )[ 'weight_fitted' ]
        .values.dot( abundance_age1_length.loc[ 'male' , : ]).sum( )
    )
    # ---- Calculate the estimated biomass proportions/contributions
    female_ratio = fitted_biomass_female / fitted_biomass_all
    male_ratio = fitted_biomass_male / fitted_biomass_all
    # ---- Get the age-1 biomass from `population_dict`
    biomass_aged_length = population_dict[ 'tables' ][ 'biomass' ][ 'aged_biomass_df' ]
    # ---- Reindex and extract just the age-1 biomass estimates
    biomass_age1_length = (
        biomass_aged_length.T
        .reset_index( )
        .set_index( [ 'age_bin' , stratum_col ] ).loc[ 1 ].T
    )
    # ---- Calculate the age-1 biomass among male fish
    # -------- Sum the total age-1 biomass
    biomass_age1_total = biomass_age1_length.loc[ 'all' , ].sum( ).sum( )
    # -------- Sum the total female age-1 biomass
    biomass_age1_female = biomass_age1_total * female_ratio
    # -------- Sum the total male age-1 biomass
    biomass_age1_male = biomass_age1_total * male_ratio
    # -------- Sum the mixed age-1 biomass
    biomass_age1_mixed = biomass_age1_total - biomass_age1_female - biomass_age1_male

    # Return outputs
    # ---- Convert biomass estimates into a summary table
    # -------- Initialize dataframe
    biomass_summary_df = pd.DataFrame( { 'sex': [ 'all' , 'female' , 'male' , 'unsexed' , 'mixed' ] } )
    # -------- Calculate biomass estimates for age-2+ fish
    if settings_dict[ 'transect' ][ 'exclude_age1' ]:
        biomass_summary_df[ 'biomass_adult' ] = (
            np.array( [ biomass_aged_length.loc[ 'all' , : ].sum( ).sum( ) - biomass_age1_total ,
                        biomass_aged_length.loc[ 'female' , : ].sum( ).sum( ) - biomass_age1_female ,
                        biomass_aged_length.loc[ 'male' , : ].sum( ).sum( ) - biomass_age1_male ,
                        nasc_biology_df[ 'biomass_unsexed' ].sum( ) ,
                        nasc_biology_df[ nasc_biology_df[ 'fraction_hake' ] < 1.0 ][ 'biomass' ].sum( )
                        - biomass_age1_mixed ] )
        )
    else: 
        biomass_summary_df[ 'biomass_adult' ] = (
            np.array( [ adult_data[ 'biomass' ].sum( ) ,
                        adult_data[ 'biomass_female' ].sum( ) ,
                        adult_data[ 'biomass_male' ].sum( ) ,
                        nasc_biology_df[ 'biomass_unsexed' ].sum( ) ,
                        nasc_biology_df[ nasc_biology_df[ 'fraction_hake' ] < 1.0 ][ 'biomass' ].sum( ) ] )
        )
    # -------- Calculate biomass estimates for age-1 fish
    biomass_summary_df[ 'biomass_age1' ] = (
        np.array( [ biomass_age1_total ,
                    biomass_age1_female ,
                    biomass_age1_male ,
                    0.0 ,
                    biomass_age1_mixed ] )
    )
    # -------- Calculate biomass estimates for all fish
    biomass_summary_df[ 'biomass_all' ] = (
        biomass_summary_df[ 'biomass_adult' ] + biomass_summary_df[ 'biomass_age1' ]
    )    
    # ---- Generate outputs
    return adult_data , biomass_summary_df , abundance_unaged_age1_tbl


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