import numpy as np
import pandas as pd

from ..spatial.transect import correct_transect_intervals
from ..biology import age1_metric_proportions

def ts_length_regression( x, slope, intercept ):
    return slope * np.log10(x) + intercept

def to_linear( x ):
    return 10.0 ** ( x / 10.0 )

def to_dB( x ):
    return 10 * np.log10(x)

def impute_missing_sigma_bs( strata_dict: dict , sigma_bs_stratum: pd.DataFrame ):
    """
    Imputes sigma_bs for strata without measurements or values

    Parameters
    ----------
    strata_dict: dict
    sigma_bs_stratum: pd.DataFrame    

    Notes
    -----
    This function iterates through all stratum layers to impute either the
    nearest neighbor or mean sigma_bs for strata that are missing values.
    """    

    # Collect stratum column name
    stratum_col = [ col for col in sigma_bs_stratum.columns if 'stratum' in col.lower( ) ][ 0 ]

    # Collect all possible strata
    strata_options = np.unique( strata_dict[ 'strata_df' ][ stratum_col ] )

    # Collect present strata
    present_strata = np.unique(sigma_bs_stratum[ stratum_col ]).astype(int)

    # Invert to retrieve missing strata
    # ---- KS
    missing_strata = strata_options[ ~ np.isin( strata_options , present_strata ) ]
    # -------- Concatenate the existing data with the missing strata
    sigma_bs_stratum_impute = (
            pd.concat( [ sigma_bs_stratum , 
                            pd.DataFrame( {
                                f"{stratum_col}": missing_strata , 
                                'species_id': np.repeat( np.unique( sigma_bs_stratum.species_id ) ,
                                                        len( missing_strata ) ) ,
                                'sigma_bs_mean': np.repeat( np.nan ,
                                                            len( missing_strata ) )
                            } ) ] )
            .sort_values( stratum_col )        
        )

    # Impute values for missing strata
    # ---- KS 
    if len( missing_strata ) > 0 :
        # ---- Find strata intervals to impute over        
        for i in missing_strata:
            strata_floor = present_strata[present_strata < i]
            strata_ceil = present_strata[present_strata > i]

            new_stratum_below = np.max(strata_floor) if strata_floor.size > 0 else None
            new_stratum_above = np.min(strata_ceil) if strata_ceil.size > 0 else None      
            
            sigma_bs_indexed = sigma_bs_stratum_impute[sigma_bs_stratum_impute[stratum_col].isin([new_stratum_below, new_stratum_above])]
            
            sigma_bs_stratum_impute.loc[sigma_bs_stratum_impute[stratum_col]==i , 'sigma_bs_mean' ] = sigma_bs_indexed[ 'sigma_bs_mean' ].mean()

    # Return output
    return sigma_bs_stratum_impute         

def summarize_sigma_bs( length_data: pd.DataFrame , specimen_data: pd.DataFrame , strata_dict: dict , configuration_dict: dict , settings_dict: dict ) :
    """
    Calculates the stratified mean sigma_bs for each stratum

    Parameters
    ----------
    length_data: pd.DataFrame
    specimen_data: pd.DataFrame
    strata_dict: dict
    configuration_dict: dict
    settings_dict: dict

    Notes
    -----
    This function iterates through each stratum to fit acoustic target 
    strength (TS, dB re. 1 m^2) values based on length distributions recorded for each
    stratum. These fitted values are then convereted from the logarithmic
    to linear domain (sigma_bs, m^2) and subsequently averaged. These are required for 
    later functions that will convert vertically integrated backscatter (e.g. the nautical
    area scattering coefficient, or NASC, m^2 nmi^-2) to estimates of areal density (animals nmi^-2).
    """

    # Select the appropriate stratum column name
    stratum_col = settings_dict[ 'transect' ][ 'stratum_name' ]

    # Meld the specimen and length dataframes together for downstream calculations
    aggregate_lengths = (
        specimen_data
        .meld( length_data , contrasts = [ 'haul_num' , stratum_col , 
                                           'species_id' , 'length' , 'group_sex' ] )
    )

    # Re-group the length-specific length counts
    aggregate_lengths = (
        aggregate_lengths
        .groupby( [ 'haul_num' , stratum_col , 'species_id' , 'length' , 'group_sex' ] )
        [ 'length_count' ].sum( )
        .reset_index( )
    )    

    # Extract TS-length regression coefficients
    # ---- Pull in TS-length regression dictionary
    ts_length_parameters = configuration_dict[ 'TS_length_regression_parameters' ]
    # ---- Sub-sample regression coefficients for target species
    ts_length_parameters_spp = (
        [ spp for spp in ts_length_parameters.values( ) if spp[ 'number_code' ] in np.unique( aggregate_lengths.species_id ).astype( int ) ]
    )
    # ---- Create DataFrame containing all necessary regression coefficients
    ts_length_parameters_df = pd.DataFrame.from_dict( ts_length_parameters_spp )
    # ---- Merge with the biological data
    ts_lengths_df = (
        aggregate_lengths
        .merge( ts_length_parameters_df.drop( 'length_units' , axis = 1 ) , 
                left_on = [ 'species_id' ] , right_on = [ 'number_code' ] )
    )

    # Calculate predicted TS from the length values
    ts_lengths_df[ 'TS' ] = ts_length_regression( ts_lengths_df[ 'length' ] , 
                                                  ts_lengths_df[ 'TS_L_slope' ] ,
                                                  ts_lengths_df[ 'TS_L_intercept' ] )

    # Convert to the linear domain ('sigma_bs')
    ts_lengths_df[ 'sigma_bs' ] = to_linear( ts_lengths_df[ 'TS' ] )

    # Calculate mean sigma_bs for all hauls, KS-strata, and INPFC strata
    # ---- By haul
    sigma_bs_haul = (
        ts_lengths_df
        .groupby( [ 'species_id' , 'haul_num' , stratum_col ] )[ [ 'sigma_bs' , 'length_count' ] ]
        .apply( lambda df: np.average( df.sigma_bs , weights = df.length_count ) )
        .to_frame( 'sigma_bs_mean' )
        .reset_index( )
    ) 
    # ---- By stratum
    sigma_bs_stratum = (
        sigma_bs_haul
        .groupby( [ stratum_col , 'species_id' ] )[ 'sigma_bs_mean' ]
        .mean( )
        .reset_index( )
    )

    # Impute sigma_bs values, if necessary, for missing strata
    sigma_bs_stratum_updated = impute_missing_sigma_bs( strata_dict ,
                                                        sigma_bs_stratum )

    # Calculate mean TS for all hauls, KS-strata, and INPFC strata
    # ---- By haul
    sigma_bs_haul[ 'TS_mean' ] = to_dB( sigma_bs_haul[ 'sigma_bs_mean' ] )
    # ---- By KS stratum
    sigma_bs_stratum_updated[ 'TS_mean' ] = to_dB( sigma_bs_stratum_updated[ 'sigma_bs_mean' ] )

    # Return output
    return (
        {
            'haul_mean_df': sigma_bs_haul ,
            'strata_mean_df': sigma_bs_stratum_updated
        }
    )

def nasc_to_biomass( input_dict: dict ,
                     analysis_dict: dict ,
                     configuration_dict: dict ,
                     settings_dict: dict ) :
    """
    Converts integrated acoustic backscatter (NASC) into estimates of 
    areal number/biomass densities, total abundance, and total biomass

    Parameters
    ----------

    Notes
    -----
    This function converts NASC into estimates of population-level metrics 
    (abundance, biomass, areal densities) stratified by transects, sex, 
    length-based strata, and age.
    """ 

    # Extract the necessary correct strata mean sigma_bs
    sigma_bs_strata = analysis_dict[ 'acoustics' ][ 'sigma_bs' ][ 'strata_mean_df' ]

    # Pull out the length-weight conversion for each stratum
    length_weight_strata = analysis_dict[ 'biology' ][ 'weight' ]['weight_stratum_df']

    # Get the name of the stratum column
    stratum_col = settings_dict[ 'transect' ][ 'stratum_name' ]

    # Correct the acoustic survey transect intervals
    nasc_interval_df = correct_transect_intervals( input_dict[ 'acoustics' ][ 'nasc_df' ] )

    # Select the appropriate NASC column based on the inclusion or exclusion of age-1 fish
    if settings_dict[ 'transect' ][ 'exclude_age1' ] :
        # ---- Drop the unused column
        nasc_interval_df.drop( 'NASC_all_ages' , axis = 1 , inplace = True )
        # ---- Rename the used column
        nasc_interval_df.rename( columns = { 'NASC_no_age1': 'nasc' } , inplace = True )
        # ---- Calculate age-1 NASC and weight proportions
        age1_proportions = (
            age1_metric_proportions( 
                input_dict[ 'biology' ][ 'distributions' ] ,
                analysis_dict[ 'biology' ][ 'proportions' ] ,
                configuration_dict[ 'TS_length_regression_parameters' ][ 'pacific_hake' ] ,
                settings_dict 
            )
        )
        # ---- Calculate adult proportions
        # -------- Initialize dataframe
        adult_proportions = age1_proportions.copy( )
        # -------- Invert proportions
        adult_proportions.iloc[ : , np.arange( 1 , 4 ) ] = (
            1.0 - adult_proportions.iloc[ : , np.arange( 1 , 4 ) ]
        )

    else:
        # ---- Drop the unused column
        nasc_interval_df.drop( 'NASC_no_age1' , axis = 1 , inplace = True )
        # ---- Rename the used column
        nasc_interval_df.rename( columns = { 'NASC_all_ages': 'nasc' } , inplace = True )
        # ---- Assign filled adult proportions dataframe
        adult_proportions = pd.DataFrame( { f"{stratum_col}": np.unique( nasc_interval_df[ stratum_col ] ) ,
                                            "number_proportion": np.ones( len( np.unique( nasc_interval_df[ stratum_col ] ) ) ) ,
                                            "number_proportion": np.ones( len( np.unique( nasc_interval_df[ stratum_col ] ) ) ) ,
                                            "number_proportion": np.ones( len( np.unique( nasc_interval_df[ stratum_col ] ) ) ) } )

    # Merge hake fraction data into `nasc_interval_df`
    # ---- Initial merge
    nasc_interval_df = (
        nasc_interval_df
        .merge( input_dict[ 'spatial' ][ 'strata_df' ] , 
                on = [ stratum_col , 'haul_num' ] , 
                how = 'outer' )
    )
    # ---- Replace `fraction_hake` where NaN occurs
    nasc_interval_df[ 'fraction_hake' ] = nasc_interval_df[ 'fraction_hake' ].fillna( 0.0 )
    # ---- Drop NaN
    nasc_interval_df.dropna( subset = [ 'transect_num' ] , inplace = True )

    # Calculate the along-transect number density (animals per nmi^2)
    # ---- Merge NASC measurements with mean sigma_bs for each stratum
    nasc_biology = nasc_interval_df.merge( sigma_bs_strata , on = [ stratum_col ] )
    # ---- Calculate the number densities
    nasc_biology[ 'number_density' ] = (
        nasc_biology[ 'fraction_hake' ]
        * nasc_biology[ 'nasc' ]
        / ( 4.0 * np.pi * nasc_biology[ 'sigma_bs_mean' ] ) 
    )
    # ---- Round to the nearest whole number (we don't want fractions of a fish)
    nasc_biology[ 'number_density' ] = np.round( nasc_biology[ 'number_density' ] )
    
    # Calculate the along-transect abundance (# animals)
    nasc_biology[ 'abundance' ] = nasc_biology[ 'interval_area' ] * nasc_biology[ 'number_density' ]

    # Calculate the along-transect abundance (#, N) and biomass for each sex (kg)
    # ---- Extract the sex-specific number proportions for each stratum
    sex_stratum_proportions = analysis_dict[ 'biology' ][ 'proportions' ][ 'number' ][ 'sex_proportions_df' ].copy( )
    # ---- Filter out sexes besides 'male' and 'female'
    sex_stratum_proportions = sex_stratum_proportions[ sex_stratum_proportions.sex.isin( [ 'male' , 'female' ] ) ]
    # ---- Merge with the NASC measurements
    nasc_biology_sex = nasc_biology.merge( sex_stratum_proportions , on = [ stratum_col , 'species_id' ] )
    # ---- Apportion number density by sex (animals/nmi^2)
    nasc_biology_sex[ 'number_density_sex' ] = (
        np.round( nasc_biology_sex[ 'number_density' ] * nasc_biology_sex[ 'proportion_number_overall' ] )
    )
    # ---- Calculate/update the sex-specific abundance estimates
    nasc_biology_sex[ 'abundance_sex' ] = (
        nasc_biology_sex[ 'abundance' ] * nasc_biology_sex[ 'proportion_number_overall' ]        
    )
    # ---- Merge with sex-specific average weights per stratum
    nasc_biology_sex = nasc_biology_sex.merge( length_weight_strata  , 
                                               on = [ stratum_col , 'sex' ] )
    # ---- Calculate biomass density (kg/nmi^2)
    nasc_biology_sex[ 'biomass_density_sex' ] = (
        nasc_biology_sex[ 'number_density_sex' ] * nasc_biology_sex[ 'average_weight' ]
    )    
    # ---- Calculate biomass for each sex
    nasc_biology_sex[ 'biomass_sex' ] = (
        nasc_biology_sex[ 'abundance_sex' ] * nasc_biology_sex[ 'average_weight' ]
    )

    # Parse out the biomass estimates from unsexed fish 
    # ---- Pivot from long-to-wide
    nasc_biology_grp = nasc_biology_sex.pivot( index = [ stratum_col , 'transect_num' , 'longitude' , 
                                                         'fraction_hake' , 'latitude' , 'nasc' , 
                                                         'abundance' , 'number_density' ] ,
                                                columns = [ 'sex' ] ,
                                                values = [ 'abundance_sex' , 'biomass_sex' , 
                                                           'number_density_sex' , 'biomass_density_sex' ] )
    # ---- Collapse the column names
    nasc_biology_grp.columns = nasc_biology_grp.columns.map('_'.join).str.strip('|')
    # ---- Reset the index
    nasc_biology_grp.reset_index( inplace = True )
    # ---- Merge with the average weights per strata for all fish
    nasc_biology_grp = (
        nasc_biology_grp
        .merge( length_weight_strata[ length_weight_strata.sex == 'all' ].drop( 'sex' , axis = 1 ) ,
                on = [ stratum_col ] )
    ) 
    # ---- Calculate unsexed number density
    nasc_biology_grp[ 'number_density_unsexed' ] = (
        nasc_biology_grp[ 'number_density' ]
        - nasc_biology_grp[ 'number_density_sex_male' ] 
        - nasc_biology_grp[ 'number_density_sex_female' ]
    )
    # ---- Calculate unsexed biomass density
    nasc_biology_grp[ 'biomass_density_unsexed' ] = (
        nasc_biology_grp[ 'number_density_unsexed' ] * nasc_biology_grp[ 'average_weight' ]
    )    
    # ---- Calculate unsexed abundance
    nasc_biology_grp[ 'abundance_unsexed' ] = (
        nasc_biology_grp[ 'abundance' ] 
        - nasc_biology_grp[ 'abundance_sex_male' ] - nasc_biology_grp[ 'abundance_sex_female' ]
    )
    # ---- Calculate unsexed biomass
    nasc_biology_grp[ 'biomass_unsexed' ] = (
        nasc_biology_grp[ 'abundance_unsexed' ] * nasc_biology_grp[ 'average_weight' ]
    ) 
    # ---- Calculate sexed biomass
    nasc_biology_grp[ 'biomass_sexed' ] = (
        nasc_biology_grp[ 'biomass_sex_male' ] + nasc_biology_grp[ 'biomass_sex_female' ]
    )     
    # ---- Sum the total biomass density
    nasc_biology_grp[ 'biomass_density' ] = (
        nasc_biology_grp[ 'biomass_density_sex_female' ]
        + nasc_biology_grp[ 'biomass_density_sex_male' ]
        + nasc_biology_grp[ 'biomass_density_unsexed' ]
    )
    # ---- Sum the total biomass
    nasc_biology_grp[ 'biomass' ] = (
        nasc_biology_grp[ 'biomass_sex_female' ] + nasc_biology_grp[ 'biomass_sex_male' ]
        + nasc_biology_grp[ 'biomass_unsexed' ]
    )

    # Edit the output dataframe format
    # ---- Overall biological distributions
    # -------- Edit column names
    nasc_biology_grp.columns = nasc_biology_grp.columns.str.replace('_sex_', '_')
    # -------- Explicitly drop unnecessary columns
    nasc_biology_grp.drop( [ 'average_weight' , 'biomass_sexed' ] , axis = 1 , inplace = True )
    # -------- Rearrange column order
    nasc_biology_grp = nasc_biology_grp.reindex( columns = [ stratum_col , 'transect_num' , 
                                                            'longitude' , 'latitude' ,
                                                            'nasc' , 'fraction_hake' , 
                                                            'number_density' , 
                                                            'number_density_female' , 
                                                            'number_density_male' , 
                                                            'number_density_unsexed' , 'abundance' , 
                                                            'abundance_female' , 'abundance_male' ,
                                                            'abundance_unsexed' , 'biomass_density' , 
                                                            'biomass_density_female' , 
                                                            'biomass_density_male' ,
                                                            'biomass_density_unsexed' , 'biomass' , 
                                                            'biomass_female' , 'biomass_male' , 
                                                            'biomass_unsexed' ] )
    # ---- Return output
    return adult_proportions , nasc_biology_grp
