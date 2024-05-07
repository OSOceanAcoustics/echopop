import numpy as np
import pandas as pd

def ts_length_regression( x, slope, intercept ):
    return slope * np.log10(x) + intercept

def to_linear( x ):
    return 10.0 ** ( x / 10.0 )

def to_dB( x ):
    return 10 * np.log10(x)

def impute_missing_sigma_bs( strata_dict: dict , sigma_bs_ks: pd.DataFrame , sigma_bs_inpfc: pd.DataFrame ):
    """
    Imputes sigma_bs for strata without measurements or values

    Parameters
    ----------
    strata_dict: dict
    sigma_bs_ks: pd.DataFrame
    sigma_bs_inpfc: pd.DataFrame            

    Notes
    -----
    This function iterates through all stratum layers to impute either the
    nearest neighbor or mean sigma_bs for strata that are missing values.
    """    

    # Collect all possible strata
    # ---- KS
    strata_options_ks = np.unique( strata_dict[ 'strata_df' ].stratum_num )
    # ---- INPFC
    strata_options_inpfc = np.unique( strata_dict[ 'inpfc_strata_df' ].stratum_inpfc )

    # Collect present strata
    # ---- KS
    present_strata_ks = np.unique(sigma_bs_ks[ 'stratum_num' ]).astype(int)
    # ---- INPFC
    present_strata_inpfc = np.unique(sigma_bs_inpfc[ 'stratum_inpfc' ]).astype(int)

    # Invert to retrieve missing strata
    # ---- KS
    missing_strata_ks = strata_options_ks[ ~ np.isin( strata_options_ks , present_strata_ks ) ]
    # -------- Concatenate the existing data with the missing strata
    sigma_bs_ks_impute = (
            pd.concat( [ sigma_bs_ks , 
                            pd.DataFrame( {
                                'stratum_num': missing_strata_ks , 
                                'species_id': np.repeat( np.unique( sigma_bs_ks.species_id ) ,
                                                        len( missing_strata_ks ) ) ,
                                'sigma_bs_mean': np.repeat( np.nan ,
                                                            len( missing_strata_ks ) )
                            } ) ] )
            .sort_values( 'stratum_num' )        
        )
    # ---- INPFC
    missing_strata_inpfc = strata_options_inpfc[ ~ np.isin( strata_options_inpfc , present_strata_inpfc ) ]
    # -------- Concatenate the existing data with the missing strata
    sigma_bs_inpfc_impute = (
            pd.concat( [ sigma_bs_inpfc , 
                            pd.DataFrame( {
                                'stratum_inpfc': missing_strata_inpfc , 
                                'species_id': np.repeat( np.unique( sigma_bs_inpfc.species_id ) ,
                                                        len( missing_strata_inpfc ) ) ,
                                'sigma_bs_mean': np.repeat( np.nan ,
                                                            len( missing_strata_inpfc ) )
                            } ) ] )
            .sort_values( 'stratum_inpfc' )        
        )

    # Impute values for missing strata
    # ---- KS 
    if len( missing_strata_ks ) > 0 :
        # ---- Find strata intervals to impute over        
        for i in missing_strata_ks:
            strata_floor = present_strata_ks[present_strata_ks < i]
            strata_ceil = present_strata_ks[present_strata_ks > i]

            new_stratum_below = np.max(strata_floor) if strata_floor.size > 0 else None
            new_stratum_above = np.min(strata_ceil) if strata_ceil.size > 0 else None      
            
            sigma_bs_indexed = sigma_bs_ks_impute[sigma_bs_ks_impute['stratum_num'].isin([new_stratum_below, new_stratum_above])]
            
            sigma_bs_ks_impute.loc[sigma_bs_ks_impute.stratum_num==i , 'sigma_bs_mean' ] = sigma_bs_indexed[ 'sigma_bs_mean' ].mean()
    # ---- INPFC
    if len( missing_strata_inpfc ) > 0 :
        # ---- Find strata intervals to impute over        
        for i in missing_strata_inpfc:
            strata_floor = present_strata_inpfc[present_strata_inpfc < i]
            strata_ceil = present_strata_inpfc[present_strata_inpfc > i]

            new_stratum_below = np.max(strata_floor) if strata_floor.size > 0 else None
            new_stratum_above = np.min(strata_ceil) if strata_ceil.size > 0 else None      
            
            sigma_bs_indexed = sigma_bs_inpfc_impute[sigma_bs_inpfc_impute['stratum_inpfc'].isin([new_stratum_below, new_stratum_above])]
            
            sigma_bs_inpfc_impute.loc[sigma_bs_inpfc_impute.stratum_inpfc==i , 'sigma_bs_mean' ] = sigma_bs_indexed[ 'sigma_bs_mean' ].mean()

    # Return output
    return sigma_bs_ks_impute , sigma_bs_inpfc_impute               

def summarize_sigma_bs( length_data: pd.DataFrame , specimen_data: pd.DataFrame , strata_dict: dict , configuration_dict: dict ) :
    """
    Calculates the stratified mean sigma_bs for each stratum

    Parameters
    ----------
    length_data: pd.DataFrame
    specimen_data: pd.DataFrame
    strata_dict: dict
    configuration_dict: dict

    Notes
    -----
    This function iterates through each stratum to fit acoustic target 
    strength (TS, dB re. 1 m^2) values based on length distributions recorded for each
    stratum. These fitted values are then convereted from the logarithmic
    to linear domain (sigma_bs, m^2) and subsequently averaged. These are required for 
    later functions that will convert vertically integrated backscatter (e.g. the nautical
    area scattering coefficient, or NASC, m^2 nmi^-2) to estimates of areal density (animals nmi^-2).
    """

    # Meld the specimen and length dataframes together for downstream calculations
    aggregate_lengths = (
        specimen_data
        .meld( length_data , contrasts = [ 'haul_num' , 'stratum_num' , 'stratum_inpfc' , 
                                           'species_id' , 'length' , 'group_sex' ] )
    )

    # Re-group the length-specific length counts
    aggregate_lengths = (
        aggregate_lengths
        .groupby( [ 'haul_num' , 'stratum_num' , 'stratum_inpfc' , 'species_id' , 'length' , 'group_sex' ] )
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
        .groupby( [ 'species_id' , 'haul_num' , 'stratum_num' , 'stratum_inpfc' ] )[ [ 'sigma_bs' , 'length_count' ] ]
        .apply( lambda df: np.average( df.sigma_bs , weights = df.length_count ) )
        .to_frame( 'sigma_bs_mean' )
        .reset_index( )
    ) 
    # ---- By KS stratum
    sigma_bs_ks = (
        sigma_bs_haul
        .groupby( [ 'stratum_num' , 'species_id' ] )[ 'sigma_bs_mean' ]
        .mean( )
        .reset_index( )
    )
    # ---- By INPFC stratum
    sigma_bs_inpfc = (
        sigma_bs_haul
        .groupby( [ 'stratum_inpfc' , 'species_id' ] )[ 'sigma_bs_mean' ]
        .mean( )
        .reset_index( )
    )

    # Impute sigma_bs values, if necessary, for missing strata
    sigma_bs_ks_updated , sigma_bs_inpfc_updated = impute_missing_sigma_bs( strata_dict ,
                                                                            sigma_bs_ks ,
                                                                            sigma_bs_inpfc )

    # Calculate mean TS for all hauls, KS-strata, and INPFC strata
    # ---- By haul
    sigma_bs_haul[ 'TS_mean' ] = to_dB( sigma_bs_haul[ 'sigma_bs_mean' ] )
    # ---- By KS stratum
    sigma_bs_ks_updated[ 'TS_mean' ] = to_dB( sigma_bs_ks_updated[ 'sigma_bs_mean' ] )
    # ---- By INPFC stratum
    sigma_bs_inpfc_updated[ 'TS_mean' ] = to_dB( sigma_bs_inpfc_updated[ 'sigma_bs_mean' ] )

    # Return output
    return (
        {
            'haul_mean_df': sigma_bs_haul ,
            'inpfc_mean_df': sigma_bs_inpfc_updated ,
            'strata_mean_df': sigma_bs_ks_updated ,
        }
    )
