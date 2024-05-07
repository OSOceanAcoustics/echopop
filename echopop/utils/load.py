from pathlib import Path
from openpyxl import load_workbook
import yaml
import numpy as np
import pandas as pd

def load_configuration( init_config_path: Path , 
                        survey_year_config_path: Path ):
    """
    Loads the biological, NASC, and stratification
    data using parameters obtained from the configuration
    files.

    Parameters
    ----------
    init_config_path : pathlib.Path
        A string specifying the path to the initialization YAML file
    survey_year_config_path : pathlib.Path
        A string specifying the path to the survey year YAML file

    Notes
    -----
    This function parses the configuration files and incorporates them into
    the Survey class object. This initializes the `config` attribute that 
    becomes available for future reference and functions.
    """
    ### Validate configuration files
    # Retreive the module directory to begin mapping the configuration file location
    #current_directory = os.path.dirname(os.path.abspath(__file__))

    # Build the full configuration file paths and verify they exist
    config_files = [init_config_path, survey_year_config_path]
    config_existence = [init_config_path.exists(), survey_year_config_path.exists()] 

    # Error evaluation and print message (if applicable)
    if not all(config_existence):
        missing_config = [ files for files, exists in zip( config_files, config_existence ) if not exists ]
        raise FileNotFoundError(f"The following configuration files do not exist: {missing_config}")

    ### Read configuration files
    # If configuration file existence is confirmed, proceed to reading in the actual files
    ## !!! TODO: Incorporate a configuration file validator that enforces required variables and formatting
    init_config_params = yaml.safe_load(init_config_path.read_text())
    survey_year_config_params = yaml.safe_load(survey_year_config_path.read_text())        

    # Validate that initialization and survey year configuration parameters do not intersect
    config_intersect = set(init_config_params.keys()).intersection(set(survey_year_config_params.keys()))
    
    # Error evaluation, if applicable
    if config_intersect:
        raise RuntimeError(
            f"The initialization and survey year configuration files comprise the following intersecting variables: {config_intersect}"
        )

    ### Format dictionary that will parameterize the `config` class attribute
    # Join the initialization and survey year parameters into a single dictionary
    config_to_add = { **init_config_params , **survey_year_config_params }
    
    # Amend length/age distribution locations within the configuration attribute
    config_to_add[ 'biometrics' ] = { 
        'bio_hake_len_bin': init_config_params[ 'bio_hake_len_bin' ] ,
        'bio_hake_age_bin': init_config_params[ 'bio_hake_age_bin' ]
    }
    
    del config_to_add['bio_hake_len_bin'] , config_to_add['bio_hake_age_bin']
    
    # Pass 'full_params' to the class instance
    return config_to_add

def validate_data_columns( file_name: Path ,
                           sheet_name: str ,
                           config_map: list ,
                           validation_settings: dict ):
    """
    Opens a virtual instance of each .xlsx file to validate the presence 
    of require data column/variable names

    Parameters
    ----------
    file_name: Path
        File path of data
    sheet_name: str
        Name of Excel sheet containing data
    config_map: list
        A list parsed from the file name that indicates how data attributes
        within `self` are organized
    validation_settings: dict
        The subset CONFIG_MAP settings that contain the target column names
    """
    
    # Open connection with the workbook and specific sheet 
    # This is useful for not calling the workbook into memory and allows for parsing
    # only the necessary rows/column names 
    try:
        workbook = load_workbook(file_name, read_only=True)

        # If multiple sheets, iterate through
        sheet_name = [ sheet_name ] if isinstance( sheet_name , str ) else sheet_name

        for sheets in sheet_name:
            sheet = workbook[ sheets ]
            
            # Validate that the expected columns are contained within the parsed
            # column names of the workbook   
            if 'vario_krig_para' in config_map:
                data_columns = [list(row) for row in zip(*sheet.iter_rows(values_only=True))][0]
            else:   
                data_columns = {col.value for col in sheet[1]}

            # Error evaluation and print message (if applicable)
            if not set(validation_settings.keys()).issubset(set(data_columns)):
                missing_columns = set(validation_settings.keys()) - set(data_columns)
                raise ValueError(f"Missing columns in the Excel file: {missing_columns}")
            
        # Close connection to the work book
        workbook.close()

    except Exception as e:
        print(f"Error reading file '{str(file_name)}': {e}")

def prepare_input_data( input_dict: dict , configuration_dict: dict ) :
    """
    Rearranges and organizes data formats of the initial file inputs

    Parameters
    ----------
    input_dict: dict
        Dictionary corresponding to the `input` attribute belonging to `Survey`-class 
    configuration_dict: dict
        Dictionary corresponding to the `config` attribute belonging to `Survey`-class 
    """

    # Generate length and age vectors
    # ---- Length vector
    length_bins = np.linspace( configuration_dict[ 'biometrics' ]['bio_hake_len_bin'][0] ,
                               configuration_dict[ 'biometrics' ]['bio_hake_len_bin'][1] ,
                               configuration_dict[ 'biometrics' ]['bio_hake_len_bin'][2] ,
                               dtype = np.float64 )
    # ---- Age vector
    age_bins = np.linspace( configuration_dict[ 'biometrics' ]['bio_hake_age_bin'][0] ,
                            configuration_dict[ 'biometrics' ]['bio_hake_age_bin'][1] , 
                            configuration_dict[ 'biometrics' ]['bio_hake_age_bin'][2] ,
                            dtype = np.float64 )   
    
    # Discretize these values into discrete intervals 
    # ---- Calculate binwidths
    # -------- Length
    length_binwidth = np.mean( np.diff( length_bins / 2.0 ) )
    # -------- Age
    age_binwidth = np.mean( np.diff( age_bins / 2.0 ) )
    # ---- Center the bins within the binwidths
    # -------- Length
    length_centered_bins = np.concatenate( ( [ length_bins[0] - length_binwidth ] ,
                                                length_bins + length_binwidth ) )
    # -------- Age
    age_centered_bins = np.concatenate( ( [ age_bins[0] - age_binwidth ] ,
                                            age_bins + age_binwidth ) ) 
    
    # Merge the vector and centered bins into dataframes that will be added into the `input` 
    # attribute
    # ---- Generate DataFrame for length
    length_bins_df = pd.DataFrame( { 'length_bins': length_bins } )
    # -------- Discretize the bins as categorical intervals
    length_bins_df[ 'length_intervals' ] = pd.cut( length_bins_df[ 'length_bins' ] ,
                                                   length_centered_bins )
    # ---- Generate DataFrame for age
    age_bins_df = pd.DataFrame( { 'age_bins': age_bins } )
    # -------- Discretize the bins as categorical intervals
    age_bins_df[ 'age_intervals' ] = pd.cut( age_bins_df[ 'age_bins' ] ,
                                             age_centered_bins )
    # ---- Update `input` attribute
    # -------- Length
    input_dict[ 'biology' ][ 'distributions' ][ 'length_bins_df' ] = length_bins_df
    # -------- Age 
    input_dict[ 'biology' ][ 'distributions' ][ 'age_bins_df' ] = age_bins_df
    # -------- Delete the duplicate configuration keys
    del configuration_dict[ 'biometrics' ]

    # Update `geo_strata` column names
    input_dict[ 'spatial' ][ 'geo_strata_df' ].rename( columns = { 'haul start': 'haul_start' ,
                                                                   'haul end': 'haul_end' } ,
                                                       inplace = True )
    
    # Create INPFC stratum key with correct latitude bins/intervals
    # ---- Rename stratum column name to avoid conflicts
    input_dict[ 'spatial' ][ 'inpfc_strata_df' ].rename( columns = { 'stratum_num': 'stratum_inpfc' } , 
                                                         inplace = True )
    # ---- Create latitude intervals to bin the strata
    latitude_bins = np.concatenate( [ [ -90 ] , 
                                      input_dict[ 'spatial' ][ 'inpfc_strata_df' ][ 'northlimit_latitude' ] ,
                                      [ 90 ] ] )
    # ---- Add categorical intervals
    input_dict[ 'spatial' ][ 'inpfc_strata_df' ][ 'latitude_interval' ] = (
        pd.cut( input_dict[ 'spatial' ][ 'inpfc_strata_df' ][ 'northlimit_latitude' ] * 0.99 , 
                latitude_bins )
    )

    # Merge haul numbers across biological variables
    # ---- Consolidate information linking haul-transect-stratum indices
    input_dict[ 'biology' ][ 'haul_to_transect_df' ] = (
        input_dict[ 'biology' ][ 'haul_to_transect_df' ]
        .merge( input_dict[ 'spatial' ][ 'strata_df' ] , on = 'haul_num' , how = 'outer' )
    )
    # ---- Create interval key for haul numbers to assign INPFC stratum
    haul_bins = np.sort( np.unique( np.concatenate( 
        [ input_dict[ 'spatial' ][ 'inpfc_strata_df' ][ 'haul_start' ] - int( 1 ) ,
          input_dict[ 'spatial' ][ 'inpfc_strata_df' ][ 'haul_end' ] ]
    ) ) )
    # ---- Quantize the INPFC dataframe hauls based on strata 
    input_dict[ 'spatial' ][ 'inpfc_strata_df' ][ 'haul_bin' ] = (
        pd.cut( ( input_dict[ 'spatial' ][ 'inpfc_strata_df' ][ 'haul_start' ]
                + input_dict[ 'spatial' ][ 'inpfc_strata_df' ][ 'haul_end' ] ) / 2 ,
                haul_bins )
    )
    # ---- Rename `stratum_num` column
    input_dict[ 'spatial' ][ 'inpfc_strata_df' ].rename( columns = { 'stratum_num': 'stratum_inpfc' } ,
                                                         inplace = True )
    # ---- Define haul bins with `haul_to_transect_df`
    input_dict[ 'biology' ][ 'haul_to_transect_df' ][ 'haul_bin' ] = (
        pd.cut( input_dict[ 'biology' ][ 'haul_to_transect_df' ][ 'haul_num' ] ,
                haul_bins )
    )
    # ---- Define INPFC stratum for `haul_to_transect_df`
    input_dict[ 'biology' ][ 'haul_to_transect_df' ] = (
        input_dict[ 'biology' ][ 'haul_to_transect_df' ]
        .merge( input_dict[ 'spatial' ][ 'inpfc_strata_df' ][ [ 'stratum_inpfc' , 'haul_bin' ] ] , how = 'left' )
        # .filter( regex = '^((?!_bin).)*$')
    )    
    # ---- Distribute this information to other biological variables
    # -------- Specimen
    input_dict[ 'biology' ][ 'specimen_df' ] = (
        input_dict[ 'biology' ][ 'specimen_df' ]
        .merge( input_dict[ 'biology' ][ 'haul_to_transect_df' ] , how = 'left' ) 
    )
    # -------- Length
    input_dict[ 'biology' ][ 'length_df' ] = (
        input_dict[ 'biology' ][ 'length_df' ] 
        .merge( input_dict[ 'biology' ][ 'haul_to_transect_df' ] , how = 'left' ) 
    )
    # -------- Catch
    input_dict[ 'biology' ][ 'catch_df' ] = (
        input_dict[ 'biology' ][ 'catch_df' ] 
        .merge( input_dict[ 'biology' ][ 'haul_to_transect_df' ] , how = 'left' ) 
    )

    # Relabel sex to literal words among biological data
    # ---- Specimen
    input_dict[ 'biology' ][ 'specimen_df' ][ 'sex' ] = (
        np.where( input_dict[ 'biology' ][ 'specimen_df' ][ 'sex' ] == int( 1 ) ,
                  'male' ,
                   np.where( input_dict[ 'biology' ][ 'specimen_df' ][ 'sex' ] == int( 2 ) ,
                             'female' , 'unsexed' ) )
    )
    # -------- Sex group
    input_dict[ 'biology' ][ 'specimen_df' ][ 'group_sex' ] = (
        np.where( input_dict[ 'biology' ][ 'specimen_df' ][ 'sex' ] != 'unsexed' , 
                  'sexed' , 'unsexed' )
    )
    # ---- Length
    input_dict[ 'biology' ][ 'length_df' ][ 'sex' ] = (
        np.where( input_dict[ 'biology' ][ 'length_df' ][ 'sex' ] == int( 1 ) ,
                  'male' ,
                   np.where( input_dict[ 'biology' ][ 'length_df' ][ 'sex' ] == int( 2 ) ,
                             'female' , 'unsexed' ) )
    )
    # -------- Sex group
    input_dict[ 'biology' ][ 'length_df' ][ 'group_sex' ] = (
        np.where( input_dict[ 'biology' ][ 'length_df' ][ 'sex' ] != 'unsexed' , 
                  'sexed' , 'unsexed' )
    )

    # Discretize the age and length bins of appropriate biological data
    # ---- Specimen
    input_dict[ 'biology' ][ 'specimen_df' ] = (
        input_dict[ 'biology' ][ 'specimen_df' ]
        .bin_variable( [ length_centered_bins , age_centered_bins ] , [ 'length' , 'age' ] )
    )
    # ---- Length
    input_dict[ 'biology' ][ 'length_df' ] = (
        input_dict[ 'biology' ][ 'length_df' ]
        .bin_variable( length_centered_bins , 'length' )
    )

    # Reorganize kriging/variogram parameters
    # ---- Kriging
    # -------- Generate dictionary comprising kriging model configuration
    kriging_params = (
        input_dict[ 'statistics' ][ 'kriging' ][ 'vario_krig_para_df' ]
        .filter( regex = 'krig[.]' )
        .rename( columns = lambda x: x.replace( 'krig.' , '' ) )
        .rename( columns = { 'ratio': 'anisotropy' ,
                                'srad': 'search_radius' } )
        .to_dict( orient = 'records' )[ 0 ]
    )
    # -------- Concatenate configuration settings for kriging
    kriging_params.update( configuration_dict[ 'kriging_parameters' ] )
    # ---- Variogram
    # -------- Generate dictionary comprising variogram model configuration
    variogram_params = (
        input_dict[ 'statistics' ][ 'kriging' ][ 'vario_krig_para_df' ]
        .filter( regex = 'vario[.]' )
        .rename( columns = lambda x: x.replace( 'vario.' , '' ) )
        .rename( columns = { 'lscl': 'correlation_range' ,
                                'powr': 'decay_power' ,
                                'hole': 'hole_effect_range' ,
                                'res': 'lag_resolution' ,
                                'nugt': 'nugget' } )
        .to_dict( orient = 'records' )[ 0 ]
    )
    # ---- Update the input attribute with the reorganized parameters
    input_dict[ 'statistics' ][ 'variogram' ].update( { 'model_config': variogram_params } )
    input_dict[ 'statistics' ][ 'kriging' ].update( { 'model_config': kriging_params } )  
    # -------- Delete the duplicate dataframe
    del input_dict[ 'statistics' ][ 'kriging' ][ 'vario_krig_para_df' ]
    # -------- Delete the duplicate configuration keys
    del configuration_dict[ 'kriging_parameters' ]
