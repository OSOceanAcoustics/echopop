from typing import List, Union, Optional
from pathlib import Path
import pandas as pd
import numpy as np
import copy
from .core import CONFIG_MAP, LAYER_NAME_MAP , DATA_STRUCTURE
from .utils.load import load_configuration , validate_data_columns , prepare_input_data

from .computation.acoustics import summarize_sigma_bs , nasc_to_biomass
from .computation.statistics import stratified_transect_statistic

from .spatial.transect import (
    prepare_transect_strata , 
    transect_distance , 
    summarize_transect_strata , 
    save_transect_coordinates
)

from .computation.biology import (
    filter_species , 
    fit_length_weight_relationship , 
    number_proportions , 
    quantize_number_counts , 
    fit_length_weights , 
    weight_proportions ,
    distribute_length_age , 
    partition_transect_age ,
    sum_strata_weight , 
    calculate_aged_unaged_proportions ,
    calculate_aged_biomass , 
    calculate_unaged_biomass ,
    apply_age_bins
)

from .utils import message as em
from .spatial.mesh import crop_mesh
from .spatial.projection import transform_geometry
from .spatial.krige import kriging

class Survey:
    """
    echopop base class that imports and prepares parameters for
    a survey. Additionally, it includes functions for accessing
    the modules associated with the transect and Kriging variable
    calculations, CV analysis, semi-variogram algorithm, and Kriging.

    Parameters
    ----------
    init_config_path : str or pathlib.Path
        A string specifying the path to the initialization YAML file
    survey_year_config_path : str or pathlib.Path
        A string specifying the path to the survey year YAML file

    Attributes
    ----------
    meta : dict 
        Metadata variable that provides summary information concerning the
        data contained within the class object (e.g. 'self.summary').
    config : dict 
        Configuration settings and parameters that can be referenced for
        various downstream and internal functions.
    data : dict
        Various dictionaries are incorporated into the Survey class object that 
        are directly referenced for various downstream and internal functions. This
        includes attributes such as 'biology', 'acoustics', and 'spatial' that represent
        various nested biological, acoustic, and spatial/stratification datasets imported
        based on the input files defined via the configuration settings.
    
    """
    def __init__(
        self,
        init_config_path: Union[str, Path] ,
        survey_year_config_path: Union[str, Path] ,
    ):
        # Loading the configuration settings and definitions that are used to 
        # initialize the Survey class object
        # ---- ATTRIBUTE ADDITIONS: `config`
        self.config = load_configuration( Path( init_config_path ) , Path( survey_year_config_path ) )

        # ---- Initialize the `input` attribute
        self.input = copy.deepcopy( DATA_STRUCTURE[ 'input' ] )

        # Loading the datasets defined in the configuration files
        self.load_survey_data()

    def load_survey_data( self ):
        """
        Loads the biological, NASC, and stratification
        data using parameters obtained from the configuration
        files. This will generate data attributes associated with the tags
        defined in both the configuration yml files and the reference CONFIG_MAP
        and LAYER_NAME_MAP dictionaries.
        """

        # Check whether data files defined from the configuration file exists
        # ---- Generate flat JSON table comprising all configuration parameter names
        flat_configuration_table = pd.json_normalize(self.config).filter(regex="filename")

        # ---- Parse the flattened configuration table to identify data file names and paths
        parsed_filenames = flat_configuration_table.values.flatten()

        # ---- Evaluate whether either file is missing
        data_existence = [(Path(self.config['data_root_dir']) / file).exists() for file in parsed_filenames]

        # Assign the existence status to each configuration file for error evaluation
        # ---- Error evaluation and print message (if applicable)
        if not all(data_existence):
            missing_data = parsed_filenames[ ~ np.array( data_existence ) ]
            raise FileNotFoundError(f"The following data files do not exist: {missing_data}")
        
        # Data validation and import
        # ---- Iterate through known datasets and datalayers
        for dataset in [*CONFIG_MAP.keys()]:

            for datalayer in [*self.config[dataset].keys()]:

                # Define validation settings from CONFIG_MAP
                validation_settings = CONFIG_MAP[dataset][datalayer]

                # Define configuration settings w/ file + sheet names
                config_settings = self.config[dataset][datalayer]

                # Create reference index of the dictionary path
                config_map = [dataset, datalayer]

                # Define the data layer name 
                # ---- Based on the lattermost portion of the file path string
                # Create list for parsing the hard-coded API dictionary
                if dataset == 'biological':
                    for region_id in [*self.config[dataset][datalayer].keys()]:

                        # Get file and sheet name
                        file_name = Path(self.config['data_root_dir']) / config_settings[region_id]['filename']
                        sheet_name = config_settings[region_id]['sheetname']
                        
                        # Update `config_map`
                        if len( config_map ) == 2:
                            config_map = config_map + [region_id]
                        else:
                            config_map[2] = region_id

                        # Validate column names of this iterated file
                        validate_data_columns( file_name , sheet_name , config_map , validation_settings )

                        # Validate datatypes within dataset and make appropriate changes to dtypes (if necessary)
                        # -- This first enforces the correct dtype for each imported column
                        # -- This then assigns the imported data to the correct class attribute
                        self.read_validated_data( file_name , sheet_name , config_map , validation_settings )
                else:
                    file_name = Path(self.config['data_root_dir']) / config_settings['filename']
                    sheet_name = config_settings['sheetname']

                    # If multiple sheets, iterate through
                    # ---- If multiple sheets, then this needs to be converted accordingly
                    sheet_name = [ sheet_name ] if isinstance( sheet_name , str ) else sheet_name

                    for sheets in sheet_name:
                        # Update if INPFC
                        if sheets.lower( ) == 'inpfc' :
                            # Update validation settings from CONFIG_MAP
                            validation_settings = CONFIG_MAP[dataset]['inpfc_strata']

                            # Update configuration key map
                            config_map = [dataset,  'inpfc_strata']

                        elif datalayer == 'geo_strata' :
                            # Update validation settings from CONFIG_MAP
                            validation_settings = CONFIG_MAP[dataset][datalayer]
                           
                            # Update configuration key map
                            config_map = [dataset,  datalayer]

                        # Validate datatypes within dataset and make appropriate changes to dtypes (if necessary)
                        # ---- This first enforces the correct dtype for each imported column
                        # ---- This then assigns the imported data to the correct class attribute
                        validate_data_columns( file_name , sheets , config_map , validation_settings )
                        
                        # Read in data and add to `Survey` object
                        self.read_validated_data( file_name , sheets , config_map , validation_settings ) 

        # Update the data format of various inputs within `Survey`
        self.input , self.config = prepare_input_data( self.input , self.config )

    def read_validated_data( self ,
                             file_name: Path ,
                             sheet_name: str ,
                             config_map: list ,
                             validation_settings: dict ):
        """
        Reads in data and validates the data type of each column/variable

        Parameters
        ----------
        file_name: Path
            The file name without the prepended file path
        sheet_name: str
            The Excel sheet name containing the target data
        config_map: list
            A list parsed from the file name that indicates how data attributes
            within `self` are organized
        validation_settings: dict
            The subset CONFIG_MAP settings that contain the target column names
        """
    
        # Based on the configuration settings, read the Excel files into memory. A format
        # exception is made for 'kriging.vario_krig_para' since it requires additional
        # data wrangling (i.e. transposing) to resemble the same dataframe format applied
        # to all other data attributes.
        if 'vario_krig_para' in config_map:
            # Read Excel file into memory and then transpose
            df_initial = pd.read_excel(file_name, header=None).T

            # Take the values from the first row and redfine them as the column headers
            df_initial.columns = df_initial.iloc[0]
            df_initial = df_initial.drop(0)

            # Slice only the columns that are relevant to the echopop module functionality
            valid_columns = list(set(validation_settings.keys()).intersection(set(df_initial.columns)))
            df_filtered = df_initial[valid_columns]

            # Ensure the order of columns in df_filtered matches df_initial
            df_filtered = df_filtered[df_initial.columns]

            # Apply data types from validation_settings to the filtered DataFrame
            df = df_filtered.apply(lambda col: col.astype(validation_settings.get(col.name, type(df_filtered.iloc[0][col.name]))))
        
        else:
            # Read Excel file into memory -- this only reads in the required columns
            df = pd.read_excel(file_name, sheet_name=sheet_name, usecols=validation_settings.keys())

            # Apply data types from validation_settings to the filtered DataFrame
            df = df.apply(lambda col: col.astype(validation_settings.get(col.name, type(col[0])))) 

        # Assign the data to their correct data attributes/keys
        if LAYER_NAME_MAP[config_map[0]]['superlayer'] == []:
            sub_attribute  = LAYER_NAME_MAP[config_map[0]]['name']
        else:
            sub_attribute = LAYER_NAME_MAP[config_map[0]]['superlayer'][0]
            
        # Step 2: Determine whether the dataframe already exists
        if sub_attribute in ['biology' , 'statistics' , 'spatial']:
            if sub_attribute == 'biology':
                # Add US / CAN as a region index 
                df['region'] = config_map[2] 

                # Apply CAN haul number offset 
                if config_map[2] == 'CAN':
                    df['haul_num'] += self.config['CAN_haul_offset']      

            # A single dataframe per entry is expected, so no other fancy operations are needed
            if sheet_name.lower() == 'inpfc':
                df_list = [ self.input[ sub_attribute ][ 'inpfc_strata_df' ] , df ]
                self.input[ sub_attribute ][ 'inpfc_strata_df' ] = pd.concat( df_list )
            else: 
                if config_map[ 0 ] == 'kriging' :
                    df_list = [ self.input[ sub_attribute ][ 'kriging' ][ config_map[1] + '_df' ] , df ]
                    self.input[ sub_attribute ][ 'kriging' ][ config_map[1] + '_df' ] = pd.concat(df_list)
                else :
                    df_list = [ self.input[ sub_attribute ][ config_map[1] + '_df' ] , df ]
                    self.input[ sub_attribute ][ config_map[1] + '_df' ] = pd.concat(df_list)
        elif sub_attribute == 'acoustics':
            
            # Toggle through including and excluding age-1
            if config_map[1] == 'no_age1':
                df = df.rename(columns={'NASC': 'NASC_no_age1'})
            else:
                df = df.rename(columns={'NASC': 'NASC_all_ages'})
            
            column_to_add = df.columns.difference(self.input['acoustics']['nasc_df'].columns).tolist()
            self.input['acoustics']['nasc_df'][column_to_add] = df[column_to_add]      
        else:
            raise ValueError('Unexpected data attribute structure. Check API settings located in the configuration YAML and core.py')              
            
    def transect_analysis(self ,
                          species_id: np.float64 = 22500 ,
                          exclude_age1: bool = True ,
                          stratum: str = 'ks' ,
                          verbose: bool = True ):
        """
        Calculate population-level metrics from acoustic transect measurements
        """    

        # Initialize the `analysis` data structure
        self.analysis = copy.deepcopy( DATA_STRUCTURE[ 'analysis' ] )

        # Initialize the `results` data structure
        self.results = copy.deepcopy( DATA_STRUCTURE[ 'results' ] )

        # Update settings to reflect the stratum definition
        self.analysis[ 'settings' ].update(
            {
                'transect': {
                    'stratum': stratum.lower( ) ,
                    'stratum_name': 'stratum_num' if stratum == 'ks' else 'inpfc' ,
                    'exclude_age1': exclude_age1
                }
            }
        )

        # Save transect coordinate information
        self.analysis[ 'transect' ].update(
            {
                'coordinates': save_transect_coordinates( self.input[ 'acoustics' ][ 'nasc_df' ] )
            }
        )

        # Filter out non-target species
        length_spp , specimen_spp , catch_spp = filter_species( [ self.input[ 'biology' ][ 'length_df' ] ,
                                                                  self.input[ 'biology' ][ 'specimen_df' ] ,
                                                                  self.input[ 'biology' ][ 'catch_df' ] ] ,
                                                                species_id )
        
        # Calculate mean sigma_bs per individual haul, KS stratum, and INPFC stratum
        self.analysis[ 'acoustics' ][ 'sigma_bs' ].update(
            summarize_sigma_bs( length_spp , 
                                specimen_spp , 
                                self.input[ 'spatial' ] , 
                                self.config ,
                                self.analysis[ 'settings' ] )
        )
        
        # Fit length-weight regression required for biomass calculation
        self.analysis[ 'biology' ][ 'weight' ].update(
            fit_length_weight_relationship( specimen_spp , 
                                            self.input[ 'biology' ][ 'distributions' ][ 'length_bins_df' ] )            
        )

        # Count the number of specimens across age and length bins
        self.analysis[ 'biology' ][ 'distributions' ].update(
            quantize_number_counts( specimen_spp , length_spp , stratum = stratum )
        )
        
        # Calculate the number proportions
        self.analysis[ 'biology' ][ 'proportions' ].update(
            {
                'number': number_proportions( self.analysis[ 'biology' ][ 'distributions' ] )
            }
        )

        # Calculate the average weights among male, female, and all fish across strata
        self.analysis[ 'biology' ][ 'weight' ].update(
            {
                'weight_stratum_df': fit_length_weights( self.analysis[ 'biology' ][ 'proportions' ][ 'number' ] ,
                                                         self.analysis[ 'biology' ][ 'weight' ] )
            }
        )

        # Calculate the weight proportions
        self.analysis[ 'biology' ][ 'proportions' ].update(
            {
                'weight': weight_proportions( specimen_spp , 
                                              length_spp , 
                                              catch_spp , 
                                              self.analysis['biology']['weight']['length_weight_regression']['weight_fitted_df'] )
            }
        )

        # Convert NASC into number density (animals/nmi^2), biomass density (kg/nmi^2), abundance
        # (# animals), and biomass (kg) for all fish, sexed (male/female) fish, and unsexed fish
        strata_adult_proportions , nasc_to_biology = (
            nasc_to_biomass( nasc_df = self.input[ 'acoustics' ][ 'nasc_df' ] ,
                             acoustics_dict = self.analysis[ 'acoustics' ] ,
                             distributions_dict = self.input[ 'biology' ][ 'distributions' ] ,
                             proportions_dict = self.analysis[ 'biology' ][ 'proportions' ] ,
                             TS_L_parameters = self.config[ 'TS_length_regression_parameters' ][ 'pacific_hake' ] ,
                             haul_hake_fractions_df = self.input[ 'spatial' ][ 'strata_df' ] ,
                             length_weight_strata = self.analysis[ 'biology' ][ 'weight' ]['weight_stratum_df'] ,
                             settings_dict = self.analysis[ 'settings' ] )
        ) 

        # Distribute abundance and biomass over length and age bins
        self.analysis[ 'biology' ][ 'population' ].update(
            {
                'tables': distribute_length_age( nasc_to_biology ,                                                 
                                                 self.analysis[ 'biology' ][ 'proportions' ] ,
                                                 self.analysis[ 'settings' ] )
            }
        )

        # Reapportion transect results to separate age-1 and age-2+ fish, generate age-1
        # abundance distributions for unaged fish, and generate biomass summary
        adult_transect , biomass_summary , unaged_age1_abundance = partition_transect_age( nasc_to_biology , 
                                                                                            self.analysis[ 'biology' ][ 'weight' ] ,
                                                                                            self.analysis[ 'settings' ] ,
                                                                                            self.analysis[ 'biology' ][ 'population' ] ,
                                                                                            strata_adult_proportions )
        # ---- Update analysis attributes
        # -------- Adult transect data
        self.analysis[ 'transect' ].update( { 'adult_transect_df' : adult_transect } )
        # -------- Unaged age-1 abundances
        self.analysis[ 'biology' ][ 'population' ][ 'tables' ][ 'abundance' ].update(
            {
                'unaged_age1_abundance_df': unaged_age1_abundance
            }
        )
        # ---- Update results (biomass summary)
        self.results[ 'transect' ].update( { 'biomass_summary_df': biomass_summary } )

        # Print result if `verbose == True`
        if verbose:
            em.transect_results_msg( biomass_summary )

    def stratified_summary( self ,
                            stratum: str = 'inpfc' ,
                            variable: str = 'biomass' ,
                            transect_sample: Optional[ float ] = None ,
                            transect_replicates: Optional[ float ] = None ,
                            verbose = True ):
        """
        Calculates the stratified summary statistics for biomass

        Notes
        -----
        This function calculates estimates and confidence intervals (95%) for biomass mean,
        variance, and coefficients of variation (CVs). This currently only calculates this
        metric for adult animals (age-2+) and is not calculated for other contrasts such as 
        age-class and sex. This also only applies to the transect results and is not currently 
        designed to be compatible with other derived population-level statistics (e.g. kriging).
        """ 
        
        # Error message for `stratum == 'ks'`
        if stratum == 'ks': 
            raise ValueError( """The Jolly and Hampton (1990) stratified analysis is not"""
        """ currently comaptible for calculating over KS strata. Please change `stratum` to"""
        """ 'inpfc'.""")

        # Parameterize analysis settings that will be applied to the stratified analysis
        self.analysis[ 'settings' ].update(
            {
                'stratified': {
                    'stratum': stratum.lower( ) ,
                    'stratum_name': 'stratum_num' if stratum == 'ks' else 'stratum_inpfc' ,
                    'transect_sample': (
                        self.config[ 'stratified_survey_mean_parameters' ][ 'strata_transect_proportion' ]
                        if transect_sample is None else transect_sample
                    ) ,
                    'transect_replicates': (
                        self.config[ 'stratified_survey_mean_parameters' ][ 'num_replicates' ]
                        if transect_replicates is None else transect_replicates
                    ) ,
                    'variable': variable ,
                    'exclude_age1': self.analysis[ 'settings' ][ 'transect' ][ 'exclude_age1' ]
                }
            }
        )

        # Pull in the updated age-2+ transect data        
        nasc_df = prepare_transect_strata( self.analysis[ 'transect' ] ,
                                           self.analysis[ 'settings' ][ 'stratified' ] )
        
        # Summarize transect spatial information
        transect_summary = transect_distance( nasc_df )

        # Summarize strata spatial information
        strata_summary = summarize_transect_strata( transect_summary )

        # Calculate the stratified mean, variance, and coefficient of variation
        replicates , stratified_results = stratified_transect_statistic( nasc_df ,
                                                                         transect_summary ,
                                                                         strata_summary ,
                                                                         self.analysis[ 'settings' ] )
        # ---- Add replicates to the analysis attribute 
        self.analysis[ 'stratified' ].update(
            { 'stratified_replicates_df': replicates }
        )
        # ---- Add the stratified results
        self.results[ 'stratified' ].update( stratified_results )

        # Print result if `verbose == True`
        if verbose:
            em.stratified_results_msg( stratified_results , 
                                       self.analysis[ 'settings' ][ 'stratified' ] )
            
    def krige( self ,
               coordinate_transform: bool = True ,
               extrapolate: bool = False ,
               kriging_parameters: Optional[ dict ] = None ,
               mesh_buffer_distance: float = 1.25 ,
               num_nearest_transects: int = 4 ,
               projection: Optional[ str ] = None ,
               stratum: str = 'ks' ,               
               variable: str = 'biomass_density' ,
               variogram_parameters: Optional[ dict ] = None ,
               verbose: bool = True ):
        """
        Interpolates biomass data using ordinary kriging
        
        Parameters
        ----------
        variable
            Biological variable that will be interpolated via kriging
        """  

        # Parameterize analysis settings that will be applied to the stratified analysis
        self.analysis[ 'settings' ].update(
            {
                'kriging': {
                    'exclude_age1': self.analysis[ 'settings' ][ 'transect' ][ 'exclude_age1' ] ,
                    'extrapolate': extrapolate ,
                    'kriging_parameters': (
                        self.input[ 'statistics' ][ 'kriging' ][ 'model_config' ]
                        if kriging_parameters is None else kriging_parameters
                    ) ,
                    'projection': (
                        self.config[ 'geospatial' ][ 'init' ] if projection is None else projection
                    ) ,
                    'standardize_coordinates': coordinate_transform ,
                    'stratum': stratum.lower( ) ,
                    'stratum_name': 'stratum_num' if stratum == 'ks' else 'inpfc' ,
                    'variable': variable ,
                    'variogram_parameters': (
                        self.input[ 'statistics' ][ 'variogram' ][ 'model_config' ]
                        if variogram_parameters is None else variogram_parameters
                    ) ,
                    'verbose': verbose
                }
            }
        )

        # Prepare temporary message concerning coordinate transformation = False
        if not coordinate_transform:
            raise ValueError( """Kriging without coordinate standardization is currently """
        """unavailable due to the kriging parameter `search_radius` being only defined for """
        """transformed x- and y-coordinates.""")

        # Pull in the biological transect data        
        nasc_df = prepare_transect_strata( self.analysis[ 'transect' ] ,
                                           self.analysis[ 'settings' ][ 'kriging' ] )
        
        # Crop the mesh grid if the kriged data will not be extrapolated
        if extrapolate is False:
            # ---- Update the analysis settings
            self.analysis[ 'settings' ][ 'kriging' ].update(
                {
                    'mesh_buffer_distance': mesh_buffer_distance ,
                    'num_nearest_transect': num_nearest_transects 
                }
            )
            # ---- Compute the cropped mesh
            mesh_full = crop_mesh( nasc_df ,
                                   self.input[ 'statistics' ][ 'kriging' ][ 'mesh_df' ] ,
                                   self.analysis[ 'settings' ][ 'kriging' ] )
            # ---- Print alert 
            print( "Kriging mesh cropped to prevent extrapolation beyond the defined "
            f"""`mesh_buffer_distance` value ({mesh_buffer_distance} nmi).""")
            
        else:
            # ---- Else, extract original mesh dataframe
            mesh_df = self.input[ 'statistics' ][ 'kriging' ][ 'mesh_df' ].copy( )
            # ---- Extract longitude column name
            mesh_longitude = (
                [ col for col in mesh_df.columns if 'lon' in col.lower( ) ][ 0 ]
            )  
            # ---- Latitude
            mesh_latitude = (
                [ col for col in mesh_df.columns if 'lat' in col.lower( ) ][ 0 ]
            )
            # ---- Rename the dataframe
            mesh_full = mesh_df.copy( ).rename( columns = { f"{mesh_longitude}": 'longitude' ,
                                                            f"{mesh_latitude}": 'latitude' } ) 

        # Standardize the x- and y-coordinates, if necessary
        if coordinate_transform:
            # ---- Transform transect data geometry (generate standardized x- and y-coordinates)
            nasc_df , d_x , d_y = (
                transform_geometry( nasc_df ,
                                    self.input[ 'statistics' ][ 'kriging' ][ 'isobath_200m_df' ] ,
                                    self.analysis[ 'settings' ][ 'kriging' ] )
            ) 
            # ---- Transform mesh grid geometry (generate standardized x- and y-coordinates)
            mesh_full , _ , _ = (
                transform_geometry( mesh_full ,
                                    self.input[ 'statistics' ][ 'kriging' ][ 'isobath_200m_df' ] ,
                                    self.analysis[ 'settings' ][ 'kriging' ] ,
                                    d_x , d_y )
            ) 
            # ---- Print alert 
            print( """Longitude and latitude coordinates (WGS84) converted to standardized """
            """coordinates (x and y).""")
        else:
            # ---- Else, duplicate the transect longitude and latitude coordinates as 'x' and 'y'
            # -------- x
            nasc_df[ 'x' ] = nasc_df[ 'longitude' ]
            # -------- y
            nasc_df[ 'y' ] = nasc_df[ 'latitude' ]
            # ---- Duplicate the mesh grid longitude and latitude coordinates as 'x' and 'y'
            # -------- x
            mesh_full[ 'x' ] = mesh_full[ 'longitude' ]
            # -------- y
            mesh_full[ 'y' ] = mesh_full[ 'latitude' ]
        # ---- Save to the analysis attribute
        self.analysis.update( { 'kriging': { 'mesh_df': mesh_full ,
                                             'transect_df': nasc_df } } )  
        
        # Run the kriging algorithm and procedure
        self.results[ 'kriging' ].update(
            {
                'interpolation': kriging( self.analysis[ 'kriging' ][ 'transect_df' ] ,
                                          self.analysis[ 'kriging' ][ 'mesh_df' ] ,
                                          self.analysis[ 'settings' ][ 'kriging' ] )
            }
        )

        # Print result if `verbose == True`
        if verbose:
            em.kriging_mesh_results_msg( self.results[ 'kriging' ] , 
                                         self.analysis[ 'settings' ][ 'kriging' ] )

    def apportion_kriged_biomass( self ,
                                  species_id ):
        """
        Apportion kriged biomass based on sex, age, and other proportions
        
        Notes
        -----
        This function apportions the kriged biomass values (interpolated along
        a user-defined mesh) based on sex and age. This provides biomass estimates
        for aged and unaged fish, total biomass inferred from the biological and 
        kriged bioamss estimates, the coefficient of variation within each cell, 
        and georeferenced biomass estimates for specific paired values of sex 
        and age-class. Binned length data (`length_df`) represent sexed, but
        unaged fish measured at Station 1. Specimen data (`specimen_data`) represen
        sexed and aged fish measured at Station 2. Total haul weights (`catch_df`) 
        represent the bulk weights of unaged and unsexed fish measured in Station 1.
        """  
        
        ### Import biological data and filter out non-target species
        # ---- Species filter
        length_spp , specimen_spp , haul_spp = filter_species( [ self.biology[ 'length_df' ] ,
                                                                 self.biology[ 'specimen_df' ] ,
                                                                 self.biology[ 'catch_df' ] ] ,
                                                               species_id )
        
        # ---- Remove 'bad' values 
        # ---- `specimen_spp`
        specimen_spp_filtered = specimen_spp[ specimen_spp.sex != 'unsexed' ].dropna( how = 'any' , subset = 'age' )

        # ---- `length_spp`
        length_spp_filtered = length_spp[ length_spp.sex != 'unsexed' ]

        ### Import discrete distribution bins
        # ---- Length
        length_intervals = self.biology[ 'distributions' ][ 'length' ][ 'length_interval_arr' ]
        
        # ---- Age
        age_intervals = self.biology[ 'distributions' ][ 'age' ][ 'age_interval_arr' ]

        ### Import length-weight regression parameters
        regression_parameters = self.statistics[ 'length_weight' ][ 'regression_parameters' ]
        regression_parameters = regression_parameters[ regression_parameters.sex != 'all' ]

        ### Import length-weight relationship
        length_weight_df = self.statistics[ 'length_weight' ][ 'length_weight_df' ]
        
        # ---- Construct the complete dataframe containing all possible
        # ---- strata, sexes, lengths, and ages
        # full_biological_indices = pd.DataFrame(
        #     list( product( np.unique( specimen_spp.stratum_num ) ,
        #                    np.unique( specimen_spp.species_id ) ,
        #                    np.unique( specimen_spp.sex ) ,
        #                    pd.cut( length_bins , length_intervals ) ,
        #                    pd.cut( age_bins , age_intervals ) ) ) ,
        #     columns = [ 'stratum_num' , 'species_id' , 'sex' , 
        #                 'length_bin' , 'age_bin' ]
        # )

        ### Process `haul_spp`
        # Remove haul numbers not found within `length_spp`
        haul_spp_matched = haul_spp[ haul_spp.haul_num.isin( length_spp.haul_num ) ]

        ### Sum weights for aged/unaged and all within each stratum
        weight_strata = sum_strata_weight( haul_spp_matched ,
                                           specimen_spp )
        
        ### Calculate the summed aged proportions for age-1+ (*_all) and age-2+ (*_adult) fish      
        aged_unaged_weight_proportions = calculate_aged_unaged_proportions( specimen_spp_filtered ,
                                                                            weight_strata )

        # !!! TODO: This does end up chewing up * a ton * of memory since the output dataframes are quite large
        aged_sex_biomass , aged_biomass = calculate_aged_biomass( self.statistics[ 'kriging' ][ 'kriged_biomass_df' ] ,
                                                                  specimen_spp_filtered ,
                                                                  length_intervals ,
                                                                  age_intervals ,
                                                                  aged_unaged_weight_proportions )

        ### Calculate unaged biomass for each sex and all animals
        unaged_sex_biomass = calculate_unaged_biomass( self.statistics[ 'kriging' ][ 'kriged_biomass_df' ] ,
                                                       length_spp_filtered ,
                                                       length_intervals ,
                                                       length_weight_df ,
                                                       regression_parameters ,
                                                       aged_unaged_weight_proportions )
        
        ### Re-distribute unaged biomass so it is compatible with aged biomass to calculate the overall summed biomass
        redistributed_unaged_sex_biomass , redistributed_unaged_biomass = apply_age_bins( aged_sex_biomass , 
                                                                                          unaged_sex_biomass )
        
        ### Sum the grand total by combining the aged and unaged biomass estimates post-apportionment
        # ---- Merge sexed
        overall_sexed_biomass = redistributed_unaged_sex_biomass.merge( aged_sex_biomass ,
                                                                        on = [ 'length_bin' , 'age_bin' , 'sex' , 'species_id' ] ,
                                                                        how = 'left' )
        
        # ---- Aggregate (sum)
        overall_sexed_biomass[ 'total_sexed_biomass_all' ] = (
            overall_sexed_biomass.biomass_sexed_unaged_all + 
            overall_sexed_biomass.biomass_sexed_aged_all
        )
        overall_sexed_biomass[ 'total_sexed_biomass_adult' ] = (
            overall_sexed_biomass.biomass_sexed_unaged_adult + 
            overall_sexed_biomass.biomass_sexed_aged_adult
        )
        
        # ---- Merge total
        overall_biomass = redistributed_unaged_biomass.merge( aged_biomass ,
                                                              on = [ 'length_bin' , 'age_bin' , 'species_id' ] ,
                                                              how = 'left' )
        
        # ---- Aggregate (sum)
        overall_biomass[ 'total_biomass_all' ] = (
            overall_biomass.biomass_unaged_all + 
            overall_biomass.biomass_aged_all
        )
        overall_biomass[ 'total_biomass_adult' ] = (
            overall_biomass.biomass_unaged_adult + 
            overall_biomass.biomass_aged_adult
        )

        ### Drop unnecessary columns
        # ---- Overall sexed biomass 
        overall_sexed_biomass.drop( [ 'summed_aged_biomass_all' ,
                                      'summed_aged_biomass_adult' ] ,
                                    axis = 1 ,
                                    inplace = True )

        ### Assign results to an attribute
        self.statistics[ 'kriging' ].update(
            {
                'apportioned_kriged_total_biomass_df': overall_biomass ,
                'apportioned_kriged_sexed_biomass_df': overall_sexed_biomass ,
            }
        )