from typing import List, Union
from pathlib import Path
import pandas as pd
import numpy as np
import copy
from .core import CONFIG_MAP, LAYER_NAME_MAP
### !!! TODO : This is a temporary import call -- this will need to be changed to 
# the correct relative structure (i.e. '.core' instead of 'EchoPro.core' at a future testing step)
from .computation.operations import bin_variable , bin_stats , count_variable , stretch
from .utils.data_file_validation import load_configuration , validate_data_columns
from .computation.acoustics import to_linear , ts_length_regression
from .computation.spatial import calculate_transect_distance
from .computation.statistics import stratified_transect_statistic

### !!! TODO : This is a temporary import call -- this will need to be changed to 
# the correct relative structure (i.e. '.utils.data_structure_utils' instead of 
# 'EchoPro.utils.data_structure_utils' at a future testing step)

class Survey:
    """
    EchoPro base class that imports and prepares parameters for
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
        ### Loading the configuration settings and definitions that are used to 
        # initialize the Survey class object
        # ATTRIBUTE ADDITIONS: `config`
        self.config = load_configuration( Path( init_config_path ) , Path( survey_year_config_path ) )

        # Initialize data attributes ! 
        self.acoustics = copy.deepcopy(LAYER_NAME_MAP['NASC']['data_tree'])
        self.biology = copy.deepcopy(LAYER_NAME_MAP['biological']['data_tree'])
        self.spatial = copy.deepcopy(LAYER_NAME_MAP['stratification']['data_tree'])
        self.statistics = copy.deepcopy(LAYER_NAME_MAP['kriging']['data_tree'])

        ### Loading the datasets defined in the configuration files
        # EXAMPLE ATTRIBUTE ADDITIONS: `biology`, `spatial`, `acoustics`
        self.load_survey_data()

        # Define length and age distributions
        self.biometric_distributions()

        ### !!! THIS IS TEMPORARY FOR DEBUGGING / TRACKING DATA ATTRIBUTE ASSIGNMENT 
        ### A utility function that helps to map the datasets currently present
        # within Survey object via the `___.summary` property. This also initializes
        # the `meta` attribute
        # ATTRIBUTE ADDITIONS: `meta`
        ##
        # self.populate_tree()

    def load_survey_data( self ):
        """
        Loads the biological, NASC, and stratification
        data using parameters obtained from the configuration
        files. This will generate data attributes associated with the tags
        defined in both the configuration yml files and the reference CONFIG_MAP
        and LAYER_NAME_MAP dictionaries.
        """

        ### Check whether data files defined from the configuration file exists
        # Generate flat JSON table comprising all configuration parameter names
        flat_configuration_table = pd.json_normalize(self.config).filter(regex="filename")

        # Parse the flattened configuration table to identify data file names and paths
        parsed_filenames = flat_configuration_table.values.flatten()

        # Evaluate whether either file is missing
        data_existence = [(Path(self.config['data_root_dir']) / file).exists() for file in parsed_filenames]

        # Assign the existence status to each configuration file for error evaluation
        # Error evaluation and print message (if applicable)
        if not all(data_existence):
            missing_data = parsed_filenames[ ~ np.array( data_existence ) ]
            raise FileNotFoundError(f"The following data files do not exist: {missing_data}")
        
        ### Data validation and import
        # Iterate through known datasets and datalayers
        for dataset in [*CONFIG_MAP.keys()]:

            for datalayer in [*self.config[dataset].keys()]:

                # Define validation settings from CONFIG_MAP
                validation_settings = CONFIG_MAP[dataset][datalayer]

                # Define configuration settings w/ file + sheet names
                config_settings = self.config[dataset][datalayer]

                # Create reference index of the dictionary path
                config_map = [dataset, datalayer]

                # Define the data layer name 
                # -- Based on the lattermost portion of the file path string
                # Create list for parsing the hard-coded API dictionary
                if dataset == 'biological':
                    for region_id in [*self.config[dataset][datalayer].keys()]:

                        # Get file and sheet name
                        file_name = Path(self.config['data_root_dir']) / config_settings[region_id]['filename']
                        sheet_name = config_settings[region_id]['sheetname']
                        config_map = config_map + ['region']
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

                    # Validate column names of this iterated file
                    validate_data_columns( file_name , sheet_name , config_map , validation_settings )

                    # Validate datatypes within dataset and make appropriate changes to dtypes (if necessary)
                    # -- This first enforces the correct dtype for each imported column
                    # -- This then assigns the imported data to the correct class attribute
                    
                    ## If multiple sheets, iterate through
                    # If multiple sheets, then this needs to be converted accordingly
                    sheet_name = [ sheet_name ] if isinstance( sheet_name , str ) else sheet_name

                    for sheets in sheet_name:
                        self.read_validated_data( file_name , sheets , config_map , validation_settings )     
        

        ### Merge haul numbers and regional indices across biological variables
        # Also add strata values/indices here alongside transect numbers 
        # -- Step 1: Consolidate information linking haul-transect-stratum
        self.biology['haul_to_transect_df'] = (
            self.biology['haul_to_transect_df']
            .merge(self.spatial['strata_df'], on ='haul_num' , how = 'outer' )
        )
        # -- Step 2: Distribute this information to the other biological variables
        # ---- Specimen 
        self.biology['specimen_df'] = (
            self.biology['specimen_df']
            .merge( self.biology['haul_to_transect_df'] , on = ['haul_num' , 'region' ] ) 
        )
        # ---- Length
        self.biology['length_df'] = (
            self.biology['length_df']
            .merge( self.biology['haul_to_transect_df'] , on = ['haul_num' , 'region'] )
        )
        # ---- Catch
        self.biology['catch_df'] = (
            self.biology['catch_df']
            .merge( self.biology['haul_to_transect_df'] , on = ['haul_num' , 'region'] )
        )

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
        # TODO : REVISIT THIS LATER
        if 'vario_krig_para' in config_map:
            # Read Excel file into memory and then transpose
            df_initial = pd.read_excel(file_name, header=None).T

            # Take the values from the first row and redfine them as the column headers
            df_initial.columns = df_initial.iloc[0]
            df_initial = df_initial.drop(0)

            # Slice only the columns that are relevant to the EchoPro module functionality
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

        # Assign the data to their correct data attributes
        # As of now this entails:
        # -- biology --> biology
        # -- stratification --> spatial
        # -- kriging --> statistics
        # -- NASC --> acoustics
        # Step 1: Step into the data attribute 
        if LAYER_NAME_MAP[config_map[0]]['superlayer'] == []:
            attribute_name  = LAYER_NAME_MAP[config_map[0]]['name']
            internal = getattr( self , attribute_name )
        else:
            attribute_name = LAYER_NAME_MAP[config_map[0]]['superlayer'][0]
            internal = getattr( self , attribute_name )
        # ------------------------------------------------------------------------------------------------
        # Step 2: Determine whether the dataframe already exists -- this only applies to some datasets
        # such as length that comprise multiple region indices (i.e. 'US', 'CAN')
        if attribute_name in ['biology' , 'statistics' , 'spatial']:
            if attribute_name == 'biology':
                # Add US / CAN as a region index 
                df['region'] = config_map[2] 

                # Apply CAN haul number offset 
                if config_map[2] == 'CAN':
                    df['haul_num'] += self.config['CAN_haul_offset']
            
            # If kriging dataset, then step one layer deeper into dictionary
            elif config_map[0] == 'kriging':
                internal = internal['kriging']    

            # A single dataframe per entry is expected, so no other fancy operations are needed
            # If geo_strata, differentiate between inpfc and stratification1 sheets
            ### TODO: Temporary approach for incorporating the inpfc dataset. An improved approach
            ### can be incorporated later on.
            if sheet_name.lower() == 'inpfc':
                df_list = [ internal[ 'inpfc_strata_df' ] , df ]
                internal[ 'inpfc_strata_df' ] = pd.concat( df_list )
            else: 
                df_list = [ internal[ config_map[1] + '_df' ] , df ]
                internal[ config_map[1] + '_df' ] = pd.concat(df_list)

        elif attribute_name == 'acoustics':
            
            # Step forward into 'acoustics' attribute
            internal = internal['nasc']

            # Toggle through including and excluding age-1
            # -- This is required for merging the NASC dataframes together
            if config_map[1] == 'no_age1':
                df = df.rename(columns={'NASC': 'NASC_no_age1'})
            else:
                df = df.rename(columns={'NASC': 'NASC_all_ages'})
            
            column_to_add = df.columns.difference(internal['nasc_df'].columns).tolist()
            internal['nasc_df'][column_to_add] = df[column_to_add]
        
        else:
            raise ValueError('Unexpected data attribute structure. Check API settings located in the configuration YAML and core.py')
        
    def biometric_distributions( self ):
        """
        Expand bin parameters into actual bins for length and age distributions
        """

        # Pull the relevant age and length bins and output a dictionary
        length_bins = np.linspace( self.config[ 'biometrics' ]['bio_hake_len_bin'][0] ,
                                   self.config[ 'biometrics' ]['bio_hake_len_bin'][1] ,
                                   self.config[ 'biometrics' ]['bio_hake_len_bin'][2] ,
                                   dtype = np.float64 )

        age_bins = np.linspace( self.config[ 'biometrics' ]['bio_hake_age_bin'][0] ,
                                self.config[ 'biometrics' ]['bio_hake_age_bin'][1] , 
                                self.config[ 'biometrics' ]['bio_hake_age_bin'][2] ,
                                dtype = np.float64 )

        ### Discretize the age and length arrays into user-defined bins that will be used later on
        ### to calculate various length- and age-weighted statistics
        # Determine bin widths
        length_binwidth = np.mean( np.diff( length_bins / 2.0 ) )
        age_binwidth = np.mean( np.diff( age_bins / 2.0 ) )

        # Now the bins are centered with the first and last elements properly appended
        # These create an along-array interval such that values can be cut/discretized into
        # bins that fall between each value/element of the array
        length_centered_bins = np.concatenate( ( [ length_bins[0] - length_binwidth ] ,
                                                length_bins + length_binwidth ) )
        age_centered_bins = np.concatenate( ( [ age_bins[0] - age_binwidth ] ,
                                            age_bins + age_binwidth ) ) 

        # Add to the biological data attribute so it can be accessed downstream
        self.biology['distributions'] = {
            'length': {
                'length_bins_arr': length_bins ,
                'length_interval_arr': length_centered_bins ,
            } ,
            'age': {
                'age_bins_arr': age_bins ,
                'age_interval_arr': age_centered_bins ,
            } ,
        }

    def transect_analysis(self ,
                          species_id: np.float64 = 22500 ):
    #     # INPUTS
    #     # This is where the users can designate specific transect numbers,
    #     # stratum numbers, species, etc. These would be applied to the functions
    #     # below
        
    #     # Initialize new attribute
    #     self.results = {}
    #     #### TODO: THIS SHOULD BE ADDED TO THE ORIGINAL SURVEY OBJECT CREATION
    #     #### THIS IS INCLUDED HERE FOR NOW FOR TESTING PURPOSES -- ALL CAPS IS CRUISE
    #     #### CONTROL FOR "REMEMBER TO MAKE THIS CHANGE BRANDYN !!!"
    
        # Initialize major data structures that will be added (**tentative names**)
        self.acoustics['sigma_bs'] = {}
        self.biology['weight'] = {}
        self.biology['population'] = {}
        self.statistics['length_weight'] = {}
                
        # Calculate sigma_bs per stratum 
        ### This will also provide dataframes for the length-binned, mean haul, and mean strata sigma_bs       
        self.strata_mean_sigma_bs( species_id )        
        
        # Fit length-weight regression required for biomass calculation
        self.fit_binned_length_weight_relationship( species_id )

        # Calculate the average sex-distributed weight and proportion per stratum
        self.strata_sex_weight_proportions( species_id )

        ### TODO : REMOVE AGE-0 !! -- OR separate into a separete bin
        ### Keep in mind -- NASC exports are age-2+
        # Calculate the age-binned weight per sex per stratum when both considering and ignoring age-0 and age-1 fish
        self.strata_age_binned_weight_proportions( species_id )
        
        # Synthesize all of the above steps to begin the conversion from 
        # integrated acoustic backscatter (ala NASC) to estimates of biological
        # relevance 
        self.nasc_to_biomass_conversion( species_id )
        
    def strata_mean_sigma_bs( self ,
                              species_id: np.float64 ):
        """
        Calculates the stratified mean sigma_bs for each stratum

        Parameters
        ----------
        species_id : np.float64
            Numeric code representing a particular species of interest

        Notes
        -----
        This function iterates through each stratum to fit acoustic target 
        strength (TS, dB re. 1 m^2) values based on length distributions recorded for each
        stratum. These fitted values are then convereted from the logarithmic
        to linear domain (sigma_bs, m^2) and subsequently averaged. These are required for 
        later functions that will convert vertically integrated backscatter (e.g. the nautical
        area scattering coefficient, or NASC, m^2 nmi^-2) to estimates of areal density (animals nmi^-2).
        """
        
        # Reformat 'specimen_df' to match the same format as 'len_df'
        ### First make copies of each
        specimen_df_copy = self.biology['specimen_df'].copy()
        specimen_df_copy = specimen_df_copy[ specimen_df_copy.species_id == species_id ]
        length_df_copy = self.biology['length_df'].copy()
        length_df_copy = length_df_copy[ length_df_copy.species_id == species_id ]
        
        ### Iterate through 'specimen_df_copy' to grab 'length' and the number of values in that bin
        ### Indexed by 'haul_num' , 'stratum_num' , 'species_id' , 'region' , 'length'
        spec_df_reframed = (
            specimen_df_copy
            .groupby(['haul_num', 'stratum_num' , 'species_id', 'length'])
            .apply(lambda x: len(x['length']))
            .reset_index(name= 'length_count' )
            )
        
        ### Concatenate the two dataframes
        all_length_df = pd.concat( [ spec_df_reframed , length_df_copy ] , join = 'inner' )
        
        # Import parameters from configuration
        ts_length_parameters = self.config[ 'TS_length_regression_parameters' ]['pacific_hake']
        slope = ts_length_parameters[ 'TS_L_slope' ]
        intercept = ts_length_parameters[ 'TS_L_intercept' ]
        
        # Convert length values into TS
        ### ??? TODO: Not necessary for this operation, but may be useful to store for future use ?
        ### ??? TODO: May need functions later on that estimate TS based on length using other methods,
        ### ??? TODO: so that would need to be tested/triaged at this step of the code
        all_length_df[ 'TS' ] = ts_length_regression( all_length_df[ 'length' ] , slope , intercept )
        
        # Convert TS into sigma_bs
        all_length_df[ 'sigma_bs' ] = to_linear( all_length_df[ 'TS' ] )
        
        # Calculate the weighted mean sigma_bs per haul
        ### This will track both the mean sigma_bs and sample size since this will propagate as a
        ### grouped mean contained with a shared stratum
        mean_haul_sigma_bs = (
            all_length_df
            .groupby(['haul_num' , 'stratum_num' , 'species_id' ])[['sigma_bs' , 'length_count']]
            .apply(lambda x: np.average( x[ 'sigma_bs' ] , weights=x[ 'length_count' ]))
            .to_frame( 'sigma_bs_mean' )
            .reset_index()
        )
                
        # Now these values can be re-merged with stratum information and averaged over strata
        mean_strata_sigma_bs = (
            mean_haul_sigma_bs
            .groupby(['stratum_num' , 'species_id'])[ 'sigma_bs_mean' ]
            .mean()
            .reset_index()
        )
        
        # Add back into object
        self.acoustics['sigma_bs'] = {
            'length_binned': all_length_df ,
            'haul_mean': mean_haul_sigma_bs ,
            'strata_mean': mean_strata_sigma_bs
        }
        
        # Fill in missing sigma_bs values
        self.impute_missing_sigma_bs( species_id )
    
    def impute_missing_sigma_bs( self ,
                                 species_id: np.float64 ):
        """
        Imputes sigma_bs for strata without measurements or values

        Parameters
        ----------
        species_id : np.float64
            Numeric code representing a particular species of interest

        Notes
        -----
        This function iterates through all stratum layers to impute either the
        nearest neighbor or mean sigma_bs for strata that are missing values.
        """    
        #### TODO: CURRENTLY : species_id is unused since only hake are being processed, but this will need
        ### to actually be parameterized in the future
        # Collect all possible strata values
        strata_options = np.unique( self.spatial[ 'strata_df' ].copy().stratum_num )
        
        #
        strata_mean = self.acoustics[ 'sigma_bs' ][ 'strata_mean' ].copy()
        
        # impute missing strata values
        present_strata = np.unique(strata_mean[ 'stratum_num' ]).astype(int)
        missing_strata = strata_options[~(np.isin(strata_options, present_strata))]
        
        if len(missing_strata) > 0:
            
            # Concatenate the existing data with a DataFrame including the missing strata 
            # with NaN placeholders for 'mean_sigma_bs'            
            sigma_bs_impute = (
                pd.concat( [ strata_mean , 
                             pd.DataFrame( {
                                 'stratum_num': missing_strata , 
                                 'species_id': np.repeat( np.unique( strata_mean.species_id ) ,
                                                         len( missing_strata ) ) ,
                                 'sigma_bs_mean': np.repeat( np.nan ,
                                                             len( missing_strata ) )
                             } ) ] )
                .sort_values( 'stratum_num' )        
            )
            
            # Find strata intervals to impute over        
            for i in missing_strata:
                strata_floor = present_strata[present_strata < i]
                strata_ceil = present_strata[present_strata > i]

                new_stratum_below = np.max(strata_floor) if strata_floor.size > 0 else None
                new_stratum_above = np.min(strata_ceil) if strata_ceil.size > 0 else None      
                
                sigma_bs_indexed = sigma_bs_impute[sigma_bs_impute['stratum_num'].isin([new_stratum_below, new_stratum_above])]
                
                sigma_bs_impute.loc[sigma_bs_impute.stratum_num==i , 'sigma_bs_mean' ] = sigma_bs_indexed[ 'sigma_bs_mean' ].mean()
                
            self.acoustics[ 'sigma_bs' ][ 'strata_mean' ] = sigma_bs_impute        
   
    
    def fit_binned_length_weight_relationship( self ,
                                               species_id: np.float64 ):
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
        
        ### First make copies of each
        specimen_df_spp = self.biology['specimen_df'].copy().pipe( lambda df: df.loc[ df.species_id == species_id ] )
        length_weight_df = specimen_df_spp[['length', 'weight']].dropna(how='any')
        
        # pull distribution values
        length_bins = self.biology['distributions']['length']['length_bins_arr']
        length_intervals = self.biology['distributions']['length']['length_interval_arr']
        
        # length-weight regression
        [ rate , initial ] = np.polyfit(np.log10(length_weight_df['length']), np.log10(length_weight_df['weight']), 1)
        
        # predict weight (fit equation)
        fitted_weight = 10 ** initial * length_bins ** rate
               
        # add to object        
        self.statistics['length_weight']['regression'] = {
            'rate': rate ,
            'initial': initial ,
        }
        
        # summarize mean/samples of length and weight values (binned statistics)
        length_bin_stats = length_weight_df.bin_stats( bin_variable = 'length' ,
                                                       bin_values = length_intervals )
        
        # fill bins where n_length < 5 w/ regressed weight values
        # TODO : the `bin_stats` function defaults to prepending 'mean' and 'n' -- this will
        # TODO : need to be changed to comport with the same name formatting as the rest of 
        # TODO : the module
        length_bin_stats[ 'weight_modeled' ] = length_bin_stats[ 'mean_weight' ].copy()
        low_n_indices = np.where(length_bin_stats[ 'n_weight' ].values < 5)[0].copy()   
        length_bin_stats.loc[low_n_indices, 'weight_modeled'] = fitted_weight[low_n_indices].copy()
            
        self.statistics[ 'length_weight' ][ 'length_weight_df' ] = length_bin_stats
        
    def strata_sex_weight_proportions( self ,
                                       species_id: np.float64 ):
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
        
        ### Reformat 'specimen_df' to match the same format as 'len_df'
        # First make copies of each
        specimen_df_copy = self.biology['specimen_df'].copy().pipe( lambda df: df.loc[ df.species_id == species_id ] )
        length_df_copy = self.biology['length_df'].copy().pipe( lambda df: df.loc[ df.species_id == species_id ] )
        
        ### Pull length distribution values
        length_intervals = self.biology['distributions']['length']['length_interval_arr']   
        
        ### Calculate the sex proportions/frequencies for station 1 (length_df) across all strata        
        length_grouped = (   
            length_df_copy
            .bin_variable( length_intervals , 'length' ) # appends `length_bin` column
            .assign( group = lambda x: np.where( x['sex'] == int(1) , 'male' , 'female' ) ) # assigns str variable for comprehension
            .pipe( lambda df: pd.concat( [ df.loc[ df[ 'sex' ] != 3 ] , df.assign( group = 'all' ) ] ) ) # appends male-female to an 'all' dataframe
            .assign( station = 1 ) # assign station number for later functions
            )
        
        ### Calculate the sex proportions/frequencies for station 2 (specimen_df) across all strata
        specimen_grouped = (
            specimen_df_copy
            .bin_variable( length_intervals , 'length' ) # appends `length_bin` column
            .assign( group = lambda x: np.where( x['sex'] == int(1) , 'male' , 'female' ) ) # assigns str variable for comprehension
            .pipe( lambda df: pd.concat( [ df.loc[ df[ 'sex' ] != 3 ] , df.assign( group = 'all' ) ] ) ) # appends male-female to an 'all' dataframe
            .assign( station = 2 ) # assign station number for later functions
        )
        
        ### "Meld" the two datasets together to keep downstream code tidy
        station_sex_length = (
            specimen_grouped # begin reformatting specimen_grouped so it resembles same structure as length_grouped
            .meld( length_grouped ) # "meld": reformats specimen_grouped and then concatenates length_grouped 
            # sum up all of the length values with respect to each length bin
            .count_variable( contrasts = [ 'group' , 'station' , 'stratum_num' , 'length_bin' ] , # grouping variables
                             variable = 'length_count' , # target value to apply a function
                             fun = 'sum' ) # function to apply
        )
        
        ### Calculate total sample size ('group' == 'all') that will convert 'counts' into 'frequency/proportion'
        total_n = (
            station_sex_length 
            .loc[ station_sex_length.group.isin( [ 'all' ] ) ] # filter out to get to just 'all': this is where sex = [ 1 , 2 , 3]
            .groupby( [ 'stratum_num' ] )[ 'count' ] # group by each stratum with 'count' being the target variable
            .sum( ) # sum up counts per stratum
            .reset_index( name = 'n_total' ) # rename index
        )
        
        ### Convert counts in station_sex_length into frequency/proportion
        station_length_aggregate = (
            station_sex_length
            # calculate the within-sample sum and proportions (necessary for the downstream dot product calculation)
            .pipe( lambda x: x.assign( within_station_n = x.groupby( [ 'group' , 'station' , 'stratum_num' ] )[ 'count' ].transform( sum ) ,
                                       within_station_p = lambda x: x[ 'count' ] / x[ 'within_station_n' ] ) )
            .replace( np.nan, 0 ) # remove erroneous NaN (divide by 0 or invalid values)
            .merge( total_n , on = 'stratum_num' ) # merge station_sex_length with total_n
            # proportion of each count indexed by each length_bin relative to the stratum total
            .assign( overall_length_p = lambda x: x[ 'count' ] / x[ 'n_total' ] ,
                     overall_station_p = lambda x: x[ 'within_station_n' ] / x[ 'n_total' ] )
            .replace( np.nan, 0 ) # remove erroneous NaN (divide by 0 or invalid values)
        )
        
        ### Calculate the sex distribution across strata
        sex_proportions = (
            station_length_aggregate
            .loc[ station_length_aggregate.group.isin( ['male' , 'female'] ) ] # only parse 'male' and 'female'
            # create a pivot that will reorient data to the desired shape
            .pivot_table( index = [ 'group' , 'station' ] , 
                          columns = [ 'stratum_num' ] , 
                          values = [ 'overall_station_p' ] )
            .groupby( 'group' )
            .sum()
        )
        
        ### Calculate the proportions each dataset / station contributed within each stratum
        station_proportions = (
            station_length_aggregate
            .loc[ station_length_aggregate.group.isin( [ 'all' ] ) ] # only parse 'all'
            # create a pivot that will reorient data to the desired shape
            .pivot_table( index = [ 'group' , 'station' ] , 
                          columns = 'stratum_num' , 
                          values = 'overall_station_p' )
            .groupby( 'station' )
            .sum()
        )
        
        ### Calculate the sex distribution within each station across strata
        sex_station_proportions = (
            station_length_aggregate
            .loc[ station_length_aggregate.group.isin( ['male' , 'female'] ) ] # only parse 'male' and 'female'
            # create a pivot that will reorient data to the desired shape
            .pivot_table( index = [ 'group' , 'station' ] , 
                          columns = 'stratum_num' , 
                          values = 'overall_station_p' )
            .groupby( [ 'group' , 'station' ] )
            .sum()
        )
        
        ### Merge the station and sex-station indexed proportions to calculate the combined dataset mean fractions
        sex_stn_prop_merged = (
            sex_station_proportions
            .stack( )
            .reset_index( name = 'sex_stn_p')
            .merge( station_proportions
                .stack()
                .reset_index( name = 'stn_p' ) , on = [ 'stratum_num' , 'station' ] )
            .pivot_table( columns = 'stratum_num' ,
                          index = [ 'station' , 'group' ] ,
                          values = [ 'stn_p' , 'sex_stn_p' ] )    
        )
        
        ### Format the length bin proportions so they resemble a similar table/matrix shape as the above metrics
        # Indexed / organized by sex , station, and stratum
        length_proportion_table = (
            station_length_aggregate
            .pivot_table( columns = [ 'group' , 'station' , 'stratum_num' ] , 
                          index = [ 'length_bin' ] ,
                          values = [ 'within_station_p' ] )[ 'within_station_p' ]
        )
        
        ### Calculate combined station fraction means
        # Station 1 
        stn_1_fraction = ( sex_stn_prop_merged.loc[ 1 , ( 'stn_p' ) ]                           
                          / ( sex_stn_prop_merged.loc[ 1 , ( 'stn_p' ) ] 
                            + sex_stn_prop_merged.loc[ 2 , ( 'sex_stn_p' ) ] ) )
        
        # Station 2
        stn_2_fraction = ( sex_stn_prop_merged.loc[ 2 , ( 'sex_stn_p' ) ]                          
                          / ( stn_1_fraction
                            + sex_stn_prop_merged.loc[ 2 , ( 'sex_stn_p' ) ] ) )
        
        ### Calculate the average weight across all animals, males, and females
        # Pull fitted weight values
        fitted_weight = self.statistics['length_weight']['length_weight_df']
        
        # Total
        total_weighted_values = ( length_proportion_table.loc[ : , ( 'all' , 1 ) ] * station_proportions.loc[ 1 , ] +
                                  length_proportion_table.loc[ : , ( 'all' , 2 ) ] * station_proportions.loc[ 2 , ] )
        total_weight = fitted_weight[ 'weight_modeled' ].dot( total_weighted_values.reset_index( drop = True ) )
        
        # Male
        male_weighted_values = ( length_proportion_table.loc[ : , ( 'male' , 1 ) ] * stn_1_fraction.loc[ 'male' ] +
                                 length_proportion_table.loc[ : , ( 'male' , 2 ) ] * stn_2_fraction.loc[ 'male' ] )
        male_weight = fitted_weight[ 'weight_modeled' ].dot( male_weighted_values.reset_index( drop = True ) )
        
        # Female
        female_weighted_values = ( length_proportion_table.loc[ : , ( 'female' , 1 ) ] * stn_1_fraction.loc[ 'female' ] +
                                   length_proportion_table.loc[ : , ( 'female' , 2 ) ] * stn_2_fraction.loc[ 'female' ] )
        female_weight = fitted_weight[ 'weight_modeled' ].dot( female_weighted_values.reset_index( drop = True ) )
        
        ### Store the data frame in an accessible location
        self.biology[ 'weight' ][ 'weight_strata_df' ] = pd.DataFrame( {
            'stratum_num': total_weight.index ,
            'proportion_female': sex_proportions.loc[ 'female' , : ][ 'overall_station_p' ].reset_index(drop=True) ,
            'proportion_male': sex_proportions.loc[ 'male' , : ][ 'overall_station_p' ].reset_index(drop=True) ,
            'average_weight_female': female_weight ,            
            'average_weight_male': male_weight ,
            'average_weight_total': total_weight
        } )

    #!!! TODO : Provide argument that will exclude age-0 and age-1 fish when flagged
    def strata_age_binned_weight_proportions( self , 
                                              species_id: np.float64 ):
        """
        Calculates the age- and sex-binned proportions across all strata
        with respect to specimen counts and weights

        Parameters
        ----------
        species_id : np.float64
            Numeric code representing a particular species of interest

        Notes
        -----
        The sex-stratified proportions for both counts and weights currently 
        incorporate age-0 and age-1 fish without an option to exclude. 
        
        """ 

        ### Reformat 'specimen_df' to match the same format as 'len_df'
        # First make copies of each
        specimen_df_copy = self.biology['specimen_df'].copy().pipe( lambda df: df.loc[ df.species_id == species_id ] )

        # Import length bins
        length_intervals = self.biology[ 'distributions' ][ 'length' ][ 'length_interval_arr' ]

        ### Calculate age bin proportions when explicitly excluding age-0 and age-1 fish
        # Calculate age proportions across all strata and age-bins
        age_proportions = (
            specimen_df_copy
            .dropna( how = 'any' )
            .count_variable( variable = 'length' ,
                             contrasts = [ 'stratum_num' , 'age' ] ,
                             fun = 'size' )
            .pipe( lambda x: x.assign( stratum_count = x.groupby( [ 'stratum_num' ] )[ 'count' ].transform( sum ) ,
                                       stratum_proportion = lambda x: x[ 'count' ] / x[ 'stratum_count' ] ) )               
        )

        # Calculate adult proportions/contributions (in terms of summed presence) for each stratum
        age_2_and_over_proportions = (
            age_proportions
            .pipe( lambda df: df
                .groupby( 'stratum_num' )
                .apply( lambda x: 1 - x.loc[ x[ 'age' ] <= 1 ][ 'count' ].sum() / x[ 'count' ].sum( ) ) ) 
            .reset_index( name = 'number_proportion' )
        )

        ### Calculate proportional contributions of each age-bin to the summed weight of each stratum
        ### when explicitly excluding age-0 and age-1 fish
        # Calculate the weight proportions across all strata and age-bins
        age_weight_proportions = (
            specimen_df_copy
            .dropna( how = 'any' )
            .pipe( lambda x: x.assign( weight_age = x.groupby( [ 'stratum_num' , 'age' ] )[ 'weight' ].transform( sum ) ,
                                    weight_stratum = x.groupby( [ 'stratum_num' ] )[ 'weight' ].transform( sum ) ) )
            .groupby( [ 'stratum_num' , 'age' ] )
            .apply( lambda x: x[ 'weight_age' ].sum( ) / x[ 'weight_stratum' ].sum ( ) )
            .reset_index( name = 'weight_stratum_proportion' )
        )
        
        # Calculate adult proportions/contributions (in terms of summed weight) for each stratum
        age_2_and_over_weight_proportions = (
            age_weight_proportions
            .pipe( lambda df: df
                .groupby( 'stratum_num' )
                .apply( lambda x: 1 - x.loc[ x[ 'age' ] <= 1 ][ 'weight_stratum_proportion' ].sum()) ) 
                .reset_index( name = 'weight_proportion' ) 
        )

        ### Now this process will be repeated but adding "sex" as an additional contrast and considering
        ### all age-bins
        age_sex_weight_proportions = (
            specimen_df_copy
            .dropna( how = 'any' )
            .bin_variable( bin_values = length_intervals ,
                        bin_variable = 'length' )
            .count_variable( contrasts = [ 'stratum_num' , 'age' , 'length_bin' , 'sex' ] ,
                            variable = 'weight' ,
                            fun = 'sum' )
            .pipe( lambda df: df.assign( total = df.groupby( [ 'stratum_num' , 'sex' ] )[ 'count' ].transform( sum ) ) )
            .groupby( [ 'stratum_num' , 'age' , 'sex' ] )
            .apply( lambda x: ( x[ 'count' ] / x[ 'total' ] ).sum() ) 
            .reset_index( name = 'weight_sex_stratum_proportion' )
            .assign( sex = lambda x: np.where( x[ 'sex' ] == 1 , 'male' , np.where( x[ 'sex' ] == 2 , 'female' , 'unsexed' ) ) )
        )

        ### Add these dataframes to the appropriate data attribute
        self.biology[ 'weight' ].update( {
            'sex_stratified': {
                'weight_proportions': age_sex_weight_proportions ,
            } ,
            'age_stratified': {
                'age_1_excluded': {
                    'number_proportions': age_2_and_over_proportions ,
                    'weight_proportions': age_2_and_over_weight_proportions ,
                } ,
                'age_1_included': {
                    'number_proportions': age_proportions ,
                    'weight_proportions': age_weight_proportions
                }
            }
        } )

    #!!! TODO : Provide argument that will exclude age-0 and age-1 fish when flagged
    def nasc_to_biomass_conversion( self , 
                                    species_id: np.float64  ):
        """
        Converts integrated acoustic backscatter (NASC) into estimates of 
        areal number/biomass densities, total abundance, and total biomass

        Parameters
        ----------
        species_id : np.float64
            Numeric code representing a particular species of interest

        Notes
        -----
        This function converts NASC into estimates of population-level metrics 
        (abundance, biomass, areal densities) stratified by transects, sex, 
        length-based strata, and age.
        """ 
        ### Calculate sex-indexed weight proportions for each stratum
        sex_indexed_weight_proportions = index_sex_weight_proportions( copy.deepcopy( self.biology ) )

        ### Join acoustic and biological dataframes that incorporate the fractions of integrated 
        ### acoustic backscatter specific to target organisms with adult-specific number and weight
        ### proportions needed to estimate the convers from NASC to biomass
        nasc_fraction_total_df = index_transect_age_sex_proportions( copy.deepcopy( self.acoustics ) ,
                                                                     copy.deepcopy( self.biology ) , 
                                                                     self.spatial[ 'strata_df' ].copy() )

        ### Calculate the areal number densities (rho_a)
        # rho_a = animals / nmi^-2
        nasc_to_areal_number_density_df = (
            nasc_fraction_total_df
            # Calculate areal number density (rho_a) for the total sample
            .assign( rho_a_total = lambda x: np.round(
                x[ 'fraction_hake' ] * x[ 'NASC_no_age1' ] / ( 4 * np.pi * x[ 'sigma_bs_mean' ] ) ) )
            # Apportion rho_a based on sex and general age categorization (adult versus not adult)
            .pipe( lambda df: 
                   df.assign( 
                       rho_a_male = np.round( df[ 'rho_a_total' ] * df[ 'proportion_male' ] ) ,
                       rho_a_female = np.round( df[ 'rho_a_total' ] * df[ 'proportion_female' ] ) )
                   .assign( rho_a_unsexed = lambda x: np.round( x[ 'rho_a_total' ] - x[ 'rho_a_male' ] - x[ 'rho_a_female' ] ) ) )             
        )

        # Create dataframe to save it to Survey object
        areal_number_density_df = (
           nasc_to_areal_number_density_df
           .stretch( variable = 'rho_a' )
           .merge( nasc_adult_number_proportions , on = [ 'stratum_num' ] )
           .assign( rho_a_adult = lambda x: x[ 'rho_a' ] * x[ 'adult_number_proportion' ] )
        )

        ### Calculate the areal biomass densitiies (B_a)
        # B_a = kg / nmi^-2
        nasc_to_areal_biomass_density_df = (
            nasc_to_areal_number_density_df
            # Convert rho_a into biomass densities
            .assign( B_a_total = lambda x: x[ 'rho_a_total' ] * x[ 'average_weight_total' ] ,
                     B_a_male = lambda x: x[ 'rho_a_male' ] * x[ 'average_weight_male' ] ,
                     B_a_female = lambda x: x[ 'rho_a_female' ] * x[ 'average_weight_female' ] ,
                     B_a_unsexed = lambda x: x[ 'rho_a_unsexed' ] * x[ 'average_weight_total' ] )           
        )

        # Create dataframe to save it to Survey object
        areal_biomass_density_df = (
            nasc_to_areal_biomass_density_df
            .stretch( variable = 'B_a' )
            .merge( nasc_adult_number_proportions , on = [ 'stratum_num' ] )
            .assign( B_a_adult = lambda x: x[ 'B_a' ] * x[ 'adult_number_proportion' ] )
        )

        ### Calculate abundances (N)
        # N = # animals
        nasc_to_abundance_df = (
            nasc_to_areal_biomass_density_df
            # Convert to abundances (N)
            .assign( N_total = lambda x: x[ 'rho_a_total' ] * x[ 'interval_area' ] ,
                     N_male = lambda x: x[ 'rho_a_male' ] * x[ 'interval_area' ] ,
                     N_female = lambda x: x[ 'rho_a_female' ] * x[ 'interval_area' ] ,
                     N_unsexed = lambda x: x[ 'rho_a_unsexed' ] * x[ 'interval_area' ] )                     
        )

        # Create dataframe to save it to Survey object
        abundance_df = (
            nasc_to_abundance_df
            .stretch( variable = 'N' )
            .merge( nasc_adult_number_proportions , on = [ 'stratum_num' ] )
            .assign( N_adult = lambda x: x[ 'N' ] * x[ 'adult_number_proportion' ] )
        )
     
        ### Calculate biomass (B)
        # B = kg animals
        nasc_to_biomass_df = (
            nasc_to_areal_biomass_density_df
            # Convert to biomasses
            .assign( B_total = lambda x: x[ 'B_a_total' ] * x[ 'interval_area' ] ,
                     B_male = lambda x: x[ 'B_a_male' ] * x[ 'interval_area' ] ,
                     B_female = lambda x: x[ 'B_a_female' ] * x[ 'interval_area' ] ,
                     B_unsexed = lambda x: x[ 'B_a_unsexed' ] * x[ 'interval_area' ] )
        ) 

        # Create datafrane to save it to Survey object
        biomass_df = (
            nasc_to_biomass_df
            .stretch( variable = 'B' )
            .merge( nasc_adult_number_proportions , on = [ 'stratum_num' ] )
            .assign( B_adult = lambda x: x[ 'B' ] * x[ 'adult_number_proportion' ] )
        )

        ### Expand biomass by age ! 
        age_biomass_df = (
            biomass_df
            .merge( sex_indexed_weight_proportions , on = [ 'stratum_num' , 'sex' ] )
            .assign( B_age = lambda x: x[ 'B_adult' ] * x[ 'weight_proportion' ] )[ [ 'transect_num' , 'latitude' , 'longitude' , 'stratum_num' , 'age' , 'sex' , 'weight_proportion' , 'B_age' ] ]
            .rename( columns = { 'weight_proportion': 'age_proportion' } )
        )

        ### Update Survey object 
        self.biology[ 'population' ].update(
            {
                'areal_density': {
                    'number_density_df': areal_number_density_df ,
                    'biomass_density_df': areal_biomass_density_df ,
                } ,
                'abundance': {
                    'abundance_df': abundance_df ,
                } ,
                'biomass': {
                    'biomass_df': biomass_df ,
                    'biomass_age_df': age_biomass_df ,
                } ,
            }
        )

    #!!! TODO : Provide argument that will exclude age-0 and age-1 fish when flagged
    def stratified_summary( self ):
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
        ### Call in NASC dataframe
        nasc_df = self.acoustics[ 'nasc' ][ 'nasc_df' ].copy( deep = True )

        ### Call in other parameters for stratified statistic calculation
        # Transect sample faction
        transect_sample = self.config[ 'stratified_survey_mean_parameters' ][ 'strata_transect_proportion' ]

        # Number of replicates
        transect_replicates = self.config[ 'stratified_survey_mean_parameters' ][ 'replicates' ]

        ### Call in biomass
        adult_biomass = (
            self.biology[ 'population' ][ 'biomass' ][ 'biomass_df' ]
            .pipe( lambda df: df.loc[ df.sex == 'total' ] )
            .groupby('transect_num')[ [ 'B_adult' ] ]
            .agg( sum )
            .reset_index()
        )        

        ### Call in INPFC stratification
        inpfc_df = self.spatial[ 'inpfc_strata_df' ].rename( { 'stratum_num': 'stratum_inpfc' } , axis = 1 )

        # Bin the strata based on latitude limits 
        latitude_bins = np.concatenate( [ [ -90 ] , inpfc_df.northlimit_latitude ] )

        ### Summarize transect features
        transect_summary = (
                    nasc_df
                    .pipe( lambda df: calculate_transect_distance( df ) )
                    .merge( adult_biomass , on = [ 'transect_num' ] )
                    .assign( **{ 'stratum_inpfc': lambda x: pd.cut( x[ 'center_latitude' ] , 
                                                                    latitude_bins , 
                                                                    right = False , 
                                                                    labels=range( len( latitude_bins ) - 1 ) ) } ) # assign bin
                    .assign( stratum_inpfc = lambda x: x[ 'stratum_inpfc' ].astype( int ) + 1 )
                )

        ### Summarize strata features
        strata_summary = (
                    transect_summary
                    .loc[ : , [ 'stratum_inpfc' , 'transect_area' ] ]
                    .groupby( 'stratum_inpfc' )
                    .agg( [ 'count' , 'sum' ] )
                    .rename( columns= { "count": "num_transects" , "sum": "total_transect_area" } )
                    .droplevel( 0 , axis = 1 )
                    .reset_index( )
                )
        
        ### Combine the above parameters and data to calculate the stratified statistics
        stratified_results = stratified_transect_statistic( transect_summary , strata_summary , transect_sample , transect_replicates )
        
        ###
        self.biology[ 'population' ].update(
            {
                'stratified_results': stratified_results
            }
        )