from typing import Union, Optional
from pathlib import Path
import numpy as np
import copy
from .core import DATA_STRUCTURE
from .analysis import (
    acoustics_to_biology , 
    krige ,
    process_transect_data ,     
    stratified_summary
)

from .biology import (
    filter_species , 
    sum_strata_weight , 
    calculate_aged_unaged_proportions ,
    calculate_aged_biomass , 
    calculate_unaged_biomass ,
    apply_age_bins
)

from .utils import message as em
from .utils import load as el

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
        self.config = el.load_configuration( Path( init_config_path ) , 
                                             Path( survey_year_config_path ) )

        # Loading the datasets defined in the configuration files
        self.input = el.load_survey_data( self.config )

        # Initialize the `analysis` data attribute
        self.analysis = copy.deepcopy( DATA_STRUCTURE[ 'analysis' ] )

        # Initialize the `results` data attribute
        self.results = copy.deepcopy( DATA_STRUCTURE[ 'results' ] )

    def transect_analysis(self ,
                          species_id: np.float64 = 22500 ,
                          exclude_age1: bool = True ,
                          stratum: str = 'ks' ,
                          verbose: bool = True ):
        """
        Calculate population-level metrics from acoustic transect measurements
        """    

        # Update settings to reflect the stratum definition
        self.analysis[ 'settings' ].update(
            {
                'transect': {
                    'species_id': species_id ,
                    'stratum': stratum.lower( ) ,
                    'stratum_name': 'stratum_num' if stratum == 'ks' else 'inpfc' ,
                    'exclude_age1': exclude_age1
                }
            }
        )

        # Initial data processing of the transect biological and acoustic data
        self.analysis[ 'transect' ] = process_transect_data( self.input ,
                                                             self.analysis[ 'transect' ] ,
                                                             self.analysis[ 'settings' ] ,
                                                             self.config )

        # Convert NASC into number density (animals/nmi^2), biomass density (kg/nmi^2), abundance
        # (# animals), and biomass (kg) for all fish, sexed (male/female) fish, and unsexed fish
        # ---- This further provides the resulting distributions of biomass and abundance over 
        # ---- length and age for each sex across the entire survey
        biomass_summary , self.analysis[ 'transect' ] = (
            acoustics_to_biology( self.input ,
                                  self.analysis[ 'transect' ] ,
                                  self.config ,
                                  self.analysis[ 'settings' ] )
        )

        # Add the biomass summary table to the results attribute
        # ---- Update results (biomass summary)
        self.results[ 'transect' ].update( { 'biomass_summary_df': biomass_summary } )

        # Print result if `verbose == True`
        if verbose:
            em.transect_results_msg( biomass_summary )

    def stratified_analysis( self ,
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

        # Calculate the stratified mean, variance, and coefficient of variation
        stratified_results , self.analysis = (
            stratified_summary( self.analysis ,
                                self.analysis[ 'settings' ][ 'stratified' ] )
        )

        # Add the stratified statistics dictionary results to the `results` attribute
        # ---- Update results (stratified results)
        self.results[ 'stratified' ].update( stratified_results )

        # Print result if `verbose == True`
        if verbose:
            em.stratified_results_msg( stratified_results , 
                                       self.analysis[ 'settings' ][ 'stratified' ] )
            
    def kriging_analysis( self ,
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

        # Run kriging analysis
        # ----> Generates a georeferenced dataframe for the entire mesh grid, summary statistics,
        # ----> and adds intermediate data products to the analysis attribute
        # ---- If kriging results are not extrapolated beyond the survey region:
        if not extrapolate: 
            # ---- Update the analysis settings
            self.analysis['settings']['kriging'].update(
                {
                    'mesh_buffer_distance': mesh_buffer_distance ,
                    'num_nearest_transect': num_nearest_transects 
                }
            )
            # ---- Run kriging algorithm
            kriged_results , self.analysis = (
                krige( self.input , self.analysis , self.analysis[ 'settings' ][ 'kriging' ] )
            ) 
        # ---- If kriging results are extrapolated beyond the survey region:
        else: 
            # ---- Run kriging algorithm
            kriged_results , self.analysis = (
                krige( self.input , self.analysis , self.analysis[ 'settings' ][ 'kriging' ] )
            ) 
        
        # Save the results to the `results` attribute
        self.results.update( { 'kriging': kriged_results } )

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
        length_spp , specimen_spp , haul_spp = filter_species( [ self.input['biology'][ 'length_df' ] ,
                                                                 self.input['biology'][ 'specimen_df' ] ,
                                                                 self.input['biology'][ 'catch_df' ] ] ,
                                                               species_id )
        
        # ---- Remove 'bad' values 
        # ---- `specimen_spp`
        specimen_spp_filtered = specimen_spp[ specimen_spp.sex != 'unsexed' ].dropna( how = 'any' , subset = 'age' )

        # ---- `length_spp`
        length_spp_filtered = length_spp[ length_spp.sex != 'unsexed' ]

        ### Import discrete distribution bins
        # ---- Length
        length_intervals = self.input['biology'][ 'distributions' ][ 'length' ][ 'length_interval_arr' ]
        
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