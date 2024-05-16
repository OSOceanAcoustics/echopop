"""
General analysis orchestration functions that bundle related functions and procedures
"""

import pandas as pd
import numpy as np
import copy

from echopop.computation.acoustics import summarize_sigma_bs , nasc_to_biomass
from echopop.computation.statistics import stratified_transect_statistic

from echopop.spatial.mesh import crop_mesh
from echopop.spatial.projection import transform_geometry
from echopop.spatial.krige import kriging

from echopop.biology import (
    filter_species , 
    fit_length_weight_relationship , 
    number_proportions , 
    quantize_number_counts , 
    quantize_weights ,
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

from echopop.spatial.transect import (
    edit_transect_columns , 
    transect_distance , 
    summarize_transect_strata , 
    save_transect_coordinates
)

def process_transect_data( input_dict: dict ,
                           analysis_dict: dict ,
                           settings_dict: dict ,
                           configuration_dict: dict ):
    
    # Filter out non-target species
    length_data , specimen_data , catch_data = (
        filter_species( [ input_dict[ 'biology' ][ 'length_df' ] ,
                        input_dict[ 'biology' ][ 'specimen_df' ] ,
                        input_dict[ 'biology' ][ 'catch_df' ] ] ,
                        settings_dict[ 'transect' ][ 'species_id' ] )
    )
    # ---- For cases where all samples were aged (i.e. in `specimen_data` and absent from
    # ---- `length_data`), these hauls are removed from `catch_data`
    catch_data = catch_data[ catch_data.haul_num.isin( length_data.haul_num ) ]
    
    # Save the transect coordinate information
    analysis_dict.update(
        { 'coordinates': save_transect_coordinates( input_dict[ 'acoustics' ][ 'nasc_df' ] ) }
    )

    # Calculate mean sigma_bs per individual haul, KS stratum, and INPFC stratum
    analysis_dict[ 'acoustics' ][ 'sigma_bs' ].update(
        summarize_sigma_bs( length_data , 
                            specimen_data , 
                            input_dict[ 'spatial' ] , 
                            configuration_dict ,
                            settings_dict )        
    )

    # Fit length-weight regression required for biomass calculation
    analysis_dict[ 'biology' ][ 'weight' ].update(
        fit_length_weight_relationship( 
            specimen_data , 
            input_dict[ 'biology' ][ 'distributions' ][ 'length_bins_df' ] )            
    )

    # Count the number of specimens across age and length bins
    analysis_dict[ 'biology' ][ 'distributions' ].update(
        quantize_number_counts( specimen_data , 
                                length_data , 
                                stratum = settings_dict[ 'transect' ][ 'stratum' ] )
    )

    # Calculate the number proportions
    analysis_dict[ 'biology' ][ 'proportions' ].update(
        {
            'number': number_proportions( analysis_dict[ 'biology' ][ 'distributions' ] )
        }
    )

    # Sum the weights of both aged and unaged specimens across length and weight bins
    # ---- Extract the length-weight fit
    length_weight_df = (
        analysis_dict['biology']['weight']['length_weight_regression']['weight_fitted_df']
    )
    # ---- Quantize the weights
    analysis_dict[ 'biology' ][ 'proportions' ].update(
        {
            'weight': quantize_weights( specimen_data ,
                                        length_data ,
                                        length_weight_df )
        }
    )

    # Calculate the average weights among male, female, and all fish across strata
    analysis_dict[ 'biology' ][ 'weight' ].update(
        {
            'weight_stratum_df': fit_length_weights( analysis_dict[ 'biology' ][ 'proportions' ][ 'number' ] ,
                                                     analysis_dict[ 'biology' ][ 'weight' ] )
        }
    )

    # Calculate the weight proportions
    analysis_dict[ 'biology' ][ 'proportions' ].update(
        {
            'weight': weight_proportions( specimen_data , 
                                          length_data , 
                                          catch_data , 
                                          analysis_dict['biology']['weight']['length_weight_regression']['weight_fitted_df'] )
        }
    )

    # Calculate the weight proportions
    analysis_dict[ 'biology' ][ 'proportions' ].update(
        {
            'weight': weight_proportions( specimen_data , 
                                          length_data , 
                                          catch_data , 
                                          analysis_dict['biology']['weight']['length_weight_regression']['weight_fitted_df'] )
        }
    )

    # Return the analysis dictionary
    return analysis_dict 

def acoustics_to_biology( input_dict: dict ,
                          analysis_dict: dict ,
                          configuration_dict: dict ,
                          settings_dict: dict ):
    
    # Convert NASC into number density (animals/nmi^2), biomass density (kg/nmi^2), abundance
    # (# animals), and biomass (kg) for all fish, sexed (male/female) fish, and unsexed fish
    strata_adult_proportions , nasc_to_biology = (
        nasc_to_biomass( input_dict ,
                            analysis_dict ,
                            configuration_dict ,
                            settings_dict )
    ) 

    # Distribute abundance and biomass estimates over length and age bins
    analysis_dict[ 'biology' ][ 'population' ].update(
        {
            'tables': distribute_length_age( nasc_to_biology ,                                                 
                                             analysis_dict[ 'biology' ][ 'proportions' ] ,
                                             settings_dict )
        }    
    )

    # Reapportion transect results to separate age-1 and age-2+ fish, generate age-1
    # abundance distributions for unaged fish, and generate biomass summary
    adult_transect , biomass_summary , unaged_age1_abundance = (
        partition_transect_age( nasc_to_biology , 
                                analysis_dict[ 'biology' ][ 'weight' ] ,
                                settings_dict ,
                                analysis_dict[ 'biology' ][ 'population' ] ,
                                strata_adult_proportions )
    )   

    # Append outputs to their appropriate positions within the analysis attribute
    # ---- Adult transect data
    analysis_dict[ 'acoustics' ].update( { 'adult_transect_df': adult_transect } )
    # ---- Age-1 abundances (unaged fish)
    analysis_dict[ 'biology' ][ 'population' ][ 'tables' ][ 'abundance' ].update(
        { 'unaged_age1_abundance_df': unaged_age1_abundance }
    )

    # Return the updated analysis attribute and `biomass_summary`
    return biomass_summary , analysis_dict

def stratified_summary( analysis_dict: dict ,
                        settings_dict: dict ) :
    
    # Define the and prepare the processed and georeferenced transect data
    transect_data = edit_transect_columns( analysis_dict[ 'transect' ] , settings_dict )

    # Summarize transect spatial information
    # ---- Transect distances
    transect_summary = transect_distance( transect_data )
    # ---- Counts and area coverage per stratum
    strata_summary = summarize_transect_strata( transect_summary )

    # Compute the stratified mean, variance, and coefficient of variation (CV)
    # ---- This includes the statistical (Gaussian) estimates (mean) and 95% confidence interval
    # ---- for each statistic
    replicates , stratified_results = stratified_transect_statistic( transect_data  ,
                                                                     transect_summary ,
                                                                     strata_summary ,
                                                                     settings_dict )
    
    # Update the analysis attribute with the resampled/bootstrapped replicates
    analysis_dict[ 'stratified' ].update( { 'stratified_replicates_df': replicates } )
    
    # Return the outputs
    return stratified_results , analysis_dict

def krige( input_dict: dict ,
           analysis_dict: dict ,
           settings_dict: dict ):
    
    # Extract kriging mesh data
    mesh_data = input_dict[ 'statistics' ][ 'kriging' ][ 'mesh_df' ]

    # Extract the reference grid (200 m isobath)
    isobath_data = input_dict[ 'statistics' ][ 'kriging' ][ 'isobath_200m_df' ] 

    # Define the and prepare the processed and georeferenced transect data
    transect_data = edit_transect_columns( analysis_dict[ 'transect' ] , 
                                           settings_dict )
    
    # Crop the mesh grid if the kriged data will not be extrapolated
    if not settings_dict[ 'extrapolate' ]:
        # ---- Compute the cropped mesh
        mesh_full = crop_mesh( transect_data ,
                               mesh_data ,
                               settings_dict )
        if settings_dict[ 'verbose' ]:
        # ---- Print alert 
            print( "Kriging mesh cropped to prevent extrapolation beyond the defined "
            f"""`mesh_buffer_distance` value ({settings_dict['mesh_buffer_distance']} nmi).""")
        
    else:
        # ---- Else, extract original mesh dataframe
        mesh_df = mesh_data.copy( )
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
    if settings_dict['standardize_coordinates' ]:
        # ---- Transform transect data geometry (generate standardized x- and y-coordinates)
        transect_data , d_x , d_y = (
            transform_geometry( transect_data ,
                                isobath_data  ,
                                settings_dict )
        ) 
        # ---- Transform mesh grid geometry (generate standardized x- and y-coordinates)
        mesh_full , _ , _ = (
            transform_geometry( mesh_full ,
                                isobath_data  ,
                                settings_dict ,
                                d_x , d_y )
        ) 
        if settings_dict[ 'verbose' ]:
            # ---- Print alert 
            print( """Longitude and latitude coordinates (WGS84) converted to standardized """
            """coordinates (x and y).""")
    else:
        # ---- Else, duplicate the transect longitude and latitude coordinates as 'x' and 'y'
        # -------- x
        transect_data[ 'x' ] = transect_data[ 'longitude' ]
        # -------- y
        transect_data[ 'y' ] = transect_data[ 'latitude' ]
        # ---- Duplicate the mesh grid longitude and latitude coordinates as 'x' and 'y'
        # -------- x
        mesh_full[ 'x' ] = mesh_full[ 'longitude' ]
        # -------- y
        mesh_full[ 'y' ] = mesh_full[ 'latitude' ]
    # --- Append to the analysis attribute
    analysis_dict.update( { 'kriging': { 'mesh_df': mesh_full , 'transect_df': transect_data } } )

    # Kriged results
    kriged_results = kriging( analysis_dict[ 'kriging' ][ 'transect_df' ] ,
                              analysis_dict[ 'kriging' ][ 'mesh_df' ] ,
                              settings_dict )
    
    # Return kriged (interpolated) results
    return kriged_results , analysis_dict
