from echopop.survey import Survey
import numpy as np 
import pandas as pd

import geopandas as gpd
import geopy.distance

from typing import Union

from echopop.spatial.projection import wgs84_to_utm
from echopop.spatial.transect import transect_extent

survey = Survey( init_config_path = "./config_files/initialization_config.yml" ,
                 survey_year_config_path = "./config_files/survey_year_2019_config.yml" )
survey.transect_analysis( )

survey.stratified_analysis( )

survey.kriging_analysis( )

survey.stratified_analysis( dataset = 'kriging' )

Number of { 'virtual transects' if settings_dict[ 'dataset' ] == 'kriging' else 'transects' } : { stratified_results[ 'num_transects' ] }
am = transect_distances
bm = transect_data.groupby( [ 'transect_num' ] )[ 'biomass' ].sum( ) * 1e-6
km = analysis_dict[ 'transect' ]['acoustics']['adult_transect_df']

np.mean( km[ 'biomass' ] ) * 1e-6

survey.transect_analysis( )

survey.stratified_analysis( )

survey.kriging_analysis( )

survey.stratified_analysis( dataset = 'kriging' )
self = survey

survey.results[ 'stratified' ][ 'kriging' ][ 'mean' ][ 'unweighted_estimate' ] * 1e-6 / 0.75
survey.results[ 'kriging' ][ 'mesh_results_df' ][ 'biomass' ].sum() / survey.results[ 'kriging' ][ 'mesh_results_df' ][ 'area' ].sum()

from echopop.spatial.transect import transect_distance

dlat = 1.25 / 60.0
transect_coords = self.analysis[ 'transect' ]['coordinates' ].copy( )
lat_t = transect_coords[ 'latitude' ].to_numpy() ; lon_t = transect_coords[ 'longitude' ].to_numpy()
US_CAN_Mesh = self.input[ 'statistics' ][ 'kriging' ][ 'mesh_df' ].copy( )
lat_k = US_CAN_Mesh[ 'centroid_latitude' ].to_numpy() ; lon_k = US_CAN_Mesh[ 'centroid_longitude' ].to_numpy()
tx = transect_coords[ 'transect_num' ].to_numpy()
uniq_tx = np.unique( tx )
nt = len( uniq_tx )
lat_b1=49 ; lat_b1=34.3
i0=1

transect_data = self.input[ 'acoustics' ][ 'nasc_df' ].copy()
transect_headings = transect_bearing( transect_data )
# ---- Roudn to the nearest 10.0 degree increment
transect_headings[ 'heading' ] = np.round( transect_headings[ 'heading' ] )
# ---- Find N-S 
region_2_bearings = transect_headings[ transect_headings[ 'heading' ] < 15 ]
# filter
transect_data_2 = transect_data[ transect_data[ 'transect_num' ].isin( region_2_bearings[ 'transect_num' ] )]
# ---- region 1
region_1 = transect_headings.copy()[ transect_headings[ 'transect_num' ] < region_2_bearings[ 'transect_num' ].min() ]
transect_data_1 = transect_data.copy()[ transect_data[ 'transect_num' ].isin( region_1[ 'transect_num' ] )]
td_1 = transect_distance( transect_data_1 )
td_1.set_index( 'transect_num' , inplace = True )
# ---- region 3
region_3 = transect_headings[ transect_headings[ 'transect_num' ] > region_2_bearings[ 'transect_num' ].max() ]
transect_data_3 = transect_data[ transect_data[ 'transect_num' ].isin( region_3[ 'transect_num' ] )]

# increments -- region 1
region_1_steps = (
    np.arange( transect_data_1[ 'latitude' ].min( ) , transect_data_1[ 'latitude' ].max( ) , dlat )
)

transect_data_1[ 'longitude_west' ] = (
    transect_data_1.groupby( [ 'transect_num' ] )[ 'longitude' ].transform( 'min' )
)
transect_data_1[ 'longitude_east' ] = (
    transect_data_1.groupby( [ 'transect_num' ] )[ 'longitude' ].transform( 'max' )
)

transect_1_west = transect_data_1[ transect_data_1[ 'longitude' ] == transect_data_1[ 'longitude_west' ] ].sort_values( [ 'latitude' ] )
transect_1_east = transect_data_1[ transect_data_1[ 'longitude' ] == transect_data_1[ 'longitude_east' ] ].sort_values( [ 'latitude' ] )

from scipy import interpolate

interpolator_west = interpolate.interp1d(
                    transect_1_west[ 'latitude' ] ,
                    transect_1_west[ 'longitude' ] ,
                    kind = 'linear' ,
                    bounds_error = False
                )
interpolator_east = interpolate.interp1d(
                    transect_1_east[ 'latitude' ] ,
                    transect_1_east[ 'longitude' ] ,
                    kind = 'linear' ,
                    bounds_error = False
                )

# new_lon_1w = np.interp( region_1_steps , transect_1_west[ 'latitude' ] , transect_1_west[ 'longitude' ] )
# new_lon_1e = np.interp( region_1_steps , transect_1_east[ 'latitude' ] , transect_1_east[ 'longitude' ] )
new_lon_1w = interpolator_west( region_1_steps )
new_lon_1e = interpolator_east( region_1_steps )
dlon = dlat * np.cos( np.radians( region_1_steps ) )
US_CAN_Mesh.rename( columns = { 'centroid_latitude': 'latitude' , 'centroid_longitude': 'longitude' } , inplace = True )
dlon = np.zeros( len( new_lon_1w ) )

indx_1 = []
for i in range( len( new_lon_1w ) ) :
    dlon[ i ] = dlat * np.cos( np.radians( region_1_steps[i] ) )
    idx = (
        np.where( ( US_CAN_Mesh[ 'longitude' ] >= new_lon_1w[i] - dlon[ i ] ) & 
                  ( ( US_CAN_Mesh[ 'longitude' ] <= new_lon_1e[i] + dlon[ i ] ) ) & 
                  ( US_CAN_Mesh[ 'latitude' ] >= region_1_steps[ i ] - dlat ) & 
                  ( US_CAN_Mesh[ 'latitude' ] < region_1_steps[ i ] + dlat ) )
    )
    indx_1.append( idx[ 0 ] )

indx1 = np.unique( np.concatenate( indx_1 ) )
len(indx1)
transect_data_2[ 'latitude_south' ] = (
    transect_data_2.groupby( [ 'transect_num' ] )[ 'latitude' ].transform( 'min' )
)
transect_data_2[ 'latitude_mean' ] = (
    transect_data_2.groupby( [ 'transect_num' ] )[ 'latitude' ].transform( 'mean' )
)
transect_data_2[ 'latitude_north' ] = (
    transect_data_2.groupby( [ 'transect_num' ] )[ 'latitude' ].transform( 'max' )
)

transect_2_south = transect_data_2[ transect_data_2[ 'latitude' ] == transect_data_2[ 'latitude_south' ] ].sort_values( [ 'longitude' ] )
transect_2_north = transect_data_2[ transect_data_2[ 'latitude' ] == transect_data_2[ 'latitude_north' ] ].sort_values( [ 'longitude' ] )
lat_m = transect_data_2[ 'latitude' ].mean( )
dlon = ( dlat * np.cos( np.radians( transect_2_north.loc[ transect_2_north[ 'transect_num' ] == transect_2_north[ 'transect_num' ].max( ) , 'latitude_mean' ] ) ) ).values
region_2_steps = (
    np.arange( transect_data_2[ 'longitude' ].min( ) , transect_data_2[ 'longitude' ].max( ) , dlon )
)
# find the lat/lon of the southern bounds / interpolate
interpolator_south = interpolate.interp1d(
                        transect_2_south[ 'longitude' ] ,
                        transect_2_south[ 'latitude' ] ,
                        kind = 'linear' ,
                        bounds_error = False
                    )

interpolator_north = interpolate.interp1d(
                        transect_2_north[ 'longitude' ] ,
                        transect_2_north[ 'latitude' ] ,
                        kind = 'linear' ,
                        bounds_error = False
                    )

new_lat_2s = interpolator_south( region_2_steps )
# new_lat_2s = np.interp( region_2_steps , transect_2_south[ 'longitude' ] , transect_2_south[ 'latitude' ] )
new_lat_2n = interpolator_north( region_2_steps )
# new_lat_2n = np.interp( region_2_steps , transect_2_north[ 'longitude' ] , transect_2_north[ 'latitude' ] )
dlon = ( dlat * np.cos( np.radians( transect_2_north.loc[ transect_2_north[ 'transect_num' ] == transect_2_north[ 'transect_num' ].max( ) , 'latitude_mean' ] ) ) ).values

indx_2 = []

for i in range( len( region_2_steps ) ) :
    if ( ~ np.isnan( new_lat_2s[i] ) ) | ( ~ np.isnan( new_lat_2n[i] ) ):
        if ( np.isnan( new_lat_2s[i] ) ) :
            lon_n_min = np.argmin( np.abs( region_2_steps[ i ] - transect_2_north[ 'longitude' ] ) )
            lon_s_min = np.argmin( np.abs( region_2_steps[ i ] - transect_2_south[ 'longitude' ] ) )
            slp = ( 
                ( transect_2_north[ 'latitude' ].iloc[ lon_n_min ] 
                - transect_2_south[ 'latitude' ].iloc[ lon_s_min ] )
                / ( transect_2_north[ 'longitude' ].iloc[ lon_n_min ] 
                - transect_2_south[ 'longitude' ].iloc[ lon_s_min ] )
            )
            lat_s2_i = (
                slp * ( region_2_steps[ i ] - transect_2_south[ 'longitude' ].iloc[ lon_s_min ] ) 
                + transect_2_south[ 'latitude' ].iloc[ lon_s_min ]
            )

            idx = (
                np.where( ( US_CAN_Mesh[ 'longitude' ] >= region_2_steps[i] - dlon[0] ) & 
                        ( US_CAN_Mesh[ 'longitude' ] <= region_2_steps[i] + dlon[0]  ) & 
                        ( US_CAN_Mesh[ 'latitude' ] >= lat_s2_i - dlat ) & 
                        ( US_CAN_Mesh[ 'latitude' ] < new_lat_2n[ i ] + dlat ) )
            )
        elif ( np.isnan( new_lat_2n[i] ) ):
            lon_n_min = np.argmin( np.abs( region_2_steps[ i ] - transect_2_north[ 'longitude' ] ) )
            lon_s_min = np.argmin( np.abs( region_2_steps[ i ] - transect_2_south[ 'longitude' ] ) )
            slp = ( 
                ( transect_2_north[ 'latitude' ].iloc[ lon_n_min ] 
                - transect_2_south[ 'latitude' ].iloc[ lon_s_min ] )
                / ( transect_2_north[ 'longitude' ].iloc[ lon_n_min ] 
                - transect_2_south[ 'longitude' ].iloc[ lon_s_min ] )
            )
            lat_n2_i = (
                slp * ( region_2_steps[ i ] - transect_2_south[ 'longitude' ].iloc[ lon_s_min ] ) 
                + transect_2_south[ 'latitude' ].iloc[ lon_s_min ]
            )   
            idx = (
                np.where( ( US_CAN_Mesh[ 'longitude' ] >= region_2_steps[i] - dlon[0] ) & 
                        ( US_CAN_Mesh[ 'longitude' ] <= region_2_steps[i] + dlon[0]  ) & 
                        ( US_CAN_Mesh[ 'latitude' ] >= new_lat_2s[ i ] - dlat ) & 
                        ( US_CAN_Mesh[ 'latitude' ] < lat_n2_i + dlat ) )
            )  
        else:
            idx = (
                np.where( ( US_CAN_Mesh[ 'longitude' ] >= region_2_steps[i] - dlon[0] ) & 
                        ( US_CAN_Mesh[ 'longitude' ] <= region_2_steps[i] + dlon[0]  ) & 
                        ( US_CAN_Mesh[ 'latitude' ] >= new_lat_2s[ i ] - dlat ) & 
                        ( US_CAN_Mesh[ 'latitude' ] < new_lat_2n[ i ] + dlat ) )
            )

    indx_2.append( idx[ 0 ] )

indx2 = np.unique( np.concatenate( indx_2 ) )
len( indx2 )

# find the lat/lon of the region 3 transects

transect_data_3[ 'latitude_south' ] = (
    transect_data_3.groupby( [ 'transect_num' ] )[ 'latitude' ].transform( 'min' )
)
transect_data_3[ 'latitude_north' ] = (
    transect_data_3.groupby( [ 'transect_num' ] )[ 'latitude' ].transform( 'max' )
)
transect_data_3[ 'longitude_west' ] = (
    transect_data_3.groupby( [ 'transect_num' ] )[ 'longitude' ].transform( 'min' )
)
transect_data_3[ 'longitude_east' ] = (
    transect_data_3.groupby( [ 'transect_num' ] )[ 'longitude' ].transform( 'max' )
)
region_3_steps = (
    np.arange( transect_data_3[ 'latitude' ].min( ) , transect_data_3[ 'latitude' ].max( ) , dlat )
)

transect_3_west = transect_data_3[ transect_data_3[ 'longitude' ] == transect_data_3[ 'longitude_west' ] ].sort_values( [ 'longitude' ] )
transect_3_east = transect_data_3[ transect_data_3[ 'longitude' ] == transect_data_3[ 'longitude_east' ] ].sort_values( [ 'longitude' ] )

interpolator_west3 = interpolate.interp1d(
                        transect_3_west[ 'latitude' ] ,
                        transect_3_west[ 'longitude' ] ,
                        kind = 'linear' ,
                        bounds_error = False
                    )
interpolator_east3 = interpolate.interp1d(
                        transect_3_east[ 'latitude' ] ,
                        transect_3_east[ 'longitude' ] ,
                        kind = 'linear' ,
                        bounds_error = False
                    )

new_lon_3w = interpolator_west3( region_3_steps )
new_lon_3e = interpolator_east3( region_3_steps )

indx_3 = []
dlon = np.zeros( len( region_3_steps ) )

for i in range( len( region_3_steps ) ) :
    dlon[ i ] = dlat * np.cos( np.radians( region_3_steps[ i ] ) )

    if np.isnan( new_lon_3w[ i ] ):
        lon_w_max = np.argmax( transect_3_west[ 'latitude' ] )
        lon_e_max = np.argmax( transect_3_east[ 'latitude' ] )

        slp = ( 
            ( transect_3_west[ 'longitude' ].iloc[ lon_w_max ] 
            - transect_3_east[ 'longitude' ].iloc[ lon_e_max ] ) 
            / ( transect_3_west[ 'latitude' ].max( ) 
            - transect_3_east[ 'latitude' ].max( ) )
        )

        lon_w3_i = ( 
            slp * ( region_3_steps[ i ] - transect_3_east[ 'latitude' ].max( ) ) 
            + transect_3_east[ 'longitude' ].iloc[ lon_e_max ]
        )

        idx = (
            np.where( 
                ( US_CAN_Mesh[ 'longitude' ] >= lon_w3_i - dlon[i] ) & 
                ( US_CAN_Mesh[ 'longitude' ] <= new_lon_3e[i] + dlon[i] ) & 
                ( US_CAN_Mesh[ 'latitude' ] >= region_3_steps[i] - dlat ) & 
                ( US_CAN_Mesh[ 'latitude' ] < region_3_steps[i] + dlat ) )
        )
    elif np.isnan( new_lon_3e[i] ) :

        lon_w_max = np.argmax( transect_3_west[ 'latitude' ] )
        lon_e_max = np.argmax( transect_3_east[ 'latitude' ] )

        slp = ( 
            ( transect_3_west[ 'longitude' ].iloc[ lon_w_max ] 
            - transect_3_east[ 'longitude' ].iloc[ lon_e_max ] ) 
            / ( transect_3_west[ 'latitude' ].max( ) 
            - transect_3_east[ 'latitude' ].max( ) )
        )

        lon_e3_i = ( 
            slp * ( region_3_steps[ i ] - transect_3_east[ 'latitude' ].max( ) ) 
            + transect_3_east[ 'longitude' ].iloc[ lon_e_max ]
        )

        idx = (
            np.where(
                ( US_CAN_Mesh[ 'longitude' ] >= new_lon_3w[ i ] - dlon[i] )
                ( US_CAN_Mesh[ 'longitude' ] <= lon_e3_i + dlon[i] ) & 
                ( US_CAN_Mesh[ 'latitude' ] >= region_3_steps[i] - dlat ) & 
                ( US_CAN_Mesh[ 'latitude' ] < region_3_steps[i] + dlat )
            )
        )
    
    else:
        idx = (
            np.where(
                ( US_CAN_Mesh[ 'longitude' ] >= new_lon_3w[ i ] - dlon[i] ) &
                ( US_CAN_Mesh[ 'longitude' ] <= new_lon_3e[i] + dlon[i] ) & 
                ( US_CAN_Mesh[ 'latitude' ] >= region_3_steps[i] - dlat ) & 
                ( US_CAN_Mesh[ 'latitude' ] < region_3_steps[i] + dlat )
            )
        )

    indx_3.append( idx[ 0 ] )
        
indx3 = np.unique( np.concatenate( indx_3 ) )

len( indx1 ) + len( indx2 ) + len( indx3 )

indx_full = np.concatenate( [ indx1 , indx2 , indx3 ] )
mesh_cropped = US_CAN_Mesh.loc[ indx_full ]

import matplotlib.pyplot as plt

plt.scatter( mesh_data_cropped[ 'longitude' ] , mesh_data_cropped[ 'latitude' ] )
plt.show()


transect_3_east = transect_data_2[ transect_data_2[ 'latitude' ] == transect_data_2[ 'latitude_north' ] ].sort_values( [ 'longitude' ] )
lat_m = transect_data_2[ 'latitude' ].mean( )
dlon = ( dlat * np.cos( np.radians( transect_2_north.loc[ transect_2_north[ 'transect_num' ] == transect_2_north[ 'transect_num' ].max( ) , 'latitude_mean' ] ) ) ).values
region_2_steps = (
    np.arange( transect_data_2[ 'longitude' ].min( ) , transect_data_2[ 'longitude' ].max( ) , dlon )
)

interpolator_south = interpolate.interp1d(
                        transect_2_south[ 'longitude' ] ,
                        transect_2_south[ 'latitude' ] ,
                        kind = 'linear' ,
                        bounds_error = False
                    )

interpolator_north = interpolate.interp1d(
                        transect_2_north[ 'longitude' ] ,
                        transect_2_north[ 'latitude' ] ,
                        kind = 'linear' ,
                        bounds_error = False
                    )

new_lat_2s = interpolator_south( region_2_steps )
new_lat_2n = interpolator_north( region_2_steps )

indx_2 = []

for i in range( len( region_2_steps ) ) :
    if ( np.isnan( new_lat_2s[i] ) ) | ( np.isnan( new_lat_2n[i] ) ):
        lon_n_min = np.argmin( np.abs( region_2_steps[ i ] - transect_2_north[ 'longitude' ] ) )
        lon_s_min = np.argmin( np.abs( region_2_steps[ i ] - transect_2_south[ 'longitude' ] ) )
        slp = ( 
            ( transect_2_north[ 'latitude' ].iloc[ lon_n_min ] 
            - transect_2_south[ 'latitude' ].iloc[ lon_s_min ] )
            / ( transect_2_north[ 'longitude' ].iloc[ lon_n_min ] 
            - transect_2_south[ 'longitude' ].iloc[ lon_s_min ] )
        )
        lat_s2_i = (
            slp * ( region_2_steps[ i ] - transect_2_south[ 'longitude' ].iloc[ lon_s_min ] ) 
            + transect_2_south[ 'latitude' ].iloc[ lon_s_min ]
        )
        idx = (
            np.where( ( US_CAN_Mesh[ 'longitude' ] >= region_2_steps[i] - dlon[0] ) & 
                    ( US_CAN_Mesh[ 'longitude' ] <= region_2_steps[i] + dlon[0]  ) & 
                    ( US_CAN_Mesh[ 'latitude' ] >= lat_s2_i - dlat ) & 
                    ( US_CAN_Mesh[ 'latitude' ] < new_lat_2n[ i ] + dlat ) )
        )
    elif ( np.isnan( new_lat_2n[i] ) ):
        lon_n_min = np.argmin( np.abs( region_2_steps[ i ] - transect_2_north[ 'longitude' ] ) )
        lon_s_min = np.argmin( np.abs( region_2_steps[ i ] - transect_2_south[ 'longitude' ] ) )
        slp = ( 
            ( transect_2_north[ 'latitude' ].iloc[ lon_n_min ] 
            - transect_2_south[ 'latitude' ].iloc[ lon_s_min ] )
            / ( transect_2_north[ 'longitude' ].iloc[ lon_n_min ] 
            - transect_2_south[ 'longitude' ].iloc[ lon_s_min ] )
        )
        lat_n2_i = (
            slp * ( region_2_steps[ i ] - transect_2_south[ 'longitude' ].iloc[ lon_s_min ] ) 
            + transect_2_south[ 'latitude' ].iloc[ lon_s_min ]
        )   
        idx = (
            np.where( ( US_CAN_Mesh[ 'longitude' ] >= region_2_steps[i] - dlon[0] ) & 
                    ( US_CAN_Mesh[ 'longitude' ] <= region_2_steps[i] + dlon[0]  ) & 
                    ( US_CAN_Mesh[ 'latitude' ] >= new_lat_2s[ i ] - dlat ) & 
                    ( US_CAN_Mesh[ 'latitude' ] < lat_n2_i + dlat ) )
        )  
    else:
        idx = (
            np.where( ( US_CAN_Mesh[ 'longitude' ] >= region_2_steps[i] - dlon[0] ) & 
                    ( US_CAN_Mesh[ 'longitude' ] <= region_2_steps[i] + dlon[0]  ) & 
                    ( US_CAN_Mesh[ 'latitude' ] >= new_lat_2s[ i ] - dlat ) & 
                    ( US_CAN_Mesh[ 'latitude' ] < new_lat_2n[ i ] + dlat ) )
        )

    indx_2.append( idx[ 0 ] )

indx2 = np.unique( np.concatenate( indx_2 ) )



len( indx1 ) + len( indx2 )






transect_data.loc[3][ 'longitude' ].min( ) == region_1.loc[3, 'longitude_min']
region_1.loc[3, 'longitude_min']
region_1.set_index( 'transect_num' , inplace = True )
transect_data.set_index( 'transect_num' , inplace = True )
transect_data_1.set_index( 'transect_num' , inplace = True )
transect_data_1.loc[ : , 'lon_w' ] = td_1.copy()['longitude_min']

from scipy import interpolate

region_1[ 'lat_w' ] = transect_data_1[ transect_data_1[ 'longitude' ] == transect_data_1[ 'lon_w' ] ]['latitude']



transect_data


np.degrees( region_1[ 'longitude_min' ] )
np.degrees( region_1[ 'longitude1' ] )



transect_data_1[ 'latitude' ].min( )

transect_coords = self.analysis[ 'transect' ]['coordinates' ].copy( )
lat_t = transect_coords[ 'latitude' ].to_numpy() ; lon_t = transect_coords[ 'longitude' ].to_numpy()
US_CAN_Mesh = self.input[ 'statistics' ][ 'kriging' ][ 'mesh_df' ].copy( )
lat_k = US_CAN_Mesh[ 'centroid_latitude' ].to_numpy() ; lon_k = US_CAN_Mesh[ 'centroid_longitude' ].to_numpy()
tx = transect_coords[ 'transect_num' ].to_numpy()
uniq_tx = np.unique( tx )
nt = len( uniq_tx )
lat_b1=49 ; lat_b1=34.3
i0=1
lon1 = np.zeros( nt )
lon2 = np.zeros( nt )
lat1 = np.zeros( nt )
lat2 = np.zeros( nt )

for i in np.arange( i0 , nt , 1 ) :
    ind=np.where(tx == uniq_tx[i-1])[0]
    lat=lat_t[ind]
    lon=lon_t[ind]

    if np.max(lat)-np.min(lat) < np.max(lon)-np.min(lon) :
        indx1 = np.where( lon == np.min( lon ) )[0]
        val = np.min(lon) 

        lon1[i-1]=val
        lat1[i-1]=lat[indx1][0]

        indx2 = np.where( lon == np.max( lon ) )[0]
        val = np.max(lon)
        lon1[i-1]=val
        lat1[i-1]=lat[indx2][0]   

    else: 
        indx1 = np.where( lat == np.min( lat ) )[0]
        val = np.min(lat)    

        lat1[i-1]=val
        lon1[i-1]=lon[indx1][0]

        indx2 = np.where( lat == np.max( lat ) )[0]
        val = np.max(lat)

        lat2[i-1]=val
        lon2[i-1]=lon[indx2][0]

dlat=1.25/60
lat_mean = transect_coords.groupby( [ 'transect_num' ] )[ 'latitude' ].mean( ).reset_index( )
region1 = lat_mean[ lat_mean[ 'latitude' ] < 51.9 ]
tx0 = np.min( region1[ 'transect_num' ] )
tx1 = np.max( region1[ 'transect_num' ] )
tt = transect_coords.groupby( [ 'transect_num' ] )[ 'latitude' ].mean( ).sort_values().reset_index()
tt[ tt[ 'transect_num' ] == 119 ]
tt[ tt[ 'latitude' ] <= 58 ].sort_values( by = 'latitude' , ascending = False )


def crop_mesh( transect_data: pd.DataFrame ,
               mesh_data: pd.DataFrame ,               
               settings_dict: dict ):
    
    # Extract the analysis settings
    # ---- Number of nearest transects
    num_nearest_transects = settings_dict[ 'num_nearest_transect' ]
    # ---- Grid buffer distance (nmi)
    mesh_buffer_distance = settings_dict[ 'mesh_buffer_distance' ]

    # Rename the mesh coordinate names, if necessary
    # ---- Longitude
    mesh_longitude = (
        [ col for col in mesh_data.columns if 'lon' in col.lower( ) ][ 0 ]
    )  
    # ---- Latitude
    mesh_latitude = (
        [ col for col in mesh_data.columns if 'lat' in col.lower( ) ][ 0 ]
    )
    # ---- Rename the dataframe
    mesh = mesh_data.copy( ).rename( columns = { f"{mesh_longitude}": 'longitude' ,
                                                 f"{mesh_latitude}": 'latitude' } )
    
    # Convert the mesh dataframes to a geodataframe
    mesh_gdf = gpd.GeoDataFrame( mesh ,
                                 geometry = gpd.points_from_xy( mesh[ 'longitude' ] ,
                                                                mesh[ 'latitude' ] ) ,
                                 crs = settings_dict[ 'projection' ] )
    
    # Convert the mesh projection to UTM (m) -- In place
    wgs84_to_utm( mesh_gdf )
    
    # Determine the survey extent by generating the border polygon    
    survey_polygon = transect_extent( transect_data , 
                                      settings_dict[ 'projection' ] , 
                                      num_nearest_transects )
    
    # Find the mesh coordinates that fall within the buffered polygon
    # ---- Convert `grid_buffer` (nmi) to m and add buffer to polygon
    survey_polygon_buffered = survey_polygon.buffer( mesh_buffer_distance * 1852 ) 
    # ---- Inclusion/union filter mask
    within_polygon_mask = mesh_gdf.geometry.within( survey_polygon_buffered )
    # ---- Apply mask to the mesh grid
    mesh_gdf_masked = mesh_gdf[ within_polygon_mask ]

    # Return the masked mesh dataframe
    return mesh_gdf_masked.drop( columns = 'geometry' )

# ----
kriged_mesh = survey.results[ 'kriging' ][ 'mesh_results_df' ].copy( )
spatial_dict = survey.input[ 'spatial' ]
settings_dict = survey.analysis[ 'settings' ][ 'stratified' ]
n_nodes = kriged_mesh.shape[0]
lon = kriged_mesh[ 'longitude' ]
lat = kriged_mesh[ 'latitude' ]
n_transect_per_lat = 5
lat_eq_inc = np.round( lat * n_transect_per_lat + 0.5 ) / n_transect_per_lat
lon_inc=0.05
uniq_lat_eq_inc=np.unique(lat_eq_inc)  # unique equal-spacing transects 
nlat=len(uniq_lat_eq_inc)

lat_inc=round(lat*5+0.5)/5
lon_inc=0.02
uniq_lat_inc=np.unique(lat_inc)
nlat=len(uniq_lat_inc)
m2nmi = 1852

# ---- Partition the mesh into virtual/synthetic transects
kriged_mesh[ 'transect_num' ] = (
    np.round( kriged_mesh[ 'latitude' ] * n_transect_per_lat + 0.5 ) / n_transect_per_lat
)

# ----> reduce
virtual_transects = pd.DataFrame( { 'transect_num': np.arange( 0 , nlat , 1 ) ,
                                    'latitude': np.unique( kriged_mesh[ 'transect_num' ] ) } )
virtual_transects.set_index( 'latitude' , inplace = True )
# Calculate the minimum and maximum longitude, mean latitude, and area of each transect
virtual_transects[ 'longitude_mean' ] = (
    kriged_mesh.groupby( [ 'transect_num' ] )[ 'longitude' ].mean( )
) 
# ---- Min longitude
virtual_transects[ 'longitude_min' ] = (
    kriged_mesh.groupby( [ 'transect_num' ] )[ 'longitude' ].min( )
)
# ---- Max longitude
virtual_transects[ 'longitude_max' ] = (
    kriged_mesh.groupby( [ 'transect_num' ] )[ 'longitude' ].max( )
)

virtual_transects.reset_index( inplace = True )
# ---- Calculate mean distance (nmi)
virtual_transects[ 'transect_distance' ] = virtual_transects.apply( lambda row:
    geopy.distance.distance(
        ( row[ 'latitude' ] , row[ 'longitude_min' ] ) ,
        ( row[ 'latitude' ] , row[ 'longitude_max' ] )
    ).nm , axis = 1 )
# ---- Calculate the area
# --------- Calculate the mean latitudinal increment
d_latitude = np.concatenate( [ [ np.nan ] , np.diff( uniq_lat_inc ) ] )

virtual_transects[ 'd_latitude' ] = np.where(
    ~ virtual_transects[ 'transect_num' ].isin( [ 0 , nlat ] ) ,
    d_latitude ,
    np.nanmean( d_latitude )
)

virtual_transects[ 'mean_spacing' ] = virtual_transects.apply( lambda row:
    geopy.distance.distance(
        ( row[ 'latitude' ] , row[ 'longitude_mean'] ) ,
        ( row[ 'latitude' ] + row[ 'd_latitude' ] , row[ 'longitude_mean'] )
    ).nm , axis = 1 )

# --------- Append area
virtual_transects[ 'transect_area' ] = np.where(
    ~ virtual_transects[ 'transect_num' ].isin( [ 0 , nlat ] ) ,
    virtual_transects[ 'transect_distance' ] * virtual_transects[ 'mean_spacing' ] / 2.0 ,
    virtual_transects[ 'transect_distance' ] * virtual_transects[ 'mean_spacing' ]
)

# Stratify 
if settings_dict[ 'stratum' ].lower() == 'inpfc': 
    # ---- 
    strata = spatial_dict[ 'inpfc_strata_df' ]
elif settings_dict[ 'stratum' ].tolower() == 'ks': 
    # ----
    strata = spatial_dict[ 'geo_strata_df' ]

# stratum column
stratum_col = settings_dict[ 'stratum_name' ]
latitude_bins = latitude_bins = np.concatenate( [ [ -90.0 ] , strata[ 'northlimit_latitude' ] , [ 90.0 ] ] )
virtual_transects[ f"{stratum_col}" ] = pd.cut( virtual_transects[ 'latitude' ] ,
                                                latitude_bins ,
                                                labels = list( strata[ f"{stratum_col}" ] ) + [ 1 ] ,
                                                ordered = False )

virtual_strata_summary = virtual_transects.groupby( [ f"{stratum_col}" ] )[ 'transect_area' ].sum( )

self = survey


stratified_results_dict = self.results['stratified']['kriging']

settings_dict = self.analysis[ 'settings' ][ 'stratified' ]
stratified_results = copy.deepcopy( stratified_results_dict )
stratified_results.update( 
    { 'back_transform': copy.deepcopy( stratified_results['total'] ) } 
    )
stratified_message['total:estimate']