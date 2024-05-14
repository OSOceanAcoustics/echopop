import numpy as np
import pandas as pd

from echopop.survey import Survey
from echopop.computation.kriging_methods import kriging_interpolation , compute_kriging_weights , range_index_threshold , compute_kriging_statistics , ordinary_kriging
from echopop.computation.spatial_old import local_search_index , griddify_lag_distances
from echopop.computation.variogram_models import variogram

ss = self.input[ 'biology' ][ 'length_df' ]
ss = ss [ ( ss.species_id == 22500 ) ]
ss[ ss.length < 15 ]
ss[ ss.age == 1 ]
ss.groupby( [ 'stratum_num' ] )[ 'haul_weight' ].sum( )

at = specimen_data[ specimen_data.stratum_num == 7 ]
at[ ( at.sex != 'unsexed' ) ]
distributions_dict = self.input[ 'biology' ][ 'distributions' ]
proportions_dict = self.analysis[ 'biology' ][ 'proportions' ]
TS_L_parameters = self.config[ 'TS_length_regression_parameters' ][ 'pacific_hake' ]
settings_dict[ 'transect' ][ 'stratum_name' ]
sigma_bs_strata

# Get stratum column name
stratum_col = settings_dict[ 'transect' ][ 'stratum_name' ]

# Calculate the age1 number proportion
# ---- Initialize the dataframe
age_proportions = proportions_dict[ 'number' ][ 'aged_length_proportions_df' ]
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
age1_proportions = age_proportions_table.sum( )[ : , 1 ]

# Calculate the new length-averaged sigma_bs for each stratum
# ---- Extract the length-binned values
length_bins = distributions_dict[ 'length_bins_df' ]
# ---- Square the length values
length_bins[ 'length_sq' ] = length_bins[ 'length_bins' ] ** 2.0
# ---- Multiply by the TS-length regression coefficient (in the linear domain)
length_bins[ 'length_sq' ] = length_bins[ 'length_sq' ] * to_linear( TS_L_parameters[ 'TS_L_intercept' ] )
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
age1_nasc_proportion = age1_sigma_bs / updated_sigma_bs

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
    unage_proportions_table.iloc[ np.arange( min_index , max_index ) , : ] * 
    age_proportions_age1_table.iloc[ np.arange( min_index , max_index ) , : ]
).sum( )
# ---- Just aged length distribution (age-1)
age1_specific_length_distribution_threshold = age_proportions_age1_table.sum( )


age_weight_proportions = proportions_dict[ 'weight' ][ 'aged_weight_proportions_df' ]
age_weight_proportions = age_weight_proportions[ age_weight_proportions[ stratum_col ] != 0 ]
age_weight_proportions_table = age_weight_proportions.pivot_table( index = [ 'length_bin' ] ,
                                                                   columns = [ 'age_bin' , stratum_col ] ,
                                                                   values = 'weight_proportion_aged' ,
                                                                   aggfunc = 'sum' ,
                                                                   observed = False )
age_weight_proportions_table[ 1 ].sum( ) 

age_weight_proportions_repivot = age_weight_proportions.pivot_table( index = [ 'age_bin' ] ,
                                                                     columns = [ stratum_col ] ,
                                                                     values = 'weight_proportion_aged' ,
                                                                     aggfunc = 'sum' ,
                                                                     observed = False )

age1_weight_proportions = np.where( ( age1_length_distribution_threshold <= 1e-10 ) & 
                                    ( age1_specific_length_distribution_threshold <= 1e-10 ) ,
                                    0.0 ,
                                    age_weight_proportions_table[ 1 ].sum( )
                                    / age_weight_proportions_repivot.sum( ) )


age_weight_proportions_table[ 1 ].sum( ) / age_weight_proportions_repivot.sum( )

age_weight_proportions_table.T.sum( axis = 1 )
age_proportions_table.sum( ).loc[ 5 ].sum( )

age_proportions_age1_table.sum( axis = 0 )

unage_proportions_table.loc[ : , : ] * age_proportions_age1_table.T.loc[ : , : ]

unage_proportions_table.iloc[ np.arange( min_index , max_index ) , : ] * age_proportions_age1_table.T.iloc[ : , np.arange( min_index , max_index ) ]
age_proportions_age1_table.T.iloc[ 5 , np.arange( min_index , max_index ) ]

unage_proportions_table.iloc[ np.arange( min_index , max_index ) , 5 ]
updated_index = age_proportions_age1_table.index.union( age_proportions_age1_table.columns )
age_proportions_age1_table.reindex( index = updated_index , columns = [ 'length_bins' ] , fill_value = 0)
unage_proportions_table.iloc[ 4 , np.arange( min_index , max_index ) ] * age_proportions_age1_table.iloc[ np.arange( min_index , max_index ) , 5 ]

age_proportions_table.loc[ : , 5 ]
age_weight_proportions = proportions_dict[ 'weight' ][ 'aged_weight_proportions_df' ]
age_weight_proportions_table = age_weight_proportions.pivot_table( index = [ 'length_bin' ] ,
                                                                   columns = [ stratum_col , 'age_bin' ] ,
                                                                   values = 'weight_proportion_aged' ,
                                                                   aggfunc = 'sum' ,
                                                                   observed = False )

min_index = np.where( length_bins[ 'length_bins' ] == 10.0 )[ 0 ]
max_index = length_bins[ 'length_bins' ].size
np.arange( min_index , max_index ).tolist()
age_proportions_table.loc[ : , ]
np.arange( min_index , max_index )[1]
age_weight_proportions_table.iloc[ np.arange( min_index , max_index ) , 5 ] * age_proportions_age1_table.iloc[ np.arange( min_index , max_index ) , 5 ]
age_weight_proportions.groupby( [ stratum_col , 'length_bin' , 'age_bin' ] )[ 'weight_proportion_aged' ].sum( )

to_dB( age1_sigma_bs )
age1_proportions
( length_bins[ 'length_sq' ] * age1_proportions[ 1 ] ).sum( )
# Return output
# ---- Initialize dataframe
new_sigma_bs = pd.DataFrame( { f"{stratum_col}": np.unique( age_proportions[ stratum_col ] ) } )
# ---- Add age-1 proportions column
new_sigma_bs[ 'age1_proportion' ] = age1_proportions
# ---- Add updated sigma_bs
new_sigma_bs[ 'sigma_bs_mean' ] = updated_sigma_bs
# ---- Return
# return new_sigma_bs

age_proportions_table.loc[ : , 7 ]

length_bins[ 'length_sq' ].values.dot( age_proportions_alt_table )
age_proportions_alt_table.sum( )

number_dist_lambda = lambda x: x.droplevel( 0 ,axis = 1 ).dot( b )

tt_t.groupby( level = 0 , axis = 1 )

tt_t.droplevel( 0 , axis = 1 )

tt_t.sum )

length_bins[ 'length_sq' ].values.dot( tt_t.sum( ) )

a = to_linear( TS_L_parameters[ 'TS_L_intercept' ] ) * np.sort( np.unique( age_proportions[ 'length_mid' ] ) ) ** 2.0
a.dot( tt_t.loc[ : , 1 ] ).sum( )

a.T * tt_t.loc[ : , 1 ].sum( )

to_linear( TS_L_parameters[ 'TS_L_intercept' ] )

tt = tt[ tt.sex == 'all' ]
tt_t = tt.pivot_table( index = [ 'length_bin' ] , columns = [ 'stratum_num' , 'age_bin' ] , values = 'proportion_number_aged' , aggfunc = 'sum' , observed = False )
tt_t.sum( )

tt_t.loc[ : , 3 ].sum()

denom = aged_number_proportion[ aged_number_proportion.sex == 'all' ].groupby( [ 'stratum_num' ] )[ 'proportion_number_aged' ].sum( )
numer = aged_number_proportion[ ( aged_number_proportion.sex == 'all' ) & ( aged_number_proportion.age_bin == pd.Interval( left = 0.5 , right = 1.5 ) ) ].groupby( [ 'stratum_num' ] )[ 'proportion_number_aged' ].sum( )
age1_wgt_prop = numer / denom
sigma_bs = self.analysis[ 'acoustics' ][ 'sigma_bs' ][0]['strata_mean_df'].loc[ lambda x: x.stratum_num == 5.0 ]
ll = self.input[ 'biology' ][ 'distributions' ][ 'length_bins_df' ][ 'length_bins' ]
vals = pd.DataFrame( { 'ss': to_linear( -68.0 ) * ll**2 } )
pp = aged_number_proportion[ aged_number_proportion.sex == 'all' ]

pt = pp[ pp.stratum_num == 5 ].pivot_table( index = [ 'length_bin' ] , columns = [ 'age_bin' ] , values = [ 'proportion_number_aged' ] , observed = False )

vals.T.values.dot( pt[ 'proportion_number_aged' ] ).sum( )

vals.merge( pp[ pp.stratum_num == 5 ] , how = 'cross' ).assign( ps = lambda x: x.ss * x.proportion_number_aged )[ 'ps' ]

self.input['acoustics']['nasc_df'][ self.input['acoustics']['nasc_df'].haul_num == 205 ]

self.input['acoustics']['nasc_df'].iloc[ 2358 ][ 'NASC_all_ages' ] * ( 1 - age1_wgt_prop.loc[ 5 ] )
survey = Survey( "./config_files/initialization_config.yml" ,
                 "./config_files/survey_year_2019_config.yml" )
species_id = 22500
self = survey
self.acoustics['sigma_bs'] = {}
self.biology['weight'] = {}
self.biology['population'] = {}
self.statistics['length_weight'] = {}
# self.transect_analysis( )
self.strata_mean_sigma_bs( species_id )  
self.fit_binned_length_weight_relationship( species_id )
# self.strata_sex_weight_proportions( species_id )
# self.strata_age_binned_weight_proportions( species_id )
# self.nasc_to_biomass_conversion( species_id )
##################################################
# spec = self.biology[ 'specimen_df' ][ ( self.biology[ 'specimen_df' ].species_id == species_id ) ]
spec = self.input[ 'biology' ][ 'specimen_df' ]
spec = spec[ ( spec.species_id == species_id ) ]
spec_filter = spec[ ( spec.sex != 'unsexed' ) ].dropna( subset = [ 'length' , 'weight' ] )
leng = self.biology[ 'length_df' ][ self.biology[ 'length_df' ].species_id == species_id ]

length_distribution = self.biology[ 'distributions' ][ 'length' ][ 'length_interval_arr' ]
length_intervals = length_distribution
age_distribution = self.biology[ 'distributions' ][ 'age' ][ 'age_interval_arr' ]
age_intervals = age_distribution
### Binned counts 
# Specimen

specimen_len_age = spec.bin_variable( length_distribution , 'length' ).bin_variable( age_distribution , 'age' ).count_variable( contrasts = [ 'stratum_num' , 'length_bin' , 'age_bin' ] , variable = 'length' , fun = 'size' )
specimen_len_age_filter = spec_filter.bin_variable( length_distribution , 'length' ).bin_variable( age_distribution , 'age' ).count_variable( contrasts = [ 'stratum_num' , 'sex' , 'length_bin' , 'age_bin' ] , variable = 'length' , fun = 'size' )
specimen_len_age_filter_full = pd.concat( [ specimen_len_age_filter ,
                                            specimen_len_age_filter.assign( sex = 'all' ).groupby( [ 'stratum_num' , 'sex' , 'length_bin' , 'age_bin' ] , observed = False )[ 'count' ].sum( ).reset_index( ) ] )
specimen_total = specimen_len_age.groupby( [ 'stratum_num' ] )[ 'count' ].sum( ).reset_index( name = 'specimen_total' )
specimen_filter_total = specimen_len_age_filter_full.groupby( [ 'stratum_num' , 'sex' ] )[ 'count' ].sum( ).reset_index( name = 'specimen_total' )
# Length 
length_len = leng.bin_variable( length_distribution , 'length' ).count_variable( contrasts = [ 'stratum_num' , 'sex' , 'length_bin'] , variable = 'length_count' , fun = 'sum' )
length_len_agg = pd.concat( [ length_len ,
                              length_len.assign( sex = 'all' ).groupby( [ 'stratum_num' , 'sex' , 'length_bin' ] , observed = False )[ 'count' ].sum( ).reset_index( ) ] )
length_sexed_total = length_len_agg.groupby( [ 'sex' , 'stratum_num' ] )[ 'count' ].sum( ).reset_index( name = 'length_total' )

aged_data = spec
unaged_data = leng

(
    aged_length_proportions ,
    unaged_length_proportions ,
    sex_proportions ,
    age_proportions
) = transect_number_proportions( spec , leng )

### Import length-weight relationship
length_weight_df = self.statistics[ 'length_weight' ][ 'length_weight_df' ]
### Extract length-weight regression results calculated for all animals
length_weight_df = length_weight_df[ length_weight_df.sex != 'all' ][ [ 'length_bin' , 'sex' , 'weight_modeled' ] ]
fitted_weight = self.statistics[ 'length_weight' ][ 'length_weight_df' ][[ 'length_bin' , 'sex' , 'weight_modeled' ] ]
unaged_catch_weights = self.biology[ 'catch_df' ][ self.biology[ 'catch_df' ].species_id == species_id ]

stratum_weights = transect_stratum_weights( aged_length_proportions ,
                                            unaged_length_proportions ,
                                            sex_proportions , 
                                            fitted_weight )

(
    aged_weight_proportions ,
    unaged_sex_weight_proportions
) = transect_weight_proportions( spec , unaged_catch_weights , unaged_length_proportions , length_weight_df )

acoustics_dict = self.input['acoustics']
biology_dict = self.biology
from echopop.computation.spatial_old import correct_transect_intervals
nasc_interval_df = correct_transect_intervals( acoustics_dict[ 'nasc_df' ] )
sigma_bs_strata = acoustics_dict[ 'sigma_bs' ][ 'strata_mean' ]
info_strata = self.spatial[ 'strata_df' ].copy( )
nasc_fraction_total_df = nasc_interval_df.merge( info_strata , on = [ 'stratum_num' , 'haul_num' ] , how = 'outer' ).merge( sigma_bs_strata , on = [ 'stratum_num' ] , how = 'outer' ).dropna( subset = 'transect_num' )
nasc_fraction_total_df[ 'fraction_hake' ] = nasc_fraction_total_df[ 'fraction_hake' ].fillna( 0.0 )
nasc_fraction_total_df[ 'rho_number_all' ] = (
   nasc_fraction_total_df.fraction_hake
   * nasc_fraction_total_df.NASC_all_ages
   / ( 4.0 * np.pi * nasc_fraction_total_df.sigma_bs_mean )
)
nasc_fraction_total_df[ 'rho_number_adult' ] = (
   nasc_fraction_total_df.fraction_hake
   * nasc_fraction_total_df.NASC_no_age1
   / ( 4.0 * np.pi * nasc_fraction_total_df.sigma_bs_mean )
)
# ---- Abundance
nasc_fraction_total_df[ 'abundance_adult' ] = (
    nasc_fraction_total_df.rho_number_adult
    * nasc_fraction_total_df.interval_area
)

# Apportioned abundance
abundance_conv_df = nasc_fraction_total_df[ [ 'latitude' , 'longitude' , 'transect_num' , 'stratum_num' , 'interval_area' , 'rho_number_adult' , 'abundance_adult' ] ]
abundance_conv_df = abundance_conv_df.merge( sex_proportions[ sex_proportions.sex.isin( [ 'male' , 'female' ] ) ][ [ 'stratum_num' , 'sex' , 'proportion_number_overall' ] ] , on = [ 'stratum_num' ] )
abundance_conv_df[ 'abundance_sex_adult' ] = abundance_conv_df[ 'abundance_adult' ] * abundance_conv_df[ 'proportion_number_overall' ]
biomass_conv_df = abundance_conv_df.merge( stratum_weights[ stratum_weights.sex != 'all' ] )
biomass_conv_df[ 'rho_biomass' ] = biomass_conv_df.average_weight * biomass_conv_df.rho_number_adult * biomass_conv_df.proportion_number_overall
biomass_conv_df[ 'biomass' ] = biomass_conv_df.average_weight * biomass_conv_df.abundance_sex_adult
cross_tabs = biomass_conv_df.merge( aged_weight_proportions[ aged_weight_proportions.sex != 'all' ] , on = [ 'stratum_num' , 'sex' ] )
cross_tabs[ 'sub_biomass' ] = cross_tabs[ 'biomass' ] * cross_tabs[ 'weight_aged_sex_proportion' ]

cross_tabs_table = cross_tabs.pivot_table( columns = 'age_bin' , index = 'length_bin' , values = 'sub_biomass' , aggfunc = 'sum' )
cross_tabs_table.iloc[ : , 1 : 25 ].unstack().reset_index( name = 'biom' )['biom'].sum() * 1e-6

cross_tabs_table.iloc[ : , 1 : 25 ].sum( ).sum( ) * 1e-6
cross_tabs_table.iloc[ : , 0 ].sum( ).sum( ) * 1e-6
cross_tabs[ 'sub_biomass' ].sum( ) * 1e-6


biomass_conv_df[ 'biomass' ].sum( ) * 1e-6

biomass_conv_df[ 'rho_biomass_adult' ].max( ) * 1e-6
biomass_conv_df[ 'rho_biomass_all' ].max( ) * 1e-6
biomass_conv_df[ 'biomass_adult' ] = biomass_conv_df.average_weight * biomass_conv_df.abundance_sex_adult
biomass_conv_df[ 'biomass_all' ] = biomass_conv_df.average_weight * biomass_conv_df.abundance_sex_all
biomass_conv_df[ 'biomass_adult' ].sum( ) * 1e-6
biomass_conv_df[ 'biomass_all' ].sum( ) * 1e-6


biomass_conv_df[ 'biomass_all' ].sum( ) * 1e-6
biomass_conv_df[ 'biomass_adult' ].sum( ) * 1e-6
biomass_conv_df[ biomass_conv_df.sex == 'male' ][ 'biomass_all' ].sum( ) * 1e-6
biomass_conv_df[ biomass_conv_df.sex == 'female' ][ 'biomass_adult' ].sum( ) * 1e-6

test = wgt_proportion[ wgt_proportion.sex == 'all' ].drop( 'sex' , axis = 1 )
test_table = test.pivot_table( columns = [ 'age_bin' , 'stratum_num' ] , index = [ 'length_bin' ] , values = 'proportion' , observed = False )
new_out = biomass_conv_df.merge( test , on = [ 'stratum_num' ] )
new_out[ 'apportioned_biomass_all' ] = (
    new_out.biomass_all * new_out.proportion
)
new_out[ 'apportioned_biomass_adult' ] = (
    new_out.biomass_adult * new_out.proportion
)
new_out[ 'rho_biomass_all' ] = (
    new_out[ 'apportioned_biomass_all' ] 
    / new_out[ 'interval_area' ]
)
new_out[ 'rho_biomass_adult' ] = (
    new_out[ 'apportioned_biomass_adult' ] 
    / new_out[ 'interval_area' ]
)

a = new_out.groupby( [ 'latitude' , 'longitude' , 'transect_num' , 'stratum_num' ] )[ 'apportioned_biomass_adult' ].sum( ).reset_index( )
a[ a.transect_num == 33.0 ]
a[ np.round( a.longitude , 4 ) == -1.242025e+02 ]
full_out = new_out
full_out[ 'age' ] = full_out[ 'age_bin' ].apply( lambda x: x.mid )

test = full_out[ ( full_out.age != 1 ) ]
test['apportioned_biomass_adult' ].sum( ) * 1e-6
out = test.groupby( [ 'transect_num' , 'stratum_num' , 'longitude' , 'latitude' ] )[ 'apportioned_biomass_adult' ].sum( ).reset_index( name = 'biomass' )
out.to_csv( "C:/Users/Brandyn/Documents/echopop_biomass.csv" , index = False )



nasc_interval_df.iloc[ 2358 ]
out.iloc[ 2360 ]
"""For comparison, my code yields 1643.1868 kmt, while Chu's code yields 1643.1898 kmt."""


full_out = new_out.groupby( [ 'length_bin' , 'age' ] , observed = False )[ [ 'apportioned_biomass_adult' ] ].sum( ).reset_index( )
full_out_table = full_out.pivot_table( columns = [ 'age' ] , index = [ 'length_bin' ] , values = 'apportioned_biomass_adult' , observed = False )
full_out_table.sum( )
full_out_table.iloc[ : , 1 : 22 ].sum( ).sum( ) * 1e-6

tt1 = wgt_proportion.pivot_table( index = 'length_bin' , columns = [ 'stratum_num' , 'sex' , 'age_bin' ] , values = 'proportion' )
tt2 = aged_weight_proportions.pivot_table( index = 'age_bin' , columns = [ 'stratum_num' , 'sex' , 'length_bin' ] , values = 'proportion_weight' , aggfunc = 'sum' )

tt1[ 3 , 'all' ].sum( )
tt2[ 3 , 'male' ].sum( )

weight_male = fitted_weight_table.loc[ 'all' ][ 'weight_modeled' ].values.dot( ( length_props_table.loc[ 1 , 'male' ] * station_props_table.loc[ 1 , 'male' ] + length_props_table.loc[ 2 , 'male' ] * station_props_table.loc[ 2 , 'male' ] ) )
weight_female = fitted_weight_table.loc[ 'all' ][ 'weight_modeled' ].values.dot( ( length_props_table.loc[ 1 , 'female' ] * station_props_table.loc[ 1 , 'female' ] + length_props_table.loc[ 2 , 'female' ] * station_props_table.loc[ 2 , 'female' ] ) )
weight_all = fitted_weight_table.loc[ 'all' ][ 'weight_modeled' ].values.dot( ( length_props_table.loc[ 1 , 'all' ] * station_props_table.loc[ 1 , 'all' ] + length_props_table.loc[ 2 , 'all' ] * station_props_table.loc[ 2 , 'all' ] ) )

stratum_weight_df = pd.DataFrame(
    {
        'stratum_num': np.tile( np.unique( station_props.stratum_num ) , 3 ) ,
        'sex': np.repeat( [ 'all' , 'male' , 'female' ] , len( np.unique( station_props.stratum_num ) ) ) ,
        'average_weight': np.concatenate( [ weight_all , weight_male , weight_female ] ) ,
    } ,
)
    
    length_proportions_table.loc[ ( 'aged' , 'male' ) , : ]
    aged_data[ [ 'stratum_num' , 'sex' , 'length_bin' ]].assign( group = 'aged' )
    unaged_data
    length_props = pd.concat( [ specimen_full_proportion.assign( station = 2 ) ,
                            length_len_proportion[ length_len_proportion.sex != 'unsexed' ][ [ 'stratum_num' , 'sex' , 'length_bin' , 'number_proportion_bin' ] ].assign( station = 1 ) ] )

length_props_table = length_props.pivot_table( index = [ 'station' , 'sex' , 'length_bin' ] ,
                                               columns = [ 'stratum_num' ] ,
                                               values = 'number_proportion_bin' ,
                                               observed = False ).fillna( 0.0 )
length_props_table.loc[ ( 1 , 'all' ) , : ]
tt = aged_length_proportions.pivot_table( index = [ 'length_bin' ] , columns = [ 'stratum_num' , 'sex' , 'age_bin' ] , values = 'proportion_number_overall_aged' )
tt[ ( 8 , 'female' ) ]

age_length_proportions_table.sum( )
specimen_full_proportion_table.sum( )

age_length_proportions = aged_length_proportions.groupby( [ 'stratum_num' , 'sex' , 'length_bin' ] , observed = False )[ 'proportion_number_sex_aged' ].sum( ).reset_index( )
aged_length_proportions.groupby( [ 'stratum_num' , 'sex' , 'length_bin' ] )[ 'proportion_number_sex_aged' ].sum( )
specimen_full_proportion = specimen_len_age_proportion[ specimen_len_age_proportion.sex != 'unsexed' ].groupby( [ 'stratum_num' , 'sex' , 'length_bin' ] , observed = False )[ 'number_proportion_bin' ].sum( ).reset_index( )

# specimen_full_proportion = specimen_len_age_proportion[ specimen_len_age_proportion.sex == 'all' ].groupby( [ 'stratum_num' , 'length_bin' ] , observed = False )[ 'number_proportion_bin' ].sum( ).reset_index( )

specimen_full_proportion_table = specimen_full_proportion.pivot_table( index = [ 'length_bin' ] ,
                                                                       columns = [ 'stratum_num' ] ,
                                                                       values = 'number_proportion_bin' ,
                                                                       observed = False )
station_props = pd.concat( [ length_proportion_all[ [ 'stratum_num' , 'fac1' ] ].rename( columns = { 'fac1': 'proportion' } ).assign( station = 1 ) ,
                             length_proportion_all[ [ 'stratum_num' , 'fac2' ] ].rename( columns = { 'fac2': 'proportion' } ).assign( station = 2 ) ] )

length_props = pd.concat( [ specimen_full_proportion.assign( station = 2 ) ,
                            length_len_proportion[ length_len_proportion.sex != 'unsexed' ][ [ 'stratum_num' , 'sex' , 'length_bin' , 'number_proportion_bin' ] ].assign( station = 1 ) ] )

length_props_table = length_props.pivot_table( index = [ 'station' , 'sex' , 'length_bin' ] ,
                                               columns = [ 'stratum_num' ] ,
                                               values = 'number_proportion_bin' ,
                                               observed = False ).fillna( 0.0 )

fitted_weight = self.statistics[ 'length_weight' ][ 'length_weight_df' ][[ 'length_bin' , 'sex' , 'weight_modeled' ] ]
fitted_weight_table = fitted_weight.pivot_table( index = [ 'sex' , 'length_bin' ] ,
                                                 values = 'weight_modeled' ,
                                                 observed = False )

fitted_weight_table.loc[ 'all' ][ 'weight_modeled' ].values.dot( length_props_table.loc[ 1 , 'male' ] * station_props_table.loc[ 1 , 'male' ] + length_props_table.loc[ 2 , 'male' ] * station_props_table.loc[ 2 , 'male' ] )

weight_male = fitted_weight_table.loc[ 'all' ][ 'weight_modeled' ].values.dot( ( length_props_table.loc[ 1 , 'male' ] * station_props_table.loc[ 1 , 'male' ] + length_props_table.loc[ 2 , 'male' ] * station_props_table.loc[ 2 , 'male' ] ) )
weight_female = fitted_weight_table.loc[ 'all' ][ 'weight_modeled' ].values.dot( ( length_props_table.loc[ 1 , 'female' ] * station_props_table.loc[ 1 , 'female' ] + length_props_table.loc[ 2 , 'female' ] * station_props_table.loc[ 2 , 'female' ] ) )
weight_all = fitted_weight_table.loc[ 'all' ][ 'weight_modeled' ].values.dot( ( length_props_table.loc[ 1 , 'all' ] * station_props_table.loc[ 1 , 'all' ] + length_props_table.loc[ 2 , 'all' ] * station_props_table.loc[ 2 , 'all' ] ) )

stratum_weight_df = pd.DataFrame(
    {
        'stratum_num': np.tile( np.unique( station_props.stratum_num ) , 3 ) ,
        'sex': np.repeat( [ 'all' , 'male' , 'female' ] , len( np.unique( station_props.stratum_num ) ) ) ,
        'average_weight': np.concatenate( [ weight_all , weight_male , weight_female ] ) ,
    } ,
)

station_props = pd.concat( [ length_proportion_all[ [ 'stratum_num' , 'fac1' , 'fac2' ] ].assign( sex = 'all' ) ,
                            length_proportion_full[ [ 'stratum_num' , 'sex' , 'fac1' , 'fac2' ] ] ] )

station_props = pd.wide_to_long( station_props , stubnames = 'fac' , i = [ 'stratum_num' , 'sex' ] , j = 'station' ).reset_index( ).rename( columns = { 'fac': 'proportion' } )

station_props_table = station_props.pivot_table( index = [ 'station' , 'sex' ] ,
                                                 columns = [ 'stratum_num' ] ,
                                                 values = 'proportion' ).fillna( 0.0 )
length_proportion_all = specimen_proportion_all[ specimen_proportion_all.sex == 'all' ].assign( number_proportion = lambda x: 1.0 - x.specimen_number_proportion_overall ).drop( [ 'sex' ] , axis = 1 )
summed_length_bin_proportion = length_len_proportion[ length_len_proportion.sex.isin( [ 'male' , 'female' ] ) ].groupby( [ 'stratum_num' , 'sex' ] )[ 'number_proportion_bin' ].sum( ).reset_index( )
length_proportion_full = length_proportion_all[ [ 'stratum_num' , 'number_proportion' ] ].merge( summed_length_bin_proportion , on = [ 'stratum_num' ] , how = 'outer' ).dropna( subset = [ 'number_proportion' , 'number_proportion_bin' ] )
length_proportion_full[ 'fac1' ] = (
   length_proportion_full.number_proportion * length_proportion_full.number_proportion_bin 
)
length_proportion_full = length_proportion_full.merge( specimen_sex_proportion_agg , on = [ 'stratum_num' , 'sex' ] )
length_proportion_full[ 'fac1' ] = (
    length_proportion_full.fac1  / (
        length_proportion_full.specimen_number_proportion_overall
        + length_proportion_full.fac1
    )
)
length_proportion_all[ 'fac1' ] = (
    length_proportion_all.number_proportion / (
        length_proportion_all.number_proportion
        + length_proportion_all.specimen_number_proportion_overall 
    )
) 
length_proportion_full[ 'fac2' ] = (
    length_proportion_full.specimen_number_proportion_overall  / (
        length_proportion_full.specimen_number_proportion_overall
        + length_proportion_full.fac1
    )
)
length_proportion_all[ 'fac2' ] = (
    length_proportion_all.specimen_number_proportion_overall  / (
        length_proportion_all.fac1
        + length_proportion_all.specimen_number_proportion_overall 
    )
) 

station_props = pd.concat( [ length_proportion_all[ [ 'stratum_num' , 'fac1' , 'fac2' ] ].assign( sex = 'all' ) ,
                             length_proportion_full[ [ 'stratum_num' , 'sex' , 'fac1' , 'fac2' ] ] ] )

station_props = pd.wide_to_long( station_props , stubnames = 'fac' , i = [ 'stratum_num' , 'sex' ] , j = 'station' ).reset_index( ).rename( columns = { 'fac': 'proportion' } )

station_props_table = station_props.pivot_table( index = [ 'station' , 'sex' ] ,
                                                 columns = [ 'stratum_num' ] ,
                                                 values = 'proportion' ).fillna( 0.0 )

aged_number_proportion.groupby( [ 'stratum_num' , 'age_bin' ] )[ 'proportion_number_aged' ].sum( )
aged_number_distribution_filtered.groupby( 'stratum_num' )[ 'count' ].sum( )


### Totals
grand_total = specimen_filter_total[ specimen_filter_total.sex == 'all' ].merge( length_sexed_total[ length_sexed_total.sex == 'all' ] , on = [ 'stratum_num' , 'sex' ] , how = 'outer' ).fillna( 0.0 )
grand_total[ 'overall_total' ] = grand_total.specimen_total + grand_total.length_total
sub_total = specimen_filter_total.merge( length_sexed_total[ length_sexed_total.sex == 'all' ] , on = [ 'stratum_num' ] , how = 'outer' ).fillna( 0.0 )
sub_total[ 'filtered_total' ] = sub_total.specimen_total + sub_total.length_total

### Proportion key
# Specimen
specimen_len_age_proportion = specimen_len_age_filter_full.merge( specimen_filter_total , on = [ 'stratum_num' , 'sex' ] , how = 'outer' )
specimen_len_age_proportion[ 'number_proportion_bin' ] = (
    specimen_len_age_proportion[ 'count' ] / specimen_len_age_proportion[ 'specimen_total' ]
)

specimen_len_age_proportion_all = specimen_len_age_filter_full[ specimen_len_age_filter_full.sex != 'all' ].drop( 'sex' , axis = 1 ).merge( specimen_filter_total[ specimen_filter_total.sex == 'all' ] , on = [ 'stratum_num'  ] )
specimen_len_age_proportion_all[ 'number_proportion_all' ] = (
    specimen_len_age_proportion_all[ 'count' ]
    / specimen_len_age_proportion_all[ 'specimen_total' ]
)

specimen_sexed_proportions = specimen_filter_total[ specimen_filter_total.sex != 'all' ].merge( grand_total.drop( 'sex' , axis = 1 ) , on = [ 'stratum_num' ] , how = 'outer' )
specimen_sexed_proportions[ 'sex_proportion' ] = (
    specimen_sexed_proportions[ 'specimen_total_x' ] / specimen_sexed_proportions.specimen_total_y
)

# Length
length_len_proportion = length_len_agg.merge( length_sexed_total[ [ 'stratum_num' , 'sex' , 'length_total' ] ] , on = [ 'stratum_num' , 'sex' ] , how = 'outer' ).dropna( subset = [ 'count' ] )
length_len_proportion[ 'number_proportion_bin' ] = (
    length_len_proportion[ 'count' ] / length_len_proportion[ 'length_total' ]
).fillna( 0.0 )
length_len_proportion_overall = length_len_agg.merge( grand_total[ [ 'stratum_num' , 'overall_total' ] ] , on = [ 'stratum_num' ] , how = 'outer' ).dropna( subset = [ 'count' ] )
length_len_proportion_overall[ 'number_proportion' ] = (
    length_len_proportion_overall[ 'count' ]
    / length_len_proportion_overall[ 'overall_total' ]
)
### Sex number proportions
# Specimen
specimen_sex_proportion = specimen_len_age_filter_full[ specimen_len_age_filter_full.sex != 'all' ].merge( grand_total[ [ 'stratum_num' , 'overall_total' ] ] , on = [ 'stratum_num' ] , how = 'outer' )
specimen_sex_proportion[ 'number_proportion' ] = (
    specimen_sex_proportion[ 'count' ] / specimen_sex_proportion[ 'overall_total' ]
).fillna( 0.0 )
# ----> Aggregate by stratum
specimen_sex_proportion_agg = specimen_sex_proportion.groupby( [ 'stratum_num' , 'sex' ] )[ 'number_proportion' ].sum( ).reset_index( name = 'specimen_number_proportion_overall' )
### Aggregate
specimen_proportion_all = pd.concat( [ specimen_sex_proportion_agg ,
                                       specimen_sex_proportion_agg.assign( sex = 'all' ).groupby( [ 'stratum_num' , 'sex' ] )[ 'specimen_number_proportion_overall' ].sum( ).reset_index( name = 'specimen_number_proportion_overall' ) ] )
length_sex_proportion_agg = length_len_proportion_overall.groupby( [ 'stratum_num' , 'sex' ] )[ 'number_proportion' ].sum( ).reset_index( name = 'length_number_proportion_overall' )
sexed_number_proportions = specimen_sex_proportion_agg.merge( length_sex_proportion_agg , on = [ 'stratum_num' , 'sex' ] )
sexed_number_proportions[ 'number_proportions' ] = (
    sexed_number_proportions.specimen_number_proportion_overall
    + sexed_number_proportions.length_number_proportion_overall
)
sexed_number_proportions = sexed_number_proportions[ [ 'stratum_num' , 'sex' , 'number_proportions' ] ]

# Station 1
length_proportion_all = specimen_proportion_all[ specimen_proportion_all.sex == 'all' ].assign( number_proportion = lambda x: 1.0 - x.specimen_number_proportion_overall ).drop( [ 'sex' ] , axis = 1 )
summed_length_bin_proportion = length_len_proportion[ length_len_proportion.sex.isin( [ 'male' , 'female' ] ) ].groupby( [ 'stratum_num' , 'sex' ] )[ 'number_proportion_bin' ].sum( ).reset_index( )
length_proportion_full = length_proportion_all[ [ 'stratum_num' , 'number_proportion' ] ].merge( summed_length_bin_proportion , on = [ 'stratum_num' ] , how = 'outer' ).dropna( subset = [ 'number_proportion' , 'number_proportion_bin' ] )
length_proportion_full[ 'fac1' ] = (
   length_proportion_full.number_proportion * length_proportion_full.number_proportion_bin 
)
length_proportion_full = length_proportion_full.merge( specimen_sex_proportion_agg , on = [ 'stratum_num' , 'sex' ] )
length_proportion_full[ 'fac1' ] = (
    length_proportion_full.fac1  / (
        length_proportion_full.specimen_number_proportion_overall
        + length_proportion_full.fac1
    )
)
length_proportion_all[ 'fac1' ] = (
    length_proportion_all.number_proportion / (
        length_proportion_all.number_proportion
        + length_proportion_all.specimen_number_proportion_overall 
    )
) 
length_proportion_full[ 'fac2' ] = (
    length_proportion_full.specimen_number_proportion_overall  / (
        length_proportion_full.specimen_number_proportion_overall
        + length_proportion_full.fac1
    )
)
length_proportion_all[ 'fac2' ] = (
    length_proportion_all.specimen_number_proportion_overall  / (
        length_proportion_all.fac1
        + length_proportion_all.specimen_number_proportion_overall 
    )
) 

station_props = pd.concat( [ length_proportion_all[ [ 'stratum_num' , 'fac1' , 'fac2' ] ].assign( sex = 'all' ) ,
                             length_proportion_full[ [ 'stratum_num' , 'sex' , 'fac1' , 'fac2' ] ] ] )

station_props = pd.wide_to_long( station_props , stubnames = 'fac' , i = [ 'stratum_num' , 'sex' ] , j = 'station' ).reset_index( ).rename( columns = { 'fac': 'proportion' } )

station_props_table = station_props.pivot_table( index = [ 'station' , 'sex' ] ,
                                                 columns = [ 'stratum_num' ] ,
                                                 values = 'proportion' ).fillna( 0.0 )

specimen_full_proportion = specimen_len_age_proportion[ specimen_len_age_proportion.sex != 'unsexed' ].groupby( [ 'stratum_num' , 'sex' , 'length_bin' ] , observed = False )[ 'number_proportion_bin' ].sum( ).reset_index( )

# specimen_full_proportion = specimen_len_age_proportion[ specimen_len_age_proportion.sex == 'all' ].groupby( [ 'stratum_num' , 'length_bin' ] , observed = False )[ 'number_proportion_bin' ].sum( ).reset_index( )

specimen_full_proportion_table = specimen_full_proportion.pivot_table( index = [ 'length_bin' ] ,
                                                                       columns = [ 'stratum_num' ] ,
                                                                       values = 'number_proportion_bin' ,
                                                                       observed = False )
station_props = pd.concat( [ length_proportion_all[ [ 'stratum_num' , 'fac1' ] ].rename( columns = { 'fac1': 'proportion' } ).assign( station = 1 ) ,
                             length_proportion_all[ [ 'stratum_num' , 'fac2' ] ].rename( columns = { 'fac2': 'proportion' } ).assign( station = 2 ) ] )

length_props = pd.concat( [ specimen_full_proportion.assign( station = 2 ) ,
                            length_len_proportion[ length_len_proportion.sex != 'unsexed' ][ [ 'stratum_num' , 'sex' , 'length_bin' , 'number_proportion_bin' ] ].assign( station = 1 ) ] )

length_props_table = length_props.pivot_table( index = [ 'station' , 'sex' , 'length_bin' ] ,
                                               columns = [ 'stratum_num' ] ,
                                               values = 'number_proportion_bin' ,
                                               observed = False ).fillna( 0.0 )

fitted_weight = self.statistics[ 'length_weight' ][ 'length_weight_df' ][[ 'length_bin' , 'sex' , 'weight_modeled' ] ]
fitted_weight_table = fitted_weight.pivot_table( index = [ 'sex' , 'length_bin' ] ,
                                                 values = 'weight_modeled' ,
                                                 observed = False )

fitted_weight_table.loc[ 'all' ][ 'weight_modeled' ].values.dot( length_props_table.loc[ 1 , 'male' ] * station_props_table.loc[ 1 , 'male' ] + length_props_table.loc[ 2 , 'male' ] * station_props_table.loc[ 2 , 'male' ] )

weight_male = fitted_weight_table.loc[ 'all' ][ 'weight_modeled' ].values.dot( ( length_props_table.loc[ 1 , 'male' ] * station_props_table.loc[ 1 , 'male' ] + length_props_table.loc[ 2 , 'male' ] * station_props_table.loc[ 2 , 'male' ] ) )
weight_female = fitted_weight_table.loc[ 'all' ][ 'weight_modeled' ].values.dot( ( length_props_table.loc[ 1 , 'female' ] * station_props_table.loc[ 1 , 'female' ] + length_props_table.loc[ 2 , 'female' ] * station_props_table.loc[ 2 , 'female' ] ) )
weight_all = fitted_weight_table.loc[ 'all' ][ 'weight_modeled' ].values.dot( ( length_props_table.loc[ 1 , 'all' ] * station_props_table.loc[ 1 , 'all' ] + length_props_table.loc[ 2 , 'all' ] * station_props_table.loc[ 2 , 'all' ] ) )

stratum_weight_df = pd.DataFrame(
    {
        'stratum_num': np.tile( np.unique( station_props.stratum_num ) , 3 ) ,
        'sex': np.repeat( [ 'all' , 'male' , 'female' ] , len( np.unique( station_props.stratum_num ) ) ) ,
        'average_weight': np.concatenate( [ weight_all , weight_male , weight_female ] ) ,
    } ,
)
### Weight proportions
# ---- Age-count proportions
specimen_data = spec
full_wgt_sex = spec.bin_variable( length_intervals , 'length' ).bin_variable( age_intervals , 'age' ).groupby( [ 'sex' , 'stratum_num' , 'length_bin' , 'age_bin' ] )[ 'weight' ].sum( ).reset_index( )
full_wgt = full_wgt_sex[ full_wgt_sex.sex != 'unsexed' ].groupby( [ 'stratum_num' , 'sex' ] )[ 'weight' ].sum( ).reset_index( name = 'total_weight' )
full_wgt_all = full_wgt.groupby( [ 'stratum_num' ] )[ 'total_weight' ].sum( ).reset_index( )
full_wgt_full = pd.concat( [ full_wgt , full_wgt_all.assign( sex = 'all' ) ] )
full_wgt_sex_all = full_wgt_sex[ full_wgt_sex.sex != 'unsexed' ].assign( sex = 'all' ).groupby( [ 'stratum_num' , 'length_bin' , 'age_bin' , 'sex' ] )[ 'weight' ].sum( ).reset_index( name = 'weight' )
full_wgt_sex_full = pd.concat( [ full_wgt_sex , full_wgt_sex_all ] )
wgt_proportion = full_wgt_sex_full.merge( full_wgt_full ).assign( proportion = lambda x: x.weight / x.total_weight )[ [ 'stratum_num' , 'sex' , 'length_bin' , 'age_bin' , 'proportion' ] ]

acoustics_dict = self.acoustics
biology_dict = self.biology
from echopop.computation.spatial_old import correct_transect_intervals
np.unique( self.acoustics[ 'nasc' ][ 'nasc_df' ][ 'transect_num' ] )
nasc_interval_df = correct_transect_intervals( acoustics_dict[ 'nasc' ][ 'nasc_df' ] )
sigma_bs_strata = acoustics_dict[ 'sigma_bs' ][ 'strata_mean' ]
info_strata = self.spatial[ 'strata_df' ].copy( )

nasc_fraction_total_df = nasc_interval_df.merge( info_strata , on = [ 'stratum_num' , 'haul_num' ] , how = 'outer' ).merge( sigma_bs_strata , on = [ 'stratum_num' ] , how = 'outer' ).dropna( subset = 'transect_num' )
nasc_fraction_total_df[ 'fraction_hake' ] = nasc_fraction_total_df[ 'fraction_hake' ].fillna( 0.0 )

nasc_fraction_total_df[ 'rho_number_all' ] = (
   nasc_fraction_total_df.fraction_hake
   * nasc_fraction_total_df.NASC_all_ages
   / ( 4.0 * np.pi * nasc_fraction_total_df.sigma_bs_mean )
)
nasc_fraction_total_df[ 'rho_number_adult' ] = (
   nasc_fraction_total_df.fraction_hake
   * nasc_fraction_total_df.NASC_no_age1
   / ( 4.0 * np.pi * nasc_fraction_total_df.sigma_bs_mean )
)

# ---- Abundance
nasc_fraction_total_df[ 'abundance_all' ] = (
    nasc_fraction_total_df.rho_number_all
    * nasc_fraction_total_df.interval_area
)

nasc_fraction_total_df[ 'abundance_adult' ] = (
    nasc_fraction_total_df.rho_number_adult
    * nasc_fraction_total_df.interval_area
)

# Apportioned abundance

abundance_conv_df = nasc_fraction_total_df[ [ 'latitude' , 'longitude' , 'transect_num' , 'stratum_num' , 'interval_area' , 'NASC_no_age1' , 'rho_number_all' , 'rho_number_adult' , 'abundance_all' , 'abundance_adult' ] ]
abundance_conv_df = abundance_conv_df.merge( sexed_number_proportions , on = [ 'stratum_num' ] )
abundance_conv_df[ 'abundance_sex_all' ] = abundance_conv_df[ 'abundance_all' ] * abundance_conv_df[ 'number_proportions' ]
abundance_conv_df[ 'abundance_sex_adult' ] = abundance_conv_df[ 'abundance_adult' ] * abundance_conv_df[ 'number_proportions' ]

biomass_conv_df = abundance_conv_df.merge( stratum_weight_df[ stratum_weight_df.sex != 'all' ] )
biomass_conv_df[ 'rho_biomass_all' ] = biomass_conv_df.average_weight * biomass_conv_df.rho_number_all * biomass_conv_df.number_proportions
biomass_conv_df[ 'rho_biomass_adult' ] = biomass_conv_df.average_weight * biomass_conv_df.rho_number_adult * biomass_conv_df.number_proportions
biomass_conv_df[ 'rho_biomass_adult' ].max( ) * 1e-6
biomass_conv_df[ 'rho_biomass_all' ].max( ) * 1e-6
biomass_conv_df[ 'biomass_adult' ] = biomass_conv_df.average_weight * biomass_conv_df.abundance_sex_adult
biomass_conv_df[ 'biomass_all' ] = biomass_conv_df.average_weight * biomass_conv_df.abundance_sex_all
biomass_conv_df[ 'biomass_adult' ].sum( ) * 1e-6
biomass_conv_df[ 'biomass_all' ].sum( ) * 1e-6

biomass_conv_df[ ( np.round( biomass_conv_df.longitude , 4 ) == -1.242025e+02 ) ]
a = biomass_conv_df.groupby( [ 'latitude' , 'longitude' , 'transect_num' , 'stratum_num' , 'NASC_no_age1' ] )[ 'biomass_adult' ].sum( ).reset_index( )
a[ ( a.transect_num == 33 ) & ( np.round( a.longitude , 4 ) == -1.242025e+02 ) ]

import geopandas as gpd
from shapely.geometry import Polygon
from shapely.ops import unary_union
from echopop.computation.spatial_old import to_utm
from echopop.computation.spatial_old import transform_geometry

dataset = 'biomass_density'
variable = 'rho_biomass_adult'
standardize_coordinates = True 
extrapolate = False
grid_buffer = 1.00
nearest_transects = 4
projection: str = 'epsg:4326'
spatial_data = biomass_conv_df.groupby( [ 'latitude' , 'longitude' , 'transect_num' , 'stratum_num' , 'interval_area' ] )[ 'rho_biomass_adult' ].sum( ).reset_index( )
kriging_parameters = self.statistics[ 'kriging' ][ 'model_config' ].copy( )
variogram_parameters = self.statistics[ 'variogram' ][ 'model_config' ].copy( )
dataframe_mesh = self.statistics[ 'kriging' ][ 'mesh_df' ].copy( )
dataframe_geostrata = self.spatial[ 'geo_strata_df' ].copy( )
reference_grid = self.statistics[ 'kriging' ][ 'isobath_200m_df' ]


lat_b1 = 34.3
dlat = 1.25 / 60.0

inpfc_df = self.spatial[ 'inpfc_strata_df' ].rename( { 'stratum_num': 'stratum_inpfc' } , axis = 1 )
# Bin the strata based on latitude limits 
latitude_bins = np.concatenate( [ [ -90 ] , inpfc_df.northlimit_latitude ] )
tes = self.acoustics[ 'nasc' ][ 'nasc_df' ]
tes[ 'stratum_inpfc' ] = pd.cut( tes[ 'latitude' ] , latitude_bins , right = False , labels = inpfc_df[ 'stratum_inpfc' ] )
tes.groupby( [ 'stratum_inpfc' ] )[ 'latitude' ].max( )
tes[ ( tes.stratum_inpfc == 6 ) ]

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


test[ test.haul_num ]
spatial_data[ spatial_data.transect_num == 119 ][ 'latitude' ].max( )
transect_data = nasc_df
mesh_data = self.input[ 'statistics' ][ 'kriging' ][ 'mesh_df' ]
settings_dict = self.analysis[ 'settings' ][ 'kriging' ]

if extrapolate is False :
    
    lat_col = [ col for col in spatial_data.columns if 'latitude' in col.lower( ) ][ 0 ]
    lon_col = [ col for col in spatial_data.columns if 'longitude' in col.lower( ) ][ 0 ]
    grid_lat_col = [ col for col in dataframe_mesh.columns if 'latitude' in col.lower( ) ][ 0 ]
    grid_lon_col = [ col for col in dataframe_mesh.columns if 'longitude' in col.lower( ) ][ 0 ]
    
    spatial_data_geometry = spatial_data[ [ 'transect_num' , lon_col , lat_col ] ].drop_duplicates( [ lon_col , lat_col ] )

    spatial_data_gdf = gpd.GeoDataFrame(
        spatial_data_geometry ,
        geometry = gpd.points_from_xy( spatial_data_geometry[ lon_col ] ,
                                       spatial_data_geometry[ lat_col ] ) ,
        crs = projection ,      
    )
    
    dataframe_mesh_gdf = gpd.GeoDataFrame(
        dataframe_mesh ,
        geometry = gpd.points_from_xy( dataframe_mesh[ grid_lon_col ] ,
                                       dataframe_mesh[ grid_lat_col ] ) ,
        crs = projection ,       
    )
    
    # Convert to UTM
    utm_code = to_utm( np.median( spatial_data_gdf.longitude ) , np.median( spatial_data_gdf.latitude ) )
    spatial_data_gdf.to_crs( f"epsg:{utm_code}" , inplace = True )
    dataframe_mesh_gdf.to_crs( f"epsg:{utm_code}" , inplace = True )
    
    # Calculate centroid of each transect line 
    coord_centroid = spatial_data_gdf.dissolve( by = 'transect_num' , aggfunc = 'mean' ).reset_index( )
    
    #
    transect_polygons = [ ]

    for transect in coord_centroid[ 'transect_num' ].unique( ) :
        
        # Extract coordinates of specific transect
        transect_centroid = coord_centroid[ coord_centroid.transect_num == transect ][ 'geometry' ]
        other_centroids = coord_centroid[ coord_centroid.transect_num != transect ].copy( )
        
        # Calculate distance
        other_centroids[ 'distance_centroid' ] = other_centroids.geometry.apply( lambda g: transect_centroid.distance( g ) )
        
        # Find closest transects
        closest_centroids = other_centroids[ other_centroids.distance_centroid.isin( other_centroids.distance_centroid.nsmallest( nearest_transects ) ) ]
        
        # Collect the nearest transects
        unique_transects = np.unique( closest_centroids.transect_num.values )
        unique_transects = np.append( unique_transects , transect )
        # Create polygon of transects
        polygon = Polygon( list( spatial_data_gdf[ spatial_data_gdf.transect_num.isin( unique_transects ) ][ 'geometry' ] ) )
        
        # Append the convex hull to list
        transect_polygons.append( polygon.convex_hull )

    # Polygon union to collapse 
    survey_polygon = unary_union( transect_polygons )
    
    # Create distance buffer that limits extrapolation
    # ---- Convert `grid_buffer` (nmi) to m
    survey_polygon_buffered = survey_polygon.buffer( grid_buffer * 1852 )
    
    # Filter mesh grid
    dataframe_mesh_gdf_within = dataframe_mesh_gdf[ dataframe_mesh_gdf[ 'geometry' ].within( survey_polygon_buffered ) ]
    # dataframe_mesh_gdf_outside = dataframe_mesh_gdf[ ~ dataframe_mesh_gdf[ 'geometry' ].within( survey_polygon_buffered ) ]
    dataframe_mesh = dataframe_mesh_gdf_within[ [ 'fraction_cell_in_polygon' , 'centroid_longitude' , 'centroid_latitude' ] ].reset_index( drop = True )
    
    print( "Kriging mesh reduced to limit extrapolation." )

from shapely.ops import cascaded_union
from shapely.ops import nearest_points, snap

from pyproj import Transformer
from shapely.ops import transform


### West flank
spatial_data_gdf = spatial_data
unique_coords = spatial_data_gdf
y_resolution = 0.5
# unique_coords[ 'x' ] = spatial_data_gdf.geometry.x
# unique_coords[ 'y' ] = spatial_data_gdf.geometry.y
unique_coords[ 'x' ] = spatial_data_gdf.longitude
unique_coords[ 'y' ] = spatial_data_gdf.latitude

### GET CENTER OF EACH (X)
unique_coords[ 'x_mean' ]

transect_centers_x = unique_coords.groupby( [ 'transect_num' ] )[ [ 'x' ] ].mean( ).reset_index( )
transect_max_y = unique_coords.groupby( [ 'transect_num' ] )[ [ 'y' ] ].max( ).reset_index( )

transect_y_threshold = transect_max_y[ transect_max_y.y > 51 ].reset_index( )
transect_x_threshold = transect_centers_x[ ( transect_centers_x.transect_num.isin( transect_y_threshold.transect_num ) ) ].reset_index( )
transect_west = transect_x_threshold[ transect_x_threshold.x == transect_x_threshold[ transect_x_threshold.x < -132.5 ].x.min( ) ].reset_index( )
transect_west[ 'y' ] = transect_y_threshold[ transect_y_threshold.transect_num.values == transect_west.transect_num.values ][ 'y' ]
transect_east = transect_x_threshold[ transect_x_threshold.x == transect_x_threshold[ transect_x_threshold.x > -132.5 ].x.max( ) ]
transect_east[ 'y' ] = transect_y_threshold.loc[ transect_y_threshold.transect_num.values == transect_east.transect_num.values ].loc[ : , 'y' ]

### GET NORTH
rat = unique_coords.groupby( [ 'transect_num' ] )[ [ 'x' , 'y' ] ].apply( lambda a: np.abs( np.diff( a.x ) ).mean( ) / np.abs( np.diff( a.y ) ).mean( ) ).reset_index( name = 'ratio' )
trans = rat[ rat.ratio < 1 ]
uu = unique_coords[ unique_coords.transect_num.isin( trans.transect_num ) ]
# ---- Get 
plt.scatter( unique_coords.x , unique_coords.y , s = 0.5 , color = 'black' )
plt.scatter( uu.x , uu.y , s = 0.5 , color = 'red' )
plt.scatter( transect_east.x , transect_east.y , s = 0.5 , color = 'red' )
plt.show( )

# Universal (global) -- Eastern border
east_x = unique_coords.groupby( [ 'transect_num' ] )[ 'x' ].max( ).reset_index( )
west_x = unique_coords.groupby( [ 'transect_num' ] )[ 'x' ].min( ).reset_index( )

pos = unique_coords.set_index( 'transect_num' )
pos[ 'x_min' ] = west_x.set_index( 'transect_num' )
pos[ 'x_max' ] = east_x.set_index( 'transect_num' )
pos = pos.reset_index()

east_pos = pos.copy( )[ pos.x == pos.x_max ]
west_pos = pos.copy( )[ pos.x == pos.x_min ]
inter = np.linspace( pos[ 'y' ].min( ) - y_resolution ,
                     pos[ 'y' ].max( ) + y_resolution ,
                     ( ( ( pos[ 'y' ].max( ) + y_resolution ) - ( pos[ 'y' ].min( ) - y_resolution ) ) // ( y_resolution ) ).astype( int ) )

east_pos[ 'y_rounded' ] = pd.cut( east_pos[ 'y' ] , inter )
west_pos[ 'y_rounded' ] = pd.cut( west_pos[ 'y' ] , inter )

east_pos[ 'x_east' ] = east_pos.groupby( [ 'y_rounded' ] , observed = True )[ 'x_max' ].transform( 'max' )
west_pos[ 'x_west' ] = west_pos.groupby( [ 'y_rounded' ] , observed = True )[ 'x_min' ].transform( 'min' )

east_pos = east_pos[ np.abs( east_pos.x - east_pos.x_east ) < 1 ]
# west_pos = west_pos[ ~ west_pos.transect_num.isin( east_pos.transect_num ) ]
west_pos = west_pos[ np.abs( west_pos.x - west_pos.x_west ) < 0.5 ]
west_pos = west_pos.sort_values( [ 'y' ] )
y_array_e = np.arange( unique_coords.y.min( ) , unique_coords.y.max( ) , y_resolution / 60.0 )
y_array_w = np.arange( unique_coords.y.min( ) , unique_coords.y.max( ) , y_resolution / 60.0 )
# y_array_w = np.arange( west_pos.y.min( ) , west_pos.y.max( ) , y_resolution / 60.0 )
east_pos_interp = np.interp( y_array_e , east_pos.y , east_pos.x_max )
west_pos_interp = np.interp( y_array_w  , west_pos.y , west_pos.x_min )

plt.scatter( dataframe_mesh.centroid_longitude , dataframe_mesh.centroid_latitude , color = 'gray' , s = 0.5 )
for transect_num, group in tes.groupby('transect_num'):
    plt.scatter(group.sort_values( [ 'longitude' , 'latitude' ] ).longitude, group.sort_values( [ 'longitude' , 'latitude' ] ).latitude, label=f'Transect {int(transect_num)}', linestyle='-', s = 1 , color = 'blue' , zorder = 100 )
plt.scatter( east_pos.x , east_pos.y , color = 'green' )
plt.plot( east_pos_interp , y_array_e , color = 'black' ,linestyle = '--' )
plt.plot( west_pos_interp , y_array_w , color = 'black' ,linestyle = '--' )
plt.show( )


pos[ 'y_rounded' ] = pd.cut( pos[ 'y' ] , inter )
pos[ 'x_east' ] = pos.groupby( [ 'y_rounded' ] , observed = True )[ 'x_max' ].transform( 'max' )
pos[ 'x_dist' ] = np.abs( pos[ 'x_east' ] - pos[ 'x_max' ] )

east_pos = pos[ ( pos.x == pos.x_max ) & ( pos.x_dist < 1 ) ].copy( )
y_array = np.arange( unique_coords.y.min( ) , unique_coords.y.max( ) , y_resolution )
east_pos_interp = np.interp( y_array , east_pos.y , east_pos.x_max )

west_pos[ ~ west_pos.transect_num.isin( east_pos.transect_num ) ]

plt.scatter( dataframe_mesh.centroid_longitude , dataframe_mesh.centroid_latitude , color = 'gray' , s = 0.5 )
for transect_num, group in tes.groupby('transect_num'):
    plt.scatter(group.sort_values( [ 'longitude' , 'latitude' ] ).longitude, group.sort_values( [ 'longitude' , 'latitude' ] ).latitude, label=f'Transect {int(transect_num)}', linestyle='-', s = 1 , color = 'black' , zorder = 100 )
plt.scatter( east_pos.x_max , east_pos.y , color = 'red' )
plt.scatter( west_pos.x_min , west_pos.y , color = 'blue' )
plt.plot( east_pos_interp , y_array_e , color = 'black' ,linestyle = '--' )
plt.plot( west_pos_interp , y_array_w , color = 'black' ,linestyle = '--' )
# plt.scatter( west_pos.x_west , west_pos.y , color = 'blue' )
plt.show( )

plt.plot( east_pos.x , east_pos.y )
plt.show()

bb = west_pos.groupby( [ 'transect_num' ] )[ 'x' ].size( )
bb[ bb > 1 ]

y_array = np.arange( unique_coords.y.min( ) , unique_coords.y.max( ) , y_resolution )
east_pos_interp = np.interp( y_array , east_pos.y , east_pos.x )

for transect_num, group in tes.groupby('transect_num'):
    plt.scatter(group.sort_values( [ 'longitude' , 'latitude' ] ).longitude, group.sort_values( [ 'longitude' , 'latitude' ] ).latitude, label=f'Transect {int(transect_num)}', linestyle='-', s = 1 , color = 'blue' , zorder = 100 )
plt.scatter( mesh_grid.centroid_longitude , mesh_grid.centroid_latitude )
plt.plot( east_pos_interp , y_array , color = 'black' ,linestyle = '--' )
plt.show( )
# BELOW --> 49.0 degrees
region_1 = unique_coords[ unique_coords.y <= 51.90 ]
west_x = region_1.groupby( [ 'transect_num' ] )[ 'x' ].min( ).reset_index( )
west_pos = region_1[ region_1.x.isin( west_x.x ) ].reset_index( drop = True )
east_x = region_1.groupby( [ 'transect_num' ] )[ 'x' ].max( ).reset_index( )
east_pos = region_1[ region_1.x.isin( east_x.x ) ].reset_index( drop = True )
max_y_per_x = region_1.groupby(['transect_num'])['y'].max()
mean_x_per_y = pd.DataFrame( region_1.groupby(['transect_num'])['x'].min() )
mean_x_per_y[ 'y' ] = max_y_per_x
mean_x_per_y = mean_x_per_y.reset_index( ).set_index( [ 'x' ] )
sorted_indices = np.argsort(mean_x_per_y.index)
x_values_sorted = mean_x_per_y.index.values[sorted_indices]
y_values_sorted = mean_x_per_y[ 'y' ].values[sorted_indices]
y_resolution = 1.25 / 60
y_array = np.arange( region_1.y.min( ) , region_1.y.max( ) , y_resolution )
west_pos = west_pos.sort_values( by = [ 'y' ] , ascending = True )
west_pos_interp = np.interp( y_array , west_pos.y , west_pos.x )
east_pos_interp = np.interp( y_array , east_pos.y , east_pos.x )

for transect_num, group in tes.groupby('transect_num'):
    plt.scatter(group.sort_values( [ 'longitude' , 'latitude' ] ).longitude, group.sort_values( [ 'longitude' , 'latitude' ] ).latitude, label=f'Transect {int(transect_num)}', linestyle='-', s = 1 , color = 'blue' , zorder = 100 )
plt.scatter( mesh_grid.centroid_longitude , mesh_grid.centroid_latitude )
plt.plot( east_pos_interp , y_array , color = 'black' ,linestyle = '--' )
plt.show( )



west_x = unique_coords.groupby( [ 'transect_num' ] )[ 'x' ].min( ).reset_index( )
west_pos = unique_coords[ unique_coords.x.isin( west_x.x ) ].reset_index( drop = True )
east_x = unique_coords.groupby( [ 'transect_num' ] )[ 'x' ].max( ).reset_index( )
east_pos = unique_coords[ unique_coords.x.isin( east_x.x ) ].reset_index( drop = True )
# Interpolate the upper boundary
max_y_per_x = unique_coords.groupby(['transect_num'])['y'].max()
mean_x_per_y = pd.DataFrame( unique_coords.groupby(['transect_num'])['x'].min() )
mean_x_per_y[ 'y' ] = max_y_per_x
mean_x_per_y = mean_x_per_y.reset_index( ).set_index( [ 'x' ] )
# Sort values by x-coordinate for interpolation
sorted_indices = np.argsort(mean_x_per_y.index)
x_values_sorted = mean_x_per_y.index.values[sorted_indices]
y_values_sorted = mean_x_per_y[ 'y' ].values[sorted_indices]

y_resolution = 1.25 # nmi
y_resolution_m = 1.25 / 60
y_array = np.arange( unique_coords.y.min( ) , unique_coords.y.max( ) , y_resolution_m )

west_pos = west_pos.sort_values( by = [ 'y' ] , ascending = True )
west_pos_interp = np.interp( y_array , west_pos.y , west_pos.x )
east_pos_interp = np.interp( y_array , east_pos.y , east_pos.x )

import matplotlib.pyplot as plt

plt.scatter(west_pos_interp, y_array, color='black', label='Convex Shape', linestyle='--')
plt.scatter(east_pos_interp, y_array, color='red', label='Convex Shape', linestyle='--')
plt.show( )

x = np.arange(-5.01, 5.01, 0.25)
y = np.arange(-5.01, 5.01, 0.25)
xx, yy = np.meshgrid(x, y)
z = np.sin(xx**2+yy**2)
f = interpolate.interp2d(x, y, z, kind='cubic')

# Interpolate upper boundary (assuming linear interpolation)
upper_boundary_interpolated = np.interp(unique_coords['x'], x_values_sorted, y_values_sorted)
north_pos = pd.DataFrame(
    {
        'x': x_values_sorted  ,
        'y': y_values_sorted ,
    } ,
)

plt.plot( north_pos.sort_values( [ "x" , "y" ] ).x , north_pos.sort_values( [ "x" , "y" ] ).y )
plt.show( )

tas = transect_extent( spatial_data , projection , nearest_neighbors = 4 )
transformer = Transformer.from_crs('epsg:32610', 'epsg:4326', always_xy=True)
new_tas = transform(transformer.transform, tas)
mesh_grid = mesh_crop( spatial_data , dataframe_mesh , nearest_neighbors = 4 , grid_buffer = 0 )
lat_col = [ col for col in mesh_grid.columns if 'lat' in col.lower( ) ][ 0 ]
# ---- Longitude
lon_col = [ col for col in mesh_grid.columns if 'long' in col.lower( ) ][ 0 ]  

# Convert to GeoDataFrame
mesh_grid_gdf = gpd.GeoDataFrame(
    mesh_grid ,
    geometry = gpd.points_from_xy( mesh_grid[ lon_col ] ,
                                    mesh_grid[ lat_col ] ) ,
    crs = projection ,       
)
mesh_grid_points = mesh_grid_gdf[ [ 'centroid_longitude' , 'centroid_latitude' ] ].to_numpy( )
mesh_grid_multipoint = MultiPoint(mesh_grid_points)
hole_polygon = mesh_grid_multipoint.intersection( new_tas )
# hole_polygon = new_tas.intersection( mesh_grid_multipoint )
hole_polygon = mesh_grid_multipoint.union( new_tas )

for geom in hole_polygon.geoms:
    xs , ys = geom.exterior.xy 
    plt.plot(xs, ys, color='red', label='Convex Shape', linestyle='--' , zorder = 5 )
plt.show()
x , y = new_tas1.exterior.xy
plt.plot(x, y, color='red', label='Convex Shape', linestyle='--' , zorder = 5 )

x, y = [], []
for point in hole_polygon.geoms:
    x.append(point.x)
    y.append(point.y)

plt.scatter(x, y, color='blue', label='Points')

x , y = hole_polygon.exterior.xy
plt.plot(x, y, color='black', label='Convex Shape', linestyle='--' , zorder = 5 )



x , y = new_tas.exterior.xy
plt.plot(x, y, color='red', label='Convex Shape', linestyle='--' , zorder = 5 )
plt.show()

hole_polygon = new_tas.intersects(mesh_grid_multipoint)




x , y = new_tas.exterior.xy
plt.plot(x, y, color='red', label='Convex Shape', linestyle='--' , zorder = 5 )
plt.show()
from shapely.geometry import Polygon, MultiPoint , Point

mesh_grid_points = dataframe_mesh[ [ 'centroid_longitude' , 'centroid_latitude' ] ].to_numpy( )
mesh_grid_multipoint = MultiPoint(mesh_grid_points)
from geopandas import GeoSeries

hole_polygon = mesh_grid_multipoint.difference( new_tas )
x, y = [], []
for point in hole_polygon.geoms:
    x.append(point.x)
    y.append(point.y)

plt.scatter(x, y, color='blue', label='Points')
plt.show()
hole_polygon = new_tas.intersect(mesh_grid_multipoint)

x , y = hole_polygon.exterior.xy
plt.plot( x , y )
plt.show()

# Create distance buffer that limits extrapolation
# ---- Convert `grid_buffer` (nmi) to m

lat_col = [ col for col in mesh_grid.columns if 'lat' in col.lower( ) ][ 0 ]
# ---- Longitude
lon_col = [ col for col in mesh_grid.columns if 'long' in col.lower( ) ][ 0 ]  

# Convert to GeoDataFrame
mesh_grid_gdf = gpd.GeoDataFrame(
    mesh_grid ,
    geometry = gpd.points_from_xy( mesh_grid[ lon_col ] ,
                                    mesh_grid[ lat_col ] ) ,
    crs = projection ,       
)

tas1 = tas.buffer( 0 )
new_tas1 = transform(transformer.transform, unary_union( tas1 ))
tt = mesh_crop( spatial_data , dataframe_mesh , nearest_neighbors = 4 , grid_buffer = 0 )
for geom in new_tas.geoms:
    xs , ys = geom.exterior.xy 
    plt.plot(xs, ys, color='black', label='Convex Shape', linestyle='--' , zorder = 0 )
x , y = new_tas.exterior.xy
plt.plot(x, y, color='black', label='Convex Shape', linestyle='--' , zorder = 0 )


for geom in new_tas1.geoms:
    xs , ys = geom.exterior.xy 
    plt.plot(xs, ys, color='red', label='Convex Shape', linestyle='--' , zorder = 5 )
x , y = new_tas1.exterior.xy
plt.plot(x, y, color='red', label='Convex Shape', linestyle='--' , zorder = 5 )

plt.scatter( mesh_grid.centroid_longitude , mesh_grid.centroid_latitude , s = 0.5 , color = 'gray' , zorder = 0 )
plt.scatter( tt.centroid_longitude , tt.centroid_latitude , s = 0.5 , color = 'orange' , zorder = 1 )

for transect_num, group in tes[ tes.stratum_inpfc == 6 ].groupby('transect_num'):
    plt.scatter(group.sort_values( [ 'longitude' , 'latitude' ] ).longitude, group.sort_values( [ 'longitude' , 'latitude' ] ).latitude, label=f'Transect {int(transect_num)}', linestyle='-', s = 1 , color = 'blue' , zorder = 100 )
plt.show( )

ttet = new_tas.difference( new_tas )
for geom in new_tas1.geoms:
    xs , ys = geom.exterior.xy 
    plt.plot(xs, ys, color='red', label='Convex Shape', linestyle='--' , zorder = 5 )
plt.show( )

if standardize_coordinates is True :
    spatial_data , d_x , d_y = transform_geometry( spatial_data , 
                                                   reference_grid , 
                                                   kriging_parameters ,
                                                   projection )

    transformed_mesh , _ , _ = transform_geometry( dataframe_mesh , 
                                                   reference_grid , 
                                                   kriging_parameters ,
                                                   projection ,
                                                   d_longitude = d_x , d_latitude = d_y )
    
    print( f"Georeferenced dataset ({dataset}) and mesh coordinates transformed." )

spatial_data = spatial_data
transformed_mesh = transformed_mesh
dataframe_mesh = dataframe_mesh
dataframe_geostrata = dataframe_geostrata 
kriging_parmaeters = kriging_parameters
variogram_parameters = variogram_parameters
variable = variable

### -------------------------
### START KRIGING_INTERPOLATION
### -------------------------
### Discretize latitudinal bins
latitude_bins = np.concatenate( [ [ -90.0 ] , dataframe_geostrata.northlimit_latitude , [ 90.0 ] ] )

### Discretize mesh data into the same strata
dataframe_mesh[ 'stratum_num' ] = pd.cut( dataframe_mesh.centroid_latitude ,
                                          latitude_bins , 
                                          labels = list( dataframe_geostrata.stratum_num ) + [ 1 ] ,
                                          ordered = False )
### -------------------------
### ----> START ORDINARY_KRIGING !!!
### -------------------------
### --------> START LOCAL_SEARCH_INDEX !!!
mesh_grid = transformed_mesh
k_max = kriging_parameters[ 'kmax' ]
### ------------> START GRIDDIFY_LAG_DISTANCES !!!
if isinstance( mesh_grid , pd.DataFrame ) and isinstance( spatial_data , pd.DataFrame ):
    ### 
    x_distance = np.subtract.outer( mesh_grid.x_transformed.values ,
                                    spatial_data.x_transformed.values )
    y_distance = np.subtract.outer( mesh_grid.y_transformed.values ,
                                    spatial_data.y_transformed.values )
elif isinstance( mesh_grid , np.ndarray ) and isinstance( spatial_data , np.ndarray ):
    ###
    x_distance = np.subtract.outer( mesh_grid , mesh_grid )
    y_distance = np.subtract.outer( spatial_data , spatial_data )  

distance_matrix = np.sqrt( x_distance * x_distance + y_distance * y_distance )
### ------------> END GRIDDIFY_LAG_DISTANCES !!!
sorted_distance_matrix = np.argpartition( distance_matrix , k_max , axis = 1 )
distance_matrix = distance_matrix
local_point_grid = sorted_distance_matrix[ : , : k_max ]
### --------> END LOCAL_SEARCH_INDEX !!!
### Initialize kriging grid, weights, and results prior to loop    
kriging_prediction_variance = np.empty( local_point_grid.shape[ 0 ] ) # Prediction variance    
kriging_sample_variance = np.empty( local_point_grid.shape[ 0 ] ) # Sample variance     
kriging_mean = np.empty( local_point_grid.shape[ 0 ] ) # Mean

### Pull target data variable
variable_data = spatial_data[ variable ].values
variable_x = spatial_data.x_transformed.values
variable_y = spatial_data.y_transformed.values

nmax_var = 500
rR = np.linspace( 0.0 , np.sqrt( 2 ) + 0.01 , nmax_var )
dr = np.diff( rR ).mean( )

r = rR
Nugt = variogram_parameters[ 'nugget' ]
Sill = variogram_parameters[ 'sill' ]
C = Sill - Nugt
L = variogram_parameters[ 'correlation_range' ]
p = variogram_parameters[ 'decay_power' ]
b = variogram_parameters[ 'hole_effect_range' ]
rL = r / L
R = 0.021

vgram = C * ( 1.0 - np.exp( - ( r / L ) ** p ) ) + Nugt
transformed_mesh.x_transformed.max( )
mesh_idx = 4
xp = transformed_mesh.x_transformed
yp = transformed_mesh.y_transformed
xn = variable_x
yn = variable_y

Vp = np.zeros( len( yp ) )
Ep = np.zeros( len( yp ) )
Eps = np.zeros( len( yp ) )
k_min = kriging_parameters[ 'kmin' ]
k_max = kriging_parameters[ 'kmax' ]
variogram_parameters[ 'range' ] = 0.006

# distance_matrix , local_point_grid  = local_search_index( transformed_mesh , 
#                                                           spatial_data , 
#                                                           kriging_parameters[ 'kmax' ] )

# local_point_grid[ i ]
# np.sort( distance_matrix[ i ][ local_point_grid[ i ] ] )
# r_sort[ 0 : 10 ]

# test_ind_vals = distance_matrix[ i ][ local_point_grid[ i ] ]
# np.sort( test_ind_vals[ test_ind_vals <= R ] )

# distance_matrix[ test_ind ]
# distance_matrix[ distance_matrix[ local_point_grid[ i ] ] <= R ]

import time

distance_matrix , local_point_grid  = local_search_index( transformed_mesh , 
                                                          spatial_data , 
                                                          kriging_parameters[ 'kmax' ] )


local_matrix = distance_matrix[ local_point_grid ]
windowed_matrix = local_matrix[ local_matrix <= R ]
x_wb = spatial_data[ [ 'transect_num' , 'x_transformed' , 'y_transformed' ] ].groupby( [ 'transect_num' ] )[ [ 'x_transformed' ] ].idxmin( )
xx_wb = spatial_data[ [ 'transect_num' , 'x_transformed' , 'y_transformed' ] ].loc[ x_wb.x_transformed ]
xx_wb.rename( columns = { 'x_transformed': 'x_west' , 'y_transformed': 'y_west' } , inplace = True )
spatial_data = spatial_data.merge( xx_wb )

for i in range( len( yp ) ):
# for i in range( 1000 ):
    st = time.time( )
    xpp = xp[ i ]
    ypp = yp[ i ]
    Dx = ( xpp - xn )
    Dx2 = Dx * Dx
    Dy = ( ypp - yn )
    Dy2 = Dy * Dy
    r = np.sqrt( Dx2 + Dy2 )
    M2_unity = 1
    r_sort = np.sort( r )
    indx_sort = np.argsort( r )

    indx = r_sort[ r_sort <= R ]

    if len( indx ) > k_max:
        indx1 = indx_sort[ 0 : ( k_max ) ]
    elif len( indx ) < k_min:
        indx1 = indx_sort[ 0 : ( k_min ) ]
    else:
        indx1 = indx_sort[ 0 : len( indx ) ]
    nk = len( indx1 )

    variogram( distance_lags = np.array( [ 0.002 , 0.003 , 0.004 , 0.005 , np.nan , np.nan ] ) , variogram_parameters = variogram_parameters , decay_power = 1.5 )
    M20 = variogram( distance_lags = r[ indx1 ] , variogram_parameters = variogram_parameters , decay_power = 1.5 )
    M2 = np.concatenate( [ M20 , [ 1.0 ] ] )

    ### Only need to calculate once
    x1 = xn[ indx1 ]
    y1 = yn[ indx1 ]
    var1 = variable_data[ indx1 ]
    ind0 = var1[ var1 > 1e-5 ]
    nvar_pts = len( ind0 )
    npts = len( variable_data )

    mat_x = np.array( [ x1 ] * nk )
    dx1 = mat_x.T - mat_x
    mat_y = np.array( [ y1 ] * nk )
    dy1 = mat_y.T - mat_y
    rs = np.sqrt( dx1 * dx1 + dy1 * dy1 )
    K0 = variogram( distance_lags = rs , variogram_parameters = variogram_parameters , decay_power = 1.5 )
    K = np.concatenate( [ K0 , np.ones( ( nk , 1 ) ) ] , axis = 1 )
    K = np.concatenate( [ K , np.ones( ( 1 , nk + 1 ) ) ] , axis = 0 )
    np.fill_diagonal( K , 0.0 )

    U , S , V = np.linalg.svd( K )
    Sigma_mask = np.abs( S / S[ 0 ] ) > kriging_parameters[ 'anisotropy' ]
    K_inv = np.matmul(
        np.matmul( V.T[:, Sigma_mask] , np.diag( 1.0 / S[ Sigma_mask ] ) ) ,                    
        U[ : , Sigma_mask ].T
    )
    lamb = np.dot( K_inv , M2 )

    # lamb = compute_kriging_weights( kriging_parameters[ 'anisotropy' ] , M2 , K )
    var1[ r[ indx1 ] > R ] = 0.0
    Vp[ i ] = np.nansum( lamb[ 0 : nk ] * var1 ) * M2_unity
    Vp_taper = np.nansum( lamb[ 0 : nk ] * var1 )
    Vpm = np.nanmean( var1 )
    eps = 2.2204e-16
    ind_non0 = var1[ abs( var1 ) > 2 * eps ]
    ind_non0 = np.linspace( 0 , nk - 1 , nk ).astype( int )
    CVs = np.sqrt( np.nanvar( var1[ ind_non0 ] , ddof = 1 ) ) / np.abs( np.nanmean( var1[ ind_non0 ] ) )
    Ep[ i ] = np.nansum( lamb * M2 )
    Eps[ i ] = np.sqrt( Ep[ i ] * np.nanvar( var1[ ind_non0 ] , ddof = 1 ) ) / np.abs( Vp[ i ] )
    ed = time.time( )
    print( f"{ed - st} seconds elapsed [{i}]" )

def within_radius_distance( array , radius ) :
    output = array[ array < radius ]
    return np.concatenate( [ output , np.full( len( array ) - len( output ) , np.nan ) ] )

def within_radius_index( array , radius ) :
    return np.sort( tt[ tt < R ] )

pad = 9278
from echopop.computation.spatial_old import griddify_lag_distances
distance_matrix = griddify_lag_distances( mesh_grid , spatial_data )
R = self.statistics[ 'kriging' ][ 'model_config' ][ 'search_radius' ]
dd = np.apply_along_axis( within_radius_distance , 1 , distance_matrix , R )

def within_radius_mask( matrix , radius ) :
    matrix_copy = matrix.copy( )
    mask = matrix_copy > radius 
    matrix_copy[ mask ] = np.nan
    return matrix_copy

def collect_nan( matrix ) :
    matrix_copy = matrix.copy( )
    nan_mask = np.isnan( matrix_copy )
    nan_counts = np.sum( ~ nan_mask , axis = 1 )
    return nan_counts

dd = within_radius_mask( distance_matrix , R )
adc = collect_nan( dd )
k_min = self.statistics[ 'kriging' ][ 'model_config' ][ 'kmin' ]
k_max = self.statistics[ 'kriging' ][ 'model_config' ][ 'kmax' ]
x_wb = spatial_data[ spatial_data.latitude < 51 ][ [ 'transect_num' , 'x_transformed' , 'y_transformed' ] ].groupby( [ 'transect_num' ] )[ [ 'x_transformed' ] ].idxmin( )
xx_wb = spatial_data[ [ 'transect_num' , 'x_transformed' , 'y_transformed' ] ].loc[ x_wb.x_transformed ]
transect_extent = xx_wb

def define_western_extent( transect_data: pd.DataFrame , 
                           latitude_threshold: float = 51.0 ):
    """
    Find the western extent of the survey transect data required for extrapolating 
    semivariogram ranges during ordinary kriging

    Parameters
    ----------
    transect_data: pd.DataFrame 
        A DataFrame that includes the georeferenced coordinates of survey transect intervals
    latitude_threshold: float
        A threshold that is applied to the georeferenced coordinates that further constrains
        any extrapolation that occurs during the kriging analysis. 
    """

    # Apply latitude threshold
    transect_thresholded = transect_data[ transect_data.latitude < latitude_threshold ]

    # Reduce the dataframe to only the necessary columns
    transect_thresholded = transect_thresholded[ [ 'transect_num' , 'x_transformed' , 'y_transformed' ] ]

    # Parse the western-most coordinates of each transect
    western_extent_idx = (
        transect_thresholded.groupby( [ 'transect_num' ] )[ [ 'x_transformed' ] ].idxmin( )
    )

    # Use indices from `western_extent_idx` to extract the actual dataframe rows
    transect_western_extent = transect_thresholded.loc[ western_extent_idx.x_transformed ]

    # Return output
    return transect_western_extent

def kriging_grid_search( distance_matrix : np.ndarray ,
                         transect_extent: pd.DataFrame , 
                         k_min = int , 
                         k_max = int , 
                         R = float , ):
    """
    Find the indices of the k-th nearest points (relative to a reference coordinate) required
    for computing the lagged semivariogram   

    Parameters
    ----------
    distance_matrix: np.ndarray
        An array/matrix that includes the distances of each mesh points from every 
        georeferenced along-transect interval
    tranect_extent: pd.DataFrame
        A DataFrame that includes the x- and y-coordinates of the western-most coordinates for 
        each survey transect
    k_min: int 
        The minimum number of nearest neighbors. Mesh points with fewer than `k_min` valid distances
        are supplemented with extrapolated ranges numbering up to `k_min`.
    k_max: int
        The maximum number of nearest neighbors.
    R: float 
        The search radius used to identify between `k_min` and `k_max` nearest neighbors.
    """ 

    # Apply a NaN mask to distance matrix values outside of the search radius R
    valid_distance_matrix = within_radius_mask( distance_matrix , R )

    # Calculate the row-wise number of NaN values within `valid_distance_matrix` to determine whether 
    # extrapolation to fulfill k[min] is required
    invalid_distances = collect_nan( valid_distance_matrix )

    # Find rows where total_counts is less than k_min
    sparse_radius = np.hstack( np.where( invalid_distances < k_min ) )

    # Calculate the closest grid points for all rows
    local_indices = distance_matrix.argsort( axis = 1 )
    local_points = np.take_along_axis( distance_matrix , local_indices , axis = 1 )

    # Initialize out-of-sample weights
    out_of_sample_weights = np.ones( len( invalid_distances ) )

    # Initialize out-of-sample (extrapolated) indices
    out_of_sample_indices = np.full( ( len( invalid_distances ), k_min ) , np.nan )

    # Initialize within-sample (interpolated, not extrapolated) indices
    # ---- This assumes that all of the local point values are valid between 
    # ---- k[min] and k[max]
    within_sample_indices = local_indices[ : , : k_max ].astype( float )

    # Define the within- and out-of-sample indices for mesh points where 
    # fewer than k[min] values are present within the search radius (R)
    if np.size( sparse_radius ) > 0 :

        # Identify the k[min]th closest points -- values up to k[min] beyond R are considered
        # to be extrapolated
        closest_indices = local_indices[ sparse_radius ][ : , : k_min ]

        # Index the mesh grid coordinates needed for bounding the search radius expansion
        # ---- y-coordinates with array transformation to access matrix operations
        mesh_y = np.array( mesh_grid.loc[ sparse_radius ].y_transformed).reshape(-1, 1)
        # ---- x-coordinates
        mesh_x = mesh_grid.loc[ sparse_radius ].x_transformed.values

        # Update local points
        # ---- Fill out NaN values
        local_points[ sparse_radius , k_min : ] = np.nan

        # Calculate the mesh distance matrix from the western boundary of the survey transect coordinates
        # ---- This further finds the index of the closest western boundary y-coordinate
        mesh_distance_matrix = np.abs(mesh_y - np.array( transect_extent.y_transformed ) ).argmin( axis = 1 )

        # Identify the western-most x-coordinate based on `mesh_distance_matrix`
        western_limit = transect_extent.x_transformed.iloc[ np.ravel( mesh_distance_matrix ) ]

        # Construct new threshold that limits the spatial extent of the extrapolation beyond the survey transect region
        western_threshold = western_limit - R 

        # Develop the tapered function for extrapolating values outside the search radius
        # ---- Based on the western-most transect values
        # -------- Evaluate index for closest western-most latitude
        if np.any( mesh_x < western_threshold ) :

            # ---- Create a thresholded mask for lazy-coding downstream in this operation
            western_threshold_mask = mesh_x < western_threshold
            # ---- Index these values
            extrapolation_index = sparse_radius[ western_threshold_mask ]
            # ---- Calculate the out-of-sample kriging weights
            # -------- Out-of-sample mean
            oos_mean = np.apply_along_axis( np.nanmean , 1 , local_points[ extrapolation_index , : k_min ] )
            # -------- Exponentiate 
            oos_exp = np.exp( - oos_mean / R )
            # -------- Update `out_of_sample_weights`
            out_of_sample_weights[ extrapolation_index ] = oos_exp
            # ---- Get outside indices that correspond to this extrapolation
            # -------- Get sparse index associated with extrapolation
            sparse_extrapolation_index = closest_indices[ western_threshold_mask ].astype( float )
            # -------- Apply indices as a mask to the NaN-included distance matrix
            extrapolated_distance = np.take_along_axis( valid_distance_matrix[ sparse_radius ][ western_threshold_mask ] , 
                                                        sparse_extrapolation_index.astype( int ) , axis = 1 )
            # -------- Create NaN mask
            extrapolated_nan_mask = ~ np.isnan( extrapolated_distance )
            # -------- Apply mask to indices
            sparse_extrapolation_index_nan = sparse_extrapolation_index.copy( )
            sparse_extrapolation_index_nan[ extrapolated_nan_mask ] = np.nan
            # -------- Update `out_of_sample_indices` matrix
            out_of_sample_indices[ extrapolation_index ] = np.sort( sparse_extrapolation_index_nan )
            # ---- Get inside indices that apply to these points
            # -------- Create NaN mask for within-sample values
            interpolated_nan_mask = np.isnan( extrapolated_distance )
            # -------- Apply mask to indices
            sparse_interpolation_index_nan = sparse_extrapolation_index.copy( )
            sparse_interpolation_index_nan[ interpolated_nan_mask ] = np.nan
            # -------- Pad NaN to match `within_sample_indices` matrix
            sparse_interpolation_pad = np.pad( sparse_interpolation_index_nan ,
                                               [ ( 0 , 0 ) , ( 0 , k_max - k_min ) ] ,
                                               mode = 'constant' ,
                                               constant_values = np.nan )
            # -------- Updated `within_sample_indices` matrix
            within_sample_indices[ extrapolation_index ] = np.sort( sparse_interpolation_pad )

            print( f"""Extrapolation applied to kriging mesh points ({len( sparse_radius )} of {within_sample_indices.shape[ 0 ]}):
                    * {len( invalid_distances[ invalid_distances == 0 ] )} points had 0 valid range estimates without extrapolation
                    * {len( invalid_distances[ ( invalid_distances != 0 ) & ( invalid_distances< k_min ) ] )} points had at least 1 valid point but fewer than {k_min} valid neighbors""")

    # Return output
    return local_points[ : , : k_max ] , within_sample_indices , out_of_sample_indices , out_of_sample_weights

# Get western extent of all transects
transect_extent = define_western_extent( spatial_data )
range_grid , inside_indices , outside_indices, outside_weights = kriging_grid_search( distance_matrix , transect_extent , kriging_parameters[ 'kmin' ] , k_max , kriging_parameters[ 'search_radius' ] )

# Calculate the M2 lagged semivariogram
local_variograms_M20 = np.apply_along_axis( variogram , 1 , range_grid , variogram_parameters , decay_power = 1.5 )

# Append 1.0 for ordinary kriging assumption
local_variograms_M2 = np.sort( np.hstack( ( local_variograms_M20 , np.ones( ( local_variograms_M20.shape[ 0 ] , 1 ) ) ) ) )

# Construct variable matrix
def spatial_data_matrix( indices , variable_data , length ) :

    # Index the NaN values
    # ---- Interior 
    interior_non_nan = indices[ : length ][ ~ np.isnan( indices[ : length ] ) ].astype( int ) 
    # ---- Exterior
    exterior_non_nan = indices[ length : ][ ~ np.isnan( indices[ length : ] ) ].astype( int ) 
    # ---- Concatenate indices
    full_indices = np.concatenate( [ interior_non_nan , exterior_non_nan ] )

    # Map the indexed variable data for both interpolated and extrapolated indices
    # ---- Pad interior values
    interior_values = np.pad( variable_data[ interior_non_nan ] ,
                              [ ( 0 , length - len( interior_non_nan ) ) ] ,
                              mode = 'constant' ,
                              constant_values = np.nan )
    # ---- Set exterior (or extrapolated) values
    if len( exterior_non_nan ) > 0 :
        exterior_values = variable_data[ exterior_non_nan ]
        # -------- Pad exterior values
        exterior_values = np.pad( np.array( exterior_values ) ,
                                  [ ( 0 , len( indices[ length : ] ) - len( exterior_values ) ) ] ,
                                  mode = 'constant' ,
                                  constant_values = np.nan )
        # -------- Set extrapolated values to 0.0
        exterior_values[ ~ np.isnan( exterior_values ) ] = 0.0
    else :
        exterior_values = indices[ length : ] 
    
    # ---- Concatenate all values
    full_values = np.sort( np.concatenate( [ interior_values , exterior_values ] ) )[ : length ]

    # Return output
    return full_values

# Stack (horizontal) the interpolated and extrapolated indices
stacked_indices = np.hstack( ( inside_indices.tolist( ) , outside_indices.tolist( ) ) )

# Map the spatial data variables based on the horizontally stacked indices
indexed_values = np.apply_along_axis( spatial_data_matrix , 1 , stacked_indices , variable_data , k_max )

def local_kriging_matrix( indices , x_coordinates , y_coordinates , variogram_parameters , length ) :

    # Find the correct non-NaN indices
    non_na_indices = indices[ ~ np.isnan( indices ) ].astype( int )

    # Gather the coordinates belonging to the appropriate indices
    # ---- x-coordinates
    x_local = x_coordinates[ non_na_indices ]
    # ---- y-coordinates
    y_local = y_coordinates[ non_na_indices ]

    # Compute the kriging matrix for the local grid
    indexed_kriging_matrix = kriging_matrix( x_local , y_local , variogram_parameters )

    # Pad the matrix so that it conforms to a uniform size for the entire array
    # ---- Calculate dimension differences (i.e. rows and columns)
    pad_length = length + 1 - len( indexed_kriging_matrix )
    # ---- NaN pad
    padded_indexed_kriging_matrix = np.pad( indexed_kriging_matrix ,
                                            [ ( 0 , pad_length ) , ( 0 , pad_length ) ] ,
                                            mode = 'constant' ,
                                            constant_values = np.nan )

    # Return output
    return padded_indexed_kriging_matrix

# Define coordinates
# ---- x
xx = spatial_data.x_transformed.values
# ---- y
yy = spatial_data.y_transformed.values

from echopop.computation.kriging_methods import kriging_matrix
# Calculate the local kriging matrix for each mesh point
point_kriging_matrix = local_kriging_matrix( stacked_indices , xx , yy , variogram_parameters , k_max )
tt = local_kriging_matrix( stacked_indices[ 2394 ] , xx , yy , variogram_parameters , k_max )
test = np.apply_along_axis( local_kriging_matrix , 1 , stacked_indices[ [ 0 , 2394 ] ] , xx , yy , variogram_parameters , k_max )

from echopop.computation.kriging_methods import kriging_lambda
lam = kriging_lambda( 0.001 , local_variograms_M2[ 0 ] , km )


np.vstack( inside_indices.to_list( ) , la)

def kriging_estimate()

# Apply the lambda weights to calculate the local point estimates
local_weighted_values = lam[ : k_max ] * var1[ inside_indices[ 0 ].astype( int ) ]

# Sum over the weighted point estimates to get the kriged mesh node point estimate
local_value = np.nansum( local_weighted_values ) * outside_weights[ 0 ]

# Redefine local values
var2 = var1[ inside_indices[ 0 ].astype( int ) ]

# Calculate the sample and prediction variances
# ---- Local mean
local_mean = np.nanmean( var1[ inside_indices[ 0 ].astype( int ) ] )
# ---- Define precision
eps = 2.2204e-16
# ---- Find all values that are above the machine precision
index_non_zero = var2[ abs( var2 ) > 2 * eps ]
index_non_zero = np.linspace( 0 , k_max - 1 , k_max ).astype( int )
# ---- Calculate the coefficient of variation (CV)
if local_mean < 2 * eps :
    CV_arr = np.nan
else :
    CV_arr = np.sqrt( np.nanvar( var2[ index_non_zero ] , ddof = 1 ) ) / np.abs( np.nanmean( var2[ index_non_zero ] ) )
# ---- Calculate the kriged prediction variance
prediction_variance = np.nansum( lam * local_variograms_M2[ 0 ] )
# ---- Calculate kriged sample variance
if local_value < 2 * eps :
    sample_variance = np.nan
else :
    sample_variance = np.sqrt( prediction_variance * np.nanvar( var2[ index_non_zero ] , ddof = 1 ) ) / np.abs( local_value )


local_variograms_M2[ 0 ] * lam
var1[ r[ indx1.astype( int ) ] > R ] = 0.0
Vp[ i ] = np.nansum( lamb[ 0 : nk ] * var1 ) * M2_unity
Vp_taper = np.nansum( lamb[ 0 : nk ] * var1 )
Vpm = np.nanmean( var1 )
eps = 2.2204e-16
ind_non0 = var1[ abs( var1 ) > 2 * eps ]
ind_non0 = np.linspace( 0 , nk - 1 , nk ).astype( int )
CVs = np.sqrt( np.nanvar( var1[ ind_non0 ] , ddof = 1 ) ) / np.abs( np.nanmean( var1[ ind_non0 ] ) )
Ep[ i ] = np.nansum( lamb * M2 )
Eps[ i ] = np.sqrt( Ep[ i ] * np.nanvar( var1[ ind_non0 ] , ddof = 1 ) ) / np.abs( Vp[ i ] )
ed = time.time( )


exist_index = ~ np.isnan( interpolated_indices )

np.apply_along_axis( variogram , 1 , range_grid , variogram_parameters , decay_power = 1.5 )


local_ranges[ 0 ].size
len( local_ranges[ 0 ] )



mat_mask = np.isnan( inside_indices )
inside_indices[ ~ mat_mask ].min()
outside_indices.min()
np.argpartition( distance_matrix , inside_indices[ ~ mat_mask ].astype( int ) )
np.take_along_axis( distance_matrix[ ~ mat_mask ] , inside_indices[ ~ mat_mask ].astype( int ) , axis = 1 )

distance_within_indices = distance_matrix[ row , : ][ inside_indices ]
### --------> START VARIOGRAM !!!
    modeled_semivariogram = variogram( distance_lags = distance_matrix[ row , distance_within_indices ] , 
                                        variogram_parameters = variogram_parameters )
### --------> END VARIOGRAM !!!
    inside_indices = inside_indices
    outside_indices = outside_indices
    out_of_sample_weights = out_of_sample_weights
    ### For Ordinary Kriging, we shift the lags by 1 to include a constant term (1.0)
    # TODO: Should we put in statements for Objective mapping and Universal Kriging w/ Linear drift?  # noqa
    lagged_semivariogram = np.concatenate( [ modeled_semivariogram , [ 1.0 ] ] )



# Calculate the absolute difference between each value in mesh_y and all values in transect_extent.y_transformed
diff_matrix = np.abs(mesh_y - transect_extent.y_transformed)
        western_limit = np.abs( mesh_grid.y_transformed - transect_extent.y_transformed ).argmin( axis = 1 )

        if mesh_point.x_transformed < transect_extent.iloc[ western_limit ].x_transformed - R :
            # ---- Collect closest coordinate indices from the total distance matrix up to k[min]
            closest_global_points = distance_matrix.argsort( )[ : k_min ]
            # ---- Calculate the out-of-sample kriging weights
            out_of_sample_weights = np.exp( -np.nanmean( np.sort( distance_matrix )[ closest_global_points ] / R ) )

    dd[ rows_to_process ].argsort( axis = 1 )
    dd[ rows_to_process ].shape
    local_indices[0].shape
    # Initialize out-of-sample weights
    out_of_sample_weights = np.ones_like(adc)

    # Expand search radius if the number of local points are insufficient
    insufficient_rows = rows_to_process    
    local_indices[ insufficient_rows ][ : , : k_min ]
    dd[ insufficient_rows ][ 0 , inside_indices[ 0 , 1 ] ]
    if np.size( insufficient_rows ) > 0:

        # Sort the closest within-sample points for insufficient rows


        too_low = dd.copy( )[ insufficient_rows ]
        np.argsort( too_low , axis = 1 )
        np.sort( too_low[ 0 ] )
        local_indices[ : , too_low , : ]

        local_indices = np.argsort( dd , axis = 1 )
        dd[ 0 , local_indices[ 1 ] ]
        ins = local_indices[ insufficient_rows , : k_min ]
        dd[ local_indices[ : , : k_min ] , : ]

        mask = np.zeros_like(dd, dtype=bool)
        mask[np.arange(len(insufficient_rows))[:, None], local_indices[insufficient_rows, :k_min]] = True
        local_indices[ 0 : 10 , : k_min ]
        selected_indices = [local_indices[row, :k_min] for row in insufficient_rows]
        local_indices[ np.ravel( insufficient_rows ) ]
        ins = np.argpartition( dd[ insufficient_rows ] , k_min , axis = 1 )[ : , : k_min ]
        dd[ ins ]
        inside_indices = local_indices[insufficient_rows, :k_min]
        np.argpartition( local_indices , k_min , axis = 1 )

        np.argpartition( dd[ local_indices ] , k_min , axis = 1 )
        # Calculate the out-of-sample kriging weights for insufficient rows
        closest_global_points = np.argsort(distance_matrix.ravel())[:k_min]
        out_of_sample_weights[insufficient_rows] = np.exp(-np.nanmean(np.sort(distance_matrix.ravel()[closest_global_points]) / R))

        # Pad the output rows with NaN values for insufficient rows
        num_cols = distance_matrix.shape[1]
        for row in insufficient_rows:
            num_values = total_counts[row]
            distance_matrix[row, num_values:] = np.nan

    # Apply the function here, adjust according to your specific requirements
    # Your function goes here

    return distance_matrix

# Example usage:
k_min = 10  # Adjust as needed
R = 0.5  # Adjust as needed
distance_matrix = apply_function_for_small_total_counts(distance_matrix, total_counts_per_row, k_min, R)


dd = within_radius_mask( distance_matrix , R )
adc = collect_nan( dd )
rows_to_process = np.where(adc < k_min)
local_indices = np.argsort(dd, axis=1)
out_of_sample_weights = np.ones_like(distance_matrix)

rows_to_process[ adc[ rows_to_process ] < k_min ]

    # Calculate the closest grid points for all rows
    local_indices = np.argsort(distance_matrix, axis=1)

    # Initialize out-of-sample weights
    out_of_sample_weights = np.ones_like(distance_matrix)

    # Expand search radius if the number of local points are insufficient
    insufficient_rows = rows_to_process[total_counts[rows_to_process] < k_min]
    if len(insufficient_rows) > 0:
        # Sort the closest within-sample points for insufficient rows
        inside_indices = local_indices[insufficient_rows, :k_min]

        # Calculate the out-of-sample kriging weights for insufficient rows
        closest_global_points = np.argsort(distance_matrix.ravel())[:k_min]
        out_of_sample_weights[insufficient_rows] = np.exp(-np.nanmean(np.sort(distance_matrix.ravel()[closest_global_points]) / R))

    # Apply the function here, adjust according to your specific requirements
    # Your function goes here

    #return distance_matrix



###
sorted_distance_matrix = np.argpartition( distance_matrix , k_max , axis = 1 )

distance_matrix , sorted_distance_matrix[ : , : k_max ]

local_matrix = distance_matrix[ local_point_grid ]
windowed_matrix = local_matrix[ local_matrix <= R ]
x_wb = spatial_data[ spatial_data.latitude < 51 ][ [ 'transect_num' , 'x_transformed' , 'y_transformed' ] ].groupby( [ 'transect_num' ] )[ [ 'x_transformed' ] ].idxmin( )
xx_wb = spatial_data[ [ 'transect_num' , 'x_transformed' , 'y_transformed' ] ].loc[ x_wb.x_transformed ]
transect_extent = xx_wb

local_point_grid = local_point_grid[ row , : ]
distance_matrix = distance_matrix[ row , : ]
sort_indx = np.argmin( np.abs( transformed_mesh.y_transformed[ row ] - xx_wb.y_transformed ) )
xx_wb.iloc[79]
np.abs( transformed_mesh.y_transformed[ row ] - xx_wb.y_transformed ).argmin( )
np.sort( distance_matrix )[ indcs[ : k_min ] ]
r_sort = np.sort( distance_matrix )
indx_sort = np.argsort( distance_matrix )
indx1 = indx_sort[ : k_min ]
mesh_point = transformed_mesh.loc[ row , : ]
indcs = distance_matrix.argsort( )[ : k_min ]distance_matrix.sort( )
if transformed_mesh.x_transformed[ row ] < xx_wb.iloc[ sort_indx ].x_transformed - R :
    M2_unity = np.exp( -np.nanmean( np.sort( distance_matrix )[ indcs ] ) / R )


a , b =np.min()

row = 169

np.where( np.isclose( transformed_mesh.y_transformed , 0.468896424973205 ) )
transformed_mesh.y_transformed[ np.isclose( transformed_mesh.y_transformed , 0.468896424973205 ) ]
y_wb = spatial_data[ spatial_data.x_transformed.isin( x_wb.x_transformed.values ) ][ 'y_transformed' ]
vpp = Vp
vpp[ vpp < 0.0 ] = 0.0
vpp[ np.isnan( vpp ) ] = 0.0
dataframe_mesh[ 'kriged_rho_a' ] = vpp
dataframe_mesh[ 'kriged_biomass' ] = (
    dataframe_mesh.fraction_cell_in_polygon
    * kriging_parameters[ 'A0' ]
    * dataframe_mesh.kriged_rho_a
)
dataframe_mesh.kriged_biomass.sum( ) * 1e-6
dataframe_mesh.kriged_biomass.sum( ) / len( dataframe_mesh.kriged_biomass )

### Iterate through the local point grid to evaluate the kriged/interpolated results
for row in range( local_point_grid.shape[ 0 ] ):
    ### Parse the k-closest points based on the defined range-threshold
### --------> START RANGE_INDEX_THRESHOLD !!!    
    local_point_grid_copy = local_point_grid[ row , : ]
    distance_matrix_copy = distance_matrix[ row , : ]
    # R = variogram_parameters[ 'range' ]
    R = 0.021
    k_min = kriging_parameters[ 'kmin' ]

    ### Calculate the within-sample grid indices
    inside_indices = np.where( distance_matrix_copy[ local_point_grid_copy ] <= R )[ 0 ]

    ### Expand search radius if the number of local points are insufficient
    if len( inside_indices ) < k_min:

        ### Sort the closest within-sample points
        inside_indices = np.argsort( distance_matrix_copy[ local_point_grid_copy ] )[ : k_min ]

        ### Look beyond the sample to collect indices beyond the range threshold
        outside_indices = np.where( distance_matrix_copy[ local_point_grid_copy[ inside_indices ] ] > R )[ 0 ]

        ### Calculate the appropraite weights for these values beyond the threshold
        # TODO: should we change this to how Chu does it?
        # tapered function to handle extrapolation
        out_of_sample_weights = np.exp( - np.nanmean( distance_matrix_copy[ local_point_grid_copy[ inside_indices ] ] ) / R )

    else:
        ### If the sample size is sufficient
        outside_indices = [ ]
        out_of_sample_weights = 1.0
        
    inside_indices = inside_indices
    outside_indices = outside_indices
    out_of_sample_weights = out_of_sample_weights
### --------> END RANGE_INDEX_THRESHOLD !!!  
    ### Index the within-sample distance indices
    distance_within_indices = local_point_grid[ row , : ][ inside_indices ]
### --------> START VARIOGRAM !!!
    modeled_semivariogram = variogram( distance_lags = distance_matrix[ row , distance_within_indices ] , 
                                        variogram_parameters = variogram_parameters )
### --------> END VARIOGRAM !!!
    inside_indices = inside_indices
    outside_indices = outside_indices
    out_of_sample_weights = out_of_sample_weights
    ### For Ordinary Kriging, we shift the lags by 1 to include a constant term (1.0)
    # TODO: Should we put in statements for Objective mapping and Universal Kriging w/ Linear drift?  # noqa
    lagged_semivariogram = np.concatenate( [ modeled_semivariogram , [ 1.0 ] ] )
    ### Calculate the covariance/kriging matrix
### --------> START COMPUTE_KRIGIGN_MATRIX !!!
    x_coordinates = variable_x[ distance_within_indices ]    
    y_coordinates = variable_y[ distance_within_indices ]
    ### Calculate local distance matrix of within-range samples
### ------------> START GRIDDIFY_LAG_DISTANCES !!!
    mesh_grid_copy = x_coordinates
    spatial_data_copy = y_coordinates
    if isinstance( mesh_grid_copy , pd.DataFrame ) and isinstance( spatial_data_copy , pd.DataFrame ):
        ### 
        x_distance = np.subtract.outer( mesh_grid_copy.x_transformed.values ,
                                        spatial_data_copy.x_transformed.values )
        y_distance = np.subtract.outer( mesh_grid.y_transformed.values ,
                                        spatial_data_copy.y_transformed.values )
    elif isinstance( mesh_grid_copy , np.ndarray ) and isinstance( spatial_data_copy , np.ndarray ):
        ###
        x_distance = np.subtract.outer( mesh_grid_copy , mesh_grid_copy )
        y_distance = np.subtract.outer( spatial_data_copy , spatial_data_copy )  

    local_distance_matrix = np.sqrt( x_distance * x_distance + y_distance * y_distance )
### ------------> END GRIDDIFY_LAG_DISTANCES !!!
    ### Calculate the covariance/kriging matrix (without the constant term)
    kriging_matrix_initial = variogram( distance_lags = local_distance_matrix , 
                                        variogram_parameters = variogram_parameters )
    ### Expand the covariance/kriging matrix with a constant
    # In Ordinary Kriging, this should be '1'
    # Columns
    kriging_matrix = np.concatenate( [ kriging_matrix_initial , 
                                       np.ones( ( x_coordinates.size , 1 ) ) ] , 
                                    axis = 1 )
    
    # Rows
    kriging_matrix = np.concatenate( [ kriging_matrix , 
                                       np.ones( ( 1 , x_coordinates.size + 1 ) ) ] , 
                                    axis = 0 )
    
    ### Diagonal fill (0.0)
    # TODO: Should we put in statements for Objective mapping and Universal Kriging w/ Linear drift?  # noqa
    # Add column and row of ones for Ordinary Kriging
    np.fill_diagonal( kriging_matrix , 0.0 )
### --------> END COMPUTE_KRIGIGN_MATRIX !!!    
    ### Use singular value decomposition (SVD) to solve for the 
    ### kriging weights (lambda)
### --------> START COMPUTE_KRIGING_WEIGHTS !!! 
    ratio = kriging_parameters[ 'anisotropy' ]
    M2 = lagged_semivariogram
    K = kriging_matrix
    ### Singular value decomposition (SVD)
    # U: left singular vectors (directions of maximum variance)
    # Sigma: singular values (amount of variance captured by each singular vector, U)
    # VH: conjugate transpose of the right singular vectors
    U , Sigma , VH = np.linalg.svd( K , full_matrices = True )
    
    ### Create Sigma mask informed by the ratio-threshold
    # The ratio between all singular values and their respective
    # maximum is used to apply a mask that informs which values
    # are used to tabulate the kriging weights (aka lambda)
    Sigma_mask = np.abs( Sigma / Sigma[ 0 ] ) > ratio

    ### Inverse masked semivariogram (K)
    K_inv = np.matmul(
        np.matmul( VH.T[:, Sigma_mask] , np.diag( 1.0 / Sigma[ Sigma_mask ] ) ) ,                    
        U[ : , Sigma_mask ].T
    )
    
    kriging_weights = np.dot( K_inv , M2 )
### --------> END COMPUTE_KRIGING_WEIGHTS !!! 
    ### Compute Kriging value and variance
    # Index the appropriate data values
    point_values = variable_data[ distance_within_indices ]
    ### Calculate the kriged mean, predication variance, and sample variance
### --------> START COMPUTE_KRIGING_STATISTICS !!!     
    ### Remove any extrapolation
    if len( outside_indices ) > 0:
        point_values[ outside_indices ] = 0.0

    ### Calculate locally weighted kriging mean
    local_mean = np.nansum( kriging_weights[: len( inside_indices ) ] * point_values ) * out_of_sample_weights

    ### Calculate locally weighted kriging prediction variance
    local_prediction_variance = np.nansum( kriging_weights * lagged_semivariogram )

    ### Calculate locally weighted point sample variance
    if abs( local_mean ) < np.finfo( float ).eps:
        local_sample_variance = np.nan
    else:
        local_arithmetic_variance = np.nanvar( point_values , ddof = 1 ) # Non-spatial, arithmetic sample variance
        local_sample_variance = np.sqrt( local_prediction_variance * local_arithmetic_variance ) / abs( local_mean )

    ### Carriage return
    estimate = local_mean
    pred_variance = local_prediction_variance 
    samp_variance = local_sample_variance
### --------> END COMPUTE_KRIGING_STATISTICS !!!  
    ### Fill in the iterable values from the output (tuple) of 
    ### `compute_kriging_statistics()`
    kriging_mean[ row ] = estimate
    kriging_prediction_variance[ row ] = pred_variance
    kriging_sample_variance[ row ] = samp_variance

### Remove nonsense values (NaN, negative)
kriging_mean = np.where( ( kriging_mean < 0 ) | np.isnan( kriging_mean ) , 0.0 , kriging_mean )

### Carriage return
estimate = kriging_mean 
pred_variance = kriging_prediction_variance 
samp_variance = kriging_sample_variance
### -------------------------
### ----> END ORDINARY_KRIGING !!!
### -------------------------
dataframe_mesh[ 'B_a_adult_mean' ] = estimate
dataframe_mesh[ 'B_a_adult_prediction_variance' ] = pred_variance
dataframe_mesh[ 'B_a_adult_sample_variance' ] = samp_variance

### Calculate cell area
dataframe_mesh[ 'cell_area_nmi2' ] = (
    dataframe_mesh.fraction_cell_in_polygon * kriging_parameters[ 'A0' ]
)

### Calculate the kriged biomass estimate
dataframe_mesh[ 'B_adult_kriged' ] = (
    dataframe_mesh.B_a_adult_mean * dataframe_mesh.cell_area_nmi2
)

dataframe_mesh[ 'B_adult_kriged' ].sum( ) * 1e-6

### Calculate the kriged biomass estimate
dataframe_mesh[ 'B_adult_kriged' ] = (
    dataframe_mesh.B_a_adult_mean * dataframe_mesh.cell_area_nmi2
)


biomass_conv_df[ 'biomass_all' ].sum( ) * 1e-6
biomass_conv_df[ 'biomass_adult' ].sum( ) * 1e-6
biomass_conv_df[ biomass_conv_df.sex == 'male' ][ 'biomass_all' ].sum( ) * 1e-6
biomass_conv_df[ biomass_conv_df.sex == 'female' ][ 'biomass_adult' ].sum( ) * 1e-6

test = wgt_proportion[ wgt_proportion.sex == 'all' ].drop( 'sex' , axis = 1 )
test_table = test.pivot_table( columns = [ 'age_bin' , 'stratum_num' ] , index = [ 'length_bin' ] , values = 'proportion' , observed = False )
new_out = biomass_conv_df.merge( test , on = [ 'stratum_num' ] )
new_out[ 'apportioned_biomass_all' ] = (
    new_out.biomass_all * new_out.proportion
)
new_out[ 'apportioned_biomass_adult' ] = (
    new_out.biomass_adult * new_out.proportion
)
new_out[ 'rho_biomass_all' ] = (
    new_out[ 'apportioned_biomass_all' ] 
    / new_out[ 'interval_area' ]
)
new_out[ 'rho_biomass_adult' ] = (
    new_out[ 'apportioned_biomass_adult' ] 
    / new_out[ 'interval_area' ]
)

a = new_out.groupby( [ 'latitude' , 'longitude' , 'transect_num' , 'stratum_num' ] )[ 'apportioned_biomass_adult' ].sum( ).reset_index( )
a[ a.transect_num == 33.0 ]
a[ np.round( a.longitude , 4 ) == -1.242025e+02 ]
full_out = new_out
full_out[ 'age' ] = full_out[ 'age_bin' ].apply( lambda x: x.mid )

test = full_out[ ( full_out.age != 1 ) ]
test['apportioned_biomass_all' ].sum( ) * 1e-6
out = test.groupby( [ 'transect_num' , 'stratum_num' , 'longitude' , 'latitude' ] )[ 'apportioned_biomass_adult' ].sum( ).reset_index( name = 'biomass' )
# out.to_csv( "C:/Users/Brandyn/Documents/echopop_biomass.csv" , index = False )



nasc_interval_df.iloc[ 2358 ]
out.iloc[ 2360 ]
"""For comparison, my code yields 1643.1868 kmt, while Chu's code yields 1643.1898 kmt."""


full_out = new_out.groupby( [ 'length_bin' , 'age' ] , observed = False )[ [ 'apportioned_biomass_adult' ] ].sum( ).reset_index( )
full_out_table = full_out.pivot_table( columns = [ 'age' ] , index = [ 'length_bin' ] , values = 'apportioned_biomass_adult' , observed = False )
full_out_table.sum( )
full_out_table.iloc[ : , 1 : 23 ].sum( ).sum( ) * 1e-6





mask = nasc_fraction_total_df[ np.round( nasc_fraction_total_df.NASC_no_age1 ) == 774 ]
mask[ 'rho_number_all' ]
mask[ 'abundance_all' ]
specimen_len_age_proportion[ ( specimen_len_age_proportion.sex == 'all' ) & ( specimen_len_age_proportion.age_bin == pd.Interval( left = 0.5 , right = 1.5 ) ) ]
test = specimen_len_age_proportion[ ( specimen_len_age_proportion.sex == 'all' ) & ( specimen_len_age_proportion.stratum_num == 5 ) ].pivot_table( columns = [ 'stratum_num' , 'age_bin' ] , index = 'length_bin' , values = 'number_proportion_bin' )
test
age1_proportions = ( 1.0 - specimen_len_age_proportion[ ( specimen_len_age_proportion.sex == 'all' ) & ( specimen_len_age_proportion.age_bin != pd.Interval( left = 0.5 , right = 1.5 ) ) ].groupby( [ 'stratum_num' ] )[ 'number_proportion_bin' ].sum( ) ).reset_index( name = 'nonadult_proportion' )
sig_b_coef = 10**(-68.0/10)

new_number_proportions = age1_proportions.merge( specimen_len_age_proportion[ ( specimen_len_age_proportion.sex == 'all' ) ] , on = [ 'stratum_num' ] )
length_bins = self.biology[ 'distributions' ][ 'length' ][ 'length_bins_arr' ]
age_bins = self.biology[ 'distributions' ][ 'age' ][ 'age_bins_arr' ]

full_intervals = pd.DataFrame( { 'length': length_bins , 'length_bin': pd.cut( length_bins , length_intervals ) } ).merge( pd.DataFrame( { 'age': age_bins , 'age_bin': pd.cut( age_bins , age_intervals ) } ) , how = 'cross' )
full_intervals[ 'sigma_bs_scaled' ] = sig_b_coef * full_intervals.length ** 2
full_intervals = pd.concat( [ full_intervals ] * len( np.unique( new_number_proportions.stratum_num ) ) )
full_intervals[ 'stratum_num' ] = np.repeat( np.unique( new_number_proportions.stratum_num ) , len( length_bins ) * len( age_bins ) )

full_intervals_tbl = full_intervals.pivot_table( columns = [ 'stratum_num' , 'age_bin' ] ,
                                                 index = [ 'length_bin' ] ,
                                                 values = 'sigma_bs_scaled' ,
                                                 observed = False )


rescaled_number_densities = pd.DataFrame(
    {
        'length': length_bins ,
        'length_bin': pd.cut( length_bins , length_intervals ) ,
        'sigma_bs_coefficient_scaled': sig_b_coef * length_bins ** 2 ,
    } ,
)

rescaled_sigma_bs = pd.concat( [ rescaled_number_densities ] * len( np.unique( new_number_proportions.stratum_num ) ) )
rescaled_sigma_bs[ 'stratum_num' ] = np.repeat( np.unique( new_number_proportions.stratum_num ) , len( length_bins ) )
rescaled_sigma_bs_tbl = rescaled_sigma_bs.pivot_table( index = [ 'length_bin' ] , columns = [ 'stratum_num' ] , values = 'sigma_bs_coefficient_scaled' , observed = False )

new_number_proportions = new_number_proportions.merge( rescaled_number_densities , on = [ 'length_bin' ] )
# Convert to table
specimen_len_age_proportion_table = new_number_proportions.pivot_table( columns = [ 'stratum_num' , 'age_bin' ] ,
                                                                        index = [ 'length_bin' ] ,
                                                                        values = 'number_proportion_bin' ,
                                                                        aggfunc = 'sum' ,
                                                                        observed = False )

full_intervals_tbl.dot()
specimen_len_age_proportion_table.values
specimen_len_age_proportion_table.values.dot( full_intervals_tbl )
specimen_len_age_proportion_table.dot( rescaled_sigma_bs_tbl )
rescaled_sigma_bs_tbl.dot( specimen_len_age_proportion_table )

all_ages = ( full_intervals_tbl * specimen_len_age_proportion_table ).sum( ).reset_index( name = 'sigma_bs_average' )
all_sigma_bs_average = all_ages.groupby( [ 'stratum_num' ] )[ 'sigma_bs_average' ].sum( )

nonadult_sigma_bs_average = all_ages[ all_ages.age_bin == pd.Interval( left = 0.5 , right = 1.5 ) ].set_index( 'stratum_num' )
nonadult_sigma_bs_average[ 'sigma_bs' ] = all_sigma_bs_average
nonadult_sigma_bs_average[ 'nonadult_nasc_proportion' ] = nonadult_sigma_bs_average[ 'sigma_bs_average' ] / nonadult_sigma_bs_average[ 'sigma_bs' ]
nonadult_sigma_bs_average = nonadult_sigma_bs_average.reset_index( )




adult_sigma_bs_average == all_sigma_bs_average
sigma_bs_average = ( ( rescaled_sigma_bs_tbl * specimen_len_age_proportion_table ).sum( ) )

rescaled_sigma_bs_tbl * specimen_len_age_proportion_table / sigma_bs_average

rescaled_sigma_bs_tbl.dot( specimen_len_age_proportion_table )
rescaled_sigma_bs_tbl.shape
specimen_len_age_proportion_table.shape
specimen_len_age_proportion_table.index
specimen_len_age_proportion_table.columns
rescaled_sigma_bs_tbl.index == specimen_len_age_proportion_table.index
rescaled_sigma_bs_tbl.columns == specimen_len_age_proportion_table.columns
scaled_sigma_bs = sig_b_coef * converted_lengths ** 2

wgt_proportion[ wgt_proportion.sex == 'all' ].groupby( [ 'stratum_num' ] )[ 'proportion' ].sum( )
1 - wgt_proportion[ ( wgt_proportion.sex == 'all' ) & ( wgt_proportion.age_bin != pd.Interval( left = 0.5 , right = 1.5 ) ) ].groupby( [ 'stratum_num' ] )[ 'proportion' ].sum( )

nasc_fraction_total_df.columns
nasc_fraction_total_df[ ( nasc_fraction_total_df.interval > 0.20 ) & ( nasc_fraction_total_df.interval < 0.25 ) & ( nasc_fraction_total_df.NASC_no_age1 > ) ]
nasc_fraction_total_df.loc[ np.round( nasc_fraction_total_df.interval , 4 ) == 15.0163 ][ 'NASC_no_age1' ]
nasc_adult_number_proportions.groupby( [ 'stratum_num' ] )[ 'count_age_proportion_all' ].sum( )
nasc_fraction_total_df.loc[ 2366 ]
### Apportion
full_props = sexed_number_proportions[ sexed_number_proportions.sex != 'all' ]# .merge( specimen_len_age_proportion_all.drop( 'sex' , axis = 1 ) , on = [ 'stratum_num'  ] , how = 'outer' )
full_props.groupby( [ 'stratum_num' ] )[ 'number_proportions' ].sum( )
nasc_fraction_total_mrg_df = nasc_fraction_total_df.merge( full_props , on = [ 'stratum_num' ] )
nasc_fraction_total_mrg_df[ 'sexed_rho_number_all' ] = (
    nasc_fraction_total_mrg_df.rho_number_all * nasc_fraction_total_mrg_df.number_proportions
)
nasc_fraction_total_mrg_df[ 'sexed_rho_number_adult' ] = (
    nasc_fraction_total_mrg_df.rho_number_adult * nasc_fraction_total_mrg_df.number_proportions
)


### Convert to biomass
nasc_areal_density = nasc_fraction_total_mrg_df.merge( stratum_weight_df , on = [ 'stratum_num' , 'sex' ] )
nasc_areal_density[ 'sexed_rho_biomass_all' ] = (
    nasc_areal_density.sexed_rho_number_all
    * nasc_areal_density.average_weight
)
nasc_areal_density[ 'sexed_rho_biomass_adult' ] = (
    nasc_areal_density.sexed_rho_number_adult
    * nasc_areal_density.average_weight
)

nasc_areal_density[ 'biomass_all' ] = (
    nasc_areal_density.interval_area
    * nasc_areal_density.sexed_rho_number_all
)
nasc_areal_density[ 'biomass_adult' ] = (
    nasc_areal_density.interval_area
    * nasc_areal_density.sexed_rho_number_adult
)

test = specimen_sex_proportion
test.groupby( [ 'stratum_num' ] )[ 'number_proportion' ].sum( )
nasc_areal_density[ nasc_areal_density.sex == 'male' ][ 'biomass_adult' ].sum( ) * 1e-6
nasc_areal_density[ nasc_areal_density.sex == 'female' ].biomass_all.max( ) * 1e-6
nasc_areal_density[ nasc_areal_density.sex == 'female' ].biomass_adult.sum( ) * 1e-6


biomass_df[ biomass_df.sex == 'all' ][ 'B' ].sum( ) * 1e-6

sexed_number_proportions
test =  index_transect_age_sex_proportions( copy.deepcopy( self.acoustics ) ,
                                                                     copy.deepcopy( self.biology ) , 
                                                                     self.spatial[ 'strata_df' ].copy() )

test[ ( test.stratum_num == 2 ) ][ 'proportion_female' ]

acoustics_dict[ 'nasc' ][ 'nasc_df' ]
nasc_interval_df[ nasc_interval_df.haul_num == 0 ]
nasc_fraction_total_df

nasc_number_proportions = (
    biology_dict[ 'weight' ][ 'proportions' ][ 'age_proportions_df' ]
)
nasc_weight_proportions = (
        biology_dict[ 'weight' ][ 'proportions' ][ 'age_weight_proportions_df' ]
)



sex_indexed_weight_proportions = copy.deepcopy( self.biology )[ 'weight' ][ 'proportions' ][ 'sex_age_weight_proportions_df' ]

tt = wgt_proportion.pivot_table( columns = [ 'age_bin' ] ,
                                 index = [ 'stratum_num' , 'sex' ] ,
                                 values = 'proportion' ,
                                 aggfunc = 'sum' )

tt.pivot_table( index = [ 'sex' ] , aggfunc = 'sum')

tt = full_wgt_sex.pivot_table( columns = [ 'stratum_num' , 'age_bin' ] ,
                               index = [ 'sex' , 'length_bin' ] ,
                               values = 'weight' )
tt.loc[ ( 'all' ) ]

tt = length_sex_age_weight_proportions[ length_sex_age_weight_proportions.stratum_num == 2 ].pivot_table( columns = 'age' ,
                                                                                                     index = [ 'sex' , 'length_bin' ] ,
                                                                                                     values = 'weight_length_sex_proportion_all' )
tt.loc[ 'male' ]
fitted_weight_table * length_props_table
fitted_weight_table.index
length_props_table.loc[ 1 , 1 ].index
fitted_weight_table * length_props_table.loc[ 1 , 1 ]
( ( length_props_table * station_props_table ).pivot_table( index = 'sex' , aggfunc = 'sum' ).reset_index( ) )

fitted_weight = fitted_weight[ fitted_weight.sex == 'all' ]
specimen_full_proportion_table * station_props_table.loc[ 2 , ]
length_props_table.loc[ 2 , : ]
( length_props_table * station_props_table ).sum( )
fitted_weight_table.dot( length_props_table * station_props_table + station_props_table * length_props_table )

fitted_weight[ 'weight_modeled' ].values.dot( station_props_table.loc[ 1 , : ] * length_props_table.loc[ 1 , ] + station_props_table.loc[ 2 , : ] * length_props_table.loc[ 2 , ] )

length_proportion_all.merge( fitted_weight )

summed_length_bin_proportion = length_len_proportion.groupby( [ 'stratum_num' , 'sex' ] )[ 'number_proportion_bin' ].sum( ).reset_index( )
length_proportion_full = length_proportion_all.merge( summed_length_bin_proportion , on = [ 'stratum_num' , 'sex' ] , how = 'outer' ).dropna( subset = [ 'number_proportion' , 'number_proportion_bin' ] )
length_proportion_full[ 'number_proportion_overall' ] = (
    length_proportion_full.number_proportion * length_proportion_full.number_proportion_bin
)

test = length_proportion_full

test[ 'fac1' ] = test.number_proportion_overall / ( test.number_proportion_overall + test.specimen_number_proportion_overall )

specimen_proportion_all

length_proportion_full[ [ 'stratum_num' , 'sex' , 'specime']]
test = length_proportion_full.merge( specimen_proportion_all , on = [ 'stratum_num' , 'sex' ] )
test[ 'fac1' ] = test.number_proportion_overall / ( test.number_proportion_overall + test.number_specimen_proportion )
test[ 'fac2' ] = test.number_specimen_proportion / ( test.fac1 + test.number_specimen_proportion )

[ 'number_proportion_overall' ] = (
    length_proportion_all.number_proportion_bin * length_proportion_all.number_proportion
)
####
fac1_M = length_proportion_all[ length_proportion_all.sex == 'male' ] / ( length_proportion_all[ length_proportion_all.sex == 'male' ] + )

    fac1_M=fac1_M/(fac1_M+fac2_M);
    fac1_F=fac1_F/(fac1_F+fac2_F);
    fac1_ALL=fac1_ALL/(fac1_ALL+fac2_ALL);
    fac2_M=fac2_M/(fac1_M+fac2_M);
    fac2_F=fac2_F/(fac1_F+fac2_F);
    fac2_ALL=fac2_ALL/(fac1_ALL+fac2_ALL);

length_proportion_all.merge( length_len_proportion , on = [ 'stratum_num' ] , how = 'outer' ).dropna( subset = [ 'number_proportion_bin' ] )
fac1_MFN = fac1_ALL.merge( length_len_proportion , on = [ 'stratum_num' , 'sex' ] , how = 'outer' ).dropna( subset = [ 'count' ] )
fac1_MFN
length_len_proportion[ ( length_len_proportion.sex == 'male' ) & ( length_len_proportion.stratum_num == 2 ) ].pivot_table( index = [ 'length_bin' ] , columns = [ ] , values = 'number_proportion_bin' )
specimen_sex_proportion.groupby( [ 'stratum_num' , 'sex' ] )[ 'count' ].sum( ).reset_index( )
specimen_sex_proportion
specimen_sex_proportion_agg
grand_total

specimen_sex_proportion[ ( specimen_sex_proportion.sex == 'male' ) & ( specimen_sex_proportion.stratum_num == 2.0 ) ].pivot_table( index = [ 'length_bin' ] , columns = [ 'age_bin' ] , values = 'count' )
specimen_sex_proportion[ ( specimen_sex_proportion.sex == 'female' ) & ( specimen_sex_proportion.stratum_num == 2.0 ) ].pivot_table( index = [ 'length_bin' ] , columns = [ 'age_bin' ] , values = 'count' )


specimen_sex_len_age_proportion = (
    specimen_len_age_proportion.groupby( [ 'stratum_num' , 'sex' ] )[ 'number_proportion_bin' ].sum( )
).reset_index( name = 'number_proportion' )
# Length
length_sex_len_proportion = (
    length_len_proportion.groupby( [ 'stratum_num' , 'sex' ] )[ 'number_proportion_bin' ].sum( )
).reset_index( name = 'number_proportion' )

### Proportion stations
# Station 2
fac2_ALL = specimen_sex_len_age_proportion[ specimen_sex_len_age_proportion.sex != 'unsexed' ]
# Station 1
fac1_ALL = fac2_ALL.assign( number_proportion = lambda x: 1.0 - x.number_proportion )
fac1_MFN = fac1_ALL.merge( length_len_proportion , on = [ 'stratum_num' , 'sex' ] , how = 'outer' ).dropna( subset = [ 'count' ] )
fac1_MFN

length_all_total = length_len[ ( length_len.sex == 'all' ) ].groupby( [ 'stratum_num' ] )[ 'count' ].sum( ).reset_index( name = 'length_total' )
length_len_proportion = length_len_proportion.merge( length_all_total[ [ 'stratum_num' , 'true_total' ] ] , on = [ 'stratum_num' ] , how = 'outer' ).dropna( subset = [ 'count' ] )
length_len_proportion[ 'number_proportion_bin' ] = np.where(
    length_len_proportion[ 'sex' ] == 'all' ,
    length_len_proportion[ 'count' ] / length_len_proportion[ 'length_total' ] ,
    length_len_proportion[ 'count' ] / length_len_proportion[ 'sexed_total' ]
)
### Totals
# ---- Unfiltered
specimen_len_age.groupby( [ 'stratum_num' , 'sex' , 'length_bin' ] )

specimen_len_age.assign( station = 2 ).meld( length_len.assign( station = 1 ) )

### Total sexed length-age number proportions per stratum (number key)
### Binned counts 
# Specimen
specimen_len_age = spec.bin_variable( length_distribution , 'length' ).bin_variable( age_distribution , 'age' ).count_variable( contrasts = [ 'stratum_num' , 'sex' , 'length_bin' , 'age_bin' ] , variable = 'length' , fun = 'size' )
# Length 
length_len = leng.bin_variable( length_distribution , 'length' ).count_variable( contrasts = [ 'stratum_num' , 'sex' , 'length_bin'] , variable = 'length_count' , fun = 'sum' )

specimen_total = spec.groupby( [ 'stratum_num' ] )[ 'length' ].size( ).reset_index( name = 'spec_ct' )

### Total sexed length number per stratum
length_total = leng.groupby( [ 'stratum_num' ] )[ 'length_count' ].sum( ).reset_index( name = 'leng_ct' )

### Grand total
grand_total = specimen_total.merge( length_total , how = 'outer' ).fillna( 0.0 ).assign( sexed_total = lambda x: x.spec_ct + x.leng_ct )

### True total
spec_true = self.biology[ 'specimen_df' ][ ( self.biology[ 'specimen_df' ].species_id == species_id ) ]
specimen_true_total = spec_true.groupby( [ 'stratum_num' ] )[ 'length' ].size( ).reset_index( name = 'spec_ct' )
true_total = specimen_true_total.merge( length_total , how = 'outer' ).fillna( 0.0 ).assign( true_total = lambda x: x.spec_ct + x.leng_ct )

### Number key
# Specimen
specimen_len_age = spec.bin_variable( length_distribution , 'length' ).bin_variable( age_distribution , 'age' ).count_variable( contrasts = [ 'stratum_num' , 'sex' , 'length_bin' , 'age_bin' ] , variable = 'length' , fun = 'size' )
# Length
length_len = leng.bin_variable( length_distribution , 'length' ).count_variable( contrasts = [ 'stratum_num' , 'sex' , 'length_bin'] , variable = 'length_count' , fun = 'sum' )
length_len = pd.concat( [ length_len ,
                          length_len.groupby( [ 'stratum_num' , 'length_bin' ] )[ 'count' ].sum( ).reset_index( ).assign( sex = 'all' ) ] )
### Proportion key
# Specimen
specimen_len_age_proportion = specimen_len_age.merge( grand_total[ [ 'stratum_num' , 'sexed_total' ] ] , on = [ 'stratum_num' ] , how = 'outer' ).dropna( subset = [ 'count' ] )
specimen_len_age_proportion[ 'number_proportion_bin' ] = (
    specimen_len_age_proportion[ 'count' ] / specimen_len_age_proportion[ 'sexed_total' ]
)

# Length
length_len_proportion = length_len.merge( grand_total[ [ 'stratum_num' , 'sexed_total' ] ] , on = [ 'stratum_num' ] , how = 'outer' ).dropna( subset = [ 'count' ] )
length_all_total = length_len[ ( length_len.sex == 'all' ) ].groupby( [ 'stratum_num' ] )[ 'count' ].sum( ).reset_index( name = 'length_total' )
length_len_proportion = length_len_proportion.merge( length_all_total[ [ 'stratum_num' , 'true_total' ] ] , on = [ 'stratum_num' ] , how = 'outer' ).dropna( subset = [ 'count' ] )
length_len_proportion[ 'number_proportion_bin' ] = np.where(
    length_len_proportion[ 'sex' ] == 'all' ,
    length_len_proportion[ 'count' ] / length_len_proportion[ 'length_total' ] ,
    length_len_proportion[ 'count' ] / length_len_proportion[ 'sexed_total' ]
)

length_len_proportion[ ( length_len_proportion.sex == 'all' ) & ( length_len_proportion.stratum_num == 2 ) ]

### Sex number proportions
# Specimen
specimen_sex_len_age_proportion = (
    specimen_len_age_proportion.groupby( [ 'stratum_num' , 'sex' ] )[ 'number_proportion_bin' ].sum( )
).reset_index( name = 'number_proportion' )
# Length
length_sex_len_proportion = (
    length_len_proportion.groupby( [ 'stratum_num' , 'sex' ] )[ 'number_proportion_bin' ].sum( )
).reset_index( name = 'number_proportion' )

### Proportion stations
# Station 2
fac2_ALL = specimen_sex_len_age_proportion[ specimen_sex_len_age_proportion.sex != 'unsexed' ]
# Station 1
fac1_ALL = fac2_ALL.assign( number_proportion = lambda x: 1.0 - x.number_proportion )
fac1_MFN = fac1_ALL.merge( length_len_proportion , on = [ 'stratum_num' , 'sex' ] , how = 'outer' ).dropna( subset = [ 'count' ] )
fac1_MFN

specimen_sex_len_age_proportion.reset_index( ).loc[ lambda x: x.sex == 'male' ].assign( number_proportion = lambda x: np.round( x.number_proportion , 4 ) )
length_sex_len_proportion.reset_index( ).loc[ lambda x: x.sex == 'female' ].assign( number_proportion = lambda x: np.round( x.number_proportion , 4 ) )

length_sex_len_proportion.loc[ ( 0 , 1 ): "male" ]

# Grand total
a = spec.groupby( [ 'sex' , 'stratum_num' ] )[ 'weight' ].size( ).reset_index( name = 'spec_ct' )
a_tot = a.groupby( [ 'stratum_num' ] )[ 'spec_ct' ].sum( ).reset_index( name = 'spec_tot' )
b = leng.groupby( [ 'sex' , 'stratum_num' ] )[ 'length_count' ].sum( ).reset_index( name = 'len_ct' )
b_tot = b.groupby( [ 'stratum_num' ] )[ 'len_ct' ].sum( ).reset_index( name = 'len_tot' )
overall_total = a_tot.merge( b_tot , on = [ 'stratum_num' ] , how = 'outer' ).fillna( 0.0 ).assign( total = lambda x: x.len_tot + x.spec_tot )
# Sexed total
aa = spec.groupby( [ 'sex' , 'stratum_num' ] )[ 'length' ].size( ).reset_index( name = 'spec_ct' )
bb = leng.groupby( [ 'sex' , 'stratum_num' ] )[ 'length_count' ].sum( ).reset_index( name = 'len_ct' )
sexed_total = aa.merge( bb , on = [ 'stratum_num' , 'sex' ] , how = 'outer' ).fillna( 0.0 ).assign( total_sex = lambda x: x.len_ct + x.spec_ct )


just_sexed_overall = sexed_total[ sexed_total.sex != 'unsexed' ].groupby( [ 'stratum_num' ] )[ 'total_sex' ].sum( ).reset_index( name = 'sexed_total' )
dd = aa.merge( overall_total , on = 'stratum_num', how = 'outer' ).fillna( 0.0 ).assign( prop = lambda x: x.spec_ct / x.total )
dd[ dd.sex == 'male' ]
just_sexed_overall.drop( [ 'spec_tot' , 'len_tot' ] , axis = 1 ).merge( sexed_total , on = 'stratum_num' ).assign( prop = lambda x: x.total_sex / x.sexed_total )

c = a.merge( b , on = [ 'sex' , 'stratum_num' ] )
total = c.groupby( [ 'stratum_num' ] )[ [ 'len_ct' , 'spec_ct' ] ].apply( lambda x: x.len_ct.sum( ) + x.spec_ct.sum( ) ).reset_index( name = 'total' )

stratum_filter = [ 0 ] 
age_filter = [ 0 , 1 ]
region = [ 'US' , 'CAN' ]
species_id = 22500

dataframe_list = [ self.biology[ 'specimen_df' ] ,
                   self.biology[ 'length_df' ] ,
                   self.biology[ 'catch_df' ] ]

if isinstance( dataframe_list , pd.DataFrame ):
    dataframe_list = [ dataframe_list ]

# Stratum removal
aa = unaged_weights_df[ ( unaged_weights_df.stratum_num == 5 ) ]
aa.haul_weight[ 11 : 20 ]
unaged_weights_df[ ( unaged_weights_df.stratum_num == 5 ) & np.logical_not( unaged_weights_df.haul_num.isin( [ 204 ] ) ) ][ 'haul_weight' ].sum( )
13045.579999999998 - 147.7600
filtered_dataframes = tuple( df[ df.species_id == species_id ] for df in dataframe_list )

self.acoustics['sigma_bs'] = { }
self.biology['weight'] = { }
self.biology['population'] = { }
self.statistics['length_weight'] = { }

survey.transect_analysis( )
survey.stratified_summary( )
survey.standardize_coordinates( )
survey.krige( )

data_regrid = self.biology[ 'population' ][ 'areal_density' ][ 'biomass_density_df' ]
overall = data_regrid[ [ 'transect_num' , 'geometry' ] ]
unique_coords = overall.drop_duplicates( subset = [ 'geometry' ] )
# Convert to UTM (m)
mesh_grid = self.statistics[ 'kriging' ][ 'mesh_df' ]
mesh_grid = mesh_grid[ [ 'geometry' ] ]

import math

def convert_wgs_to_utm(lon: float, lat: float):
    """Based on lat and lng, return best utm epsg-code"""
    utm_band = str((math.floor((lon + 180) / 6 ) % 60) + 1)
    if len(utm_band) == 1:
        utm_band = '0'+utm_band
    if lat >= 0:
        epsg_code = '326' + utm_band
        return epsg_code
    epsg_code = '327' + utm_band
    return epsg_code

utm_code = convert_wgs_to_utm( np.median( unique_coords.geometry.x ) , np.median( unique_coords.geometry.y.min( ) ) )
unique_coords_utm = unique_coords.to_crs( f"epsg:{utm_code}" )
unique_coords_utm[ 'x' ] = unique_coords_utm.geometry.x
unique_coords_utm[ 'y' ] = unique_coords_utm.geometry.y
mesh_grid_utm = mesh_grid.to_crs( f"epsg:{utm_code}" )
buffer = 1 # nmi
coord_centroid = unique_coords_utm.dissolve( by = 'transect_num' , aggfunc = 'mean' )
coord_centroid[ 'geometry' ] = coord_centroid.geometry.centroid
coord_centroid = coord_centroid.reset_index( )
coord_min = unique_coords_utm.dissolve( by = 'transect_num' , aggfunc = 'min' )
coord_min[ 'geometry' ] = coord_min.geometry.values
coord_min = coord_min.reset_index( )
# Max
coord_max = unique_coords_utm.dissolve( by = 'transect_num' , aggfunc = 'max' )
coord_max[ 'geometry' ] = coord_max.geometry.centroid
coord_max = coord_max.reset_index( )

n_close = 2

n_close = 3
polygons = []
for transect in np.unique( unique_coords_utm[ 'transect_num' ] ):
    # Extract coordinates and find westernmost and easternmost points
    # pos_min = coord_min[ coord_min.transect_num == transect ][ 'geometry' ]
    pos_centroid = coord_centroid[ coord_centroid.transect_num == transect ][ 'geometry' ]
    # pos_max = coord_max[ coord_max.transect_num == transect ][ 'geometry' ]

    # pos_min_out = coord_min[ coord_min.transect_num != transect ].copy( )
    pos_centroid_out = coord_centroid[ coord_centroid.transect_num != transect ].copy( )
    # pos_max_out = coord_max[ coord_max.transect_num != transect ].copy( )

    # pos_min_out[ 'distance_min' ] = pos_min_out.geometry.apply( lambda g: pos_min.distance( g ) )
    # pos_min_out[ 'distance_centroid' ] = pos_min_out.geometry.apply( lambda g: pos_centroid.distance( g ) )
    # pos_min_out[ 'distance_max' ] = pos_min_out.geometry.apply( lambda g: pos_max.distance( g ) )

    # pos_centroid_out[ 'distance_min' ] = pos_centroid_out.geometry.apply( lambda g: pos_min.distance( g ) )
    pos_centroid_out[ 'distance_centroid' ] = pos_centroid_out.geometry.apply( lambda g: pos_centroid.distance( g ) )
    # pos_centroid_out[ 'distance_max' ] = pos_centroid_out.geometry.apply( lambda g: pos_max.distance( g ) )

    # pos_max_out[ 'distance_min' ] = pos_max_out.geometry.apply( lambda g: pos_min.distance( g ) )
    # pos_max_out[ 'distance_centroid' ] = pos_max_out.geometry.apply( lambda g: pos_centroid.distance( g ) )
    # pos_max_out[ 'distance_max' ] = pos_max_out.geometry.apply( lambda g: pos_max.distance( g ) )

    proximity_coords = [ #pos_min_out[ pos_min_out.distance_min.isin( pos_min_out.distance_min.nsmallest( n_close ) ) ] ,
                        #pos_min_out[ pos_min_out.distance_centroid.isin( pos_min_out.distance_centroid.nsmallest( n_close ) ) ] ,
                        #pos_min_out[ pos_min_out.distance_max.isin( pos_min_out.distance_max.nsmallest( n_close ) ) ] ,
                        #pos_centroid_out[ pos_centroid_out.distance_min.isin( pos_centroid_out.distance_min.nsmallest( n_close ) ) ] ,
                        pos_centroid_out[ pos_centroid_out.distance_centroid.isin( pos_centroid_out.distance_centroid.nsmallest( n_close ) ) ] ]
                        # pos_centroid_out[ pos_centroid_out.distance_max.isin( pos_centroid_out.distance_max.nsmallest( n_close ) ) ] ,
                        # pos_max_out[ pos_max_out.distance_min.isin( pos_max_out.distance_min.nsmallest( n_close ) ) ] ,
                        # pos_max_out[ pos_max_out.distance_centroid.isin( pos_max_out.distance_centroid.nsmallest( n_close ) ) ] ,
                        # pos_max_out[ pos_max_out.distance_max.isin( pos_max_out.distance_max.nsmallest( n_close ) ) ] ]
    
    # stacked_coords = pd.concat( proximity_coords )
    stacked_coords = pd.concat( proximity_coords )
    unique_transects = np.unique( stacked_coords.transect_num )
    full_pol = Polygon( list( unique_coords_utm[ unique_coords_utm.transect_num.isin( unique_transects ) ][ 'geometry' ] ) )
    # polygons.append( sp.concave_hull( full_pol ).convex_hull )
    polygons.append( full_pol.convex_hull )

combined_polygon = unary_union(polygons)
buffer_nmi = 1852 * 1.25
buffered_polygon = combined_polygon.buffer( buffer_nmi )
buffered_polygon = unary_union(buffered_polygon)
mesh_grid_utm_updated = mesh_grid_utm[ mesh_grid_utm[ 'geometry' ].within( buffered_polygon ) ]
mesh_grid_utm_outside = mesh_grid_utm[ ~ mesh_grid_utm[ 'geometry' ].within( buffered_polygon ) ]

x , y = buffered_polygon.exterior.xy
plt.plot(x, y, color='black', label='Convex Shape', linestyle='-', zorder = 400 , linewidth = 1 )

for transect_num, group in unique_coords_utm.groupby('transect_num'):
    plt.scatter(group.sort_values( [ 'x' , 'y' ] ).geometry.x, group.sort_values( [ 'x' , 'y' ] ).geometry.y, label=f'Transect {int(transect_num)}', linestyle='-', s = 1 , color = 'black' , zorder = 100 )

plt.scatter(mesh_grid_utm_updated.geometry.x, mesh_grid_utm_updated.geometry.y, color='red', label='Mesh Grid Points', s = 0.5 , zorder = 0 )

plt.scatter(mesh_grid_utm_outside.geometry.x, mesh_grid_utm_outside.geometry.y, label='Mesh Grid Points', s = 0.5 , zorder = 0 , color = 'gray' )

plt.show( )




    # Centroid
    coord_centroid = coord_tmp.dissolve( by = 'transect_num' , aggfunc = 'mean' )
    coord_centroid[ 'geometry' ] = coord_centroid.geometry.centroid
    coord_centroid = coord_centroid.reset_index( )
    # Min
    coord_min = coord_tmp.dissolve( by = 'transect_num' , aggfunc = 'min' )
    coord_min[ 'geometry' ] = coord_min.geometry.values[ 0 ]
    coord_min = coord_min.reset_index( )
    # Max
    coord_max = coord_tmp.dissolve( by = 'transect_num' , aggfunc = 'max' )
    coord_max[ 'geometry' ] = coord_max.geometry.centroid
    coord_max = coord_max.reset_index( )

    coord_otmp[ 'distance_max' ] = ccoord_otmp.geometry.apply( lambda g: coord_max[ 'geometry' ].distance( g ) )
    coord_otmp[ 'distance_centroid' ] = coord_otmp.geometry.apply( lambda g: coord_centroid[ 'geometry' ].distance( g ) )
    coord_otmp[ 'distance_min' ] = coord_otmp.geometry.apply( lambda g: coord_min[ 'geometry' ].distance( g ) )

    dist_threshold = [ coord_otmp[ coord_otmp[ 'distance_min' ].isin( coord_otmp[ 'distance_min' ].nsmallest( n_close ) ) ] ,
                       coord_otmp[ coord_otmp[ 'distance_centroid' ].isin( coord_otmp[ 'distance_max' ].nsmallest( n_close ) ) ] ,
                       coord_otmp[ coord_otmp[ 'distance_max' ].isin( coord_otmp[ 'distance_max' ].nsmallest( n_close ) ) ] ]


    pos = coord_tmp[ coord_tmp.transect_num == transect_num ][ 'geometry' ]
    coord_tmp[ 'distance' ] = coord_centroid.geometry.apply( lambda g: pos.distance( g ) )
    closest_pos = coord_tmp[ coord_tmp[ 'distance' ].isin( coord_tmp[ 'distance' ].nsmallest( n_close ) ) ]
    full_pol = Polygon( list( unique_coords_utm[ unique_coords_utm.transect_num.isin( closest_pos.transect_num ) ][ 'geometry' ] ) )
    polygons.append( full_pol.convex_hull )


from shapely.geometry import Point, Polygon
from shapely.ops import cascaded_union , unary_union
from shapely.geometry import box
unique_coords_utm.groupby( [ 'transect_num' ] )

### West flank
west_x = unique_coords_utm.groupby( [ 'transect_num' ] )[ 'x' ].min( ).reset_index( )
west_pos = unique_coords_utm[ unique_coords_utm.x.isin( west_x.x ) ].reset_index( drop = True )
east_x = unique_coords_utm.groupby( [ 'transect_num' ] )[ 'x' ].max( ).reset_index( )
east_pos = unique_coords_utm[ unique_coords_utm.x.isin( east_x.x ) ].reset_index( drop = True )
# Interpolate the upper boundary
max_y_per_x = unique_coords_utm.groupby(['transect_num'])['y'].max()
mean_x_per_y = pd.DataFrame( unique_coords_utm.groupby(['transect_num'])['x'].min() )
mean_x_per_y[ 'y' ] = max_y_per_x
mean_x_per_y = mean_x_per_y.reset_index( ).set_index( [ 'x' ] )
# Sort values by x-coordinate for interpolation
sorted_indices = np.argsort(mean_x_per_y.index)
x_values_sorted = mean_x_per_y.index.values[sorted_indices]
y_values_sorted = mean_x_per_y[ 'y' ].values[sorted_indices]

# Interpolate upper boundary (assuming linear interpolation)
upper_boundary_interpolated = np.interp(unique_coords_utm['x'], x_values_sorted, y_values_sorted)
north_pos = pd.DataFrame(
    {
        'x': x_values_sorted  ,
        'y': y_values_sorted ,
    } ,
)

plt.plot( north_pos.sort_values( [ "x" , "y" ] ).x , north_pos.sort_values( [ "x" , "y" ] ).y )
plt.show( )


test[ 'diff' ].max( )
west_pos[ west_pos.geometry.x >= east_pos.geometry.x ]

y_resolution = 1.25 # nmi
y_resolution_m = 1.25 * 1852
y_array = np.arange( unique_coords_utm.geometry.y.min( ) , unique_coords_utm.geometry.y.max( ) , y_resolution_m )

west_pos = west_pos.sort_values( by = [ 'y' ] , ascending = True )
np.all( np.diff( west_pos.y ) > 0 )
west_pos_interp = np.interp( y_array , west_pos.y , west_pos.x )

east_pos_interp = np.interp( y_array , east_pos.y , east_pos.x )

np.arange( east_pos.geometry.y.values )

plt.plot(west_pos.geometry.x, west_pos.geometry.y, color='black', label='Convex Shape', linestyle='--')
plt.plot(east_pos.geometry.x, east_pos.geometry.y, color='red', label='Convex Shape', linestyle='--')

xs = [point.x for point in west_pos.geometry]
ys = [point.y for point in west_pos.geometry]
plt.scatter(xs, ys , zorder = 100 )

xs = [point.x for point in east_pos.geometry]
ys = [point.y for point in east_pos.geometry]
plt.scatter( xs , ys , color = 'red' , zorder = 200 )

plt.show( )

import shapely as sp
from shapely.ops import triangulate
poly1 = Polygon( [ [ p.x , p.y ] for p in west_pos.sort_values( by = 'transect_num' ).geometry ] )
poly2 = Polygon( [ [ p.x , p.y ] for p in east_pos.sort_values( by = 'transect_num' ).geometry ] )
poly_new_df = pd.concat( [ west_pos.sort_values( by = 'transect_num' ) ,
                           east_pos.sort_values( by = 'transect_num' , ascending = False ) ,
                           west_pos.sort_values( by = 'transect_num' ).head( 1 ) ] )
poly3 = Polygon( [ [ p.x , p.y ] for p in poly_new_df.geometry ] )
from scipy.spatial import Delaunay
import geopandas as gpd
import matplotlib.pyplot as plt

# Assuming 'unique_coords_utm' is your GeoDataFrame containing the points

# Extract x and y coordinates from the 'geometry' column
x_coords = unique_coords_utm.geometry.apply(lambda point: point.x)
y_coords = unique_coords_utm.geometry.apply(lambda point: point.y)

# Perform Delaunay triangulation
tri = Delaunay(list(zip(x_coords, y_coords)))

# Plot the triangulation
gdf = gpd.GeoDataFrame(geometry=[unique_coords_utm.geometry], crs=unique_coords_utm.crs)
gdf.plot()
plt.triplot(x_coords, y_coords, tri.simplices)
plt.show()

import geopandas as gpd
import matplotlib.pyplot as plt
import triangle

# Assuming 'unique_coords_utm' is your GeoDataFrame containing the points

# Extract x and y coordinates from the 'geometry' column
x_coords = unique_coords_utm.geometry.apply(lambda point: point.x)
y_coords = unique_coords_utm.geometry.apply(lambda point: point.y)

# Perform Delaunay triangulation
points = dict(vertices=list(zip(x_coords, y_coords)))
triangles = triangle.triangulate(points, 'q')

# Plot the triangulation
gdf = gpd.GeoDataFrame(geometry=[unique_coords_utm.geometry], crs=unique_coords_utm.crs)
gdf.plot()
plt.triplot(x_coords, y_coords, triangles['vertices'], triangles['triangles'])
plt.show()


from shapely.ops import voronoi_diagram
from shapely.geometry import MultiPoint
multi_point = MultiPoint(np.array(unique_coords_utm['geometry'].values))
regions = voronoi_diagram( multi_point )
multi_point.buffer(0.5)
regions = sp.delaunay_triangles( multi_point )

x , y = regions.exterior.xy
plt.plot(x, y, color='black', label='Convex Shape', linestyle='--', zorder = 400 , linewidth = 3 )
plt.show( )

for geom in regions.geoms:
    xs , ys = geom.exterior.xy 
    plt.plot(xs, ys, color='black', label='Convex Shape', linestyle='--', zorder = 400)
plt.show( )

unary_union(list(shapely.ops.polygonize(lines))
sp.ops.triangulate( unique_coords_utm.geometry.values[ 0 ] )
x, y = poly1.exterior.xy
plt.plot(x, y, color='red', label='Convex Shape', linestyle='--', zorder = 300)

x, y = poly2.exterior.xy
plt.plot(x, y, color='orange', label='Convex Shape', linestyle='--', zorder = 300)

x, y = poly3.exterior.xy
plt.plot(x, y, color='green', label='Convex Shape', linestyle='--', zorder = 300)

x , y = poly3.exterior.xy
plt.plot(x, y, color='black', label='Convex Shape', linestyle='--', zorder = 400 , linewidth = 3 )
plt.show( )

for geom in sp.union( sp.make_valid( poly3 ) , sp.make_valid( poly2 ) ).geoms:
    xs , ys = geom.exterior.xy 
    plt.plot(xs, ys, color='black', label='Convex Shape', linestyle='--', zorder = 400)



poly_temp = poly2.buffer(0).union( poly1.buffer(0) )
poly_all = unary_union( [ poly1.buffer(0) , poly2.buffer(0) ] )
poly_new_df = pd.concat( [ west_pos.sort_values( by = 'transect_num' ) ,
                           east_pos.sort_values( by = 'transect_num' , ascending = False ) ,
                           west_pos.sort_values( by = 'transect_num' ).head( 1 ) ] )
poly3 = Polygon( [ [ p.x , p.y ] for p in poly_new_df.geometry ] )
poly5 = sp.polygonize( [ [ p.x , p.y ] for p in poly_new_df.geometry ] )

n_close = 5
polygons = []
for transect_num in coord_centroid[ 'transect_num' ]:
    # Extract coordinates and find westernmost and easternmost points
    coord_tmp = coord_centroid.copy( )
    pos = coord_tmp[ coord_tmp.transect_num == transect_num ][ 'geometry' ]
    coord_tmp[ 'distance' ] = coord_centroid.geometry.apply( lambda g: pos.distance( g ) )
    closest_pos = coord_tmp[ coord_tmp[ 'distance' ].isin( coord_tmp[ 'distance' ].nsmallest( n_close ) ) ]
    full_pol = Polygon( list( unique_coords_utm[ unique_coords_utm.transect_num.isin( closest_pos.transect_num ) ][ 'geometry' ] ) )
    polygons.append( full_pol.convex_hull )

combined = unary_union( polygons )

poly4 = sp.union( combined , sp.make_valid( poly3 ) )

unique_coords[ unique_coords.transect_num == 143 ].geometry.y.min()

for transect_num, group in unique_coords_utm.groupby('transect_num'):
    plt.scatter(group.sort_values( by = [ 'x' , 'y' ] ).geometry.x, group.sort_values( by = [ 'x' , 'y' ] ).geometry.y, label=f'Transect {int(transect_num)}', linestyle='-', zorder = 0 , color = 'black' , s = 1 )
plt.scatter( unique_coords_utm[ unique_coords_utm.transect_num == 1 ].sort_values( by = [ 'x' , 'y' ] ).geometry.x , unique_coords_utm[ unique_coords_utm.transect_num == 1 ].sort_values( by = [ 'x' , 'y' ] ).geometry.y , color = 'red' , linewidth = 3.0 )
plt.scatter( west_pos.x , west_pos.y , color = 'red' , zorder = 300 )
plt.scatter( east_pos.x , east_pos.y , color = 'gray' , zorder = 300 )
x, y = unary_union( polygons ).exterior.xy
plt.plot(x, y, color='orange', label='Convex Shape', linestyle='--', zorder = 500)
x, y = poly4.exterior.xy
plt.plot(x, y, color='black', label='Convex Shape', linestyle='--', zorder = 500)
for geom in poly4.geoms:
    xs , ys = geom.exterior.xy 
    plt.plot(xs, ys, color='black', label='Convex Shape', linestyle='--')

plt.show( )

x, y = poly1.exterior.xy
plt.plot(x, y, color='black', label='Convex Shape', linestyle='--', zorder = 300)

x, y = poly2.exterior.xy
plt.plot(x, y, color='black', label='Convex Shape', linestyle='--', zorder = 300)



for geom in poly2.geoms:
    xs , ys = geom.exterior.xy 
    plt.plot(xs, ys, color='red', label='Convex Shape', linestyle='--')

xs, ys = poly_all.exterior.xy
plt.plot(xs, ys, color='red', label='Convex Shape', linestyle='--', zorder = 300)

for geom in poly_all.geoms:
    xs , ys = geom.exterior.xy 
    plt.plot(xs, ys, color='red', label='Convex Shape', linestyle='--')

xx, yy = poly_temp.exterior.xy
plt.plot(xx, yy, color='red', label='Convex Shape', linestyle='--', zorder = 300)

for geom in poly_temp.geoms:
    xs , ys = geom.exterior.xy 
    plt.plot(xs, ys, color='red', label='Convex Shape', linestyle='--')

plt.show( )

plt.plot( west_pos_interp , y_array , color = 'gray' , zorder = 300 )
plt.plot( east_pos_interp , y_array , color = 'gray' , zorder = 300 )

x, y = buffered_polygon.exterior.xy
plt.plot(x, y, color='black', label='Convex Shape', linestyle='--', zorder = 300)



plt.show( )

plt.plot(gpd.GeoSeries(east_pos.geometry))

west_pos.geometry.x , west_pos.geometry.y

coord_centroid = unique_coords_utm.dissolve( by = 'transect_num' , aggfunc = 'mean' )
coord_centroid[ 'geometry' ] = coord_centroid.geometry.centroid
coord_centroid = coord_centroid.reset_index( )
n_close = 2


polygons = []
for transect_num, group in unique_coords_utm.groupby('transect_num'):
    min_lon = group.geometry.x.min()
    max_lon = group.geometry.x.max()
    min_lat = group.geometry.y.min()
    max_lat = group.geometry.y.max()

    buffer_nmi = 1852 * buffer

    polygon = Polygon(list([(min_lon - buffer_nmi, min_lat - buffer_nmi),
                      (max_lon + buffer_nmi, min_lat - buffer_nmi),
                      (max_lon + buffer_nmi, max_lat + buffer_nmi),
                      (min_lon - buffer_nmi, max_lat + buffer_nmi)]))
    polygons.append(polygon.convex_hull)

combined_polygon = unary_union(polygons)
buffer_nmi = 1852 * 0
buffered_polygon = combined_polygon.buffer( buffer_nmi )
mesh_grid_utm_updated = mesh_grid_utm[ mesh_grid_utm[ 'geometry' ].within( combined_polygon ) ]

import matplotlib.pyplot as plt

# Plot transect lines
for transect_num, group in unique_coords_utm.groupby('transect_num'):
    plt.plot(group.geometry.x, group.geometry.y, label=f'Transect {int(transect_num)}', linestyle='-', marker='o')

plt.scatter(mesh_grid_utm_updated.geometry.x, mesh_grid_utm_updated.geometry.y, color='red', label='Mesh Grid Points')


x, y = combined_polygon.exterior.xy
plt.plot(x, y, color='black', label='Convex Shape', linestyle='--')
plt.show( )

from shapely import concave_hull

n_close = 5
polygons = []
for transect_num in coord_centroid[ 'transect_num' ]:
    # Extract coordinates and find westernmost and easternmost points
    coord_tmp = coord_centroid.copy( )
    pos = coord_tmp[ coord_tmp.transect_num == transect_num ][ 'geometry' ]
    coord_tmp[ 'distance' ] = coord_centroid.geometry.apply( lambda g: pos.distance( g ) )
    closest_pos = coord_tmp[ coord_tmp[ 'distance' ].isin( coord_tmp[ 'distance' ].nsmallest( n_close ) ) ]
    full_pol = Polygon( list( concave_hull( unique_coords_utm[ unique_coords_utm.transect_num.isin( closest_pos.transect_num ) ][ 'geometry' ] ) ) )
    polygons.append( full_pol.convex_hull )

combined_polygon = unary_union(polygons)
buffer_nmi = 1852 * 1.25
buffered_polygon = combined_polygon.buffer( buffer_nmi )
temp = concave_hull( buffered_polygon.union( buffered_polygon ).convex_hull )
mesh_grid_utm_updated = mesh_grid_utm[ mesh_grid_utm[ 'geometry' ].within( buffered_polygon ) ]

import matplotlib.pyplot as plt

# Plot transect lines
for transect_num, group in unique_coords_utm.groupby('transect_num'):
    plt.plot(group.geometry.x, group.geometry.y, label=f'Transect {int(transect_num)}', linestyle='-', marker='o')

plt.scatter(mesh_grid_utm_updated.geometry.x, mesh_grid_utm_updated.geometry.y, color='red', label='Mesh Grid Points')

for geom in buffered_polygon.geoms:
    xs , ys = geom.exterior.xy 
    plt.plot(xs, ys, color='black', label='Convex Shape', linestyle='--')

x , y = buffered_polygon.exterior.xy
plt.plot(x, y, color='black', label='Convex Shape', linestyle='--')

# for geom in temp.geoms:
#     xs , ys = geom.exterior.xy 
#     plt.plot(xs, ys, color='red', label='Convex Shape', linestyle='--')

x , y = temp.exterior.xy
plt.plot(x, y, color='blue', label='Convex Shape', linestyle='--')

plt.show( )

bbox = combined_polygon.bounds
bbox_polygon = box(*bbox)

# Filter mesh grid points using the bounding box
mesh_grid_filtered_bbox = mesh_grid_utm[mesh_grid_utm.geometry.within(bbox_polygon)]

# Filter mesh grid points using spatial indexing and the combined polygon
spatial_index = mesh_grid_filtered_bbox.sindex
possible_matches_index = list(spatial_index.intersection(combined_polygon.bounds))
possible_matches = mesh_grid_filtered_bbox.iloc[possible_matches_index]
final_matches = possible_matches[possible_matches.intersects(combined_polygon)]
mesh_grid_filtered = mesh_grid_utm[mesh_grid_utm.geometry.within(combined_polygon)]

import matplotlib.pyplot as plt

# Plot transect lines
for transect_num, group in unique_coords_utm.groupby('transect_num'):
    plt.plot(group.geometry.x, group.geometry.y, label=f'Transect {int(transect_num)}', linestyle='-', marker='o')

# Plot mesh grid points
plt.scatter(mesh_grid.geometry.x, mesh_grid.geometry.y, color='red', label='Mesh Grid Points')

# Plot the outline of the convex shape
x, y = combined_polygon.exterior.xy
plt.plot(x, y, color='black', label='Convex Shape', linestyle='--')

# Add legend and labels
plt.legend()
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Transect Lines, Mesh Grid Points, and Convex Shape')

# Show plot
plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()
# Filter mesh grid points that fall within the combined polygon
mesh_grid_filtered = mesh_grid[mesh_grid.geometry.within(combined_polygon)]

transect_trans_df , d_x , d_y = transform_geometry( data_regrid , 
                                                    self.statistics[ 'kriging' ][ 'isobath_200m_df' ] , 
                                                    self.config[ 'kriging_parameters' ] ,
                                                    self.config[ 'geospatial' ][ 'init' ] )



overall = self.biology[ 'population' ][ 'biomass' ][ 'biomass_age_df' ]
overall_nonstrat = overall[ overall.stratum_num != 0 ]
overall_nonstrat[ overall.sex == 'female' ][ 'B_age' ].sum( ) * 1e-6

spatial_data = spatial_data_transformed
transformed_mesh = self.statistics[ 'kriging' ][ 'mesh_df' ].copy( )
dataframe_mesh = self.statistics[ 'kriging' ][ 'mesh_df' ].copy( )
dataframe_geostrata = self.spatial[ 'geo_strata_df' ].copy( )
variogram_parameters = variogram_parameters
kriging_parameters = kriging_parameters
self.statistics[ 'kriging' ][ 'isobath_200m_df' ]
transformed_mesh[ 'x_transformed' ].min( )
transformed_mesh[ 'x_transformed' ].max( )
spatial_data[ 'longitude' ].min( )
transformed_mesh[ 'centroid_longitude' ].min( )
transformed_mesh[ 'centroid_longitude' ][ transformed_mesh[ 'centroid_longitude' ] > spatial_data[ 'longitude' ].min( ) ]
spatial_data.longitude.min( )
latitude_bins = np.concatenate( [ [ -90.0 ] , dataframe_geostrata.northlimit_latitude , [ 90.0 ] ] )
dataframe_mesh[ 'stratum_num' ] = pd.cut( dataframe_mesh.centroid_latitude ,
                                            latitude_bins ,
                                            labels = list( dataframe_geostrata.stratum_num ) + [ 1 ] ,
                                            ordered = False )

variable = 'B_a_adult'

mesh_grid = transformed_mesh
transformed_mesh[ 'centroid_latitude' ].min( )
k_max = kriging_parameters[ 'kmax' ]

distance_matrix = griddify_lag_distances( mesh_grid , spatial_data )

### Calculate the kriging distance matrix and corresponding indices
distance_matrix , local_point_grid  = local_search_index( transformed_mesh , 
                                                            spatial_data , 
                                                            kriging_parameters[ 'kmax' ] )