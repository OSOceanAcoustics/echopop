from EchoPro.survey import Survey
import numpy as np
import pandas as pd
import geopandas as gpd
import copy
from pathlib import Path
import os
from EchoPro.core import CONFIG_MAP, LAYER_NAME_MAP
from EchoPro.utils.data_file_validation import load_configuration , validate_data_columns
from EchoPro.computation.operations import bin_variable , bin_stats , count_variable
from EchoPro.utils.data_file_validation import validate_data_columns
from EchoPro.computation.acoustics import to_linear , ts_length_regression
from EchoPro.computation.operations import group_merge
os.getcwd()

init_config_path=Path('./config_files/initialization_config.yml')
survey_year_config_path=Path('./config_files/survey_year_2019_config.yml')
obj = Survey( init_config_path='./config_files/initialization_config.yml' , survey_year_config_path='./config_files/survey_year_2019_config.yml' )






obj.transect_analysis()
obj.stratified_summary()
obj.standardize_coordinates()
self = obj
species_id = 22500

################# Distribute biological variables to kriging data
### Reformat 'specimen_df' to match the same format as 'len_df'
# First make copies of each
specimen_df_copy = self.biology['specimen_df'].copy().pipe( lambda df: df.loc[ df.species_id == species_id ] )
length_df_copy = self.biology['length_df'].copy().pipe( lambda df: df.loc[ df.species_id == species_id ] )

# Import length bins
length_intervals = self.biology[ 'distributions' ][ 'length' ][ 'length_interval_arr' ]
length_bins = self.biology[ 'distributions' ][ 'length' ][ 'length_bins_arr' ]
age_intervals = self.biology[ 'distributions' ][ 'age' ][ 'age_interval_arr' ]

# Import haul catch numbers
haul_catch_df = self.biology[ 'catch_df' ].copy().pipe( lambda df: df.loc[ df.species_id == species_id ] )

# Remove hauls not found in `length_df_copy`
haul_catch_filtered = (
    haul_catch_df
    .loc[ lambda x: x.haul_num.isin( length_df_copy.haul_num ) ]
)

### Calculate station weights -- stratify by stratum
weight_strata_station = (
    pd.concat( [
        haul_catch_filtered 
        .groupby( 'stratum_num' )
        .apply( lambda df: pd.Series( { 'stratum_weight': df.haul_weight.sum( ) } ) )
        .reset_index( )
        .assign( station = 1 ) , 
        specimen_df_copy
        .groupby( [ 'stratum_num' ] )
        .apply( lambda df: pd.Series( { 'stratum_weight': df.weight.sum( ) } ) )
        .reset_index()
        .assign( station = 2 )
    ] )
)

weight_strata = (
    weight_strata_station
    .groupby( [ 'stratum_num' ] )
    .apply( lambda x: pd.Series( { 'weight_stratum_total': x.stratum_weight.sum( ) } ) )            
    .reset_index()
)

### Calculate summed length-age-sex-stratum weight bins
length_age_sex_stratum_weight = (
    specimen_df_copy
    .dropna( how = 'any' )
    .loc[ lambda x: x.sex != 'unsexed' ]
    .bin_variable( age_intervals , 'age' )
    .bin_variable( length_intervals , 'length' )      
    .groupby( [ 'stratum_num' , 'age_bin' , 'length_bin' , 'sex' ] )
    .apply( lambda df: pd.Series( { 'summed_weight_all': df.weight.sum( ) ,
                                    'summed_weight_adult': df.loc[ df.age > 1 ].weight.sum( ) } ) ) 
    .reset_index( )
    .replace( np.nan , 0 )
    .assign( total_weight_sex_all = lambda df: df.groupby( [ 'stratum_num' , 'sex' ] )[ 'summed_weight_all' ].transform( sum ) ,
                total_weight_sex_adult = lambda df: df.groupby( [ 'stratum_num' , 'sex' ] )[ 'summed_weight_adult' ].transform( sum ) ,
                proportion_weight_all = lambda df: df.summed_weight_all / df.total_weight_sex_all ,
                proportion_weight_adult = lambda df: df.summed_weight_adult / df.total_weight_sex_adult )
)

# Append the total weight to calculate the aged-to-unaged weight ratio/proportion
aged_proportions = (
    length_age_sex_stratum_weight
    .merge( weight_strata , on = [ 'stratum_num' ] )     
    .assign( proportion_weight_aged_all = lambda x: x.summed_weight_all / x.weight_stratum_total ,
                proportion_weight_unaged_all = lambda x: 1.0 - x.proportion_weight_aged_all ,
                proportion_weight_aged_adult = lambda x: x.summed_weight_adult / x.weight_stratum_total ,
                proportion_weight_unaged_adult = lambda x: 1.0 - x.proportion_weight_aged_adult )
)

### Calculate interpolated weights based on length bins for each sex per haul
# Extract length from binned values 
length_weight_fit = (
    self.statistics[ 'length_weight' ][ 'length_weight_df' ]
    .assign( length_bin_value  = lambda x: x[ 'length_bin' ].apply( lambda y: y.mid ) )
)

haul_weights = (
    length_df_copy
    .bin_variable( length_intervals , 'length' )
    .loc[ lambda x: x.group == 'sexed' ]
    .pipe( lambda df: pd.concat( [ df , df.assign( sex = 'all' ) ] ) )  
    .merge( length_weight_fit , on = [ 'sex' , 'length_bin' ] )
    .groupby( [ 'stratum_num' , 'haul_num' , 'sex' ] )
    .apply( lambda df: pd.Series( { 'weight_interp': ( np.interp( df.length_bin_value , 
                                                                  length_weight_fit.loc[ lambda x: x.sex.isin( df.sex ) ][ 'length_bin_value' ] , 
                                                                  length_weight_fit.loc[ lambda x: x.sex.isin( df.sex ) ][ 'weight_modeled' ] ) *
                                                          df.length_count ).sum( ) } ) )
    .groupby( [ 'stratum_num' , 'sex' ] )
    .apply( lambda df: pd.Series( { 'summed_haul_weights': df.weight_interp.sum( ) } ) )
    .reset_index( )
)

### Some normalization
(
    haul_weights
    .merge( weight_strata_station.loc[ lambda x: x.station == 1 ] , on = 'stratum_num' )
    .assign( normalized_weight = lambda x: )
)
haul_weights.merge( weight_strata_station.loc[ lambda x: x.station == 1 ] , on = 'stratum_num' )



### Kriged data input
dataframe_input = self.statistics[ 'kriging' ][ 'kriged_biomass_df' ]

### Need aged/unaged biomass for males/females 
# First need the length-age-sex proportions
self.biology[ 'weight' ][ 'proportions' ][ 'length_sex_age_weight_proportions_df' ].loc[ lambda x: x[ 'count' ] > 0.0 ]
self.biology[ 'weight' ][ 'weight_strata_df' ]
(
    self.biology[ 'weight' ][ 'proportions' ][ 'length_sex_age_weight_proportions_df' ]
    .assign( prop_sex_sum = lambda x: x.groupby( [ 'stratum_num' , 'sex' ] )[ 'count' ].transform( sum ) ,
             prop_holder = lambda x: np.where( x.sex == 'all' , x[ 'count' ] , 0.0 ) ,
             prop_strat_sum = lambda x: x.groupby( [ 'stratum_num' ] )[ 'prop_holder' ].transform( sum ) ,
             prop_sex_norm = lambda x: x[ 'count' ] / x.prop_sex_sum ,
             prop_strat_norm = lambda x: x[ 'count' ] / x.prop_strat_sum )
    .groupby( [ 'stratum_num' , 'sex' ] )
    .apply( lambda df: pd.Series( { 'weight': ( df.prop_sex_norm * df.prop_strat_norm ).sum( ) } ) )
)
x.groupby( [ 'stratum_num' ] )[ 'count' ]
self.biology[ 'weight' ][ 'proportions' ][ 'sex_age_weight_proportions_df' ].loc[ lambda x: x.weight_sex_proportion_adult > 0 ]
(
    self.biology[ 'weight' ][ 'proportions' ][ 'sex_age_weight_proportions_df' ]
    .loc[ lambda df: ( df.age > 1 ) & ( df.sex == 'all' ) ]
)


total_weights = (
    self.biology[ 'weight' ][ 'proportions' ][ 'length_sex_age_weight_proportions_df' ]
    .loc[ lambda df: ( df.age > 1 ) & ( df.sex == 'all' ) ]
    .drop_duplicates( subset = [ 'stratum_num' , 'sex' ] )
)
    .groupby( [ 'stratum_num' ] )
    .apply( lambda df: pd.Series( { 'stratum_total_weight': np.nansum( df.weight_total_adult ) } ) )
)

dist_weight_sum = (
    self.biology[ 'weight' ][ 'proportions' ][ 'length_sex_age_weight_proportions_df' ]
    .merge( total_weights , on = 'stratum_num' )
)

self.biology[ 'weight' ][ 'proportions' ][ 'length_sex_age_weight_proportions_df' ]




.apply( lambda df: pd.Series( {
    'count_age_proportion_all': ( df[ 'count' ] / df.stratum_count_all ).sum() ,
    'count_age_proportion_adult': ( df.loc[ df.age > 1 ][ 'count' ] / df.stratum_count_total ).sum( )
} ) )


import matplotlib.pyplot as plt


l1 = length_df_copy.groupby('stratum_num').apply( lambda x: x[ 'length_count' ].sum()).reset_index(name='count')
s1 = specimen_df_copy.groupby('stratum_num').size().reset_index(name='count')

ma = l1.merge( s1 , on = 'stratum_num' , how = 'outer' ).replace(np.nan, 0).assign( total = lambda x: x.count_x + x.count_y )

plt.bar( ma.stratum_num , ma.total )
plt.xlabel( 'KS stratum' )
plt.ylabel( 'Frequency' )
plt.show( )

l2 = length_df_copy.groupby(['stratum_num' , 'transect_num']).apply( lambda x: x[ 'length_count' ].sum()).reset_index(name='count')
s2 = specimen_df_copy.groupby(['stratum_num' , 'transect_num']).size().reset_index(name='count')
m2 = l2.merge( s2 , on = ['stratum_num' , 'transect_num'] , how = 'outer' ).replace(np.nan, 0).assign( total = lambda x: x.count_x + x.count_y )
m2 = m2.sort_values( 'transect_num' )
m2.loc[ m2.transect_num > 100 ]
for stratum_num in np.unique( m2.stratum_num ):
    subset = m2[m2['stratum_num'] == stratum_num]
    plt.plot(subset['transect_num'], subset['total'], label=f'Stratum {int(stratum_num)}')
    plt.scatter(subset['transect_num'], subset['total'], marker='o')
plt.xlabel('Transect Number')
plt.ylabel('Total')
plt.title('Total vs Transect Number for Different Stratum')
plt.legend()
plt.grid(True)
plt.show()

(
    nasc_interval_df
    # Merge stratum information ( join = 'outer' since missing values will be filled later on)
    .merge( info_strata , on = [ 'stratum_num' , 'haul_num' ] , how = 'outer' )
    # Drop unused hauls
    .dropna( subset = 'transect_num' )
    # Fill NaN w/ 0's for 'fraction_hake'
    .assign( fraction_hake = lambda x: x[ 'fraction_hake' ].fillna( 0 ) )
    # Group merge
    .group_merge( dataframes_to_add = dataframes_to_add , on = 'stratum_num' )
    )
)



a = age_sex_weight_proportions.pivot_table( index = [ 'age' , 'sex' ] , columns = [ 'stratum_num' ] )
a.loc[ : , ( 'weight_age_sex_proportion_all' ) ]
age_sex_weight_proportions.pivot_table( index = [ 'sex' , 'station' ] , 
                          columns = [ 'stratum_num' ] , 
                          values = [ 'overall_station_p' ] )


            age_proportions
            .pipe( lambda df: df
                .groupby( 'stratum_num' )
                .apply( lambda x: 1 - x.loc[ x[ 'age' ] < 2 ][ 'count' ].sum() / x[ 'count' ].sum( ) ) ) 
            .reset_index( name = 'adult_proportion' )
        )

 age_proportions = (
            specimen_df_copy
            .assign( age_group = lambda x: np.where( x.age < 2 , 'nonadult' , 'adult' ) )
            .dropna( how = 'any' )
            .count_variable( variable = 'length' ,
                             contrasts = [ 'stratum_num' , 'age' , 'age_group' ] ,
                             fun = 'size' )
            .pipe( lambda x: x.assign( stratum_count_all = x.groupby( [ 'stratum_num' ] )[ 'count' ].transform( sum ) ,
                                       stratum_proportion_all = lambda x: x[ 'count' ] / x[ 'stratum_count_all' ] ,
                                       stratum_count_adult = x.loc[ x.age_group == 'adult' ].groupby( [ 'stratum_num' ] )[ 'count' ].transform( sum ).astype( int ) ,
                                       stratum_proportion_adult = lambda x: x[ 'count' ] / x[ 'stratum_count_adult' ] ) ) 
        )

plt.bar( ma.stratum_num , ma.total )
plt.xlabel( 'KS stratum' )
plt.ylabel( 'Frequency' )
plt.show( )

plt.boxplot( specimen_df_copy.groupby('stratum_num').size() , specimen_df_copy.groupby('stratum_num').size().values )
plt.show()
### Reformat 'specimen_df' to match the same format as 'len_df'
# First make copies of each
specimen_df_copy = self.biology['specimen_df'].copy().pipe( lambda df: df.loc[ df.species_id == species_id ] )
length_df_copy = self.biology['length_df'].copy().pipe( lambda df: df.loc[ df.species_id == species_id ] )

self.biology['specimen_df'].copy().loc[ lambda x: x.age == 0]
self.biology['length_df'].copy().loc[ lambda x: x.age == 0]
### Pull length distribution values
length_intervals = self.biology['distributions']['length']['length_interval_arr']   

### Calculate the sex proportions/frequencies for station 1 (length_df) across all strata
specimen_df_copy.loc[ lambda x: x.sex == 'male' ]
        
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
    .assign( group = lambda x: np.where( x[ 'sex' ] == int( 1 ) , 'male' , np.where( x[ 'sex' ] == int( 2 ) , 'female' , 'unsexed' ) ) ) # assigns str variable for comprehension
    .pipe( lambda df: pd.concat( [ df.loc[ df[ 'sex' ] != 3 ] , df.assign( group = 'all' ) ] ) ) # appends male-female to an 'all' dataframe
    .assign( station = 2 ) # assign station number for later functions
)

aa = (
    specimen_df_spp
    .dropna( subset = [ 'weight' , 'length' ] )
    .assign( group = lambda x: np.where( x[ 'sex' ] == int( 1 ) , 'male' , 
                                np.where( x[ 'sex' ] == int( 2 ) , 'female' , 'unsexed' ) ) ) # assigns str variable for comprehension
    .pipe( lambda df: pd.concat( [ df.loc[ df[ 'group' ] != 'unsexed'  ] , df.assign( group = 'all' ) ] ) ) # appends male-female to an 'all' dataframe
)

np.unique( aa.group )

aa.loc[ aa.group == 'all' ].shape
specimen_sex_grouped_df.bin_stats( bin_variable = 'length' , bin_values = length_intervals , contrasts = 'group' ).loc[ lambda x: x.group == 'female' ]

len_wgt_all = survey.bio_calc._generate_length_val_conversion(
    len_name="length",
    val_name="weight",
    df=survey.bio_calc.specimen_df,
)

(
    fitted_weights
    .rename( columns = { 'length_bin': 'length' } )
    .bin_variable( bin_values = length_intervals ,
                   bin_variable = 'length' )
    .merge( length_bin_stats , on = [ 'group' , 'length_bin' ] )
    .assign( weight_modeled = lambda x: np.where( x.n_weight < 5 ,
                                                  x.fitted_weight ,
                                                  x.mean_weight ) )
)

specimen_df_copy = self.biology['specimen_df'].copy().pipe( lambda df: df.loc[ df.species_id == species_id ] )

# Import length bins
length_intervals = self.biology[ 'distributions' ][ 'length' ][ 'length_interval_arr' ]
age_intervals = self.biology[ 'distributions' ][ 'age' ][ 'age_interval_arr' ]
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

age_proportions.groupby( [ 'age' , 'stratum_num' ] ).apply( lambda x: x[ 'stratum_proportion' ].sum() ).reset_index().loc[ lambda x: x.stratum_num == 1.0 ]

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
    .dropna( subset = [ 'length' , 'weight' ] )
    .pipe( lambda x: x.assign( weight_age = x.groupby( [ 'stratum_num' , 'age' ] )[ 'weight' ].transform( sum ) ,
                               weight_stratum = x.groupby( [ 'stratum_num' ] )[ 'weight' ].transform( sum ) ) )
    .groupby( [ 'stratum_num' , 'age' ] )
    .apply( lambda x: x[ 'weight_age' ].sum( ) / x[ 'weight_stratum' ].sum ( ) )
    .reset_index( name = 'weight_stratum_proportion' )
)

a = (
    specimen_df_copy
    .dropna( how = "any" )
    .bin_variable( bin_values = length_intervals ,
                   bin_variable = 'length' )
    .bin_variable( bin_values = age_intervals ,
                   bin_variable = 'age' )
    .assign( group = lambda x: np.where( x[ 'sex' ] == int( 1 ) , 'male' , 
                                        np.where( x[ 'sex' ] == int( 2 ) , 'female' , 'unsexed' ) ) ) # assigns str variable for comprehension
    .pipe( lambda df: pd.concat( [ df.loc[ df[ 'group' ] != 'unsexed'  ] , df.loc[ df[ 'group' ] != 'unsexed'  ].assign( group = 'all' ) ] ) ) # appends male-female to an 'all' dataframe
    .bin_stats( bin_variable = 'age' ,
                bin_values = age_intervals ,
                variables = 'weight' ,
                contrasts = [ 'group' , 'stratum_num' , 'length_bin' ] ,
                functions = ['mean', 'size' , 'sum' ] )
    # .groupby( [ 'group' ] )
    # .apply( lambda x: x.mean)
)

a.loc[ ( a.stratum_num == 1 ) & ( a.group == 'male' ) & ( a.sum_weight > 0 ) ][ 'sum_weight' ]
a.loc[ ( a.stratum_num == 1 ) & ( a.group == 'all' ) & ( a.sum_weight > 0 ) ][ 'sum_weight' ].sum( )
a.loc[ ( a.stratum_num == 1 ) & ( a.group == 'all' ) ].groupby( 'age_bin' ).apply( lambda x: x[ 'sum_weight' ].sum())

(
    specimen_df_copy
    .dropna( subset = [ 'weight' , 'length' ] )
    .pipe( lambda x: x.assign( weight_age = x.groupby( [ 'stratum_num' , 'age' ] )[ 'weight' ].transform( sum ) ,
                               weight_stratum = x.groupby( [ 'stratum_num' ] )[ 'weight' ].transform( sum ) ) )
).loc[ lambda x: ( x.stratum_num == 1 ) ]

b = age_weight_proportions.groupby( [ 'age' , 'stratum_num' ] ).apply( lambda x: x[ 'weight_stratum_proportion' ].sum() ).reset_index().loc[ lambda x: x.stratum_num == 1.0 ]
c = b.rename( columns = { 0: 'prop' } ).loc[ b.age != 0 ][ 'prop' ].sum()
b.rename( columns = { 0: 'prop' } ).loc[ b.age != 0 ].assign( prop_new = lambda x: x.prop / c )
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


fitted_weights.bin_variable( bin_values = length_intervals , 
                             bin_variable = 'length_bin' )
length_bin_stats.merge( fitted_weights , on = [ 'group' , 'length_bin' ] )

length_bin_stats.assign( weight_modeled = lambda x: np.wherE( x.n_weight < 5 , ))

length_bin_stats[ 'weight_modeled' ] = length_bin_stats[ 'mean_weight' ].copy()
low_n_indices = np.where(length_bin_stats[ 'n_weight' ].values < 5)[0].copy()   
length_bin_stats.loc[low_n_indices, 'weight_modeled'] = fitted_weight[low_n_indices].copy()

a1 = (
    specimen_df_spp    
    .assign( group = lambda x: np.where( x[ 'sex' ] == int( 1 ) , 'male' , 
                               np.where( x[ 'sex' ] == int( 2 ) , 'female' , 'unsexed' ) ) ) # assigns str variable for comprehension
    .pipe( lambda df: pd.concat( [ df.loc[ df[ 'group' ] != 'unsexed'  ] , df.assign( group = 'all' ) ] ) ) # appends male-female to an 'all' dataframe
    .dropna( subset = [ 'weight' , 'length' ] )
    .groupby( 'group' )
    .apply( lambda x: pd.Series( np.polyfit( np.log10( x[ 'length' ].values ) , np.log10( x[ 'weight' ].values ) , 1 ) ,
                                index = [ 'rate' , 'initial' ] ) )
    .reset_index()
)

b = a.merge( pd.DataFrame({'length_bins': length_bins}) , how = 'cross' ).assign( fitted_weight = lambda x: 10 ** x.initial * x.length_bins.astype( int ) ** x.rate )
b.loc[ b.group == 'all' ][ 'fitted_weight' ][ 5 ]
a.groupby( 'group' ).apply( lambda x: 10 ** x.initial * int( length_bins ) ** x.rate , axis = 1 )
a.rate[ 0 ]
a.reset_index()

specimen_df_spp.bin_stats( bin_variable = 'length' , bin_values = length_intervals )


specimen_df_spp.assign(group = lambda x: x.sex.astype(str)).loc[lambda x: ( x.group == '1' ) ].apply(lambda x: pd.Series( np.polyfit( np.log10( x[ 'length' ] ) , np.log10( x[ 'weight' ] ) , 1 ) , index=['slope', 'intercept'] ) )
df = specimen_df_spp.assign(group = lambda x: x.sex.astype(str)).loc[lambda x: ( x.group == '1' ) ].dropna(how='any')
np.polyfit( np.log10( df[ 'length' ] ) , np.log10( df[ 'weight' ] ) , 1 )


specimen_df_spp.assign(group = lambda x: str(x.sex)).group

length_weight_df.bin_stats( bin_variable = 'length' , bin_values = length_intervals , contrasts = [ 'sex' ] )


### Calculate the sex proportions/frequencies for station 2 (specimen_df) across all strata
specimen_grouped = (
    specimen_df_copy
    .bin_variable( length_intervals , 'length' ) # appends `length_bin` column
    .assign( group = lambda x: np.where( x[ 'sex' ] == int( 1 ) , 'male' , np.where( x[ 'sex' ] == int( 2 ) , 'female' , 'unsexed' ) ) # assigns str variable for comprehension
    .pipe( lambda df: pd.concat( [ df.loc[ df[ 'sex' ] != 3 ] , df.assign( group = 'all' ) ] ) ) # appends male-female to an 'all' dataframe
    .assign( station = 2 ) # assign station number for later functions
)

specimen_df_copy.dropna( subset = [ 'weight' , 'length' ] )

(
    specimen_df_copy
    .dropna( subset = 'weight' )
    .bin_variable( length_intervals , 'length' ) # appends `length_bin` column
    .count_variable( contrasts = [  'length_bin' ] , # grouping variables
                        variable = 'length' , # target value to apply a function
                        fun = 'count' ) # function to apply
    . sum( )
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

station_sex_length.loc[ lambda x: ( x.station == 1 ) & ( x.group == 'all' ) ] 
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
fitted_weight = self.statistics['length_weight']['length_weight_df'].loc[ lambda x: x.group == 'all' ]

# Total
total_weighted_values = ( length_proportion_table.loc[ : , ( 'all' , 1 ) ] * station_proportions.loc[ 1 , ] +
                            length_proportion_table.loc[ : , ( 'all' , 2 ) ] * station_proportions.loc[ 2 , ] )
total_weight = fitted_weight[ 'weight_modeled' ].dot( total_weighted_values.reset_index( drop = True ) )
total_weighted_values.loc[ : , 1.0 ].sort_values( ascending = False ).head(3).sum()
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



self.biology[ 'weight' ][ 'weight_strata_df' ]


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

age_sex_weight_proportions.loc[ lambda x: ( x.sex == 'female' ) & ( x.stratum_num == 1.0 ) ]

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
