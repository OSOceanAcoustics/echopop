import pandas as pd
import numpy as np
from EchoPro.computation.acoustics import to_linear , ts_length_regression

### Dummy dataframes
specimen_dataframe = pd.DataFrame(
    {
        'stratum_num': np.repeat( [ 0 , 1 , 2 , 4 , 5 ] , 4 ) ,
        'haul_num': np.repeat( [ 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10 ] , 2 ) ,
        'species_id': np.append( np.repeat( [ 19350 ] , 19 ) , 43130 ) ,
        'length': np.linspace( 10 , 100 , 20 ) ,
        'weight': np.linspace( 1 , 5 , 20 ) ,     
    }
)

length_dataframe = pd.DataFrame(
    {
        'stratum_num': np.repeat( [ 0 , 1 , 2 , 4 , 5 ] , 4 ) ,
        'haul_num': np.repeat( [ 1 , 2 , 3 , 4 , 5 ] , 4 ) ,
        'species_id': np.append( np.repeat( [ 19350 ] , 19 ) , 43130 ) ,
        'length': np.linspace( 10 , 100 , 20 ) ,
        'length_count': np.linspace( 10 , 100 , 20 ) ,     
    }
)

strata_dataframe = pd.DataFrame(
    {
        'stratum_num': [ 0 , 1 , 2 , 3 , 4 , 5 , 6 ]
    }
)

### Dummy parameters
ts_length_parameters = { 'species_code': 19350 ,
                         'TS_L_slope': 20.0 ,
                         'TS_L_intercept': -84.1 ,
                         'length_units': 'cm' }

### Evaluate function
# temp
specimen_df_copy = specimen_dataframe
length_df_copy = length_dataframe

### First step: species_id
species_id = 19350
# Run
specimen_df_copy = specimen_df_copy[ specimen_df_copy.species_id == species_id ]
length_df_copy = length_df_copy[ length_df_copy.species_id == species_id ]
   
### Next step: Expected to combine the two datasets 
spec_df_reframed = (
    specimen_df_copy
    .groupby(['haul_num', 'stratum_num' , 'species_id', 'length'])
    .apply(lambda x: len(x['length']))
    .reset_index(name= 'length_count' )
    )

all_length_df = pd.concat( [ spec_df_reframed , length_df_copy ] , join = 'inner' )

### Pull TS-length parameter values
slope = ts_length_parameters[ 'TS_L_slope' ]
intercept = ts_length_parameters[ 'TS_L_intercept' ]

### Generate TS -- append as column to `all_length_df`
all_length_df[ 'TS' ] = ts_length_regression( all_length_df[ 'length' ] , slope , intercept )

# Convert TS into sigma_bs
all_length_df[ 'sigma_bs' ] = to_linear( all_length_df[ 'TS' ] )

### Mean sigma_bs per haul
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

### Expected output dictionary 
output_dictionary = {
    'length_binned': all_length_df ,
    'haul_mean': mean_haul_sigma_bs ,
    'strata_mean': mean_strata_sigma_bs
}

### Imputation
# temp
strata_options = np.unique( strata_dataframe.stratum_num )

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
   