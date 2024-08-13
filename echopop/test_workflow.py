from echopop.live.live_survey import LiveSurvey
from echopop.live.sql_methods import SQL
####################################################################################################
# TEST: Set up `LiveSurvey` object
# NOTE: General initialization parameter configuration
live_init_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/live_initialization_config.yml"
# NOTE: File configuration
live_file_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/live_survey_year_2019_config.yml"
# NOTE: Create object
realtime_survey = LiveSurvey(live_file_config_path, live_init_config_path, verbose=True)
# NOTE: String-representation via `LiveSurvey.__repr__`: 
# NOTE: Lists current files being processed and linked databases (WIP)
realtime_survey
####################################################################################################
# TEST: TRIGGER --> NEW ACOUSTIC DATA
# NOTE: Load new acoustic data (Either glob file search or `input_filenames Optional[List[str]]`)
realtime_survey.load_acoustic_data()
# NOTE: Process new acoustic data
# NOTE: This will update linked database tables
realtime_survey.process_acoustic_data()
# NOTE: Generate population estimates (or pass if there are no biological data)
# NOTE: `working_dataset = Literal["acoustic", "biology"]`
realtime_survey.estimate_population(working_dataset="acoustic")
# NOTE: String-representation via `LiveSurvey.__repr__`: 
# NOTE: Lists current files being processed and linked databases (WIP)
realtime_survey
####################################################################################################
# TEST: TRIGGER --> NEW BIOLOGY DATA
# NOTE: Load new biological data (Either glob file search or `input_filenames Optional[List[str]]`)
realtime_survey.load_biology_data()
# NOTE: Process new biological data
# NOTE: This will update linked database tables
realtime_survey.process_biology_data()
# NOTE: Generate population estimates (or pass if there are no acoustic data)
# NOTE: `working_dataset = Literal["acoustic", "biology"]`
realtime_survey.estimate_population(working_dataset="biology")
# NOTE: String-representation via `LiveSurvey.__repr__`: 
# NOTE: Lists current files being processed and linked databases (WIP)
realtime_survey
####################################################################################################
# TEST: `LiveSurvey` --[`files_processed`]--> `Echodataflow`
# NOTE: `LiveSurvey.meta` attribute
# ---- ACOUSTIC
realtime_survey.meta["provenance"]["acoustic_files"]
# ---- BIOLOGICAL
realtime_survey.meta["provenance"]["biology_files"]
# NOTE: SQL function query from database file [cumulative list]
# ---- ACOUSTIC
SQL(db_file=realtime_survey.config["database"]["acoustics"],
    command="select", table_name="files_processed")
# ---- BIOLOGICAL
SQL(db_file=realtime_survey.config["database"]["biology"],
    command="select", table_name="files_processed")
####################################################################################################
# TEST: `LiveSurvey` --[(key) SQL tables]--> Users
# !!! The SQL functions will fail if the tables have not yet been created/initialized
# ---- ACOUSTICS
# NOTE: Mean linear backscatter coefficient (`sigma_bs`) keyed for each haul and stratum
SQL(realtime_survey.config["database"]["acoustics"], "select", table_name="sigma_bs_mean_df")
# NOTE: Along-track acoustically-derived number/biomass densities and NASC 
SQL(realtime_survey.config["database"]["acoustics"], "select", table_name="survey_data_df")
# ---- BIOLOGICAL
# NOTE: Fitted (discretized) length-weight relationship
SQL(realtime_survey.config["database"]["biology"], "select", table_name="weight_fitted_df")
# NOTE: Quantized length-binned weights (summed)
SQL(realtime_survey.config["database"]["biology"], "select", table_name="length_weight_df")
# NOTE: Average weights per stratum
SQL(realtime_survey.config["database"]["biology"], "select", table_name="weight_stratum_df")

dat = realtime_survey.input["acoustics"]["prc_nasc_df"].copy()
dat = dat[dat.latitude > 40]
dat = dat[dat.depth > 20]

import matplotlib.pyplot as plt
import seaborn as sns
from geopy.distance import geodesic
import pandas as pd
from datetime import datetime
import matplotlib.dates as mdates
def calculate_distances(df):
    distances = [0]  # Start with 0 for the first point
    for i in range(1, len(df)):
        point1 = (df.iloc[i - 1]['latitude'], df.iloc[i - 1]['longitude'])
        point2 = (df.iloc[i]['latitude'], df.iloc[i]['longitude'])
        distances.append(geodesic(point1, point2).meters)
    return distances

def parse_datetime(date_str):
    # List of possible formats
    formats = [
        '%Y-%m-%d %H:%M:%S.%f',  # With fractional seconds
        '%Y-%m-%d %H:%M:%S',     # Without fractional seconds
        '%Y-%m-%dT%H:%M:%S.%f',  # ISO 8601 format with fractional seconds
        '%Y-%m-%dT%H:%M:%S'      # ISO 8601 format without fractional seconds
    ]
    
    for fmt in formats:
        try:
            return pd.to_datetime(date_str, format=fmt)
        except (ValueError, TypeError):
            continue  # Try the next format
    
    return pd.NaT  # Return NaT if no formats match

dat["ping_time"] = dat["ping_time"].apply(parse_datetime)

pivot_table = dat.pivot_table(index=["depth"], columns=["ping_time"], values=["NASC"], aggfunc="mean")
# Get the unique distance and depth values for plotting
plt.figure(figsize=(10, 8))
ax = sns.heatmap(pivot_table, cmap="viridis", cbar_kws={'label': 'NASC'})
plt.gca().xaxis.set_major_locator(mdates.MinuteLocator(interval=30))  # Major ticks every 30 minutes
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))  # Format as hour:minute
plt.gcf().autofmt_xdate()
ax.set_xticks(ax.get_xticks()[::max(len(ax.get_xticks()) // 10, 1)])  # Show fewer ticks if necessary
plt.xlabel('Ping time')
plt.ylabel('Depth')
# plt.gca().invert_yaxis()  # To have depth increasing downwards like in a typical depth plot
plt.show()


dat.groupby(["ping_time"]).size()
unique_pairs = dat.drop_duplicates(subset=['latitude', 'longitude']).sort_values("ping_time")

unique_pairs["d"] = calculate_distances(dat)
df['cumulative_distance'] = df['distance'].cumsum()

unique_distances = dat.groupby('source')[['latitude', 'longitude']].unique().reset_index()
unique_distances = unique_distances.explode('distance')


# Create a pivot table to reshape the dataframe suitable for a heatmap
dat['source_id'] = dat['source'].astype('category').cat.codes
pivot_table = dat.pivot(index=["depth"], columns=["distance"], values=["NASC"])
dat.groupby('source')['distance'].cumsum()
plt.plot(index="depth", columns="distance", values="NASC")
plt.show()

data = {
    'distance': [1, 1, 2, 2, 1, 1, 3, 3],
    'depth': [1, 2, 1, 2, 1, 2, 1, 2],
    'source': ['A', 'A', 'A', 'A', 'B', 'B', 'B', 'B']
}
dat = pd.DataFrame(data)
dat = dat.sort_values(by=['source', 'distance'])
unique_distances = dat.groupby('source')['distance'].unique().reset_index()
unique_distances = unique_distances.explode('distance')

unique_distances['distance'] = pd.to_numeric(unique_distances['distance'], errors='coerce')
unique_distances['distance_diff'] = unique_distances.groupby('source')['distance'].diff().fillna(0)
unique_distances['cumsum_diff'] = unique_distances.groupby('source')['distance_diff'].cumsum()
unique_distances['Cumsum_dist'] = unique_distances['cumsum_diff'].cumsum()
unique_distances['Cumsum_dist'] = pd.to_numeric(unique_distances['Cumsum_dist'], errors='coerce')
dat = dat.merge(unique_distances[['source', 'distance', 'Cumsum_dist']], on=['source', 'distance'], how='left')


# Calculate cumulative sum of distances for each source
dat['Cumsum_dist'] = dat.groupby('source')['distance'].transform(lambda x: x.cumsum())
dat['Cumsum_dist_within_source'] = dat.groupby('source')['distance'].cumsum()

dat['Cumsum_dist'] = dat.groupby('source')['Cumsum_dist_within_source'].transform(lambda x: x + x.shift(1).fillna(0).cumsum())
