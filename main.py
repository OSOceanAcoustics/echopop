

## Proposed Structure of the API ##

# This is a class that creates all preliminary files for EchoPro.
# I think this class should be disconnected from EchoPro i.e. can be run independent of EchoPro
# 1. One function could be to mimic the Construct NASC button
# 2. We could also use this to construct the configuration files
from EchoPro import CreateFiles
CreateFiles()

# This is the base class of EchoPro
# 1. Will read the configuration files and store these variables into a dictionary
# 2. Load in biological, acoustic, and stratification files (not associated with
# bootstrapping) into a DataFrame or Xarray objects that can be accessed downstream
from EchoPro import EchoPro
epro_2019 = EchoPro(init_file_path='./config_files/initialization_config.yml',
                    survey_year_file_path='./config_files/survey_year_2019_config.yml',
                    source=3,
                    bio_data_type=1,
                    age_data_status=1)

print(epro_2019.gear_df.head())

# This is a class that runs bootstrapping, the main input will be
# an EchoPro object and bootstrapping/Kriging inputs
# This class will call several other classes:
# 1. A class to Process un-kriged Acoustic & Biological Data
# 2. A class that performs the CV analysis
# 3. A class that performs kriging

# TODO: 1. How do we want to reorganize this api (or do we want to?), now that we know that bootstrapping is not
# TODO: always necessary? If odd or even transects are chosen, then no bootstrapping occurs.

# TODO: Previous idea
# This class will then modify epro_2019 so that it can be used downstream
# epro_2019.run_bootstrapping()

# TODO: new idea
# This class will return important info that can be used downstream (i.e. generate reports and plots)
# epro_2019.run_process(start_transect=1, end_transect=200, kriging=True)

# epro_2019.down_sample()

# TODO: or do we want finer control?
# epro_2019.process_un_kriged()
# epro_2019.run_CV_analysis()
# epro_2019.run_kriging_no_gui()
# epro_2019.run_kriging_gui()
# epro_2019.construct_final_tables()

# This is a class that generates all the reports associated with
# the EchoPro object.
# This class can only be ran after bootstrapping is performed.
# 1. Use epro_2019 to write variables to consolidated xlsx files
# from EchoPro import GenerateReports
# GenerateReports(epro_2019)

# This is a class that creates all visualization plots using
# the EchoPro object.
# This class can only be ran after bootstrapping is performed.
# 1. Use epro_2019 to create plots.
# from EchoPro import Visualize
# Visualize(epro_2019)

