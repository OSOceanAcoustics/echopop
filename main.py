# import sys
# sys.path.insert(1, './parameters/')
from EchoPro import EchoPro

epro_2019 = EchoPro(init_file_path='./config_files/initialization_config.yml',
                    survey_year_file_path='./config_files/survey_year_2019_config.yml')

# print(epro_2019.get_initialized_parameters())
# print(epro_2019.get_survey_year_parameters())

epro_2019.process()

