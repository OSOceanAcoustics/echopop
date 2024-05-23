from echopop.survey import Survey

survey = Survey( init_config_path = "./config_files/initialization_config.yml" ,
                 survey_year_config_path = "./config_files/survey_year_2019_config.yml" )
survey.transect_analysis( )
survey.stratified_analysis( )
survey.kriging_analysis( )
self = survey

import copy 
import pandas as pd
import numpy as np

import echopop.utils.message as em

stratified_results_dict = self.results[ 'stratified' ][ 'transect' ]
settings_dict = self.analysis[ 'settings' ][ 'stratified' ]

kriging_results_dict = self.results[ 'kriging' ]
settings_dict = self.analysis[ 'settings' ][ 'kriging' ]

input = 'stratified:transect'
":" in input
input_1 , input_2 = input.split( ":" )
f"{input}_results_msg(self.results['{input_1}']['{input_2}'], self.analysis['settings']['{input_1}'])"
string = f"{input_1}_results_msg(self.results['{input_1}']['{input_2}'], self.analysis['settings']['{input_1}'])"
eval(string)
self.results[ input ].keys()

self.analysis['settings']['stratified']
transect_results_dict = self.results['transect']
settings_dict = self.analysis['settings']['transect']
transect_results_dict.keys()

example_string = "print('Hello, user!')"
eval(example_string)