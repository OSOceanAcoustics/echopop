from echopop.survey import Survey

survey = Survey( init_config_path = "./config_files/initialization_config.yml" ,
                 survey_year_config_path = "./config_files/survey_year_2019_config.yml" )

survey.transect_analysis( )

survey.stratified_analysis( )

survey.kriging_analysis( )
