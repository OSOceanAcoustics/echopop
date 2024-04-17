from echopop.survey import Survey
survey = Survey( "./config_files/initialization_config.yml" ,
                 "./config_files/survey_year_2019_config.yml" )
survey.transect_analysis( )
survey.stratified_summary( )