# from echopop import Survey
# from pathlib import Path
# from typing import Any, Dict, List, Literal, Optional, Union
# import re
# from collections import OrderedDict
# from typing import Any, Dict, List, Literal, Optional, Union
# import numpy as np
# import pytest
# from pydantic import RootModel, ValidationError
# from echopop.utils.validate_dict import (
#     CONFIG_DATA_MODEL,
#     CONFIG_INIT_MODEL,
#     BiologicalFiles,
#     FileSettings,
#     Geospatial,
#     HaulTransectMap,
#     InputModel,
#     KrigingFiles,
#     KrigingParameters,
#     NASCExports,
#     PatternParts,
#     StratificationFiles,
#     StratifiedSurveyMeanParameters,
#     TransectRegionMap,
#     TSLRegressionParameters,
#     XLSXFile,
# )
# from echopop.utils.validate import posfloat, posint, realcircle, realposfloat

# init_config = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config.yml"
# file_config = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_config.yml"
# survey = Survey(init_config, file_config)

# survey.load_survey_data()
# survey.load_acoustic_data()
# survey.transect_analysis()
# survey.fit_variogram()
# survey.kriging_analysis(best_fit_variogram=True, variogram_parameters={"n_lags": 30})

# input = {
#     "length": {"filename": "blurgh", "sheetname": "sheet1"},
#     "specimen": {"filename": "blargh", "sheetname": "sheet2"},
#     "catch": {"filename": "blorgh", "sheetname": "sheet3"},
#     "haul_to_transect": {"filename": "blorgh", "sheetname": "sheet4"},
# }
# exception = None
# BiologicalFiles.create(**input)
# param = "haul_to_transect"

# input = {
#     "save_file_template": "blurgh_{YEAR}_{COUNTRY}.xlsx",
#     "country_code": ["A", "B"],
#     "file_settings": {
#         "A": {"directory": "blurgh/blargh", "sheetname": "sheet1"},
#         "B": {"directory": "blergh/blorgh", "sheetname": "sheet1"},
#     },
#     "excess": "erroneous",
# }

# HaulTransectMap.create(**input)
# InputModel.judge(invalid_field="invalid_value")
# InputModel.judge(2)

# InputModel.judge()

# print(CONFIG_DATA_MODEL_fields[param]["annotation"].__module__)
# print(CONFIG_DATA_MODEL.model_fields[param]._attributes_set["annotation"].__module__)

# KrigingParameterInputs.validate_realposfloat(None)
