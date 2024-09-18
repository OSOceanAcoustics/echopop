from echopop.survey import Survey

survey = Survey( init_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config.yml" ,
                 survey_year_config_path = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_config.yml" )
survey.input
survey.load_acoustic_data()
survey.load_survey_data()

import pandera as pa
import pandas as pd
from echopop.utils.validate import posfloat, posint, realcircle, realposfloat
import re

df = pd.DataFrame({
    "haul_num": [0, 1, 2, 3, 4, 5, "a", 7, 8],
    "stratum_num": [0, 1, 2, 3, 4, 5, 6, 7, 8],
    "sex": [1, 2, 3, "male", "female", "unsexed", "m", "f", "u"],
    "species_id": [22500, 154.2, 22500, 154.2, 22500, 154.2, "a", "b", "c"],
    "transect_num": [1, 2, 3, 4, 5, 6.1, 7.2, 8.3, 9.4],
    "length_count": 1
})

import copy
from pathlib import Path
from typing import List, Optional, Union

import numpy as np
import pandas as pd
import yaml
from openpyxl import load_workbook

from echopop.core import BIODATA_HAUL_MAP, CONFIG_MAP, DATA_STRUCTURE, LAYER_NAME_MAP
from echopop.utils.data_structure_utils import map_imported_datasets
from echopop.utils.load import read_validated_data, validate_config_structure, validate_data_columns
from echopop.utils.validate import CONFIG_DATA_MODEL, CONFIG_INIT_MODEL

print(pa.infer_schema(df).to_script())

from pandera import DataFrameModel, Field
from pandera.typing import Series
import numpy as np



df = self.input["biology"]["length_df"].copy().reset_index()
cls = LengthBiodata
LengthBiodata.validate_df(df)



df = self.input["biology"]["catch_df"].copy().reset_index()
CatchBiodata.validate_df(df)



df = self.input["biology"]["specimen_df"].copy().reset_index(drop=True)
SpecimenBiodata.validate_df(df)



df = self.input["biology"]["haul_to_transect_df"].copy().reset_index(drop=True)
HaulTransect.validate_df(df)



df = self.input["spatial"]["strata_df"].copy().reset_index(drop=True)
KSStrata.validate_df(df)



df = self.input["spatial"]["geo_strata_df"].copy().reset_index(drop=True)
GeoStrata.validate_df(df)

df = self.input["spatial"]["inpfc_strata_df"].copy().reset_index(drop=True)
GeoStrata.validate_df(df)



df = self.input["acoustics"]["nasc_df"].copy().reset_index(drop=True)
AcousticData.validate_df(df)



df = self.input["statistics"]["kriging"]["isobath_200m_df"].copy().reset_index(drop=True)
IsobathData.validate_df(df)

df = self.input["statistics"]["kriging"]["mesh_df"].copy().reset_index(drop=True)
KrigingMesh.validate_df(df)

class KrigingParameters(BaseDataFrame):
    y_offset: Series[float] = Field(ge=-90.0, le=90.0, nullable=False, alias="dataprep.y_offset")
    corr: Series[float] = Field(ge=0.0, nullable=False, alias="vario.corr")
    hole: Series[float] = Field(ge=0.0, nullable=False, alias="vario.hole")
    lscl: Series[float] = Field(ge=0.0, nullable=False, alias="vario.lscl")
    model: Series[int] = Field(nullable=False, alias="vario.model")
    nugt: Series[float] = Field(ge=0.0, nullable=False, alias="vario.nugt")
    powr: Series[float] = Field(ge=0.0, nullable=False, alias="vario.powr")
    range: Series[float] = Field(ge=0.0, nullable=False, alias="vario.range")
    res: Series[float] = Field(gt=0.0, nullable=False, alias="vario.res")
    sill: Series[float] = Field(ge=0.0, nullable=False, alias="vario.sill")  
    ytox_ratio: Series[float] = Field(nullable=False, alias="vario.ytox_ratio")
    ztox_ratio: Series[float] = Field(nullable=False, alias="vario.ztox_ratio")
    blk_nx: Series[int] = Field(gt=0, nullable=False, alias="krig.blk_nx")
    blk_ny: Series[int] = Field(gt=0, nullable=False, alias="krig.blk_ny")
    blk_nz: Series[int] = Field(gt=0, nullable=False, alias="krig.blk_nz")
    dx0: Series[float] = Field(ge=-180.0, le=180.0, nullable=False, alias="krig.dx0")
    dx: Series[float] = Field(nullable=False, alias="krig.dx")
    dy0: Series[float] = Field(ge=-90.0, le=90.0, nullable=False, alias="krig.dy0")
    dy: Series[float] = Field(nullable=False, alias="krig.dy")
    dz: Series[float] = Field(nullable=False, alias="krig.dz")
    elim: Series[float] = Field(nullable=False, alias="krig.elim")
    eps: Series[float] = Field(nullable=False, alias="krig.eps")
    kmax: Series[int] = Field(gt=0, nullable=False, alias="krig.kmax")
    kmin: Series[int] = Field(gt=0, nullable=False, alias="krig.kmin")
    nx: Series[int] = Field(gt=0, nullable=False, alias="krig.nx")
    ny: Series[int] = Field(gt=0, nullable=False, alias="krig.ny")
    nz: Series[int] = Field(gt=0, nullable=False, alias="krig.nz")
    ratio: Series[float] = Field(nullable=False, alias="krig.ratio")
    srad: Series[float] = Field(gt=0.0, nullable=False, alias="krig.srad")
    x_res: Series[float] = Field(nullable=False, alias="krig.x_res")
    xmin: Series[float] = Field(nullable=False, alias="krig.xmin")
    xmax: Series[float] = Field(nullable=False, alias="krig.xmax")
    xmin0: Series[float] = Field(ge=-180.0, le=180.0, nullable=False, alias="krig.xmin0")
    xmax0: Series[float] = Field(ge=-180.0, le=180.0, nullable=False, alias="krig.xmax0")    
    y_res: Series[float] = Field(nullable=False, alias="krig.y_res")
    ymin: Series[float] = Field(nullable=False, alias="krig.ymin")
    ymax: Series[float] = Field(nullable=False, alias="krig.ymax")
    ymin0: Series[float] = Field(ge=-90.0, le=90.0, nullable=False, alias="krig.ymin0")
    ymax0: Series[float] = Field(ge=-90.0, le=90.0, nullable=False, alias="krig.ymax0")

    # @classmethod
    # def validate_df(cls, df):
    #     # --- Format the column names
    #     df.columns = [col.lower() for col in df.columns]
    #     # ---- Get column types
    #     column_types = cls.get_column_types()

    #     # Apply type coercion based on column types and alias patterns
    #     for column_name, dtype in column_types.items():
    #         pandas_dtype = str(dtype)  # Convert dtype to string directly
            
    #         # Apply coercion based on column patterns
    #         for col in df.columns:
    #             if re.match(column_name, col):
    #                 # Coerce the column to the appropriate dtype
    #                 if df[col].dtype == object and dtype.type.kind in ["i", "u", "f"]:
    #                     df[col] = df[col].astype(float).astype(pandas_dtype)
    #                 else:
    #                     df[col] = df[col].astype(pandas_dtype, errors="ignore")        
    #     # ---- Pass to validation
    #     return cls.judge(df.filter(column_types.keys()))

df = self.input["statistics"]["kriging"]["vario_krig_para_df"].copy()
KrigingParameters.validate_df(df)
