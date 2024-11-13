from openpyxl import Workbook, load_workbook
from openpyxl.styles import Font, Alignment, Border, Side
import numpy as np
import pandas as pd
from pathlib import Path 

# SIM
LENGTH_BINS = np.linspace(2, 80, 40)
AGE_BINS = np.linspace(1, 22, 22)
SEX_OPTS = ["male", "female", "unsexed"]

def cartesian_coord(*arrays):
    grid = np.meshgrid(*arrays)
    coord_list = [entry.ravel() for entry in grid]
    points = np.vstack(coord_list).T
    return points

data = pd.DataFrame(cartesian_coord(LENGTH_BINS, AGE_BINS, SEX_OPTS), 
                    columns=["length", "age", "sex"])
data["biomass"] = 0.0
data.loc[(data.sex == "male") &
         (
             ((data.age == "1.0") & (data.length == "24.0"))
             | ((data.age == "2.0") & (data.length.isin(["26.0", "28.0"])))
             | ((data.age == "3.0") & (data.length.isin(["28.0", "30.0", "32.0"])))
             | ((data.age == "3.0") & (data.length.isin(["30.0", "32.0", "34.0", "36.0"])))
         ),
         "biomass"
    ] = [10, 20, 20, 40, 50, 40, 30, 10]
data.loc[(data.sex == "female") &
         (
             ((data.age == "1.0") & (data.length.isin(["22.0", "24.0"]))
             | ((data.age == "2.0") & (data.length.isin(["24.0", "26.0"])))
             | ((data.age == "3.0") & (data.length.isin(["26.0", "30.0", "32.0", "34.0"])))
             | ((data.age == "3.0") & (data.length.isin(["34.0", "36.0", "38.0", "40.0"])))
             | ((data.age == "4.0") & (data.length.isin(["40.0", "42.0"]))))
         ),
         "biomass"
    ] = [10, 10, 20, 20, 40, 50, 60, 70, 60, 60, 50, 30, 10]
data.loc[(data.sex == "unsexed") &
         (
             ((data.age == "1.0") & (data.length.isin(["22.0", "24.0"]))
             | ((data.age == "2.0") & (data.length.isin(["26.0"]))))
         ),
         "biomass"
    ] = [10, 20, 10]
data_all = (
    data
    .groupby(["age", "length"])
    .sum()
    .reset_index()
    .assign(sex="all")
)

data = pd.concat([data, data_all], ignore_index=True)
data = data.astype({"age": float, "length": float, "biomass": float, "sex": str})

###################
FILENAME = "C:/Users/Brandyn Lucca/Documents/custom_format.xlsx"
sexes = ["all", "male", "female"]
###################
# Create a workbook and select the active sheet
if Path(FILENAME).exists():
    wb = load_workbook(FILENAME)
else:
    wb = Workbook()
    wb.remove(wb["Sheet"])

for sex in sexes:
    data_pvt = data.loc[data.sex == sex].pivot_table(
        columns=["age"],
        index=["length"], 
        values="biomass", 
        margins_name="Subtotal",
        margins=True, 
        aggfunc="sum"
    )
    # data_pvt

    # 
    # data.loc[data.sex == "male", "biomass"].sum() + data.loc[data.sex == "female", "biomass"].sum()
    #
    data_proc = data_pvt.rename_axis(columns=None, index=None)
    # data_proc
    DTSHAPE = data_pvt.shape

    # ----
    if sex.capitalize() in wb.sheetnames: 
        wb.remove(wb[sex.capitalize()])
    # ----
    _ = wb.create_sheet(sex.capitalize())
    # ---- Activate sheet
    ws = wb[sex.capitalize()]

    # Write data with custom headers and gaps
    # ---- Get midway column
    midway_col_idx = int(np.ceil((DTSHAPE[1] - 1) / 2))

    # Create headers
    headers = ["Length (cm)"] + [""] * (midway_col_idx - 1) + ["Age"]
    # ---- Complete the headers
    headers = headers + [""] * (DTSHAPE[1] - len(headers))

    # Append to the worksheet
    ws.append(headers)

    # Add actual column headers
    ws.append([""] + data_proc.columns.values.tolist())

    # Populate the file with all rows except for the last margins row
    # ---- Iterate along the length bins
    for row in LENGTH_BINS:
        # ---- Get the row 
        row_data = data_proc.loc[row, :]
        # ---- Concatenate row name with the rest of the data
        row_data = [row_data.name] + list(row_data)
        # ---- Append
        ws.append(row_data)

    # Get the final row
    ws.append(
        ["Subtotal"] + data_proc.loc["Subtotal"].values[:-1].tolist()    
    )

    # Get the length sums
    sum_length = data_proc.iloc[:-1, -1]

    # Get the aged sums
    sum_age = data_proc.loc["Subtotal"].values[:-1]

    # Total number across all fish
    total_aged = sum_age.sum()

    # Calculate the age-1 proportion
    # ---- Proportion
    age1_proportion = sum_age[0] / total_aged
    # ---- Compute age-2+ values
    age2_aged = total_aged*(1-age1_proportion)

    # Add next row 
    # ---- Age 1+
    age_1_row = [
        "Total (age1+)",
        total_aged
    ]
    # ---- Add the "Over Age" value
    age_1_row = age_1_row + [
        "",
        "Over Age:",
        total_aged,
        ""
    ]
    # ---- Add sexed 
    if sex == "all":
        age_1_row = age_1_row + [
            "Male+Female:",
            data.loc[data.sex.isin(["male", "female"]), "biomass"].sum()
        ]    
    # ---- Append
    ws.append(age_1_row)

    # Add next row
    # ---- Age-2+
    age_all_row = [
        "Total (age2+)",
        age2_aged
    ]
    # ---- Add sexed 
    if sex == "all":
        age_all_row = (
            age_all_row 
            + [""] * 4
            + [
                "Male+Female:",
                data.loc[data.sex.isin(["male", "female"]), "biomass"].sum()
            ]
        )
    # ---- Append
    ws.append(age_all_row)

    # Add empty row
    ws.append([""])

    # Add sheet label
    sheet_label = (
        f"Un-kriged Acoustically Weighted Biomass (mmt) ({sex.capitalize()})"
    )
    # X 1e-9
    # ---- Append
    ws.append(
        [""] * 9 + [sheet_label]
    )

# Save the workbook
wb.save(FILENAME)
# ---- Close
wb.close()

####################################################################################################
from typing import Any, Dict, Literal, List, Union
from echopop import Survey

data_pvt = data.pivot_table(
    columns=["sex", "age"],
    index=["length"], 
    values="biomass", 
    # margins_name="Subtotal",
    # margins=True, 
    aggfunc="sum"
)

data_pvt.unstack().reset_index(name="biomass")
test_survey = Survey
test_survey.analysis = {}
test_survey.analysis.update({
    "transect": dict(population=dict(tables=dict(aged_abundance_df=data_pvt.copy())))
})
survey = test_survey

test_cls = FEATReports(test_survey, reports)
self = test_cls

class FEATReports: 
       
    def __init__(self, 
                 survey: Survey,
                 reports: List[
                     Literal[
                        "kriged_aged_biomass_mesh",
                        "kriged_biomass_mesh",
                        "kriged_length_age_abundance",
                        "kriged_length_age_biomass",
                        "transect_aged_biomass",
                        "transect_length_age_abundance",
                        "transect_length_age_biomass",    
                        "transect_population_estimates",                         
                     ]
                 ],
                 ):
    
        # Assign table-type(s)
        self.reports = reports
        
        # Assign the data
        self.data = survey
   
    def generate(
        self,
        save_directory: Union[str, Path],
    ):
        
        # Validate the reports
        report_methods = self.get_report_type(
            self.reports
        )
        
        #
        report_files = [report_methods[k](self, save_directory) for k in report_methods.keys()]
        
        # Map methods to the appropriate report-type
        # ---- Join over lines
        report_str = "\n".join(report_files)
        # ---- Print
        print(f"The following report tables were generated: \n '{report_str}'")

    @classmethod
    def get_report_type(
        cls,
        requested_reports = List[str],
    ):
        
        # Find matching methods
        mapped_report_methods = {
            k: cls.__METHOD_REFERENCE__.get(k) for k in requested_reports
            if k in cls.__METHOD_REFERENCE__
        }
        
        # Return
        return mapped_report_methods        

    def transect_length_age_biomass_report(
        self,
        save_directory = Union[str, Path],
    ):
        
        # Get the data used for generating the report
        dataset = self.data.analysis["transect"]["population"]["tables"]
        
        # Create the *.xlsx filepath/filename
        filepath = Path(save_directory) / "transect_length_age_biomass_table.xlsx"
        
        # Return the filepath name
        return filepath.as_posix()
    
    # Reference dictionary: report methods
    __METHOD_REFERENCE__ = {
        "transect_length_age_biomass": transect_length_age_biomass_report,
    }
    
    # Reference dictionary: report titles
    __REPORT_TITLES__ = {
        "transect_length_age_biomass":(
            "Un-kriged Acoustically Weighted Biomass (mmt) ({SEX})"
        )
    }
        