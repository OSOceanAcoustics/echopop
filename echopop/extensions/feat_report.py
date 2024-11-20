from pathlib import Path
from typing import Any, Callable, Dict, List, Literal, Optional, Tuple, Union

import numpy as np
import pandas as pd
import pandas.io.formats.excel as pdif
from openpyxl import Workbook, load_workbook
from openpyxl.worksheet.worksheet import Worksheet

from .. import Survey
from ..biology import filter_species

####################################################################################################
# FILE WRITING UTILITY
# --------------------------------------------------------------------------------------------------


def initialize_workbook(filepath: Path) -> Workbook:
    """
    Initialize Excel workbook
    """

    # Check if workbook already exists
    if filepath.exists():
        wb = load_workbook(filepath)
    # ---- If not, generate it and remove the default "Sheet" that is created
    else:
        wb = Workbook()
        wb.remove(wb["Sheet"])

    # Return
    return wb


def format_file_sheet(
    sheetname: str,
    workbook: Workbook,
) -> Worksheet:
    """
    Initialize and format Excel worksheet(s)
    """

    # Check whether the sheetname already exists, and remove if needed
    if sheetname in workbook.sheetnames:
        workbook.remove(workbook[sheetname])

    # Create the new sheet
    _ = workbook.create_sheet(sheetname)

    # Activate and return the new sheet
    return workbook[sheetname]


def format_table_headers(
    dataframe: pd.DataFrame, col_name: str = "Age", row_name: str = "Length (cm)"
) -> List[str]:
    """
    Format worksheet headers independent of the input data
    """

    # Get the shape of the table
    DTSHAPE = dataframe.shape

    # Write data with custom headers and gaps
    # ---- Get midway column
    midway_col_idx = int(np.ceil((DTSHAPE[1] - 1) / 2))

    # Create headers
    headers = [row_name] + [""] * (midway_col_idx - 1) + [col_name]
    # ---- Complete the headers
    return headers + [""] * (DTSHAPE[1] - len(headers))


def append_datatable_rows(worksheet: Worksheet, dataframe: pd.DataFrame) -> None:
    """
    Write data rows into Excel sheet
    """

    # Add column names (i.e. age)
    worksheet.append([""] + dataframe.columns.values.tolist())

    # Get the dataset row indices
    ROW_INDICES = dataframe.index.values

    # Populate the file with all rows except for the last margins row
    # ---- Iterate along the length bins
    for row in ROW_INDICES[:-1]:
        # ---- Get the row
        row_data = dataframe.loc[row, :]
        # ---- Concatenate row name with the rest of the data
        row_data = [row_data.name] + list(row_data)
        # ---- Append
        worksheet.append(row_data)

    # Append the final row
    # ---- Only apply if 'Subtotal' margin included
    if "Subtotal" in ROW_INDICES:
        worksheet.append(["Subtotal"] + dataframe.loc["Subtotal"].values[:-1].tolist())


def append_table_aggregates(
    worksheet: Worksheet,
    dataframe: pd.DataFrame,
    sex: Literal["all", "female", "male"],
    tables_dict: Dict[str, pd.DataFrame],
) -> None:
    """
    Add table aggregate calculations to bottom-left margin of worksheet
    """

    # Get the aged sums
    if "Unaged" in dataframe.columns:
        sum_age = dataframe.loc["Subtotal"].values[:-2]
        # ---- Unaged sum
        sum_unaged = dataframe.loc["Subtotal", "Unaged"]
    else:
        sum_age = dataframe.loc["Subtotal"].values[:-1]
        # ---- Unaged sum
        sum_unaged = 0.0

    # Total number across all fish
    total_aged = sum_age.sum() + sum_unaged

    # Calculate the age-1 proportion
    # ---- Proportion
    age1_proportion = dataframe.loc["Subtotal", 1] / sum_age.sum()
    # ---- Compute age-2+ values
    age2_aged = total_aged * (1 - age1_proportion)

    # Add next row
    # ---- Age 1+
    age_1_row = ["Total (age1+)", total_aged]
    # ---- Add the "Over Age" value
    age_1_row = age_1_row + ["", "Over Age:", sum_age.sum(), ""]
    # ---- Add sexed
    if sex == "all":
        age_1_row = age_1_row + [
            "Male+Female:",
            tables_dict["female"].iloc[:-1, :-1].sum().sum()
            + tables_dict["male"].iloc[:-1, :-1].sum().sum(),
        ]
    # ---- Append
    worksheet.append(age_1_row)

    # Add next row
    # ---- Age-2+
    age_all_row = ["Total (age2+)", age2_aged]
    # ---- Add sexed

    if sex == "all":
        age_all_row = (
            age_all_row
            + [""] * 4
            + [
                "Male+Female:",
                tables_dict["female"].iloc[:-1, 1:-1].sum().sum()
                + tables_dict["male"].iloc[:-1, 1:-1].sum().sum(),
            ]
        )
    # ---- Append
    worksheet.append(age_all_row)


def append_sheet_label(
    worksheet: Worksheet,
    title: str,
    sex: str,
) -> None:
    """
    Append sheet title to tail of Excel worksheet
    """

    # Add empty row
    worksheet.append([""])

    # Add sheet label
    sheet_label = title.replace("{SEX}", sex.capitalize())
    # ---- Append
    worksheet.append([""] * 9 + [sheet_label])


####################################################################################################
# DATA WRANGLING
# --------------------------------------------------------------------------------------------------


def pivot_dataframe_data(
    transect_df: pd.DataFrame,
    tables_dict: Dict[str, pd.DataFrame],
    sex: str,
    stratum_name: str,
    contrasts: List[str] = ["transect_num", "longitude", "latitude"],
):
    """
    Pivot DataFrame object for row and column indexing
    """

    # Filter the correct sex for the aged weight proportions
    aged_weight_proportions_df = tables_dict[sex]

    # Assign the correct biomass column name (based on sex)
    if sex != "all":
        biomass_name = f"biomass_{sex}" if f"biomass_{sex}" in transect_df.columns else "biomass"
    else:
        biomass_name = "biomass"

    # Filter the transect dataframe
    transect_data = transect_df.filter(contrasts + [stratum_name, biomass_name])

    # Set the index to the strata
    transect_data.set_index([stratum_name], inplace=True)

    # Merge the transect and weight proportions dataframes
    aged_transect_data = pd.merge(
        transect_data, aged_weight_proportions_df, left_index=True, right_index=True
    )

    # Propagate the biomass column over the age column proportions
    aged_transect_data.loc[:, 1.0:] = aged_transect_data.loc[:, 1.0:].mul(
        aged_transect_data[biomass_name], axis=0
    )

    # Sort the final dataframe and return
    return aged_transect_data.sort_values(contrasts).reset_index()


def repivot_table(
    age_length_dataframe: pd.DataFrame,
    variable: str,
    length_dataframe: Optional[pd.DataFrame] = None,
):
    """
    Repivot existing pivot table(s)
    """

    # Restack the data table
    data_stk = age_length_dataframe.stack(future_stack=True).sum(axis=1).reset_index(name=variable)

    # Concatenate unaged data, if needed
    if length_dataframe is not None:
        # ---- Restack
        unaged_data_stk = length_dataframe.sum(axis=1).reset_index(name=variable)
        # ---- Get the right-most age-bin and add 1
        shifted_bin = data_stk["age_bin"].max() + 1
        # ---- Get the mid-point for later renaming
        shifted_bin_mid = shifted_bin.mid
        # ---- Add age column
        unaged_data_stk["age_bin"] = shifted_bin
        # ---- Concatenate with the aged data
        data_stk = pd.concat([data_stk, unaged_data_stk], ignore_index=True)
    else:
        shifted_bin_mid = None

    # Convert the length and age columns from intervals into numerics
    # ---- Length
    data_stk["length"] = data_stk["length_bin"].apply(lambda x: x.mid).astype(float)
    # ---- Age
    data_stk["age"] = data_stk["age_bin"].apply(lambda x: x.mid).astype(float)

    # Pivot
    data_pvt = data_stk.pivot_table(
        columns=["age"],
        index=["length"],
        values=variable,
        margins_name="Subtotal",
        margins=True,
        aggfunc="sum",
    )

    # Drop the extra axes
    data_proc = data_pvt.rename_axis(columns=None, index=None)

    # Rename the unaged placeholder to "Unaged" and return
    return data_proc.rename(columns={shifted_bin_mid: "Unaged"})


def pivot_aged_weight_proportions(
    weight_distribution: pd.DataFrame,
    age1_exclusion: bool,
    stratum_name: str,
) -> Dict[str, pd.DataFrame]:
    """
    Pivot distributions of aged weight proportions
    """

    # Stack the age distributions
    weight_stk = weight_distribution.sum().reset_index(name="weight")

    # Concatenate a copy representing "all" fish
    weight_stk_all = (
        weight_stk.groupby([stratum_name, "age_bin"], observed=False)["weight"].sum().reset_index()
    )
    # ---- Add to the original stacked dataframe
    weight_stk_full = pd.concat([weight_stk, weight_stk_all.assign(sex="all")], ignore_index=True)

    # Convert age bins into number ages
    weight_stk_full["age"] = weight_stk_full["age_bin"].apply(lambda x: x.mid).astype(float)

    # Pivot the table
    weight_pvt = weight_stk_full.pivot_table(
        columns=["sex", stratum_name], index=["age"], values="weight"
    )

    # Zero out the age-1 category, if they are being excluded
    if age1_exclusion:
        weight_pvt.loc[1] = 0.0

    # Convert to proportions
    weight_prop_pvt = (weight_pvt / weight_pvt.sum()).fillna(0.0)

    # Create table dictionary
    tables = {sex: weight_prop_pvt.loc[:, sex].T for sex in ["all", "male", "female"]}

    # Return the output
    return tables


def pivot_haul_tables(
    dataframe: pd.DataFrame,
    sex: str,
) -> pd.DataFrame:
    """
    Pivot haul count datatable(s)
    """

    # Subset the data into the specific sex
    haul_sex = dataframe.loc[sex, :]

    # Stack the table
    haul_stk = haul_sex.stack(future_stack=True).reset_index(name="count")

    # Convert the length bins into number lengths
    haul_stk.loc[:, "length"] = haul_stk["length_bin"].apply(lambda x: x.mid).astype(float)

    # Repivot the table with margin subtotals
    haul_agg = haul_stk.pivot_table(
        index=["length"],
        columns=["haul_num"],
        values="count",
        margins=True,
        margins_name="Subtotal",
        aggfunc="sum",
    )

    # Return the new pivot table with column and index names removed
    return haul_agg.rename_axis(columns=None, index=None)


####################################################################################################
# EXCEL REPORT WRITING
# --------------------------------------------------------------------------------------------------


def write_aged_dataframe_report(
    tables_dict: Dict[str, pd.DataFrame],
    filepath: Path,
):
    """
    Write georeferenced population estimates apportioned across age-classes
    """

    # Create the workbook, if needed
    wb = initialize_workbook(filepath)

    # Initialize
    for sex in ["all", "female", "male"]:

        # Add/overwrite the sheet
        ws = format_file_sheet(sex.capitalize(), wb)

        # Set up superheader representing 'sex'
        if sex != "all":
            ws.append([sex.capitalize()])
        else:
            ws.append([f"{sex.capitalize()} (Male + Female)"])

    # Save the workbook
    wb.save(filepath)
    # ---- Close
    wb.close()

    # Remove header formatting
    pdif.ExcelFormatter.header_style = None

    # Open up a writer engine
    with pd.ExcelWriter(filepath, engine="openpyxl", mode="a", if_sheet_exists="overlay") as writer:

        # Iterate through all sexes
        for sex in ["all", "female", "male"]:

            # Subset the data dictionary for the particular sheet
            sheet_data = tables_dict[sex]

            # Add actual column names (i.e. age) and datatable rows
            sheet_data.to_excel(writer, sheet_name=sex.capitalize(), startrow=1, index=None)


def write_age_length_table_report(
    tables_dict: Dict[str, pd.DataFrame],
    table_title: str,
    filepath: Path,
) -> None:
    """
    Write pivot table comprising population estimates apportioned across age and length bins
    """

    # Create the workbook, if needed
    wb = initialize_workbook(filepath)

    # Iterate through all sexes
    for sex in ["all", "female", "male"]:

        # Subset the data dictionary for the particular sheet
        sheet_data = tables_dict[sex]

        # Add/overwrite the sheet
        ws = format_file_sheet(sex.capitalize(), wb)

        # Set up the headers for the file sheet
        ws.append(format_table_headers(sheet_data))

        # Add actual column names (i.e. age) and datatable rows
        append_datatable_rows(ws, sheet_data)

        # Append the data table aggregates
        append_table_aggregates(ws, sheet_data, sex, tables_dict)

        # Append the trailing sheet label along the bottom-most row
        append_sheet_label(ws, table_title, sex)

    # Save the workbook
    wb.save(filepath)
    # ---- Close
    wb.close()


def write_haul_report(
    tables_dict: Dict[str, pd.DataFrame],
    table_title: str,
    filepath: Path,
):
    """
    Write haul count table(s)
    """

    # Create the workbook, if needed
    wb = initialize_workbook(filepath)

    # Iterate through all sexes
    for sex in ["all", "female", "male"]:

        # Subset the data dictionary for the particular sheet
        sheet_data = tables_dict[sex]

        # Add/overwrite the sheet
        ws = format_file_sheet(sex.capitalize(), wb)

        # Set up the headers for the file sheet
        ws.append(format_table_headers(sheet_data, col_name=f"Haul Number ({sex.capitalize()})"))

        # Add actual column names (i.e. age) and datatable rows
        append_datatable_rows(ws, sheet_data)

        # Add grand total
        total_row = ["Total", sheet_data.iloc[:-1, :-1].sum().sum()]
        # ---- Extend, if needed
        if sex == "all":
            total_row = total_row + [
                "",
                "Male+Female",
                tables_dict["female"].iloc[:-1, :-1].sum().sum()
                + tables_dict["male"].iloc[:-1, :-1].sum().sum(),
            ]
        # ---- Append
        ws.append(total_row)

        # Add two new rows
        ws.append([""])

        # Append the trailing sheet label along the bottom-most row
        append_sheet_label(ws, table_title, sex)

    # Save the workbook
    wb.save(filepath)
    # ---- Close
    wb.close()


####################################################################################################
# FEATReports class
# --------------------------------------------------------------------------------------------------


class FEATReports:
    """
    Echopop FEAT report generation function

    This class includes methods processing and writing reports used by the FEAT team

    Parameters
    ----------
    survey: echopop.survey.Survey
        An `echopop.survey.Survey` class object containing all of the relevant results required for
        producing various `*.xlsx`-formatted reports

    reports: Optional[List[Literal["aged_length_haul_counts", "kriged_aged_biomass_mesh", \
        "kriged_biomass_mesh", "kriged_length_age_biomass", "kriging_input", \
            "total_length_haul_counts", "transect_aged_biomass", "transect_length_age_abundance", \
                "transect_length_age_biomass", "transect_population_results"]]]]]
        A list of compatible reports that can be generated. When no list is supplied, then all
        reports are generated. Available options include:

            - *"aged_length_haul_counts"* \n
            Table comprising distributions of counts from aged fish for each haul

            - *"kriged_biomass_mesh"* \n
            Dataframe comprising georeferenced kriged mesh results with biomass estimates

            - *"kriged_aged_biomass_mesh"* \n
            Dataframe comprising georeferenced kriged mesh results with biomass distributed across
            all age bins

            - *"kriging_input"* \n
            Dataframe comprising georeferenced abundance, biomass, and NASC values that can be used
            for the kriging analysis. Note that only the 'biomass' column is currently used within
            `Echopop`

            - *"kriged_length_age_biomass"* \n
            Table comprising kriged biomass estimates distributed over age and length bins

            - *"transect_aged_biomass"* \n
            Dataframe comprising georeferenced along-transect biomass estimates with biomass
            distributed across age bins

            - *"transect_length_age_abundance"* \n
            Table comprising along-transect abundance estimates distributed over age and length bins

            - *"transect_length_age_biomass"* \n
            Table comprising along-transect biomass estimates distributed over age and length bins

            - *"transect_population_estimates"* \n
            Dataframe comprising georeferenced along-transect population estimates with that
            includes abundance, biomass, biomass density, number density, and values specific to
            each sex

    save_directory: Optional[Union[str, pathlib.Path]]
        A string or Path object specifying the directory where generated reports will be saved.
        When no directory path is supplied for this argument, the `FEATReports` class will produce
        reports in the folder called `reports` contained within `config["data_root_dir"]` stored
        in the input `echopop.survey.Survey` object

    Attributes
    ----------
    data: echopop.survey.Survey
        The `echopop.survey.Survey`-class object containing all of the data used to generate the
        report `*.xlsx` files

    reports: List[str]
        A list of the report-types that will be produced

    save_directory: pathlib.Path
        A Path object comprising the save directory path where report files will be written

    Methods
    -------
    *generate()*
        Report factory method that generates the `*.xlsx` files used by the FEAT team

    Notes
    -----
    The entire workflow including results from `Survey.transect_analysis()` and
    `Survey.kriging_analysis()` are required to produce all of the report files generated by the
    `FEATreports` object

    Examples
    --------
    >>> from echopop.survey import Survey
    >>> survey = Survey(init_config_path, survey_year_config_path)
    >>> survey.load_survey_data()
    >>> survey.load_acoustic_data()
    >>> survey.transect_analysis(verbose=False)
    >>> survey.kriging_analysis(variogram_parameters=dict(model=["exponential", "bessel"],
    ... n_lags=30), variable="biomass_density", verbose=False)
    >>> from echopop.extensions.FEATreports import FEATreports
    >>> reports = ["aged_length_haul_counts", "kriged_biomass_mesh"]
    >>> FEATReports(survey, reports, save_directory="C:/Users/Data").generate()
    The following report tables were generated:
       -'C:/Users/Data/aged_length_haul_counts_table.xlsx'
       -'C:/Users/Data/kriged_biomass_mesh_dataframe.xlsx'
       -'C:/Users/Data/kriged_biomass_mesh_reduced_dataframe.xlsx'
    >>> FEATReports(survey, save_directory="C:/Users/Data").generate()
    The following report tables were generated:
       -'C:/Users/Data/aged_length_haul_counts_table.xlsx'
       -'C:/Users/Data/kriged_aged_biomass_mesh_dataframe.xlsx'
       -'C:/Users/Data/kriged_aged_biomass_reduced_mesh_dataframe.xlsx'
       -'C:/Users/Data/kriged_biomass_mesh_dataframe.xlsx'
       -'C:/Users/Data/kriged_biomass_mesh_reduced_dataframe.xlsx'
       -'C:/Users/Data/kriged_length_age_biomass_table.xlsx'
       -'C:/Users/Data/kriging_input_dataframe.xlsx'
       -'C:/Users/Data/total_length_haul_counts_table.xlsx'
       -'C:/Users/Data/transect_aged_biomass_dataframe.xlsx'
       -'C:/Users/Data/transect_aged_biomass_reduced_dataframe.xlsx'
       -'C:/Users/Data/transect_length_age_abundance_table.xlsx'
       -'C:/Users/Data/transect_length_age_biomass_table.xlsx'
       -'C:/Users/Data/transect_population_results_dataframe.xlsx'

    """

    def __init__(
        self,
        survey: Survey,
        reports: Optional[
            List[
                Literal[
                    "aged_length_haul_counts",
                    "kriged_aged_biomass_mesh",
                    "kriged_biomass_mesh",
                    "kriged_length_age_biomass",
                    "kriging_input",
                    "total_length_haul_counts",
                    "transect_aged_biomass",
                    "transect_length_age_abundance",
                    "transect_length_age_biomass",
                    "transect_population_results",
                ]
            ]
        ] = None,
        save_directory: Optional[Union[str, Path]] = None,
    ):

        # Assign table-type(s)
        if reports is None:
            reports = [
                "aged_length_haul_counts",
                "kriged_aged_biomass_mesh",
                "kriged_biomass_mesh",
                "kriged_length_age_biomass",
                "kriging_input",
                "total_length_haul_counts",
                "transect_aged_biomass",
                "transect_length_age_abundance",
                "transect_length_age_biomass",
                "transect_population_results",
            ]
        # ---- Assign to attribute
        self.reports = reports

        # Assign the data
        self.data = survey

        # Save directory
        if not save_directory:
            # ---- Get the save directory configured within the `Survey` configuration
            survey_directory = Path(survey.config["data_root_dir"])
            # ---- Create the save directory
            save_directory = survey_directory / "reports"
            # ---- Create the directory, if needed
            save_directory.mkdir(parents=True, exist_ok=True)
            # ---- Assign the save directory
            self.save_directory = save_directory
        else:
            # ---- Validate existence
            if not Path(save_directory).exists():
                raise FileNotFoundError(
                    f"The reports save directory '{Path(save_directory).as_posix()}' does not "
                    f"exist!"
                )
            # ---- Assign the save directory
            self.save_directory = save_directory

    def generate(self):
        """
        Report factory method
        """

        # Validate the reports
        report_methods = self.get_report_type(self.reports)

        # Validate the reports
        unknown_requests, report_methods = self.get_report_type(self.reports)

        # Print unexpected reports
        if len(unknown_requests) > 0:
            # ---- Format join
            unknown_str = "\n   ".join(f"-'{request}'" for request in unknown_requests)
            # ---- Print
            print(
                f"The following requested reports do not match available report-types (use "
                f"`FEATReports.report_options()` for a complete list/print-out of available "
                f"reports):\n   {unknown_str}"
            )

        # Run all input report generation functions
        report_files = [v.get("function")(self, **v) for _, v in report_methods.items()]

        # Map methods to the appropriate report-type
        # ---- Initialize an empty list
        report_list = sum(
            ([item] if isinstance(item, str) else list(item) for item in report_files), []
        )
        # ---- Join over lines
        report_str = "\n   ".join(f"-'{file}'" for file in report_list)
        # ---- Print
        print(f"The following report tables were generated:\n   {report_str}")

    @classmethod
    def get_report_type(
        cls,
        requested_reports=List[str],
    ) -> Tuple[List[str], Dict[str, Any]]:

        # Identify superfluous requested reports
        unknown_reports = [report for report in requested_reports
                           if report not in cls.__METHOD_REFERENCE__]

        # Find matching methods
        mapped_report_methods = {
            k: cls.__METHOD_REFERENCE__.get(k)
            for k in requested_reports
            if k in cls.__METHOD_REFERENCE__
        }

        # Return
        return unknown_reports, mapped_report_methods

    # REPORT METHODS
    def aged_length_haul_counts_report(
        self,
        data: Callable[[pd.DataFrame], pd.DataFrame],
        title: str,
        filename: Callable[[Union[str, Path]], Path],
        **kwargs,
    ):

        # Get the dataset
        (specimen_data,) = filter_species(
            data(self).input["biology"]["specimen_df"],
            data(self).analysis["settings"]["transect"]["species_id"],
        )

        # Consolidate the dataset to also include "all"
        full_data = pd.concat(
            [
                specimen_data.loc[specimen_data.group_sex == "sexed", :],
                specimen_data.loc[specimen_data.group_sex == "sexed", :].assign(sex="all"),
            ]
        )

        # Remove missing data
        full_data = full_data.loc[~np.isnan(full_data.age)]

        # Initial pivot table
        full_pvt = full_data.pivot_table(
            index=["sex", "length_bin"],
            columns=["haul_num"],
            values="length",
            aggfunc="count",
            observed=False,
        )

        # Create data dictionary containing the sex-specific tables
        haul_pvt_tables = {
            sex: pivot_haul_tables(full_pvt, sex) for sex in ["all", "female", "male"]
        }

        # Write the *.xlsx sheet
        write_haul_report(haul_pvt_tables, title, filename(self))

        # Return the filepath name
        return filename(self).as_posix()

    def kriged_aged_biomass_mesh_report(
        self,
        data: Callable[[pd.DataFrame], pd.DataFrame],
        filename: Callable[[Union[str, Path]], Path],
        **kwargs,
    ):

        # Get the data dictionary
        data_dict = data(self)

        # Get 'exclude_age1' argument value
        age1_exclusion = data_dict.analysis["settings"]["transect"]["exclude_age1"]

        # Get the stratum name
        stratum_name = data_dict.analysis["settings"]["transect"]["stratum_name"]

        # Get the weight distributions
        weight_distribution = data_dict.analysis["transect"]["biology"]["distributions"]["weight"][
            "aged_length_weight_tbl"
        ].copy()

        # Reformat the aged weight proportions tables
        tables = pivot_aged_weight_proportions(weight_distribution, age1_exclusion, stratum_name)

        # Get mesh dataframe
        mesh_df = data_dict.results["kriging"]["mesh_results_df"].copy()
        # ---- Filter
        mesh_df_copy = mesh_df.filter(["latitude", "longitude", stratum_name, "biomass"]).copy()

        # Combine the pivot tables with the transect tdata
        mesh_pvt_tables = {
            sex: pivot_dataframe_data(
                mesh_df_copy, tables, sex, stratum_name, contrasts=["longitude", "latitude"]
            )
            for sex in ["all", "female", "male"]
        }

        # Write the *.xlsx sheet (full dataset)
        write_aged_dataframe_report(mesh_pvt_tables, filename(self)[0])

        # Remove the empty cells
        mesh_pvt_reduced_tables = {
            sex: mesh_pvt_tables[sex].loc[lambda x: x.biomass > 0.0]
            for sex in ["all", "female", "male"]
        }

        # Write the *.xlsx sheet (reduced dataset)
        write_aged_dataframe_report(mesh_pvt_reduced_tables, filename(self)[1])

        # Return the filepath name
        return filename(self)[0].as_posix(), filename(self)[1].as_posix()

    def kriged_biomass_mesh_report(
        self,
        data: Callable[[pd.DataFrame], pd.DataFrame],
        filename: Callable[[Union[str, Path]], Path],
        **kwargs,
    ):

        # Get the stratum name
        stratum_name = self.data.analysis["settings"]["transect"]["stratum_name"]

        # Get apportionment information
        apportion_df = (
            data(self)
            .analysis["transect"]["biology"]["proportions"]["weight"][
                "aged_unaged_sex_weight_proportions_df"
            ]
            .copy()
        )
        # ---- Pivot
        apportion_pvt = apportion_df.pivot_table(index=["stratum_num"], columns=["sex"])
        # ---- Compute the proportions across strata
        stratum_prop = (
            apportion_pvt["weight_proportion_overall_aged"]
            + apportion_pvt["weight_proportion_overall_unaged"]
        )

        # Get mesh results
        mesh_df = data(self).results["kriging"]["mesh_results_df"].copy()
        # ---- Set the index
        mesh_df.set_index(stratum_name, inplace=True)

        # Merge the mesh and apportionment dataframes
        merged_mesh = pd.merge(mesh_df, stratum_prop, left_index=True, right_index=True)
        # ---- Propagate the proportions across biomass estimates
        merged_mesh.loc[:, "female":"male"] = merged_mesh.loc[:, "female":"male"].mul(
            merged_mesh["biomass"], axis=0
        )
        # ---- Rename the columns
        merged_mesh.rename(
            columns={"female": "biomass_female", "male": "biomass_male"}, inplace=True
        )

        # Back-calculate the kriged standard deviation
        merged_mesh.loc[:, "kriged_sd"] = (
            merged_mesh.loc[:, "kriged_mean"] * merged_mesh.loc[:, "sample_cv"]
        )

        # Update and filter the columns that will be exported
        output_df = merged_mesh.reset_index().filter(
            [
                "latitude",
                "longitude",
                stratum_name,
                "biomass",
                "biomass_female",
                "biomass_male",
                "sample_cv",
                "kriged_sd",
            ]
        )

        # Save the *.xlsx sheet (full dataset)
        output_df.to_excel(filename(self)[0], sheet_name="Sheet1", index=None)

        # Reduce the dataframe and then save
        output_df.loc[output_df.biomass > 0.0].to_excel(
            filename(self)[1], sheet_name="Sheet1", index=None
        )

        # Return the filepath name
        return filename(self)[0].as_posix(), filename(self)[1].as_posix()

    def kriged_length_age_biomass_report(
        self,
        data: Callable[[pd.DataFrame], pd.DataFrame],
        title: str,
        filename: Callable[[Union[str, Path]], Path],
        **kwargs,
    ):

        # Initialize a pivot table from the dataset
        dataset = data(self).pivot_table(
            index=["sex", "length_bin"],
            columns=["age_bin"],
            values=["biomass_apportioned"],
            aggfunc="sum",
            observed=False,
        )

        # Repivot the datasets for each sex
        tables = {
            sex: repivot_table(dataset.loc[sex, :] * 1e-9, "biomass")
            for sex in ["all", "male", "female"]
        }

        # Write the *.xlsx sheet
        write_age_length_table_report(tables, title, filename(self))

        # Return the filepath name
        return filename(self).as_posix()

    def kriging_input_report(
        self,
        data: Callable[[pd.DataFrame], pd.DataFrame],
        filename: Callable[[Union[str, Path]], Path],
        **kwargs,
    ):

        # Get the dataset
        dataset = data(self).sort_values(["transect_num", "longitude", "latitude"])

        # Filter the dataset columns
        dataset_filt = dataset.filter(
            ["latitude", "longitude", "biomass_density", "nasc", "number_density"]
        )

        # and save the *.xlsx sheet
        dataset_filt.to_excel(filename(self), sheet_name="Sheet1", index=None)

        # Return the filepath name
        return filename(self).as_posix()

    def total_length_haul_counts_report(
        self,
        data: Callable[[pd.DataFrame], pd.DataFrame],
        title: str,
        filename: Callable[[Union[str, Path]], Path],
        **kwargs,
    ):

        # Get the dataset
        length_data, specimen_data = filter_species(
            [
                data(self).input["biology"]["length_df"],
                data(self).input["biology"]["specimen_df"],
            ],
            data(self).analysis["settings"]["transect"]["species_id"],
        )

        # Consolidate the dataset to also include "all"
        # ---- Length data
        full_length_data = pd.concat(
            [length_data.loc[length_data.group_sex == "sexed", :], length_data.assign(sex="all")]
        )
        # ---- Specimen data
        full_specimen_data = pd.concat(
            [
                specimen_data.loc[specimen_data.group_sex == "sexed", :],
                specimen_data.assign(sex="all"),
            ]
        )

        # Remove missing data
        # ---- Length
        full_length_data = full_length_data.loc[~np.isnan(full_length_data.length)]
        # ---- Specimen
        full_specimen_data = full_specimen_data.loc[
            (~np.isnan(full_specimen_data.length)) & (~np.isnan(full_specimen_data.age))
        ]

        # Initial pivot table
        # ---- Length
        full_length_pvt = full_length_data.pivot_table(
            index=["sex", "length_bin"],
            columns=["haul_num"],
            values="length_count",
            aggfunc="sum",
            observed=False,
        )
        # ---- Specimen
        full_specimen_pvt = full_specimen_data.pivot_table(
            index=["sex", "length_bin"],
            columns=["haul_num"],
            values="length",
            aggfunc="count",
            observed=False,
        )

        # Combine the datasets
        full_pvt = full_length_pvt.add(full_specimen_pvt, fill_value=0).astype(int)

        # Create data dictionary containing the sex-specific tables
        haul_pvt_tables = {
            sex: pivot_haul_tables(full_pvt, sex) for sex in ["all", "female", "male"]
        }

        # Write the *.xlsx sheet
        write_haul_report(haul_pvt_tables, title, filename(self))

        # Return the filepath name
        return filename(self).as_posix()

    def transect_aged_biomass_report(
        self,
        data: Callable[[pd.DataFrame], pd.DataFrame],
        filename: Callable[[Union[str, Path]], Path],
        **kwargs,
    ):

        # Get the data dictionary
        data_dict = data(self)

        # Get 'exclude_age1' argument value
        age1_exclusion = data_dict["settings"]["transect"]["exclude_age1"]

        # Get the stratum name
        stratum_name = data_dict["settings"]["transect"]["stratum_name"]

        # Get the weight distributions
        weight_distribution = data_dict["transect"]["biology"]["distributions"]["weight"][
            "aged_length_weight_tbl"
        ].copy()

        # Get transect dataframe
        transect_df = data_dict["transect"]["acoustics"]["adult_transect_df"].copy()

        # Reformat the aged weight proportions tables
        tables = pivot_aged_weight_proportions(weight_distribution, age1_exclusion, stratum_name)

        # Combine the pivot tables with the transect data
        transect_pvt_tables = {
            sex: pivot_dataframe_data(transect_df, tables, sex, stratum_name)
            for sex in ["all", "female", "male"]
        }

        # Write the *.xlsx sheet (full dataset)
        write_aged_dataframe_report(transect_pvt_tables, filename(self)[0])

        # Reduce the datasets
        transect_pvt_reduced_tables = {
            sex: pivot_dataframe_data(
                transect_df.loc[transect_df.biomass > 0.0], tables, sex, stratum_name
            )
            for sex in ["all", "female", "male"]
        }

        # Write the *.xlsx sheet (reduced dataset)
        write_aged_dataframe_report(transect_pvt_reduced_tables, filename(self)[1])

        # Return the filepath name
        return filename(self)[0].as_posix(), filename(self)[1].as_posix()

    def transect_length_age_abundance_report(
        self,
        data: Callable[[pd.DataFrame], pd.DataFrame],
        title: str,
        filename: Callable[[Union[str, Path]], Path],
        **kwargs,
    ):

        # Repivot the datasets for each sex
        # ---- First get the complete aged table
        tables = {
            sex: repivot_table(
                age_length_dataframe=data(self)["aged_abundance_df"].loc[sex, :],
                length_dataframe=data(self)["unaged_abundance_df"].loc[sex, :],
                variable="abundance",
            )
            for sex in ["all", "male", "female"]
        }

        # Write the *.xlsx sheet
        write_age_length_table_report(tables, title, filename(self))

        # Return the filepath name
        return filename(self).as_posix()

    def transect_length_age_biomass_report(
        self,
        data: Callable[[pd.DataFrame], pd.DataFrame],
        title: str,
        filename: Callable[[Union[str, Path]], Path],
        **kwargs,
    ) -> str:

        # Repivot the datasets for each sex
        tables = {
            sex: repivot_table(data(self).loc[sex, :] * 1e-9, "biomass")
            for sex in ["all", "male", "female"]
        }

        # Write the *.xlsx sheet
        write_age_length_table_report(tables, title, filename(self))

        # Return the filepath name
        return filename(self).as_posix()

    def transect_population_results_report(
        self,
        data: Callable[[pd.DataFrame], pd.DataFrame],
        filename: Callable[[Union[str, Path]], Path],
        **kwargs,
    ):

        # Save the *.xlsx sheet
        data(self).to_excel(filename(self), sheet_name="Sheet1", index=None)

        # Return the filepath name
        return filename(self).as_posix()

    @staticmethod
    def report_options():
        """
        Helper method for listing viable report-types
        """

        # Get reference
        options = list(FEATReports.__METHOD_REFERENCE__.keys())

        # Join
        options_str = "\n   ".join(f"-'{report}'" for report in options)

        # Print
        print(f"The following report-types are available for export:\n   {options_str}")

        # Return the list
        return list(options)

    # CLASS REFERENCE DICTIONARIES
    # Reference dictionary: report methods
    __METHOD_REFERENCE__ = {
        "aged_length_haul_counts": {
            "function": aged_length_haul_counts_report,
            "data": lambda v: v.data,
            "title": "Aged Length-Haul Counts ({SEX})",
            "filename": lambda v: Path(v.save_directory) / "aged_length_haul_counts_table.xlsx",
        },
        "kriged_aged_biomass_mesh": {
            "function": kriged_aged_biomass_mesh_report,
            "data": lambda v: v.data,
            "title": None,
            "filename": lambda v: (
                Path(v.save_directory) / "kriged_aged_biomass_mesh_dataframe.xlsx",
                Path(v.save_directory) / "kriged_aged_biomass_reduced_mesh_dataframe.xlsx",
            ),
        },
        "kriged_biomass_mesh": {
            "function": kriged_biomass_mesh_report,
            "data": lambda v: v.data,
            "title": None,
            "filename": lambda v: (
                Path(v.save_directory) / "kriged_biomass_mesh_dataframe.xlsx",
                Path(v.save_directory) / "kriged_biomass_mesh_reduced_dataframe.xlsx",
            ),
        },
        "kriged_length_age_biomass": {
            "function": kriged_length_age_biomass_report,
            "data": (lambda v: v.data.results["kriging"]["tables"]["overall_apportionment_df"]),
            "title": "Kriged Acoustically Weighted Biomass (mmt) ({SEX})",
            "filename": lambda v: Path(v.save_directory) / "kriged_length_age_biomass_table.xlsx",
        },
        "kriging_input": {
            "function": kriging_input_report,
            "data": lambda v: v.data.analysis["transect"]["acoustics"]["adult_transect_df"],
            "title": None,
            "filename": lambda v: Path(v.save_directory) / "kriging_input_dataframe.xlsx",
        },
        "total_length_haul_counts": {
            "function": total_length_haul_counts_report,
            "data": lambda v: v.data,
            "title": "Un-Aged Length-Haul Counts ({SEX})",
            "filename": lambda v: Path(v.save_directory) / "total_length_haul_counts_table.xlsx",
        },
        "transect_aged_biomass": {
            "function": transect_aged_biomass_report,
            "data": lambda v: v.data.analysis,
            "title": None,
            "filename": lambda v: (
                Path(v.save_directory) / "transect_aged_biomass_dataframe.xlsx",
                Path(v.save_directory) / "transect_aged_biomass_reduced_dataframe.xlsx",
            ),
        },
        "transect_length_age_abundance": {
            "function": transect_length_age_abundance_report,
            "data": (
                lambda v: v.data.analysis["transect"]["biology"]["population"]["tables"][
                    "abundance"
                ]
            ),
            "title": "Transect-based Acoustically Weighted Abundance ({SEX})",
            "filename": lambda v: Path(v.save_directory)
            / "transect_length_age_abundance_table.xlsx",
        },
        "transect_length_age_biomass": {
            "function": transect_length_age_biomass_report,
            "data": lambda v: v.data.analysis["transect"]["biology"]["population"]["tables"][
                "biomass"
            ]["aged_biomass_df"],
            "title": "Transect-based Acoustically Weighted Biomass (mmt) ({SEX})",
            "filename": lambda v: Path(v.save_directory) / "transect_length_age_biomass_table.xlsx",
        },
        "transect_population_results": {
            "function": transect_population_results_report,
            "data": lambda v: v.data.analysis["transect"]["acoustics"]["adult_transect_df"],
            "title": None,
            "filename": lambda v: Path(v.save_directory)
            / "transect_population_results_dataframe.xlsx",
        },
    }
