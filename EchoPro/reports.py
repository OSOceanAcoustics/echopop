import pathlib
from typing import List, Union

import geopandas as gpd
import numpy as np
import pandas as pd

from .data_loader.nasc_data import _process_nasc_data


class Reports:
    """
    A Class that writes requested variables to
    consolidated files.

    Parameters
    ----------
    survey : Survey
        Survey object that contains all necessary parameters

    """

    def __init__(self, survey):

        self.survey = survey

    def __process_gear_data(self, df):
        """
        Parameters
        ----------
        df : Pandas Dataframe
            Dataframe holding the gear data
        Returns
        -------
        Processed Dataframe
        """

        df = df[
            [
                "Haul",
                "Average_Footrope_Depth",
                "Surface_Temperature",
                "Gear_Temperature",
                "Average_Wireout",
                "Transect",
                "Net_Height",
            ]
        ].copy()

        # set data types of dataframe
        df = df.astype(
            {
                "Haul": int,
                "Average_Footrope_Depth": np.float64,
                "Surface_Temperature": np.float64,
                "Gear_Temperature": np.float64,
                "Average_Wireout": np.float64,
                "Transect": np.float64,
            }
        )

        if self.EPro.params["exclude_age1"] is False:
            raise NotImplementedError("Including age 1 data has not been implemented!")

        df["Net_Height"] = (
            df["Net_Height"]
            .apply(lambda x: np.nan if type(x) == str else x)
            .astype(float)
        )

        df.set_index("Haul", inplace=True)
        df.sort_index(inplace=True)

        return df

    def __load_gear_data(self):

        # TODO: if this function ends up being used, we need to document and create a check for df

        if self.EPro.params["source"] == 3:

            if self.EPro.params["filename_gear_US"]:
                gear_us_df = pd.read_excel(
                    self.EPro.params["data_root_dir"]
                    + self.EPro.params["filename_gear_US"],
                    sheet_name="biodata_gear",
                )
                gear_us_df = self.__process_gear_data(gear_us_df)
            else:
                gear_us_df = None

            if self.EPro.params["filename_gear_CAN"]:
                gear_can_df = pd.read_excel(
                    self.EPro.params["data_root_dir"]
                    + self.EPro.params["filename_gear_CAN"],
                    sheet_name="biodata_gear_CAN",
                )
                gear_can_df = self.__process_gear_data(gear_can_df)

                # transect offset is set to zero since we will not have overlap transects
                # TODO: Is this a necessary variable? If so, we should make it an input to the function.  # noqa
                CAN_Transect_offset = 0

                # combine US & CAN gear files
                gear_can_df.index = (
                    gear_can_df.index + self.EPro.params["CAN_haul_offset"]
                )  # add haul_offset
                gear_can_df["Transect"] = gear_can_df["Transect"] + CAN_Transect_offset
            else:
                gear_can_df = None

            # combine US & CAN trawl files
            if isinstance(gear_us_df, pd.DataFrame) and isinstance(
                gear_can_df, pd.DataFrame
            ):
                self.EPro.gear_df = pd.concat([gear_us_df, gear_can_df])
            else:
                raise SystemError(
                    "Cannot construct gear_df for source = 3, "
                    "since either US or CAN is not available."
                )

        else:
            raise NotImplementedError(
                f"Source of {self.EPro.params['source']} not implemented yet."
            )

    def __process_trawl_data(self, df):
        """
        Parameters
        ----------
        df : Pandas Dataframe
            Dataframe holding the trawl data
        Returns
        -------
        Processed Dataframe
        """

        df = df[
            [
                "Haul",
                "Haul_Type",
                "Performance_Code",
                "Duration",
                "Distance_Fished",
                "Stratum",
                "EQ_Latitude",
                "EQ_Longitude",
                "Average_Bottom_Depth",
                "Vessel_Log_Start",
                "Vessel_Log_Stop",
            ]
        ].copy()

        if df.dtypes["Average_Bottom_Depth"] == object:
            # remove any instances of >, <, or m from the column Average_Bottom_Depth
            df["Average_Bottom_Depth"] = (
                df["Average_Bottom_Depth"]
                .apply(lambda x: x.strip("><m") if type(x) == str else x)
                .astype(float)
            )

        df["Duration"] = pd.to_timedelta(
            df["Duration"]
        )  # change datatype from string to timedelta

        # set data types of dataframe
        df = df.astype(
            {
                "Haul": int,
                "Haul_Type": int,
                "Performance_Code": np.float64,
                "Distance_Fished": np.float64,
                "Stratum": int,
                "EQ_Latitude": np.float64,
                "EQ_Longitude": np.float64,
                "Average_Bottom_Depth": np.float64,
                "Vessel_Log_Start": np.float64,
                "Vessel_Log_Stop": np.float64,
            }
        )

        df.set_index("Haul", inplace=True)
        df.sort_index(inplace=True)

        if self.EPro.params["exclude_age1"] is False:
            raise NotImplementedError("Including age 1 data has not been implemented!")

        if self.EPro.params["hemisphere"][0] == "N":
            df["EQ_Latitude"] = df["EQ_Latitude"].abs()
        elif self.EPro.params["hemisphere"][0] == "S":
            df["EQ_Latitude"] = -1.0 * (df["EQ_Latitude"].abs())
        else:
            raise ValueError("Wrong N/S Hemisphere provided! \n")

        if self.EPro.params["hemisphere"][1] == "W":
            df["EQ_Longitude"] = -1.0 * (df["EQ_Longitude"].abs())
        elif self.EPro.params["hemisphere"][1] == "E":
            df["EQ_Longitude"] = df["EQ_Longitude"].abs()
        else:
            raise ValueError("Wrong E/W Hemisphere provided! \n")

        return df

    def __load_trawl_data(self):

        # TODO: if this function ends up being used, we need to document and create a check for df

        if self.EPro.params["source"] == 3:

            if self.EPro.params["filename_trawl_US"]:
                trawl_us_df = pd.read_excel(
                    self.EPro.params["data_root_dir"]
                    + self.EPro.params["filename_trawl_US"],
                    converters={"Duration": str},
                    sheet_name="biodata_haul",
                )
                trawl_us_df = self.__process_trawl_data(trawl_us_df)
            else:
                trawl_us_df = None

            if self.EPro.params["filename_trawl_CAN"]:
                trawl_can_df = pd.read_excel(
                    self.EPro.params["data_root_dir"]
                    + self.EPro.params["filename_trawl_CAN"],
                    converters={"Duration": str},
                    sheet_name="biodata_haul_CAN",
                )
                trawl_can_df = self.__process_trawl_data(trawl_can_df)

                trawl_can_df.index = (
                    trawl_can_df.index + self.EPro.params["CAN_haul_offset"]
                )  # add haul_offset

                # change lon to -lon
                trawl_can_df["EQ_Longitude"] = trawl_can_df["EQ_Longitude"].apply(
                    lambda x: -1.0 * x if x > 0 else x
                )
            else:
                trawl_can_df = None

            # combine US & CAN trawl files
            if isinstance(trawl_us_df, pd.DataFrame) and isinstance(
                trawl_can_df, pd.DataFrame
            ):
                self.EPro.trawl_df = pd.concat([trawl_us_df, trawl_can_df])
            else:
                raise SystemError(
                    "Cannot construct trawl_df for source = 3, "
                    "since either US or CAN is not available."
                )

        else:
            raise NotImplementedError(
                f"Source of {self.EPro.params['source']} not implemented yet."
            )

    @staticmethod
    def _write_dfs_to_excel(
        df_list: List[pd.DataFrame],
        sheet_name_list: List[str],
        excel_path: pathlib.Path,
    ) -> None:
        """
        Writes a list of DataFrames to an Excel file using a
        provided output path and associated sheet name.

        Parameters
        ----------
        df_list: list of pd.DataFrame
            A list of DataFrames to write to the Excel file
        sheet_name_list: list of str
            A list of sheet names corresponding to ``df_list``
        excel_path: pathlib.Path
            The path to the Excel file where data should be written
        """

        # write DataFrames to Excel sheet
        with pd.ExcelWriter(excel_path) as writer:

            for i in range(len(df_list)):
                df_list[i].to_excel(writer, sheet_name=sheet_name_list[i])

    def _aged_len_haul_counts_report(self):
        """
        Creates aged length-haul-counts table, which specifies the aged hake
        counts at each length bin from all hauls for male, female, and all.
        Additionally, writes this report to an Excel file.

        Returns
        -------

        """

        pass

    def _write_biomass_ages_report(
        self,
        output_excel_path_all: pathlib.Path,
        output_excel_path_non_zero: pathlib.Path,
        results: gpd.GeoDataFrame,
        results_male: gpd.GeoDataFrame,
        results_female: gpd.GeoDataFrame,
        krig_result: bool,
    ) -> None:
        """
        Generates a biomass at different age bins report based off of the input
        DataFrames and writes it to an Excel file.

        Parameters
        ----------
        output_excel_path_all: pathlib.Path
            The output Excel file path where all reports that include all the data
            should be saved
        output_excel_path_non_zero: pathlib.Path
            The output Excel file path where all reports that include non-zero
            biomass should be saved
        results: gpd.GeoDataFrame
            A GeoDataFrame containing data that includes all genders
        results_male: gpd.GeoDataFrame
            A GeoDataFrame containing data that only includes males
        results_female: gpd.GeoDataFrame
            A GeoDataFrame containing data that only includes females
        krig_result: bool
            It True the data being written corresponds to Kriging based variables,
            else corresponds to Transect based variables
        """

        # specify columns that should be included in the report
        if krig_result:
            lat_lon_names = ["centroid_latitude", "centroid_longitude"]
        else:
            lat_lon_names = ["latitude", "longitude"]
        wanted_columns = (
            lat_lon_names
            + ["stratum_num", "biomass_adult"]
            + [
                "biomass_age_bin_" + str(i + 1)
                for i in range(len(self.survey.params["bio_hake_age_bin"]))
            ]
        )

        # write all results to Excel file
        df_list = [
            results[wanted_columns],
            results_male[wanted_columns],
            results_female[wanted_columns],
        ]
        sheet_names = ["all genders", "male", "female"]
        self._write_dfs_to_excel(df_list, sheet_names, output_excel_path_all)

        # write only output corresponding to non-zero biomass values
        df_list = [
            results[results["biomass_adult"] != 0.0][wanted_columns],
            results_male[results["biomass_adult"] != 0.0][wanted_columns],
            results_female[results["biomass_adult"] != 0.0][wanted_columns],
        ]
        self._write_dfs_to_excel(df_list, sheet_names, output_excel_path_non_zero)

    def _transect_based_core_variables_report(
        self,
        output_excel_path_all: pathlib.Path,
        output_excel_path_non_zero: pathlib.Path,
    ) -> None:
        """
        Generates a report containing core transect based variables
        and writes it to an Excel file.

        Parameters
        ----------
        output_excel_path_all: pathlib.Path
            The output Excel file path where all reports that include all the data
            should be saved
        output_excel_path_non_zero: pathlib.Path
            The output Excel file path where all reports that include non-zero
            NASC should be saved
        """

        # specify column names to grab and their corresponding type
        # TODO: these columns names have not been reviewed, should they be changed?
        nasc_var_types = {
            "Region ID": int,
            "Bottom depth": np.float64,
            "Layer mean depth": np.float64,
            "Layer height": np.float64,
        }

        # obtain NASC data data that was not necessary for core routines
        extra_nasc_df = _process_nasc_data(self.survey, nasc_var_types).set_index(
            self.survey.bio_calc.nasc_df.index
        )

        # set variables to improve readability
        transect_results = self.survey.bio_calc.transect_results_gdf
        transect_results_male = self.survey.bio_calc.transect_results_male_gdf
        transect_results_female = self.survey.bio_calc.transect_results_female_gdf
        stratum_vals = self.survey.bio_calc.nasc_df.stratum_num

        # TODO: change this so that it uses the NASC proportion
        fraction_adult_stratum_df = self.survey.bio_calc.num_fraction_adult_df.loc[
            stratum_vals
        ].values.flatten()

        # define columns to grab from the produced results
        wanted_columns = [
            "latitude",
            "longitude",
            "stratum_num",
            "biomass_adult",
            "biomass_density_adult",
            "numerical_density_adult",
            "abundance_adult",
            "transect_spacing",
            "interval",
        ]
        gender_wanted_columns = [
            "biomass_adult",
            "biomass_density",
            "abundance_adult",
            "numerical_density",
        ]

        # select columns from transect data and rename gender specific columns
        df = transect_results[wanted_columns]
        male_df = transect_results_male[gender_wanted_columns].rename(
            columns={
                "biomass_adult": "biomass_male_adult",
                "biomass_density": "biomass_density_male",
                "abundance_adult": "abundance_male_adult",
                "numerical_density": "numerical_density_male",
            }
        )
        female_df = transect_results_female[gender_wanted_columns].rename(
            columns={
                "biomass_adult": "biomass_female_adult",
                "biomass_density": "biomass_density_female",
                "abundance_adult": "abundance_female_adult",
                "numerical_density": "numerical_density_female",
            }
        )

        # get sig_b over the NASC points
        sig_b = (
            self.survey.bio_calc.strata_sig_b.loc[stratum_vals]
            .reset_index(drop=True)
            .set_axis(self.survey.bio_calc.nasc_df.index)
        )

        # get average weight over the NASC points
        avg_wgt = (
            self.survey.bio_calc.bio_param_df.averaged_weight.loc[stratum_vals]
            .reset_index(drop=True)
            .set_axis(self.survey.bio_calc.nasc_df.index)
        )

        # obtain the NASC for adults
        NASC_adult = self.survey.bio_calc.nasc_df["NASC"] * fraction_adult_stratum_df

        # put all results into one DataFrame
        final_df = pd.concat(
            [
                df,
                male_df,
                female_df,
                extra_nasc_df,
                self.survey.bio_calc.nasc_df[["vessel_log_start", "vessel_log_end"]],
                NASC_adult,
                self.survey.bio_calc.mix_sa_ratio,
                sig_b,
                avg_wgt,
            ],
            axis=1,
        )

        # define column names in the same order as the defined report
        ordered_columns = [
            "Region ID",
            "vessel_log_start",
            "vessel_log_end",
            "latitude",
            "longitude",
            "stratum_num",
            "Bottom depth",
            "NASC",
            "abundance_male_adult",
            "abundance_female_adult",
            "abundance_adult",
            "biomass_male_adult",
            "biomass_female_adult",
            "biomass_adult",
            "numerical_density_male",
            "numerical_density_female",
            "numerical_density_adult",
            "biomass_density_male",
            "biomass_density_female",
            "biomass_density_adult",
            "Layer mean depth",
            "Layer height",
            "transect_spacing",
            "interval",
            "hake_mix_coefficient",
            "sig_bs_haul",
            "averaged_weight",
        ]

        # write all results to Excel file
        df_list = [final_df[ordered_columns]]
        sheet_names = ["Sheet1"]
        self._write_dfs_to_excel(df_list, sheet_names, output_excel_path_all)

        # write only output corresponding to non-zero NASC values
        df_list = [final_df[final_df["NASC"] != 0.0][ordered_columns]]
        self._write_dfs_to_excel(df_list, sheet_names, output_excel_path_non_zero)

    def _total_len_haul_counts_report(self):
        """
        Generates a report for total hake counts at length with
        actually measured length (sampled at both biological sampling
        stations for US/CAN) for all hauls. Additionally, writes this
        report to an Excel file.

        Returns
        -------

        """

        pass

    def _transect_based_len_age_abundance_report(self):
        """
        Generates a report that is a 40 x 21 matrix, with 40 length bins (1st column)
        (L =2, 4, 6, … 80 cm, and ith length bin include length i-1 and i cm) and 21
        age bins (1st row), 1, 2, 3, … 20 year old, and the 21st age bin is un-aged
        abundance. The abundance values are based on the transect results. There are
        three tabs with male, female, and all (male+female). Additionally, writes this
        report to an Excel file.

        Returns
        -------

        """

        pass

    def _transect_based_len_age_biomass_report(self):
        """
        Generates a report that is a 40 x 21 matrix, with 40 length bins (1st column)
        (L =2, 4, 6, … 80 cm, and ith length bin include length i-1 and i cm) and 21
        age bins (1st row), 1, 2, 3, … 20 year old. The biomass values are based on the
        transect results. There are three tabs with male, female, and all (male+female).
        Additionally, writes this report to an Excel file.

        Returns
        -------

        """

        pass

    def create_and_write_reports(self, output_path: Union[str, pathlib.Path]) -> None:
        """
        Constructs Kriging mesh and Transect report DataFrames and writes
        them to Excel files.

        Parameters
        ----------
        output_path: str or pathlib.Path
            The output path where all Excel files should be saved

        # TODO: maybe include an option to overwrite files (default to False)

        # TODO: should we include an option that allows you to generate a specific report?
        """

        if not isinstance(output_path, (str, pathlib.Path)):
            raise TypeError("output_path must be a string or pathlib.Path!")

        if isinstance(output_path, str):
            # convert string path to pathlib.Path
            output_path = pathlib.Path(output_path)

        # check if path exists, if it doesn't create it
        output_path.mkdir(parents=True, exist_ok=True)

        # TODO: perform a check that all necessary DataFrames have been constructed

        self._aged_len_haul_counts_report()

        # self._write_biomass_ages_report(
        #     output_excel_path_all=output_path / "transect_based_aged_output_all.xlsx",
        #     output_excel_path_non_zero=output_path
        #     / "transect_based_aged_output_non_zero.xlsx",
        #     results=self.survey.bio_calc.transect_results_gdf,
        #     results_male=self.survey.bio_calc.transect_results_male_gdf,
        #     results_female=self.survey.bio_calc.transect_results_female_gdf,
        #     krig_result=False,
        # )
        #
        # self._write_biomass_ages_report(
        #     output_excel_path_all=output_path / "kriging_based_aged_output_all.xlsx",
        #     output_excel_path_non_zero=output_path
        #     / "kriging_based_aged_output_non_zero.xlsx",
        #     results=self.survey.bio_calc.kriging_results_gdf,
        #     results_male=self.survey.bio_calc.kriging_results_male_gdf,
        #     results_female=self.survey.bio_calc.kriging_results_female_gdf,
        #     krig_result=True,
        # )

        self._transect_based_core_variables_report(
            output_excel_path_all=output_path / "transect_based_core_output_all.xlsx",
            output_excel_path_non_zero=output_path
            / "transect_based_core_output_non_zero.xlsx",
        )

        self._total_len_haul_counts_report()

        self._transect_based_len_age_abundance_report()

        self._transect_based_len_age_biomass_report()
