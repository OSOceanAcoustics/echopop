import pathlib
from typing import List, Tuple, Union

import geopandas as gpd
import numpy as np
import pandas as pd

from .data_loader.nasc_data import _process_nasc_data
from .utils.binning import get_bin_ind


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

    # def __process_gear_data(self, df):
    #     """
    #     Parameters
    #     ----------
    #     df : Pandas Dataframe
    #         Dataframe holding the gear data
    #     Returns
    #     -------
    #     Processed Dataframe
    #     """
    #
    #     df = df[
    #         [
    #             "Haul",
    #             "Average_Footrope_Depth",
    #             "Surface_Temperature",
    #             "Gear_Temperature",
    #             "Average_Wireout",
    #             "Transect",
    #             "Net_Height",
    #         ]
    #     ].copy()
    #
    #     # set data types of dataframe
    #     df = df.astype(
    #         {
    #             "Haul": int,
    #             "Average_Footrope_Depth": np.float64,
    #             "Surface_Temperature": np.float64,
    #             "Gear_Temperature": np.float64,
    #             "Average_Wireout": np.float64,
    #             "Transect": np.float64,
    #         }
    #     )
    #
    #     if self.EPro.params["exclude_age1"] is False:
    #         raise NotImplementedError("Including age 1 data has not been implemented!")
    #
    #     df["Net_Height"] = (
    #         df["Net_Height"]
    #         .apply(lambda x: np.nan if type(x) == str else x)
    #         .astype(float)
    #     )
    #
    #     df.set_index("Haul", inplace=True)
    #     df.sort_index(inplace=True)
    #
    #     return df
    #
    # def __load_gear_data(self):
    #
    #     # TODO: if this function ends up being used, we need to document and create a check for df
    #
    #     if self.EPro.params["source"] == 3:
    #
    #         if self.EPro.params["filename_gear_US"]:
    #             gear_us_df = pd.read_excel(
    #                 self.EPro.params["data_root_dir"]
    #                 + self.EPro.params["filename_gear_US"],
    #                 sheet_name="biodata_gear",
    #             )
    #             gear_us_df = self.__process_gear_data(gear_us_df)
    #         else:
    #             gear_us_df = None
    #
    #         if self.EPro.params["filename_gear_CAN"]:
    #             gear_can_df = pd.read_excel(
    #                 self.EPro.params["data_root_dir"]
    #                 + self.EPro.params["filename_gear_CAN"],
    #                 sheet_name="biodata_gear_CAN",
    #             )
    #             gear_can_df = self.__process_gear_data(gear_can_df)
    #
    #             # transect offset is set to zero since we will not have overlap transects
    #             # TODO: Is this a necessary variable? If so, we should make it an input to the function.  # noqa
    #             CAN_Transect_offset = 0
    #
    #             # combine US & CAN gear files
    #             gear_can_df.index = (
    #                 gear_can_df.index + self.EPro.params["CAN_haul_offset"]
    #             )  # add haul_offset
    #             gear_can_df["Transect"] = gear_can_df["Transect"] + CAN_Transect_offset
    #         else:
    #             gear_can_df = None
    #
    #         # combine US & CAN trawl files
    #         if isinstance(gear_us_df, pd.DataFrame) and isinstance(
    #             gear_can_df, pd.DataFrame
    #         ):
    #             self.EPro.gear_df = pd.concat([gear_us_df, gear_can_df])
    #         else:
    #             raise SystemError(
    #                 "Cannot construct gear_df for source = 3, "
    #                 "since either US or CAN is not available."
    #             )
    #
    #     else:
    #         raise NotImplementedError(
    #             f"Source of {self.EPro.params['source']} not implemented yet."
    #         )
    #
    # def __process_trawl_data(self, df):
    #     """
    #     Parameters
    #     ----------
    #     df : Pandas Dataframe
    #         Dataframe holding the trawl data
    #     Returns
    #     -------
    #     Processed Dataframe
    #     """
    #
    #     df = df[
    #         [
    #             "Haul",
    #             "Haul_Type",
    #             "Performance_Code",
    #             "Duration",
    #             "Distance_Fished",
    #             "Stratum",
    #             "EQ_Latitude",
    #             "EQ_Longitude",
    #             "Average_Bottom_Depth",
    #             "Vessel_Log_Start",
    #             "Vessel_Log_Stop",
    #         ]
    #     ].copy()
    #
    #     if df.dtypes["Average_Bottom_Depth"] == object:
    #         # remove any instances of >, <, or m from the column Average_Bottom_Depth
    #         df["Average_Bottom_Depth"] = (
    #             df["Average_Bottom_Depth"]
    #             .apply(lambda x: x.strip("><m") if type(x) == str else x)
    #             .astype(float)
    #         )
    #
    #     df["Duration"] = pd.to_timedelta(
    #         df["Duration"]
    #     )  # change datatype from string to timedelta
    #
    #     # set data types of dataframe
    #     df = df.astype(
    #         {
    #             "Haul": int,
    #             "Haul_Type": int,
    #             "Performance_Code": np.float64,
    #             "Distance_Fished": np.float64,
    #             "Stratum": int,
    #             "EQ_Latitude": np.float64,
    #             "EQ_Longitude": np.float64,
    #             "Average_Bottom_Depth": np.float64,
    #             "Vessel_Log_Start": np.float64,
    #             "Vessel_Log_Stop": np.float64,
    #         }
    #     )
    #
    #     df.set_index("Haul", inplace=True)
    #     df.sort_index(inplace=True)
    #
    #     if self.EPro.params["exclude_age1"] is False:
    #         raise NotImplementedError("Including age 1 data has not been implemented!")
    #
    #     if self.EPro.params["hemisphere"][0] == "N":
    #         df["EQ_Latitude"] = df["EQ_Latitude"].abs()
    #     elif self.EPro.params["hemisphere"][0] == "S":
    #         df["EQ_Latitude"] = -1.0 * (df["EQ_Latitude"].abs())
    #     else:
    #         raise ValueError("Wrong N/S Hemisphere provided! \n")
    #
    #     if self.EPro.params["hemisphere"][1] == "W":
    #         df["EQ_Longitude"] = -1.0 * (df["EQ_Longitude"].abs())
    #     elif self.EPro.params["hemisphere"][1] == "E":
    #         df["EQ_Longitude"] = df["EQ_Longitude"].abs()
    #     else:
    #         raise ValueError("Wrong E/W Hemisphere provided! \n")
    #
    #     return df
    #
    # def __load_trawl_data(self):
    #
    #     # TODO: if this function ends up being used, we need to document and create a check for df
    #
    #     if self.EPro.params["source"] == 3:
    #
    #         if self.EPro.params["filename_trawl_US"]:
    #             trawl_us_df = pd.read_excel(
    #                 self.EPro.params["data_root_dir"]
    #                 + self.EPro.params["filename_trawl_US"],
    #                 converters={"Duration": str},
    #                 sheet_name="biodata_haul",
    #             )
    #             trawl_us_df = self.__process_trawl_data(trawl_us_df)
    #         else:
    #             trawl_us_df = None
    #
    #         if self.EPro.params["filename_trawl_CAN"]:
    #             trawl_can_df = pd.read_excel(
    #                 self.EPro.params["data_root_dir"]
    #                 + self.EPro.params["filename_trawl_CAN"],
    #                 converters={"Duration": str},
    #                 sheet_name="biodata_haul_CAN",
    #             )
    #             trawl_can_df = self.__process_trawl_data(trawl_can_df)
    #
    #             trawl_can_df.index = (
    #                 trawl_can_df.index + self.EPro.params["CAN_haul_offset"]
    #             )  # add haul_offset
    #
    #             # change lon to -lon
    #             trawl_can_df["EQ_Longitude"] = trawl_can_df["EQ_Longitude"].apply(
    #                 lambda x: -1.0 * x if x > 0 else x
    #             )
    #         else:
    #             trawl_can_df = None
    #
    #         # combine US & CAN trawl files
    #         if isinstance(trawl_us_df, pd.DataFrame) and isinstance(
    #             trawl_can_df, pd.DataFrame
    #         ):
    #             self.EPro.trawl_df = pd.concat([trawl_us_df, trawl_can_df])
    #         else:
    #             raise SystemError(
    #                 "Cannot construct trawl_df for source = 3, "
    #                 "since either US or CAN is not available."
    #             )
    #
    #     else:
    #         raise NotImplementedError(
    #             f"Source of {self.EPro.params['source']} not implemented yet."
    #         )

    def _get_bin_count(
        self, df: pd.DataFrame, haul: int, len_cnt_exists: bool
    ) -> np.ndarray:
        """
        Obtains the number of animals in each length bin for the
        given length data in ``df`` and provided ``haul``.

        Parameters
        ----------
        df: pd.DataFrame
            A DataFrame corresponding to either the length or specimen data, which
            has the haul as its index
        haul: int
            The haul number
        len_cnt_exists: bool
            If True, ``df`` contains ``length_count``, otherwise it does not

        Returns
        -------
        bin_cnt: np.ndarray
            The number of animals in each bin

        Notes
        -----
        The length bins are determined by the parameter
        ``self.survey.params["bio_hake_len_bin"]``.
        """

        # obtain length data
        length_data = df.loc[haul]["length"]

        if len_cnt_exists:
            # obtain length count data
            length_count_data = df.loc[haul]["length_count"]

        # get numpy arrays of data accounting for single values
        if not isinstance(length_data, pd.Series):
            length_data = np.array([length_data])
            if len_cnt_exists:
                length_count_data = np.array([length_count_data])

        else:
            length_data = length_data.values.flatten()
            if len_cnt_exists:
                length_count_data = length_count_data.values.flatten()

        # get bin indices of length data
        len_bin_ind = get_bin_ind(length_data, self.survey.params["bio_hake_len_bin"])

        # get total number of lengths in a bin
        if len_cnt_exists:
            bin_cnt = np.array([length_count_data[i].sum() for i in len_bin_ind])
        else:
            bin_cnt = np.array([i.shape[0] for i in len_bin_ind])

        return bin_cnt

    def _bin_len_by_haul_df(
        self,
        all_hauls: np.ndarray,
        df_len: pd.DataFrame,
        len_uniq_haul: np.ndarray,
        df_spec: pd.DataFrame,
        spec_uniq_haul: np.ndarray,
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        For each haul in ``haul_nums`` bin the length data
        contained the provided DataFrames.

        Parameters
        ----------
        all_hauls: np.ndarray
            All haul numbers that should be included (will be columns
            of returned DataFrames)
        df_len: pd.DataFrame
            A DataFrame corresponding to the length data (must contain
            ``length`` and ``length_counts`` columns)
        len_uniq_haul: np.ndarray
            All unique haul numbers in ``df_len``
        df_spec: pd.DataFrame
            A DataFrame corresponding to the specimen data (must contain
            the ``length`` column)
        spec_uniq_haul: np.ndarray
            All unique haul numbers in ``df_spec``

        Returns
        -------
        len_bin_haul_df: pd.DataFrame
            A DataFrame containing the number of animals in each length bin
            and haul for the length data
        spec_bin_haul_df: pd.DataFrame
            A DataFrame containing the number of animals in each length bin
            and haul for the specimen data
        """

        # DataFrame containing the length bin data for each haul
        len_bin_haul_df = pd.DataFrame(
            data=0,
            columns=list(all_hauls),
            index=self.survey.params["bio_hake_len_bin"],
            dtype=np.int64,
        )

        # DataFrame containing the specimen bin data for each haul
        spec_bin_haul_df = pd.DataFrame(
            data=0,
            columns=list(all_hauls),
            index=self.survey.params["bio_hake_len_bin"],
            dtype=np.int64,
        )

        # name index
        len_bin_haul_df.index.name = "length_bin"
        spec_bin_haul_df.index.name = "length_bin"

        for haul in all_hauls:

            if haul in len_uniq_haul:
                # get the number of animals in each length bin
                len_bin_cnt = self._get_bin_count(df_len, haul, len_cnt_exists=True)

                # store data
                len_bin_haul_df[haul] = len_bin_cnt

            if haul in spec_uniq_haul:
                # get the number of animals in each length bin
                len_bin_cnt = self._get_bin_count(df_spec, haul, len_cnt_exists=False)

                # store data
                spec_bin_haul_df[haul] = len_bin_cnt

        return len_bin_haul_df, spec_bin_haul_df

    def _bin_len_by_haul_all_dfs(
        self,
    ) -> Tuple[
        pd.DataFrame,
        pd.DataFrame,
        pd.DataFrame,
        pd.DataFrame,
        pd.DataFrame,
        pd.DataFrame,
    ]:
        """
        For each haul in ``haul_nums`` bin the length data
        for the length and specimen data (completed for males,
        females, and all genders).

        Returns
        -------
        len_haul_all_df: pd.DataFrame
            A DataFrame containing the number of animals in each length bin
            and haul for the length data and all genders
        spec_haul_all_df: pd.DataFrame
            A DataFrame containing the number of animals in each length bin
            and haul for the specimen data and all genders
        len_haul_all_df_male: pd.DataFrame
            A DataFrame containing the number of animals in each length bin
            and haul for the length data and males
        spec_haul_all_df_male: pd.DataFrame
            A DataFrame containing the number of animals in each length bin
            and haul for the specimen data and males
        len_haul_all_df_female: pd.DataFrame
            A DataFrame containing the number of animals in each length bin
            and haul for the length data and females
        spec_haul_all_df_female: pd.DataFrame
            A DataFrame containing the number of animals in each length bin
            and haul for the specimen data and females
        """

        # create variables to improve readability
        len_df = self.survey.length_df.dropna(how="all")
        spec_df = self.survey.specimen_df.dropna(subset=["age"])

        # unique hauls in length data
        len_uniq_haul = len_df.index.unique().values

        # unique hauls in specimen data
        spec_uniq_haul = spec_df.index.unique().values

        # get all haul numbers
        all_hauls = np.union1d(len_uniq_haul, spec_uniq_haul)

        # get binned lengths at each haul for all data
        len_haul_all_df, spec_haul_all_df = self._bin_len_by_haul_df(
            all_hauls, len_df, len_uniq_haul, spec_df, spec_uniq_haul
        )

        # get DataFrame values corresponding to males
        len_df_male = len_df[len_df["sex"] == 1]
        spec_df_male = spec_df[spec_df["sex"] == 1]

        # get binned lengths at each haul for data corresponding to males
        len_haul_all_df_male, spec_haul_all_df_male = self._bin_len_by_haul_df(
            all_hauls,
            len_df_male,
            len_df_male.index.unique().values,
            spec_df_male,
            spec_df_male.index.unique().values,
        )

        # get DataFrame values corresponding to males
        len_df_female = len_df[len_df["sex"] == 2]
        spec_df_female = spec_df[spec_df["sex"] == 2]

        # get binned lengths at each haul for data corresponding to females
        len_haul_all_df_female, spec_haul_all_df_female = self._bin_len_by_haul_df(
            all_hauls,
            len_df_female,
            len_df_female.index.unique().values,
            spec_df_female,
            spec_df_female.index.unique().values,
        )

        return (
            len_haul_all_df,
            spec_haul_all_df,
            len_haul_all_df_male,
            spec_haul_all_df_male,
            len_haul_all_df_female,
            spec_haul_all_df_female,
        )

    @staticmethod
    def add_len_bin_haul_total(df_list: List[pd.DataFrame]) -> None:
        """
        A helper function that adds ``length_bin_total`` and ``length_haul_total``,
        which are the number of animals for each length bin and haul, respectively
        to each DataFrame in ``df_list``.

        Parameters
        ----------
        df_list: list of pd.DataFrame
            All DataFrames where totals should be calculated and assigned

        Notes
        -----
        The totals are directly applied to the DataFrames in ``df_list``.
        """

        for df in df_list:

            # get totals for the length bins and hauls
            bin_total = df.sum(axis=1)
            haul_total = df.sum(axis=0)

            # assign totals to DataFrame
            df["length_bin_total"] = bin_total
            df.loc["length_haul_total"] = haul_total

    def _get_adult_NASC(self, stratum_vals: np.ndarray) -> pd.Series:
        """
        Computes NASC values corresponding to the adult animal population.

        Parameters
        ----------
        stratum_vals: np.ndarray
            The ``stratum_num`` column of the NASC DataFrame

        Returns
        -------
        pd.Series
            A Series representing the NASC for the adult animal population
        """

        # get the normalized length-age distribution
        len_age_dist_all_norm = (
            self.survey.bio_calc.len_age_dist_all
            / self.survey.bio_calc.len_age_dist_all.sum(dim=["len_bin", "age_bin"])
        )

        # create adult NASC proportion coefficient
        nasc_fraction_adult_df = pd.DataFrame(
            columns=["val"], index=len_age_dist_all_norm.stratum, dtype=np.float64
        )

        for i in len_age_dist_all_norm.stratum.values:
            sig_bs_aged_ave = np.sum(
                self.survey.params["sig_b_coef"]
                * np.matmul(
                    (self.survey.params["bio_hake_len_bin"] ** 2),
                    len_age_dist_all_norm.sel(stratum=i).values,
                )
            )

            temp = self.survey.params["sig_b_coef"] * np.matmul(
                (self.survey.params["bio_hake_len_bin"] ** 2),
                len_age_dist_all_norm.sel(stratum=i).isel(age_bin=0).values,
            )

            age1_nasc_proportion = temp / sig_bs_aged_ave

            nasc_fraction_adult_df.loc[i] = abs(1.0 - age1_nasc_proportion)

        # obtain the adult NASC proportion coefficient for each stratum value
        fraction_adult_stratum_df = nasc_fraction_adult_df.loc[
            stratum_vals
        ].values.flatten()

        # obtain the NASC for adults
        NASC_adult = self.survey.bio_calc.nasc_df["NASC"] * fraction_adult_stratum_df
        NASC_adult.name = "NASC_adult"
        return NASC_adult

    @staticmethod
    def _write_dfs_to_excel(
        df_list: List[pd.DataFrame],
        sheet_name_list: List[str],
        excel_path: pathlib.Path,
        include_index: bool = True,
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
        include_index: bool, default=True
            If True, the index will be included in the Excel sheet, else it won't be
        """

        # write DataFrames to Excel sheet
        with pd.ExcelWriter(excel_path) as writer:

            for i in range(len(df_list)):
                df_list[i].to_excel(
                    writer, sheet_name=sheet_name_list[i], index=include_index
                )

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
        NASC_adult: pd.Series,
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
        NASC_adult: pd.Series
            A Series defining the NASC values corresponding to the adult population
            at each transect number defined in ``self.survey.bio_calc.nasc_df``
        """

        # specify column names to grab and their corresponding type
        # TODO: these columns names have not been reviewed, should they be changed?
        nasc_var_types = {
            "Region ID": int,
            "Bottom depth": np.float64,
            "Layer mean depth": np.float64,
            "Layer height": np.float64,
        }

        # obtain NASC data that was not necessary for core routines
        extra_nasc_df = _process_nasc_data(self.survey, nasc_var_types).set_index(
            self.survey.bio_calc.nasc_df.index
        )

        # set variables to improve readability
        transect_results = self.survey.bio_calc.transect_results_gdf
        transect_results_male = self.survey.bio_calc.transect_results_male_gdf
        transect_results_female = self.survey.bio_calc.transect_results_female_gdf
        stratum_vals = self.survey.bio_calc.nasc_df.stratum_num

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
            "NASC_adult",
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
        df_list = [final_df[final_df["NASC_adult"] != 0.0][ordered_columns]]
        self._write_dfs_to_excel(df_list, sheet_names, output_excel_path_non_zero)

    def _len_haul_count_reports(
        self,
        output_excel_path_specimen: pathlib.Path,
        output_excel_path_total: pathlib.Path,
    ) -> None:
        """
        Generates two reports for hake counts at each length bin for all hauls.
        One report uses only the specimen data and the other uses both the
        length and specimen data. Additionally, writes these reports to an Excel file.

        Parameters
        ----------
        output_excel_path_specimen: pathlib.Path
            The output Excel file path where the report corresponding to the specimen
            data should be saved
        output_excel_path_total
            The output Excel file path where the report corresponding to the specimen
            and length data should be saved
        """

        # obtain length counts at each haul for the length and specimen data
        (
            len_haul_all_df,
            spec_haul_all_df,
            len_haul_all_df_male,
            spec_haul_all_df_male,
            len_haul_all_df_female,
            spec_haul_all_df_female,
        ) = self._bin_len_by_haul_all_dfs()

        # add the length bin total and length haul total values to DataFrames
        self.add_len_bin_haul_total(
            [
                len_haul_all_df,
                spec_haul_all_df,
                len_haul_all_df_male,
                spec_haul_all_df_male,
                len_haul_all_df_female,
                spec_haul_all_df_female,
            ]
        )

        # write the reports corresponding to the specimen data to Excel file
        df_list = [spec_haul_all_df, spec_haul_all_df_male, spec_haul_all_df_female]
        sheet_names = ["all genders", "male", "female"]
        self._write_dfs_to_excel(df_list, sheet_names, output_excel_path_specimen)

        # combine results for length and specimen data
        total_haul_all_df = len_haul_all_df + spec_haul_all_df
        total_haul_male_df = len_haul_all_df_male + spec_haul_all_df_male
        total_haul_female_df = len_haul_all_df_female + spec_haul_all_df_female

        # write the reports corresponding to the specimen and length data to Excel file
        df_list = [total_haul_all_df, total_haul_male_df, total_haul_female_df]
        self._write_dfs_to_excel(df_list, sheet_names, output_excel_path_total)

    def _redistribute_age1_data(self, Len_Age_Matrix_Acoust):
        # TODO: document!

        # get the sum of the values for each age bin
        age_bin_sum = Len_Age_Matrix_Acoust.sum("len_bin").values

        # fill in age 1 data with zero
        Len_Age_Matrix_Acoust.values[:, 0] = 0.0

        # redistribute age 1 data to the rest of the age bins
        Len_Age_Matrix_Acoust.values[:, 1:] += (
            age_bin_sum[0] * Len_Age_Matrix_Acoust.values[:, 1:] / age_bin_sum[1:].sum()
        )

    def _get_len_age_abundance(self, abundance_df, ds, sex, kriging_vals: bool):
        # TODO: document!

        # obtain only those strata that are defined in abundance_df
        defined_stratum = abundance_df.index.unique().values

        # compute the total number of animals in both stations
        total_N = ds.num_M + ds.num_F + ds.station_1_N

        # compute the proportion of the sex using station 2 data
        Len_Age_sex_proportion = ds[f"num_{sex}"] / total_N

        # compute the proportion of the sex using station 1 data
        # TODO: add this in to get the unaged bin
        Len_sex_proportion = ds[f"station_1_N_{sex}"] / total_N

        # get the normalized distribution of data in station 2
        Len_Age_key_sex_norm = ds[f"len_age_dist_{sex}"] / ds[
            f"len_age_dist_{sex}"
        ].sum(["len_bin", "age_bin"])

        # sum together all abundance values in each stratum
        # TODO: we should probably rename ds coordinate to stratum_num
        if kriging_vals:
            abundance_sum_stratum = (
                abundance_df.groupby(level=0)
                .sum()["abundance_adult"]
                .to_xarray()
                .rename({"stratum_num": "stratum"})
            )
        else:
            abundance_sum_stratum = (
                abundance_df.groupby(level=0)
                .sum()["abundance"]
                .to_xarray()
                .rename({"stratum_num": "stratum"})
            )

        # get the abundance for the sex for each stratum (station 1)
        N_len = abundance_sum_stratum * Len_sex_proportion.sel(stratum=defined_stratum)

        # get the abundance for the sex for each stratum (station 2)
        N_len_age = abundance_sum_stratum * Len_Age_sex_proportion.sel(
            stratum=defined_stratum
        )

        Len_Age_Matrix_Acoust_unaged = (
            N_len
            * ds[f"len_dist_station1_normalized_{sex}"].sel(stratum=defined_stratum)
        ).sum("stratum")

        # get the abundance for the sex at each length and age bin
        Len_Age_Matrix_Acoust = (
            N_len_age * Len_Age_key_sex_norm.sel(stratum=defined_stratum)
        ).sum("stratum")

        # redistribute the age 1 data if Kriging values are being used
        if self.survey.params["exclude_age1"] and kriging_vals:

            self._redistribute_age1_data(Len_Age_Matrix_Acoust)

        return Len_Age_Matrix_Acoust, Len_Age_Matrix_Acoust_unaged

    def _get_len_age_biomass(self, biomass_df, ds, kriging_vals: bool):
        # TODO: document!

        # obtain only those strata that are defined in biomass_df
        defined_stratum = biomass_df.index.unique().values

        # obtain the total biomass for each stratum
        if kriging_vals:
            biomass_sum_stratum = (
                biomass_df.groupby(level=0)
                .sum()["biomass_adult"]
                .to_xarray()
                .rename({"stratum_num": "stratum"})
            )
        else:
            biomass_sum_stratum = (
                biomass_df.groupby(level=0)
                .sum()["biomass"]
                .to_xarray()
                .rename({"stratum_num": "stratum"})
            )

        # get the abundance for the sex at each length and age bin
        Len_Age_Matrix_biomass = (
            biomass_sum_stratum
            * ds.len_age_weight_dist_all_normalized.sel(stratum=defined_stratum)
        ).sum("stratum")

        # redistribute the age 1 data if Kriging values are being used
        # if self.survey.params["exclude_age1"] and kriging_vals:
        #     self._redistribute_age1_data(Len_Age_Matrix_biomass)
        # TODO: need to redistribute the data in a special way for biomass

        return Len_Age_Matrix_biomass

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

        # abundance_df = survey_2019.bio_calc.transect_results_gdf[["abundance", "stratum_num"]]

        # TODO: for Kriging
        # abundance_df = survey_2019.bio_calc.kriging_results_gdf[
        # ["abundance_adult", "stratum_num"]]

        # temp_M, temp_unaged_M = self._get_len_age_abundance(abundance_df, krig.ds, sex="M",
        # kriging_vals=False)
        # temp_F, temp_unaged_F = self._get_len_age_abundance(abundance_df, krig.ds, sex="F",
        # kriging_vals=False)

        # temp_F.to_pandas()

        # temp_total = temp_M + temp_F

        # TODO: Take note that EchoPro Matlab does not correctly account for age bins, thus
        #  the output produced for ``Un-aged`` column is not correct when an age bin end is not 20!

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

        # biomass_df = survey_2019.bio_calc.transect_results_male_gdf[["biomass", "stratum_num"]]
        # biomass_df = biomass_df.reset_index(drop=True).set_index("stratum_num")
        #
        # temp_biomass_M = reports._get_len_age_biomass(biomass_df, krig.ds, kriging_vals=False)
        # temp_biomass_M = temp_biomass_M * 1e-9

        # biomass_df = survey_2019.bio_calc.transect_results_female_gdf[["biomass", "stratum_num"]]
        # biomass_df = biomass_df.reset_index(drop=True).set_index("stratum_num")
        #
        # temp_biomass_F = reports._get_len_age_biomass(biomass_df, krig.ds, kriging_vals=False)
        # temp_biomass_F = temp_biomass_F * 1e-9
        #
        # biomass_df = survey_2019.bio_calc.transect_results_gdf[["biomass", "stratum_num"]]
        # biomass_df = biomass_df.reset_index(drop=True).set_index("stratum_num")
        #
        # temp_biomass_all = reports._get_len_age_biomass(biomass_df, krig.ds, kriging_vals=False)
        # temp_biomass_all = temp_biomass_all * 1e-9

        pass

    def _kriging_based_core_variables_report(
        self,
        output_excel_path_all: pathlib.Path,
        output_excel_path_non_zero: pathlib.Path,
    ) -> None:
        """
        Generates a report containing core Kriging based variables
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

        # set variables to improve readability
        krig_results = self.survey.bio_calc.kriging_results_gdf
        krig_results_male = self.survey.bio_calc.kriging_results_male_gdf
        krig_results_female = self.survey.bio_calc.kriging_results_female_gdf

        # define columns grab from the produced results
        wanted_columns = [
            "centroid_latitude",
            "centroid_longitude",
            "stratum_num",
            "NASC",
            "biomass_adult",
            "abundance_adult",
            "sig_b",
            "biomass_adult_cell_CV",
        ]
        gender_wanted_columns = ["biomass_adult", "abundance_adult"]

        # collect wanted columns from results and rename gender based results
        df = krig_results[wanted_columns]
        male_df = krig_results_male[gender_wanted_columns].rename(
            columns={
                "biomass_adult": "biomass_male_adult",
                "abundance_adult": "abundance_male_adult",
            }
        )
        female_df = krig_results_female[gender_wanted_columns].rename(
            columns={
                "biomass_adult": "biomass_female_adult",
                "abundance_adult": "abundance_female_adult",
            }
        )

        # put together all wanted results
        final_df = pd.concat([df, male_df, female_df], axis=1)

        # define column names in the same order as the defined report
        ordered_columns = [
            "centroid_latitude",
            "centroid_longitude",
            "stratum_num",
            "NASC",
            "abundance_male_adult",
            "abundance_female_adult",
            "abundance_adult",
            "biomass_male_adult",
            "biomass_female_adult",
            "biomass_adult",
            "sig_b",
            "biomass_adult_cell_CV",
        ]

        # write all results to Excel file
        df_list = [final_df[ordered_columns]]
        sheet_names = ["Sheet1"]
        self._write_dfs_to_excel(
            df_list, sheet_names, output_excel_path_all, include_index=False
        )

        # write only output corresponding to non-zero NASC values
        df_list = [final_df[final_df["NASC"] != 0.0][ordered_columns]]
        self._write_dfs_to_excel(
            df_list, sheet_names, output_excel_path_non_zero, include_index=False
        )

    def _kriging_input_report(
        self, NASC_adult: pd.Series, output_excel_path: pathlib.Path
    ) -> None:
        """
        Creates a report that contains input used for Kriging.

        Parameters
        ----------
        NASC_adult: pd.Series
            A Series defining the NASC values corresponding to the adult population
            at each transect number defined in ``self.survey.bio_calc.nasc_df``
        output_excel_path: pathlib.Path
            The output Excel file path where the report should be saved
        """

        # put together all wanted results
        final_df = pd.concat(
            [
                self.survey.bio_calc.transect_results_gdf[["latitude", "longitude"]],
                self.survey.bio_calc.transect_results_gdf["biomass_density_adult"],
                NASC_adult,
                self.survey.bio_calc.transect_results_gdf["numerical_density_adult"],
            ],
            axis=1,
        )

        # write all results to Excel file
        df_list = [final_df]
        sheet_names = ["Sheet1"]
        self._write_dfs_to_excel(
            df_list, sheet_names, output_excel_path, include_index=False
        )

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

        NASC_adult = self._get_adult_NASC(self.survey.bio_calc.nasc_df.stratum_num)

        # self._transect_based_core_variables_report(
        #     output_excel_path_all=output_path / "transect_based_core_output_all.xlsx",
        #     output_excel_path_non_zero=output_path
        #     / "transect_based_core_output_non_zero.xlsx", NASC_adult=NASC_adult
        # )

        # self._transect_based_len_age_abundance_report()
        # self._transect_based_len_age_biomass_report()

        # self._kriging_based_core_variables_report(
        #     output_excel_path_all=output_path / "kriging_based_core_output_all.xlsx",
        #     output_excel_path_non_zero=output_path
        #     / "kriging_based_core_output_non_zero.xlsx",
        # )

        self._kriging_input_report(
            NASC_adult=NASC_adult, output_excel_path=output_path / "kriging_input.xlsx"
        )

        # self._len_haul_count_reports(
        #     output_excel_path_specimen=output_path / "specimen_length_counts_haul.xlsx",
        #     output_excel_path_total=output_path / "total_length_counts_haul.xlsx",
        # )
