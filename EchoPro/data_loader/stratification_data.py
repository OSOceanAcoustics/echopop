from pathlib import Path

import numpy as np
import pandas as pd

from ..utils.input_checks import check_column_names, check_existence_of_file


class LoadStrataData:  # TODO: Does it make sense for this to be a class?
    """
    A Class that loads and processes all
    stratification data. Additionally, it
    calculates the associated mean backscattering
    cross-section for each stratum.

    Parameters
    ----------
    survey : Survey
        An initialized Survey object. Note that any change to
        self.survey will also change this object.
    """

    def __init__(self, survey=None):
        self.survey = survey

        # expected columns for strata Dataframe
        self.strata_cols_types = {
            "stratum_num": int,
            "haul_num": int,
            "fraction_hake": np.float64
        }

        # expected columns for geo strata Dataframe
        self.geo_strata_cols_types = {"stratum_num": int, "Latitude (upper limit)": np.float64}

        self._load_stratification_file()
        self._load_geographic_stratification()

    def _load_stratification_file(self) -> None:
        """
        Loads and checks the stratification file associated
        with relating the stratification to the Haul.

        Returns
        -------
        strata_df : pd.Dataframe
            Dataframe representation of stratification file.
        """

        if self.survey.params["strata_sheetname"] in ["Base KS", "INPFC"]:

            # check existence of the file
            file_path = (
                self.survey.params["data_root_dir"]
                / self.survey.params["strata_filename"]
            )
            check_existence_of_file(file_path)

            # read and check stratification file
            strata_df = pd.read_excel(
                file_path, sheet_name=self.survey.params["strata_sheetname"]
            )
            check_column_names(
                df=strata_df,
                expected_names=set(self.strata_cols_types.keys()),
                path_for_df=file_path
            )

            # extract only those columns that are necessary
            strata_df = strata_df[list(self.strata_cols_types.keys())].copy()

            # set data types of dataframe
            strata_df = strata_df.astype(self.strata_cols_types)

            # set index of dataframe
            strata_df.set_index(["haul_num", "stratum_num"], inplace=True)
            strata_df.sort_index(inplace=True)

            self.survey.strata_df = strata_df
        else:
            raise NotImplementedError("strata_filename has unknown sheet name!")

    def _load_geographic_stratification(self) -> None:
        """
        Loads and checks the geographic stratification
        file defining the stratum and the associated
        Latitude upper limit.

        Returns
        -------
        geo_strata_df : pd.Dataframe
            Dataframe representation of geographic stratification file.
        """

        if self.survey.params["geo_strata_sheetname"] in ["stratification1", "INPFC"]:

            # check existence of the file
            file_path = (
                self.survey.params["data_root_dir"]
                / self.survey.params["geo_strata_filename"]
            )
            check_existence_of_file(file_path)

            # read and check geographic stratification file
            geo_strata_df = pd.read_excel(
                file_path, sheet_name=self.survey.params["geo_strata_sheetname"]
            )
            check_column_names(
                df=geo_strata_df,
                expected_names=set(self.geo_strata_cols_types.keys()),
                path_for_df=file_path
            )

            # extract only those columns that are necessary
            geo_strata_df = geo_strata_df[list(self.geo_strata_cols_types.keys())].copy()

            # set data types of dataframe
            geo_strata_df = geo_strata_df.astype(self.geo_strata_cols_types)

            self.survey.geo_strata_df = geo_strata_df
        else:
            raise NotImplementedError("geo_strata_filename has unknown sheet name!")
