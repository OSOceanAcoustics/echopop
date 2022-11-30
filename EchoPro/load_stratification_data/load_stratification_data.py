import numpy as np
import pandas as pd
from pathlib import Path
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
        self.strata_cols = {'stratum_num', 'haul_num', 'fraction_hake'}

        # expected columns for geo strata Dataframe
        self.geo_strata_cols = {'stratum_num', 'Latitude (upper limit)'}

        self._load_stratification_file()
        self._load_geographic_stratification()

    def _check_strata_df(self, strata_df: pd.DataFrame, df_path: Path) -> None:
        """
        Ensures that the appropriate columns are
        contained in the stratification Dataframe.

        Parameters
        ----------
        strata_df: pd.DataFrame
            The constructed Strata DataFrame
        df_path: Path
            The path to the Excel file used to construct the DataFrame
        """

        # TODO: should we add more in-depth checks here?

        check_column_names(df=strata_df, expected_names=self.strata_cols, path_for_df=df_path)

    def _load_stratification_file(self) -> None:
        """
        Loads and checks the stratification file associated
        with relating the stratification to the Haul.

        Returns
        -------
        strata_df : pd.Dataframe
            Dataframe representation of stratification file.
        """

        if self.survey.params['strata_sheetname'] in ['Base KS', 'INPFC']:

            # check existence of the file
            file_path = self.survey.params['data_root_dir'] / self.survey.params['strata_filename']
            check_existence_of_file(file_path)

            # read and check stratification file
            strata_df = pd.read_excel(file_path, sheet_name=self.survey.params['strata_sheetname'])
            self._check_strata_df(strata_df, file_path)

            # extract only those columns that are necessary
            strata_df = strata_df[['stratum_num', 'haul_num', 'fraction_hake']].copy()

            # set data types of dataframe
            strata_df = strata_df.astype({'stratum_num': int, 'haul_num': int,
                                          'fraction_hake': np.float64})

            # set index of dataframe
            strata_df.set_index(['haul_num', 'stratum_num'], inplace=True)
            strata_df.sort_index(inplace=True)

            self.survey.strata_df = strata_df

        else:
            raise NotImplementedError(f"strata_filename has unknown sheet name!")

    def _check_geo_strata_df(self, geo_strata_df: pd.DataFrame, df_path: Path) -> None:
        """
        Ensures that the appropriate columns are
        contained in the geographic stratification Dataframe.

        Parameters
        ----------
        geo_strata_df: pd.DataFrame
            The constructed Geo Strata DataFrame
        df_path: Path
            The path to the Excel file used to construct the DataFrame
        """

        # TODO: should we add more in-depth checks here?

        check_column_names(df=geo_strata_df, expected_names=self.geo_strata_cols, path_for_df=df_path)

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

        if self.survey.params['geo_strata_sheetname'] in ['stratification1', 'INPFC']:

            # check existence of the file
            file_path = self.survey.params['data_root_dir'] / self.survey.params['geo_strata_filename']
            check_existence_of_file(file_path)

            # read and check geographic stratification file
            geo_strata_df = pd.read_excel(file_path, sheet_name=self.survey.params['geo_strata_sheetname'])
            self._check_geo_strata_df(geo_strata_df, file_path)

            # extract only those columns that are necessary
            geo_strata_df = geo_strata_df[['stratum_num', 'Latitude (upper limit)']].copy()

            # set data types of dataframe
            geo_strata_df = geo_strata_df.astype({'stratum_num': int,
                                                  'Latitude (upper limit)': np.float64})

            self.survey.geo_strata_df = geo_strata_df

        else:
            raise NotImplementedError(f"geo_strata_filename has unknown sheet name!")
