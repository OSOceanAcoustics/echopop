import numpy as np
import pandas as pd
from pathlib import Path
from ..utils.input_checks import check_column_names, check_existence_of_file


class LoadBioData:  # TODO: Does it make sense for this to be a class?
    """
    This class loads and checks those files
    associated with the biological data.

    Parameters
    ----------
    survey : Survey
        An initialized Survey object. Note that any change to
        self.survey will also change this object.
    """

    def __init__(self, survey=None):

        self.survey = survey

        # expected columns for length Dataframe
        self.len_cols = {'Haul', 'Species_Code', 'Sex', 'Length', 'Frequency'}

        # expected columns for specimen Dataframe
        self.spec_cols = {'Haul', 'Species_Code', 'Sex', 'Length', 'Weight',
                          'Age'}

        # expected columns for gear Dataframe
        self.gear_cols = {'Haul', 'Transect'}

        self._load_length_data()
        self._load_specimen_data()
        self._load_gear_data()

    def _check_length_df(self, len_df: pd.DataFrame, df_path: Path) -> None:
        """
        Ensures that the appropriate columns are
        contained in the length Dataframe.

        Parameters
        ----------
        len_df: pd.DataFrame
            The constructed Length DataFrame
        df_path: Path
            The path to the Excel file used to construct the DataFrame
        """

        # TODO: should we add more in-depth checks here?

        check_column_names(df=len_df, expected_names=self.len_cols, path_for_df=df_path)

    def _check_specimen_df(self, spec_df: pd.DataFrame, df_path: Path) -> None:
        """
        Ensures that the appropriate columns are
        contained in the specimen Dataframe.

        Parameters
        ----------
        spec_df: pd.DataFrame
            The constructed Specimen DataFrame
        df_path: Path
            The path to the Excel file used to construct the DataFrame
        """

        # TODO: should we add more in-depth checks here?

        check_column_names(df=spec_df, expected_names=self.spec_cols, path_for_df=df_path)

    def _check_gear_df(self, gear_df: pd.DataFrame, df_path: Path) -> None:
        """
        Ensures that the appropriate columns are
        contained in the gear Dataframe.

        Parameters
        ----------
        gear_df: pd.DataFrame
            The constructed Gear DataFrame
        df_path: Path
            The path to the Excel file used to construct the DataFrame
        """

        # TODO: should we add more in-depth checks here?

        check_column_names(df=gear_df, expected_names=self.gear_cols, path_for_df=df_path)

    def _process_length_data_df(self, df: pd.DataFrame,
                                haul_num_offset: int) -> pd.DataFrame:
        """
        Processes the length dataframe by:
        * Obtaining the required columns from the dataframe
        * Extracting only the target species
        * Setting the data type of each column
        * Applying a haul offset, if necessary
        * Replacing the length and sex columns with an array
        of length frequency and dropping the frequency column
        * Setting the index required for downstream processes

        Parameters
        ----------
        df : Pandas Dataframe
            Dataframe holding the length data
        haul_num_offset : int
            The offset that should be applied to the Haul column

        Returns
        -------
        Processed Dataframe
        """

        # obtaining those columns that are required
        df = df[['Haul', 'Species_Code', 'Sex', 'Length', 'Frequency']].copy()

        # set data types of dataframe
        df = df.astype({'Haul': int, 'Species_Code': int, 'Sex': int,
                        'Length': np.float64, 'Frequency': np.float64})

        # extract target species
        df = df.loc[df['Species_Code'] == self.survey.params['species_code_ID']]

        # Apply haul offset
        df['Haul'] = df['Haul'] + haul_num_offset

        if self.survey.params['exclude_age1'] is False:
            raise NotImplementedError("Including age 1 data has not been implemented!")

        # remove species code column
        df.drop(columns=['Species_Code'], inplace=True)

        df.set_index('Haul', inplace=True)

        return df

    def _process_specimen_data(self, df: pd.DataFrame,
                               haul_num_offset: int) -> pd.DataFrame:
        """
        Processes the specimen dataframe by:
        * Obtaining the required columns from the dataframe
        * Extracting only the target species
        * Setting the data type of each column
        * Applying a haul offset, if necessary
        * Setting the index required for downstream processes

        Parameters
        ----------
        df : Pandas Dataframe
            Dataframe holding the specimen data
        haul_num_offset : int

        Returns
        -------
        Processed Dataframe
        """

        # obtaining those columns that are required
        df = df[['Haul', 'Species_Code', 'Sex', 'Length', 'Weight', 'Age']].copy()

        # set data types of dataframe
        df = df.astype({'Haul': int, 'Species_Code': int, 'Sex': int,
                        'Length': np.float64, 'Weight': np.float64,
                        'Age': np.float64})

        # extract target species
        df = df.loc[df['Species_Code'] == self.survey.params['species_code_ID']]

        # Apply haul_num_offset
        df['Haul'] = df['Haul'] + haul_num_offset

        if self.survey.params['exclude_age1'] is False:
            raise NotImplementedError("Including age 1 data has not been implemented!")

        # remove species code column
        df.drop(columns=['Species_Code'], inplace=True)

        # set and organize index
        df.set_index('Haul', inplace=True)
        df.sort_index(inplace=True)

        # perform check on data
        if len(df['Age']) - df['Age'].isna().sum() < 0.1 * len(df['Age']):
            raise RuntimeWarning('Aged data are less than 10%!\n')

        return df

    def _load_length_data(self) -> None:
        """
        Loads and prepares data associated with a station
        that records the length and sex of the animal.
        Additionally, it sets survey.length_df using the
        final processed dataframe.
        """

        if self.survey.params['source'] == 3:

            # check existence of the files
            file_path_us = self.survey.params['data_root_dir'] / self.survey.params['length_US_filename']
            file_path_can = self.survey.params['data_root_dir'] / self.survey.params['length_CAN_filename']
            check_existence_of_file(file_path_us)
            check_existence_of_file(file_path_can)

            # read in and check US and Canada Excel files
            df_us = pd.read_excel(file_path_us, sheet_name=self.survey.params['length_US_sheet'])
            self._check_length_df(df_us, file_path_us)

            df_can = pd.read_excel(file_path_can, sheet_name=self.survey.params['length_CAN_sheet'])
            self._check_length_df(df_can, file_path_can)

            # process US and Canada dataframes
            length_us_df = self._process_length_data_df(df_us, 0)
            length_can_df = self._process_length_data_df(df_can, self.survey.params['CAN_haul_offset'])

            # Construct full length dataframe from US and Canada sections
            self.survey.length_df = pd.concat([length_us_df, length_can_df])

        else:
            raise NotImplementedError(f"Source of {self.survey.params['source']} not implemented yet.")

    def _load_specimen_data(self) -> None:
        """
        Loads and prepares data associated with a station
        that records the length, weight, age, and sex of
        the animal. Additionally, it sets survey.specimen_df
        using the final processed dataframe.
        """

        if self.survey.params['source'] == 3:

            # check existence of the files
            file_path_us = self.survey.params['data_root_dir'] / self.survey.params['specimen_US_filename']
            file_path_can = self.survey.params['data_root_dir'] / self.survey.params['specimen_CAN_filename']
            check_existence_of_file(file_path_us)
            check_existence_of_file(file_path_can)

            # read in and check US and Canada Excel files
            specimen_us_df = pd.read_excel(file_path_us, sheet_name=self.survey.params['specimen_US_sheet'])
            self._check_specimen_df(specimen_us_df, file_path_us)

            specimen_can_df = pd.read_excel(file_path_can, sheet_name=self.survey.params['specimen_CAN_sheet'])
            self._check_specimen_df(specimen_can_df, file_path_can)

            # process US and Canada dataframes
            specimen_us_df = self._process_specimen_data(specimen_us_df, 0)
            specimen_can_df = self._process_specimen_data(specimen_can_df, self.survey.params['CAN_haul_offset'])

            # Construct full specimen dataframe from US and Canada sections
            self.survey.specimen_df = pd.concat([specimen_us_df, specimen_can_df])

        else:
            raise NotImplementedError(f"Source of {self.survey.params['source']} not implemented yet.")

    def _process_gear_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Processes the gear data by
        * selecting the haul and transect columns
        * ensuring the dataframe has the appropriate data types
        * setting the ``Haul`` column as the Dataframe index
        * sorting the ``Haul`` index in ascending order

        Parameters
        ----------
        df : pd. Dataframe
            Dataframe holding the gear data

        Returns
        -------
        pd.DataFrame
            Processed gear Dataframe
        """

        # obtain those columns necessary for core downstream processes
        df = df[['Haul', 'Transect']].copy()

        # set data types of dataframe
        df = df.astype({'Haul': int, 'Transect': np.float64})

        if self.survey.params['exclude_age1'] is False:
            raise NotImplementedError("Including age 1 data has not been implemented!")

        # set Haul as index and sort it
        df.set_index('Haul', inplace=True)
        df.sort_index(inplace=True)

        return df

    def _load_gear_data(self) -> None:
        """
        Loads and prepares the gear data ``Haul`` and ``Transects``. Additionally,
        it sets survey.gear_df using the final processed dataframe.

        Notes
        -----
        This data is currently only being used as a mapping between hauls and
        transects, which is necessary for obtaining subsets of data during
        transect selection.
        """
        # TODO: replace the gear file with a more purposeful mapping between hauls and transects
        # TODO: make sheet_name part of configuration file, if we need it

        if self.survey.params['source'] == 3:

            # check existence of the files
            file_path_us = self.survey.params['data_root_dir'] / self.survey.params['filename_gear_US']
            file_path_can = self.survey.params['data_root_dir'] / self.survey.params['filename_gear_CAN']
            check_existence_of_file(file_path_us)
            check_existence_of_file(file_path_can)

            # read in, check, and process the US gear file
            gear_us_df = pd.read_excel(file_path_us, sheet_name='biodata_gear')
            self._check_gear_df(gear_us_df, file_path_us)
            gear_us_df = self._process_gear_data(gear_us_df)

            # read in, check, and process the Canada gear file
            gear_can_df = pd.read_excel(file_path_can, sheet_name='biodata_gear_CAN')
            self._check_gear_df(gear_can_df, file_path_can)
            gear_can_df = self._process_gear_data(gear_can_df)

            # transect offset is set to zero since we will not have overlapping transects
            # TODO: Is this a necessary variable? If so, we should make it an input to the function.
            CAN_Transect_offset = 0

            # add Canada transect and haul offset
            gear_can_df.index = gear_can_df.index + self.survey.params['CAN_haul_offset']
            gear_can_df['Transect'] = gear_can_df['Transect'] + CAN_Transect_offset

            # combine US & CAN trawl files
            self.survey.gear_df = pd.concat([gear_us_df, gear_can_df])

        else:
            raise NotImplementedError(f"Source of {self.survey.params['source']} not implemented yet.")
