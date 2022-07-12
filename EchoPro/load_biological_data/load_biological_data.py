import numpy as np
import pandas as pd


class LoadBioData:  # TODO: Does it make sense for this to be a class?
    """
    This class loads and checks those files
    associated with the biological data.

    Parameters
    ----------
    epro : EchoPro object
        An initialized EchoPro object. Note that any change to
        self.epro will also change this object.
    """

    def __init__(self, epro=None):

        self.epro = epro

        # expected columns for length Dataframe
        self.len_cols = {'Haul', 'Species_Code', 'Sex', 'Length', 'Frequency'}

        # expected columns for specimen Dataframe
        self.spec_cols = {'Haul', 'Species_Code', 'Sex', 'Length', 'Weight',
                          'Age'}

        self._load_length_data()
        self._load_specimen_data()

    def _check_length_df(self, len_df: pd.DataFrame):
        """
        Ensures that the appropriate columns are
        contained in the length Dataframe.

        TODO: should we add more in-depth checks here?
        """

        if not set(len_df.columns).intersection(self.len_cols):
            raise NameError("Length dataframe does not contain all expected columns!")

    def _check_specimen_df(self, spec_df: pd.DataFrame):
        """
        Ensures that the appropriate columns are
        contained in the specimen Dataframe.

        TODO: should we add more in-depth checks here?
        """

        if not set(spec_df.columns).intersection(self.spec_cols):
            raise NameError("Specimen dataframe does not contain all expected columns!")

    def _process_length_data_df(self, df: pd.DataFrame, haul_num_offset: int):
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
        df = df.loc[df['Species_Code'] == self.epro.params['species_code_ID']]

        # Apply haul offset
        df['Haul'] = df['Haul'] + haul_num_offset

        if self.epro.params['exclude_age1'] is False:
            raise NotImplementedError("Including age 1 data has not been implemented!")

        # remove species code column
        df.drop(columns=['Species_Code'], inplace=True)

        df.set_index('Haul', inplace=True)

        return df

    def _process_specimen_data(self, df: pd.DataFrame, haul_num_offset: int):
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
        df = df.loc[df['Species_Code'] == self.epro.params['species_code_ID']]

        # Apply haul_num_offset
        df['Haul'] = df['Haul'] + haul_num_offset

        if self.epro.params['exclude_age1'] is False:
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

    def _load_length_data(self):
        """
        Loads and prepares data associated with a station
        that records the length and sex of the animal.
        Additionally, it sets epro.length_df using the
        final processed dataframe.
        """

        if self.epro.params['source'] == 3:

            # read in and check US and Canada Excel files
            df_us = pd.read_excel(self.epro.params['data_root_dir'] + self.epro.params['length_US_filename'],
                                  sheet_name=self.epro.params['length_US_sheet'])
            self._check_length_df(df_us)

            df_can = pd.read_excel(self.epro.params['data_root_dir'] + self.epro.params['length_CAN_filename'],
                                   sheet_name=self.epro.params['length_CAN_sheet'])
            self._check_length_df(df_can)

            # process US and Canada dataframes
            length_us_df = self._process_length_data_df(df_us, 0)
            length_can_df = self._process_length_data_df(df_can, self.epro.params['CAN_haul_offset'])

            # Construct full length dataframe from US and Canada sections
            self.epro.length_df = pd.concat([length_us_df, length_can_df])

        else:
            raise NotImplementedError(f"Source of {self.epro.params['source']} not implemented yet.")

    def _load_specimen_data(self):
        """
        Loads and prepares data associated with a station
        that records the length, weight, age, and sex of
        the animal. Additionally, it sets epro.specimen_df
        using the final processed dataframe.
        """

        if self.epro.params['source'] == 3:

            # read in and check US and Canada Excel files
            specimen_us_df = pd.read_excel(self.epro.params['data_root_dir'] + self.epro.params['specimen_US_filename'],
                                           sheet_name=self.epro.params['specimen_US_sheet'])
            self._check_specimen_df(specimen_us_df)

            specimen_can_df = pd.read_excel(self.epro.params['data_root_dir'] + self.epro.params['specimen_CAN_filename']
                                            , sheet_name=self.epro.params['specimen_CAN_sheet'])
            self._check_specimen_df(specimen_can_df)

            # process US and Canada dataframes
            specimen_us_df = self._process_specimen_data(specimen_us_df, 0)
            specimen_can_df = self._process_specimen_data(specimen_can_df, self.epro.params['CAN_haul_offset'])

            # Construct full specimen dataframe from US and Canada sections
            self.epro.specimen_df = pd.concat([specimen_us_df, specimen_can_df])

        else:
            raise NotImplementedError(f"Source of {self.epro.params['source']} not implemented yet.")
