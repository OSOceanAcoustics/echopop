import numpy as np
import pandas as pd
import xarray as xr
import warnings
import sys


class LoadBioData:
    """
    A Class that loads and processes all
    biological data

    Parameters
    ----------
    EPro : EchoPro object
        An initialized EchoPro object. Note that any change to
        self.EPro will also change this object.
    """

    def __init__(self, EPro = None):

        self.EPro = EPro

        self.__load_biological_data()

    def __process_length_data_ds_wu_jung(self, df, haul_num_offset):
        """
        Process and turn the length data file into a xarray Dataset, where the frequency
        column has not been expanded. The Dataset has coordinates Haul, Length, and Sex.
        Parameters
        ----------
        df : Pandas Dataframe
            Dataframe holding the length data
        haul_num_offset : int
            The value that should be added to the haul index to differentiate it from other ships

        Returns
        -------
        xarray Dataset

        Notes
        -----
        The produced Dataset is a compressed version of the true data.
        """

        # obtaining those columns that are required
        df = df[['Haul', 'Species_Code', 'Sex', 'Length', 'Frequency']].copy()

        # extract target species
        df = df.loc[df['Species_Code'] == self.EPro.params['species_code_ID']]

        # set data types of dataframe
        df = df.astype({'Haul': int, 'Species_Code': int, 'Sex': int, 'Length': np.float64, 'Frequency': int})

        # Apply haul_num_offset
        df['Haul'] = df['Haul'] + haul_num_offset

        if self.EPro.params['exclude_age1'] is False:
            raise NotImplementedError("Including age 1 data has not been implemented!")

        df.drop(columns=['Species_Code'], inplace=True)

        df.set_index(['Haul', 'Length', 'Sex'], inplace=True)

        ds = df.groupby(level=[0, 1, 2]).sum().to_xarray()

        ds['n'] = ds.Frequency.sum(['Length', 'Sex'])

        ds['meanlen'] = (ds.Length * ds.Frequency).sum(['Length', 'Sex']) / ds.n

        ds['stdlen'] = np.sqrt((((np.abs(ds.Length - ds.meanlen) ** 2) * ds.Frequency).sum(['Length', 'Sex'])) / (ds.n - 1.0))

        TS0 = 20.0 * np.log10(ds.Length) - 68.0

        ds['TS_lin'] = 10.0 * np.log10((10.0 ** (TS0 / 10.0) * ds.Frequency).sum(['Length', 'Sex']) / (ds.n))

        ds['TS_log'] = (TS0 * ds.Frequency).sum(['Length', 'Sex']) / (ds.n)

        ds['TS_sd'] = np.sqrt((((np.abs(TS0 - ds.TS_log) ** 2) * ds.Frequency).sum(['Length', 'Sex'])) / (ds.n - 1.0))

        ds['nM'] = ds.sel(Sex=1).Frequency.sum('Length')
        ds['nF'] = ds.sel(Sex=2).Frequency.sum('Length')

        return ds

    def __process_length_data_df(self, df, haul_num_offset):
        """
        Process the length data file, which includes expanding the frequencies column.
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

        # extract target species
        df = df.loc[df['Species_Code'] == self.EPro.params['species_code_ID']]

        # set data types of dataframe
        df = df.astype({'Haul': int, 'Species_Code': int, 'Sex': int, 'Length': np.float64, 'Frequency': np.float64})

        # Apply haul_num_offset
        df['Haul'] = df['Haul'] + haul_num_offset

        if self.EPro.params['exclude_age1'] is False:
            raise NotImplementedError("Including age 1 data has not been implemented!")

        df.drop(columns=['Species_Code'], inplace=True)

        # expand length and sex columns
        df['Frequency'] = df['Frequency'].apply(lambda x: np.ones(int(x)))
        df['Length'] = df['Length'] * df['Frequency']
        df['Sex'] = df['Sex'] * df['Frequency']
        df.drop(columns='Frequency', inplace=True)

        # Group by Haul and calculate necessary statistics
        df = df.groupby('Haul').agg(lambda x: np.concatenate(x.values))
        df['TS0'] = 20.0 * df['Length'].apply(lambda x: np.log10(x)) - 68.0
        df = pd.concat([df, df.apply(lambda x: pd.Series([np.where(x.Sex == 1)[0], np.where(x.Sex == 2)[0],
                                                          len(x.Length), np.nanmean(x.Length),
                                                          np.nanstd(x.Length, ddof=1),
                                                          10.0 * np.log10(np.nanmean(10.0 ** (x.TS0 / 10.0))),
                                                          np.std(x.TS0, ddof=1), np.mean(x.TS0)],
                                                         index=['Male_ind', 'Female_ind', 'n', 'meanlen', 'stdlen',
                                                                'TS_lin', 'TS_sd', 'TS_log']),
                                     result_type='expand', axis=1)], axis=1)
        df.drop(columns='TS0', inplace=True)

        # get number of female and males
        df[['nM', 'nF']] = df.apply(lambda x: pd.Series([len(x.Male_ind), len(x.Female_ind)]), result_type='expand',
                                    axis=1)

        return df

    def __process_specimen_data(self, df, haul_num_offset):
        """
        Parameters
        ----------
        df : Pandas Dataframe
            Dataframe holding the specimen data
        haul_num_offset : int
            # TODO: what does this refer to? Should we account for this in the code?
        Returns
        -------
        Processed Dataframe
        """

        df = df[['Haul', 'Species_Code', 'Sex', 'Length', 'Weight', 'Specimen_Number', 'Age']].copy()

        # extract target species
        df = df.loc[df['Species_Code'] == self.EPro.params['species_code_ID']]

        # set data types of dataframe
        df = df.astype({'Haul': int, 'Species_Code': int, 'Sex': int, 'Length': np.float64, 'Weight': np.float64,
                        'Specimen_Number': np.float64, 'Age': np.float64})

        # Apply haul_num_offset
        df['Haul'] = df['Haul'] + haul_num_offset

        if self.EPro.params['exclude_age1'] is False:
            raise NotImplementedError("Including age 1 data has not been implemented!")

        df.drop(columns=['Species_Code'], inplace=True)

        df.set_index('Haul', inplace=True)
        df.sort_index(inplace=True)

        if len(df['Age']) - df['Age'].isna().sum() < 0.1 * len(df['Age']):
            raise RuntimeWarning('Aged data are less than 10%!\n')

        return df

    def __load_length_data(self):

        if self.EPro.params['source'] == 3:

            df_us = pd.read_excel(self.EPro.params['data_root_dir'] + self.EPro.params['filename_length_US'],
                                  sheet_name='biodata_length')
            df_can = pd.read_excel(self.EPro.params['data_root_dir'] + self.EPro.params['filename_length_CAN'],
                                   sheet_name='biodata_length_CAN')

            # TODO: Pick an option and stick with it!!
            # Option 1, Dataframe with array elements
            length_us_df = self.__process_length_data_df(df_us, 0)
            length_can_df = self.__process_length_data_df(df_can, self.EPro.params['CAN_haul_offset'])
            self.EPro.length_df = pd.concat([length_us_df, length_can_df])

            # Option 3, Dataset with Haul, Length, and Sex as coordinates
            length_us_ds_wu_jung = self.__process_length_data_ds_wu_jung(df_us, 0)
            length_can_ds_wu_jung = self.__process_length_data_ds_wu_jung(df_can, self.EPro.params['CAN_haul_offset'])
            self.EPro.length_ds = length_us_ds_wu_jung.merge(length_can_ds_wu_jung)

        else:
            raise NotImplementedError(f"Source of {self.EPro.params['source']} not implemented yet.")

    def __load_specimen_data(self):

        if self.EPro.params['source'] == 3:

            specimen_us_df = pd.read_excel(self.EPro.params['data_root_dir'] + self.EPro.params['filename_specimen_US'],
                                           sheet_name='biodata_specimen')
            specimen_us_df = self.__process_specimen_data(specimen_us_df, 0)

            specimen_can_df = pd.read_excel(self.EPro.params['data_root_dir'] + self.EPro.params['filename_specimen_CAN']
                                            , sheet_name='biodata_specimen_CAN')
            specimen_can_df = self.__process_specimen_data(specimen_can_df, self.EPro.params['CAN_haul_offset'])

            self.EPro.specimen_df = pd.concat([specimen_us_df, specimen_can_df])

        else:
            raise NotImplementedError(f"Source of {self.EPro.params['source']} not implemented yet.")

    def __load_biological_data(self):
        """
        A function that loads the appropriate data based on the
        source provided.

        Returns
        -------
        self.length_df : Dataframe
            A dataframe containing the length information to
            be used downstream. If source=3 it contains US and Canada.
        self.specimen_df : Dataframe
            A dataframe containing the specimen information to be used
            downstream. If source=3 it contains US and Canada.
        """

        self.__load_length_data()

        self.__load_specimen_data()
