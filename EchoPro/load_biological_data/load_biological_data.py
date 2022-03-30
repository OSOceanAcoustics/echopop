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
    params : dict
        The parameter dictionary constructed in the EchoPro class. Note that any
        modification to this dictionary made within this class will be also
        done to params in the EchoPro class.
    """

    def __init__(self, EPro = None):

        # self.params = params

        self.EPro = EPro

        print("Loading biological data ...")

        self.__load_biological_data()

    def __process_length_data_ds(self, df, haul_num_offset):
        """
        Process and turn the length data file into a xarray Dataset, where the frequency
        column has not been expanded. This one has 'haul' and 'haul_index' as coordinates.
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
        df = df.astype({'Haul': int, 'Species_Code': int, 'Sex': int, 'Length': np.float64, 'Frequency': np.float64})

        # check to make sure that the Sex column is composed of 1's and 2's
        unknown_sex_values = set(df.Sex.unique()) - {1, 2}
        if unknown_sex_values:
            warnings.warn(
                f"The Sex column contains values {unknown_sex_values}, converting values greater than 2 to 2 and values less than 2 to 1")

            df.loc[df.Sex > 2, 'Sex'] = int(2)
            df.loc[df.Sex < 1, 'Sex'] = int(1)

        # Apply haul_num_offset
        df['Haul'] = df['Haul'] + haul_num_offset

        if self.EPro.params['exclude_age1'] is False:
            raise NotImplementedError("Including age 1 data has not been implemented!")

        df.drop(columns=['Species_Code'], inplace=True)

        df.set_index('Haul', inplace=True)
        max_haul_size = df.index.value_counts().max()
        series_haul = df.groupby('Haul').apply(np.array)

        # pad the numpy arrays of each haul to match the maximum haul size
        series_haul_pad = series_haul.apply(
            lambda x: np.vstack((x, np.nan * np.ones((max_haul_size - x.shape[0], 3))))
            if x.shape[0] != max_haul_size else x)

        # collect all hauls into a 2D array
        np_haul_pad = series_haul_pad.agg(lambda x: np.vstack(x.values))

        ds = xr.Dataset(
            data_vars={
                'sex': (['haul', 'haul_index'], np_haul_pad[:, 0].reshape((series_haul_pad.shape[0], max_haul_size))),
                'bio_length': (
                ['haul', 'haul_index'], np_haul_pad[:, 1].reshape((series_haul_pad.shape[0], max_haul_size))),
                'frequency': (
                ['haul', 'haul_index'], np_haul_pad[:, 2].reshape((series_haul_pad.shape[0], max_haul_size)))
            },
            coords={
                'haul': (['haul'], series_haul_pad.index),
                'haul_index': (['haul_index'], range(max_haul_size))
            })

        # add variables that are helpful for downstream processes
        ds['n'] = ds.frequency.sum('haul_index', skipna=True)
        ds['meanlen'] = (ds.bio_length * ds.frequency).sum('haul_index', skipna=True) / (ds.n)
        ds['stdlen'] = np.sqrt(
            (((abs(ds.bio_length - ds.meanlen) ** 2) * ds.frequency).sum('haul_index', skipna=True)) / (ds.n - 1.0))

        TS0 = 20.0 * np.log10(ds.bio_length) - 68.0
        ds['TS_lin'] = 10.0 * np.log10((10.0 ** (TS0 / 10.0) * ds.frequency).sum('haul_index', skipna=True) / (ds.n))
        ds['TS_log'] = (TS0 * ds.frequency).sum('haul_index', skipna=True) / (ds.n)
        ds['TS_sd'] = np.sqrt(
            (((abs(TS0 - ds.TS_log) ** 2) * ds.frequency).sum('haul_index', skipna=True)) / (ds.n - 1.0))

        ds['nM'] = (ds[['sex', 'frequency']].where(ds.sex == 1.0, drop=True)).frequency.sum('haul_index', skipna=True)
        ds['nF'] = (ds[['sex', 'frequency']].where(ds.sex == 2.0, drop=True)).frequency.sum('haul_index', skipna=True)

        return ds

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

        # check to make sure that the Sex column is composed of 1's and 2's
        unknown_sex_values = set(df.Sex.unique()) - {1, 2}
        if unknown_sex_values:
            warnings.warn(
                f"The Sex column contains values {unknown_sex_values}, converting values greater than 2 to 2 and values less than 2 to 1")

            df.loc[df.Sex > 2, 'Sex'] = int(2)
            df.loc[df.Sex < 1, 'Sex'] = int(1)

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

        # check to make sure that the Sex column is composed of 1's and 2's
        unknown_sex_values = set(df.Sex.unique()) - {1, 2}
        if unknown_sex_values:
            warnings.warn(
                f"The Sex column contains values {unknown_sex_values}, converting values greater than 2 to 2 and values less than 2 to 1")

            df.loc[df.Sex > 2, 'Sex'] = int(2)
            df.loc[df.Sex < 1, 'Sex'] = int(1)

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

    @staticmethod
    def __process_catch_data_df(df):
        """
        Parameters
        ----------
        df : Pandas Dataframe
            Dataframe holding the catch data
        Returns
        -------
        Processed Dataframe
        """

        # obtaining those columns that are required
        df = df[['Haul', 'Species_Code', 'Species_Name', 'Number_In_Haul', 'Weight_In_Haul']].copy()

        # set data types of dataframe
        df = df.astype({'Haul': int, 'Species_Code': int, 'Species_Name': str, 'Number_In_Haul': np.float64,
                        'Weight_In_Haul': np.float64})

        df.set_index('Haul', inplace=True)

        df.sort_index(inplace=True)

        return df

    @staticmethod
    def __process_catch_data_ds(df):
        """
        Parameters
        ----------
        df : Pandas Dataframe
            Dataframe holding the catch data
        Returns
        -------
        Processed xarray Dataset
        """

        # obtaining those columns that are required
        df = df[['Haul', 'Species_Code', 'Species_Name', 'Number_In_Haul', 'Weight_In_Haul']].copy()

        # set data types of dataframe
        df = df.astype({'Haul': int, 'Species_Code': int, 'Species_Name': str, 'Number_In_Haul': np.float64,
                        'Weight_In_Haul': np.float64})

        df.set_index('Haul', inplace=True)

        df.sort_index(inplace=True)

        max_haul_size = df.index.value_counts().max()
        series_haul = df.groupby('Haul').apply(np.array)

        # pad the numpy arrays of each haul to match the maximum haul size
        series_haul_pad = series_haul.apply(
            lambda x: np.vstack((x, np.nan * np.ones((max_haul_size - x.shape[0], 4)))) if x.shape[
                                                                                               0] != max_haul_size else x)

        # collect all hauls into a 2D array
        np_haul_pad = series_haul_pad.agg(lambda x: np.vstack(x.values))

        ds = xr.Dataset(
            data_vars={
                'species_code': (
                ['haul', 'haul_index'], np_haul_pad[:, 0].reshape((series_haul_pad.shape[0], max_haul_size))),
                'species_name': (
                ['haul', 'haul_index'], np_haul_pad[:, 1].reshape((series_haul_pad.shape[0], max_haul_size))),
                'number_in_haul': (
                ['haul', 'haul_index'], np_haul_pad[:, 2].reshape((series_haul_pad.shape[0], max_haul_size))),
                'weight_in_haul': (
                ['haul', 'haul_index'], np_haul_pad[:, 3].reshape((series_haul_pad.shape[0], max_haul_size))),
                'no_null_length': (['haul'], df.index.value_counts().sort_index().to_numpy())
            },
            coords={
                'haul': (['haul'], series_haul_pad.index),
                'haul_index': (['haul_index'], range(max_haul_size))
            })

        return ds

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

        df = df[['Haul', 'Haul_Type', 'Performance_Code', 'Duration', 'Distance_Fished', 'Stratum', 'EQ_Latitude',
                 'EQ_Longitude', 'Average_Bottom_Depth', 'Vessel_Log_Start', 'Vessel_Log_Stop']].copy()

        if df.dtypes['Average_Bottom_Depth'] == object:
            # remove any instances of >, <, or m from the column Average_Bottom_Depth
            df['Average_Bottom_Depth'] = df['Average_Bottom_Depth'].apply(
                lambda x: x.strip('><m') if type(x) == str else x).astype(float)

        df['Duration'] = pd.to_timedelta(df['Duration'])  # change datatype from string to timedelta

        # set data types of dataframe
        df = df.astype({'Haul': int, 'Haul_Type': int, 'Performance_Code': np.float64, 'Distance_Fished': np.float64,
                        'Stratum': int, 'EQ_Latitude': np.float64, 'EQ_Longitude': np.float64,
                        'Average_Bottom_Depth': np.float64, 'Vessel_Log_Start': np.float64,
                        'Vessel_Log_Stop': np.float64})

        df.set_index('Haul', inplace=True)
        df.sort_index(inplace=True)

        if self.EPro.params['exclude_age1'] is False:
            raise NotImplementedError("Including age 1 data has not been implemented!")

        if self.EPro.params['hemisphere'][0] == 'N':
            df['EQ_Latitude'] = df['EQ_Latitude'].abs()
        elif self.EPro.params['hemisphere'][0] == 'S':
            df['EQ_Latitude'] = -1.0 * (df['EQ_Latitude'].abs())
        else:
            raise ValueError('Wrong N/S Hemisphere provided! \n')

        if self.EPro.params['hemisphere'][1] == 'W':
            df['EQ_Longitude'] = -1.0 * (df['EQ_Longitude'].abs())
        elif self.EPro.params['hemisphere'][1] == 'E':
            df['EQ_Longitude'] = df['EQ_Longitude'].abs()
        else:
            raise ValueError('Wrong E/W Hemisphere provided! \n')

        return df

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

        df = df[['Haul', 'Average_Footrope_Depth', 'Surface_Temperature', 'Gear_Temperature', 'Average_Wireout',
                 'Transect', 'Net_Height']].copy()

        # set data types of dataframe
        df = df.astype({'Haul': int, 'Average_Footrope_Depth': np.float64, 'Surface_Temperature': np.float64,
                        'Gear_Temperature': np.float64, 'Average_Wireout': np.float64, 'Transect': np.float64})

        if self.EPro.params['exclude_age1'] is False:
            raise NotImplementedError("Including age 1 data has not been implemented!")

        df['Net_Height'] = df['Net_Height'].apply(lambda x: np.nan if type(x) == str else x).astype(float)

        df.set_index('Haul', inplace=True)
        df.sort_index(inplace=True)

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

        # check to make sure that the Sex column is composed of 1's and 2's
        unknown_sex_values = set(df.Sex.unique()) - {1, 2}
        if unknown_sex_values:
            warnings.warn(
                f"The Sex column contains values {unknown_sex_values}, converting values greater than 2 to 2 and values less than 2 to 1")

            df.loc[df.Sex > 2, 'Sex'] = int(2)
            df.loc[df.Sex < 1, 'Sex'] = int(1)

        # Apply haul_num_offset
        df['Haul'] = df['Haul'] + haul_num_offset

        if self.EPro.params['exclude_age1'] is False:
            raise NotImplementedError("Including age 1 data has not been implemented!")

        df.drop(columns=['Species_Code'], inplace=True)

        if self.EPro.params['age_data_status'] == 2:
            raise NotImplementedError("age_data_status = 2 has not been implemented")
        elif self.EPro.params['age_data_status'] != 1:
            raise NotImplementedError(
                f"age_data_status = {self.EPro.params['age_data_status']} has not been implemented")

        df.set_index('Haul', inplace=True)
        df.sort_index(inplace=True)

        if len(df['Age']) - df['Age'].isna().sum() < 0.1 * len(df['Age']):
            raise RuntimeWarning('Aged data are less than 10%!\n')

        return df

    def __load_catch_data(self):

        if self.EPro.params['source'] == 3:

            if self.EPro.params['database_type'] == 'Oracle':
                catch_us_df = pd.read_excel(self.EPro.params['data_root_dir'] + self.EPro.params['filename_catch_US'])
                catch_us_df = self.__process_catch_data_df(catch_us_df)

                catch_can_df = pd.read_excel(self.EPro.params['data_root_dir'] + self.EPro.params['filename_catch_CAN'])
                catch_can_df = self.__process_catch_data_df(catch_can_df)

                # combine US & CAN catch files
                catch_can_df.index = catch_can_df.index + self.EPro.params["CAN_haul_offset"]  # add haul_offset
                self.EPro.catch_df = pd.concat([catch_us_df, catch_can_df])

            else:
                raise NotImplementedError("Loading data from a non-Oracle database has not been implemented!")
        else:

            raise NotImplementedError(f"Source of {self.EPro.params['source']} not implemented yet.")

    def __load_length_data(self):

        if self.EPro.params['source'] == 3:

            df_us = pd.read_excel(self.EPro.params['data_root_dir'] + self.EPro.params['filename_length_US'])
            df_can = pd.read_excel(self.EPro.params['data_root_dir'] + self.EPro.params['filename_length_CAN'])

            # Option 1, Dataframe with array elements
            length_us_df = self.__process_length_data_df(df_us, 0)
            length_can_df = self.__process_length_data_df(df_can, self.EPro.params['CAN_haul_offset'])
            self.EPro.length_df = pd.concat([length_us_df, length_can_df])

            # Option 2, Dataset with haul and haul_index as coordinates
            length_us_ds_bran = self.__process_length_data_ds(df_us, 0)
            length_can_ds_bran = self.__process_length_data_ds(df_can, self.EPro.params['CAN_haul_offset'])
            self.EPro.length_ds_bran = length_us_ds_bran.merge(length_can_ds_bran)

            # Option 3, Dataset with Haul, Length, and Sex as coordinates
            length_us_ds_wu_jung = self.__process_length_data_ds_wu_jung(df_us, 0)
            length_can_ds_wu_jung = self.__process_length_data_ds_wu_jung(df_can, self.EPro.params['CAN_haul_offset'])
            self.EPro.length_ds_wu_jung = length_us_ds_wu_jung.merge(length_can_ds_wu_jung)

        else:
            raise NotImplementedError(f"Source of {self.EPro.params['source']} not implemented yet.")

    def __load_trawl_data(self):

        if self.EPro.params['source'] == 3:

            if self.EPro.params['filename_trawl_US']:
                trawl_us_df = pd.read_excel(self.EPro.params['data_root_dir'] + self.EPro.params['filename_trawl_US'],
                                            converters={'Duration': str})
                trawl_us_df = self.__process_trawl_data(trawl_us_df)
            else:
                trawl_us_df = None

            if self.EPro.params['filename_trawl_CAN']:
                trawl_can_df = pd.read_excel(self.EPro.params['data_root_dir'] + self.EPro.params['filename_trawl_CAN'],
                                             converters={'Duration': str})
                trawl_can_df = self.__process_trawl_data(trawl_can_df)

                trawl_can_df.index = trawl_can_df.index + self.EPro.params['CAN_haul_offset']  # add haul_offset

                # change lon to -lon
                trawl_can_df['EQ_Longitude'] = trawl_can_df['EQ_Longitude'].apply(lambda x: -1.0 * x if x > 0 else x)
            else:
                trawl_can_df = None

            # combine US & CAN trawl files
            if isinstance(trawl_us_df, pd.DataFrame) and isinstance(trawl_can_df, pd.DataFrame):
                self.EPro.trawl_df = pd.concat([trawl_us_df, trawl_can_df])
            else:
                raise SystemError(f"Cannot construct trawl_df for source = 3, since either US or CAN is not available.")

        else:
            raise NotImplementedError(f"Source of {self.EPro.params['source']} not implemented yet.")

    def __load_gear_data(self):

        if self.EPro.params['source'] == 3:

            if self.EPro.params['filename_gear_US']:
                gear_us_df = pd.read_excel(self.EPro.params['data_root_dir'] + self.EPro.params['filename_gear_US'])
                gear_us_df = self.__process_gear_data(gear_us_df)
            else:
                gear_us_df = None

            if self.EPro.params['filename_gear_CAN']:
                gear_can_df = pd.read_excel(self.EPro.params['data_root_dir'] + self.EPro.params['filename_gear_CAN'])
                gear_can_df = self.__process_gear_data(gear_can_df)

                # transect offset is set to zero since we will not have overlap transects
                # TODO: Is this a necessary variable? If so, we should make it an input to the function.
                CAN_Transect_offset = 0

                # combine US & CAN gear files
                gear_can_df.index = gear_can_df.index + self.EPro.params['CAN_haul_offset']  # add haul_offset
                gear_can_df['Transect'] = gear_can_df['Transect'] + CAN_Transect_offset
            else:
                gear_can_df = None

            # combine US & CAN trawl files
            if isinstance(gear_us_df, pd.DataFrame) and isinstance(gear_can_df, pd.DataFrame):
                self.EPro.gear_df = pd.concat([gear_us_df, gear_can_df])
            else:
                raise SystemError(f"Cannot construct gear_df for source = 3, since either US or CAN is not available.")

        else:
            raise NotImplementedError(f"Source of {self.EPro.params['source']} not implemented yet.")

    def __load_specimen_data(self):

        if self.EPro.params['source'] == 3:

            specimen_us_df = pd.read_excel(self.EPro.params['data_root_dir'] + self.EPro.params['filename_specimen_US'])
            specimen_us_df = self.__process_specimen_data(specimen_us_df, 0)

            specimen_can_df = pd.read_excel(self.EPro.params['data_root_dir'] + self.EPro.params['filename_specimen_CAN'])
            specimen_can_df = self.__process_specimen_data(specimen_can_df, self.EPro.params['CAN_haul_offset'])

            self.EPro.specimen_df = pd.concat([specimen_us_df, specimen_can_df])

        else:
            raise NotImplementedError(f"Source of {self.EPro.params['source']} not implemented yet.")

    def __load_biological_data(self):
        """
        A function that loads the appropriate data based on the source provided.

        Returns
        -------
        self.length_df : Dataframe
            A dataframe containing the length information to be used downstream. If source=3 it contains US and Canada.
        self.specimen_df : Dataframe
            A dataframe containing the specimen information to be used downstream. If source=3 it contains US and Canada.
        self.trawl_df : Dataframe
            A dataframe containing the trawl information to be used downstream. If source=3 it contains US and Canada.
        self.gear_df : Dataframe
            A dataframe containing the gear information to be used downstream. If source=3 it contains US and Canada.
        self.catch_df : Dataframe
            A dataframe containing the catch information to be used downstream. If source=3 it contains US and Canada.
        """

        self.__load_catch_data()

        self.__load_length_data()

        self.__load_trawl_data()

        self.__load_gear_data()

        self.__load_specimen_data()

    @staticmethod
    def __get_bin_counts(input_data: np.array, centered_bins: np.array):
        """
        This function manually computes bin counts given ``input_data``. This
        function is computing the histogram of ``input_data`` using
        bins that are centered, rather than bins that are on the edge.
        The first value is between negative infinity and the first bin
        center plus the bin width divided by two. The last value is
        between the second to last bin center plus the bin width
        divided by two to infinity.


        Parameters
        ----------
        input_data: numpy array
            The data to create a histogram of.
        centered_bins: numpy array
            A function that specifies the bin centers.

        Returns
        -------
        hist: numpy array
            The values of the histogram
        """

        bin_diff = np.diff(centered_bins) / 2.0

        hist = [(input_data <= centered_bins[0] + bin_diff[0]).sum()]

        for i in range(len(centered_bins) - 2):
            hist.append(((centered_bins[i] + bin_diff[i] < input_data)
                         & (input_data <= centered_bins[i + 1] + bin_diff[i + 1])).sum())

        hist.append((input_data > centered_bins[-2] + bin_diff[-1]).sum())

        return np.array(hist)

    def process_length_weight_data(self, specimen_df: pd.DataFrame = None):
        """
        process length weight data (all hake trawls) to obtain
        (1) length-weight regression or length-weight-key
        (2) length-age keys

        Parameters
        ----------
        specimen_df: Pandas Dataframe
            Pandas Dataframe describing the length weight data

        """

        # select the indices that do not have nan in either Length or Weight
        len_wgt_nonull = np.logical_and(specimen_df.Length.notnull(), specimen_df.Weight.notnull())

        df_no_null = specimen_df.loc[len_wgt_nonull]

        # all valid length - weight measurements (no nan's)
        L = df_no_null.Length.values
        W = df_no_null.Weight.values
        Lm = df_no_null.Length[df_no_null.Sex == 1].values
        Wm = df_no_null.Weight[df_no_null.Sex == 1].values
        Lf = df_no_null.Length[df_no_null.Sex == 2].values
        Wf = df_no_null.Weight[df_no_null.Sex == 2].values

        # length-weight regression for all trawls (male & female)
        x = np.log10(L)
        y = np.log10(W)

        p = np.polyfit(x, y, 1)  # linear regression

        self.EPro.params['reg_w0'] = 10.0 ** p[1]
        self.EPro.params['reg_p'] = p[0]

        # length-weight regression for all trawls for male
        xM = np.log10(Lm)
        yM = np.log10(Wm)

        pM = np.polyfit(xM, yM, 1)  # linear regression

        self.EPro.params['reg_w0M'] = 10.0 ** pM[1]
        self.EPro.params['reg_pM'] = pM[0]

        # length-weight regression for all trawls for female
        xF = np.log10(Lf)
        yF = np.log10(Wf)

        pF = np.polyfit(xF, yF, 1)  # linear regression

        self.EPro.params['reg_w0F'] = 10.0 ** pF[1]
        self.EPro.params['reg_pF'] = pF[0]

        # total number of fish individuals at length specified by bio_hake_len_bin
        self.EPro.params['len_nM'] = self.__get_bin_counts(Lm, self.EPro.params['bio_hake_len_bin'])
        self.EPro.params['len_nF'] = self.__get_bin_counts(Lf, self.EPro.params['bio_hake_len_bin'])
        self.EPro.params['len_nALL'] = self.__get_bin_counts(L, self.EPro.params['bio_hake_len_bin'])

        # length-key
        self.EPro.params['len_key_M'] = self.EPro.params['len_nM'] / sum(self.EPro.params['len_nM'])
        self.EPro.params['len_key_F'] = self.EPro.params['len_nF'] / sum(self.EPro.params['len_nF'])
        self.EPro.params['len_key_ALL'] = self.EPro.params['len_nALL'] / sum(self.EPro.params['len_nALL'])

        # weight at length or length-weight-key
        # length-weight-key per fish over entire survey region (an array)
        self.EPro.params['len_wgt_M'] = self.EPro.params['reg_w0M'] * self.EPro.params['bio_hake_len_bin'] ** \
                                        self.EPro.params['reg_pM']
        self.EPro.params['len_wgt_F'] = self.EPro.params['reg_w0F'] * self.EPro.params['bio_hake_len_bin'] ** \
                                        self.EPro.params['reg_pF']
        self.EPro.params['len_wgt_ALL'] = self.EPro.params['reg_w0'] * self.EPro.params['bio_hake_len_bin'] ** \
                                          self.EPro.params['reg_p']

        # create length-weight sex structured relations
        for i in range(len(self.EPro.params['len_nM'])):

            # bins with less than 5 samples will be replaced by the regression curve
            # this is done for all Male, Female, and both Females and Males.
            if self.EPro.params['len_nM'][i] >= 5:
                indm = (self.EPro.params['bio_hake_len_bin'][i] - 1 < Lm) & (
                            Lm <= self.EPro.params['bio_hake_len_bin'][i] + 1)
                self.EPro.params['len_wgt_M'][i] = np.mean(Wm[indm])

            if self.EPro.params['len_nF'][i] >= 5:
                indf = (self.EPro.params['bio_hake_len_bin'][i] - 1 < Lf) & (
                            Lf <= self.EPro.params['bio_hake_len_bin'][i] + 1)
                self.EPro.params['len_wgt_F'][i] = np.mean(Wf[indf])

            if self.EPro.params['len_nALL'][i] >= 5:
                ind = (self.EPro.params['bio_hake_len_bin'][i] - 1 < L) & (
                            L <= self.EPro.params['bio_hake_len_bin'][i] + 1)
                self.EPro.params['len_wgt_ALL'][i] = np.mean(W[ind])

        # average length-weight-key per fish over entire survey region (a scalar)
        self.EPro.params['ave_len_wgt_M'] = np.dot(self.EPro.params['len_wgt_M'], self.EPro.params['len_key_M'])
        self.EPro.params['ave_len_wgt_F'] = np.dot(self.EPro.params['len_wgt_F'], self.EPro.params['len_key_F'])
        self.EPro.params['ave_len_wgt_ALL'] = np.dot(self.EPro.params['len_wgt_ALL'], self.EPro.params['len_key_ALL'])

    def construct_catch_trawl_output_matrices(self, KS_stratification, stratification_index):
        """
        construct the hake (or taget species) catch output matrix for
        visualization & generate report tables

        Parameters
        ----------
        KS_stratification : int
            Specifies the type of stratification to be used.
            0 = Pre-Stratification or customized stratification (geographically defined)
            1 = Post-Stratification (KS-based or trawl-based)
        stratification_index : int  # TODO: this looks to mirror KS_stratification, it seems to be unnecessary to use this
                                    # TODO: Ask Chu if it is ok to remove this.
            Index for the chosen stratification
            0 = INPFC strata
            1 = KS (trawl)-based
            2-6 = geographically based but close to trawl-based stratification
            7 = mix-proportion, rather than 85% & 20% hake/hake-mix rules
            10 = one stratum for the whole survey
        """

        print("constructing catch trawl")
        # Create catch table with
        # Columns   1-3:   'Trawl number'   'Abundance'     'Biomass'
        self.final_table_catch = self.EPro.catch_df[self.EPro.catch_df['Species_Code'] == self.EPro.params['species_code_ID']]

        # # TODO: make strata_df a variable of epro_2019
        # # load stratification file
        # stratification_index = 1
        # if stratification_index != 1 and stratification_index != 0:
        #     raise NotImplementedError(f"stratification_index of {stratification_index} has not been implemented!")
        # else:
        #
        #     if stratification_index == 1:
        #         strata_df = pd.read_excel(epro_2019.params['data_root_dir'] + epro_2019.params['filename_strata'],
        #                                   sheet_name='Base KS')
        #         strata_df = strata_df[['Year', 'Cluster name', 'Haul', 'wt']].copy()
        #
        #         # set data types of dataframe
        #         strata_df = strata_df.astype({'Year': int, 'Cluster name': int, 'Haul': int, 'wt': np.float64})
        #
        #         strata_df.set_index('Haul', inplace=True)
        #         strata_df.sort_index(inplace=True)
        #
        #     else:
        #         strata_df = pd.read_excel(epro_2019.params['data_root_dir'] + epro_2019.params['filename_strata'],
        #                                   sheet_name='INPFC')
        #         strata_df = strata_df[['Year', 'INPFC', 'Haul', 'wt']].copy()
        #
        #         # set data types of dataframe
        #         strata_df = strata_df.astype({'Year': int, 'INPFC': int, 'Haul': int, 'wt': np.float64})
        #
        #         strata_df.set_index('Haul', inplace=True)
        #         strata_df.sort_index(inplace=True)
        #
        # strata_df.head()

        # # TODO: make final_table_trawl a variable of epro_2019
        # if stratification_index == 1:
        #     selected_columns = ['Haul', 'Transect', 'EQ_Latitude', 'EQ_Longitude', 'Cluster name',
        #                         'Average_Bottom_Depth',
        #                         'Surface_Temperature', 'Gear_Temperature', 'Length', 'Sex', 'Age', 'Weight',
        #                         'Frequency',
        #                         'Average_Footrope_Depth', 'Weight_In_Haul']
        #     final_table_trawl = pd.DataFrame(columns=selected_columns)
        # elif stratification_index == 0:
        #     selected_columns = ['Haul', 'Transect', 'EQ_Latitude', 'EQ_Longitude', 'INPFC', 'Average_Bottom_Depth',
        #                         'Surface_Temperature', 'Gear_Temperature', 'Length', 'Sex', 'Age', 'Weight',
        #                         'Frequency',
        #                         'Average_Footrope_Depth', 'Weight_In_Haul']
        #     final_table_trawl = pd.DataFrame(columns=selected_columns)
        # else:
        #     raise NotImplementedError(f"stratification_index of {stratification_index} has not been implemented!")
        #
        # catch_specific = epro_2019.catch_df[epro_2019.catch_df['Species_Code'] == epro_2019.params['species_code_ID']]
        #
        # # get all unique Haul values accross all of the data
        # full_haul = np.unique(np.concatenate([strata_df.index.values, epro_2019.gear_df.index.values,
        #                                       catch_specific.index.values, epro_2019.length_ds_wu_jung.Haul.values,
        #                                       epro_2019.trawl_df.index.values,
        #                                       epro_2019.specimen_df.index.unique().values]))
        #
        # full_haul


        # from functools import reduce
        #
        # temp = epro_2019.length_ds_wu_jung.Frequency.to_dataframe().reset_index().dropna(subset='Frequency').set_index(
        #     'Haul')
        #
        # # TODO: change strata_df['Cluster name'] to be more inclusive
        # data_frames_length = [temp, catch_specific, epro_2019.gear_df, epro_2019.trawl_df, strata_df['Cluster name']]
        # df_merged_length = reduce(lambda left, right: pd.merge(left, right, on=['Haul'], how='outer'),
        #                           data_frames_length)
        #
        # df_merged_length['Age'] = np.nan
        # df_merged_length['Weight'] = np.nan
        #
        # # TODO: change strata_df['Cluster name'] to be more inclusive
        # data_frames_specimen = [epro_2019.specimen_df, catch_specific, epro_2019.gear_df, epro_2019.trawl_df,
        #                         strata_df['Cluster name']]
        # df_merged_specimen = reduce(lambda left, right: pd.merge(left, right, on=['Haul'], how='outer'),
        #                             data_frames_specimen)
        #
        # df_merged_specimen['Frequency'] = np.int64(
        #     1)  # TODO: check to make sure if adding a frequency column makes sense with Chu
        #
        # df_merged_length[selected_columns[1:]]
        # # df_merged_specimen[selected_columns[1:]]
        #
        # # TODO: combine df_merged_length and df_merged_specimen


