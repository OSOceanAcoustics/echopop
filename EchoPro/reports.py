from .survey import Survey
import pandas as pd
import numpy as np
import xarray as xr


class Reports:
    """
    A Class that writes requested variables to
    consolidated files.

    Parameters
    ----------
    echopro : EchoPro
        EchoPro object that contains all necessary parameters

    """

    def __init__(self,
                 survey: Survey):

        if survey.bootstrapping_performed:
            print("Generate all consolidated files")
        else:
            raise RuntimeError("Bootstrapping must be performed first!")

        self.EPro = survey

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

    def __load_gear_data(self):

        # TODO: if this function ends up being used, we need to document and create a check for df

        if self.EPro.params['source'] == 3:

            if self.EPro.params['filename_gear_US']:
                gear_us_df = pd.read_excel(self.EPro.params['data_root_dir'] + self.EPro.params['filename_gear_US'],
                                           sheet_name='biodata_gear')
                gear_us_df = self.__process_gear_data(gear_us_df)
            else:
                gear_us_df = None

            if self.EPro.params['filename_gear_CAN']:
                gear_can_df = pd.read_excel(self.EPro.params['data_root_dir'] + self.EPro.params['filename_gear_CAN'],
                                            sheet_name='biodata_gear_CAN')
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

    def __load_catch_data(self):

        # TODO: if this function ends up being used, we need to document and create a check for df

        if self.EPro.params['source'] == 3:

            if self.EPro.params['database_type'] == 'Oracle':
                catch_us_df = pd.read_excel(
                    self.EPro.params['data_root_dir'] + self.EPro.params['filename_catch_US'],
                    sheet_name='biodata_catch')
                catch_us_df = self.__process_catch_data_df(catch_us_df)

                catch_can_df = pd.read_excel(
                    self.EPro.params['data_root_dir'] + self.EPro.params['filename_catch_CAN'],
                    sheet_name='biodata_catch_CAN')
                catch_can_df = self.__process_catch_data_df(catch_can_df)

                # combine US & CAN catch files
                catch_can_df.index = catch_can_df.index + self.EPro.params["CAN_haul_offset"]  # add haul_offset
                self.EPro.catch_df = pd.concat([catch_us_df, catch_can_df])

            else:
                raise NotImplementedError("Loading data from a non-Oracle database has not been implemented!")
        else:

            raise NotImplementedError(f"Source of {self.EPro.params['source']} not implemented yet.")

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
        df = df.astype(
            {'Haul': int, 'Haul_Type': int, 'Performance_Code': np.float64, 'Distance_Fished': np.float64,
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

    def __load_trawl_data(self):

        # TODO: if this function ends up being used, we need to document and create a check for df

        if self.EPro.params['source'] == 3:

            if self.EPro.params['filename_trawl_US']:
                trawl_us_df = pd.read_excel(
                    self.EPro.params['data_root_dir'] + self.EPro.params['filename_trawl_US'],
                    converters={'Duration': str}, sheet_name='biodata_haul')
                trawl_us_df = self.__process_trawl_data(trawl_us_df)
            else:
                trawl_us_df = None

            if self.EPro.params['filename_trawl_CAN']:
                trawl_can_df = pd.read_excel(
                    self.EPro.params['data_root_dir'] + self.EPro.params['filename_trawl_CAN'],
                    converters={'Duration': str}, sheet_name='biodata_haul_CAN')
                trawl_can_df = self.__process_trawl_data(trawl_can_df)

                trawl_can_df.index = trawl_can_df.index + self.EPro.params['CAN_haul_offset']  # add haul_offset

                # change lon to -lon
                trawl_can_df['EQ_Longitude'] = trawl_can_df['EQ_Longitude'].apply(
                    lambda x: -1.0 * x if x > 0 else x)
            else:
                trawl_can_df = None

            # combine US & CAN trawl files
            if isinstance(trawl_us_df, pd.DataFrame) and isinstance(trawl_can_df, pd.DataFrame):
                self.EPro.trawl_df = pd.concat([trawl_us_df, trawl_can_df])
            else:
                raise SystemError(
                    f"Cannot construct trawl_df for source = 3, since either US or CAN is not available.")

        else:
            raise NotImplementedError(f"Source of {self.EPro.params['source']} not implemented yet.")