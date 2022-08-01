import numpy as np
import pandas as pd


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
        self.strata_cols = {'Year', 'stratum', 'Haul', 'wt'}

        # expected columns for geo strata Dataframe
        self.geo_strata_cols = {'Strata index', 'Latitude (upper limit)'}

        self._load_stratification_file()
        self._load_geographic_stratification()
        self._get_strata_sig_b()

    def _check_strata_df(self, strata_df: pd.DataFrame) -> None:
        """
        Ensures that the appropriate columns are
        contained in the stratification Dataframe.

        TODO: should we add more in-depth checks here?
        """

        if not set(strata_df.columns).intersection(self.strata_cols):
            raise NameError("Strata dataframe does not contain all expected columns!")

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

            # read and check stratification file
            strata_df = pd.read_excel(self.survey.params['data_root_dir'] + self.survey.params['strata_filename'],
                                      sheet_name=self.survey.params['strata_sheetname'])
            self._check_strata_df(strata_df)

            # extract only those columns that are necessary
            strata_df = strata_df[['Year', 'stratum', 'Haul', 'wt']].copy()

            # set data types of dataframe
            strata_df = strata_df.astype({'Year': int, 'stratum': int,
                                          'Haul': int, 'wt': np.float64})

            # set index of dataframe
            strata_df.set_index(['Haul', 'stratum'], inplace=True)
            strata_df.sort_index(inplace=True)

            self.survey.strata_df = strata_df

        else:
            raise NotImplementedError(f"strata_filename has unknown sheet name!")

    def _check_geo_strata_df(self, geo_strata_df: pd.DataFrame) -> None:
        """
        Ensures that the appropriate columns are
        contained in the geographic stratification Dataframe.

        TODO: should we add more in-depth checks here?
        """

        if not set(geo_strata_df.columns).intersection(self.geo_strata_cols):
            raise NameError("geo_strata dataframe does not contain all expected columns!")

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

            # read and check geographic stratification file
            geo_strata_df = pd.read_excel(self.survey.params['data_root_dir'] + self.survey.params['geo_strata_filename'],
                                          sheet_name=self.survey.params['geo_strata_sheetname'])
            self._check_geo_strata_df(geo_strata_df)

            # extract only those columns that are necessary
            geo_strata_df = geo_strata_df[['Strata index', 'Latitude (upper limit)']].copy()

            # set data types of dataframe
            geo_strata_df = geo_strata_df.astype({'Strata index': int,
                                                  'Latitude (upper limit)': np.float64})

            self.survey.geo_strata_df = geo_strata_df

        else:
            raise NotImplementedError(f"geo_strata_filename has unknown sheet name!")

    def _get_strata_sig_b(self) -> None:
        """
        Computes the backscattering cross-section (sigma_b),
        using the strata, specimen, and length dataframes.
        These values are then stored in self.survey.strata_sig_b
        as a Pandas series with index "stratum".
        """

        # TODO: the target strength functions are specific to Hake, replace with input in the future

        # initialize sig_bs_haul column in strata_df
        self.survey.strata_df["sig_bs_haul"] = np.nan

        # select the indices that do not have nan in either Length or Weight
        spec_df = self.survey.specimen_df[['Length', 'Weight']].copy()
        spec_df = spec_df.dropna(how='any')

        for haul_num in spec_df.index.unique():

            # lengths from specimen file associated with index haul_num
            spec_len = spec_df.loc[haul_num]['Length']

            if haul_num in self.survey.length_df.index:

                # add lengths from length file associated with index haul_num
                length_len = self.survey.length_df.loc[haul_num]['Length'].values
                length_freq = self.survey.length_df.loc[haul_num]['Frequency'].values

                # empirical relation for target strength
                TS0j_length = 20.0 * np.log10(length_len) - 68.0

                # sum of target strengths
                sum_TS0j_length = np.nansum((10.0 ** (TS0j_length / 10.0)) * length_freq)

                # total number of values used to calculate sum_TS0j_length
                num_length = np.nansum(length_freq)

            else:

                # sum of target strengths
                sum_TS0j_length = 0.0

                # total number of values used to calculate sum_TS0j_length
                num_length = 0.0

            # empirical relation for target strength
            TS0j_spec = 20.0 * np.log10(spec_len) - 68.0

            # sum of target strengths
            sum_TS0j_spec = np.nansum(10.0 ** (TS0j_spec / 10.0))

            # mean differential backscattering cross-section for each haul
            self.survey.strata_df.loc[haul_num,
                                    "sig_bs_haul"] = (sum_TS0j_spec + sum_TS0j_length)/(num_length + TS0j_spec.size)

        # mean backscattering cross-section for each stratum
        self.survey.strata_sig_b = 4.0 * np.pi * self.survey.strata_df['sig_bs_haul'].groupby('stratum').mean()
