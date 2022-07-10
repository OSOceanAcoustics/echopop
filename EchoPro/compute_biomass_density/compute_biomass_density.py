import numpy as np
import pandas as pd
import xarray as xr


class ComputeBiomassDensity:
    """
    A class that computes the biomass density
    of an animal population based off of NASC
    values and other associated data.

    Parameters
    ----------
    EPro : EchoPro object
        An initialized EchoPro object. Note that any change to
        self.EPro will also change this object.
    """

    def __init__(self, EPro = None):

        self.EPro = EPro

        self.bio_hake_len_bin = EPro.params['bio_hake_len_bin']
        self.bio_hake_age_bin = EPro.params['bio_hake_age_bin']

    @staticmethod
    def get_bin_ind(input_data: np.array, centered_bins: np.array):
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
            An array that specifies the bin centers.

        Returns
        -------
        hist_ind: numpy array
            The index values of input_data corresponding to the histogram

        """

        bin_diff = np.diff(centered_bins) / 2.0

        hist_ind = [np.argwhere(input_data <= centered_bins[0] + bin_diff[0]).flatten()]

        for i in range(len(centered_bins) - 2):
            # get values greater than lower bound
            g_lb = centered_bins[i] + bin_diff[i] < input_data

            # get values less than or equal to the upper bound
            le_ub = input_data <= centered_bins[i + 1] + bin_diff[i + 1]

            hist_ind.append(np.argwhere(g_lb & le_ub).flatten())

        hist_ind.append(np.argwhere(input_data > centered_bins[-2] + bin_diff[-1]).flatten())

        return hist_ind

    def get_age_key_das(self, spec_w_strata, bins_len, bins_age):
        """
        Constructs the DataArrays corresponding to the keys
        associated with age.

        Parameters
        ----------
        spec_w_strata : Pandas Dataframe
            Species Dataframe with the strata as an index and
            columns "Length", "Age", "Weight".
        bins_len : Numpy array
            1D numpy array specifies the bin centers for lengths.
        bins_age : Numpy array
            1D numpy array specifies the bin centers for ages.

        Returns
        -------
        The Xarray DataArrays: age_len_key_da, age_len_key_wgt_da,
        and age_len_key_norm_da with coordinates
        (strata, bins_len, bins_age).

        """

        unique_strata = np.sort(spec_w_strata.index.unique())

        age_len_key = np.empty((unique_strata.shape[0],
                                bins_len.shape[0], bins_age.shape[0]), dtype=int)

        age_len_key_wgt = np.empty((unique_strata.shape[0],
                                    bins_len.shape[0], bins_age.shape[0]), dtype=np.float64)

        age_len_key_norm = np.empty((unique_strata.shape[0],
                                     bins_len.shape[0], bins_age.shape[0]), dtype=np.float64)

        stratum_i = 0
        for stratum in unique_strata:

            input_arr_len = spec_w_strata.loc[stratum].Length.values
            input_arr_age = spec_w_strata.loc[stratum].Age.values
            input_arr_wgt = spec_w_strata.loc[stratum].Weight.values

            age_bins_ind = self.get_bin_ind(input_arr_age, bins_age)

            i = 0
            for arr_ind in age_bins_ind:
                temp = self.get_bin_ind(input_arr_len[arr_ind], bins_len)
                age_len_key[stratum_i, :, i] = np.array([i.shape[0] for i in temp])

                age_len_key_wgt[stratum_i, :, i] = np.array([np.sum(input_arr_wgt[arr_ind[i]]) for i in temp])
                i += 1

            # normalized age_len_ky
            age_len_key_norm[stratum_i, :, :] = age_len_key[stratum_i, :, :] / np.sum(age_len_key[stratum_i, :, :])

            stratum_i += 1

        age_len_key_da = xr.DataArray(data=age_len_key,
                                      coords={'strata': unique_strata,
                                              'len_bins': bins_len, 'age_bins': bins_age})

        age_len_key_wgt_da = xr.DataArray(data=age_len_key_wgt,
                                          coords={'strata': unique_strata,
                                                  'len_bins': bins_len, 'age_bins': bins_age})

        age_len_key_norm_da = xr.DataArray(data=age_len_key_norm,
                                           coords={'strata': unique_strata,
                                                   'len_bins': bins_len, 'age_bins': bins_age})

        return age_len_key_da, age_len_key_wgt_da, age_len_key_norm_da

    @staticmethod
    def get_length_val_reg_vals(len_name, val_name, df: pd.DataFrame = None):
        """
        Obtains the regression values for a provided length and value.

        Parameters
        ----------
        len_name : str
            Name of column in df corresponding to the length
        val_name : str
            Name of column in df corresponding to the value
        df : Pandas Dataframe
            Dataframe containing the length and value columns


        Returns
        -------
        reg_w0 and reg_p used in the calculation
        reg_w0 * bio_hake_len_bin ** reg_p
        """

        # select the indices that do not have nan in either Length or Weight
        len_wgt_nonull = np.logical_and(df[len_name].notnull(), df[val_name].notnull())

        df_no_null = df.loc[len_wgt_nonull]

        L = df_no_null[len_name].values
        V = df_no_null[val_name].values

        # length-value regression for all trawls (male & female)
        x = np.log10(L)
        y = np.log10(V)

        p = np.polyfit(x, y, 1)  # linear regression

        reg_w0 = 10.0 ** p[1]
        reg_p = p[0]

        return reg_w0, reg_p

    def generate_length_val_key(self, bio_hake_len_bin, reg_w0, reg_p,
                                len_name, val_name, df: pd.DataFrame = None):
        """
        process length weight data (all hake trawls) to obtain
        (1) length-weight regression or length-weight-key
        (2) length-age keys

        Parameters
        ----------
        df: Pandas Dataframe
            Pandas Dataframe describing the length weight data

        """

        # select the indices that do not have nan in either Length or Weight
        len_wgt_nonull = np.logical_and(df[len_name].notnull(), df[val_name].notnull())

        df_no_null = df.loc[len_wgt_nonull]

        L = df_no_null[len_name].values
        V = df_no_null[val_name].values

        # total number of fish individuals at length specified by bio_hake_len_bin
        len_nALL = self.get_bin_ind(L, bio_hake_len_bin)
        len_nALL = np.array([i.shape[0] for i in len_nALL])

        # normalized length-key
        norm_len_key_ALL = len_nALL / sum(len_nALL)

        if (reg_p is None) or (reg_w0 is None):

            # length-value regression for all trawls (male & female)
            x = np.log10(L)
            y = np.log10(V)

            p = np.polyfit(x, y, 1)  # linear regression

            reg_w0 = 10.0 ** p[1]
            reg_p = p[0]

        # value at length or length-value-key
        # length-weight-key per fish over entire survey region (an array)
        len_val_ALL = reg_w0 * bio_hake_len_bin ** reg_p

        # create length-value structured relations
        for i in range(len(len_nALL)):

            # bins with less than 5 samples will be replaced by the regression curve
            if len_nALL[i] >= 5:
                ind = (bio_hake_len_bin[i] - 1 < L) & (
                            L <= bio_hake_len_bin[i] + 1)
                len_val_ALL[i] = np.mean(V[ind])

        return len_val_ALL, len_nALL, norm_len_key_ALL

    def get_weight_key_das(self, spec_w_strata, bins_len, reg_w0, reg_p, len_name='Length',
                           val_name='Weight'):

        unique_strata = np.sort(spec_w_strata.index.unique())

        len_wgt_key = np.empty((unique_strata.shape[0], bins_len.shape[0]))
        len_wgt_key[:, :] = np.nan

        len_key = np.empty((unique_strata.shape[0], bins_len.shape[0]))
        len_key[:, :] = np.nan

        len_key_norm = np.empty((unique_strata.shape[0], bins_len.shape[0]))
        len_key_norm[:, :] = np.nan

        stratum_i = 0
        for stratum in unique_strata:
            l_w_k, l_k, l_k_n = self.generate_length_val_key(bins_len, reg_w0, reg_p,
                                                             len_name=len_name,
                                                             val_name=val_name,
                                                             df=spec_w_strata.loc[stratum])

            len_wgt_key[stratum_i, :] = l_w_k
            len_key[stratum_i, :] = l_k
            len_key_norm[stratum_i, :] = l_k_n

            stratum_i += 1

        len_wgt_key_da = xr.DataArray(data=len_wgt_key,
                                      coords={'strata': unique_strata, 'len_bins': bins_len})

        len_key_da = xr.DataArray(data=len_key,
                                  coords={'strata': unique_strata, 'len_bins': bins_len})

        len_key_norm_da = xr.DataArray(data=len_key_norm,
                                       coords={'strata': unique_strata, 'len_bins': bins_len})

        return len_wgt_key_da, len_key_da, len_key_norm_da

    def get_len_keys_strata(self, len_strata, bins_len):

        input_data = len_strata['Length'].values
        len_ind = self.get_bin_ind(input_data, bins_len)
        len_key = np.array([i.shape[0] for i in len_ind])
        len_key_n = len_key / np.sum(len_key)

        input_data = len_strata[len_strata['Sex'] == 1]['Length'].values
        len_ind = self.get_bin_ind(input_data, bins_len)
        len_key = np.array([i.shape[0] for i in len_ind])
        len_key_M = len_key / np.sum(len_key)

        input_data = len_strata[len_strata['Sex'] == 2]['Length'].values
        len_ind = self.get_bin_ind(input_data, bins_len)
        len_key = np.array([i.shape[0] for i in len_ind])
        len_key_F = len_key / np.sum(len_key)

        return len_key_n, len_key_M, len_key_F

    def get_age_related_sums(self, spec_strata_M, spec_strata_F, bins_len, bins_age):

        input_df = pd.concat([spec_strata_M, spec_strata_F], axis=0)

        # TODO: clean up outputs of this function call
        _, _, age_len_key_norm_da = self.get_age_key_das(input_df, bins_len, bins_age)
        _, _, age_len_key_norm_M_da = self.get_age_key_das(spec_strata_M, bins_len, bins_age)
        _, _, age_len_key_norm_F_da = self.get_age_key_das(spec_strata_F, bins_len, bins_age)

        len_age_key_sum = age_len_key_norm_da.isel(strata=0).sum(dim='age_bins').values

        len_age_key_M_sum = age_len_key_norm_M_da.isel(strata=0).sum(dim='age_bins').values

        len_age_key_F_sum = age_len_key_norm_F_da.isel(strata=0).sum(dim='age_bins').values

        return len_age_key_sum, len_age_key_M_sum, len_age_key_F_sum

    def get_biomass_constants(self, spec_w_strata, length_explode_df, bins_len, bins_age):
        """
        Obtains the constants associated with each stratum,
        which are used in the biomass density calculation
        """

        # TODO: clean up this code! There may be things saved here
        #  that do not need to be saved...

        spec_strata_ind = spec_w_strata.index.unique()
        len_strata_ind = length_explode_df.index.unique()
        strata_ind = spec_strata_ind.intersection(len_strata_ind).values

        total_N = np.empty(strata_ind.shape[0], dtype=np.float64)

        spec_M_prop = np.empty(strata_ind.shape[0], dtype=np.float64)
        spec_F_prop = np.empty(strata_ind.shape[0], dtype=np.float64)

        len_M_prop = np.empty(strata_ind.shape[0], dtype=np.float64)
        len_F_prop = np.empty(strata_ind.shape[0], dtype=np.float64)

        fac2_ALL = np.empty(strata_ind.shape[0], dtype=np.float64)
        fac1_ALL = np.empty(strata_ind.shape[0], dtype=np.float64)

        fac1_M = np.empty(strata_ind.shape[0], dtype=np.float64)
        fac1_F = np.empty(strata_ind.shape[0], dtype=np.float64)

        fac2_M = np.empty(strata_ind.shape[0], dtype=np.float64)
        fac2_F = np.empty(strata_ind.shape[0], dtype=np.float64)

        len_wgt_prod = np.empty(strata_ind.shape[0], dtype=np.float64)
        len_wgt_M_prod = np.empty(strata_ind.shape[0], dtype=np.float64)
        len_wgt_F_prod= np.empty(strata_ind.shape[0], dtype=np.float64)

        # TODO: clean up outputs of this function
        len_weight_ALL, _, _ = self.generate_length_val_key(bins_len, reg_w0=None, reg_p=None,
                                                            len_name='Length',
                                                            val_name='Weight', df=spec_w_strata)

        # # select the indices that do not have nan in either Length or Weight
        spec_w_strata = spec_w_strata.dropna(how='any')

        # select the indices that do not have nan in either Length or Weight
        length_explode_df = length_explode_df.dropna(how='any')

        ind = 0
        for stratum in strata_ind:

            spec_strata = spec_w_strata.loc[stratum]
            spec_strata_M = spec_strata[spec_strata['Sex'] == 1]
            spec_strata_F = spec_strata[spec_strata['Sex'] == 2]
            len_strata = length_explode_df.loc[stratum]
            len_strata_M = len_strata[len_strata['Sex'] == 1]
            len_strata_F = len_strata[len_strata['Sex'] == 2]

            len_age_key_sum, len_age_key_M_sum, len_age_key_F_sum \
                = self.get_age_related_sums(spec_strata_M, spec_strata_F, bins_len, bins_age)

            # total number of sexed fish at stations 1 and 2
            total_N[ind] = spec_strata_M.shape[0] + spec_strata_F.shape[0] + len_strata.shape[0]
            # total_N[ind] = spec_strata.shape[0] + len_strata.shape[0]  # TODO: This is what it should be

            # proportion of males/females in station 2
            spec_M_prop[ind] = spec_strata_M.shape[0] / total_N[ind]
            spec_F_prop[ind] = spec_strata_F.shape[0] / total_N[ind]

            # proportion of males/females in station 1
            len_M_prop[ind] = len_strata_M.shape[0] / total_N[ind]
            len_F_prop[ind] = len_strata_F.shape[0] / total_N[ind]

            # total proportion of sexed fish in station 2
            fac2_ALL[ind] = spec_M_prop[ind] + spec_F_prop[ind]

            # total proportion of sexed fish in station 1
            fac1_ALL[ind] = 1.0 - fac2_ALL[ind]

            fac1_M[ind] = fac1_ALL[ind] * 1.0
            fac1_F[ind] = fac1_ALL[ind] * 1.0  # TODO: something looks wrong here, ask Chu

            fac1_M[ind] = fac1_M[ind] / (fac1_M[ind] + spec_M_prop[ind])
            fac1_F[ind] = fac1_F[ind] / (fac1_F[ind] + spec_F_prop[ind])

            fac2_M[ind] = spec_M_prop[ind] / (fac1_M[ind] + spec_M_prop[ind])
            fac2_F[ind] = spec_F_prop[ind] / (fac1_F[ind] + spec_F_prop[ind])

            # average of total proportion of sexed fish pertaining to station 1
            # fac1_ALL = fac1_ALL / (fac1_ALL + fac2_ALL) # TODO: do we need to calculate this?

            # average of total proportion of sexed fish pertaining to station 2
            # fac2_ALL = fac2_ALL / (fac1_ALL + fac2_ALL) # TODO: do we need to calculate this?

            len_key_n, len_key_M, len_key_F = self.get_len_keys_strata(len_strata, bins_len)

            len_wgt_prod[ind] = np.dot(fac1_ALL[ind] * len_key_n + fac2_ALL[ind] * len_age_key_sum, len_weight_ALL)
            len_wgt_M_prod[ind] = np.dot(fac1_M[ind] * len_key_M + fac2_M[ind] * len_age_key_M_sum, len_weight_ALL)
            len_wgt_F_prod[ind] = np.dot(fac1_F[ind] * len_key_F + fac2_F[ind] * len_age_key_F_sum, len_weight_ALL)

            ind += 1

        total_N_da = xr.DataArray(data=total_N, coords={'strata': strata_ind})

        spec_M_prop_da = xr.DataArray(data=spec_M_prop, coords={'strata': strata_ind})
        spec_F_prop_da = xr.DataArray(data=spec_F_prop, coords={'strata': strata_ind})

        len_M_prop_da = xr.DataArray(data=len_M_prop, coords={'strata': strata_ind})
        len_F_prop_da = xr.DataArray(data=len_F_prop, coords={'strata': strata_ind})

        len_wgt_prod_da = xr.DataArray(data=len_wgt_prod, coords={'strata': strata_ind})
        len_wgt_M_prod_da = xr.DataArray(data=len_wgt_M_prod, coords={'strata': strata_ind})
        len_wgt_F_prod_da = xr.DataArray(data=len_wgt_F_prod, coords={'strata': strata_ind})

        bio_calc_const = xr.Dataset({'spec_M_prop': spec_M_prop_da, 'spec_F_prop': spec_F_prop_da,
                                     'len_M_prop': len_M_prop_da, 'len_F_prop': len_F_prop_da,
                                     'len_wgt_prod': len_wgt_prod_da, 'len_wgt_M_prod': len_wgt_M_prod_da,
                                     'len_wgt_F_prod': len_wgt_F_prod_da, 'total_N': total_N_da})

        return bio_calc_const

    def get_final_biomass_table(self):

        # TODO: clean up this code!! Is it necessary to create bio_dense_df?

        # minimal columns to do Jolly Hampton CV on data that has not been kriged
        self.EPro.final_biomass_table = self.EPro.nasc_df[['Latitude', 'Longitude', 'Stratum', 'Spacing']].copy()

        # get df relating the haul to the stratum
        strata_haul_df = self.EPro.strata_df.reset_index()[['Haul', 'strata']].set_index('Haul')

        # get all specimen data that is necessary for key generation  # TODO: we may be able to remove this line
        spec_w_strata = self.EPro.specimen_df.copy().reset_index()

        # add strata column
        spec_w_strata['Strata'] = spec_w_strata.apply(lambda x: strata_haul_df.loc[x[0]],
                                                      axis=1).values

        spec_w_strata.set_index('Strata', inplace=True)

        length_explode_df = self.EPro.length_df[['Sex', 'Length']].copy()
        # add strata column
        length_explode_df['Strata'] = length_explode_df.reset_index().apply(lambda x: strata_haul_df.loc[x[0]],
                                                                            axis=1).values

        length_explode_df.reset_index(inplace=True)

        length_explode_df.set_index('Strata', inplace=True)

        length_explode_df = length_explode_df.explode(['Sex', 'Length'])

        length_explode_df = length_explode_df.astype({'Haul': int,
                                                      'Sex': int,
                                                      'Length': np.float64})

        bc = self.get_biomass_constants(spec_w_strata, length_explode_df,
                                        self.bio_hake_len_bin, self.bio_hake_age_bin)

        # TODO: should we include the below code that is commented out?
        # # calculates the interval for the area calculation
        # interval = (nasc_df['VL start'].iloc[1:].values - nasc_df['VL start'].iloc[:-1].values)
        # last_interval = nasc_df['VL end'].iloc[-1] - nasc_df['VL start'].iloc[-1]
        #
        # interval = np.concatenate([interval, np.array([last_interval])])
        #
        # median_interval = np.median(interval)
        #
        # # remove outliers at the end of the transect
        # ind_outliers = np.argwhere(np.abs(interval - median_interval) > 0.05).flatten()
        # interval[ind_outliers] = nasc_df['VL end'].values[ind_outliers] - nasc_df['VL start'].values[ind_outliers]

        bio_dense_df = self.EPro.nasc_df[['Stratum', 'NASC', 'Haul']].copy()
        # bio_dense_df['interval'] = interval

        wgt_vals = self.EPro.strata_df.reset_index().set_index('Haul')['wt']
        wgt_vals_ind = wgt_vals.index

        mix_sa_ratio = self.EPro.nasc_df.apply(lambda x: wgt_vals[x.Haul] if x.Haul in wgt_vals_ind else 0.0, axis=1)

        self.EPro.nasc_df['mix_sa_ratio'] = mix_sa_ratio

        bio_dense_df['n_A'] = self.EPro.nasc_df.apply(
            lambda x: np.round((x.mix_sa_ratio * x.NASC) / float(self.EPro.strata_ds.sig_b.sel(strata=x.Stratum))),
            axis=1)
        # bio_dense_df['A'] = bio_dense_df['interval'] * nasc_df['Spacing']
        # bio_dense_df['N_A'] = bio_dense_df['n_A'] * bio_dense_df['A']

        nntk_male = bio_dense_df.apply(lambda x: np.round(
            x.n_A * float(bc.len_M_prop.sel(strata=x.Stratum) + bc.spec_M_prop.sel(strata=x.Stratum))), axis=1)
        nntk_female = bio_dense_df.apply(lambda x: np.round(
            x.n_A * float(bc.len_F_prop.sel(strata=x.Stratum) + bc.spec_F_prop.sel(strata=x.Stratum))), axis=1)

        bio_dense_df['nntk_male'] = nntk_male
        bio_dense_df['nntk_female'] = nntk_female

        nWgt_male_int = bio_dense_df.apply(lambda x: x.nntk_male * float(bc.len_wgt_M_prod.sel(strata=x.Stratum)),
                                           axis=1)
        nWgt_female_int = bio_dense_df.apply(lambda x: x.nntk_female * float(bc.len_wgt_F_prod.sel(strata=x.Stratum)),
                                             axis=1)

        bio_dense_df['nWgt_male'] = nWgt_male_int
        bio_dense_df['nWgt_female'] = nWgt_female_int

        nWgt_unsexed_int = bio_dense_df.apply(
            lambda x: (x.n_A - x.nntk_male - x.nntk_female) * float(bc.len_wgt_prod.sel(strata=x.Stratum)), axis=1)
        bio_dense_df['nWgt_unsexed'] = nWgt_unsexed_int

        bio_dense_df['nWgt_total'] = bio_dense_df['nWgt_male'] + bio_dense_df['nWgt_female'] + bio_dense_df[
            'nWgt_unsexed']

        spec_w_strata = spec_w_strata.dropna(how='any')
        age_len_key_da, age_len_key_wgt_da, age_len_key_norm_da = self.get_age_key_das(spec_w_strata,
                                                                                       self.bio_hake_len_bin,
                                                                                       self.bio_hake_age_bin)

        # TODO: it would probably be better to do an average of station 1 and 2 here... (Chu doesn't do this)
        age_len_key_wgt_norm_da = age_len_key_wgt_da / age_len_key_wgt_da.sum(dim=['len_bins', 'age_bins'])

        # each stratum's multiplier once normalized weight has been calculated
        age2_wgt_proportion_da = 1.0 - age_len_key_wgt_norm_da.isel(age_bins=0).sum(
            dim='len_bins') / age_len_key_wgt_norm_da.sum(dim=['len_bins', 'age_bins'])

        nWgt_total_2_prop = bio_dense_df.apply(
            lambda x: x.nWgt_total * float(age2_wgt_proportion_da.sel(strata=x.Stratum)),
            axis=1)

        self.EPro.final_biomass_table["normalized_biomass_density"] = nWgt_total_2_prop

        print("We are using our own biomass density calculation!")
