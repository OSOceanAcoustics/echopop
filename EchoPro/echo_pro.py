import sys

import yaml
import numpy as np
# from .run_bootstrapping import RunBootstrapping
from .load_biological_data import LoadBioData
from .load_stratification_data import LoadStrataData
from .cv_analysis import CVAnalysis
from .kriging import Kriging
from .kriging_mesh import KrigingMesh
from .semivariogram import SemiVariogram
import pandas as pd
import scipy.io
import xarray as xr
import warnings


class EchoPro:
    """
    EchoPro base class that importing and prepares parameters for
    later processes. Additionally, it includes functions for
    accessing the processing, visualization, and generating reports
    Classes.

    Parameters
    ----------
    init_file_path : str
        A string specifying the path to the initialization YAML file
    survey_year_file_path : str
        A string specifying the path to the survey year YAML file
    source : int
        Define the region of data to use.
        1 = US
        2 = Canada
        3 = US and Canada
    bio_data_type : int
        Specifies the biological data to be used.
        1 = Acoustic Survey Trawl Survey
        2 = Bottom Trawl Survey
        3 = Observer Data
    age_data_status : int
        1 = actual age data
        2 = from age_vs_len
    exclude_age1 : bool
        States whether age 1 hake should be included in analysis.
    KS_stratification : int
        Specifies the type of stratification to be used.
        0 = Pre-Stratification or customized stratification (geographically defined)
        1 = Post-Stratification (KS-based or trawl-based)
    stratification_index : int
        Index for the chosen stratification
        0 = INPFC strata
        1 = KS (trawl)-based
        2-6 = geographically based but close to trawl-based stratification
        7 = mix-proportion, rather than 85% & 20% hake/hake-mix rules
        10 = one stratum for the whole survey

    """
    def __init__(self,
                 init_file_path: str,
                 survey_year_file_path: str,
                 source: int = 3,
                 bio_data_type: int = 1,
                 age_data_status: int = 1,
                 exclude_age1: bool = True,
                 KS_stratification: int = 1,
                 stratification_index: int = 1):

        self.bootstrapping_performed = False

        # self.__init_file_path = init_file_path
        self.__check_init_file()

        # self.__survey_year_file_path = survey_year_file_path
        self.__check_survey_year_file()

        init_params = self.__read_initialization_config(init_file_path)
        init_params = self.__set_params_from_init(source, bio_data_type, init_params, age_data_status)

        survey_year_params = self.__read_survey_year_config(survey_year_file_path)

        self.params = self.__collect_parameters(init_params, survey_year_params)

        self.params['exclude_age1'] = exclude_age1

        self.catch_df = None

        self.length_df = None
        self.length_ds_bran = None
        self.length_ds = None

        self.trawl_df = None

        self.gear_df = None

        self.specimen_df = None

        self.strata_df = None
        self.geo_strata_df = None
        self.strata_ds = None

        self.__load_files(KS_stratification, stratification_index)

    def __check_init_file(self):
        # TODO: create this function that checks the contents of the initialization config file
        # TODO: it should make sure that certain variables are defined too
        print("A check of the initialization file needs to be done!")

    def __check_survey_year_file(self):
        # TODO: create this function that checks the contents of the survey year config file
        # TODO: it should make sure that certain variables are defined and all paths exist
        print("A check of the survey year file needs to be done!")

    @staticmethod
    def __read_initialization_config(init_file_path):

        with open(init_file_path) as f:

            init_params = yaml.load(f, Loader=yaml.SafeLoader)

        return init_params

    @staticmethod
    def __read_survey_year_config(survey_year_file_path):
        # TODO: This may need to be replaced in the future with a python file that way we can
        # TODO: mirror proc_parameters_XXXX.m more closely (specifically the switch case statement)

        with open(survey_year_file_path) as f:

            survey_params = yaml.load(f, Loader=yaml.SafeLoader)

        return survey_params

    @staticmethod
    def __set_params_from_init(source: int, bio_data_type: int, init_params: dict, age_data_status: int):

        # setting bio_hake_lin_bin variable to a numpy array
        init_params["bio_hake_len_bin"] = np.linspace(init_params["bio_hake_len_bin"][0],
                                                      init_params["bio_hake_len_bin"][1],
                                                      num=init_params["bio_hake_len_bin"][2],
                                                      dtype=np.int64)

        # setting bio_hake_age_bin variable to a numpy array
        init_params["bio_hake_age_bin"] = np.linspace(init_params["bio_hake_age_bin"][0],
                                                      init_params["bio_hake_age_bin"][1],
                                                      num=init_params["bio_hake_age_bin"][2],
                                                      dtype=np.int64)

        # making acoust_freq0 into a numpy array
        init_params["acoust_freq0"] = np.array(init_params["acoust_freq0"], dtype=np.float64)

        # turning beamwidth input into radians
        init_params["acoust_bw"] = np.array(init_params["acoust_bw"], dtype=np.float64)*np.pi/180.0

        # finding the acoustic freqency index, default = 1 --> 38 kHz
        # TODO: make sure zero indexing doesn't mess up downsteam items
        init_params.update({'acoust_freq_ind': np.intersect1d(np.floor(init_params["acoust_freq"]/1000.0),
                                                              init_params["acoust_freq0"],
                                                              return_indices=True)[2]})

        # Setting acoust_sig_b_coef to a float from string input
        init_params["acoust_sig_b_coef"] = eval(init_params["acoust_sig_b_coef"])

        # set variables based on the source used
        if source <= 3:
            init_params["platform_name"] = 'FSV'
        else:  # TODO look into this and see if this else statement is necessary
            init_params["platform_name"] = 'SD'
            if bio_data_type != 3:
                bio_data_type = 3
                print("Changing bio_data_type to 3 based on platform name.")

        init_params["opr_indx"] = 3  # TODO: look into this variable, might only be necessary for Matlab GUI

        # setting the species code ID based on bio data type
        if bio_data_type == 1:
            init_params["species_code_ID"] = 22500  # target species_code for acoustic survey
        elif bio_data_type == 2:
            init_params["species_code_ID"] = 22500  # target species_code for bottom trawl survey
        elif bio_data_type == 3:
            init_params["species_code_ID"] = 206  # target species_code for industry data(observer data)

        init_params["source"] = source
        init_params["bio_data_type"] = bio_data_type
        init_params["age_data_status"] = age_data_status

        return init_params

    def __set_params_from_survey_year(self):

        # TODO: This portion may need to be set for downstream items, but we might be able to avoid it
        # if para.proc.exclude_age1 == 1
        #     para.acoust.filename.processed_data = para.acoust.filename.processed_data_age2;
        #     para.bio_acoust.filename.Transect_region_haul = para.bio_acoust.filename.Transect_region_haul_age2;
        # else
        #     para.acoust.filename.processed_data = para.acoust.filename.processed_data_age1;
        #     para.bio_acoust.filename.Transect_region_haul = para.bio_acoust.filename.Transect_region_haul_age1;
        # end

        print("Do stuff!")

    def __collect_parameters(self, init_params, survey_params):

        # check to make sure no survey year and initialization parameters are the same
        param_intersect = set(init_params.keys()).intersection(set(survey_params.keys()))

        # if no parameters are the same, then run process, else return error
        if not param_intersect:
            # combine survey year and initialization parameters into one dictionary
            full_params = {}
            full_params.update(init_params)
            full_params.update(survey_params)

        else:
            raise RuntimeError('The initialization and survey year configuration files define the same variable! ' +
                               f'\n These variables are: {param_intersect}')

        return full_params

    def load_nasc_data(self):
        """
        Load VL interval-based NASC table.

        Parameters
        ----------
        exclude_age1 : bool
            States whether age 1 hake should be included in analysis.

        Returns
        -------
        Pandas Dataframe of NASC table.
        """

        if self.params['exclude_age1']:
            df = pd.read_excel(self.params['data_root_dir'] + self.params['filename_processed_data_no_age1'],
                               sheet_name='Sheet1')
        else:
            df = pd.read_excel(self.params['data_root_dir'] + self.params['filename_processed_data_all_ages'],
                               sheet_name='Sheet1')

        # obtaining those columns that are required
        df = df[['Transect', 'Region ID', 'VL start', 'VL end', 'Latitude', 'Longitude', 'Stratum', 'Spacing',
                 'Layer mean depth', 'Layer height', 'Bottom depth', 'NASC', 'Assigned haul']].copy()

        # set data types of dataframe
        df = df.astype({'Transect': int, 'Region ID': int, 'VL start': np.float64, 'VL end': np.float64,
                        'Latitude': np.float64, 'Longitude': np.float64, 'Stratum': int, 'Spacing': np.float64,
                        'Layer mean depth': np.float64, 'Layer height': np.float64, 'Bottom depth': np.float64,
                        'NASC': np.float64, 'Assigned haul': int})

        df.rename(columns={'Assigned haul': 'Haul'}, inplace=True)

        if self.params['survey_year'] < 2003:
            # TODO: it may be the case that we need to include lines 35-61 of
            #  EchoPro/general/load_files_parameters/get_NASC_data.m
            raise NotImplementedError("Loading the NASC table for survey years less than 2003 has not been implemented!")

        else:
            df.set_index('Transect', inplace=True)
            # df.sort_index(inplace=True)

        return df

    def get_final_biomass_table(self):

        nasc_df = self.load_nasc_data()

        # minimal columns to do Jolly Hampton CV on data that has not been kriged
        self.final_biomass_table = nasc_df[['Latitude', 'Longitude', 'Stratum', 'Spacing']].copy()

        mat = scipy.io.loadmat('../2019_consolidated_files/final_biomass_table_nwgt_total.mat')

        self.final_biomass_table["nwgt_total"] = mat['nwgt_total']

        warnings.warn("We are currently using nwgt_total from Matlab for CV, change this!")

    def get_strata_ds(self):

        self.strata_df["length_average_haul"] = np.nan
        self.strata_df["TS_lin_haul"] = np.nan
        self.strata_df["sig_bs_haul"] = np.nan

        for haul_num in self.specimen_df.index.unique():

            all_len = self.specimen_df.loc[haul_num]['Length']

            # TODO: replace with length_ds?
            if haul_num in self.length_df.index:
                all_len = np.concatenate([all_len,
                                          self.length_df.loc[haul_num]['Length']])

            # TODO: these two functions are specific to Hake, replace with input in the future
            TS0j = 20.0 * np.log10(all_len) - 68.0
            TSj = 10.0 * np.log10(np.nanmean(10.0 ** (TS0j / 10.0)))
            self.strata_df.loc[haul_num, "TS_lin_haul"] = TSj
            self.strata_df.loc[haul_num, "sig_bs_haul"] = 10.0 ** (TSj / 10.0)
            self.strata_df.loc[haul_num, "length_average_haul"] = np.nanmean(all_len)

        self.strata_ds = self.strata_df.to_xarray()

        self.strata_ds["sig_bs"] = self.strata_ds.sig_bs_haul.mean(dim="Haul", skipna=True)
        self.strata_ds["sig_b"] = 4.0 * np.pi * self.strata_ds["sig_bs"]

    def __load_files(self, KS_stratification, stratification_index):
        """
        Load the biological, NASC table, stratification file,
        Constructs final trawl and catch tables

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

        Returns
        -------

        """

        self.load_strata = LoadStrataData(self)  # TODO: remove self from load_strata, only necessary for testing

        self.load_bio = LoadBioData(self)  # TODO: remove self from load_bio, only necessary for testing

        if self.params['bio_data_type'] == 1:

            self.load_bio.process_length_weight_data(self.specimen_df)

            self.load_strata.load_stratafication_file(stratification_index)

            # self.load_bio.get_final_catch_trawl_tables(stratification_index) # TODO: this might not be needed make it optional

        else:
            raise NotImplementedError(f"Processing bio_data_type = {self.params['bio_data_type']} has not been implemented!")

        print("getting strata data")
        self.load_strata.get_strata_data(stratification_index, KS_stratification, transect_reduction_fraction=0.0)

        self.get_final_biomass_table()
        self.get_strata_ds()

    def run_cv_analysis(self, lat_INPFC=None, kriged_data=False, seed=None):

        if self.params["JH_fac"] == 1:
            nr = 1  # number of realizations
        else:
            nr = 10000  # number of realizations

        cva = CVAnalysis(self)

        if kriged_data:
            raise NotImplementedError("CV analysis for kriged data has not been implemented")
        else:
            # cva.run_jolly_hampton(nr, lat_INPFC, self.final_biomass_table)
            return cva.run_jolly_hampton(nr, lat_INPFC, self.final_biomass_table, seed)

    def get_kriging_mesh(self):

        krig_mesh = KrigingMesh(self)

        return krig_mesh

    def get_semi_variogram(self, x, y, field):

        semi_vario = SemiVariogram(x, y, field)

        return semi_vario

    def get_kriging(self):

        krig = Kriging(self)

        return krig

    # # import class to use it's functions
    # from EchoPro.load_stratification_data import LoadStrataData
    #
    # strata_class = LoadStrataData(epro_2019)
    #
    # # get df relating the haul to the stratum
    # strata_haul_df = epro_2019.strata_df.reset_index()[['Haul', 'strata']].set_index('Haul')
    #
    # # get the bins for the lengths
    # bins_len = epro_2019.params['bio_hake_len_bin']
    #
    # # get the bins for the ages
    # bins_age = epro_2019.params['bio_hake_age_bin']
    #
    # # get all specimen data that is necessary for key generation
    # spec_w_strata = epro_2019.specimen_df.drop('Specimen_Number', axis=1).copy().reset_index()
    #
    # # add strata column
    # spec_w_strata['Strata'] = spec_w_strata.apply(lambda x: strata_haul_df.loc[x[0]],
    #                                               axis=1).values
    #
    # spec_w_strata.set_index('Strata', inplace=True)
    #
    # # spec_w_strata = spec_w_strata[(spec_w_strata['Sex'] != 3)].copy() # TODO: this should be for all sexes
    #
    # age_len_key_da, age_len_key_wgt_da, age_len_key_norm_da = strata_class.get_age_key_das(spec_w_strata,
    #                                                                                        bins_len, bins_age)
    #
    # # TODO: it would probably be better to do an average of station 1 and 2 here... (Chu doesn't do this)
    # age_len_key_wgt_norm_da = age_len_key_wgt_da / age_len_key_wgt_da.sum(dim=['len_bins', 'age_bins'])
    #
    # # each stratum's multiplier once normalized weight has been calculated
    # age2_wgt_proportion_da = 1.0 - age_len_key_wgt_norm_da.isel(age_bins=0).sum(
    #     dim='len_bins') / age_len_key_wgt_norm_da.sum(dim=['len_bins', 'age_bins'])
    #
    # # get all specimen data that is necessary for key generation
    # spec_w_strata = epro_2019.specimen_df.drop('Specimen_Number', axis=1).copy().reset_index()
    #
    # # add strata column
    # spec_w_strata['Strata'] = spec_w_strata.apply(lambda x: strata_haul_df.loc[x[0]],
    #                                               axis=1).values
    #
    # spec_w_strata.set_index('Strata', inplace=True)
    #
    # reg_w0, reg_p = strata_class.get_length_val_reg_vals(len_name='Length', val_name="Weight", df=spec_w_strata)
    #
    # len_weight_ALL, len_nALL, norm_len_key_ALL = strata_class.generate_length_val_key(bins_len, reg_w0, reg_p,
    #                                                                                   len_name='Length',
    #                                                                                   val_name='Weight',
    #                                                                                   df=spec_w_strata)
    #
    # # spec_w_strata = spec_w_strata[(spec_w_strata['Sex'] != 3)] # TODO: this should be for all sexes
    # len_wgt_key_spec_da, len_key_spec_da, len_key_norm_spec_da = strata_class.get_weight_key_das(spec_w_strata,
    #                                                                                              bins_len, reg_w0,
    #                                                                                              reg_p,
    #                                                                                              len_name='Length',
    #                                                                                              val_name='Weight')
    #
    # length_explode_df = epro_2019.length_df[['Sex', 'Length']].copy()
    # # add strata column
    # length_explode_df['Strata'] = length_explode_df.reset_index().apply(lambda x: strata_haul_df.loc[x[0]],
    #                                                                     axis=1).values
    #
    # length_explode_df.reset_index(inplace=True)
    #
    # length_explode_df.set_index('Strata', inplace=True)
    #
    # length_explode_df = length_explode_df.explode(['Sex', 'Length'])
    #
    # length_explode_df = length_explode_df.astype({'Haul': int,
    #                                               'Sex': int,
    #                                               'Length': np.float64})
    #
    # # length_explode_df = length_explode_df[(length_explode_df['Sex'] != 3)] # TODO: this should be for all sexes
    #
    # unique_strata = length_explode_df.index.unique().values
    #
    # len_key_norm_length = np.empty((unique_strata.shape[0], bins_len.shape[0]), dtype=np.float64)
    # len_key_norm_length[:, :] = 0.0
    #
    # stratum_ind = 0
    # for stratum in unique_strata:
    #     input_data = length_explode_df.loc[stratum]['Length'].values
    #     len_ind = strata_class.get_bin_ind(input_data, bins_len)
    #
    #     len_key_n = np.array([i.shape[0] for i in len_ind])
    #     len_key_norm_length[stratum_ind, :] = len_key_n / np.sum(len_key_n)
    #
    #     stratum_ind += 1
    #
    # len_key_norm_length_da = xr.DataArray(data=len_key_norm_length,
    #                                       coords={'strata': unique_strata, 'len_bins': bins_len})
    #
    # len_key_norm_ave = (len_key_norm_length_da + len_key_norm_spec_da) / 2
    #
    # # get the nasc dataframe
    # nasc_df = epro_2019.load_nasc_data()
    #
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
    #
    # bio_dense_df = nasc_df[['Stratum', 'NASC', 'Haul']].copy()
    # bio_dense_df['interval'] = interval
    # bio_dense_df['n_A'] = nasc_df.apply(lambda x: x.NASC / epro_2019.strata_ds.sig_b.sel(strata=x.Stratum).values,
    #                                     axis=1)
    # # bio_dense_df['A'] = bio_dense_df['interval']*nasc_df['Spacing']
    # # bio_dense_df['N_A'] = bio_dense_df['n_A']*bio_dense_df['A']
    #
    # bio_density = bio_dense_df.apply(lambda x: x.n_A * np.dot(len_key_norm_ave.sel(strata=x.Stratum),
    #                                                           len_wgt_key_spec_da.sel(strata=x.Stratum)), axis=1)
    #
    # bio_density_2_prop = bio_dense_df.apply(lambda x: x.n_A * np.dot(len_key_norm_ave.sel(strata=x.Stratum),
    #                                                                  len_wgt_key_spec_da.sel(
    #                                                                      strata=x.Stratum)) * age2_wgt_proportion_da.sel(
    #     strata=x.Stratum).values,
    #                                         axis=1)


