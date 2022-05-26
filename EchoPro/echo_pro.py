import sys

import yaml
import numpy as np
# from .run_bootstrapping import RunBootstrapping
from .load_biological_data import LoadBioData
from .load_stratification_data import LoadStrataData
from .cv_analysis import CVAnalysis
from .kriging import Kriging
import pandas as pd
import scipy.io
import xarray as xr
import warnings


class EchoPro:
    """
    EchoPro base class that importing and prepares parameters for
    later processes. Additionally, it includes functions for
    acessing the processing, visualization, and generating reports
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
    stratification_index : int  # TODO: this looks to mirror KS_stratification, it seems to be unnecessary to use this
                                # TODO: Ask Chu if it is ok to remove this.
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

    def __load_nasc_data(self):
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

        if self.params['survey_year'] < 2003:
            # TODO: it may be the case that we need to include lines 35-61 of
            #  EchoPro/general/load_files_parameters/get_NASC_data.m
            raise NotImplementedError("Loading the NASC table for survey years less than 2003 has not been implemented!")

        else:
            df.set_index('Transect', inplace=True)
            df.sort_index(inplace=True)

        return df

    def get_final_biomass_table(self):

        nasc_df = self.__load_nasc_data()

        # minimal columns to do Jolly Hampton CV on data that has not been kriged
        self.final_biomass_table = nasc_df[['Latitude', 'Longitude', 'Stratum', 'Spacing']].copy()

        mat = scipy.io.loadmat('../2019_consolidated_files/final_biomass_table_wgt_total.mat')

        self.final_biomass_table["wgt_total"] = mat['wgt_total']

        warnings.warn("We are currently using wgt_total from Matlab for CV, change this!")

    # def get_strata_ds(self):
    #
    #     self.strata_df["length_average_haul"] = np.nan
    #     self.strata_df["TS_lin_haul"] = np.nan
    #     self.strata_df["sig_bs_haul"] = np.nan
    #     for haul_num in self.specimen_df.index.unique():
    #
    #         all_len = self.specimen_df.loc[haul_num]['Length']
    #
    #         # TODO: replace with length_ds?
    #         if haul_num in self.length_df.index:
    #             all_len = np.concatenate([all_len,
    #                                       self.length_df.loc[haul_num]['Length']])
    #
    #         TS0j = 20.0 * np.log10(all_len) - 68.0
    #         TSj = 10.0 * np.log10(np.nanmean(10.0 ** (TS0j / 10.0)))
    #         self.strata_df.loc[haul_num, "TS_lin_haul"] = TSj
    #         self.strata_df.loc[haul_num, "sig_bs_haul"] = 10.0 ** (TSj / 10.0)
    #         self.strata_df.loc[haul_num, "length_average_haul"] = np.nanmean(all_len)
    #
    #     self.strata_ds = self.strata_df.to_xarray()
    #
    #     self.strata_ds["sig_bs"] = self.strata_ds.sig_bs_haul.mean(dim="Haul", skipna=True)
    #     self.strata_ds["sig_b"] = 4.0 * np.pi * self.strata_ds["sig_bs"]

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

        # print("getting strata data")
        # self.load_strata.get_strata_data(stratification_index, KS_stratification, transect_reduction_fraction=0.0)

        self.get_final_biomass_table()
        # self.get_st rata_ds()

    def run_cv_analysis(self, lat_INPFC=None, kriged_data=False):

        if self.params["JH_fac"] == 1:
            nr = 1  # number of realizations
        else:
            nr = 10000  # number of realizations

        cva = CVAnalysis(self)

        if kriged_data:
            raise NotImplementedError("CV analysis for kriged data has not been implemented")
        else:
            cva.run_jolly_hampton(nr, lat_INPFC, self.final_biomass_table)
            return cva

    def run_kriging(self):

        krig = Kriging(self)

        return krig


    # def init_params(self):
    #
    #     # TODO: eventually bulk the below functions
    #
    #     # setting the stratification index based on user provided input
    #     # TODO: Might be able to take out this if else statement depending on downstream items
    #     if KS_stratification == 1:
    #         stratification_index = 1
    #     else:
    #         stratification_index = 0
    #
    #     # check that the stratification index is correctly set
    #     if bio_data_type != 1:
    #         if stratification_index != 0:
    #             stratification_index = 0  # non - acoustical and trawl survey data only use INPFC stratification
    #             print("Changing stratification_index to 0.")
    #     else:
    #         if stratification_index != 1:
    #             print("Changing stratification_index to 1")
    #             stratification_index = 1  # index for the chosen stratification
    #             # 1 = KS(trawl) - based, 2 - 7 = geographically based but close to trawl - based stratification
    #             # 0 = INPFC strata
    #             # 7 = mix - proportion, rather than 85 % & 20 % hake / hake - mix rules
    #             # 10 = one stratum for the whole survey
    #             # TODO: ask about the above comments
    #
    #     # TODO: make sure to take this into account!
    #     # if para.proc.exclude_age1 == 1
    #     #     para.acoust.filename.processed_data = para.acoust.filename.processed_data_age2;
    #     #     para.bio_acoust.filename.Transect_region_haul = para.bio_acoust.filename.Transect_region_haul_age2;
    #     # else
    #     #     para.acoust.filename.processed_data = para.acoust.filename.processed_data_age1;
    #     #     para.bio_acoust.filename.Transect_region_haul = para.bio_acoust.filename.Transect_region_haul_age1;
    #     # end
    #
    #     # check to make sure no survey year and initialization parameters are the same
    #     param_intersect = set(self.__init_params.keys()).intersection(set(self.__survey_params.keys()))
    #
    #     # if no parameters are the same, then run process, else return error
    #     if not param_intersect:
    #         # combine survey year and initialization parameters into one dictionary
    #         full_params = {}
    #         full_params.update(self.__init_params)
    #         full_params.update(self.__survey_params)
    #
    #         return ProcessData(full_params, extrapolation, age_data_status, source, bio_data_type, KS_stratification,
    #                            stratification_index, kriging, kriging_input, exclude_age1, start_transect,
    #                            end_transect, transect_reduction_fraction, transect_reduction_mode, bootstrap_limit,
    #                            default_parameters)
    #     else:
    #         raise RuntimeError('The initialization and survey year configuration files define the same variable! ' +
    #                            f'\n These variables are: {param_intersect}')

    # def run_process(self,
    #                       extrapolation: bool = False,
    #                       KS_stratification: int = 1,
    #                       kriging: bool = True,
    #                       kriging_input: int = 1,
    #                       exclude_age1: bool = False,
    #                       start_transect: int = 1,
    #                       end_transect: int = 200,
    #                       transect_reduction_fraction: float = 0.0,
    #                       transect_reduction_mode: int = 1,
    #                       bootstrap_limit: int = 1,
    #                       default_parameters: bool = True):
    #     """
    #     A function that performs bootstrapping. This involves processing
    #     acoustic and biological data, computation of CV analysis,
    #     and Kriging.
    #
    #     Parameters
    #     ----------
    #     extrapolation : bool
    #             Specifies if extrapolation should be used for Kriging
    #     KS_stratification : int
    #         Specifies the type of stratification to be used.
    #         0 = Pre-Stratification or customized stratification (geographically defined)
    #         1 = Post-Stratification (KS-based or trawl-based)
    #     kriging : bool
    #         States whether or not to perform kriging
    #     kriging_input : int
    #         Specifies TODO: ask Chu what this variable actually specifies
    #         1 = Biomass density
    #         2 = NASC
    #         3 = Number density
    #     exclude_age1 : bool
    #         States whether or not age 1 hake should be included in analysis.
    #     start_transect : int
    #         Value to start transect
    #     end_transect : int
    #         Value to end transect
    #     transect_reduction_fraction : float
    #         Reduction fraction for transect TODO: should this be 5 or 0.05?
    #     transect_reduction_mode : int
    #         1 = Regular
    #         2 = Random
    #     bootstrap_limit : int
    #         The number of bootstraping iterations to perform
    #     default_parameters : bool
    #         States whether or not to use the default parameters
    #     """
    #
    #     # Get all inputs to the function run_bootstrapping()
    #     function_args = locals()
    #
    #     # remove the self argument
    #     del function_args['self']
    #
    #     print(f"saved args = {function_args}")
    #
    #     # import copy
    #     #
    #     # # get copy of EchoPro
    #     # # TODO: This is creating a copy of the EchoPro object that
    #     # # TODO: could eat up memory if the input files are large
    #     # # epro_copy = copy.deepcopy(self)
    #     # #
    #     # # print(self)
    #     #
    #     # bootstrapping_routine = RunBootstrapping(**function_args)
    #     #
    #     # bootstrapping_routine.set_echopro_object(self)
    #
    #     return



