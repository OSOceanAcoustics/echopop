import yaml
import numpy as np
# from .run_bootstrapping import RunBootstrapping
from .load_biological_data import LoadBioData
from .load_stratification_data import LoadStrataData
from .compute_biomass_density import ComputeBiomassDensity
from .cv_analysis import CVAnalysis
from .kriging import Kriging
from .kriging_mesh import KrigingMesh
from .semivariogram import SemiVariogram
import pandas as pd


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
                 stratification_index: int = 1):

        self.bootstrapping_performed = False

        self._check_init_file()
        self._check_survey_year_file()

        init_params = self._read_initialization_config(init_file_path)
        init_params = self._set_params_from_init(source, bio_data_type, init_params, age_data_status)

        survey_year_params = self._read_survey_year_config(survey_year_file_path)

        self.params = self._collect_parameters(init_params, survey_year_params)

        self.params['exclude_age1'] = exclude_age1

        self.strata_df = None
        self.geo_strata_df = None
        self.strata_ds = None

        self.strata_class = LoadStrataData(stratification_index,
                                           self.params['data_root_dir'],
                                           self.params['filename_strata'],
                                           self.params['stratification_filename'])

        self.length_df = None
        self.length_ds = None

        self.specimen_df = None

        self.nasc_df = None

        self.final_biomass_table = None

        self._load_files(stratification_index)

        self._compute_biomass_density()

    def _check_init_file(self):
        # TODO: create this function that checks the contents of the initialization config file
        # TODO: it should make sure that certain variables are defined too
        print("A check of the initialization file needs to be done!")

    def _check_survey_year_file(self):
        # TODO: create this function that checks the contents of the survey year config file
        # TODO: it should make sure that certain variables are defined and all paths exist
        print("A check of the survey year file needs to be done!")

    @staticmethod
    def _read_initialization_config(init_file_path):

        with open(init_file_path) as f:

            init_params = yaml.load(f, Loader=yaml.SafeLoader)

        return init_params

    @staticmethod
    def _read_survey_year_config(survey_year_file_path):
        # TODO: This may need to be replaced in the future with a python file that way we can
        # TODO: mirror proc_parameters_XXXX.m more closely (specifically the switch case statement)

        with open(survey_year_file_path) as f:

            survey_params = yaml.load(f, Loader=yaml.SafeLoader)

        return survey_params

    @staticmethod
    def _set_params_from_init(source: int, bio_data_type: int, init_params: dict, age_data_status: int):

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

    def _collect_parameters(self, init_params, survey_params):

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

    def _load_nasc_data(self):
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

        # TODO: Should we create a class for this?

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

        return df

    def _load_files(self, stratification_index):
        """
        Load the biological, NASC table, stratification file,
        Constructs final trawl and catch tables

        Parameters
        ----------
        stratification_index : int
            Index for the chosen stratification
            0 = INPFC strata
            1 = KS (trawl)-based
            2-6 = geographically based but close to trawl-based stratification
            7 = mix-proportion, rather than 85% & 20% hake/hake-mix rules
            10 = one stratum for the whole survey

        Returns
        -------

        """

        # load specimen and length data using EchoPro variables
        LoadBioData(self)

        # load all associated stratification data
        self.strata_df, self.strata_ds, self.geo_strata_df = self.strata_class.get_strata_data(self.specimen_df,
                                                                                               self.length_df)

        self.nasc_df = self._load_nasc_data()

    def _compute_biomass_density(self):
        """
        Computes the biomass density estimate based on
        the nasc_df, strata_df, strata_ds, specimen_df,
        and length_df inputs. Additionally, it produces
        the final_biomass_table that is used by
        downstream processes.
        """

        bio_dense = ComputeBiomassDensity(self)

        bio_dense.get_final_biomass_table()

    def run_cv_analysis(self, lat_INPFC=None, kriged_data=False, seed=None):

        if self.params["JH_fac"] == 1:
            nr = 1  # number of realizations
        else:
            nr = 10000  # number of realizations

        cva = CVAnalysis(self)

        if kriged_data:
            raise NotImplementedError("CV analysis for kriged data has not been implemented")
        else:
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
