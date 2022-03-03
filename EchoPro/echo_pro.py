import yaml
import numpy as np
from .run_bootstrapping import RunBootstrapping


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
    """

    def __init__(self,
                 init_file_path: str,
                 survey_year_file_path: str):

        self.bootstrapping_performed = False

        self.__init_file_path = init_file_path
        self.__check_init_file()

        self.__survey_year_file_path = survey_year_file_path
        self.__check_survey_year_file()

        self.__read_initialization_config()
        self.__set_params_from_init()

        self.__read_survey_year_config()

        self.__load_files()

    def __check_init_file(self):
        # TODO: create this function that checks the contents of the initialization config file
        # TODO: it should make sure that certain variables are defined too
        print("A check of the initialization file needs to be done!")

    def __check_survey_year_file(self):
        # TODO: create this function that checks the contents of the survey year config file
        # TODO: it should make sure that certain variables are defined and all paths exist
        print("A check of the survey year file needs to be done!")

    def __read_initialization_config(self):

        with open(self.__init_file_path) as f:

            self.__init_params = yaml.load(f, Loader=yaml.SafeLoader)

    def __read_survey_year_config(self):
        # TODO: This may need to be replaced in the future with a python file that way we can
        # TODO: mirror proc_parameters_XXXX.m more closely (specifically the switch case statement)

        with open(self.__survey_year_file_path) as f:

            self.__survey_params = yaml.load(f, Loader=yaml.SafeLoader)

    def __set_params_from_init(self):

        # setting bio_hake_lin_bin variable to a numpy array
        self.__init_params["bio_hake_len_bin"] = np.linspace(self.__init_params["bio_hake_len_bin"][0],
                                                             self.__init_params["bio_hake_len_bin"][1],
                                                             num=self.__init_params["bio_hake_len_bin"][2],
                                                             dtype=np.float64)

        # setting bio_hake_age_bin variable to a numpy array
        self.__init_params["bio_hake_age_bin"] = np.linspace(self.__init_params["bio_hake_age_bin"][0],
                                                             self.__init_params["bio_hake_age_bin"][1],
                                                             num=self.__init_params["bio_hake_age_bin"][2],
                                                             dtype=np.float64)

        # making acoust_freq0 into a numpy array
        self.__init_params["acoust_freq0"] = np.array(self.__init_params["acoust_freq0"], dtype=np.float64)

        # turning beamwidth input into radians
        self.__init_params["acoust_bw"] = np.array(self.__init_params["acoust_bw"], dtype=np.float64)*np.pi/180.0

        # finding the acoustic freqency index, default = 1 --> 38 kHz
        # TODO: make sure zero indexing doesn't mess up downsteam items
        self.__init_params.update({'acoust_freq_ind': np.intersect1d(np.floor(self.__init_params["acoust_freq"]/1000.0),
                                                                     self.__init_params["acoust_freq0"],
                                                                     return_indices=True)[2]})

        # Setting acoust_sig_b_coef to a float from string input
        self.__init_params["acoust_sig_b_coef"] = eval(self.__init_params["acoust_sig_b_coef"])

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

    def __load_files(self):

        print("Loading files")

    def run_bootstrapping(self,
                          extrapolation: bool = False,
                          age_data_status: int = 1,
                          source: int = 3,
                          bio_data_type: int = 1,
                          KS_stratification: int = 1,
                          kriging: bool = True,
                          kriging_input: int = 1,
                          exclude_age1: bool = False,
                          start_transect: int = 1,
                          end_transect: int = 200,
                          transect_reduction_fraction: float = 0.0,
                          transect_reduction_mode: int = 1,
                          bootstrap_limit: int = 1,
                          default_parameters: bool = True):
        """
        A function that performs bootstrapping. This involves processing
        acoustic and biological data, computation of CV analysis,
        and Kriging.

        Parameters
        ----------
        extrapolation : bool
                Specifies if extrapolation should be used for Kriging
        age_data_status : int
            0 = fake age data (TODO: this option looks unused)
            1 = actual age data
            2 = from age_vs_len
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
        KS_stratification : int
            Specifies the type of stratification to be used.
            0 = Pre-Stratification or customized stratification (geographically defined)
            1 = Post-Stratification (KS-based or trawl-based)
        kriging : bool
            States whether or not to perform kriging
        kriging_input : int
            Specifies TODO: ask Chu what this variable actually specifies
            1 = Biomass density
            2 = NASC
            3 = Number density
        exclude_age1 : bool
            States whether or not age 1 hake should be included in analysis.
        start_transect : int
            Value to start transect
        end_transect : int
            Value to end transect
        transect_reduction_fraction : float
            Reduction fraction for transect TODO: should this be 5 or 0.05?
        transect_reduction_mode : int
            1 = Regular
            2 = Random
        bootstrap_limit : int
            The number of bootstraping iterations to perform
        default_parameters : bool
            States whether or not to use the default parameters
        """

        # Get all inputs to the function run_bootstrapping()
        function_args = locals()

        # remove the self argument
        del function_args['self']

        print(f"saved args = {function_args}")

        import copy

        # get copy of EchoPro
        # TODO: This is creating a copy of the EchoPro object that
        # TODO: could eat up memory if the input files are large
        # epro_copy = copy.deepcopy(self)
        #
        # print(self)

        bootstrapping_routine = RunBootstrapping(**function_args)

        bootstrapping_routine.set_echopro_object(self)

        return



