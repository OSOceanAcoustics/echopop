

class ProcessData:
    """
    EchoPro class for conducting the processing of the provided parameters.

    Parameters
    ----------
    params : dict
        A dictionary containing all parameters from the survey year and
        initialization configuration files
    extrapolation : bool
        Specifies if extrapolation should be used for Kriging
    age_data_status : int
        Specifies the type of age data provided
    source : int
        Define the region of data to use.
    bio_data_type : int
        Specifies the biological data to be used.
    KS_stratification : int
        Specifies the type of stratification to be used.
    stratification_index : int
        Index used for later downstream processing
    kriging : bool
            States whether or not to perform kriging
    kriging_input : int
        Specifies TODO: ask Chu what this variable actually specifies
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
    bootstrap_limi: int
        The number of bootstraping iterations to perform
    default_parameters : bool
        States whether or not to use the default parameters
    """

    def __init__(self,
                 params: dict,
                 extrapolation: bool,
                 age_data_status: int,
                 source: int,
                 bio_data_type: int,
                 KS_stratification: int,
                 stratification_index: int,
                 kriging: bool,
                 kriging_input: int,
                 exclude_age1: bool,
                 start_transect: int,
                 end_transect: int,
                 transect_reduction_fraction: float,
                 transect_reduction_mode: int,
                 bootstrap_limit: int,
                 default_parameters: bool):

        self.__params = params
        self.__extrapolation = extrapolation
        self.__age_data_status = age_data_status
        self.__source = source
        self.__bio_data_type = bio_data_type
        self.__KS_stratification = KS_stratification
        self.__stratification_index = stratification_index
        self.__kriging = kriging
        self.__kriging_input = kriging_input
        self.__exclude_age1 = exclude_age1
        self.__start_transect = start_transect
        self.__end_transect = end_transect
        self.__transect_reduction_fraction = transect_reduction_fraction
        self.__transect_reduction_mode = transect_reduction_mode
        self.__bootstrap_limit = bootstrap_limit
        self.__default_parameters = default_parameters

        self.__process_acoustic_data()


    def __process_nasc_data(self):
        """
        Process echo-integration NASC data to obtain length-age-sex structured biomass estimate and all
        required and useful tables & results
        """
        print(f'\n ================== Biomass Estimate of {self.__params["survey_year"]} Data ===================\n')





    def __process_raw_acoustic_data(self):
        """
        Process raw acoustic echogram data, either EK60 or EK500.
        """
        print("Needs to be implemented!")

    def __process_acoustic_data(self):

        if self.__params["acoust_file_type"] == 1:
            self.__process_nasc_data()
        else:
            self.__process_raw_acoustic_data()
