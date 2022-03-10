

class RunBootstrapping:
    """
    A Class that performs bootstrapping. This involves processing
    acoustic and biological data, computation of CV analysis,
    and Kriging.

    Parameters
    ----------
    echopro : EchoPro
        EchoPro object that contains all initial parameters
    extrapolation : bool
        Specifies if extrapolation should be used for Kriging
    age_data_status : int
        Specifies the type of age data provided
    KS_stratification : int
        Specifies the type of stratification to be used.
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
    bootstrap_limit: int
        The number of bootstraping iterations to perform
    default_parameters : bool
        States whether or not to use the default parameters
    """

    def __init__(self,
                 extrapolation: bool,
                 age_data_status: int,
                 source: int,
                 bio_data_type: int,
                 KS_stratification: int,
                 kriging: bool,
                 kriging_input: int,
                 exclude_age1: bool,
                 start_transect: int,
                 end_transect: int,
                 transect_reduction_fraction: float,
                 transect_reduction_mode: int,
                 bootstrap_limit: int,
                 default_parameters: bool):

        self.class_input_args = locals()

        # remove the self argument
        del self.class_input_args['self']

        print("Performing bootstrapping")

        print(self.class_input_args)

        self.epro = None

        # self.__get_init_parameters()

    def set_echopro_object(self, epro):

        self.epro = epro

        print(self.epro)


    def __get_init_parameters(self):

        print("function to initialize parameters")

