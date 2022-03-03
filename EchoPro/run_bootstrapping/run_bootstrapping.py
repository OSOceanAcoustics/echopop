

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
    source : int
        Define the region of data to use.
    bio_data_type : int
        Specifies the biological data to be used.
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

        # # TODO: eventually bulk the below functions
        #
        # # set variables based on the source used
        # if source <= 3:
        #     self.__init_params["platform_name"] = 'FSV'
        # else:  # TODO look into this and see if this else statement is necessary
        #     self.__init_params["platform_name"] = 'SD'
        #     if bio_data_type != 3:
        #         bio_data_type = 3
        #         print("Changing bio_data_type to 3 based on platform name.")
        #
        # self.__init_params["opr_indx"] = 3  # TODO: look into this variable, might only be necessary for Matlab GUI
        #
        # # setting the species code ID based on bio data type
        # if bio_data_type == 1:
        #     self.__init_params["species_code_ID"] = 22500  # target species_code for acoustic survey
        # elif bio_data_type == 2:
        #     self.__init_params["species_code_ID"] = 22500  # target species_code for bottom trawl survey
        # elif bio_data_type == 3:
        #     self.__init_params["species_code_ID"] = 206  # target species_code for industry data(observer data)
        #
        # # setting the stratification index based on user provided input
        # # TODO: Might be able to take out this if else statement depending on downstream items
        # if KS_stratification == 1:
        #     stratification_index = 1
        # else:
        #     stratification_index = 0
        #
        # # check that the stratification index is correctly set
        # if bio_data_type != 1:
        #     if stratification_index != 0:
        #         stratification_index = 0  # non - acoustical and trawl survey data only use INPFC stratification
        #         print("Changing stratification_index to 0.")
        # else:
        #     if stratification_index != 1:
        #         print("Changing stratification_index to 1")
        #         stratification_index = 1  # index for the chosen stratification
        #         # 1 = KS(trawl) - based, 2 - 7 = geographically based but close to trawl - based stratification
        #         # 0 = INPFC strata
        #         # 7 = mix - proportion, rather than 85 % & 20 % hake / hake - mix rules
        #         # 10 = one stratum for the whole survey
        #         # TODO: ask about the above comments
        #
        # # TODO: make sure to take this into account!
        # # if para.proc.exclude_age1 == 1
        # #     para.acoust.filename.processed_data = para.acoust.filename.processed_data_age2;
        # #     para.bio_acoust.filename.Transect_region_haul = para.bio_acoust.filename.Transect_region_haul_age2;
        # # else
        # #     para.acoust.filename.processed_data = para.acoust.filename.processed_data_age1;
        # #     para.bio_acoust.filename.Transect_region_haul = para.bio_acoust.filename.Transect_region_haul_age1;
        # # end
        #
        # # check to make sure no survey year and initialization parameters are the same
        # param_intersect = set(self.__init_params.keys()).intersection(set(self.__survey_params.keys()))
        #
        # # if no parameters are the same, then run process, else return error
        # if not param_intersect:
        #     # combine survey year and initialization parameters into one dictionary
        #     full_params = {}
        #     full_params.update(self.__init_params)
        #     full_params.update(self.__survey_params)
        #
        #     return ProcessData(full_params, extrapolation, age_data_status, source, bio_data_type, KS_stratification,
        #                        stratification_index, kriging, kriging_input, exclude_age1, start_transect,
        #                        end_transect, transect_reduction_fraction, transect_reduction_mode, bootstrap_limit,
        #                        default_parameters)
        # else:
        #     raise RuntimeError('The initialization and survey year configuration files define the same variable! ' +
        #                        f'\n These variables are: {param_intersect}')

