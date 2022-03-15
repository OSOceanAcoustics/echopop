import yaml
import numpy as np
# from .run_bootstrapping import RunBootstrapping
import pandas as pd


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
    """

    def __init__(self,
                 init_file_path: str,
                 survey_year_file_path: str,
                 source: int = 3,
                 bio_data_type: int = 1,
                 age_data_status: int = 1):

        self.bootstrapping_performed = False

        # self.__init_file_path = init_file_path
        self.__check_init_file()

        # self.__survey_year_file_path = survey_year_file_path
        self.__check_survey_year_file()

        init_params = self.__read_initialization_config(init_file_path)
        init_params = self.__set_params_from_init(source, bio_data_type, init_params, age_data_status)

        survey_year_params = self.__read_survey_year_config(survey_year_file_path)

        self.params = self.__collect_parameters(init_params, survey_year_params)

        self.__load_files()

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
                                                      dtype=np.float64)

        # setting bio_hake_age_bin variable to a numpy array
        init_params["bio_hake_age_bin"] = np.linspace(init_params["bio_hake_age_bin"][0],
                                                      init_params["bio_hake_age_bin"][1],
                                                      num=init_params["bio_hake_age_bin"][2],
                                                      dtype=np.float64)

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

    def process_length_data(self, df, haul_num_offset):
        """
        Parameters
        ----------
        df : Pandas Dataframe
            Dataframe holding the length data
        haul_num_offset : int
            # TODO: what does this refer to? Should we account for this in the code?
        Returns
        -------
        Processed Dataframe
        """

        # obtaining those columns that are required
        df = df[['Haul', 'Species_Code', 'Sex', 'Length', 'Frequency']].copy()

        # extract target species
        df = df.loc[df['Species_Code'] == self.params['species_code_ID']]

        # Apply haul_num_offset
        df['Haul'] = df['Haul'] + haul_num_offset

        if self.params['exclude_age1'] == 1:
            raise NotImplementedError("Excluding age 1 data has not been implemented!")

        df.drop(columns=['Species_Code'], inplace=True)

        # expand length and sex columns
        df['Frequency'] = df['Frequency'].apply(lambda x: np.ones(x))
        df['Length'] = df['Length'] * df['Frequency']
        df['Sex'] = df['Sex'] * df['Frequency']
        df.drop(columns='Frequency', inplace=True)

        # Group by Haul and calculate necessary statistics
        df = df.groupby('Haul').agg(lambda x: np.concatenate(x.values))
        df['TS0'] = 20.0 * df['Length'].apply(lambda x: np.log10(x)) - 68.0
        df = pd.concat([df, df.apply(lambda x: pd.Series([np.where(x.Sex == 1)[0], np.where(x.Sex == 2)[0],
                                                          len(x.Length), np.nanmean(x.Length),
                                                          np.nanstd(x.Length, ddof=1),
                                                          10.0 * np.log10(np.nanmean(10.0 ** (x.TS0 / 10.0))),
                                                          np.std(x.TS0, ddof=1), np.mean(x.TS0)],
                                                         index=['Male_ind', 'Female_ind', 'n', 'meanlen', 'stdlen',
                                                                'TS_lin', 'TS_sd', 'TS_log']),
                                     result_type='expand', axis=1)], axis=1)
        df.drop(columns='TS0', inplace=True)

        # get number of female and males
        df[['nM', 'nF']] = df.apply(lambda x: pd.Series([len(x.Male_ind), len(x.Female_ind)]), result_type='expand',
                                    axis=1)

        return df

    @staticmethod
    def process_catch_data(df):
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

        df.set_index('Haul', inplace=True)

        df.sort_index(inplace=True)

        return df

    def process_trawl_data(self, df):
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

        df['Duration'] = pd.to_timedelta(df['Duration']) # change datatype from string to timedelta

        if self.params['exclude_age1'] == 1:
            raise NotImplementedError("Excluding age 1 data has not been implemented!")

        # TODO: Should these two lines be added or assumed?
        # ind=find(diff(out.trawl_no)< -10);
        # out.trawl_no(ind + 1: end)=out.trawl_no(ind + 1: end)+para.bio.haul_no_offset;

        if self.params['hemisphere'][0] == 'N':
            df['EQ_Latitude'] = df['EQ_Latitude'].abs()
        elif self.params['hemisphere'][0] == 'S':
            df['EQ_Latitude'] = -1.0 * (df['EQ_Latitude'].abs())
        else:
            raise ValueError('Wrong N/S Hemisphere provided! \n')

        if self.params['hemisphere'][1] == 'W':
            df['EQ_Longitude'] = -1.0 * (df['EQ_Longitude'].abs())
        elif self.params['hemisphere'][1] == 'E':
            df['EQ_Longitude'] = df['EQ_Longitude'].abs()
        else:
            raise ValueError('Wrong E/W Hemisphere provided! \n')

        return df

    def process_gear_data(self, df):
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

        if self.params['exclude_age1'] == 1:
            raise NotImplementedError("Excluding age 1 data has not been implemented!")

        # TODO: Are these lines necessary? We should maybe take this into account during data construction.
        # if str2num(char(para.survey_year))  ~= 2011 & str2num(char(para.survey_year))  ~= 2015 & str2num(char(para.survey_year))  ~= 1995
        #     ind = find(diff(out.transect) < -10);
        #     out.transect(ind + 1: end)=out.transect(ind + 1: end)+out.transect(ind);
        # ind = 1:length(out.transect);
        # out.transect = out.transect(ind);
        # ind = find(diff(out.trawl_no) < -10);
        # out.trawl_no(ind + 1: end)=out.trawl_no(ind + 1: end)+para.bio.haul_no_offset;

        df['Net_Height'] = df['Net_Height'].apply(lambda x: np.nan if type(x) == str else x).astype(float)

        return df

    def process_specimen_data(self, df, haul_num_offset):
        """
        Parameters
        ----------
        df : Pandas Dataframe
            Dataframe holding the specimen data
        haul_num_offset : int
            # TODO: what does this refer to? Should we account for this in the code?
        Returns
        -------
        Processed Dataframe
        """

        df = df[['Haul', 'Species_Code', 'Sex', 'Length', 'Weight', 'Specimen_Number', 'Age']].copy()

        # extract target species
        df = df.loc[df['Species_Code'] == self.params['species_code_ID']]

        # Apply haul_num_offset
        df['Haul'] = df['Haul'] + haul_num_offset

        if self.params['exclude_age1'] == 1:
            raise NotImplementedError("Excluding age 1 data has not been implemented!")

        df.drop(columns=['Species_Code'], inplace=True)

        if self.params['age_data_status'] == 2:
            raise NotImplementedError("age_data_status = 2 has not been implemented")
        elif self.params['age_data_status'] != 1:
            raise NotImplementedError(
                f"age_data_status = {self.params['age_data_status']} has not been implemented")

        df.set_index('Haul', inplace=True)
        df.sort_index(inplace=True)

        if len(df['Age']) - df['Age'].isna().sum() < 0.1 * len(df['Age']):
            raise RuntimeWarning('Aged data are less than 10%!\n')

        return df

    def load_biological_data_us(self):
        """
        A file that loads all biological data provided for the US section.

        Returns
        -------
        catch_us_df : Pandas Dataframe
            Dataframe containing the catch data for the US
        length_us_df  : Pandas Dataframe
            Dataframe containing the length data for the US
        trawl_us_df : Pandas Dataframe
            Dataframe containing the trawl data for the US
        gear_us_df : Pandas Dataframe
            Dataframe containing the gear data for the US
        specimen_us_df : Pandas Dataframe
            Dataframe containing the specimen data for the US
        """

        if self.params['database_type'] == 'Oracle':
            catch_us_df = pd.read_excel(self.params['data_root_dir'] + self.params['filename_catch_US'])
            catch_us_df = self.process_catch_data(catch_us_df)

            length_us_df = pd.read_excel(self.params['data_root_dir'] + self.params['filename_length_US'])
            length_us_df = self.process_length_data(length_us_df, 0)
        else:
            catch_us_df = None
            length_us_df = None
            raise NotImplementedError("Loading data from a non-Oracle database has not been implemented!")

        if self.params['filename_trawl_US']:
            trawl_us_df = pd.read_excel(self.params['data_root_dir'] + self.params['filename_trawl_US'],
                                        converters= {'Duration': str})
            trawl_us_df = self.process_trawl_data(trawl_us_df)
        else:
            trawl_us_df = None

        if self.params['filename_gear_US']:
            gear_us_df = pd.read_excel(self.params['data_root_dir'] + self.params['filename_gear_US'])
            gear_us_df = self.process_gear_data(gear_us_df)
        else:
            gear_us_df = None

        specimen_us_df = pd.read_excel(self.params['data_root_dir'] + self.params['filename_specimen_US'])
        specimen_us_df = self.process_specimen_data(specimen_us_df, 0)

        return catch_us_df, length_us_df, trawl_us_df, gear_us_df, specimen_us_df

    def load_biological_data_canada(self):
        """
        A file that loads all biological data provided for the Canada section.

        Returns
        -------
        catch_can_df : Pandas Dataframe
            Dataframe containing the catch data for the Canada
        length_can_df  : Pandas Dataframe
            Dataframe containing the length data for the Canada
        trawl_can_df : Pandas Dataframe
            Dataframe containing the trawl data for the Canada
        gear_can_df : Pandas Dataframe
            Dataframe containing the gear data for the Canada
        specimen_can_df : Pandas Dataframe
            Dataframe containing the specimen data for the Canada
        """

        if self.params['database_type'] == 'Oracle':
            catch_can_df = pd.read_excel(self.params['data_root_dir'] + self.params['filename_catch_CAN'])
            catch_can_df = self.process_catch_data(catch_can_df)

            length_can_df = pd.read_excel(self.params['data_root_dir'] + self.params['filename_length_CAN'])
            length_can_df = self.process_length_data(length_can_df, self.params['haul_no_offset'])
        else:
            catch_can_df = None
            length_can_df = None
            raise NotImplementedError("Loading data from a non-Oracle database has not been implemented!")

        specimen_can_df = pd.read_excel(self.params['data_root_dir'] + self.params['filename_specimen_CAN'])
        specimen_can_df = self.process_specimen_data(specimen_can_df, self.params['haul_no_offset'])

        if self.params['filename_trawl_CAN']:
            trawl_can_df = pd.read_excel(self.params['data_root_dir'] + self.params['filename_trawl_CAN'],
                                         converters= {'Duration': str})
            trawl_can_df = self.process_trawl_data(trawl_can_df)
        else:
            trawl_can_df = None

        if self.params['filename_gear_CAN']:
            gear_can_df = pd.read_excel(self.params['data_root_dir'] + self.params['filename_gear_CAN'])
            gear_can_df = self.process_gear_data(gear_can_df)
        else:
            gear_can_df = None

        return catch_can_df, length_can_df, trawl_can_df, gear_can_df, specimen_can_df

    def __load_biological_data(self):
        """
        A function that loads the appropriate data based on the source provided.

        Returns
        -------
        self.length_df : Dataframe
            A dataframe containing the length information to be used downstream. If source=3 it contains US and Canada.
        self.specimen_df : Dataframe
            A dataframe containing the specimen information to be used downstream. If source=3 it contains US and Canada.
        self.trawl_df : Dataframe
            A dataframe containing the trawl information to be used downstream. If source=3 it contains US and Canada.
        self.gear_df : Dataframe
            A dataframe containing the gear information to be used downstream. If source=3 it contains US and Canada.
        self.catch_df : Dataframe
            A dataframe containing the catch information to be used downstream. If source=3 it contains US and Canada.
        """
        if self.params['source'] == 3:

            print("Loading US biological data ...")

            catch_us_df, length_us_df, trawl_us_df, gear_us_df, specimen_us_df = self.load_biological_data_us()
            catch_can_df, length_can_df, trawl_can_df, gear_can_df, specimen_can_df = self.load_biological_data_canada()

            # combine specimen and length dataframes
            self.length_df = pd.concat([length_us_df, length_can_df])
            self.specimen_df = pd.concat([specimen_us_df, specimen_can_df])

            # deal with trawl number offset
            max_US_haul_no = trawl_us_df['Haul'].max()
            if max_US_haul_no > self.params['haul_no_offset']:
                CAN_haul_no_offset = 100 * np.ceil(max_US_haul_no / 100.0)
            else:
                CAN_haul_no_offset = self.params['haul_no_offset']

            # max_US_Transect = max(data.bio.gear.transect) # TODO: This looks unused, make sure this is the case

            # combine US & CAN trawl files
            trawl_can_df['Haul'] = trawl_can_df['Haul'] + CAN_haul_no_offset  # add haul_offset
            # change lon to -lon
            trawl_can_df['EQ_Longitude'] = trawl_can_df['EQ_Longitude'].apply(lambda x: -1.0 * x if x > 0 else x)
            self.trawl_df = pd.concat([trawl_us_df, trawl_can_df])

            # transect offset is set to zero since we will not have overlap transects
            # TODO: Is this a necessary variable? If so, we should make it an input to the function.
            CAN_Transect_offset = 0

            # combine US & CAN gear files
            gear_can_df['Haul'] = gear_can_df['Haul'] + CAN_haul_no_offset  # add haul_offset
            gear_can_df['Transect'] = gear_can_df['Transect'] + CAN_Transect_offset
            self.gear_df = pd.concat([gear_us_df, gear_can_df])

            # combine US & CAN catch files
            catch_can_df.index = catch_can_df.index + CAN_haul_no_offset  # add haul_offset
            self.catch_df = pd.concat([catch_us_df, catch_can_df])

        else:

            raise NotImplementedError(f"Source of {self.params['source']} not implemented yet.")

    def __load_files(self):

        self.__load_biological_data()

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

    def run_process(self,
                          extrapolation: bool = False,
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

        # import copy
        #
        # # get copy of EchoPro
        # # TODO: This is creating a copy of the EchoPro object that
        # # TODO: could eat up memory if the input files are large
        # # epro_copy = copy.deepcopy(self)
        # #
        # # print(self)
        #
        # bootstrapping_routine = RunBootstrapping(**function_args)
        #
        # bootstrapping_routine.set_echopro_object(self)

        return



