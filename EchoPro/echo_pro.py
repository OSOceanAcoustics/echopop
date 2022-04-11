import sys

import yaml
import numpy as np
# from .run_bootstrapping import RunBootstrapping
from .load_biological_data import LoadBioData
import pandas as pd
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

            # TODO: is the below code necessary?
            # [n, m] = size(dat);
            # out(1: n, 1)=dat(:, 1); % transect number
            # ind = find(out(1:n, 1) > 1000); % 1000 is a fixed number added to the first transect of the CAN survey transect
            # if ~isempty(ind)
            #     out(ind, 1) = out(ind, 1) - 1000 + transect_offset; % modify the first transect line number to Tnum + offset
            #     out(ind + 1: end, 1)=out(ind + 1: end, 1)+transect_offset; % modify the rest CAN transect line numbers to Tnum + offset
            # end
            # out(1: n, 2)=999 * ones(n, 1); % region number - not useful
            # out(1: n, 3)=dat(:, 2); % VL start
            # out(1: n, 4)=dat(:, 2)+0.5; % VL stop
            # out(1: n, 5: m + 2)=dat(:, 3: m); %
            # ind = find(out(1:n, 6) > 0  ); % convert longitude to negative
            # out(ind, 6) = -out(ind, 6);
            # if str2num(survey_year) == 1995
            #     ind_nan = find(isnan(out(:, 7)) == 1);
            #     ind_good = find(isnan(out(:, 7)) == 0);
            #     out(ind_nan, 7) = floor(interp1(ind_good, out(ind_good, 7), ind_nan));
            #     ind7 = find(out(:, 7) == 7);
            #     ind7_lt5000 = ind7(ind7 < 5000);
            #     ind7_gt5000 = ind7(ind7 >= 5000);
            #     out(ind7_lt5000, 7) = 6;
            #     out(ind7_gt5000, 7) = 8;
            # end

            raise NotImplementedError("Loading the NASC table for survey years less than 2003 has not been implemented!")

        else:
            df.set_index('Transect', inplace=True)
            df.sort_index(inplace=True)

        return df

    def __load_stratafication_file(self, stratification_index):
        """
        Loads stratification file. Produces the variable self.strata_df,
        which is a Pandas dataframe representing stratification file.

        Parameters
        ----------
        stratification_index : int  # TODO: this looks to mirror KS_stratification, it seems to be unnecessary to use this
                                    # TODO: Ask Chu if it is ok to remove this.
            Index for the chosen stratification
            0 = INPFC strata
            1 = KS (trawl)-based
            2-6 = geographically based but close to trawl-based stratification
            7 = mix-proportion, rather than 85% & 20% hake/hake-mix rules
            10 = one stratum for the whole survey
        """

        # load stratification file
        if stratification_index != 1 and stratification_index != 0:
            raise NotImplementedError(f"stratification_index of {stratification_index} has not been implemented!")
        else:

            if stratification_index == 1:
                self.strata_df = pd.read_excel(self.params['data_root_dir'] + self.params['filename_strata'],
                                               sheet_name='Base KS')
                self.strata_df = self.strata_df[['Year', 'Cluster name', 'Haul', 'wt']].copy()

                # set data types of dataframe
                self.strata_df = self.strata_df.astype({'Year': int, 'Cluster name': int, 'Haul': int,
                                                        'wt': np.float64})

                self.strata_df.set_index('Haul', inplace=True)
                self.strata_df.sort_index(inplace=True)

            else:
                self.strata_df = pd.read_excel(self.params['data_root_dir'] + self.params['filename_strata'],
                                               sheet_name='INPFC')
                self.strata_df = self.strata_df[['Year', 'INPFC', 'Haul', 'wt']].copy()

                # set data types of dataframe
                self.strata_df = self.strata_df.astype({'Year': int, 'INPFC': int, 'Haul': int, 'wt': np.float64})

                self.strata_df.set_index('Haul', inplace=True)
                self.strata_df.sort_index(inplace=True)

    def __get_bio_strata(self, stratification_index, KS_stratification, transect_reduction_fraction):
        """
        A function that obtains a subset of the strata data corresponding
        to the biological data. This information is stored in the
        Pandas Dataframe ``self.bio_strata``.
        """

        if stratification_index == 1:

            if self.params['CAN_strata_num0']:
                raise NotImplementedError("Filled CAN_strata_num0 value has not been implemented!")

            if self.params['exclude_age1']:

                # TODO: Do we need to have a condition to check for para.proc.age1_haul?
                # TODO: According to Chu, this seems to be unused.
                # raise NotImplementedError(f"age1_haul")

                if transect_reduction_fraction != 0.0:
                    raise NotImplementedError("transect reduction has not been implemented!")
                else:

                    # select only stratum not equal to zero
                    strata_mask = self.strata_df['Cluster name'] != 0
                    self.bio_strata = self.strata_df[strata_mask][
                        ['Cluster name', 'wt']].reset_index().set_index('Cluster name')
        else:

            if self.params['CAN_strata_num0']:
                raise NotImplementedError("Filled CAN_strata_num0 value has not been implemented!")

            if self.params['exclude_age1']:

                # TODO: Do we need to have a condition to check for para.proc.age1_haul?
                # TODO: According to Chu, this seems to be unused.
                # raise NotImplementedError(f"age1_haul")

                if transect_reduction_fraction != 0.0:
                    raise NotImplementedError("transect reduction has not been implemented!")
                else:

                    # select only stratum not equal to zero
                    strata_mask = self.strata_df['INPFC'] != 0
                    self.bio_strata = self.strata_df[strata_mask][['INPFC', 'wt']].reset_index().set_index('INPFC')

        # TODO: Do we have to fill-in trawl information for stratum number < first non-empty stratum?
        #  Matlab get_historical_strata_data lines 144-169

        # TODO: Do we have to interpolate trawl information to stranum (no tralws) by using the parameters
        #  of the adjacient strata (stratum)
        #  Matlab get_historical_strata_data lines 172-194

        # TODO: determine if this code is necessary.
        # # extrpolate trawl information to stranum (no tralws) > last non-empty stratum (situation in using the Observer data or A-SHOP data)
        # # automatically determine the maximum number of stratum   11/2/2021 (SD, and/or Observer uses geographic_stratification xlsx file)
        # if para.proc.KS_stratification == 1
        #     n_strata_max =  max(dat0(:,2));     # change max(dat0(:,1)) to max(dat0(:,2)) on 5/29/2021
        #                                         # change back from max(dat0(:,2)) to max(dat0(:,1)) on 10/18/2021
        #                                         # change to if loop for determining what which column corresponding strata index
        # else
        #     n_strata_max =  max(dat0(:,1));
        #     if n_strata_max > 20
        #         n_strata_max =  max(dat0(:,2));
        #     end
        # end

        # for i = data.bio.haul_strata(end)+1:n_strata_max   # change max(dat0(:,1)) to max(dat0(:,2)) on 5/29/2021
        #                                                    # change back from max(dat0(:,2)) to max(dat0(:,1)) on 10/18/2021
        #     data.bio.strata(i).trawls=data.bio.strata(data.bio.haul_strata(end)).trawls;
        #     data.bio.strata(i).wgt=data.bio.strata(data.bio.haul_strata(end)).wgt;
        #     data.bio.haul_strata(i) = i;
        #     if length(data.bio.strata(i).trawls) < 1
        #         disp(data.bio.strata(i).trawls)
        #     end
        # end

    def __get_expanded_data(self, ind):
        """
        Collect the expanded length data from length_ds for
        male, female, and unsexed genders.

        Parameters
        ----------
        ind : int
            index of haul to select data from
        """
        # Collect expanded length data from length_ds
        if ind in self.length_ds.Haul:
            male_val_len = np.repeat(self.length_ds.Length.values,
                                     np.nan_to_num(
                                         self.length_ds.Frequency.sel(Haul=ind, Sex=1).values).astype(dtype=int),
                                     axis=0)

            female_val_len = np.repeat(self.length_ds.Length.values,
                                       np.nan_to_num(
                                           self.length_ds.Frequency.sel(Haul=ind, Sex=2).values).astype(dtype=int),
                                       axis=0)

            # TODO: uncomment below if we want unsexed to be included
            unsexed_val_len = np.repeat(self.length_ds.Length.values,
                                        np.nan_to_num(
                                            self.length_ds.Frequency.sel(Haul=ind, Sex=3).values).astype(dtype=int),
                                        axis=0)
        else:
            male_val_len = []
            female_val_len = []

            # TODO: uncomment below if we want unsexed to be included
            unsexed_val_len = []

        return male_val_len, female_val_len, unsexed_val_len

    def __get_len_data_specimen(self, gen, ind):
        """
        Collect length data from specimen_df for a specific gender.

        Parameters
        ----------
        gen : float
            Gender
        ind : int
            index of haul to select data from
        """
        selected_ind = (self.specimen_df['Sex'] == gen) & pd.notnull(self.specimen_df['Age'])
        if ind in selected_ind[selected_ind].index:
            val_specimen = self.specimen_df[selected_ind].loc[ind]['Length']

            if isinstance(val_specimen, np.float64):
                val_specimen = [val_specimen]
            else:
                val_specimen = val_specimen.values
        else:
            val_specimen = []

        return val_specimen

    def __get_bio_len_age_haul(self):

        # construct hake-haul number array
        self.params['bio_hake_trawl_num'] = np.union1d(self.length_ds.Haul.values,
                                                       self.specimen_df.index.unique().values)

        # initialize parameters with a numpy array
        len_bio_hake_len_bin = len(self.params['bio_hake_len_bin'])
        len_bio_hake_trawl_num = len(self.params['bio_hake_trawl_num'])

        # TODO: uncomment below if we want unsexed to be included
        self.params['bio_len_haul_U'] = np.zeros((len_bio_hake_len_bin, len_bio_hake_trawl_num))
        self.params['bio_aged_len_haul_U'] = np.zeros((len_bio_hake_len_bin, len_bio_hake_trawl_num))

        self.params['bio_len_haul_M'] = np.zeros((len_bio_hake_len_bin, len_bio_hake_trawl_num))
        self.params['bio_len_haul_F'] = np.zeros((len_bio_hake_len_bin, len_bio_hake_trawl_num))
        self.params['bio_aged_len_haul_M'] = np.zeros((len_bio_hake_len_bin, len_bio_hake_trawl_num))
        self.params['bio_aged_len_haul_F'] = np.zeros((len_bio_hake_len_bin, len_bio_hake_trawl_num))

        j = 0
        for i in self.params['bio_hake_trawl_num']:

            male_val_len, female_val_len, unsexed_val_len = self.__get_expanded_data(i)

            # Collect length data from specimen_df for males
            male_val_specimen = self.__get_len_data_specimen(1.0, i)

            # Collect length data from specimen_df for females
            female_val_specimen = self.__get_len_data_specimen(2.0, i)

            # TODO: uncomment below if we want unsexed to be included
            # Collect length data from specimen_df for unsexed
            unsexed_val_specimen = self.__get_len_data_specimen(3.0, i)

            # store length data from length_ds and specimen_df
            self.params['bio_len_haul_M'][:, j] = self.load_bio.get_bin_counts(
                np.concatenate([male_val_len, male_val_specimen]), self.params['bio_hake_len_bin'])

            self.params['bio_len_haul_F'][:, j] = self.load_bio.get_bin_counts(
                np.concatenate([female_val_len, female_val_specimen]), self.params['bio_hake_len_bin'])

            # TODO: uncomment below if we want unsexed to be included
            self.params['bio_len_haul_U'][:, j] = self.load_bio.get_bin_counts(
                np.concatenate([unsexed_val_len, unsexed_val_specimen]), self.params['bio_hake_len_bin'])

            self.params['bio_aged_len_haul_M'][:, j] = self.load_bio.get_bin_counts(male_val_specimen,
                                                                                    self.params['bio_hake_len_bin'])
            self.params['bio_aged_len_haul_F'][:, j] = self.load_bio.get_bin_counts(female_val_specimen,
                                                                                    self.params['bio_hake_len_bin'])

            # TODO: uncomment below if we want unsexed to be included
            self.params['bio_aged_len_haul_U'][:, j] = self.load_bio.get_bin_counts(unsexed_val_specimen,
                                                                                    self.params['bio_hake_len_bin'])
            j += 1

        # self.params['bio_len_haul_ALL'] =  self.params['bio_len_haul_M'] + self.params['bio_len_haul_F'] 
        # self.params['bio_aged_len_haul_ALL'] =  self.params['bio_aged_len_haul_M'] + self.params['bio_aged_len_haul_F'] 

        # TODO: uncomment below if we want unsexed to be included
        self.params['bio_len_haul_ALL'] = self.params['bio_len_haul_M'] + self.params['bio_len_haul_F'] + \
                                               self.params['bio_len_haul_U']

        # TODO: uncomment below if we want unsexed to be included
        self.params['bio_aged_len_haul_ALL'] = self.params['bio_aged_len_haul_M'] + self.params[
            'bio_aged_len_haul_F'] + self.params['bio_aged_len_haul_U']

    def __remove_empty_bio_len_age_haul(self):
        """
        Deletes columns in bio_len_age_haul arrays that are empty.
        """

        aged_del_ind = self.params['bio_aged_len_haul_ALL'].sum(axis=0)


    def get_historical_strata_data(self, stratification_index, KS_stratification,
                                   transect_reduction_fraction: float = 0.0):
        """
        Get quantities and keys associated with strata for:
        1. trawl information
        2. Length - key
        3. Length-weight key
        4. Length-age - key for abundance & biomass

        Parameters
        ----------
        stratification_index : int  # TODO: this looks to mirror KS_stratification, it seems to be unnecessary to use this
                                    # TODO: Ask Chu if it is ok to remove this.
            Index for the chosen stratification
            0 = INPFC strata
            1 = KS (trawl)-based
            2-6 = geographically based but close to trawl-based stratification
            7 = mix-proportion, rather than 85% & 20% hake/hake-mix rules
            10 = one stratum for the whole survey
        KS_stratification : int
            Specifies the type of stratification to be used.
            0 = Pre-Stratification or customized stratification (geographically defined)
            1 = Post-Stratification (KS-based or trawl-based)
        transect_reduction_fraction : float
            Reduction fraction for transect TODO: should this be 5 or 0.05?
        """

        # TODO: is it necessary to have this? It looks like it is not necessary.
        # if para.proc.stratification_index == 0 & para.proc.KS_stratification == 0 # INPFC with
        # if dat0(1, 1) == 2015
        #     # remove age1 hauls but not those used in transect_region_hual files, to be consistent
        #     with K - S stratification
        #     para.proc.age1_haul = [4 9 10 13 15 18 23 30 35 39 41 48 60 69 227];

        if self.params['bio_data_type'] != 3:

            # TODO: Is this code necessary? It looks like the file para.acoust.filename.strata will
            # TODO: only have one year
            # indx = find(str2num(para.survey_year ) == dat0(:,1));
            # dat = dat0(indx,[2 4 5]);                         % strata-haul-weight factor

            print("Do we need to set stratum_id or just use strata_df? Look into this!")
            # stratum_id = self.strata_df['Cluster name'] or self.strata_df['INPFC']

        else:
            # TODO: If this statement is triggered, then this is probably not the right place to
            # TODO: put this code because is could potential be ran in bootstrapping.
            raise NotImplementedError("Loading the stratification file has not been "
                                      + f"implemented for bio_data_type = {self.params['bio_data_type']}")

        # TODO: these lines look unnecessary for the Python version
        # data.bio_acoust.haul_wgt_tbl = dat(:, 2: 3);
        # data.bio.haul_strata = unique(dat(:, 1));

        # TODO: these lines look unnecessary for the Python version
        # % % remove strata index 0 - - YOY
        # ind = find(data.bio.haul_strata == 0);
        # data.bio.haul_strata(ind) = [];
        # n = length(data.bio.haul_strata); % number of strata

        self.__get_bio_strata(stratification_index, KS_stratification, transect_reduction_fraction)

        self.__get_bio_len_age_haul()

        # self.__remove_empty_bio_len_age_haul()

        return


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

        self.load_bio = LoadBioData(self)  # TODO: remove self from load_bio, only necessary for testing

        self.__load_nasc_data()

        if self.params['bio_data_type'] == 1:

            self.load_bio.process_length_weight_data(self.specimen_df)

            self.__load_stratafication_file(stratification_index)

            self.load_bio.get_final_catch_trawl_tables(stratification_index) # TODO: this might not be needed make it optional

        else:
            raise NotImplementedError(f"Processing bio_data_type = {self.params['bio_data_type']} has not been implemented!")


        # self.get_historical_strata_data(stratification_index, transect_reduction_fraction)

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



