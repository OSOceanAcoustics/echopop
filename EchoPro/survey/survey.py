import yaml
import numpy as np
# from .run_bootstrapping import RunBootstrapping
from ..load_biological_data import LoadBioData
from ..load_stratification_data import LoadStrataData
from ..compute_biomass_density import ComputeBiomassDensity
from ..cv_analysis import cv_analysis
from ..kriging import Kriging
from ..kriging_mesh import KrigingMesh
from ..semivariogram import SemiVariogram
from ..load_nasc_data import load_nasc_data
from typing import Tuple


class Survey:
    """
    EchoPro base class that imports and prepares parameters for
    a survey. Additionally, it includes functions for
    accessing the modules associated with the biomass density
    calculation, CV analysis, semi-variogram algorithm, and
    Kriging.

    Parameters
    ----------
    init_file_path : str
        A string specifying the path to the initialization YAML file
    survey_year_file_path : str
        A string specifying the path to the survey year YAML file
    source : int
        The region of data to use.
        1 = US
        2 = Canada
        3 = US and Canada
    exclude_age1 : bool
        States whether age 1 hake should be included in analysis.
    """
    def __init__(self,
                 init_file_path: str,
                 survey_year_file_path: str,
                 source: int = 3,
                 exclude_age1: bool = True):

        self._check_init_file(init_file_path)
        self._check_survey_year_file(survey_year_file_path)

        # read initialization configuration file
        init_params = self._read_config(init_file_path)
        init_params = self._set_params_from_init(source, init_params)

        # read survey year configuration file
        survey_year_params = self._read_config(survey_year_file_path)

        # assign parameters from configuration files and init params
        self.params = self._collect_parameters(init_params, survey_year_params)
        self.params['exclude_age1'] = exclude_age1

        # initialize all class variables
        self.strata_df = None
        self.geo_strata_df = None
        self.strata_sig_b = None
        self.length_df = None
        self.specimen_df = None
        self.nasc_df = None
        self.final_biomass_table = None

    @staticmethod
    def _check_init_file(init_file_path: str) -> None:
        """"""
        # TODO: create this function that checks the contents of the initialization config file
        # TODO: it should make sure that certain variables are defined too
        print("A check of the initialization file needs to be done!")

    @staticmethod
    def _check_survey_year_file(survey_year_file_path: str) -> None:
        # TODO: create this function that checks the contents of the survey year config file
        # TODO: it should make sure that certain variables are defined and all paths exist
        print("A check of the survey year file needs to be done!")

    @staticmethod
    def _read_config(file_path: str) -> dict:
        """
        Reads configuration files and returns a dictionary
        with the parameters specified in the file.

        Parameters
        ----------
        file_path: str
            Path to configuration file.
        """

        with open(file_path) as f:
            params = yaml.load(f, Loader=yaml.SafeLoader)

        return params

    @staticmethod
    def _set_params_from_init(source: int, init_params: dict) -> dict:
        """
        Constructs and assigns important variables using
        parameters from the initialization configuration file.

        Parameters
        ----------
        source : int
            The region of data to use.
            1 = US
            2 = Canada
            3 = US and Canada
        init_params : dict
            Parameters obtained from the initialization file

        Returns
        -------
        init_params : dict
            The input ``init_params`` with additional variables
        """

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

        init_params["source"] = source

        return init_params

    @staticmethod
    def _collect_parameters(init_params: dict, survey_params: dict) -> dict:
        """
        Collects all parameters defined in the initialization
        and survey year configuration files  into one variable.

        Parameters
        ----------
        init_params : dict
            Parameters obtained from the initialization file
        survey_params : dict
            Parameters obtained from the survey year file

        Returns
        -------
        full_params : dict
            All parameters obtained from both the survey year
            and initialization configuration files
        """

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

    def load_survey_data(self, file_type: str = 'all') -> None:
        """
        Loads the biological, NASC, and stratification
        data using parameters obtained from the configuration
        files.

        Parameters
        ----------
        file_type : str
            Specifies what survey data should be loaded.
            Possible options:
            - 'all' -> loads all survey data
            - 'biological' -> only loads the biological data
            - 'strata' -> only loads the stratification data
            - 'nasc' -> only loads the NASC data

        Notes
        -----
        This function assigns class variables obtained from loading the
        data. Specifically, the following class variables are created
        for the file_type:
        - ``file_type='biological'``
            - ``self.length_df``
            - ``self.specimen_df``
        - ``file_type='strata'``
            - ``self.strata_df``
            - ``self.geo_strata_df``
            - ``self.strata_sig_b``
        - `file_type='nasc'``
            - ``self.nasc_df``
        """

        if file_type not in ['all', 'biological', 'strata', 'nasc']:
            raise ValueError("file_type must be 'all', 'biological', 'strata', or 'nasc'!")

        # load specimen and length data
        if file_type in ('biological', 'all'):
            LoadBioData(self)

        # load all associated stratification data
        if file_type in ('strata', 'all'):
            LoadStrataData(self)

        if file_type in ('nasc', 'all'):
            self.nasc_df = load_nasc_data.load_nasc_df(self)

    def compute_biomass_density(self):
        """
        Computes the normalized biomass density and
        creates ``self.final_biomass_table``, which
        is a Pandas DataFrame that contains the
        normalized biomass density and associated
        useful variables.
        """

        bio_dense = ComputeBiomassDensity(self)
        bio_dense.get_final_biomass_table()

    def run_cv_analysis(self,
                        lat_inpfc: Tuple[float] = (np.NINF, 36, 40.5, 43.000, 45.7667, 48.5, 55.0000),
                        kriged_data=False, seed=None) -> float:
        """
        Performs CV analysis by running the Jolly-Hampton
        algorithm.

        Parameters
        ----------
        lat_inpfc : Tuple[float]
            Bin values which represent the latitude bounds for
            each region within a survey (established by INPFC)
        kriged_data : bool
            If True, perform CV analysis on Kriged data, otherwise
            perform CV analysis on data that has not been Kriged
        seed : int
            Seed value for the random number generator

        Returns
        -------
        The mean Jolly-Hampton CV value.

        Notes
        -----
        The format of ``lat_inpfc`` should be such that it can be
        used by Pandas.cut.

        If the initialization parameter ``JH_fac`` is 1, then only
        1 realization is run, otherwise 10,000 realizations of the
        algorithm are run.
        """

        if self.params["JH_fac"] == 1:
            nr = 1  # number of realizations
        else:
            nr = 10000  # number of realizations

        if kriged_data:
            raise NotImplementedError("CV analysis for kriged data has not been implemented")
        else:
            return cv_analysis.run_jolly_hampton(nr, lat_inpfc,
                                                 self.final_biomass_table,
                                                 self.params["JH_fac"], seed)

    def get_kriging_mesh(self) -> KrigingMesh:
        """
        Initializes a ``KrigingMesh`` object using
        parameters obtained from the configuration
        files.

        Returns
        -------
        KrigingMesh
            An initialized object that contains the mesh
            data, functions to transform the mesh, and
            functions to plot the mesh.

        Notes
        -----
        This function assigns class variables to the returned object.
        Specifically, the following class variables are created:
        - ``mesh_gdf`` a GeoPandas Dataframe representing the full mesh
        - ``smoothed_contour_gdf`` a GeoPandas Dataframe representing
        the smoothed contour (e.g. 200m isobath)
        """

        return KrigingMesh(self)

    def get_semi_variogram(self, x, y, field):

        semi_vario = SemiVariogram(x, y, field)

        return semi_vario

    def get_kriging(self):

        krig = Kriging(self)

        return krig
