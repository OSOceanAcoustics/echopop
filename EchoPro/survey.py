from pathlib import Path
from typing import List, Optional, Tuple, Union
from warnings import warn

import geopandas as gpd
import numpy as np
import xarray as xr
import yaml

from .computation import (
    Bootstrapping,
    ComputeTransectVariables,
    Kriging,
    SemiVariogram,
    generate_bin_ds,
    get_len_age_abundance,
    get_len_age_biomass,
    krig_param_type,
    krig_type_dict,
    run_jolly_hampton,
    vario_param_type,
    vario_type_dict,
)
from .data_loader import KrigingMesh, LoadBioData, LoadStrataData, load_nasc_df
from .utils.input_checks import check_existence_of_file


class Survey:
    """
    EchoPro base class that imports and prepares parameters for
    a survey. Additionally, it includes functions for accessing
    the modules associated with the transect and kriging variable
    calculations, CV analysis, semi-variogram algorithm, and Kriging.

    Parameters
    ----------
    init_file_path : str or pathlib.Path
        A string specifying the path to the initialization YAML file
    survey_year_file_path : str or pathlib.Path
        A string specifying the path to the survey year YAML file
    source : int
        The region of data to use.
        1 = US
        2 = Canada
        3 = US and Canada
    exclude_age1 : bool
        States whether age 1 hake should be included in analysis.
    """

    def __init__(
        self,
        init_file_path: Union[str, Path],
        survey_year_file_path: Union[str, Path],
        source: int = 3,
        exclude_age1: bool = True,
    ):

        # convert configuration paths to Path objects, if necessary
        init_file_path = Path(init_file_path)
        survey_year_file_path = Path(survey_year_file_path)

        self._check_init_file(init_file_path)
        self._check_survey_year_file(survey_year_file_path)

        # read initialization configuration file
        init_params = self._read_config(init_file_path)
        init_params = self._set_params_from_init(source, init_params)

        # read survey year configuration file
        survey_year_params = self._read_config(survey_year_file_path)

        # assign parameters from configuration files and init params
        self.params = self._collect_parameters(init_params, survey_year_params)
        self.params["exclude_age1"] = exclude_age1

        # convert all string paths to Path objects in params
        self._convert_str_to_path_obj()

        # initialize all class variables
        self.strata_df = None
        self.geo_strata_df = None
        self.strata_sig_b = None
        self.length_df = None
        self.specimen_df = None
        self.nasc_df = None
        self.bio_calc = None

    @staticmethod
    def _check_init_file(init_file_path: Path) -> None:
        """
        Ensures that the initialization configuration file
        exists and contains the appropriate contents.

        Parameters
        ----------
        init_file_path: Path
            The path to the initialization configuration file

        Raises
        ------
        FileNotFoundError
            If the file does not exist
        """

        # make sure the configuration file exists
        check_existence_of_file(init_file_path)

        # TODO: create this function that checks the contents of the initialization config file
        # TODO: it should make sure that certain variables are defined too
        print("A full check of the initialization file contents needs to be done!")

    @staticmethod
    def _check_survey_year_file(survey_year_file_path: Path) -> None:
        """
        Ensures that the survey year configuration file
        exists and contains the appropriate contents.

        Parameters
        ----------
        survey_year_file_path: Path
            The path to the survey year configuration file

        Raises
        ------
        FileNotFoundError
            If the file does not exist
        """

        # make sure the survey year file exists
        check_existence_of_file(survey_year_file_path)

        # TODO: create this function that checks the contents of the survey year config file
        # TODO: it should make sure that certain variables are defined and all paths exist
        print("A check of the survey year file contents needs to be done!")

    @staticmethod
    def _read_config(file_path: Path) -> dict:
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
        init_params["bio_hake_len_bin"] = np.linspace(
            init_params["bio_hake_len_bin"][0],
            init_params["bio_hake_len_bin"][1],
            num=init_params["bio_hake_len_bin"][2],
            dtype=np.int64,
        )

        # setting bio_hake_age_bin variable to a numpy array
        init_params["bio_hake_age_bin"] = np.linspace(
            init_params["bio_hake_age_bin"][0],
            init_params["bio_hake_age_bin"][1],
            num=init_params["bio_hake_age_bin"][2],
            dtype=np.int64,
        )

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
        param_intersect = set(init_params.keys()).intersection(
            set(survey_params.keys())
        )

        # if no parameters are the same, then run process, else return error
        if not param_intersect:

            # combine survey year and initialization parameters into one dictionary
            full_params = {}
            full_params.update(init_params)
            full_params.update(survey_params)

        else:
            raise RuntimeError(
                "The initialization and survey year configuration files define the same variable! "
                + f"\n These variables are: {param_intersect}"
            )

        return full_params

    def _convert_str_to_path_obj(self) -> None:
        """
        Converts all string paths to pathlib.Path objects in the
        class variable ``params``.

        Notes
        -----
        The class variable ``params`` will be directly modified.
        """

        # convert the root directory to a Path object
        self.params["data_root_dir"] = Path(self.params["data_root_dir"])

        for param_name, param_val in self.params.items():

            # convert each filename path to a Path object
            if "filename" in param_name:
                self.params[param_name] = Path(param_val)

    def load_survey_data(self, file_type: str = "all") -> None:
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

        if file_type not in ["all", "biological", "strata", "nasc"]:
            raise ValueError(
                "file_type must be 'all', 'biological', 'strata', or 'nasc'!"
            )

        # load specimen and length data
        if file_type in ("biological", "all"):
            LoadBioData(self)

        # load all associated stratification data
        if file_type in ("strata", "all"):
            LoadStrataData(self)

        if file_type in ("nasc", "all"):
            self.nasc_df = load_nasc_df(self)

    def compute_transect_results(
        self, selected_transects: Optional[List] = None
    ) -> None:
        """
        Constructs ``self.bio_calc.transect_results_gdf``, `
        `self.bio_calc.transect_results_male_gdf``, and
        ``self.bio_calc.transect_results_female_gdf``, which are
        GeoDataFrames that contain variables over the transect
        points (e.g. abundance, biomass).

        Parameters
        ----------
        selected_transects : list or None
            The subset of transects used in the calculations
        """

        self.bio_calc = None
        self.bio_calc = ComputeTransectVariables(self)
        self.bio_calc.get_transect_results_gdf(selected_transects)

        # create Dataset containing useful distributions and variables over length and age
        self.bio_calc.bin_ds = generate_bin_ds(self)

    def run_cv_analysis(
        self,
        lat_inpfc: Tuple[float] = (np.NINF, 36, 40.5, 43.000, 45.7667, 48.5, 55.0000),
        kriged_data=False,
        seed=None,
    ) -> float:
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
            if self.bio_calc.kriging_results_gdf is None:
                raise RuntimeError(
                    "Kriging must be ran before performing CV analysis on Kriged data!"
                )
        else:
            if self.bio_calc.transect_results_gdf is None:
                raise RuntimeError(
                    "The biomass density must be calculated before performing CV analysis on data!"
                )

        return run_jolly_hampton(self, nr, lat_inpfc, seed, kriged_data)

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

    def get_semi_variogram(
        self,
        krig_mesh: KrigingMesh = None,
        params: vario_param_type = {},
        warning: bool = True,
    ):
        """
        Initializes a ``SemiVariogram`` object based on the provided
        ``KrigingMesh`` object, the calculated areal biomass density,
        and user provided semi-variogram parameters.

        Parameters
        ----------
        krig_mesh : KrigingMesh
            Object representing the Kriging mesh
        params : dict
            Semi-variogram specific parameters. Contains the following
            parameters:
            - ``nlag: int`` -- The total number of lag centers
            - ``lag_res: float`` -- The spacing between lag centers
        warning : bool
            If True all warnings are printed to the terminal, otherwise
            they are silenced.

        Returns
        -------
        semi_vario : SemiVariogram
            An initialized object, which provides users with access
            to a routine that calculates the normalized
            semi-variogram and routines for obtaining the best
            semi-variogram model for the estimated semi-variogram.

        Warnings
        --------
        UserWarning
            If the final biomass table being used was created from a subset
            of the full data

        Raises
        ------
        ValueError
            If ``krig_mesh`` is not a ``KrigingMesh`` object
        ValueError
            If ``params`` is empty
        ValueError
            If ``params`` does not contain all required parameters
        TypeError
            If the values of ``params`` are not the expected type
        ValueError
            If the areal biomass density has not been calculated

        Notes
        -----
        To run this routine, one must first compute the areal biomass density
        using ``compute_transect_results``. It is standard to compute the biomass
        density from the full set of data (i.e. not from a subset of the data).
        """

        if not isinstance(krig_mesh, KrigingMesh):
            raise ValueError("You must provide a KrigingMesh object!")

        if not params:
            raise ValueError("You must provide parameters for the semi-variogram!")

        # make sure all parameters are included
        if set(vario_type_dict.keys()).difference(set(params.keys())):
            raise ValueError("Some required parameters where not provided!")

        # check that all types are correct for params
        for key, val in params.items():
            expected_type = vario_type_dict.get(key)
            if not isinstance(val, expected_type):
                raise TypeError(f"{key} is not of type {expected_type}")

        # provide a warning if the transect_results_gdf being used was
        # created from a subset of the full data
        if (len(self.bio_calc.transect_results_gdf) != len(self.nasc_df)) and warning:
            warn(
                "The biomass data being used is a subset of the full dataset. "
                "It is recommended that you use the biomass data created from the full dataset. "
                "To silence this warning set the warning argument to False."
            )

        if (not isinstance(self.bio_calc.transect_results_gdf, gpd.GeoDataFrame)) and (
            "biomass_density_adult" not in self.bio_calc.transect_results_gdf
        ):
            raise ValueError(
                "The areal biomass density must be calculated before running this routine!"
            )

        semi_vario = SemiVariogram(
            krig_mesh.transformed_transect_df.x_transect.values,
            krig_mesh.transformed_transect_df.y_transect.values,
            self.bio_calc.transect_results_gdf[
                "biomass_density_adult"
            ].values.flatten(),
            params["lag_res"],
            params["nlag"],
        )

        return semi_vario

    def get_kriging(self, params: krig_param_type) -> Kriging:
        """
        Initializes a ``Kriging`` object using the
        provided parameters

        Parameters
        ----------
        params : dict
            Kriging specific parameters. Contains the following parameters:
            - ``k_max: int`` -- the maximum number of data points within the
            search radius.
            - ``k_min: int`` -- the minimum number of data points within the
            search radius.
            - ``R: float`` -- search radius for Kriging
            - ``ratio: float`` -- acceptable ratio for the singular values
            divided by the largest singular value.
            - ``s_v_params: dict`` -- dictionary specifying the parameter values
            for the semi-variogram model.
            - ``s_v_model: Callable`` -- a Semi-variogram model from the ``SemiVariogram`` class

        Returns
        -------
        krig : Kriging
            An initialized ``Kriging`` object that provides users with
            access to routines that run Ordinary Kriging and routines
            that plot final Kriging results
        """

        if not params:
            raise ValueError("You must provide parameters for the Kriging routine!")

        # make sure all parameters are included
        if set(krig_type_dict.keys()).difference(set(params.keys())):
            raise ValueError("Some required parameters where not provided!")

        # check that all types are correct for params
        for key, val in params.items():
            expected_type = krig_type_dict.get(key)
            if not isinstance(val, expected_type):
                raise TypeError(
                    f"The Kriging parameter {key} is not of type {expected_type}"
                )

        krig = Kriging(
            self,
            params["k_max"],
            params["k_min"],
            params["R"],
            params["ratio"],
            params["s_v_params"],
            params["s_v_model"],
        )

        return krig

    def get_bootstrapping(self) -> Bootstrapping:
        """
        Initializes a ``Bootstrapping`` object.

        Returns
        -------
        boot: Bootstrapping
            An initialized ``Bootstrapping`` object that provides users with
            access to the routine ``run_bootstrapping``, which runs bootstrapping
            for data with Kriging and data without Kriging.
        """

        # initialize bootstrapping class
        boot = Bootstrapping(self)

        return boot

    def compute_length_age_variables(self, data: str = "transect") -> None:
        """
        Computes abundance and biomass over each length and age bin,
        for males, females, and all genders.

        Parameters
        ----------
        data : str
            Specifies the results produced:
            - 'all' -> Both Kriging and transect based variables
            - 'transect' -> only produces transect based variables
            - 'kriging' -> only produces Kriging variables

        Notes
        -----
        The computed DataFrames containing the specified data are assigned
        to class variables within ``self.biocalc``.
        Transect based results
            - ``self.bio_calc.transect_bin_abundance_male_df`` -> abundance at
             each length and age bin for males
            - ``self.bio_calc.transect_bin_abundance_female_df`` -> abundance at
             each length and age bin for females
            - ``self.bio_calc.transect_bin_abundance_df`` -> abundance at
             each length and age bin, when using all genders
            - A similar set of variables are created for biomass results with
             'abundance' replaced with 'biomass'. For example, biomass at each
             length and age bin when using all genders will be stored in the
             class variable ``self.bio_calc.transect_bin_biomass_df``.
        Kriging based results
            - An analogous set of variables are created for the Kriging based
             results with 'transect' replaced with 'kriging'. For example,
             biomass at each length and age bin when using all genders will
             be stored in ``self.bio_calc.kriging_bin_biomass_df``.
        """

        if not isinstance(self.bio_calc.bin_ds, xr.Dataset):
            raise RuntimeError(
                "self.bio_calc.bin_ds is not a Dataset, the routine "
                "self.compute_transect_results must be ran first."
            )

        if data in ["transect", "all"]:

            # ensure that the appropriate data exists
            if not isinstance(self.bio_calc.transect_results_gdf, gpd.GeoDataFrame):
                raise RuntimeError(
                    "self.bio_calc does not contain transect based results, "
                    "self.compute_transect_results must be ran first."
                )

            # obtain and assign abundance DataFrames for transect data
            (
                self.bio_calc.transect_bin_abundance_male_df,
                self.bio_calc.transect_bin_abundance_female_df,
                self.bio_calc.transect_bin_abundance_df,
            ) = get_len_age_abundance(
                gdf=self.bio_calc.transect_results_gdf,
                ds=self.bio_calc.bin_ds,
                kriging_vals=False,
            )

            # obtain and assign biomass DataFrames for transect data
            (
                self.bio_calc.transect_bin_biomass_male_df,
                self.bio_calc.transect_bin_biomass_female_df,
                self.bio_calc.transect_bin_biomass_df,
            ) = get_len_age_biomass(
                gdf_all=self.bio_calc.transect_results_gdf,
                gdf_male=self.bio_calc.transect_results_male_gdf,
                gdf_female=self.bio_calc.transect_results_female_gdf,
                ds=self.bio_calc.bin_ds,
                kriging_vals=False,
            )

        elif data in ["kriging", "all"]:

            # ensure that the appropriate data exists
            if not isinstance(self.bio_calc.kriging_results_gdf, gpd.GeoDataFrame):
                raise RuntimeError(
                    "self.bio_calc does not contain kriging based results, "
                    "The Kriging routine compute_kriging_variables must be "
                    "ran first."
                )

            # get_len_age_abundance(gdf=self.bio_calc.kriging_results_gdf,
            #                       ds=self.bio_calc.bin_ds, kriging=True)

            raise NotImplementedError(
                "Creating abundance and biomass over each length and "
                "age bin has not been implemented for Kriging data."
            )

        else:
            raise RuntimeError(
                "The input variable data must be 'all', 'transect', or 'kriging'!"
            )
