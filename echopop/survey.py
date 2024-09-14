import copy
from pathlib import Path
from typing import List, Literal, Optional, Union

from IPython.display import display

from .analysis import (
    acoustics_to_biology,
    apportion_kriged_values,
    krige,
    process_transect_data,
    stratified_summary,
    variogram_analysis,
)
from .core import DATA_STRUCTURE
from .graphics import variogram_interactive as egv
from .spatial.projection import transform_geometry
from .spatial.transect import edit_transect_columns
from .utils import load as el, load_nasc as eln, message as em
from .utils.validate import VariogramBase, VariogramInitial, VariogramOptimize


class Survey:
    """
    Echopop base class that imports and prepares parameters for
    a survey. Additionally, it includes functions for accessing
    the modules associated with the transect and Kriging variable
    calculations, CV analysis, semi-variogram algorithm, and Kriging.

    Parameters
    ----------
    init_config_path : str or pathlib.Path
        A string specifying the path to the initialization YAML file
    survey_year_config_path : str or pathlib.Path
        A string specifying the path to the survey year YAML file

    Attributes
    ----------
    meta : dict
        Metadata variable that provides summary information concerning the
        data contained within the class object.
    config : dict
        Configuration settings and parameters that can be referenced for
        various downstream and internal functions defined within the `init_config_path` and
        `survey_year_config_path` *.yaml files.
    input: dict
        Input data based on files included within the `survey_year_config_path` *.yaml file.
    analysis: dict
        Analysis variables and intermediate data products.
    results: dict
        Overall results produced by each of the analysis methods and workflows.
    """

    def __init__(
        self, init_config_path: Union[str, Path], survey_year_config_path: Union[str, Path]
    ):
        # Initialize `meta` attribute
        self.meta = copy.deepcopy(DATA_STRUCTURE["meta"])

        # Loading the configuration settings and definitions that are used to initialize the Survey
        # class object
        self.config = el.load_configuration(Path(init_config_path), Path(survey_year_config_path))

        # Initialize the `input` data attribute
        self.input = copy.deepcopy(DATA_STRUCTURE["input"])

        # Initialize the `analysis` data attribute
        self.analysis = copy.deepcopy(DATA_STRUCTURE["analysis"])

        # Initialize the `results` data attribute
        self.results = copy.deepcopy(DATA_STRUCTURE["results"])

    def load_acoustic_data(
        self,
        index_variable: Union[str, List[str]] = ["transect_num", "interval"],
        ingest_exports: Optional[Literal["echoview", "echopype"]] = None,
        region_class_column: str = "region_class",
        transect_pattern: str = r"T(\d+)",
        unique_region_id: str = "region_id",
        verbose: bool = True,
    ):
        """
        Loads in active acoustic backscatter survey data from .xlsx files, or processes
        Echoview `csv` file exports by consolidating them into aforementioned .xlsx files.

        Parameters
        ----------
        index_variable: Union[str, List[str]]
            Index columns used for defining discrete acoustic backscatter samples and vertical
            integration.
        ingest_exports: Literal['echoview', 'echopype']
            The type of acoustic backscatter exports required for generating the associated
            consolidated .xlsx files.
        region_class_column: str
            Dataframe column denoting the Echoview export region class (e.g. "zooplankton").
        transect_pattern: str
            A (raw) string that corresponds to the transect number embedded within the base name of
            the file path associated with each export file. Defaults to ``r'T(\\d+)'``. See a
            further description below for more details.
        unique_region_id: str
            Dataframe column that denotes region-specific names and identifiers.
        verbose: bool
            Console messages that will print various messages, updates, etc. when set to True.

        Notes
        ----------
        The string pattern for `transect_pattern` requires a consistent filename format that can be
        readily parsed by `read_echoview_exports`. The default value for `transect_pattern`
        (``r'T(\\d+)'``) enables the function to parse all numbers that trail the letter "T" in the
        base filename. For example, an example file of
        "C:/Path/User/Data/random_12_V34-T56-A78_G9.csv" would yield a transect number of '56'
        since it trails the "T" and does not accidentally use other numbers in the string. However,
        a filename like "C:/Path/User/Data/randomT_12_T34-T56-T78_G9.csv" would detect multiple
        transect numbers ('34', '56', and '78') since numbers trail the letter "T" in three places.
        Therefore, it is important to ensure that embedded transect numbers are differentiable from
        other digits that may appear in the base filename.
        """

        # Check `ingest_exports` argument
        if ingest_exports is not None and ingest_exports not in ["echoview", "echopype"]:
            raise ValueError("Argument `ingest_exports` must either be 'echoview' or 'echopype'.")

        # Compile echoview acoustic backscatter exports if `echoview_exports == True`:
        if ingest_exports is not None and ingest_exports == "echoview":
            eln.batch_read_echoview_exports(
                self.config,
                transect_pattern,
                index_variable,
                unique_region_id,
                region_class_column,
                verbose,
            )
            # ---- Update key for `export_regions`
            self.meta["provenance"]["imported_datasets"].update(["export_regions"])

        # Read in compiled `*.xlsx` acoustic backscatter data file(s) and additional validation
        # ---- Update key for `NASC`
        self.meta["provenance"]["imported_datasets"].update(["NASC"])
        # ---- Load the acoustic survey data
        el.load_dataset(self.input, self.config, dataset_type="NASC")

    def load_survey_data(self, verbose: bool = True):
        """
        Loads in biological and spatial survey data

        Parameters
        ----------
        verbose: bool
            Console messages that will print various messages, updates, etc. when set to True.
        """

        # Create haul-transect-mapping key file
        # el.write_haul_to_transect_key(self.config, verbose)

        # Get previously processed datasets
        # ---- Updated datasets
        new_datasets = ["biological", "kriging", "stratification"]
        # ---- Load in the new data
        el.load_dataset(self.input, self.config, dataset_type=new_datasets)
        # ---- Update key for the new datasets
        self.meta["provenance"]["imported_datasets"].update(new_datasets)

    def transect_analysis(
        self,
        species_id: Union[float, list[float]] = 22500,
        exclude_age1: bool = True,
        stratum: Literal["inpfc", "ks"] = "ks",
        verbose: bool = True,
    ):
        """
        Calculate population-level metrics from acoustic transect measurements
        """

        # Update settings to reflect the stratum definition
        self.analysis["settings"].update(
            {
                "transect": {
                    "age_group_columns": {
                        "haul_id": "haul_no_age1" if exclude_age1 else "haul_all_ages",
                        "nasc_id": "NASC_no_age1" if exclude_age1 else "NASC_all_ages",
                        "stratum_id": "stratum_no_age1" if exclude_age1 else "stratum_all_ages",
                    },
                    "species_id": species_id,
                    "stratum": stratum.lower(),
                    "stratum_name": "stratum_num" if stratum == "ks" else "inpfc",
                    "exclude_age1": exclude_age1,
                }
            }
        )

        # Initial data processing of the transect biological and acoustic data
        self.analysis["transect"] = process_transect_data(
            self.input, self.analysis["transect"], self.analysis["settings"], self.config
        )

        # Convert NASC into number density (animals/nmi^2), biomass density (kg/nmi^2), abundance
        # (# animals), and biomass (kg) for all fish, sexed (male/female) fish, and unsexed fish
        # ---- This further provides the resulting distributions of biomass and abundance over
        # ---- length and age for each sex across the entire survey
        biomass_summary, self.analysis["transect"] = acoustics_to_biology(
            self.input, self.analysis["transect"], self.config, self.analysis["settings"]
        )

        # Add the biomass summary table to the results attribute
        # ---- Update results (biomass summary)
        self.results["transect"].update({"biomass_summary_df": biomass_summary})

        # Update meta provenance tracking
        self.meta["provenance"].update(
            {"transect": self.analysis["settings"]["transect"]["stratum"]}
        )

        # Print result if `verbose == True`
        if verbose:
            em.transect_results_msg(self.results["transect"], self.analysis["settings"]["transect"])

    def stratified_analysis(
        self,
        dataset: Literal["transect", "kriging"] = "transect",
        stratum: Literal["inpfc", "ks"] = "inpfc",
        variable: Literal["abundance", "biomass", "nasc"] = "biomass",
        mesh_transects_per_latitude: Optional[int] = None,
        transect_sample: Optional[float] = None,
        transect_replicates: Optional[float] = None,
        bootstrap_ci: float = 0.95,
        bootstrap_ci_method: Literal[
            "BC", "BCa", "empirical", "percentile", "standard", "t-jackknife", "t-standard"
        ] = "BCa",
        bootstrap_ci_method_alt: Optional[
            Literal["empirical", "percentile", "standard", "t-jackknife", "t-standard"]
        ] = "t-jackknife",
        bootstrap_adjust_bias: bool = True,
        verbose=True,
    ):
        """
        Calculates the stratified summary statistics for biomass

        Notes
        -----
        This function calculates estimates and confidence intervals (95%) for biomass mean,
        variance, and coefficients of variation (CVs). This currently only calculates this
        metric for adult animals (age-2+) and is not calculated for other contrasts such as
        age-class and sex. This also only applies to the transect results and is not currently
        metric for adult animals (age-2+) and is not calculated for other contrasts such as
        age-class and sex. This also only applies to the transect results and is not currently
        designed to be compatible with other derived population-level statistics (e.g. kriging).
        """

        # Error message for `stratum == 'ks'`
        if stratum == "ks":
            raise ValueError(
                """The Jolly and Hampton (1990) stratified analysis is not"""
                """ currently compatible for calculating over KS strata. Please change `stratum` """
                """ to 'inpfc'."""
            )

        # Parameterize analysis settings that will be applied to the stratified analysis
        self.analysis["settings"].update(
            {
                "stratified": {
                    "dataset": f"{dataset}",
                    "stratum": stratum.lower(),
                    "stratum_name": "stratum_num" if stratum == "ks" else "stratum_inpfc",
                    "transect_sample": (
                        self.config["stratified_survey_mean_parameters"][
                            "strata_transect_proportion"
                        ]
                        if transect_sample is None
                        else transect_sample
                    ),
                    "transect_replicates": (
                        self.config["stratified_survey_mean_parameters"]["num_replicates"]
                        if transect_replicates is None
                        else transect_replicates
                    ),
                    "variable": variable,
                    "exclude_age1": self.analysis["settings"]["transect"]["exclude_age1"],
                    "verbose": verbose,
                    "bootstrap_ci_method": bootstrap_ci_method,
                    "bootstrap_ci_method_alt": bootstrap_ci_method_alt,
                    "bootstrap_ci": bootstrap_ci,
                    "bootstrap_adjust_bias": bootstrap_adjust_bias,
                }
            }
        )
        # ---- Append kriging-specific parameters if necessary
        if dataset == "kriging":
            self.analysis["settings"]["stratified"].update(
                {
                    "mesh_transects_per_latitude": (
                        self.config["stratified_survey_mean_parameters"][
                            "mesh_transects_per_latitude"
                        ]
                        if mesh_transects_per_latitude is None
                        else mesh_transects_per_latitude
                    ),
                    "variable": (
                        self.analysis["settings"]["kriging"]["variable"].replace("_density", "")
                    ),
                }
            )

        # Calculate the stratified mean, variance, and coefficient of variation
        stratified_results, self.analysis = stratified_summary(
            self.analysis,
            self.results,
            self.input["spatial"],
            self.analysis["settings"]["stratified"],
        )

        # Add the stratified statistics dictionary results to the `results` attribute
        # ---- Update results (stratified results)
        self.results["stratified"].update({f"{dataset}": stratified_results})

        # Print result if `verbose == True`
        if verbose:
            em.stratified_results_msg(stratified_results, self.analysis["settings"]["stratified"])

    def variogram_gui(self):
        """
        Semivariogram plotting and parameter optimization GUI method
        """

        # Initialize results
        self.results["variogram"] = {}

        # Initialize Survey-class object
        self.analysis.update({"variogram": {}})

        # Get the stratum name
        stratum_name = self.analysis["settings"]["transect"]["stratum_name"]

        # Get standardization config for kriging
        standardization_parameters = self.input["statistics"]["kriging"]["model_config"]
        # ---- Get isobath data
        isobath_df = self.input["statistics"]["kriging"]["isobath_200m_df"]

        # Get variogram parameters
        variogram_parameters = self.input["statistics"]["variogram"]["model_config"].copy()

        # Generate settings dictionary
        settings_dict = {
            "stratum_name": stratum_name,
            "variable": "biomass",
            "verbose": False,
            "kriging_parameters": {
                "longitude_reference": standardization_parameters["longitude_reference"],
                "longitude_offset": standardization_parameters["longitude_offset"],
                "latitude_offset": standardization_parameters["latitude_offset"],
            },
        }

        # Prepare the transect data
        # ---- Create a copy of the transect dictionary
        transect_input = copy.deepcopy(self.analysis["transect"])
        # ---- Edit the transect data
        transect_data = edit_transect_columns(transect_input, settings_dict)
        isobath_df = self.input["statistics"]["kriging"]["isobath_200m_df"]
        transect_data, _, _ = transform_geometry(transect_data, isobath_df, settings_dict)

        # Generate GUI
        SEMIVARIOGRAM_GUI = egv.variogram_widgets(
            transect_data,
            variogram_parameters,
            settings_dict,
            self.analysis["variogram"],
            self.results["variogram"],
        )

        # Run GUI
        display(SEMIVARIOGRAM_GUI)

        # Update the results
        # self.analysis["variogram"].update({
        #     "model_fit": SEMIVARIOGRAM_GUI.results["best_fit"]["model_fit"],
        #     "model": SEMIVARIOGRAM_GUI.results["variogram"]["model"]
        # })

    def fit_variogram(
        self,
        variogram_parameters: VariogramBase = {},
        optimization_parameters: VariogramOptimize = {},
        model: Union[str, List[str]] = ["bessel", "exponential"],
        n_lags: int = 30,
        azimuth_range: float = 360.0,
        standardize_coordinates: bool = True,
        force_lag_zero: bool = True,
        initialize_variogram: VariogramInitial = [
            "nugget",
            "sill",
            "correlation_range",
            "hole_effect_range",
            "decay_power",
        ],
        variable: Literal["biomass", "abundance"] = "biomass",
        verbose: bool = True,
    ):
        """
        Compute the best-fit variogram parameters for tansect data

        Parameters
        ----------
        variogram_parameters: VariogramBase
            A dictionary comprising various arguments required for computing the model variogram.
            See :fun:`echopop.utils.validate.VariogramBase` and
            :fun:`echopop.spatial.variogram.variogram` for more details on the required/default
            parameters.
        optimization_parameters: VariogramOptimize
            A dictionary comprising various arguments for optimizing the variogram fit via
            non-linear least squares. See :fun:`echopop.utils.validate.VariogramOptimize` for more
            details on the required/default parameters.
        initialize_variogram: VariogramInitial
            A dictionary or list that indicates how each variogram parameter (see
            :fun:`echopop.spatial.variogram.variogram` for more details) is configured for
            optimization. Including parameter names in a list will incorporate default initial
            values imported from the associated file in the configuration *.yaml are used instead.
            This also occurs when `initialize_variogram` is formatted as a dictionary and the
            'value' key is not present for defined parameters. Parameter names excluded from either
            the list or dictionary keys are assumed to be held as fixed values. See
            :fun:`echopop.utils.validate.VariogramInitial` and
            :fun:`echopop.utils.validate.InitialValues` for more details.
        model: Union[str, List[str]]
            A string or list of model names. A single name represents a single family model. Two
            inputs represent the desired composite model (e.g. the composite J-Bessel and
            exponential model). Defaults to: ['bessel', 'exponential']. Available models and their
            required arguments can be reviewed in the :fun:`echopop.spatial.variogram.variogram`
            function.
        azimuth_range: float
            The total azimuth angle range that is allowed for constraining
            the relative angles between spatial points, particularly for cases where a high degree
            of directionality is assumed.
        n_lags: int
            See the `variogram_parameters` argument in
            :fun:`echopop.spatial.variogram.empirical_variogram` for more details on
            `n_lags`.
        force_lag_zero: bool
            See the `variogram_parameters` argument in
            :fun:`echopop.spatial.variogram.empirical_variogram` for more details on
            `force_lag_zero`.
        standardize_coordinates: bool
            When set to `True`, transect coordinates are standardized using reference coordinates.
        variable: Literal["biomass", "abundance"]
            Transect data values used for fitting the variogram. This includes two options:
            "abundance" and "biomass", with the default being "biomass". These inputs correspond
            to fitting the empirical and theoretical variograms on "number density" and "biomass
            density", respectively.
        verbose: bool
            When set to `True`, optional console messages and reports are provided to users.

        Notes
        -----
        The variogram model fitting methods makes use of the `lmfit` library. Values included in
        the `variogram_parameters` argument, but omitted from `initialize_variogram`, use default
        values imported from `self.input["statistics"]["variogram"]["model_config"]`.
        """

        # Validate "variable" input
        if variable not in ["biomass", "abundance"]:
            raise ValueError(
                f"The user input for `variable` ({variable}) is invalid. Only `variable='biomass'` "
                f"and `variable='abundance'` are valid inputs for the `fit_variogram()` method."
            )

        # Initialize Survey-class object
        self.analysis.update({"variogram": {}})

        # Parameterize analysis settings that will be applied to the variogram fitting and analysis
        self.analysis["settings"].update(
            {
                "variogram": {
                    "azimuth_range": azimuth_range,
                    "fit_parameters": (
                        initialize_variogram.keys()
                        if isinstance(initialize_variogram, dict)
                        else initialize_variogram
                    ),
                    "force_lag_zero": force_lag_zero,
                    "model": model,
                    "standardize_coordinates": standardize_coordinates,
                    "stratum_name": self.analysis["settings"]["transect"]["stratum_name"],
                    "variable": variable,
                    "verbose": verbose,
                }
            }
        )

        # Append `kriging_parameters` to the settings dictionary
        if standardize_coordinates:
            self.analysis["settings"]["variogram"].update(
                {"kriging_parameters": self.input["statistics"]["kriging"]["model_config"]}
            )

        # Create a copy of the existing variogram settings
        default_variogram_parameters = self.input["statistics"]["variogram"]["model_config"].copy()
        # ---- Update model, n_lags
        default_variogram_parameters.update({"model": model, "n_lags": n_lags})

        # Create optimization settings dictionary
        # ---- Add to settings
        self.analysis["settings"]["variogram"].update({"optimization": optimization_parameters})

        # Find the best-fit variogram parameters
        best_fit_variogram = variogram_analysis(
            variogram_parameters,
            default_variogram_parameters,
            optimization_parameters,
            initialize_variogram,
            self.analysis["transect"],
            self.analysis["settings"]["variogram"],
            self.input["statistics"]["kriging"]["isobath_200m_df"],
        )

        # Add "partial" results to analysis attribute
        self.analysis["variogram"].update(
            {
                "model": model,
                "initial_fit": best_fit_variogram["initial_fit"],
                "optimized_fit": best_fit_variogram["optimized_fit"],
            }
        )

        # Add variogram result
        self.results.update(
            {"variogram": {"model_fit": best_fit_variogram["best_fit_parameters"], "model": model}}
        )

        # Print result if `verbose == True`
        if verbose:
            em.variogram_results_msg(self.analysis["variogram"])

    # !!! TODO: develop different name for "crop_method = 'interpolation'"
    def kriging_analysis(
        self,
        bearing_tolerance: float = 15.0,
        coordinate_transform: bool = True,
        crop_method: Literal["transect_ends", "convex_hull"] = "transect_ends",
        extrapolate: bool = False,
        best_fit_variogram: bool = True,
        kriging_parameters: Optional[dict] = None,
        latitude_resolution: float = 1.25,
        mesh_buffer_distance: float = 1.25,
        num_nearest_transects: int = 4,
        projection: Optional[str] = None,
        stratum: str = "ks",
        variable: str = "biomass_density",
        variogram_model: Union[str, List[str]] = ["bessel", "exponential"],
        variogram_parameters: Optional[dict] = None,
        verbose: bool = True,
    ):
        """
        Interpolates biomass data using ordinary kriging


        Parameters
        ----------
        variable
            Biological variable that will be interpolated via kriging
        """

        # Parameterize analysis settings that will be applied to the kriging analysis
        self.analysis["settings"].update(
            {
                "kriging": {
                    "exclude_age1": self.analysis["settings"]["transect"]["exclude_age1"],
                    "extrapolate": extrapolate,
                    "kriging_parameters": (
                        self.input["statistics"]["kriging"]["model_config"]
                        if kriging_parameters is None
                        else kriging_parameters
                    ),
                    "projection": (
                        self.config["geospatial"]["init"] if projection is None else projection
                    ),
                    "standardize_coordinates": coordinate_transform,
                    "stratum": stratum.lower(),
                    "stratum_name": "stratum_num" if stratum == "ks" else "inpfc",
                    "variable": variable,
                    "variogram_parameters": (
                        self.input["statistics"]["variogram"]["model_config"]
                        if variogram_parameters is None
                        else variogram_parameters
                    ),
                    "verbose": verbose,
                }
            }
        )

        # Update variogram model
        self.analysis["settings"]["kriging"]["variogram_parameters"]["model"] = variogram_model

        # Update variogram parameters to use fitted if the values are available
        if best_fit_variogram:
            if "variogram" in self.results:
                if "model_fit" in self.results["variogram"]:
                    # ---- Parameters
                    self.analysis["settings"]["kriging"].update(
                        {"variogram_parameters": self.results["variogram"]["model_fit"]}
                    )
                    # ---- Update model
                    self.analysis["settings"]["kriging"]["variogram_parameters"]["model"] = (
                        self.results["variogram"]["model"]
                    )
                else:
                    raise ValueError(
                        "Argument `best_fit_variogram` is invalid. No fitted variogram parameters "
                        "have been estimated via the `fit_variogram()` Survey-class method. "
                        "Either run the `fit_variogram()` method first, or set the argument "
                        "`best_fit_variogram=False`."
                    )

        # Prepare temporary message concerning coordinate transformation = False
        if not coordinate_transform:
            raise ValueError(
                """Kriging without coordinate standardization is currently """
                """unavailable due to the kriging parameter `search_radius` being only defined """
                """for transformed x- and y-coordinates."""
            )

        # Run kriging analysis
        # ----> Generates a georeferenced dataframe for the entire mesh grid, summary statistics,
        # ----> and adds intermediate data products to the analysis attribute
        # ---- If kriging results are not extrapolated beyond the survey region:
        if not extrapolate:
            # ---- Update the analysis settings
            self.analysis["settings"]["kriging"].update(
                {
                    "bearing_tolerance": bearing_tolerance,
                    "crop_method": crop_method,
                    "latitude_resolution": latitude_resolution,
                    "mesh_buffer_distance": mesh_buffer_distance,
                    "num_nearest_transect": num_nearest_transects,
                }
            )
            # ---- Run kriging algorithm
            kriged_results, self.analysis = krige(
                self.input, self.analysis, self.analysis["settings"]["kriging"]
            )
        # ---- If kriging results are extrapolated beyond the survey region:
        else:
            # ---- Run kriging algorithm
            kriged_results, self.analysis = krige(
                self.input, self.analysis, self.analysis["settings"]["kriging"]
            )

        # Save the results to the `results` attribute
        self.results.update({"kriging": kriged_results})

        # Distribute the kriged results over length and age bins
        aged_apportioned, unaged_apportioned, kriged_apportioned_table = apportion_kriged_values(
            self.analysis, kriged_results["mesh_results_df"], self.analysis["settings"]["kriging"]
        )

        # Modify the resulting tables if age-1 fish are excluded!
        # if settings_dict[ 'exclude_age1' ]:

        # ---- Update results
        self.results["kriging"].update(
            {
                "tables": {
                    "overall_apportionment_df": kriged_apportioned_table,
                    "aged_tbl": aged_apportioned,
                    "unaged_tbl": unaged_apportioned,
                }
            }
        )

        # Print result if `verbose == True`
        if verbose:
            em.kriging_results_msg(self.results["kriging"], self.analysis["settings"]["kriging"])

    def summary(self, results_name: str):
        """
        Summary property that prints out formatted results for the desired analysis

        Parameters
        ----------
        results_name: str
            The name of the results that should be printed into the console. This can either be
            formatted as a single input name (e.g. 'transect' , 'kriging') or a nested/layered
            variable (e.g. 'stratified:transect') where a colon (':') is used as the delimiter that
            separates the two result layer names.
        """
        # Break up `results_name` if it contains a ':' delimiter
        if ":" in results_name:
            # ---- Split the string
            input_1, input_2 = results_name.split(":")
            # ---- Generate the string
            eval_string = f"""em.{input_1}_results_msg(self.results['{input_1}']['{input_2}'], \
            self.analysis['settings']['{input_1}'])"""
        else:
            # ---- Generate the string
            eval_string = f"""em.{results_name}_results_msg(self.results['{results_name}'], \
            self.analysis['settings']['{results_name}'])"""

        # Print out the result
        return eval(eval_string)
