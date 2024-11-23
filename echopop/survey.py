import copy
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Union

import numpy as np
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
from .graphics import plotting as egp, variogram_interactive as egv
from .spatial.projection import transform_geometry
from .spatial.transect import edit_transect_columns
from .utils import load as el, load_nasc as eln, message as em
from .utils.load import dataset_integrity


class Survey:
    """
    Echopop survey analysis class

    This class includes methods for ingesting, processing, and visualizing spatially distributed
    acoustic backscatter and biological data. These are used to provide various population
    and survey uncertainty estimates

    Parameters
    ----------
    init_config_path : str or pathlib.Path
        A string specifying the path to the initialization ``*.yaml`` file

    survey_year_config_path : str or pathlib.Path
        A string specifying the path to the survey year ``*.yaml`` file

    Attributes
    ----------
    meta: Dict[str, Any]
        Metadata variable that provides summary information concerning the
        data contained within the class object

    config : Dict[str, Any]
        Configuration settings and parameters that can be referenced for
        various downstream and internal functions defined within the `init_config_path` and
        `survey_year_config_path` ``*.yaml`` files

    input: Dict[str, Any]
        Input data based on files included within the `survey_year_config_path` ``*.yaml`` file

    analysis: Dict[str, Any]
        Analysis variables and intermediate data products

    results: Dict[str, Any]
        Overall results produced by each of the analysis methods and workflows

    Notes
    -----

    See `Survey-class data structure \
        <https://echopop.readthedocs.io/en/latest/core_data_structure.html>`_
    for more details on how this object is organized and used.

    See `Configuration file formatting \
        <https://echopop.readthedocs.io/en/latest/implementation/preprocessing_data.html>`_
    for more details on how how the user-defined files for ``init_config_path`` and
    ``survey_year_config_path`` should be organized for successful integration into ``Echopop``


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
        Echoview ``*.csv`` file exports by consolidating them into aforementioned ``*.xlsx`` files

        Parameters
        ----------
        index_variable: Union[str, List[str]]
            Index columns used for defining discrete acoustic backscatter samples and vertical
            integration

        ingest_exports: Optional[Literal["echoview", "echopype"]]
            The type of acoustic backscatter exports required for generating the associated
            consolidated ``*.xlsx`` files

        region_class_column: str
            Dataframe column denoting the Echoview export region class (e.g. "zooplankton")

        transect_pattern: str
            A (raw) string that corresponds to the transect number embedded within the base name of
            the file path associated with each export file. The default (raw) string and regular
            expression ``r"T(\\d+)"`` matches any cases in the filename where there are numbers
            appended to a leading "T". For instance, the filepath
            ``"C:/Path/To/Folder/T01random_text.csv"`` would correspond to a transect number of
            ``56``. However, this sort of pattern would not perform well where this pattern is
            repeated multiple times throughout the filepath and filename strings. For example,
            ``"C:/Path/To/Folder/T01/survey_T999_text.csv"`` would correspond to multiple transect
            numbers: ``[1, 999]``. It is crucial to ensure that ``transect_pattern`` can
            be readily extracted from the full filepath, so it is therefore important to use a
            istinctive expression

        unique_region_id: str
            Dataframe column that denotes region-specific names and identifiers

        verbose: bool
            Console messages that will print various messages, updates, etc. when ``verbose=True``

        """

        # Check `ingest_exports` argument
        if ingest_exports is not None and ingest_exports not in ["echoview", "echopype"]:
            raise ValueError("Argument `ingest_exports` must either be 'echoview' or 'echopype'.")

        # Compile echoview acoustic backscatter exports if `echoview_exports == True`:
        if ingest_exports is not None and ingest_exports == "echoview":
            eln.ingest_echoview_exports(
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
            Console messages that will print various messages, updates, etc. when ``verbose=True``

        """

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

        Parameters
        ----------

        species_id: Union[float, list[float]]
            A number or string code for target species that are incorporated into the ingested
            datasets

        exclude_age1: bool
            When ``exclude_age1=True``, age-1 fish are excluded from the overall analysis. This
            means that age, length, and joint age-length count/weight distributions are computed
            for age-2+ fish only. However, it is assumed that some age-1 fish "leak" into these
            calculations, so they are not excluded entirely from the survey population estimates.
            The default ``species_id=22500`` corresponds to the Pacific hake (*Mercluccius
            productus*) code used by NWFSC-FEAT

        stratum: Literal["inpfc", "ks"]
            Define which stratification to use. Options include:

            - *"inpfc"* \n
            Latitude-based (INPFC)

            - *"ks"* \n
            Clusters based on fish fork length via Kolmogorov-Smirnov distribution tests

        verbose: bool
            When ``verbose=True``, optional console messages and reports are provided to users

        """

        # Check dataset integrity
        dataset_integrity(self.input, analysis="transect")

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
                    "stratum_name": "stratum_num" if stratum == "ks" else "stratum_inpfc",
                    "unique_strata": (
                        np.unique(self.input["spatial"]["strata_df"]["stratum_num"])
                        if stratum == "ks"
                        else np.unique(self.input["spatial"]["inpfc_strata_df"]["stratum_inpfc"])
                    ),
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
        transect_replicates: Optional[int] = None,
        bootstrap_ci: float = 0.95,
        bootstrap_ci_method: Literal[
            "BC", "BCa", "empirical", "percentile", "standard", "t-jackknife", "t-standard"
        ] = "t-jackknife",
        bootstrap_ci_method_alt: Optional[
            Literal["empirical", "percentile", "standard", "t-jackknife", "t-standard"]
        ] = "t-standard",
        bootstrap_adjust_bias: bool = True,
        verbose=True,
    ):
        """
        Calculates the stratified summary statistics for biomass

        Parameters
        ----------

        dataset: Literal["kriging", "transect"]
            Define whether population estimates from the kriged mesh (``dataset="kriging"``) or
            along-transect measurements (``dataset="transect"``)'

        ci_percentile: float
            The confidence interval percentile used for computing the bootstrap
            confidence/uncertainty intervals. The default value is ``ci_percentile=0.95``, which
            corresponds to the 95% confidence interval (*CI*)

        bootstrap_ci_method: Literal["BC", "BCa", "empirical", "percentile", "standard", \
            "t-jackknife", "t-standard"]'
            The method for computing the bootstrap *CI*. Each method may produce different
            estimates for the lower and upper bounds. The default algorithm  is
            ``bootstrap_ci_method="t-jackknife"``, which is a jackknife method for computing the
            *CI* assuming a *t*-distribution. However, the choice of ``bootstrap_ci_method``
            is largely a decision between preferring computation speed versus bias adjustments.
            For instance, the ``"percentile"`` and ``"standard"`` approaches will produce wider
            estimates but do not incorporate bias adjustments like the bias-controlled methods
            (i.e. ``"BC"`` and ``"BCa"``). See
            :func:`echopop.statistics.bootstrap_confidence_intervals` for more details and
            references for these calculations


        stratum: Literal["inpfc", "ks"]
            Define which stratification to use. Options include:

            - *"inpfc"* \n
            Latitude-based (INPFC)

            - *"ks"* \n
            Clusters based on fish fork length via Kolmogorov-Smirnov distribution tests

        variable: Literal["abundance", "biomass", "nasc"]
            Population estimate used for the stratification analysis including abundance
            (``variable="abundance"``), biomass (``variable="biomass"``), and acoustic
            NASC (``variable="nasc"``)

        verbose: bool
            When ``verbose=True``, optional console messages and reports are provided to users

        Other Parameters
        ----------------

        bootstrap_adjust_bias: bool
            When ``bootstrap_adjust_bias=True``, the bootstrap distribution bias is calculated and
            subtracted from the bootstrapped samples prior to the *CI* being computed

        bootstrap_ci_method_alt: Optional[Literal["empirical", "percentile", "standard", \
            "t-jackknife", "t-standard"]]
            An optional argument that provides an alternative *CI* calculation for cases where
            certain algorithms (e.g. "BCa") are incompatible with the bootstrapped estimate
            distributions (e.g. too skewed). This defaults to
            ``bootstrap_ci_method_alt="t-standard"``

        mesh_transects_per_latitude: Optional[int]
            The number of "virtual" transects per degree latitude that are created from the
            kriged mesh estimates when ``dataset="kriging"``. This argument
            can also be configured within the initialization ``*.yml`` file where:

            .. code-block:: YAML

                stratified_survey_mean_parameters:
                    mesh_transects_per_latitude: 5


        transect_sample: Optional[float]
            The proportion of transects that are resampled (without replacement)
            for the Jolly and Hampton (1990) algorithm. This argument can also be
            configured within the initialization ``*.yml`` file where:

            .. code-block:: YAML

                stratified_survey_mean_parameters:
                    strata_transect_proportion: 0.75

        transect_replicates: Optional[int]
            The number of iterations used for the Jolly and Hampton (1990) algorithm.
            This can be configured within the initialization ``*.yml`` file where:

            .. code-block:: YAML

                stratified_survey_mean_parameters:
                    num_replicates: 10000

        Notes
        -----
        This method calculates estimates and confidence intervals (95%) for biomass for each
        stratum and the entire survey using the random stratified sampling algorithm developed by
        Jolly and Hampton (1990) [1]_. This provides an approach for calculating a survey
        coefficient of variation (*CV*) that represents the overall sampling uncertainty in the
        survey. This currently only calculates this metric for adult animals (age-2+) and is not
        calculated for other contrasts such as age-class and sex. This also only applies to the
        transect results and is not currently metric for adult animals (age-2+) and is not
        calculated for other contrasts such as age-class and sex. This also only applies to the
        transect results and is not currently designed to be compatible with other derived
        population-level statistics (e.g. kriging).

        See Also
        --------
        :func:`echopop.statistics.bootstrap_confidence_intervals`
            Wrapper function for computing the bootstrap *CI* using the user-defined method

        References
        ----------

        .. [1] G. M. Jolly and I. Hampton. (1990). *A stratified random transect design for acoustic
            surveys of fish stocks*. Canadian Journal of Fisheries and Aquatic Sciences, 47(7):
            1282-1291, doi:10.1139/f90-147.

        """

        # Check dataset integrity
        dataset_integrity(self.input, analysis=f"stratified:{dataset}")

        # Error message for `stratum == 'ks'`
        if stratum == "ks":
            raise ValueError(
                "The Jolly and Hampton (1990) stratified analysis is not currently compatible for "
                "calculating over KS strata. Please change `stratum` to `'inpfc'`."
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

        *See :meth:`.fit_variogram` for more details on how the ``Survey.variogram_gui()`` method
        performs the variogram model fitting*

        See Also
        --------
        :meth:`.fit_variogram`
            Variogram fitting method

        """

        # Check dataset integrity
        dataset_integrity(self.input, analysis="variogram")

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

    def fit_variogram(
        self,
        variogram_parameters: Dict[str, Any] = {},
        optimization_parameters: Dict[str, Any] = {},
        model: Union[str, List[str]] = ["bessel", "exponential"],
        n_lags: int = 30,
        azimuth_range: float = 360.0,
        standardize_coordinates: bool = True,
        force_lag_zero: bool = True,
        initialize_variogram: Union[List[str], Dict[str, Any]] = [
            "nugget",
            "sill",
            "correlation_range",
            "hole_effect_range",
            "decay_power",
        ],
        variable: Literal["biomass"] = "biomass",
        verbose: bool = True,
    ):
        """
        Compute the best-fit variogram parameters for tansect data

        Parameters
        ----------
        variogram_parameters: VariogramBase
            A dictionary containing user-defined variogram parameters. See
            :class:`echopop.spatial.variogram.VariogramBase` for more details. Valid arguments
            include:

            - **sill: realposfloat** \n
            The asymptotic value as lags approach infinity

            - **nugget: realposfloat** \n
            The semivariogram *y*-intercept that corresponds to variability at lag distances
            shorter than the lag resolution

            - **correlation_range: realposfloat** \n
            The relative length scale, or range, at which the autocorrelation between lag distances
            no longer increases and becomes asymptotic

            - **hole_effect_range: realposfloat** \n
            The (normalized) length scale/range that holes' are observed, which represent 'null'
            (or very small) points compared to their neighboring lags

            - **decay_power: realposfloat** \n
            An exponential term that is used in certain generalized exponential (or related)
            semivariogram models that modulates the ascending rate for a semivariogram

            - **enhance_semivariance: bool** \n
            A boolean term that determines whether the correlation decay in certain  cosine-related
            variogram models are enhanced (or not) with increasing lag distances

        optimization_parameters: VariogramOptimize
            A dictionary comprising various arguments for optimizing the variogram fit via
            non-linear least squares. See :class:`echopop.utils.validate.VariogramOptimize` for more
            details on the required/default parameters.

        initialize_variogram: VariogramInitial
            A dictionary or list that indicates how each variogram parameter (see
            :func:`echopop.spatial.variogram.variogram` for more details) is configured for
            optimization. Including parameter names in a list will incorporate default initial
            values imported from the associated file in the configuration ``*.yaml`` are used
            instead. This also occurs when ``initialize_variogram`` is formatted as a dictionary
            and the `'value'` key is not present for defined parameters. Parameter names excluded
            from either the list or dictionary keys are assumed to be held as fixed values. See
            :class:`echopop.utils.validate.VariogramInitial` for more details

        model: Union[str, List[str]]
            A string or list of model names. A single name represents a single family model. Two
            inputs represent the desired composite model (e.g. the composite J-Bessel and
            exponential model). Defaults to: ``model=["bessel", "exponential"]``. Available
            models and their required arguments can be reviewed in the
            :func:`echopop.spatial.variogram.variogram` function

        azimuth_range: realcircle
            The total azimuth angle range that is allowed for constraining  the relative angles
            between spatial points, particularly for cases where a high degree of directionality is
            assumed

        n_lags: int
            See the ``variogram_parameters`` argument in
            :func:`echopop.spatial.variogram.empirical_variogram` for more details on
            ``n_lags``

        force_lag_zero: bool
            See the ``variogram_parameters`` argument in
            :func:`echopop.spatial.variogram.empirical_variogram` for more details on
            ``force_lag_zero``

        standardize_coordinates: bool
            When set to ``True``, transect coordinates are standardized using reference coordinates

        variable: Literal["biomass", "abundance"]
            Transect data values used for fitting the variogram. This includes two options:
            "abundance" and "biomass", with the default being "biomass". These inputs correspond
            to fitting the empirical and theoretical variograms on "number density" and "biomass
            density", respectively. Note that ``Echopop`` currently only accepts
            ``variable="biomass"``. Support for additional variables will be available in future
            releases

        verbose: bool
            When set to ``True``, optional console messages and reports are provided to users.

        Notes
        -----
        The variogram model fitting methods makes use of the ``lmfit`` library. Values included in
        the ``variogram_parameters`` argument, but omitted from ``initialize_variogram``, use
        default values imported from ``self.input["statistics"]["variogram"]["model_config"]``.


        See Also
        --------

        :func:`echopop.spatial.variogram.variogram` :
            Variogram model calculation
        :func:`echopop.spatial.variogram.empirical_variogram` :
            Empirical variogram calculation
        :class:`echopop.utils.validate_dict.VariogramBase` :
            Variogram model parameters
        :class:`echopop.utils.validate.VariogramOptimize` :
            Variogram model parameter optimization arguments
        :class:`echopop.utils.validate.VariogramInitial` :
            Variogram model initialization arguments
        :class:`lmfit.parameter.Parameters` :
            Variogram parameter optimization leverages the ``Parameters``
            class from ``lmfit`` for model optimization
        """

        # Check dataset integrity
        dataset_integrity(self.input, analysis="variogram")

        # Validate "variable" input
        if variable not in ["biomass"]:
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

    def kriging_analysis(
        self,
        cropping_parameters: Dict[str, Any] = {},
        kriging_parameters: Dict[str, Any] = {},
        coordinate_transform: bool = True,
        extrapolate: bool = False,
        best_fit_variogram: bool = False,
        variable: Literal["biomass"] = "biomass",
        variogram_parameters: Optional[Dict[str, Any]] = None,
        verbose: bool = True,
    ):
        """
        Interpolates biomass data using ordinary kriging

        Parameters
        ----------

        cropping_parameters: Dict[str, Any]
            An optional dictionary containing user-defined arguments for how the kriging mesh will
            be cropped. This can include the following arguments:

            - **crop_method: Literal["transect_ends", "convex_hull"]** \n
            The method used for determining the survey extent when ``extrapolate=True``. The
            options for this include:

                - *"transect_ends"* \n
                Interpolate the eastern and western extents of transect lines over a latitude
                resolution (*see `latitude_resolution`*)


                - *"convex_hull"* \n
                Define the survey extent by computing the convex hull of the
                extents of all transect lines

            - **num_nearest_transects: posint** \n
            Combines the polygonal extents of the *n*-nearest transect lines. Only used when
            ``crop_method="convex_hull"``

            - **mesh_buffer_distance: realposfloat** \n
            A buffer (in nmi) added to the survey extent polygon generated when
            ``crop_method="convex_hull"``

            - **latitude_resolution: realposfloat** \n
            The latitude resolution used for interpolated the eastern and  western transect line
            extents when ``crop_method="transect_ends"``

            - **bearing_tolerance: realcircle** \n
            An angular tolerance (in degrees) used for grouping transect lines based on their
            respective bearings for interpolating the survey extent when
            ``crop_method="transect_ends"``

        extrapolate: bool
            Project the kriged estimates over the entire kriging mesh instead of only within the
            survey extent

        kriging_parameters: Dict[str, Any]
            An optional dictionary containing user-defined kriging parameters that will be used for
            fitting and interpolation. This can include the following arguments:

            - **anisotropy: realposfloat** \n
            The relative magnitude of directionality of the spatially autocorrelated process.  It
            is assumed that variogram parameters (e.g. nugget effect, sill) are the same in  all
            directions and therefore considered to be isotropic (``anisotropy=0.0``)

            - **correlation_range: Optional[realposfloat]** \n
            The relative length scale, or range, at which the autocorrelation between lag distances
            no longer increases and becomes asymptotic

            - **kmax: posint** \n
            The maximum number of nearest neighbors required for including values for kriging
            detected within the search radius

            - **kmin: posint** \n
            The minimum number of nearest neighbors required for including values for kriging
            within the search radius

            - **search_radius: Optional[realposfloat]** \n
            The adaptive search radius that identifies the *k*-nearest neighbors around each
            georeferenced value that are subsequently kriged

        coordinate_transform: bool
            Standardize and transform the spatial coordinates from latitude and longitude to 'x'
            and 'y'

        best_fit_variogram: bool
            Use optimized variogram parameters produced from the :meth:`.fit_variogram` method

        variable: Literal["biomass"]
            Biological variable that will be interpolated via kriging

        variogram_parameters: Optional[Dict[str, Any]]
            An optional dictionary containing user-defined variogram parameters. See
            :func:`echopop.utils.validate_dict.VariogramBase` for more details

        verbose: bool
            Print console messages and results

        Notes
        -----
        Although both ``correlation_range`` and ``search_radius`` are considered to be optional, a
        value for at least one of these variables must be supplied, otherwise a ``ValueError`` will
        be raised.

        See Also
        --------
        :meth:`.fit_variogram` :
            Variogram fitting method
        :func:`echopop.utils.validate_dict.VariogramBase` :
            Variogram model parameters

        """

        # Check dataset integrity
        dataset_integrity(self.input, analysis="kriging")

        # Populate settings dictionary with input argument values/entries
        self.analysis["settings"].update(
            {
                "kriging": {
                    "best_fit_variogram": best_fit_variogram,
                    "cropping_parameters": {**cropping_parameters},
                    "extrapolate": extrapolate,
                    "kriging_parameters": {**kriging_parameters},
                    "standardize_coordinates": coordinate_transform,
                    "variable": variable,
                    "verbose": verbose,
                },
            },
        )

        # Inherited settings/configurations (contingent on previously executed methods)
        self.analysis["settings"]["kriging"].update(
            {
                # ---- From `self.config`
                "projection": self.config["geospatial"]["init"],
                # ---- From `self.transect_analysis` settings
                "exclude_age1": self.analysis["settings"]["transect"]["exclude_age1"],
                "stratum": self.analysis["settings"]["transect"]["stratum"],
            },
        )

        # Calculate additional keys for the settings
        self.analysis["settings"]["kriging"].update(
            {
                "stratum_name": (
                    "stratum_num"
                    if self.analysis["settings"]["kriging"]["stratum"] == "ks"
                    else "inpfc"
                ),
                "variogram_parameters": (
                    {
                        **self.input["statistics"]["variogram"]["model_config"],
                        **variogram_parameters,
                    }
                    if (
                        variogram_parameters
                        and "model_config" in self.input["statistics"]["variogram"]
                    )
                    else (
                        {
                            **self.input["statistics"]["variogram"]["model_config"],
                            **{"model": ["exponential", "bessel"], "n_lags": 30},
                        }
                        if (
                            not variogram_parameters
                            and "model_config" in self.input["statistics"]["variogram"]
                        )
                        else (
                            {
                                **self.results["variogram"]["model"],
                                **self.results["variogram"]["model_fit"],
                            }
                            if best_fit_variogram is True
                            else {}
                        )
                    )
                ),
            }
        )
        # ---- Further append variogram parameters if they were ran
        if "variogram" in self.analysis["settings"]:
            self.analysis["settings"]["kriging"]["variogram_parameters"].update(
                **self.analysis["settings"]["variogram"]
            )

        # Run kriging analysis
        # ----> Generates a georeferenced dataframe for the entire mesh grid, summary statistics,
        # ----> and adds intermediate data products to the analysis attribute
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

    def plot(
        self,
        kind: Literal["age_length_distribution", "mesh", "transect"],
        variable: str,
        plot_parameters: Dict[str, Any] = {},
        plot_type: Optional[Literal["heatmap", "hexbin", "scatter", "pcolormesh"]] = None,
    ):
        """
        Plotting method for visualizing results

        Note
        ----------
        *See `Plotting parameterization` section for more details*

        Parameters
        ----------
        kind: Literal["age_length_distribution", "mesh", "transect"]
            The 'kind' of plot and dataset that will be visualized. Possible plots include:

            - `age_length_distribution`: two-dimensional heatmap of population estimates
                as a function of age (x-axis) and length (y-axis).

            - `mesh`:  map displaying kriged estimates distributed across the spatial kriging
                mesh.

            - `transect`: bubbleplot of various along-transect population estimates.

        variable: str
            The variable used for plotting. *See the `Plotting Variables` section for more details*.
        plot_parameters: Dict[str, Any]
            A dictionary comprising various plotting parameters. *See the `Plotting parameters`
            for more details*.
        plot_type: Optional[Literal["heatmap", "hexbin", "scatter", "pcolormesh"]]
            The type of plot. Options are specific to each `kind`:

            - `age_length_distribution`: ["heatmap"]
            - `mesh`: ["hexbin", "scatter", "pcolormesh"]
            - `transect`: ["scatter"]

        Other Parameters
        ----------------
        Valid variables depend on the input for `kind`:

        "age_length_distribution":

            - **"abundance"**: Animal abundance.
            - **"biomass"**: Animal biomass.

        "mesh":

            - **"biomass"**: Kriged animal biomass.
            - **"biomass_density"**: Kriged animal biomass density.
            - **"kriged_variance"**: Kriging variance, which represents the spatial uncertainty in
            prediced kriged estimates.
            - **"kriged_cv"**: Coefficient of variation calculated for each mesh node.
            - **"local_variance"**: Sample variance of local all transect values within the search
            radius of each mesh node.


        "transect":

            - **"abundance"**: Animal abundance.
            - **"abundance_female"**/**"abundance_male"**: Sexed animal abundance.
            - **"biomass"**: Animal biomass.
            - **"biomass_female"**/**"biomass_male"**: Sexed animal biomass.
            - **"biomass_density"**: Animal biomass density.
            - **"biomass_density_female"**/**"biomass_density_male"**: Sexed animal biomass density.
            - **"nasc"**: Nautical area scattering coefficient (NASC).
            - **"number_density"**: Animal number density.
            - **"number_density_female"**/**"number_density_male"**: Sexed animal number density.

        The format for `plotting_parameters` is:

        axis_limits: Optional[Dict[str, Any]]

            A dictionary that contains limits for the x- and y-axes. This should contain two
            nested dictionaries for `x` and `y`:

            - `xmin`/`ymin` (float): Minimum x- and y-axis values.
            - `xmax`/`ymax` (float): Maximum x- and y-axis values.
            - `left`/`right` (float): The left- and right-most axis values. This should not be
            included if `xmin`/`xmax`/`ymin`/`ymax` are included.

        geo_config: Optional[Dict[str, Any]]
            A dictionary that contains three keyword arguments:

            - `coastline` (`cartopy.feature.Feature`): A `cartopy.feature.Feature` object that
            includes land and/or coastline information used for geospatial plots.

            - `init` (str): The initial geospatial projection code (e.g. "EPSG:4326").

            - `plot_projection` (`cartopy.crs.Projection`): A `cartopy.crs.Projection`
            projection object used for plotting geospatial data.

        grid_heatmap: Optional[bool]
            An optional boolean value for adding correctly spaced gridlines to the biological
            heatmap. It is otherwise not used for the spatial plots.
        sex: Optional[Literal["all", "female", "male"]]
            An optional argument for defining which fish sex to plot for the biological
            heatmap. It is otherwise not used for the spatial plots.
        log_base: Optional[float]
            Base for rescaling the plot colormap via logarithmic transformation.
        cmap: Optional[str]
            Plotting colormap (e.g. "viridis").
        vmin, vmax: Optional[float]
            The data range used to rescale the colormap in linear or logarithmic space.
        **kwargs: Dict[str, Any]
            Plotting functions can accept additional arguments used within
            :func:`matplotlib.pyplot.plot`

        See Also
        ----------------
        :func:`matplotlib.pyplot.plot`
            For more details on `matplotlib` plotting function keyword arguments. Keep in mind that
            these should be used in caution as they may interfere with the primary variables
            defined for `plotting_parameters`.

        """

        # Get associated plotting function information
        plot_info = egp.PLOT_MAP(self, kind)

        # Initialize 'parameters' dictionary
        parameters = plot_parameters.copy()

        # Proceed with plotting
        # ---- Type: spatial
        if plot_info["type"] == "spatial":
            # ---- Get the geospatial configuration
            geo_config = self.config["geospatial"].copy()
            # ---- Create copy of user-defined geospatial configuration
            geo_param = parameters.get("geo_config", {})
            # ---- Update the parameterization
            parameters.update({"geo_config": {**geo_config, **geo_param}})

        # Add the primary arguments into the dictionary
        parameters.update(dict(kind=kind, plot_type=plot_type, variable=variable))

        # Prepare plotting parameters
        validated_parameters = egp.validate_plot_args(**parameters)

        # Plot
        plot_info.get("function")(plot_info.get("data"), **validated_parameters)

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
