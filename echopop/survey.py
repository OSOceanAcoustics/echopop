import copy
from pathlib import Path
from typing import Literal, Optional, Union

from .analysis import (
    acoustics_to_biology,
    apportion_kriged_values,
    krige,
    process_transect_data,
    stratified_summary,
)
from .core import DATA_STRUCTURE
from .utils import load as el, message as em


class Survey:
    """
    echopop base class that imports and prepares parameters for
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
        various downstream and internal functions.
    input: dict
        Input data.
    analysis: dict
        Analysis variables and intermediate data products.
    results: dict
        Analysis results.
    """

    def __init__(
        self, init_config_path: Union[str, Path], survey_year_config_path: Union[str, Path]
    ):
        # Initialize `meta` attribute
        self.meta = copy.deepcopy(DATA_STRUCTURE["meta"])

        # Loading the configuration settings and definitions that are used to
        # initialize the Survey class object
        self.config = el.load_configuration(Path(init_config_path), Path(survey_year_config_path))

        # Loading the datasets defined in the configuration files
        self.input = el.load_survey_data(self.config)

        # Initialize the `analysis` data attribute
        self.analysis = copy.deepcopy(DATA_STRUCTURE["analysis"])

        # Initialize the `results` data attribute
        self.results = copy.deepcopy(DATA_STRUCTURE["results"])

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

    # !!! TODO: develop different name for "crop_method = 'interpolation'"
    def kriging_analysis(
        self,
        bearing_tolerance: float = 15.0,
        coordinate_transform: bool = True,
        crop_method: Literal["interpolation", "convex_hull"] = "interpolation",
        extrapolate: bool = False,
        kriging_parameters: Optional[dict] = None,
        latitude_resolution: float = 1.25,
        mesh_buffer_distance: float = 1.25,
        num_nearest_transects: int = 4,
        projection: Optional[str] = None,
        stratum: str = "ks",
        variable: str = "biomass_density",
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

        # Parameterize analysis settings that will be applied to the stratified analysis
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
