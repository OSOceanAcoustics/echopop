import copy
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from ...analysis import krige
from ...spatial.projection import transform_geometry
from ...utils import message as em
from .data_ingestion import (
    ingest_inversion_files,
    load_analysis_files,
    load_inversion_configuration,
)
from .math import reduce_dataset
from .message import inversion_kriging_results_msg
from .optimize import estimate_population, group_optimizer, prepare_minimizer
from .spatial import inversion_variogram_analysis, krige_inverted_data


class AcousticInversion:
    """
    Annotation
    """

    def __init__(
        self,
        inversion_configuration_filepath: Union[str, Path],
    ):

        # Load the configuration file
        self.inversion_config = load_inversion_configuration(inversion_configuration_filepath)

        # Copy `inversion_config` for now [TEMPORARY -- JURY-RIGGED TO ENABLE PLOTTING LATER ON]
        self.config = self.inversion_config

        # Initialize the dataset attribute
        self.input = {}

        # Load the metadata and dataset
        self.input["metadata"], self.input["sv_df"] = ingest_inversion_files(
            **self.inversion_config["inversion_file_settings"],
            **self.inversion_config["processing_parameters"]["acoustic_files"],
        )

        # Update the parameterization to include the literal center frequencies ingested
        self.inversion_config["processing_parameters"]["acoustic_models"].update(
            {
                "center_frequencies": self.input["sv_df"].columns.to_numpy() * 1e3,
            }
        )

        # Read in additional analysis files
        load_analysis_files(self.input, self.inversion_config)

        # Analysis
        self.analysis = {}

        # Initialize results
        self.results = {}

    def invert_population(
        self,
        verbose: bool = True,
        subset_dataset: Optional[Dict[str, Any]] = None,
    ):
        """
        Invert population estimates from active acoustic backscatter measurements
        """

        # Initialize the dictionary for the inversion results
        self.results["inversion"] = {}

        # Create a copy of the ingested dataset
        measurements_df = self.input["sv_df"].copy()

        # Subset data, if needed
        if subset_dataset:
            self.analysis["inversion"], measurements_df = reduce_dataset(
                self.input["metadata"], measurements_df, subset_dataset
            )
        else:
            self.analysis["inversion"] = copy.deepcopy(self.input["metadata"])

        # Get specific keyword arguments
        # ---- Center frequencies
        center_frequencies = self.inversion_config["processing_parameters"]["acoustic_models"][
            "center_frequencies"
        ]
        # ---- Minimum Sv threshold
        sv_threshold = self.inversion_config["processing_parameters"]["acoustic_files"][
            "sv_threshold"
        ]
        # ---- Aggregate
        aggregate = self.inversion_config["processing_parameters"]["acoustic_files"]["aggregate"]

        # Prepare the optimizer parameterization
        measurements_df = prepare_minimizer(
            data_df=measurements_df,
            center_frequencies=center_frequencies,
            aggregate=aggregate,
            sv_threshold=sv_threshold,
            scattering_parameters=self.inversion_config["scattering_parameters"],
            simulation_parameters=self.inversion_config["processing_parameters"]["simulation"],
            processing_parameters=self.inversion_config["processing_parameters"]["acoustic_models"],
        )

        # Run the inversion algorithm
        inversion_results = measurements_df.apply(
            group_optimizer,
            axis=1,
            args=(
                self.inversion_config["processing_parameters"]["acoustic_models"],
                self.inversion_config["optimization_parameters"],
                self.inversion_config["processing_parameters"]["simulation"],
                center_frequencies,
                sv_threshold,
                verbose,
            ),
            result_type="expand",
        )

        # Store the best-fit inverted results
        self.results["inversion"]["inverted_parameters_df"] = inversion_results

        # Generate population results
        self.results["inversion"]["transect_df"] = estimate_population(
            inversion_results,
            density_sw=(
                self.inversion_config["processing_parameters"]["acoustic_models"]["density_sw"]
            ),
            aggregate=aggregate,
            **self.analysis["inversion"],
            **self.inversion_config["processing_parameters"]["inversion"],
        )

        # Copy these results to a "transect" key to enable plotting later on
        # [TEMPORARY -- JURY-RIGGED TO ENABLE PLOTTING LATER ON]
        self.analysis["transect"] = {}
        self.analysis["transect"]["acoustics"] = {}
        self.analysis["transect"]["acoustics"]["adult_transect_df"] = self.results["inversion"][
            "transect_df"
        ]
        # ---- Adjust column name
        self.analysis["transect"]["acoustics"]["adult_transect_df"]["biomass_density"] = (
            self.results["inversion"]["transect_df"]["biomass_areal_density"]
        )

    def fit_inversion_variogram(
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
        verbose: bool = True,
    ) -> None:
        """
        Compute the best-fit variogram parameters for inverted biomass densities

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

        # Initialize Survey-class object
        self.analysis.update({"variogram": {}, "settings": {}})

        # Parameterize analysis settings that will be applied to the variogram fitting and analysis
        self.analysis["settings"].update(
            {
                "variogram": {
                    "azimuth_range": (
                        self.inversion_config["variogram_parameters"]["azimuth_range"]
                        if "azimuth_range" in self.inversion_config["variogram_parameters"]
                        else azimuth_range
                    ),
                    "fit_parameters": (
                        initialize_variogram.keys()
                        if isinstance(initialize_variogram, dict)
                        else initialize_variogram
                    ),
                    "force_lag_zero": force_lag_zero,
                    "model": (
                        self.inversion_config["variogram_parameters"]["model"]
                        if "model" in self.inversion_config["variogram_parameters"]
                        else model
                    ),
                    "standardize_coordinates": standardize_coordinates,
                    "variable": "biomass_areal_density",
                    "verbose": verbose,
                }
            }
        )

        # Create a copy of the existing variogram settings
        variogram_parameters.update({"model": model, "n_lags": n_lags})

        # Create optimization settings dictionary
        # ---- Add to settings
        self.analysis["settings"]["variogram"].update({"optimization": optimization_parameters})

        # Append `kriging_parameters` to the settings dictionary
        if standardize_coordinates:
            self.analysis["settings"]["variogram"].update(
                {"kriging_parameters": self.inversion_config["kriging_parameters"]}
            )

        # Create a copy of the existing variogram settings
        default_variogram_parameters = self.inversion_config["variogram_parameters"].copy()
        # ---- Update model, n_lags
        default_variogram_parameters.update({"model": model, "n_lags": n_lags})

        # Standardize the transect coordinates, if necessary
        # ---- Create copy of the transect data
        transect_data = self.results["inversion"]["transect_df"].copy()
        # ---- Standardize
        if standardize_coordinates:
            # ---- Transform geometry
            transect_data, _, _ = transform_geometry(
                transect_data,
                self.input["isobath_reference"],
                self.analysis["settings"]["variogram"],
            )
            # ---- Print message if verbose
            if verbose:
                # ---- Print alert
                print(
                    "Longitude and latitude coordinates (WGS84) converted to standardized "
                    "coordinates (x and y)."
                )
        else:
            # ---- x
            transect_data["x"] = "longitude"
            # ---- y
            transect_data["y"] = "latitude"

        # Store the standardized transect data
        self.analysis["variogram"]["transect_df"] = transect_data

        # Find the best-fit variogram parameters
        best_fit_variogram = inversion_variogram_analysis(
            transect_data,
            variogram_parameters,
            default_variogram_parameters,
            optimization_parameters,
            initialize_variogram,
            self.analysis["settings"]["variogram"],
        )

        # Add "partial" results to analysis attribute
        self.analysis["variogram"].update(
            {
                "model": default_variogram_parameters["model"],
                "initial_fit": best_fit_variogram["initial_fit"],
                "optimized_fit": best_fit_variogram["optimized_fit"],
            }
        )

        # Add variogram result
        self.results.update(
            {
                "variogram": {
                    "model_fit": best_fit_variogram["best_fit_parameters"],
                    "n_lags": default_variogram_parameters["n_lags"],
                    "model": default_variogram_parameters["model"],
                }
            }
        )

        # Print result if `verbose == True`
        if verbose:
            em.variogram_results_msg(self.analysis["variogram"])

    def inversion_kriging_analysis(
        self,
        cropping_parameters: Dict[str, Any] = {},
        kriging_parameters: Dict[str, Any] = {},
        extrapolate: bool = False,
        verbose: bool = True,
    ):
        """
        Interpolates biomass data from inverted population estimates using ordinary kriging

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

        # Populate settings dictionary with input argument values/entries
        self.analysis["settings"].update(
            {
                "kriging": {
                    "cropping_parameters": {"crop_method": "convex_hull", **cropping_parameters},
                    "extrapolate": extrapolate,
                    "kriging_parameters": {
                        **self.inversion_config["kriging_parameters"],
                        **kriging_parameters,
                    },
                    "projection": self.inversion_config["geospatial"]["init"],
                    "standardize_coordinates": (
                        self.analysis["settings"]["variogram"]["standardize_coordinates"]
                    ),
                    "variable": "biomass_areal_density",
                    "variogram_parameters": {
                        "model": self.results["variogram"]["model"],
                        "n_lags": self.results["variogram"]["n_lags"],
                        **self.results["variogram"]["model_fit"],
                    },
                    "verbose": verbose,
                },
            },
        )

        # Initialize the analysis attribute
        self.analysis["kriging"] = {}

        # Standardize the transect coordinates, if necessary
        # ---- Create copy of the transect data
        mesh_data = self.input["mesh"].copy()
        # ---- Standardize
        if self.analysis["settings"]["variogram"]["standardize_coordinates"]:
            # ---- Transform geometry
            mesh_data, _, _ = transform_geometry(
                mesh_data, self.input["isobath_reference"], self.analysis["settings"]["variogram"]
            )
            # ---- Print message if verbose
            if verbose:
                # ---- Print alert
                print(
                    "Longitude and latitude coordinates (WGS84) converted to standardized "
                    "coordinates (x and y)."
                )
        else:
            # ---- x
            mesh_data["x"] = mesh_data["longitude"]
            # ---- y
            mesh_data["y"] = mesh_data["latitude"]

        # Get the transect data
        transect_df = self.analysis["variogram"]["transect_df"].copy()

        # Run kriging analysis
        # ----> Generates a georeferenced dataframe for the entire mesh grid, summary statistics,
        # ----> and adds intermediate data products to the analysis attribute
        # ---- Run kriging algorithm
        kriged_results = krige_inverted_data(
            transect_df, mesh_data, self.analysis["kriging"], self.analysis["settings"]["kriging"]
        )

        # Save the results to the `results` attribute
        self.results.update({"kriging": kriged_results})

        # Print result if `verbose == True`
        if verbose:
            inversion_kriging_results_msg(
                self.results["kriging"], self.analysis["settings"]["kriging"]
            )
