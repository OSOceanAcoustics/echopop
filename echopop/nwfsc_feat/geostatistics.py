import inspect
from typing import Any, Callable, Dict, Optional, Tuple

import pandas as pd
from lmfit import Parameters

from . import spatial
from .projection import reproject_dataset


class Geostats:
    """
    Class for performing geostatistics including variogram analysis and interpolation via kriging.

    This class provides a comprehensive workflow for geostatistical analysis including:
    - Coordinate projection and standardization
    - Mesh cropping using various algorithms
    - Empirical variogram calculation
    - Theoretical variogram model fitting
    - Kriging interpolation
    - Result projection and processing

    Parameters
    ----------
    data_df : pd.DataFrame
        DataFrame containing survey data with spatial coordinates and variables.
    mesh_df : pd.DataFrame
        DataFrame containing mesh/grid points for kriging interpolation.
    kriging_params : Dict[str, Any], default={}
        Dictionary containing kriging parameters (e.g., kmin, kmax, search_radius).
    variogram_params : Dict[str, Any], default={}
        Dictionary containing variogram parameters (e.g., model, nugget, sill, range).
    projection : str, default="epsg:4326"
        Initial coordinate reference system of the input data.
    coordinate_names : Tuple[str, str], default=("longitude", "latitude")
        Names of the coordinate columns in the input data.

    Attributes
    ----------
    data_df : pd.DataFrame
        Survey data with coordinates and variables.
    mesh_df : pd.DataFrame
        Mesh/grid data for interpolation.
    projection : str
        Current coordinate reference system.
    coordinates : Tuple[str, str]
        Original coordinate column names.
    projection_coordinates : Tuple[str, str]
        Projected coordinate column names (typically ("x", "y") after projection).
    variogram_params : Dict[str, Any]
        Variogram parameters and fitted values.
    kriging_params : Dict[str, Any]
        Kriging parameters.
    variable : str
        Name of the variable used for variogram and kriging analysis.
    lags : np.ndarray
        Lag distances from empirical variogram calculation.
    gamma : np.ndarray
        Semivariance values from empirical variogram calculation.
    lag_counts : np.ndarray
        Number of point pairs at each lag distance.
    lag_covariance : np.ndarray
        Covariance values at each lag distance.
    best_fit_variogram_params : Dict[str, Any]
        Best-fit theoretical variogram parameters.
    variogram_fit_initial : object
        Initial variogram fit object.
    variogram_fit_optimized : object
        Optimized variogram fit object.
    survey_cv : float
        Coefficient of variation from kriging results.

    Examples
    --------
    >>> import pandas as pd
    >>> from echopop.nwfsc_feat.geostatistics import Geostats
    >>>
    >>> # Create sample data
    >>> data_df = pd.DataFrame({
    ...     'longitude': [-125.0, -124.8, -124.6],
    ...     'latitude': [48.0, 48.2, 48.4],
    ...     'biomass_density': [10.5, 12.3, 8.7]
    ... })
    >>> mesh_df = pd.DataFrame({
    ...     'longitude': [-125.1, -124.9, -124.7],
    ...     'latitude': [47.9, 48.1, 48.3]
    ... })
    >>>
    >>> # Initialize Geostats object
    >>> geo = Geostats(data_df, mesh_df)
    >>>
    >>> # Project coordinates
    >>> geo.project_coordinates(normalize=True)
    >>>
    >>> # Calculate empirical variogram
    >>> geo.calculate_empirical_variogram(variable="biomass_density")
    >>>
    >>> # Fit variogram model
    >>> from lmfit import Parameters
    >>> params = Parameters()
    >>> params.add('nugget', value=0.1)
    >>> params.add('sill', value=1.0)
    >>> params.add('correlation_range', value=0.5)
    >>> geo.fit_variogram_model(params)
    >>>
    >>> # Perform kriging
    >>> results = geo.krige(default_mesh_cell_area=1.0)
    """

    def __init__(
        self,
        data_df: pd.DataFrame,
        mesh_df: pd.DataFrame,
        kriging_params: Dict[str, Any],
        variogram_params: Dict[str, Any],
        projection: str = "epsg:4326",
        coordinate_names: Tuple[str, str] = ("longitude", "latitude"),
    ):

        # Survey data
        self.data_df = data_df.copy()

        # Interpolation mesh
        self.mesh_df = mesh_df.copy()

        # Kriging parameters
        self.kriging_params = kriging_params.copy()

        # Variogram parameters
        self.variogram_params = variogram_params.copy()

        # Projection
        self.projection = projection

        # Coordinate names
        self.coordinates = coordinate_names

        # Projection coordinate names (will be updated after projection)
        self.projection_coordinates = coordinate_names

    def project_coordinates(
        self,
        x_offset: float = 0.0,
        y_offset: float = 0.0,
        normalize: bool = True,
        reference_df: Optional[pd.DataFrame] = None,
        delta_x: Optional[float] = None,
        delta_y: Optional[float] = None,
        crs_out: Optional[str] = None,
    ) -> None:
        """
        Project or standardize coordinates for the stored datasets.

        Parameters
        ----------
        x_offset : float, default=0.0
            Offset to apply to the x-coordinates during affine normalization.
        y_offset : float, default=0.0
            Offset to apply to the y-coordinates during affine normalization.
        coordinate_names : Tuple[str, str], optional
            Names of the coordinate columns. If None, uses self.coordinates.
        normalize : bool, default=True
            If True, performs affine normalization. If False, performs CRS projection.
        reference_df : pd.DataFrame, optional
            Reference DataFrame with coordinates for interpolation, used as
            an additional offset during normalization.
        delta_x : float, optional
            Total x-axis distance used for standardizing coordinates.
        delta_y : float, optional
            Total y-axis distance used for standardizing coordinates.
        crs_out : str, optional
            Target coordinate reference system for projection (required when normalize=False).

        Returns
        -------
        None
            Updates self.data_df and self.mesh_df with projected coordinates.

        Notes
        -----
        If normalize=True, performs affine normalization of coordinates to survey-relative
        standardized coordinates ('x', 'y') using offsets and scaling.

        If normalize=False, projects coordinates using a formal map projection via GeoPandas
        to the CRS specified by crs_out.

        Examples
        --------
        >>> geo = Geostats(data_df, mesh_df)
        >>> # Affine normalization (survey-relative coordinates)
        >>> geo.project_coordinates(x_offset=-124.78, y_offset=45.0)
        >>> print(geo.data_df[['x', 'y']].head())
        >>> # CRS projection (e.g., to UTM)
        >>> geo.project_coordinates(normalize=False, crs_out="EPSG:32610")
        >>> print(geo.data_df[['x', 'y']].head())
        """

        # Affine normalization, if required
        if normalize:
            # ---- Standardize the input data coordinates
            self.data_df, delta_x_new, delta_y_new = spatial.standardize_coordinates(
                data_df=self.data_df,
                reference_df=reference_df,
                x_offset=x_offset,
                y_offset=y_offset,
                delta_x=delta_x,
                delta_y=delta_y,
                coordinate_names=self.coordinates,
            )
            # ---- Apply the input data x- and y-coordinate intervals to the projection mesh
            self.mesh_df, _, _ = spatial.standardize_coordinates(
                data_df=self.mesh_df,
                reference_df=reference_df,
                x_offset=x_offset,
                y_offset=y_offset,
                delta_x=delta_x_new,
                delta_y=delta_y_new,
                coordinate_names=self.coordinates,
            )
        # Standard geospatial reprojection
        else:
            if crs_out is None:
                raise ValueError("crs_out must be specified when normalize=False")
            # ---- Input data coordinates
            self.data_df = reproject_dataset(
                data_df=self.data_df,
                projection=self.projection,
                coordinate_names=self.coordinates,
                crs_out=crs_out,
            )
            # ---- Mesh
            self.mesh_df = reproject_dataset(
                data_df=self.mesh_df,
                projection=self.projection,
                coordinate_names=self.coordinates,
                crs_out=crs_out,
            )

        # Adjust the projection names
        self.projection_coordinates = ("x", "y")

    def crop_mesh(self, crop_function: Callable, **kwargs) -> None:
        """
        Crop the mesh using a specified cropping function and update self.mesh_df

        Parameters
        ----------
        crop_function : Callable
            The mesh cropping function to use (e.g., mesh.hull_crop, mesh.transect_ends_crop).
        **kwargs
            Additional keyword arguments to pass to the cropping function.

        Returns
        -------
        None
            Updates self.mesh_df with the cropped mesh.

        Notes
        -----
        If the cropping function returns a tuple, the first element is assumed to be the cropped
        mesh.

        Examples
        --------
        >>> geo.crop_mesh(
        ...     mesh.hull_crop,
        ...     transect_df=geo.data_df,
        ...     mesh_df=geo.mesh_df,
        ...     num_nearest_transects=3,
        ...     mesh_buffer_distance=2.5,
        ...     projection="epsg:4326"
        ... )
        >>> geo.crop_mesh(
        ...     mesh.transect_ends_crop,
        ...     transect_df=geo.data_df,
        ...     mesh_df=geo.mesh_df,
        ...     latitude_resolution=1.25/60,
        ...     transect_mesh_region_function=FEAT.transect_mesh_region_2019
        ... )
        """

        # Get correct arguments
        args = inspect.signature(crop_function).parameters

        # Inject coordinate names, if needed
        if "coordinate_names" in args and "coordinate_names" not in kwargs:
            kwargs["coordinate_names"] = self.projection_coordinates

        # Inject projection, if needed
        if "projection" in args and "projection" not in kwargs:
            kwargs["projection"] = self.projection

        # Call the cropping function
        result = crop_function(self.data_df, self.mesh_df, **kwargs)

        # Update the cropped mesh DataFrame
        self.mesh_df = result[0] if isinstance(result, tuple) else result

    def calculate_empirical_variogram(
        self,
        variable: str,
        azimuth_filter: bool = True,
        azimuth_angle_threshold: float = 180.0,
        force_lag_zero: bool = True,
    ) -> None:
        """
        Compute the empirical variogram from transect data

        Parameters
        ----------
        variable : str, default = 'biomass_density'
            The variable used for computing the empirical variogram (e.g. 'biomass_density'), which
            must exist as a column in self.data_df.
        azimuth_filter : bool
            When True, a 2D array of azimuth angles are generated. This subsequent array represents
            the relative azimuth angles between spatial points, and can serve as a filter for cases
            where a high degree of directionality is assumed. This accompanies the argument
            'azimuth_angle_threshold' that defines the threshold azimuth angle.
        azimuth_angle_threshold : float
            This threshold is used for filtering the azimuth angles.

        force_lag_zero : bool, default = True
            When True, the nugget effect is assumed to be 0.0 for the empirical variogram. This
            adds lag 0 to the subsequent array outputs where semivariance (or 'gamma_h') is also
            equal to 0.

        Returns
        -------
        None
            Adds or updates the attributes self.lags, self.gamma, self.lag_counts, and
            self.lag_covariance. These represent the lag intervals, semivariance, lag counts, and
            mean lag covariance between head and tail points, respectively.
        """

        # Assign variable
        self.variable = variable

        # Subset stored variogram inputs
        args = inspect.signature(spatial.empirical_variogram).parameters

        # Combine variogram parameter kwargs, if needed
        empirical_variogram_args = {
            "azimuth_filter": azimuth_filter,
            "azimuth_angle_threshold": azimuth_angle_threshold,
            "force_lag_zero": force_lag_zero,
            **{k: v for k, v in self.variogram_params.items() if k in args},
        }

        # Compute the empirical variogram
        self.lags, self.gamma, self.lag_counts, self.lag_covariance = spatial.empirical_variogram(
            transect_df=self.data_df,
            variable=variable,
            coordinates=self.projection_coordinates,
            **empirical_variogram_args
        )

    def fit_variogram_model(
        self,
        parameter_values: Parameters,
        optimizer_kwargs: Dict[str, Any] = {},
    ) -> None:
        """
        Fit a theoretical variogram model to the empirical variogram.

        Parameters
        ----------
        parameter_values : Parameters
            An lmfit Parameters object containing the initial parameter values and constraints
            for the variogram model fitting.
        optimizer_kwargs : Dict[str, Any], default={}
            Additional keyword arguments to pass to the optimizer.

        Returns
        -------
        None
            Updates self.best_fit_variogram_params, self.variogram_fit_initial,
            self.variogram_fit_optimized, and self.variogram_params attributes.

        Notes
        -----
        This method requires that calculate_empirical_variogram() has been called first
        to compute the empirical variogram (self.lags, self.gamma, self.lag_counts).
        """

        # Fit the optimized theoretical variogram model
        self.best_fit_variogram_params, self.variogram_fit_initial, self.variogram_fit_optimized = (
            spatial.fit_variogram(
                lags=self.lags,
                lag_counts=self.lag_counts,
                gamma=self.gamma,
                model=self.variogram_params["model"],
                variogram_parameters=parameter_values,
                optimizer_kwargs=optimizer_kwargs,
            )
        )

        # Update the variogram parameters
        self.variogram_params = {**self.variogram_params, **self.best_fit_variogram_params}

        # Update the kriging parameters, if needed
        if "search_radius" in self.kriging_params:
            self.kriging_params["search_radius"] = self.variogram_params["correlation_range"] * 3

    def krige(
        self,
        default_mesh_cell_area: float,
        adaptive_search_strategy: Callable = spatial.uniform_search_strategy,
    ) -> pd.DataFrame:
        """
        Perform kriging interpolation and project the results.

        Parameters
        ----------
        default_mesh_cell_area : float
            The default area of each mesh cell used for projecting kriging results.
        adaptive_search_strategy : Callable, default=spatial.uniform_search_strategy
            The search strategy function to use for adaptive kriging.

        Returns
        -------
        pd.DataFrame
            A DataFrame containing the kriged results with the variable column renamed
            to match self.variable.

        Notes
        -----
        This method requires that both calculate_empirical_variogram() and
        fit_variogram_model() have been called first. The method performs ordinary
        kriging followed by projection of the results.
        """

        # Initial kriging [assumes ordinary kriging]
        kriged_estimates = spatial.krige(
            transect_df=self.data_df,
            kriging_mesh=self.mesh_df,
            coordinate_names=self.projection_coordinates,
            variable=self.variable,
            kriging_parameters=self.kriging_params,
            variogram_parameters=self.variogram_params,
            adaptive_search_strategy=adaptive_search_strategy,
        )

        # Project the results
        kriged_results, self.survey_cv = spatial.project_kriging_results(
            kriged_estimates=kriged_estimates,
            kriging_mesh=self.mesh_df,
            transect_df=self.data_df,
            variable=self.variable,
            default_mesh_cell_area=default_mesh_cell_area,
        )

        # Rename the variable and return
        return kriged_results.rename(columns={"estimate": self.variable})
