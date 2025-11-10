import abc
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd
from lmfit import Parameters
from pydantic import ValidationError


class InversionBase(abc.ABC):
    """
    Abstract base class for handling acoustic inversion methods.

    This class provides a framework for different types of acoustic inversions
    by establishing common interfaces and shared functionality for parameter
    management and stratification handling.

    Parameters
    ----------
    model_parameters : Dict[str, Any]
        Dictionary containing model configuration parameters. Common keys include:
        - 'stratify_by': str or List[str] - columns to stratify by
        - 'strata': array-like - specific strata to process
        - 'impute_missing_strata': bool - whether to impute missing strata

    Attributes
    ----------
    model_params : Dict[str, Any]
        Processed model parameters
    inversion_method : str
        String identifier for the specific inversion method (set by subclasses)

    Examples
    --------
    >>> # Example parameters for length-TS inversion
    >>> params = {
    ...     "ts_length_regression": {"slope": 20.0, "intercept": -68.0},
    ...     "stratify_by": "stratum_ks",
    ...     "strata": [1, 2, 3, 4, 5],
    ...     "impute_missing_strata": True
    ... }
    >>>
    >>> # Create concrete inversion class (subclass)
    >>> inverter = InversionLengthTS(params)
    >>> print(inverter.inversion_method)
    length_TS_regression

    Notes
    -----
    This is an abstract base class and cannot be instantiated directly.
    Subclasses must implement the abstract `invert` method.

    The class automatically converts single-string 'stratify_by' parameters
    to lists for consistent handling across different inversion methods.
    """

    def __init__(self, model_parameters):

        # Ingest model parameters
        self.model_params = model_parameters

        # Modify "stratify_by" if needed
        if "stratify_by" in self.model_params:
            if isinstance(self.model_params["stratify_by"], str):
                self.model_params["stratify_by"] = [self.model_params["stratify_by"]]

        # Initialize method
        self.inversion_method = ""

    @abc.abstractmethod
    def invert(self, df_nasc: pd.DataFrame) -> pd.DataFrame:
        """
        Perform inversion on input dataframe to convert NASC to number density.

        This abstract method must be implemented by subclasses to define the
        specific inversion algorithm (e.g., length-TS, age-TS, species-specific).

        Parameters
        ----------
        df_nasc : pd.DataFrame
            Input dataframe with NASC (Nautical Area Scattering Coefficient) values
            to perform inversion on. Must contain 'nasc' column and any required
            stratification columns.

        Returns
        -------
        pd.DataFrame
            Input dataframe with added 'number_density' column containing the
            inverted number density estimates.

        Notes
        -----
        To perform inversion on non-stratum groupings, pre-process the dataframe
        so that each row contains the minimum unit inversion will be performed on.

        Subclasses should implement the specific inversion algorithm appropriate
        for their method (length-based, age-based, etc.).
        """
        pass


class InvParameters:
    """
    Container for acoustic scattering model inversion parameters.

    This class manages parameter sets used in acoustic inversion analysis,
    providing scaling/unscaling functionality, bounds management, and
    integration with the lmfit optimization library. It serves as the
    primary interface for parameter handling in echopop inversion workflows.

    Parameters
    ----------
    parameters : Dict[str, Any]
        Dictionary of parameter specifications with each parameter containing:
        - 'value': Current parameter value
        - 'min': Lower bound (optional, default=-inf)
        - 'max': Upper bound (optional, default=+inf)
        - 'vary': Whether parameter should be optimized (optional, default=False)

    Attributes
    ----------
    parameters : Dict[str, Dict[str, Any]]
        Validated parameter dictionary
    parameter_bounds : Dict[str, Dict[str, float]]
        Original parameter bounds for unscaling operations
    scaled : bool
        Current scaling status

    Methods
    -------
    values
        Property to get current parameter values as dictionary
    bounds
        Property to get parameter bounds as dictionary of (min, max) tuples
    is_scaled
        Property to check if parameters are currently scaled
    scale()
        Scale parameters to [0,1] range and return new scaled instance
    unscale()
        Unscale parameters to original range and return new unscaled instance
    to_lmfit()
        Convert to lmfit.Parameters object for optimization
    inverse_transform(scaled_dict)
        Convert scaled parameter dictionary back to original scale
    update_bounds(bounds)
        Update the boundaries (min/max) for stored parameters
    from_series(series)
        Class method to create instance from pandas Series

    Examples
    --------
    >>> params = {
    ...     'length_mean': {'value': 25.0, 'min': 10.0, 'max': 40.0, 'vary': True},
    ...     'g': {'value': 1.02, 'min': 0.95, 'max': 1.05, 'vary': True}
    ... }
    >>> inv_params = InvParameters(params)
    >>> inv_params.scale()  # Scale to [0,1]
    >>> lmfit_params = inv_params.to_lmfit()

    Notes
    -----
    The scaling transformation uses min-max normalization:

    .. math::
        x_{scaled} = \\frac{x - x_{min}}{x_{max} - x_{min}}

    This improves optimization convergence by normalizing parameter ranges
    and reducing numerical conditioning issues in multi-parameter problems.
    """

    def __init__(
        self,
        parameters: Dict[str, Any],
    ):

        # Run initial validation
        validated_parameters = self._validate(parameters)

        # Initialize attributes
        self._unscaled_parameters = validated_parameters
        self.parameter_bounds = self._get_parameter_limits(validated_parameters)
        self._scaled = False
        self._realizations = {}

        # Store the scaled parameters
        self._scale_parameters()

        # Default parameters to unscaled, natural values
        self.parameters = validated_parameters

    def __repr__(self):
        return f"InvParameters(scaled={self._scaled}, parameters={list(self.parameters.keys())})"

    def __str__(self):
        # Format the parameters string
        parameters_str = ", ".join([f"{k}: {v}" for k, v in self.parameters.items()])

        # Return the string output
        return f"InvParameters(scaled={self._scaled}, parameters=[{parameters_str}])"

    def _validate(self, parameters):
        """
        Internal validator
        """

        # Delay import to avoid circular import issues
        from ..validators.inversion import ModelInputParameters

        # Validation step
        if not isinstance(parameters, dict):
            raise TypeError("Parameters must be a dictionary.")
        for k, v in parameters.items():
            if not isinstance(k, (int, float, str)):
                raise TypeError(f"Parameter key '{k}' must be a numeric or string.")

        # Validate parameters
        try:
            # ---- Check
            validated_parameters = ModelInputParameters.create(**parameters)
        # Break creation
        except (ValidationError, Exception) as e:
            raise e from None

        # Return parameter
        return validated_parameters

    @property
    def values(self):
        """
        Get dictionary of parameter names and their current values.

        Returns
        -------
        dict
            Dictionary mapping parameter names to their current values.
            Values are in the current scaling state (scaled or unscaled).

        Examples
        --------
        >>> params = InvParameters({'length_mean': {'value': 25.0}})
        >>> params.values
        {'length_mean': 25.0}
        """
        return {key: param["value"] for key, param in self.parameters.items()}

    @property
    def bounds(self):
        """
        Get dictionary of parameter bounds (min/max values).

        Returns
        -------
        dict
            Dictionary mapping parameter names to their bounds.
            Each bound is a dict with 'min' and 'max' keys containing
            the original unscaled bounds regardless of current scaling state.

        Examples
        --------
        >>> params = InvParameters({
        ...     'length_mean': {'value': 25.0, 'min': 10.0, 'max': 40.0}
        ... })
        >>> params.bounds
        {'length_mean': {'min': 10.0, 'max': 40.0}}
        """
        return self.parameter_bounds

    @property
    def is_scaled(self):
        """
        Check whether parameters are currently in scaled [0,1] form.

        Returns
        -------
        bool
            True if parameters are scaled to [0,1], False if in original units.

        Notes
        -----
        This property helps track the current state of parameters during
        optimization workflows where scaling/unscaling may occur multiple times.
        """
        return self._scaled

    @property
    def realizations(self):
        """
        Access the dictionary of parameter realizations.

        Returns
        -------
        dict
            Dictionary mapping realization indices (int) to Monte Carlo parameter sets.

        Examples
        --------
        >>> params.simulate_parameter_sets(mc_realizations=5) # params is an InvParameters object
        >>> realizations = params.realizations
        >>> first_real = realizations[0]  # First realization
        """
        return self._realizations

    @staticmethod
    def _get_parameter_limits(parameters):
        """
        Extract parameter bounds from validated parameter dictionary.

        Parameters
        ----------
        parameters : dict
            Dictionary of parameter specifications containing 'min' and 'max' keys

        Returns
        -------
        dict
            Dictionary mapping parameter names to their bounds with 'min' and 'max' keys

        Notes
        -----
        This static method is used internally to preserve original parameter
        bounds for later unscaling operations. It extracts only the bound
        information, ignoring other parameter attributes.
        """

        # Get the upper and lower values in case of rescaling
        parameter_limits = {
            key: {"min": value["min"], "max": value["max"]} for key, value in parameters.items()
        }

        return parameter_limits

    def _scale_parameters(self):
        """
        Internal method for scaling parameters to [0,1] range using min-max normalization.

        This method transforms parameter values to a normalized range which
        improves optimization performance by ensuring all parameters have
        similar scales and reduces numerical conditioning problems.

        Notes
        -----
        The scaling transformation is:

        .. math::
            x_{scaled} = \\frac{x - x_{min}}{x_{max} - x_{min}}

        After scaling, parameter bounds become [0,1] for all parameters.
        The original bounds are preserved in `parameter_bounds` for unscaling.

        The scaling status flag is updated to True after this operation.
        """

        # Update parameters
        scaled_params = {
            key: (
                {
                    **kwargs,
                    "value": (
                        (kwargs["value"] - kwargs["min"]) / (kwargs["max"] - kwargs["min"])
                        if kwargs["max"] != kwargs["min"]
                        else 0.0
                    ),
                    "min": 0.0,
                    "max": 1.0,
                }
                if isinstance(kwargs, dict)
                else kwargs
            )
            for key, kwargs in self._unscaled_parameters.items()
        }

        # Update the parameters
        self._scaled_parameters = scaled_params

    def scale(self):
        """
        Scale model parameters to [0, 1] range using min-max normalization.

        This method transforms all model parameter values to a normalized range which helps
        improve optimization performance by ensuring all parameters have similar scales and reduce
        numerical conditioning problems.

        Notes
        -----
        The scaling transformation is:

        .. math::
            x_{scaled} = \\frac{x - x_{min}}{x_{max} - x_{min}}

        After scaling, parameter bounds become [0,1] for all parameters.
        The original bounds are preserved in `parameter_bounds` for unscaling.

        The scaling status flag is updated to True after this operation.
        """

        # Update the `is_scaled` property
        self._scaled = True

        # Set active parameters
        self.parameters = self._scaled_parameters

    def unscale(self):
        """
        Convert scaled [0, 1] parameters back to their original units and ranges.

        This method reverses the min-max normalization, restoring parameters to their original
        physical units and boundaries. It uses the preserved bounds from `parameter_bounds` to
        perform the inverse transformation.

        Notes
        -----
        The unscaling transformation is:

        .. math::
            x_{unscaled} = x_{scaled} \\cdot (x_{max} - x_{min}) + x_{min}

        After unscaling, parameters return to their original bounds and units.
        The scaling status flag is updated to False after this operation.
        """

        # Update the `is_scaled` property
        self._scaled = False

        # Set active parameters
        self.parameters = self._unscaled_parameters

    def to_lmfit(self):
        """
        Convert internal parameters to an lmfit.Parameters object.

        This method creates an lmfit.Parameters object compatible with
        lmfit optimization routines. All parameter attributes (value, min,
        max, vary) are transferred to the lmfit format.

        Returns
        -------
        lmfit.Parameters
            Parameters object ready for use with lmfit minimization routines.
            Contains all parameter values, bounds, and vary flags.

        Notes
        -----
        The lmfit.Parameters object provides the interface for scipy.optimize
        and other optimization backends used by lmfit. Default bounds are
        0.0 (min) and np.inf (max) if not specified.

        Examples
        --------
        >>> params = InvParameters({'g': {'value': 1.02, 'min': 1.0, 'max': 1.1}})
        >>> lmfit_params = params.to_lmfit()
        >>> # Can now use with lmfit.minimize()
        """
        params = Parameters()
        for name, kwargs in self.parameters.items():
            params.add(
                name,
                value=kwargs["value"],
                min=kwargs.get("min", 0.0),
                max=kwargs.get("max", np.inf),
                vary=kwargs.get("vary", True),
            )
        return params

    def realizations_to_lmfit(self) -> Dict[str, dict]:
        """
        Convert all parameter realizations to lmfit.Parameters format.

        This method transforms the Monte Carlo parameter ensemble into formats suitable for
        optimization routines, with options for both parallel and serialized parameter handling.

        Returns
        -------
        dict
            If serialize=False: Dictionary mapping realization indices to
            lmfit.Parameters objects
            If serialize=True: Single Parameters object with parameters named
            as {param_name}_{realization_index}

        Examples
        --------
        >>> mc_params = InvParameters(base_params)
        >>> mc_params.simulate_parameter_sets(mc_realizations=3)
        >>> # Separate parameters for parallel evaluation
        >>> param_dict = mc_params.realizations_to_lmfit()

        Notes
        -----
        The serialized format is useful when evaluating multiple realizations
        simultaneously in vectorized forward models or when using optimization
        algorithms that can handle large parameter spaces efficiently.
        """

        # Convert the dictionary of realizations into lmfit.Parameters
        return {
            i: (
                lambda param_dict: (
                    (params := Parameters())
                    or [
                        params.add(
                            name,
                            value=kwargs["value"],
                            min=kwargs.get("min", 0.0),
                            max=kwargs.get("max", np.inf),
                            vary=kwargs.get("vary", True),
                        )
                        for name, kwargs in param_dict.items()
                    ]
                )
                and params
            )(param_dict)
            for i, param_dict in self._realizations.items()
        }

    def inverse_transform(self, scaled_dict: dict) -> dict:
        """
        Convert a dictionary of scaled [0,1] parameter values to original units.

        This utility method transforms scaled parameter values back to their
        original ranges using the bounds stored during initialization. It's
        useful for converting optimization results back to physical units.

        Parameters
        ----------
        scaled_dict : dict
            Dictionary of parameter names and their scaled [0,1] values

        Returns
        -------
        dict
            Dictionary of parameter names and their unscaled values in original units

        Examples
        --------
        >>> params = InvParameters({'length_mean': {'min': 10, 'max': 40}})
        >>> scaled_vals = {'length_mean': 0.5}  # Mid-range scaled value
        >>> unscaled = params.inverse_transform(scaled_vals)
        >>> unscaled['length_mean']  # Should be 25.0 (midpoint of 10-40)
        25.0

        Notes
        -----
        This method does not modify the internal parameter state, only
        converts the provided dictionary values.
        """
        return {
            k: scaled_dict[k] * (self.parameter_bounds[k]["max"] - self.parameter_bounds[k]["min"])
            + self.parameter_bounds[k]["min"]
            for k in scaled_dict
        }

    @staticmethod
    def _set_parameter_bounds(
        parameters: Dict[str, dict], bounds: Dict[str, dict]
    ) -> Dict[str, dict]:
        """
        Internal method for setting the parameter boundaries of a parameter set.
        """
        return {
            k: {
                **v,
                "min": bounds.get(k, {}).get("min", v.get("min", -np.inf)),
                "max": bounds.get(k, {}).get("max", v.get("max", np.inf)),
            }
            for k, v in parameters.items()
        }

    def update_bounds(self, bounds: Dict[str, dict]) -> None:
        """
        Add or update interval boundaries for each designated parameter.

        Parameters
        ----------
        bounds : Dict[str, dict]
            Dictionary mapping parameters to {'min': value, 'max': value}.

        Example
        -------
        >>> bounds = {'length_mean': {'min': 0.5, 'max': 2.0}, 'g': {'min': 0.9, 'max': 1.1}}
        >>> params.update_bounds(bounds)
        """

        # Map the new min/max values to the appropriate values
        updated_parameters = self._set_parameter_bounds(self._unscaled_parameters, bounds)

        # Run the validation
        validated_parameters = self._validate(updated_parameters)

        # Update the parameter set attributes
        self._unscaled_parameters = validated_parameters
        self.parameters = validated_parameters
        self._scale_parameters()
        self.parameter_bounds = self._get_parameter_limits(validated_parameters)

    @classmethod
    def from_series(cls, series: pd.Series):
        """
        Create InvParameters instance from a pandas Series.

        This class method provides a convenient way to construct InvParameters
        from pandas Series data, commonly used when reading parameters from
        DataFrames or other pandas data structures.

        Parameters
        ----------
        series : pd.Series
            Series containing parameter names as index and values as data.
            Only parameter values are extracted; bounds default to infinite.

        Returns
        -------
        InvParameters
            New instance with parameters set from series values and
            default bounds (no min/max constraints, vary=False)

        Raises
        ------
        TypeError
            If input is not a pandas Series

        Examples
        --------
        >>> import pandas as pd
        >>> s = pd.Series({'length_mean': 25.0, 'g': 1.02})
        >>> params = InvParameters.from_series(s)
        >>> params.values
        {'length_mean': 25.0, 'g': 1.02}

        Notes
        -----
        Parameters created this way have default bounds (no constraints)
        and vary=False, making them suitable for fixed-parameter scenarios.
        """

        if not isinstance(series, pd.Series):
            raise TypeError(
                f"InvParameters expects a `Series` for `.from_series`. Got type " f"{type(series)}."
            )
        return cls({k: {"value": v} for k, v in series.to_dict().items()})

    def simulate_parameter_sets(
        self,
        mc_realizations: Optional[int] = None,
        rng: Optional[np.random.Generator] = None,
    ) -> None:
        """
        Generate Monte Carlo sampled parameter sets used for initializing the optimization methods
        that inform the acoustic inversion and for uncertainty quantification.

        This class extends `InvParameters` to generate a list of realizations through Monte Carlo
        sampling. It is used for both uncertainty propagation and initializing the optimization
        algorithm in acoustic inversion by providing ensembles of parameter sets.
        """

        # Extract the initial parameter set
        params = self.parameters

        # Create random number generator if one is not defined
        if rng is None:
            rng = np.random.default_rng()

        # Generate parameter sets
        # ---- If not designated
        if not mc_realizations:
            self._realizations = {0: self.parameters}
        # ---- Generate sets
        else:
            self._realizations = {
                i: {
                    key: {
                        "value": (
                            rng.uniform(params[key]["min"], params[key]["max"])
                            if params[key].get("vary", True)
                            else params[key]["value"]
                        ),
                        "min": params[key]["min"],
                        "max": params[key]["max"],
                        "vary": params[key].get("vary", True),
                    }
                    for key in params
                }
                for i in range(mc_realizations)
            }
