from typing import Any, Dict, Optional

import numpy as np
import pandas as pd
from lmfit import Parameters
from pydantic import ConfigDict, Field, RootModel, ValidationError, model_validator

from ..core.validators import BaseDictionary


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
    scaled : bool, default=False
        Whether parameters are currently scaled to [0,1] range

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
    unscale_dict(scaled_dict)
        Convert scaled parameter dictionary back to original scale
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
        scaled: bool = False,
    ):
        if not isinstance(parameters, dict):
            raise TypeError("Parameters must be a dictionary.")
        for k, v in parameters.items():
            if not isinstance(k, (int, float, str)):
                raise TypeError(f"Parameter key '{k}' must be a numeric or string.")

        # Validate parameters
        validated_parameters = ModelInputParameters.create(**parameters)

        self.parameters = validated_parameters
        self.parameter_bounds = self._get_parameter_limits(validated_parameters)
        self.scaled = scaled

    def __str__(self):
        return "InvParameters-object"

    def __repr__(self):
        status = "Scaled" if self.scaled else "Unscaled"
        return (
            f"Status: {status}\n"
            f"Parameters\n----------\n[{', '.join(self.parameters)}]\n\n"
            f"Parameterizations\n-----------------\n{self.parameters}"
        )

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
        return self.scaled

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

    def scale(self):
        """
        Scale all parameters to [0,1] range using min-max normalization.

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

        # Update scaling status
        self.scaled = True

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
            for key, kwargs in self.parameters.items()
        }

        # Update the parameters
        self.parameters = scaled_params

    def unscale(self):
        """
        Convert scaled [0,1] parameters back to their original units and ranges.

        This method reverses the scaling transformation, restoring parameters
        to their original physical units and bounds. It uses the preserved
        bounds from `parameter_bounds` to perform the inverse transformation.

        Notes
        -----
        The unscaling transformation is:

        .. math::
            x_{unscaled} = x_{scaled} \\cdot (x_{max} - x_{min}) + x_{min}

        After unscaling, parameters return to their original bounds and units.
        The scaling status flag is updated to False after this operation.

        Raises
        ------
        ValueError
            If called on parameters that are not currently scaled
        """

        # Update scaling status
        self.scaled = False

        # Undo the min-max normalization
        unscaled_params = {
            key: {
                **kwargs,
                "value": (
                    kwargs["value"]
                    * (self.parameter_bounds[key]["max"] - self.parameter_bounds[key]["min"])
                    + self.parameter_bounds[key]["min"]
                ),
                "min": self.parameter_bounds[key]["min"],
                "max": self.parameter_bounds[key]["max"],
            }
            for key, kwargs in self.parameters.items()
        }

        # Update the parameters
        self.parameters = unscaled_params

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

    def unscale_dict(self, scaled_dict: dict) -> dict:
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
        >>> unscaled = params.unscale_dict(scaled_vals)
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


class MCInvParameters(InvParameters):
    """
    Monte Carlo sampled parameter sets for uncertainty quantification.

    This class extends InvParameters to generate multiple realizations of
    parameter values through Monte Carlo sampling. It's used for uncertainty
    propagation in acoustic inversion by providing ensembles of parameter
    sets that can be evaluated in parallel.

    Parameters
    ----------
    invparameters : InvParameters
        Base parameter set containing bounds and vary flags for sampling
    mc_realizations : int, optional
        Number of Monte Carlo realizations to generate. If None or 0,
        only the base parameter set is stored.
    rng : np.random.Generator, optional
        Random number generator for reproducible sampling. If None,
        uses default numpy random generator.

    Attributes
    ----------
    samples : dict
        Dictionary mapping realization indices to InvParameters instances,
        each containing a different random parameter realization
    parameter_bounds : dict
        Original parameter bounds inherited from base InvParameters

    Methods
    -------
    realizations
        Property to access the samples dictionary
    to_lmfit_samples(serialize=False)
        Convert all realizations to lmfit.Parameters objects
    __getitem__(idx)
        Access individual realizations by index
    __repr__()
        Detailed string representation showing sampling info
    __str__()
        Brief string identifier

    Examples
    --------
    >>> base_params = InvParameters({
    ...     'length_mean': {'value': 25.0, 'min': 20.0, 'max': 30.0, 'vary': True}
    ... })
    >>> mc_params = MCInvParameters(base_params, mc_realizations=100)
    >>> len(mc_params.samples)
    100
    >>> lmfit_ensemble = mc_params.to_lmfit_samples()

    Notes
    -----
    Parameter values are sampled uniformly within their specified bounds
    for parameters where vary=True. Fixed parameters (vary=False) retain
    their original values across all realizations.

    The sampling strategy supports warm-start optimization where the
    best-performing realization can be selected as an initialization point.
    """

    def __init__(
        self,
        invparameters: InvParameters,
        mc_realizations: Optional[int] = None,
        rng: Optional[np.random.Generator] = None,
    ):

        super().__init__(invparameters.parameters, scaled=invparameters.scaled)
        if rng is None:
            rng = np.random.default_rng()

        params = invparameters.parameters
        self.parameter_bounds = invparameters.parameter_bounds

        if not mc_realizations:
            self.samples = {0: params}
        else:
            self.samples = {
                i: InvParameters(
                    {
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
                    },
                    scaled=invparameters.scaled,
                )
                for i in range(mc_realizations)
            }
        # Update parameter_bounds for each sample
        for sample in self.samples.values():
            sample.parameter_bounds = self.parameter_bounds

    @property
    def realizations(self):
        """
        Access the dictionary of parameter realizations.

        Returns
        -------
        dict
            Dictionary mapping realization indices (int) to InvParameters
            instances, each containing a different parameter realization.

        Examples
        --------
        >>> mc_params = MCInvParameters(base_params, mc_realizations=5)
        >>> realizations = mc_params.realizations
        >>> first_real = realizations[0]  # First realization
        """
        return self.samples

    def to_lmfit_samples(self, serialize: bool = False):
        """
        Convert all parameter realizations to lmfit.Parameters format.

        This method transforms the Monte Carlo parameter ensemble into
        formats suitable for optimization routines, with options for
        both parallel and serialized parameter handling.

        Parameters
        ----------
        serialize : bool, default False
            If False, returns dict of separate Parameters objects (one per realization).
            If True, returns single Parameters object with all realizations as
            separate parameters (useful for concurrent optimization).

        Returns
        -------
        dict or lmfit.Parameters
            If serialize=False: Dictionary mapping realization indices to
            lmfit.Parameters objects
            If serialize=True: Single Parameters object with parameters named
            as {param_name}_{realization_index}

        Examples
        --------
        >>> mc_params = MCInvParameters(base_params, mc_realizations=3)
        >>> # Separate parameters for parallel evaluation
        >>> param_dict = mc_params.to_lmfit_samples(serialize=False)
        >>> # Single combined parameters for batch optimization
        >>> param_batch = mc_params.to_lmfit_samples(serialize=True)

        Notes
        -----
        The serialized format is useful when evaluating multiple realizations
        simultaneously in vectorized forward models or when using optimization
        algorithms that can handle large parameter spaces efficiently.
        """
        # Create concurrent, serialized set of parameters
        if serialize:
            # ---- Initialize
            params = Parameters()
            # ---- Iterate through each realization
            for i, sample in self.samples.items():
                # ---- Parse the parameters
                for pname, pinfo in sample.parameters.items():
                    # ---- Add them to the Parameters object
                    params.add(
                        f"{pname}_{i + 1}",
                        value=pinfo["value"],
                        min=pinfo["min"],
                        max=pinfo["max"],
                        vary=pinfo.get("vary", True),
                    )
            return params
        else:
            return {i: sample.to_lmfit() for i, sample in self.samples.items()}

    def __getitem__(self, idx):
        """
        Access a specific parameter realization by index.

        Parameters
        ----------
        idx : int
            Index of the realization to retrieve

        Returns
        -------
        InvParameters
            Parameter realization at the specified index

        Raises
        ------
        KeyError
            If the index is not in the samples dictionary
        """
        return self.samples[idx]

    def __repr__(self):
        """
        Concise string representation showing number of samples.

        Returns
        -------
        str
            Brief description of the MCInvParameters instance
        """
        if len(self.samples) == 1:
            return "MCInvParameters(1 sample)"
        else:
            return f"MCInvParameters({len(self.samples)} samples)"

    def __str__(self):
        """
        Detailed string representation showing all sample values.

        Returns
        -------
        str
            Multi-line string showing parameter values for each realization
        """
        return f"MCInvParameters with {len(self.samples)} samples:\n" + "\n".join(
            [f"Sample {i}: {sample.values}" for i, sample in self.samples.items()]
        )


class SingleParameter(
    BaseDictionary,
    arbitrary_types_allowed=True,
    title="values for lmfit.Parameters class required for optimization",
):
    min: Optional[float] = Field(default=-np.inf, allow_inf_nan=True)
    value: float = Field(allow_inf_nan=False)
    max: Optional[float] = Field(default=np.inf, allow_inf_nan=True)
    vary: bool = Field(default=False)

    @model_validator(mode="after")
    def check_bounds(self):
        if self.min > self.max:
            raise ValueError(f"min ({self.min}) cannot be greater than max ({self.max}).")
        if not (self.min <= self.value <= self.max):
            raise ValueError(
                f"value ({self.value}) must be between min ({self.min}) and max ({self.max})."
            )
        return self


class ModelInputParameters(
    RootModel[Dict[str, "SingleParameter"]],
):

    model_config = ConfigDict(
        arbitrary_types_allowed=True,
        title="scattering model parameters",
    )

    @model_validator(mode="before")
    @classmethod
    def prevalidator_trans(cls, data):
        if isinstance(data, InvParameters):
            return data.parameters
        return data

    # Validator method
    @classmethod
    def judge(cls, **kwargs):
        """
        Validator method
        """
        try:
            return cls(**kwargs)
        except ValidationError as e:
            raise e

    # Factory method
    @classmethod
    def create(cls, **kwargs):
        """
        Factory creation method
        """
        return cls.judge(**kwargs).model_dump(exclude_none=True)
