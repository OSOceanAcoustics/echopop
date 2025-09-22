import time
import warnings
from typing import Any, Dict, Literal, Optional

import numpy as np
import pandas as pd
from lmfit import Minimizer
from pydantic import ValidationError

from echopop.core.exceptions import EchopopValidationError
from echopop.inversion.inversion_base import InversionBase
from echopop.typing import InvParameters, MCInvParameters
from echopop.validators.inversion import (
    SCATTERING_MODEL_PARAMETERS,
    ValidateBuildModelArgs,
    ValidateInversionMatrix,
)


def echopop_optim_cb(params, iter, resid, *args, **kwargs):
    """
    Callback function for optimization iterations in echopop.

    This callback function is called during each iteration of the optimization
    process to provide progress updates. It prints parameter values and residuals
    at specific iteration intervals.

    Parameters
    ----------
    params : lmfit.Parameters
        Current parameter values during optimization
    iter : int
        Current iteration number
    resid : numpy.ndarray
        Residual array from current iteration
    *args : tuple
        Additional positional arguments (unused)
    **kwargs : dict
        Additional keyword arguments (unused)

    Notes
    -----
    Progress is printed at iterations 1, 10, and every 25th iteration thereafter.
    The callback displays residual magnitude in dB and key parameter values.
    """
    if iter in (1, 10) or iter % 25 == 0:
        res_str = f"{float(resid):.4g}"
        param_str = " | ".join(
            f"{name}: {par.value:.4g}" for name, par in params.items() if par.vary
        )
        print(f"Iter: {iter} | Q abs[pred. - meas.]: {res_str} dB | {param_str}")


def monte_carlo_initialize(
    parameters_lmfit,
    center_frequencies,
    Sv_measured,
    parameters_meta: MCInvParameters,
    model_settings,
    **kwargs,
):
    """
    Initialize Monte Carlo optimization by selecting the best-fit realization.

    This function evaluates multiple parameter realizations from Monte Carlo
    sampling and selects the realization that produces the smallest residual
    as the starting point for optimization.

    Parameters
    ----------
    parameters_lmfit : dict or list
        Dictionary or list of lmfit.Parameters objects for each realization
    center_frequencies : numpy.ndarray
        Array of acoustic frequencies in Hz
    Sv_measured : numpy.ndarray or pandas.Series
        Measured volume backscattering strength values in dB re 1 m^-1
    parameters_meta : MCInvParameters
        Metadata container for Monte Carlo parameter realizations
    model_settings : dict
        Model configuration settings including model type and function
    **kwargs : dict
        Additional keyword arguments passed to fit_Sv

    Returns
    -------
    lmfit.Parameters
        The parameter set from the realization with the smallest fit residual

    Notes
    -----
    This function implements a "warm start" strategy for optimization by
    pre-selecting the most promising parameter realization before running
    the full optimization procedure.
    """

    # Evaluate fits
    fits = [
        fit_Sv(p, Sv_measured, center_frequencies, parameters_meta, model_settings)
        for p in parameters_lmfit.values()
    ]

    # Now `fits` contains the fit (objective) value for each realization
    best_idx = np.argmin(fits)

    # Return the initialize fits
    return parameters_lmfit[best_idx]


def prepare_minimizer(
    Sv_measured: pd.Series,
    scattering_params,
    model_settings,
    simulation_settings,
    verbose: bool,
    **kwargs,
):
    """
    Prepare lmfit.Minimizer objects for Monte Carlo optimization.

    This function creates multiple Minimizer instances for Monte Carlo
    optimization, with each minimizer corresponding to a different
    parameter realization. The function handles parameter sampling,
    objective function setup, and callback configuration.

    Parameters
    ----------
    Sv_measured : pandas.Series
        Measured volume backscattering strength indexed by frequency
    scattering_params : InvParameters
        Scattering model parameters with bounds and initial values
    model_settings : dict
        Model configuration including model type and function reference
    simulation_settings : dict
        Simulation parameters including Monte Carlo realizations and RNG
    verbose : bool
        Whether to enable verbose output during optimization
    **kwargs : dict
        Additional keyword arguments

    Returns
    -------
    list
        List of configured lmfit.Minimizer objects for Monte Carlo optimization

    Notes
    -----
    The function generates Monte Carlo parameter realizations and creates
    a separate Minimizer for each realization. Each Minimizer is configured
    with the same objective function but different parameter starting points.

    The objective function minimizes the sum of squared differences between
    predicted and measured Sv values across all valid frequencies.
    """
    # Generate parameter sets
    parameter_sets = MCInvParameters(
        scattering_params, simulation_settings["mc_realizations"], simulation_settings["rng"]
    )

    # Generate `lmfit.Parameters` objects for each realization
    parameters_lmfit = parameter_sets.to_lmfit_samples()

    # Find which values thresholded out
    valid_idx = np.argwhere(Sv_measured > -999.0).flatten()

    # Get the valid center frequencies
    center_frequencies = np.array(Sv_measured.index.values, dtype=float)[valid_idx]

    # Turn into a numpy array
    Sv_measured = Sv_measured.to_numpy()[valid_idx].astype(float)

    # Initialize with Monte Carlo sampling
    if simulation_settings["monte_carlo"]:
        parameters_lmfit = monte_carlo_initialize(
            parameters_lmfit, center_frequencies, Sv_measured, parameter_sets, model_settings
        )
    else:
        parameters_lmfit = parameters_lmfit[0]

    # Set `iter_cb` if verbose
    if verbose:
        if "iter_cb" in simulation_settings and simulation_settings["iter_cb"]:
            iter_cb_arg = simulation_settings["iter_cb"]
        else:
            iter_cb_arg = echopop_optim_cb
    else:
        iter_cb_arg = None
    # Generate `Minimizer` function class required for bounded optimization
    return [
        Minimizer(
            fit_Sv,
            parameters_lmfit,
            fcn_args=(
                Sv_measured,
                center_frequencies,
                parameter_sets,
                model_settings,
            ),
            iter_cb=iter_cb_arg,
            nan_policy="omit",
        )
    ]


def fit_Sv(
    parameters, Sv_measured, center_frequencies, parameters_meta: InvParameters, model_settings
):

    # Extract the `lmfit.Parameters` values
    parameter_set = parameters.valuesdict()

    # Rescale to the original scale if parameters were normalized
    if parameters_meta.is_scaled:
        parameter_set = parameters_meta.unscale_dict(parameter_set)

    # Compute Sv
    Sv_prediction = model_settings["model_function"](
        center_frequencies=center_frequencies,
        **parameter_set,
        **model_settings["environment"],
        **model_settings,
    )

    return objective(Sv_prediction, Sv_measured)


def objective(Sv_prediction: np.ndarray[float], Sv_measured: np.ndarray[float]):
    r"""
    Compute the objective function for acoustic scattering model optimization.

    This function calculates the weighted sum of absolute deviations between
    predicted and measured volume backscattering strength (Sv) values. The
    objective function is minimized during parameter optimization.

    Parameters
    ----------
    Sv_prediction : numpy.ndarray
        Predicted Sv values from the scattering model in dB re 1 m^-1
    Sv_measured : numpy.ndarray
        Measured Sv values from acoustic data in dB re 1 m^-1

    Returns
    -------
    float
        Weighted sum of absolute deviations (objective function value)

    Notes
    -----
        The objective function is defined as:

    .. math::
        Q_d = \sum\limits_{i=1}^{M} \left\| \sigma_{\text{tot}}^i(\zeta) -
        \sum\limits{j=1}^{N} n_j \sigma_j^i(\zeta) \right\|_2^2

    where :math:`\zeta` is a feature vector containing scattering model parameters, and
    :math:`\sigma_\text{tot}^i` is the total differential backscattering cross-section for the
    :math:`i^\text{th}` frequency:

    .. math::
        \sigma_{\text{tot}}^i(\zeta) =
        \sum\limits_{j=1}^N n_j \sigma_j^i \sigma_j^i(\zeta),~i=1, 2, ...M

    Currently uses L1 norm (absolute deviation) rather than L2 norm (squared deviation) to be more
    robust to outliers in the acoustic measurements. This, effectively, is reduced to:

    .. math::
        Q = \sum_{i} w_i |Sv_{pred,i} - Sv_{meas,i}|

    where :math:`w_i` are frequency-specific weights.

    References
    ----------
    .. [1] Chu, D., Lawson, G.L., Wiebe, P.H. 2016. Estimation of biological parameters of marine
           organisms using lienar and nonlinear scattering model-based inversion methods. The
           Journal of the Acoustical Society of America, 139: 2885-2895. doi: 10.1121/1.4948759
    """

    # Calculate the deviation
    deviation = Sv_prediction - Sv_measured

    # Pre-allocate weight array
    wd = np.ones_like(Sv_measured)

    # Return the summed absolute deviation (Q)
    return np.sum(np.abs(deviation) * wd, axis=-1)


def perturb_parameters(params, scale=0.05):
    r"""
    Add random perturbations to optimization parameters for restart strategies.

    This function modifies parameter values by adding random noise scaled by
    the parameter magnitude. It respects parameter bounds and only perturbs
    parameters that are marked as variable in the optimization.

    Parameters
    ----------
    params : lmfit.Parameters
        Parameter set to be perturbed
    scale : float, optional
        Relative perturbation scale as a fraction of parameter value, by default 0.05

    Returns
    -------
    lmfit.Parameters
        Modified parameter set with random perturbations applied

    Notes
    -----
    Perturbations are applied as:

    .. math::
        p_{new} = p_{old} + \delta \cdot |p_{old}|

    where :math:`\delta` is drawn from :math:`U(-scale, scale)`.

    Parameters are clamped to their defined bounds after perturbation to
    ensure validity for subsequent optimization.
    """
    # params: lmfit.Parameters
    for name, par in params.items():
        if par.vary:
            # Calculate perturbation
            delta = np.random.uniform(-scale, scale) * abs(par.value if par.value != 0 else 1)
            new_value = par.value + delta
            # Respect bounds
            if par.min is not None:
                new_value = max(new_value, par.min)
            if par.max is not None:
                new_value = min(new_value, par.max)
            par.value = new_value
    return params


def optim(
    Sv_measured: pd.Series,
    scattering_parameters,
    simulation_settings,
    optimization_kwargs,
    verbose: bool = False,
):
    """
    Optimize scattering model parameters for a single measurement row.

    This function performs parameter optimization to minimize the difference
    between predicted and measured volume backscattering strength (Sv) values.
    It includes burn-in optimization, parameter perturbation for restarts,
    and comprehensive error handling.

    Parameters
    ----------
    Sv_measured : pandas.Series
        Measured Sv data with multiindex containing frequencies and metadata.
        Must have 'sv_mean', 'minimizer', and 'label' columns
    scattering_parameters : InvParameters
        Parameter definitions with bounds, initial values, and scaling info
    simulation_settings : dict
        Simulation configuration including minimum frequency count and flags
    optimization_kwargs : dict
        Optimization settings including method, tolerances, and burn-in params
    verbose : bool, optional
        Whether to print optimization progress and results, by default False

    Returns
    -------
    pandas.Series
        Optimized parameter values with residual error (Q) appended

    Notes
    -----
    The optimization process follows these steps:

    1. Filter valid (non-thresholded) frequency measurements
    2. Check minimum frequency count requirement
    3. Perform burn-in optimization if specified
    4. Run main optimization with filtered parameters
    5. Unscale parameters if needed
    6. Return results with error metrics

    The function uses L-BFGS-B or other methods from lmfit for optimization.
    Burn-in helps avoid local minima by pre-optimizing with relaxed constraints.

    Examples
    --------
    >>> results = optim(sv_row, params, sim_settings, opt_kwargs, verbose=True)
    >>> print(f"Final error: {results['Q']:.3f} dB")
    """
    # Catch start time in case `verbose=True`
    start_time = time.time()

    # Find which values are below the defined threshold
    valid_idx = np.argwhere(Sv_measured["sv_mean"] > -999.0).flatten()

    # Get the frequencies
    center_frequencies = np.array(Sv_measured["sv_mean"].index.values[valid_idx], dtype=float)

    # Only run if the correct number of frequencies are valid
    if len(valid_idx) >= simulation_settings["minimum_frequency_count"]:

        # Assign message string
        frequency_msg = ""

        with warnings.catch_warnings():
            warnings.filterwarnings(action="ignore", category=RuntimeWarning)
            # ---- Get the minimizer
            minimizer = Sv_measured["minimizer"].iloc[0][0]
            # ---- Sort kwargs
            lmfit_method_kws = {
                k: v
                for k, v in optimization_kwargs.items()
                if k not in ["burnin", "restart_strategy"]
            }
            # ---- Check for 'burn-in' parameterization
            if "burnin" in optimization_kwargs:
                # ---- Print
                if verbose:
                    print(
                        "Beginning burn-in optimization tuning.\n"
                        "--------------------------------------"
                    )
                # ---- Run coarse burn-in optimization
                parameters_coarse = minimizer.minimize(**optimization_kwargs["burnin"])
                # ---- Update the minimizer to run the full optimization configuration
                minimizer = Minimizer(
                    minimizer.userfcn,
                    parameters_coarse.params,
                    fcn_args=tuple(minimizer.userargs),
                    iter_cb=minimizer.iter_cb,
                    nan_policy="omit",
                )
            # ---- Print
            if verbose:
                print("Beginning parameter optimization.\n" "---------------------------------")
            # ---- Run optimization
            parameters_optimized = minimizer.minimize(**lmfit_method_kws)
            # ---- Extract
            best_fit_set_dict = parameters_optimized.params.valuesdict()
            # ---- Store fit deviation (Q)
            fit_error = parameters_optimized.residual[0]

            # Check for restart strategy
            restart_strategy = optimization_kwargs.get("restart_strategy", None)
            # ---- Run if defined
            if restart_strategy:
                if verbose:
                    print("Implementing restart strategy.\n" "------------------------------")
                # ---- Get number of restart attempts
                n_attempts = restart_strategy.get("max_attempts", 1)
                # ---- Perturb the parameters
                perturbation = restart_strategy.get("scale", 0.05)
                # ---- Get the Q threshold for breaking (assuming non-convergence)
                threshold = restart_strategy.get("Q_threshold", None)
                # ---- Iterate through the allotted attempts
                for attempt in np.arange(1, n_attempts + 1):
                    if threshold and fit_error < threshold:
                        break
                    if verbose:
                        print(f"Restart attempt: {attempt}.")
                    # ---- Update parameters
                    pert_params = perturb_parameters(
                        parameters_optimized.params, scale=perturbation
                    )
                    # ---- Update minimizer
                    minimizer = Minimizer(
                        minimizer.userfcn,
                        pert_params,
                        fcn_args=tuple(minimizer.userargs),
                        iter_cb=minimizer.iter_cb,
                        nan_policy="omit",
                    )
                    # ---- Restart optimizer
                    updated_params = minimizer.minimize(**lmfit_method_kws)
                    # ---- Check against previous `fit_error` and update parameters
                    if updated_params.residual[0] < fit_error:
                        best_fit_set_dict = updated_params.params.valuesdict()
                    fit_error = min(fit_error, updated_params.residual[0])

            # Unscale, if scaled
            if scattering_parameters.is_scaled:
                best_fit_set_dict = scattering_parameters.unscale_dict(best_fit_set_dict)

            # Set
            best_fit_set = pd.Series(best_fit_set_dict)

            # Add Q column
            best_fit_set["Q"] = fit_error
    else:
        # Assign message string
        frequency_msg = (
            f"\nWARNING: The number of frequencies with non-thresholded Sv [{len(valid_idx)}] "
            f"was fewer than the minimum frequency count "
            f"[{simulation_settings["minimum_frequency_count"]}]. Values were not optimized."
        )

        # Center frequencies, if any
        center_frequencies = center_frequencies

        # Create `pandas.Series`
        best_fit_set = pd.Series(scattering_parameters.values) * np.nan

        # Add error, Q, to series
        best_fit_set["Q"] = np.nan

    # Catch end time in case `verbose=True`
    end_time = time.time()

    # Print results
    if verbose:
        # ---- Row label
        row = f"{Sv_measured['label'].iloc[0]} | "
        # ---- Get error value
        error_value = f"Sv error (Q): {np.round(best_fit_set['Q'], 3)} dB "
        # ---- Parameter values
        # parameter_values = f"\n{best_fit_set[:-1].to_frame().T}"
        # ---- Get elapsed time (s)
        elapsed_time = f"\nElapsed time: {np.round(end_time - start_time, 2)} s;"
        # ---- Number of frequencies
        valid_freq = "[" + "/".join(f"{freq * 1e-3}" for freq in center_frequencies) + " Hz" + "]"

        # Print out
        print(
            f"===========\n{row}{error_value}{valid_freq}{frequency_msg}{elapsed_time}\n==========="
        )

    # Return
    return best_fit_set


class InversionMatrix(InversionBase):
    """
    Matrix-based acoustic scattering parameter inversion for marine organisms.

    This class performs acoustic inversion to estimate biological parameters
    (size, density, abundance) from multi-frequency volume backscattering
    strength measurements. It uses advanced scattering models and nonlinear
    optimization with optional Monte Carlo initialization for robust parameter
    estimation.

    Parameters
    ----------
    data : pd.DataFrame
        MultiIndex DataFrame containing acoustic measurements with columns:
        - 'sv_mean': Volume backscattering strength (dB re 1 m^-1)
        - 'nasc': Nautical Area Scattering Coefficient (m²/nmi²)
        - 'thickness_mean': Mean layer thickness (m)
        Frequency must be specified as a column index level.
    simulation_settings : Dict[str, Any]
        Configuration dictionary containing:
        - 'monte_carlo': bool, whether to use MC initialization
        - 'mc_realizations': int, number of MC samples
        - 'scale_parameters': bool, whether to scale parameters to [0,1]
        - 'environment': dict with 'sound_speed_sw' and 'density_sw'
        - 'reference_frequency': float, frequency for population estimates
        - 'minimum_frequency_count': int, minimum frequencies required
    verbose : bool, default=True
        Whether to print progress and diagnostic information

    Attributes
    ----------
    measurements : pd.DataFrame
        Validated acoustic measurement data with added processing columns
    simulation_settings : dict
        Validated and processed simulation configuration
    inversion_method : str
        Method identifier, set to "scattering_model"
    model_params : InvParameters
        Container for biological model parameters with bounds and vary flags
    model_settings : dict
        Scattering model configuration including type and numerical settings
    parameter_bounds : dict
        Original parameter bounds for unscaling operations
    rng : np.random.Generator
        Random number generator for reproducible Monte Carlo sampling

    Methods
    -------
    build_scattering_model(model_parameters, model_settings)
        Configure and validate scattering model parameters and settings
    invert(optimization_kwargs)
        Perform parameter inversion using nonlinear optimization

    Examples
    --------
    >>> # Initialize inversion matrix
    >>> inv_matrix = InversionMatrix(sv_data, simulation_settings)
    >>>
    >>> # Configure scattering model
    >>> params = InvParameters(parameter_dict)
    >>> model_config = {'type': 'pcdwba', 'taper_order': 10.0}
    >>> inv_matrix.build_scattering_model(params, model_config)
    >>>
    >>> # Run inversion
    >>> opt_kwargs = {'max_nfev': 1000, 'method': 'least_squares'}
    >>> results = inv_matrix.invert(opt_kwargs)
    >>>
    >>> # Extract population estimates
    >>> pop_estimates = estimate_population(results, nasc_data,
    ...                                   density_sw=1026.,
    ...                                   reference_frequency=120e3,
    ...                                   aggregate_method="transect")

    Notes
    -----
    The inversion process involves several steps:

    1. **Data Validation**: Ensures acoustic data has required structure
    2. **Model Configuration**: Sets up forward scattering model and parameters
    3. **Initialization**: Optionally uses Monte Carlo warm-start strategy
    4. **Optimization**: Minimizes misfit between predicted and measured Sv
    5. **Parameter Recovery**: Converts optimized values back to physical units

    Monte Carlo initialization can significantly improve convergence for
    nonlinear problems by evaluating multiple starting points and selecting
    the most promising realization.

    Parameter scaling normalizes all parameters to [0,1] range, which
    improves numerical conditioning when parameters have very different
    scales (e.g., length in mm vs density in kg/m³).

    The class supports various scattering models through a plugin architecture
    defined in `SCATTERING_MODEL_PARAMETERS`. Currently implemented models
    include PCDWBA for elongated organisms.

    Raises
    ------
    EchopopValidationError
        If data or simulation_settings fail validation requirements
    """

    def __new__(
        cls,
        data: pd.DataFrame,
        simulation_settings: Dict[str, Any],
        verbose: bool = True,
    ):
        # Validate
        try:
            # ---- Check
            valid_args = ValidateInversionMatrix.create(
                **dict(data=data, simulation_settings=simulation_settings)
            )
        # Break creation
        except ValidationError as e:
            raise EchopopValidationError(str(e)) from None

        # Create instance
        self = super().__new__(cls)

        # Update attributes
        self.measurements = valid_args["data"].copy()
        self.simulation_settings = valid_args["simulation_settings"]
        self.verbose = verbose

        # Generate
        return self

    def __init__(
        self,
        data: pd.DataFrame,
        simulation_settings: Dict[str, Any],
        verbose: bool = True,
    ):
        """
        Initialize the InversionMatrix with acoustic data and simulation settings.

        This constructor sets up the inversion framework by storing the input
        data and simulation configuration, initializing caches, and preparing
        the random number generator for Monte Carlo sampling.

        Parameters
        ----------
        data : pd.DataFrame
            MultiIndex DataFrame containing acoustic measurements with required
            columns 'sv_mean', 'nasc', and 'thickness_mean' at top level,
            and 'frequency' as nested index level.
        simulation_settings : Dict[str, Any]
            Dictionary containing simulation configuration including:
            - environment: Environmental parameters (sound speed, density)
            - monte_carlo: Whether to use Monte Carlo initialization
            - mc_realizations: Number of MC realizations (if monte_carlo=True)
            - scale_parameters: Whether to scale parameters to [0,1]
            - minimum_frequency_count: Min frequencies required for inversion
        verbose : bool, default=True
            Whether to print progress messages during initialization and inversion

        Notes
        -----
        The constructor performs minimal initialization, deferring heavy
        computation until build_scattering_model() is called. This allows
        for efficient object creation and parameter validation.

        The random number generator is configured based on simulation_settings
        to ensure reproducible Monte Carlo sampling when mc_seed is provided.
        """

        # Set inversion method
        self.inversion_method = "scattering_model"

        # Initialize attributes
        self._sv_cache = {}  # Initialize cache
        self.parameter_bounds = {}
        self.model = None
        self.model_params = {}
        self.model_settings = {}
        self.rng = None

        # Store random number generator, if required
        self._set_rng(**self.simulation_settings)

    def _set_rng(self, monte_carlo: bool, mc_seed: Optional[int] = None, **kwargs):
        """
        Configure random number generator for Monte Carlo sampling.

        Sets up a numpy random generator with optional seed for reproducible
        Monte Carlo parameter initialization. Only creates RNG when Monte
        Carlo sampling is enabled.

        Parameters
        ----------
        monte_carlo : bool
            Whether Monte Carlo sampling is enabled
        mc_seed : Optional[int], default=None
            Random seed for reproducible sampling. If None, uses system entropy.
        **kwargs : dict
            Additional simulation settings (ignored)

        Notes
        -----
        The RNG is stored in simulation_settings["rng"] for use by Monte Carlo
        initialization functions. This ensures consistent random sampling
        across all parameter realizations.
        """
        if monte_carlo:
            self.simulation_settings["rng"] = np.random.default_rng(mc_seed)

    def _set_minimizers(self):
        """
        Create and configure lmfit.Minimizer objects for each measurement.

        This method prepares optimization objects for each row of acoustic data
        by calling prepare_minimizer(). Each minimizer is configured with
        appropriate parameter bounds, Monte Carlo initialization, and model
        settings.

        Notes
        -----
        This is a computationally intensive step that:
        1. Creates parameter realizations for Monte Carlo initialization
        2. Evaluates initial fits to select best starting points
        3. Configures lmfit.Minimizer objects with proper bounds and callbacks

        The method includes extensive commented code for potential parallel
        processing, which could be enabled for large datasets.

        Progress messages are printed if verbose=True and Monte Carlo is enabled,
        as this step can take significant time for large datasets.

        Each minimizer is stored in the 'minimizer' column of self.measurements
        along with a descriptive 'label' for progress reporting during inversion.
        """

        # Define a new column for the `lmfit.Minimizer` class
        self.measurements["minimizer"] = np.array(np.nan).astype(object)

        # Create list of `Minimizer` objects depending on number of defined realizations

        _model_params = self.model_params
        _model_settings = self.model_settings
        _simulation_settings = self.simulation_settings
        _verbose = self.verbose
        if _verbose and _simulation_settings["monte_carlo"]:
            print(
                "Initializing parameter optimizers using Monte Carlo methods.\n"
                "------------------------------------------------------------"
            )
        self.measurements["minimizer"] = self.measurements["sv_mean"].apply(
            prepare_minimizer,
            axis=1,
            args=(
                _model_params,
                _model_settings,
                _simulation_settings,
                _verbose,
            ),
        )

        # Define label for verbosity
        self.measurements["label"] = [
            "; ".join(
                f"{name if name is not None else 'index'}: {val}"
                for name, val in zip(
                    self.measurements.index.names, (idx if isinstance(idx, tuple) else (idx,))
                )
            )
            for idx in self.measurements.index
        ]

    def build_scattering_model(
        self,
        model_parameters: InvParameters,
        model_settings: Dict[str, Any],
    ) -> None:
        """
        Configure the scattering model with parameters and computational settings.

        This method validates and stores the biological parameters and model
        configuration, then prepares the optimization framework by creating
        minimizer objects for each measurement location.

        Parameters
        ----------
        model_parameters : InvParameters
            Container with biological and physical parameters including:
            - Organism characteristics (length_mean, g, h, etc.)
            - Parameter bounds and vary flags for optimization
            - Current parameter values as initial estimates
        model_settings : Dict[str, Any]
            Model configuration dictionary containing:
            - type: Scattering model identifier (e.g., "pcdwba")
            - Model-specific settings (taper_order, frequency_interval, etc.)
            - Distribution parameters for length and orientation averaging
            - Environment parameters (sound_speed_sw, density_sw)

        Raises
        ------
        EchopopValidationError
            If parameters or settings fail validation against model schema

        Notes
        -----
        This method performs several important setup tasks:
        1. Validates parameter compatibility with the specified model type
        2. Scales parameters to [0,1] if scale_parameters=True in simulation_settings
        3. Loads the appropriate scattering model function from the registry
        4. Creates minimizer objects for each measurement (computationally intensive)

        Parameter scaling is recommended for optimization performance as it
        normalizes all parameters to similar ranges, improving numerical
        conditioning of the optimization problem.

        The model function is retrieved from SCATTERING_MODEL_PARAMETERS registry
        based on the model type, ensuring the correct forward model is used.

        Examples
        --------
        >>> inv_matrix = InversionMatrix(data, simulation_settings)
        >>> parameters = InvParameters(param_dict)
        >>> settings = {"type": "pcdwba", "taper_order": 10.0}
        >>> inv_matrix.build_scattering_model(parameters, settings)
        """

        # Validate
        try:
            # ---- Check
            valid_args = ValidateBuildModelArgs.create(
                **dict(model_parameters=model_parameters, model_settings=model_settings)
            )
        # Break creation
        except (ValidationError, Exception) as e:
            raise EchopopValidationError(str(e)) from None

        # Update attributes
        self.model_params = valid_args["model_parameters"]
        self.model_settings.update(valid_args["model_settings"])

        # Update the model settings to include the model function
        self.model_settings["model_function"] = SCATTERING_MODEL_PARAMETERS.get(
            self.model_settings["type"]
        )["function"]

        # Scale the parameters, if defined
        if self.simulation_settings["scale_parameters"] and not self.model_params.is_scaled:
            self.model_params.scale()

        # Add the `lmfit.Minimizer` objects to the dataset
        self._set_minimizers()

    def invert(self, optimization_kwargs):
        """
        Execute acoustic scattering parameter inversion for all measurements.

        This method runs the optimization process for each measurement location,
        estimating biological parameters by fitting the forward scattering model
        to observed volume backscattering strength data.

        Parameters
        ----------
        optimization_kwargs : dict
            Configuration for the optimization algorithm including:
            - `lmfit` keyword arguments
            - restart_strategy: Parameters for restart attempts on failure. This can include:
                - max_attempts: the mamximum number of restart attempts. The restart optimization
                loop will conclude once this number of attempts have been made.
                - Q_threshold: a threshold criterion for determining when to exit the restart
                optimization loop.
                - scale: a proportion (of the parameters) used to perturb the parameter fits upon
                restarting.
            - burnin: Pre-optimization with simpler algorithm. This can include any valid `lmfit`
            keyword argument.

        Returns
        -------
        pd.DataFrame
            DataFrame containing original acoustic measurements plus inversion
            results in 'parameters' column. Each entry is an InvParameters
            object with optimized parameter values and metadata.

        Notes
        -----
        The inversion process for each measurement involves:
        1. Forward model evaluation at current parameter values
        2. Residual calculation between predicted and measured Sv
        3. Parameter adjustment using the specified optimization algorithm
        4. Convergence checking based on tolerance criteria

        The method uses restart strategies to handle optimization failures:
        - Multiple attempts with perturbed initial conditions
        - Different optimization algorithms for difficult cases
        - Quality thresholds to detect poor fits

        Burn-in optimization with simpler algorithms (e.g., Nelder-Mead) can
        improve convergence for challenging parameter landscapes before
        switching to more sophisticated methods.

        Progress is reported during optimization if verbose=True, showing
        iteration numbers, residual values, and current parameter estimates.

        Examples
        --------
        >>> optimization_config = {
        ...     "method": "least_squares",
        ...     "max_nfev": 1000,
        ...     "loss": "huber",
        ...     "restart_strategy": {"max_attempts": 3, "scale": 0.1}
        ... }
        >>> results = inv_matrix.invert(optimization_config)
        >>> # Access optimized parameters for first measurement
        >>> optimal_params = results.iloc[0]["parameters"]
        """

        #
        inversion_results = self.measurements[["sv_mean", "minimizer", "label"]].apply(
            optim,
            axis=1,
            args=(
                self.model_params,
                self.simulation_settings,
                optimization_kwargs,
                self.verbose,
            ),
            result_type="expand",
        )

        # Prepare the output
        output = self.measurements.copy()

        # Drop the minimizer and label columns
        output.drop(["minimizer", "label"], axis=1, inplace=True)

        # Bundle the parameter results and add to the output
        output["parameters"] = [
            InvParameters.from_series(row.drop("Q")) for _, row in inversion_results.iterrows()
        ]

        return output


def estimate_population(
    inverted_data: pd.DataFrame,
    nasc_data: pd.DataFrame,
    density_sw: float,
    reference_frequency: float,
    aggregate_method: Literal["cells", "interval", "transect"],
):
    """
    Estimate population characteristics from inverted acoustic data.

    This function computes population-level estimates including body size,
    density, and biomass from inverted acoustic scattering parameters.
    It combines acoustic data with biological parameter estimates to
    provide comprehensive population assessments.

    Parameters
    ----------
    inverted_data : pd.DataFrame
        DataFrame containing inverted scattering parameters with columns
        including 'parameters' containing biological model parameters
    nasc_data : pd.DataFrame
        DataFrame containing Nautical Area Scattering Coefficient (NASC)
        values with frequency as column level
    density_sw : float
        Seawater density in kg/m³, typically ~1026 kg/m³
    reference_frequency : float
        Reference acoustic frequency in Hz for population estimation
    aggregate_method : Literal["cells", "interval", "transect"]
        Method for spatial aggregation of population estimates:
        - "cells": Individual acoustic cells
        - "interval": Depth/range intervals
        - "transect": Survey transect lines

    Returns
    -------
    pd.DataFrame
        Combined DataFrame containing:
        - nasc: NASC values at reference frequency
        - number_density: Estimated areal number density (individuals/m²)
        - biomass_density: Estimated areal biomass density (kg/m²)
        - Biological parameters: length_mean, radius_mean, body_volume,
          body_density, body_weight, etc.

    Notes
    -----
    The function assumes organisms can be modeled as uniformly bent cylinders
    for volume estimation. Body weight calculation uses:

    .. math::
        W = \\rho_{sw} \\cdot g \\cdot \\pi r^2 L

    where W is body weight, ρ_sw is seawater density, g is the relative
    density factor, r is mean radius, and L is mean length.

    References
    ----------
    .. [1] Stanton, T.K. et al. (1993). Acoustic scattering characteristics
           of several zooplankton groups. ICES Journal of Marine Science,
           50(3), 289-299.
    """

    # Unpackage the parameters using optimized extraction
    parameters = _extract_parameters_optimized(inverted_data)

    # Compute the average body radius for a single animal
    parameters["radius_mean"] = parameters["length_mean"] / parameters["length_radius_ratio"]

    # Estimate the average volume assuming an uniformly bent cylinder [single animal]
    parameters["body_volume"] = np.pi * parameters["radius_mean"] ** 2 * parameters["length_mean"]

    # Compute the animal body density [single animal]
    parameters["body_density"] = density_sw * parameters["g"]

    # Compute average body weight [single animal]
    parameters["body_weight"] = parameters["body_volume"] * parameters["body_density"]

    # Subset the dataset further to the single frequency
    acoustic_data = (
        inverted_data.loc[
            :, inverted_data.columns.get_level_values("frequency") == reference_frequency
        ]
    ).droplevel("frequency", axis=1)

    # Repeat the same for the coordinates
    reference_nasc = (
        nasc_data.loc[:, nasc_data.columns.get_level_values("frequency") == reference_frequency]
    )[reference_frequency].to_frame("nasc")

    # Calculate the areal number density
    if aggregate_method == "transect":
        reference_nasc.loc[:, "number_density"] = invert_transect_number_density(
            acoustic_data, reference_nasc, parameters
        ).to_numpy()
    else:
        reference_nasc["number_density"] = invert_intervals_number_density(
            acoustic_data, reference_nasc, parameters
        )

    # Compute areal biomass density
    reference_nasc["biomass_density"] = (
        reference_nasc["number_density"] * parameters["body_weight"]
    ).fillna(0.0)

    # Find shared indices
    idx = list(set(reference_nasc.index.names).intersection(parameters.index.names))

    # Reset indices
    reference_nasc = reference_nasc.reset_index().set_index(idx)

    # Concatenate the parameters
    return pd.concat(
        [reference_nasc, parameters.reset_index().set_index(idx).reindex(reference_nasc.index)],
        axis=1,
    )


def invert_transect_number_density(
    acoustic_data: pd.DataFrame,
    reference_nasc: pd.DataFrame,
    parameters: pd.DataFrame,
):
    """
    Compute areal number density aggregated by survey transect.

    This function estimates the areal number density of organisms by
    weighting volumetric densities by NASC contributions and aggregating
    across survey transects. The method accounts for varying layer
    thickness and spatial sampling intensity.

    Parameters
    ----------
    acoustic_data : pd.DataFrame
        DataFrame containing acoustic measurements including:
        - nasc: Nautical Area Scattering Coefficient values
        - thickness_mean: Mean layer thickness in meters
    reference_nasc : pd.DataFrame
        DataFrame with NASC values at reference frequency for weighting
    parameters : pd.DataFrame
        DataFrame containing biological parameters including:
        - number_density: Volumetric number density (individuals/m³)

    Returns
    -------
    pd.Series
        Areal number density in individuals/m² aggregated by transect,
        converted from nautical miles to meters (factor of 1852²)

    Notes
    -----
    The areal number density calculation follows:

    .. math::
        N_A = N_V \\cdot w \\cdot h \\cdot n \\cdot 1852^2

    where N_A is areal density, N_V is volumetric density, w is NASC weight,
    h is layer thickness, n is interval count, and 1852² converts from
    nautical miles² to meters².

    NASC weighting ensures proper spatial representation:

    .. math::
        w = \\frac{\\text{NASC}_{ref}}{\\text{NASC}_{acoustic}}
    """

    # Create copies
    acoustic_df = acoustic_data.copy()
    reference_df = reference_nasc.copy()

    # Find shared indices
    idx = list(set(reference_nasc.index.names).intersection(acoustic_data.index.names))

    # Reset the indices for each
    acoustic_df = acoustic_df.reset_index().set_index(idx)
    reference_df = reference_df.reset_index().set_index(idx)

    # Compute the NASC weights
    reference_df["nasc_weight"] = reference_df["nasc"] / acoustic_df["nasc"]

    # Get the counts of each transect
    interval_counts = reference_df.groupby(level="transect_num")["nasc"].count()

    # Compute the areal number density
    areal_number_density = (
        parameters["number_density"]
        * reference_df["nasc_weight"]
        * acoustic_df["thickness_mean"]
        * interval_counts
    ).fillna(0.0) * 1852**2

    # Return
    return areal_number_density


def invert_intervals_number_density(
    acoustic_data: pd.DataFrame,
    parameters: pd.DataFrame,
):
    """
    Compute areal number density for individual depth/range intervals.

    This function calculates areal number density by converting volumetric
    density estimates to areal density using layer thickness information.
    Unlike transect aggregation, this preserves interval-level resolution.

    Parameters
    ----------
    acoustic_data : pd.DataFrame
        DataFrame containing acoustic measurements with:
        - thickness_mean: Mean layer thickness in meters for each interval
    parameters : pd.DataFrame
        DataFrame containing biological parameters with:
        - number_density: Volumetric number density (individuals/m³)

    Returns
    -------
    pd.Series
        Areal number density in individuals/m² for each interval,
        converted from nautical miles to meters (factor of 1852²)

    Notes
    -----
    The conversion from volumetric to areal density is:

    .. math::
        N_A = N_V \\cdot h \\cdot 1852^2

    where N_A is areal density (individuals/m²), N_V is volumetric density
    (individuals/m³), h is layer thickness (m), and 1852² converts from
    nautical miles² to meters².

    This method maintains spatial resolution at the interval level, making
    it suitable for detailed vertical distribution analysis.
    """

    # Calculate the layer density
    layer_density = parameters["number_density"] * acoustic_data["thickness_mean"]

    # Return
    return layer_density * 1852**2


# OPTIMIZATION: More efficient parameter extraction
def _extract_parameters_optimized(inverted_data: pd.DataFrame) -> pd.DataFrame:
    """
    Extract parameters more efficiently than using .apply().

    This optimized version avoids the overhead of pandas .apply()
    by using direct iteration and batch DataFrame construction.

    Parameters
    ----------
    inverted_data : pd.DataFrame
        DataFrame with 'parameters' column containing InvParameters objects

    Returns
    -------
    pd.DataFrame
        DataFrame with parameter values as columns
    """
    # Extract all parameter dictionaries at once
    param_dicts = [obj.values for obj in inverted_data["parameters"]]

    # Construct DataFrame in one operation
    return pd.DataFrame(param_dicts, index=inverted_data.index)
