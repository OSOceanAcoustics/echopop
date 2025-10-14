def uniformly_bent_cylinder_old(
    n_segments: Union[int, np.ndarray[int]],
    radius_of_curvature_ratio: float,
    taper_order: float,
) -> pd.DataFrame:
    """
    Generates the normalized position matrix for an uniformly bent cylinder
    """

    st = time.time()
    # Curvature coefficients (for shape normalization)
    # ---- Gamma
    gamma = 0.5 / radius_of_curvature_ratio
    # ---- Normalization
    norm_ratio = radius_of_curvature_ratio * 2

    # Create normalized horizontal increments along the anterioposterior (z) axis of the body shape
    # ---- If array of values
    if isinstance(n_segments, np.ndarray):
        # ---- Maximum dimension
        max_n = n_segments.max()
        # ---- Generate the linearly spaced indices
        indices = np.arange(max_n)
        # ---- Broadcast a map for valid segment ranges
        valid_mask = indices < n_segments[:, None]
        # ---- Create a placeholder array padded with NaN
        z = np.full((n_segments.size, max_n), np.nan)
        # ---- Fill with the valid values
        z[valid_mask] = np.concatenate([np.linspace(-1.0, 1.0, n) for n in n_segments])
    # ---- If only a single value
    else:
        # ---- Maximum dimension
        max_n = n_segments
        # ---- Create vector
        z = np.linspace(-1.0, 1.0, n_segments)

    # Compute the taper vector
    taper = np.sqrt(1 - z**taper_order)

    # Bend the cylinder
    # ---- z-axis
    z_curved = np.sin(gamma) * z
    # ---- Dorsoventral axis (x-axis)
    x_curved = 1 - np.sqrt(1 - z_curved**2)

    # Normalize the curvature
    # ---- z-axis
    z_norm = z_curved * norm_ratio
    # ---- x-axis
    x_norm = x_curved * norm_ratio

    # Calculate the slope between curved segments
    gamma_tilt = np.arctan2(z_norm, x_norm)

    # Calculate the orientation angles between curved segments
    # ---- z-axis differences
    dz = np.diff(z_norm)
    # ---- x-axis differences
    dx = np.diff(x_norm) + np.finfo(float).eps
    # ---- alpha tilt angles
    alpha_tilt = np.arctan(dz / dx)
    # ---- Get the valid number of values per row
    n_valid = np.apply_along_axis(valid_array_row_length, 1, arr=alpha_tilt)
    # ---- Preallocate array
    new_column = np.array([])
    # ---- Generate values that are appended to the last valid column
    if alpha_tilt.ndim > 1:
        for i in range(alpha_tilt.shape[0]):
            if n_valid[i] == max_n - 1:
                new_column = np.append(new_column, np.arctan(dz[i, -1] / dx[i, -1]))
            else:
                new_column = np.append(new_column, np.nan)
                alpha_tilt[i, n_valid[i]] = np.arctan(dz[i, n_valid[i] - 1] / dx[i, n_valid[i] - 1])
        # ---- Now concatenate the missing column
        alpha_tilt = np.concatenate([alpha_tilt, new_column.reshape(-1, 1)], axis=1)
    else:
        alpha_tilt = np.append(alpha_tilt, np.arctan(dz[-1] / dx[-1]))
    # ---- beta tilt angles
    beta_tilt = np.where(alpha_tilt >= 0.0, alpha_tilt - np.pi / 2, alpha_tilt + np.pi / 2)

    # Compute the along-axis Euclidean distances to construct the position vector
    # ---- Position vector
    r_pos = np.sqrt(x_norm**2 + z_norm**2)
    # ---- Generate values that are appended to the first valid column
    if r_pos.ndim > 1:
        dr_pos = np.concatenate(
            [
                np.sqrt(dx[:, 0] * dx[:, 0] + dz[:, 0] * dz[:, 0]).reshape(-1, 1),
                np.sqrt(dx * dx + dz * dz),
            ],
            axis=1,
        )
    else:
        # ---- Compute the derivative of the position vector derivative
        dr_pos = np.append(np.sqrt(dx[0] * dx[0] + dz[0] * dz[0]), np.sqrt(dx * dx + dz * dz))
    ed = time.time()
    print(f"Cylinder generation time: {ed-st:5g}")
    # Return the relevant parameters
    return taper, gamma_tilt, beta_tilt, r_pos, dr_pos

def compute_ts_old(
    taper_order: float,
    length_sd_norm: float,
    length_mean: float,
    length_radius_ratio: float,
    radius_of_curvature_ratio: float,
    theta: Union[np.ndarray[float], float],
    k: Union[np.ndarray[float], float],
    ka: Union[np.ndarray[float], float],
    g: float,
    h: float,
    n_integration: int,
    ni_wavelen: int,
    model: Literal["pcdwba"] = "pcdwba",
) -> np.ndarray[complex]:
    """
    Compute the acoustic target strength (TS, dB re. 1 m^2) of a backscattering object.
    """
    # NOTE: !!! THE IS SPECIFIC TO THE PCDWBA. THESE LINES WILL NOT BE APPLICABLE TO OTHER MODELS.
    # THIS THEREFORE WILL REQUIRE REFACTORING (e.g. SCATTERING MODEL FUNCTION VALIDATORS)

    # Calculate the appropriate number of integration points
    # ---- Compute threshold
    kL_max = np.nanmax(k * length_mean, axis=1) * (1 + 3.1 * length_sd_norm)
    # ---- Adjust number of integration points based on `kL_max`, if needed
    n_int = np.where(
        kL_max < n_integration, n_integration, np.ceil(kL_max * ni_wavelen / (2 * np.pi))
    ).astype(int)

    # Create shape to build position vector and other required arrays
    taper, gamma_tilt, beta_tilt, r_pos, dr_pos = uniformly_bent_cylinder_old(
        n_int, radius_of_curvature_ratio, taper_order
    )

    # TS modeling
    if model == "pcdwba":
        f_bs_old = pcdwba_old(
            taper, gamma_tilt, beta_tilt, r_pos, dr_pos, length_radius_ratio, g, h, ka, theta
        )

    # Return output
    return f_bs_old


def compute_Sv_old(
    number_density: float,
    theta_values: np.ndarray[float],
    theta_mean: float,
    theta_sd: float,
    length_values: np.ndarray[float],
    length_mean: float,
    length_deviation: float,
    form_function: np.ndarray[complex],
    ka: np.ndarray[float],
    ka_center: np.ndarray[float],
) -> np.ndarray[float]:
    """
    Predict the volumetric backscattering strength from theoretical scattering estimates
    """

    # Orientation-averaged linear scattering coefficient
    f_bs_orientation = orientation_average_old(theta_values, form_function, theta_mean, theta_sd)

    # Length-averaged sigma_bs (normalized to length)
    sigma_bs_length = length_average_old(
        length_values, ka, ka_center, f_bs_orientation, length_mean, length_deviation
    )
    # ---- Convert to sigma_bs (linear backscattering cross-section)
    sigma_bs = sigma_bs_length * (length_mean) ** 2

    # Switch to logarithmic domain to compute S_V (volumetric backscattering strength)
    Sv_prediction = 10 * np.log10(number_density * sigma_bs)

    # Return predicted Sv
    return Sv_prediction

def generate_frequency_interval_old(
    frequency: np.ndarray[float], length_sd_norm: float, frequency_interval: float
) -> np.ndarray[float]:
    """
    Generate frequency interval 2D array centered on an array input of center frequencies.
    """

    frequency_lst = [
        np.arange(
            freq * (1 - 3.1 * length_sd_norm),
            freq * (1 + 3.1 * length_sd_norm) + frequency_interval,
            frequency_interval,
        )
        for freq in frequency
    ]

    # Find the maximum length of the generated arrays
    max_length = max(len(arr) for arr in frequency_lst)

    # Create a padded 2D array with NaN for shorter arrays
    padded_results = np.full((len(frequency_lst), max_length), np.nan)
    for i, arr in enumerate(frequency_lst):
        padded_results[i, : len(arr)] = arr

    return padded_results

def length_average_old(
    length: np.ndarray[float],
    ka: np.ndarray[float],
    ka_center: np.ndarray[float],
    form_function: np.ndarray[complex],
    length_mean: float,
    length_deviation: float,
    distribution: Literal["gaussian", "uniform"] = "gaussian",
) -> np.ndarray[float]:
    """
    """

    # Skip weighted averaging if only a single value is present
    # if len(form_function) == 1:
    #     # ---- If complex
    #     if np.all(np.iscomplex(form_function)):
    #         return np.sqrt(form_function.real**2 + form_function.imag**2)
    #     # ---- If not complex
    #     else:
    #         return form_function

    # Normalize the length values, if needed
    length_norm = length / length_mean
    # ---- Also normalize the standard deviation
    length_sd_norm = length_deviation / length_mean

    # Weight based on distribution input
    # ---- Gaussian (Normal)
    if distribution == "gaussian":
        # ---- Get the interval
        length_interval = np.diff(length_norm).mean()
        # ---- Compute the PDF
        PDF = (
            length_interval
            * np.exp(-0.5 * (length_norm - 1) ** 2 / length_sd_norm**2)
            / (np.sqrt(2 * np.pi) * length_sd_norm)
        )
    # ---- Uniform
    elif distribution == "uniform":
        # ---- Compute the PDF
        PDF = np.ones(len(form_function)) / len(form_function)
    else:
        raise ValueError("Invalid distribution type. Choose 'gaussian' or 'uniform'.")

    # Length-weighted averaged sigma_bs
    # ---- Get valid values
    n_vals = np.apply_along_axis(valid_array_row_length, 1, arr=ka)
    # ---- Compute the length-weighted ka
    ka_weighted = length_norm * ka_center.reshape(-1, 1)
    # ---- Trim values so they fall within the valid/defined bandwidth
    ka_weighted_trim = np.where(
        (ka_weighted >= np.nanmin(ka)) & (ka_weighted <= np.nanmax(ka)), ka_weighted, np.nan
    )
    # ---- Evaluate
    # -------- One frequency
    if len(form_function) == 1:
        sigma_bs_L = np.array(
            [
                (
                    length_norm**2
                    * PDF
                    * np.interp(ka_weighted_trim[0], ka[0], form_function[0] ** 2)
                ).sum()
            ]
        )
    # -------- More than one frequency
    else:
        sigma_bs_L = np.array(
            [
                (
                    length_norm**2
                    * PDF
                    * np.interp(ka_weighted_trim[i], ka[i, : n_vals[i]], form_function[i] ** 2)
                ).sum()
                for i in range(len(form_function))
            ]
        )
    # ---- Return the weighted average
    return sigma_bs_L


def inverse_normalize_series_old(series: pd.Series, ranges_dict: Dict[str, Any]):
    """
    Apply inverse min-max normalization to a `pandas.Series`
    """

    # Create copy
    series_copy = series.copy()

    # Get the columns that will be inverse transformed
    transform_cols = list(ranges_dict.keys())

    # Convert `ranges_dict` to a `pd.DataFrame`
    ranges_df = pd.DataFrame(ranges_dict).T

    # Apply inverse transformation
    series_copy[transform_cols] = (
        series_copy[transform_cols] * (ranges_df["high"] - ranges_df["low"]) + ranges_df["low"]
    )

    # Return the inverse min-max normalized series
    return series_copy


def orientation_average_old(
    angle: np.ndarray[float],
    form_function: np.ndarray[complex],
    theta_mean: float,
    theta_sd: float,
    distribution: Literal["gaussian", "uniform"] = "gaussian",
) -> np.ndarray[float]:
    """
    """

    # # Skip weighted averaging if only a single value is present
    # if len(form_function) == 1:
    #     # ---- If complex
    #     if np.all(np.iscomplex(form_function)):
    #         return [np.sqrt(f.real**2 + f.imag**2) for f in form_function[0]]
    #     # ---- If not complex
    #     else:
    #         return form_function

    # Weight based on distribution input
    # ---- Gaussian (Normal)
    if distribution == "gaussian":
        # ---- Get interval
        orientation_interval = np.diff(angle).mean()
        # ---- Compute the PDF
        PDF = (
            orientation_interval
            * np.exp(-0.5 * (angle - theta_mean) ** 2 / theta_sd**2)
            / (np.sqrt(2 * np.pi) * theta_sd)
        )
    # ---- Uniform
    elif distribution == "uniform":
        # ---- Compute the PDF
        PDF = np.ones(len(form_function)) / len(form_function)
    else:
        raise ValueError("Invalid distribution type. Choose 'gaussian' or 'uniform'.")

    # Return the weighted form function
    # ---- If complex
    return [
        (
            np.sqrt(np.matmul((f[0].real ** 2 + f[0].imag ** 2), PDF))
            if np.all(np.iscomplex(f))
            else np.sqrt(np.matmul(f, PDF))
        )
        for f in form_function
    ]

def valid_array_row_length(arr: np.ndarray[float]) -> int:
    """
    Returns the number of valid (i.e. not NaN) length of each row within an array
    """

    return np.sum(~np.isnan(arr))

def pcdwba_old(
    taper: float,
    gamma_tilt: np.ndarray[float],
    beta_tilt: np.ndarray[float],
    r_pos: np.ndarray[float],
    dr_pos: np.ndarray[float],
    length_radius_ratio: float,
    g: Union[np.ndarray[float], float],
    h: Union[np.ndarray[float], float],
    ka: Union[np.ndarray[float], float],
    theta: Union[np.ndarray[float], float],
) -> List[np.ndarray[complex]]:
    """
    """

    st = time.time()

    # Get the array sizes for later broadcasting
    # ---- Number of frequencies/wavenumbers
    n_k = (
        np.apply_along_axis(valid_array_row_length, axis=1, arr=ka)
        if isinstance(ka, np.ndarray)
        else len(ka)
    )
    # ---- Number of segments and minimum number of integration points
    n_segments = (
        np.apply_along_axis(valid_array_row_length, axis=1, arr=r_pos)
        if r_pos.ndim > 1
        else len(r_pos)
    )
    # ---- Number of orientation values
    n_theta = len(theta)

    # Compute the reflection coefficient, `C_b`
    C_b = reflection_coefficient(g, h)

    # Pre-allocate output (this will be indexed by center frequency)
    f_bs_old = []

    # Iterate across frequencies
    for i in range(len(n_k)):
        # ---- Reindex 'ka'
        ka_i = ka[i, : n_k[i]]
        # ---- Adjust `ka` to account for body shape tapering
        ka_tapered = ka_i.reshape(-1, 1) * taper[i, : n_segments[i]] / h
        # ka_tapered = ka.reshape(-1, 1) * taper / h
        # ---- Adjust along-axis tilt angles and slopes to be relative to the incident planar wave
        # -------- Along-axis curvature slopes
        delta_gamma_cos = np.cos(gamma_tilt[i, : n_segments[i]].reshape(-1, 1) - theta)
        # -------- Along-axis tilt angles between segments
        delta_theta_cos = np.abs(np.cos(beta_tilt[i, : n_segments[i]].reshape(-1, 1) - theta))
        # ---- Generate the exponentiated matrix
        M1 = length_radius_ratio * ka_i.reshape(-1, 1) * (r_pos[i, : n_segments[i]] / h)
        # ---- Generate the matrix that accounts for the material properties and position vector
        # ---- variability
        M2 = h**2 * C_b * dr_pos[i, : n_segments[i]] / 4
        # ---- Normalize the `ka` vector
        ka_norm = np.linspace(-2 * ka_i[-1], 2 * ka_i[-1], 2 * n_segments[i])
        # ---- Pre-compute the cylindrical Bessel function of the first kind of order 1
        J1 = j1(ka_norm)
        # ---- Broadcast `ka_tapered` for multiplication with `delta_gamma_cos`
        ARG = (
            2 * ka_tapered[:, :, np.newaxis] * delta_theta_cos[np.newaxis, :, :]
            + np.finfo(float).eps
        )
        # -------- Flatten for subsequent interpolation
        ARG_flat = ARG.ravel(order="F").reshape((n_k[i] * n_segments[i], n_theta), order="F")
        # ---- Interpolate the values
        J1_interp = np.array([np.interp(ARG_flat[:, m], ka_norm, J1) for m in range(n_theta)])
        # -------- Reshape and normalize
        J1_norm = (J1_interp.T / ARG_flat).reshape((n_k[i], n_segments[i], n_theta), order="F")
        # ---- Exponentiate the terms to compute the phase
        phase = np.exp(1j * M1[:, :, np.newaxis] * delta_gamma_cos[np.newaxis, :, :])
        # ---- Combine the adjusted `ka`, interpolated Bessel function output, and phase
        M3 = ka_tapered[:, :, np.newaxis] ** 2 * J1_norm * phase
        # ---- Compute the complex form function (matrix multiplication via Einstein summation)
        f_bs_old.append([np.einsum("ijk, j->ik", M3, M2) + np.finfo(float).eps])

    ed = time.time()
    print(f"OLD PCDWBA: {ed-st}")
    # Return the form function
    return f_bs_old


def normalize_parameters_old(
    parameter_sets: Dict[str, Any],
    inverse: bool = False,
    inverse_reference: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Normalize the optimization parameters

    Parameters
    ----------
    parameter_sets: Dict[str, Any]
        Dictionary comprising acoustic scattering model parameters that are to be optimized.
    inverse: bool, default: False
        Boolean flag to indicate whether the parameters should be normalized via min-max
        normalization (`inverse=True`) or denormalized (i.e. inverse normalization,
        `inverse=False`).
    """

    # Min-max normalization
    if inverse is False:
        return {
            realization: {
                key: (
                    {
                        **value,
                        "initial": (
                            (value["initial"] - value["low"]) / (value["high"] - value["low"])
                            if value["high"] != value["low"]
                            else 0.0
                        ),
                        "low": 0.0,
                        "high": 1.0,
                    }
                    if isinstance(value, dict)
                    else value
                )  # Skip scalar entries
                for key, value in sets.items()
            }
            for realization, sets in parameter_sets.items()
        }
    # Min-max inverse normalization
    else:
        return {
            key: (
                value * (inverse_reference[key]["high"] - inverse_reference[key]["low"])
                + inverse_reference[key]["low"]
            )
            for key, value in parameter_sets.items()
        }

def simulate_Sv_old(
    scattering_parameters: Parameters,
    config: Dict[str, Any],
) -> np.ndarray[float]:
    """
    Simulate the volumetric backscattering strength (Sv, dB re. m^-1)
    """

    times = pd.DataFrame({
        "milestone": ["t1", "t2", "t3", "t4", "t5", "t6"],
        "time": np.nan
    })


    st = time.time()
    
    # Extract parameter values from dictionary for parsing
    parameters_dict = scattering_parameters.valuesdict()

    # Rescale to the original scale
    if config["scale_parameters"]:
        parameters_dict = normalize_parameters_old(
            parameters_dict, inverse=True, inverse_reference=config["parameter_bounds"]
        )
    
    # Compute acoustic property metrics
    # ---------------------------------
    # Normalize the length standard deviation
    if "length_sd_norm" not in parameters_dict and "length_deviation" in parameters_dict:
        length_sd_norm = parameters_dict["length_deviation"] / parameters_dict["length_mean"]
    elif "length_sd_norm" not in parameters_dict and "length_sd_norm" in config:
        length_sd_norm = config["length_sd_norm"]
    else:
        length_sd_norm = parameters_dict["length_sd_norm"]

    to1 = time.time()
    times.loc[0, "time"] = to1 - st
    print(f"T1: Parameter handling: {(to1-st):.5g}")

    # Generate frequency intervals centered on the central frequencies
    frequencies = generate_frequency_interval_old(
        config["center_frequencies"],
        length_sd_norm,
        config["frequency_interval"],
    )

    to2 = time.time()
    times.loc[1, "time"] = to2 - st
    print(f"T2: Frequency generation: {(to2-to1):.5g}")

    # Compute the acoustic wavenumbers weighted by target size
    # ---- Center frequencies
    k_center = wavenumber(config["center_frequencies"], config["sound_speed_sw"])
    # ---- Compute ka (center frequencies)
    ka_center = k_center * parameters_dict["length_mean"] / parameters_dict["length_radius_ratio"]
    # ---- Frequency intervals
    # -------- Just wavenumber (`k`)
    k = wavenumber(frequencies, config["sound_speed_sw"])
    # -------- Now `ka`
    ka = k * parameters_dict["length_mean"] / parameters_dict["length_radius_ratio"]

    to3 = time.time()
    times.loc[2, "time"] = to3 - st
    print(f"T3: Wavenumbers: {(to3-to2):.5g}")

    # Compute over a vector of angles (centered on 90 degrees)
    theta_values = np.linspace(
        parameters_dict["theta_mean"] - 3.1 * parameters_dict["theta_sd"],
        parameters_dict["theta_mean"] + 3.1 * parameters_dict["theta_sd"],
        config["orientation_bin_count"],
    )
    theta_radians = theta_values * np.pi / 180.0

    # Compute over vector lengths
    length_values = np.linspace(
        parameters_dict["length_mean"] - 3 * (length_sd_norm * parameters_dict["length_mean"]),
        parameters_dict["length_mean"] + 3 * (length_sd_norm * parameters_dict["length_mean"]),
        config["length_bin_count"],
    )

    to4 = time.time()
    times.loc[3, "time"] = to4 - st
    print(f"T4: Distributions: {(to4-to3):.5g}")

    # PCDWBA (TS modeling step)
    # ------
    fbs = compute_ts_old(
        config["taper_order"],
        length_sd_norm,
        parameters_dict["length_mean"],
        parameters_dict["length_radius_ratio"],
        parameters_dict["radius_of_curvature_ratio"],
        theta_radians,
        k,
        ka,
        parameters_dict["g"],
        parameters_dict["h"],
        config["n_integration"],
        config["n_wavelength"],
        model=config["ts_model"]
    )
    to5 = time.time()
    times.loc[4, "time"] = to5 - st
    print(f"T5: Linear scattering coefficient: {(to5-to4):.5g}")

    # Compute S_V
    Sv_prediction = compute_Sv_old(
        parameters_dict["number_density"],
        theta_values,
        parameters_dict["theta_mean"],
        parameters_dict["theta_sd"],
        length_values,
        parameters_dict["length_mean"],
        parameters_dict["length_mean"] * length_sd_norm,
        fbs,
        ka,
        ka_center,
    )

    to6 = time.time()
    times.loc[5, "time"] = to6 - st
    print(f"T6: Sv: {(to6-to5):.5g}")

    # Return array
    return Sv_prediction, times