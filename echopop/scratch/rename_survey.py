import numpy as np
import time
def ragged_rowwise_apply(arr, fn):
    """Apply function to each row of ragged array efficiently"""
    return awk.Array([fn(row) for row in awk.to_list(arr)])

# Interpolate the Bessel function output over all orientations
def _fast_bessel_interp(ARG_flat, ka_norm_np, J1_np, n_theta):
    """Pure numpy interpolation - no compilation overhead"""
    result = np.zeros((ARG_flat.shape[0], n_theta))
    
    # Use numpy's vectorized operations instead of explicit loops
    for th in range(n_theta):
        result[:, th] = np.interp(ARG_flat[:, th], ka_norm_np, J1_np)
    
    return result

def uniformly_bent_cylinder_optimized(
    n_segments_array: awk.Array,
    radius_of_curvature_ratio: float,
    taper_order: float,
) -> Tuple[awk.Array, awk.Array, awk.Array, awk.Array, awk.Array]:
    """
    Optimized version of uniformly_bent_cylinder with better performance for large arrays
    """
    # Pre-compute constants
    gamma = 0.5 / radius_of_curvature_ratio
    norm_ratio = radius_of_curvature_ratio * 2
    eps = np.finfo(float).eps
    
    # Process each segment count efficiently
    results = []
    for n_segments in n_segments_array:
        # Create z coordinates
        z = np.linspace(-1.0, 1.0, n_segments)
        
        # Compute taper
        taper_i = np.sqrt(1 - z**taper_order)
        
        # Bend the cylinder
        z_curved = np.sin(gamma) * z
        x_curved = 1 - np.sqrt(1 - z_curved**2)
        
        # Normalize curvature
        z_norm = z_curved * norm_ratio
        x_norm = x_curved * norm_ratio
        
        # Calculate slope
        gamma_tilt_i = np.arctan2(z_norm, x_norm)
        
        # Calculate differences efficiently
        dz = np.diff(z_norm)
        dx = np.diff(x_norm) + eps
        
        # Calculate tilt angles
        alpha_tilt = np.arctan(dz / dx)
        alpha_tilt_final = np.arctan(dz[-1] / dx[-1])
        alpha_tilts = np.concatenate([alpha_tilt, [alpha_tilt_final]])
        
        # Beta tilt angles
        beta_tilt_i = np.where(alpha_tilts >= 0.0, alpha_tilts - np.pi / 2, alpha_tilts + np.pi / 2)
        
        # Position vectors
        r_pos_i = np.sqrt(x_norm**2 + z_norm**2)
        
        # First derivatives
        dr_first = np.sqrt(dx[0]**2 + dz[0]**2)
        dr_rest = np.sqrt(dx**2 + dz**2)
        dr_pos_i = np.concatenate([[dr_first], dr_rest])
        
        results.append((taper_i, gamma_tilt_i, beta_tilt_i, r_pos_i, dr_pos_i))
    
    # Convert to awkward arrays
    taper = awk.Array([result[0] for result in results])
    gamma_tilt = awk.Array([result[1] for result in results])
    beta_tilt = awk.Array([result[2] for result in results])
    r_pos = awk.Array([result[3] for result in results])
    dr_pos = awk.Array([result[4] for result in results])
    
    return taper, gamma_tilt, beta_tilt, r_pos, dr_pos

def fast_einsum_reduction(arr_3d_ragged, weights_ragged):
    """Fast Einstein summation for ragged 3D arrays"""
    return ragged_rowwise_apply(
        awk.zip(arr_3d_ragged, weights_ragged),
        lambda args: np.einsum('ijk,j->ik', args[0], args[1]) + np.finfo(float).eps
    )

def batch_bessel_j1(ka_norm_ragged):
    """Batch compute Bessel functions for ragged arrays"""
    return ragged_rowwise_apply(ka_norm_ragged, j1)

def vectorized_bessel_interp(ARG_ragged, ka_norm_ragged, J1_ragged, n_theta, n_k):
    """Highly optimized Bessel interpolation for ragged arrays"""
    def _fast_bessel_interp(ARG_flat, ka_norm_np, J1_np, n_theta):
        """Pre-allocated array interpolation - much faster than column_stack"""
        result = np.zeros((ARG_flat.shape[0], n_theta))
        
        # Use numpy's vectorized operations with pre-allocated result
        for th in range(n_theta):
            result[:, th] = np.interp(ARG_flat[:, th], ka_norm_np, J1_np)
        
        return result
    
    return awk.Array([
        _fast_bessel_interp(
            np.array(ARG_ragged[i]).ravel(order="F").reshape(
                (n_k[i] * len(ARG_ragged[i][0]), n_theta), order="F"
            ),
            np.array(ka_norm_ragged[i]), 
            np.array(J1_ragged[i]), 
            n_theta
        )
        for i in range(len(ARG_ragged))
    ])
# ==============================================================================
# OPTIMIZED EXTERNAL PCDWBA FUNCTION
# ==============================================================================
import cProfile
import pstats
from io import StringIO

def run_copilot_1(
    center_frequencies: np.ndarray[float],
    length_mean: float,
    length_sd_norm: float,
    length_radius_ratio: float,
    taper_order: float,
    radius_of_curvature_ratio: float,
    theta_mean: float,  
    theta_sd: float,
    orientation_distribution: Dict[str, Any],
    g: float,
    h: float,
    sound_speed_sw: float,
    frequency_interval: float,
    n_integration: int,
    n_wavelength: int, 
    number_density: Optional[float] = None,
    length_distribution: Optional[Dict[str, Any]] = None,
    **kwargs
):
    t0 = time.time()
    # Pre-allocate arrays based on input size
    n_freqs = len(center_frequencies)
    n_theta = orientation_distribution["bins"]
    n_length = length_distribution["bins"]
    
    # Generate frequency intervals centered on the central frequencies
    frequencies = generate_frequency_interval(
        center_frequencies,
        length_sd_norm,
        frequency_interval
    )
    t1 = time.time()
    print(f"[run_copilot_1] Frequency generation: {t1-t0:.4f}s")
    # Compute the acoustic wavenumbers weighted by target size
    k_c = wavenumber(center_frequencies, sound_speed_sw)
    ka_c = k_c * length_mean / length_radius_ratio
    k_f = wavenumber(frequencies, sound_speed_sw)
    ka_f = k_f * length_mean / length_radius_ratio
    t2 = time.time()
    print(f"[run_copilot_1] Wavenumber computation: {t2-t1:.4f}s")
    # Compute over a vector of angles (centered on 90 degrees)
    theta_values = np.linspace(
        theta_mean - 3.1 * theta_sd,
        theta_mean + 3.1 * theta_sd,
        n_theta,
    )
    theta_radians = theta_values * np.pi / 180.0
    t3 = time.time()
    print(f"[run_copilot_1] Angle vector computation: {t3-t2:.4f}s")
    # Calculate the appropriate number of integration points
    kL_max = np.nanmax(k_f * length_mean, axis=1) * (1 + 3.1 * length_sd_norm)
    n_int = awk.values_astype(
        np.where(
            kL_max < n_integration, 
            n_integration, 
            np.ceil(kL_max * n_wavelength / (2 * np.pi))
        ),
        int
    )
    t4 = time.time()
    print(f"[run_copilot_1] Integration points: {t4-t3:.4f}s")
    # Create shape/position matrix
    taper, gamma_tilt, beta_tilt, r_pos, dr_pos = uniformly_bent_cylinder(
        n_int, radius_of_curvature_ratio, taper_order
    )
    t5 = time.time()
    print(f"[run_copilot_1] Cylinder shape: {t5-t4:.4f}s")
    # Get the array sizes for later broadcasting
    # ---- Number of frequencies/wavenumbers
    n_k = awk.num(ka_f, axis=-1)
    # ---- Number of segments and minimum number of integration points
    n_segments = awk.num(r_pos, axis=-1)
    # ---- Number of orientation values
    n_theta = len(theta_radians)

    # Compute the reflection coefficient, `C_b`
    C_b = reflection_coefficient(g, h)

    # Adjust `ka` to account for body shape tapering
    ka_f_tapered = awk.Array([
        np.outer(np.array(ka_f[i, :n_k[i]]), taper[i, :n_segments[i]]) / h
        for i in range(len(n_k))
    ])

    # Adjust the along-axis tilt angles and slopes so they are relative to the incident planar 
    # wave
    # ---- Curvature slopes
    delta_gamma_cos = awk.Array([
        np.cos(np.array(gamma_tilt[i, :n_segments[i]]).reshape(-1, 1) - theta_radians)
        for i in range(len(n_k))
    ])
    # ---- Along-axis intersegment tilt angles
    delta_theta_cos = awk.Array([
        np.abs(np.cos(np.array(beta_tilt[i, :n_segments[i]]).reshape(-1, 1) - theta_radians))
        for i in range(len(n_k))
    ])

    # Calculate the phase term
    # ---- Compute the bulk of the exponentiated matrix
    M1 = awk.Array([
        length_radius_ratio * np.array(ka_f[i, :n_k[i]]).reshape(-1, 1) * 
        (np.array(r_pos[i, :n_segments[i]]) / h)
        for i in range(len(n_k))
    ])
    # ---- Compute phase
    phase = awk.Array([
        np.exp(1j * np.array(M1[i][:, :, np.newaxis]) * 
               np.array(delta_gamma_cos[i][np.newaxis, :, :]))
        for i in range(len(n_k))
    ])

    # Calculate the effect of the material properties with respect to the position matrix derivative
    M2 = awk.Array([
        h**2 * C_b * np.array(dr_pos[i, :n_segments[i]]) / 4
        for i in range(len(n_k))
    ])

    # Prepare the ka_f values for the Bessel function
    ka_f_last = awk.Array([np.array(ka_f[i, n_k[i]-1]) for i in range(len(n_k))])
    ka_f_norm = awk.Array([
        np.linspace(-2 * ka_f_last[i], 2 * ka_f_last[i], 2 * len(r_pos[i]))
        for i in range(len(n_k))
    ])
    J1_vals = awk.Array([j1(ka_f_norm[i]) for i in range(len(ka_f_norm))])
    t6 = time.time()
    print(f"[run_copilot_1] Bessel prep: {t6-t5:.4f}s")
    ARG = 2 * awk.Array([
        np.array(ka_f_tapered[i][:, :, np.newaxis]) * np.array(delta_theta_cos[i][np.newaxis, :, :])
        for i in range(len(n_k))
    ]) + eps
    J1_interp = vectorized_bessel_interp(ARG, ka_f_norm, J1_vals, orientation_distribution["bins"], n_k)
    t7 = time.time()
    print(f"[run_copilot_1] Bessel interpolation: {t7-t6:.4f}s")
    J1_norm = awk.Array([
        (
            np.array(J1_interp[i]) / 
            np.array(ARG[i]).ravel(order="F").reshape((n_k[i] * len(r_pos[i]), orientation_distribution["bins"]), order="F")
        ).reshape((n_k[i], len(r_pos[i]), orientation_distribution["bins"]), order="F")
        for i in range(len(n_k))
    ])
    t8 = time.time()
    print(f"[run_copilot_1] Bessel normalization: {t8-t7:.4f}s")
    f_j = ka_f_tapered ** 2 * J1_norm * phase
    f_bs = awk.Array([np.einsum("ijk, j->ik", f_j[i], M2[i]) for i in range(len(n_k))])
    t9 = time.time()
    print(f"[run_copilot_1] DWBA and einsum: {t9-t8:.4f}s")
    f_bs_orientation = awk.Array([
        np.sqrt(np.matmul((np.array(f_bs[i]).real ** 2 + 
                            np.array(f_bs[i]).imag ** 2), PDF))
        for i in range(len(f_bs))
    ])
    t10 = time.time()
    print(f"[run_copilot_1] Orientation averaging: {t10-t9:.4f}s")
    sigma_bs_length = awk.Array([
        np.sum(length_norm**2 * PDF_length * 
                np.interp(length_norm * ka_c[i], np.arange(len(f_bs_orientation[i])), 
                            np.array(f_bs_orientation[i]) ** 2))
        for i in range(len(f_bs_orientation))
    ])
    t11 = time.time()
    print(f"[run_copilot_1] Length integration: {t11-t10:.4f}s")
    Sv_prediction = 10 * np.log10(number_density * sigma_bs_length)
    ed = time.time()
    print(f"[run_copilot_1] Total: {ed-t0:.4f}s")
    return Sv_prediction

def run_copilot_2(
    center_frequencies: np.ndarray[float],
    length_mean: float,
    length_sd_norm: float,
    length_radius_ratio: float,
    taper_order: float,
    radius_of_curvature_ratio: float,
    theta_mean: float,  
    theta_sd: float,
    orientation_distribution: Dict[str, Any],
    g: float,
    h: float,
    sound_speed_sw: float,
    frequency_interval: float,
    n_integration: int,
    n_wavelength: int, 
    number_density: Optional[float] = None,
    length_distribution: Optional[Dict[str, Any]] = None,
    **kwargs):
    t0 = time.time()
    # Use pre-computed constants
    sqrt_2pi = np.sqrt(2 * np.pi)
    eps = np.finfo(float).eps

    # Pre-allocate arrays based on input size
    n_freqs = len(center_frequencies)
    n_theta = orientation_distribution["bins"]
    n_length = length_distribution["bins"]
    
    # Generate frequency intervals centered on the central frequencies
    frequencies = generate_frequency_interval(
        center_frequencies,
        length_sd_norm,
        frequency_interval
    )
    t1 = time.time()
    print(f"[run_copilot_2] Frequency generation: {t1-t0:.4f}s")
    # Compute the acoustic wavenumbers weighted by target size
    # ---- Center frequencies
    k_c = wavenumber(center_frequencies, sound_speed_sw)
    # ---- Compute ka (center frequencies)
    ka_c = k_c * length_mean / length_radius_ratio
    # ---- Frequency intervals
    # -------- Just wavenumber (`k`)
    k_f = wavenumber(frequencies, sound_speed_sw)
    # -------- Now `ka`
    ka_f = k_f * length_mean / length_radius_ratio
    t2 = time.time()
    print(f"[run_copilot_2] Wavenumber computation: {t2-t1:.4f}s")
    # Compute over a vector of angles (centered on 90 degrees)
    theta_values = np.linspace(
        theta_mean - 3.1 * theta_sd,
        theta_mean + 3.1 * theta_sd,
        n_theta,
    )
    # ---- Convert to radians
    theta_radians = theta_values * np.pi / 180.0
    t3 = time.time()
    print(f"[run_copilot_2] Angle vector computation: {t3-t2:.4f}s")
    # Calculate the appropriate number of integration points
    # ---- Compute threshold
    kL_max = np.nanmax(k_f * length_mean, axis=1) * (1 + 3.1 * length_sd_norm)
    # ---- Adjust number of integration points based on `kL_max`, if needed
    n_int = awk.values_astype(
        np.where(
            kL_max < n_integration, 
            n_integration, 
            np.ceil(kL_max * n_wavelength / (2 * np.pi))
        ),
        int
    )
    t4 = time.time()
    print(f"[run_copilot_2] Integration points: {t4-t3:.4f}s")
    # Create shape/position matrix
    taper, gamma_tilt, beta_tilt, r_pos, dr_pos = uniformly_bent_cylinder(
        n_int, radius_of_curvature_ratio, taper_order
    )
    t5 = time.time()
    print(f"[run_copilot_2] Cylinder shape: {t5-t4:.4f}s")
    # Get the array sizes for later broadcasting
    # ---- Number of frequencies/wavenumbers
    n_k = awk.num(ka_f, axis=-1)
    # ---- Number of segments and minimum number of integration points
    n_segments = awk.num(r_pos, axis=-1)
    # ---- Number of orientation values
    n_theta = len(theta_radians)

    # Compute the reflection coefficient, `C_b`
    C_b = reflection_coefficient(g, h)

    # Adjust `ka` to account for body shape tapering
    ka_f_tapered = awk.Array([
        np.outer(np.array(ka_f[i, :n_k[i]]), taper[i, :n_segments[i]]) / h
        for i in range(len(n_k))
    ])

    # Adjust the along-axis tilt angles and slopes so they are relative to the incident planar 
    # wave
    # ---- Curvature slopes
    delta_gamma_cos = awk.Array([
        np.cos(np.array(gamma_tilt[i, :n_segments[i]]).reshape(-1, 1) - theta_radians)
        for i in range(len(n_k))
    ])
    # ---- Along-axis intersegment tilt angles
    delta_theta_cos = awk.Array([
        np.abs(np.cos(np.array(beta_tilt[i, :n_segments[i]]).reshape(-1, 1) - theta_radians))
        for i in range(len(n_k))
    ])

    # Calculate the phase term
    # ---- Compute the bulk of the exponentiated matrix
    M1 = awk.Array([
        length_radius_ratio * np.array(ka_f[i, :n_k[i]]).reshape(-1, 1) * 
        (np.array(r_pos[i, :n_segments[i]]) / h)
        for i in range(len(n_k))
    ])
    # ---- Compute phase
    phase = awk.Array([
        np.exp(1j * np.array(M1[i][:, :, np.newaxis]) * 
               np.array(delta_gamma_cos[i][np.newaxis, :, :]))
        for i in range(len(n_k))
    ])

    # Calculate the effect of the material properties with respect to the position matrix derivative
    M2 = awk.Array([
        h**2 * C_b * np.array(dr_pos[i, :n_segments[i]]) / 4
        for i in range(len(n_k))
    ])

    # Prepare the ka_f values for the Bessel function
    ka_f_last = awk.Array([np.array(ka_f[i, n_k[i]-1]) for i in range(len(n_k))])
    ka_f_norm = awk.Array([
        np.linspace(-2 * ka_f_last[i], 2 * ka_f_last[i], 2 * n_segments[i])
        for i in range(len(n_k))
    ])
    # ---- Compute the cylindrical Bessel function(s) of the first kind
    J1 = awk.Array([j1(ka_f_norm[i]) for i in range(len(ka_f_norm))])

    # Calculate product of tapered ka_f and along-axis tilts - OPTIMIZED
    ARG = 2 * awk.Array([
        np.array(ka_f_tapered[i][:, :, np.newaxis]) * 
        np.array(delta_theta_cos[i][np.newaxis, :, :])
        for i in range(len(n_k))
    ]) + eps

    # Use the original working Bessel interpolation approach
    J1_interp = awk.Array([
        _original_bessel_interp(
            np.array(ARG[i]).ravel(order="F").reshape((n_k[i] * n_segments[i], n_theta), order="F"),
            np.array(ka_f_norm[i]), np.array(J1[i]), n_theta
        )
        for i in range(len(n_k))
    ])

    # Normalize the interpolated Bessel function - ORIGINAL METHOD
    J1_norm = awk.Array([
        (
            np.array(J1_interp[i]) / 
            np.array(ARG[i]).ravel(order="F").reshape((n_k[i] * n_segments[i], n_theta), order="F")
        ).reshape((n_k[i], n_segments[i], n_theta), order="F")
        for i in range(len(n_k))
    ])

    # Calculate the primary DWBA equation
    f_j = ka_f_tapered ** 2 * J1_norm * phase

    # Matrix multiplication via Einstein summation to get the final complex form function result
    f_bs = awk.Array([np.einsum("ijk, j->ik", f_j[i], M2[i]) for i in range(len(n_k))])
    
    # Orientation averaging
    # ---- Get interval
    orientation_interval = np.diff(theta_values).mean()
    # ---- Compute the PDF
    PDF = (
        orientation_interval
        * np.exp(-0.5 * (theta_values - theta_mean) ** 2 / theta_sd**2)
        / (sqrt_2pi * theta_sd)
    )

    f_bs_orientation = awk.Array([
        np.sqrt(np.matmul((np.array(f_bs[i]).real ** 2 + 
                           np.array(f_bs[i]).imag ** 2), PDF))
        for i in range(len(f_bs))
    ])

    # Compute over vector lengths
    length_values = np.linspace(
        length_mean - 3 * (length_sd_norm * length_mean),
        length_mean + 3 * (length_sd_norm * length_mean),
        n_length,
    )

    # Generate length deviation statistics
    length_deviation = length_mean * length_sd_norm

    # Normalize by mean length
    length_norm = length_values / length_mean
    length_sd_norm_calc = length_deviation / length_mean

    # ---- Get interval
    length_interval = np.diff(length_norm).mean()
    # ---- Compute the PDF
    PDF_length = (
        length_interval
        * np.exp(-0.5 * (length_norm - 1) ** 2 / length_sd_norm_calc**2)
        / (sqrt_2pi * length_sd_norm_calc)
    )

    # Integrate length vector over the form function
    sigma_bs_length = awk.Array([
        np.sum(length_norm**2 * PDF_length * 
               np.interp(length_norm * ka_c[i], np.array(ka_f[i]), 
                         np.array(f_bs_orientation[i]) ** 2))
        for i in range(len(f_bs_orientation))
    ])

    # Final length and normalization
    sigma_bs = sigma_bs_length * length_mean**2
    
    # Convert to Sv
    Sv_prediction = 10 * np.log10(number_density * sigma_bs)
    return Sv_prediction

def run_my_original(
    center_frequencies,
    length_mean,
    length_sd_norm,
    length_radius_ratio,
    taper_order,
    radius_of_curvature_ratio,
    theta_mean,
    theta_sd,
    orientation_distribution,
    g,
    h,
    sound_speed_sw,
    frequency_interval,
    n_integration,
    n_wavelength,
    number_density=None,
    length_distribution=None,
    **kwargs
):
    t0 = time.time()
    # Use pre-computed constants
    sqrt_2pi = np.sqrt(2 * np.pi)
    eps = np.finfo(float).eps

    # Pre-allocate arrays based on input size
    n_freqs = len(center_frequencies)
    n_theta = orientation_distribution["bins"]
    n_length = length_distribution["bins"]
    
    # Generate frequency intervals centered on the central frequencies
    frequencies = generate_frequency_interval(
        center_frequencies,
        length_sd_norm,
        frequency_interval
    )
    t1 = time.time()
    print(f"[run_my_original] Frequency generation: {t1-t0:.4f}s")
    # Compute the acoustic wavenumbers weighted by target size
    # ---- Center frequencies
    k_c = wavenumber(center_frequencies, sound_speed_sw)
    # ---- Compute ka (center frequencies)
    ka_c = k_c * length_mean / length_radius_ratio
    # ---- Frequency intervals
    # -------- Just wavenumber (`k`)
    k_f = wavenumber(frequencies, sound_speed_sw)
    # -------- Now `ka`
    ka_f = k_f * length_mean / length_radius_ratio
    t2 = time.time()
    print(f"[run_my_original] Wavenumber computation: {t2-t1:.4f}s")
    # Compute over a vector of angles (centered on 90 degrees)
    theta_values = np.linspace(
        theta_mean - 3.1 * theta_sd,
        theta_mean + 3.1 * theta_sd,
        n_theta,
    )
    # ---- Convert to radians
    theta_radians = theta_values * np.pi / 180.0
    t3 = time.time()
    print(f"[run_my_original] Angle vector computation: {t3-t2:.4f}s")
    # Calculate the appropriate number of integration points
    # ---- Compute threshold
    kL_max = np.nanmax(k_f * length_mean, axis=1) * (1 + 3.1 * length_sd_norm)
    # ---- Adjust number of integration points based on `kL_max`, if needed
    n_int = awk.values_astype(
        np.where(
            kL_max < n_integration, 
            n_integration, 
            np.ceil(kL_max * n_wavelength / (2 * np.pi))
        ),
        int
    )
    t4 = time.time()
    print(f"[run_my_original] Integration points: {t4-t3:.4f}s")
    # Create shape/position matrix
    taper, gamma_tilt, beta_tilt, r_pos, dr_pos = uniformly_bent_cylinder(
        n_int, radius_of_curvature_ratio, taper_order
    )
    t5 = time.time()
    print(f"[run_my_original] Cylinder shape: {t5-t4:.4f}s")
    # Get the array sizes for later broadcasting
    # ---- Number of frequencies/wavenumbers
    n_k = awk.num(ka_f, axis=-1)
    # ---- Number of segments and minimum number of integration points
    n_segments = awk.num(r_pos, axis=-1)
    # ---- Number of orientation values
    n_theta = len(theta_radians)

    # Compute the reflection coefficient, `C_b`
    C_b = reflection_coefficient(g, h)

    # Adjust `ka` to account for body shape tapering
    ka_f_tapered = awk.Array([
        np.outer(np.array(ka_f[i, :n_k[i]]), taper[i, :n_segments[i]]) / h
        for i in range(len(n_k))
    ])

    # Adjust the along-axis tilt angles and slopes so they are relative to the incident planar 
    # wave
    # ---- Curvature slopes
    delta_gamma_cos = awk.Array([
        np.cos(np.array(gamma_tilt[i, :n_segments[i]]).reshape(-1, 1) - theta_radians)
        for i in range(len(n_k))
    ])
    # ---- Along-axis intersegment tilt angles
    delta_theta_cos = awk.Array([
        np.abs(np.cos(np.array(beta_tilt[i, :n_segments[i]]).reshape(-1, 1) - theta_radians))
        for i in range(len(n_k))
    ])

    # Calculate the phase term
    # ---- Compute the bulk of the exponentiated matrix
    M1 = awk.Array([
        length_radius_ratio * np.array(ka_f[i, :n_k[i]]).reshape(-1, 1) * 
        (np.array(r_pos[i, :n_segments[i]]) / h)
        for i in range(len(n_k))
    ])
    # ---- Compute phase
    phase = awk.Array([
        np.exp(1j * np.array(M1[i][:, :, np.newaxis]) * 
               np.array(delta_gamma_cos[i][np.newaxis, :, :]))
        for i in range(len(n_k))
    ])

    # Calculate the effect of the material properties with respect to the position matrix derivative
    M2 = awk.Array([
        h**2 * C_b * np.array(dr_pos[i, :n_segments[i]]) / 4
        for i in range(len(n_k))
    ])

    # Prepare the ka_f values for the Bessel function
    ka_f_last = awk.Array([np.array(ka_f[i, n_k[i]-1]) for i in range(len(n_k))])
    ka_f_norm = awk.Array([
        np.linspace(-2 * ka_f_last[i], 2 * ka_f_last[i], 2 * n_segments[i])
        for i in range(len(n_k))
    ])
    J1_vals = awk.Array([j1(ka_f_norm[i]) for i in range(len(ka_f_norm))])
    t6 = time.time()
    print(f"[run_my_original] Bessel prep: {t6-t5:.4f}s")
    ARG = 2 * awk.Array([
        np.array(ka_f_tapered[i][:, :, np.newaxis]) * np.array(delta_theta_cos[i][np.newaxis, :, :])
        for i in range(len(n_k))
    ]) + eps

    # Interpolate the Bessel function output over all orientations
    def _fast_bessel_interp(ARG_flat, ka_norm_np, J1_np, n_theta):
        """Pure numpy interpolation - no compilation overhead"""
        result = np.zeros((ARG_flat.shape[0], n_theta))
        
        # Use numpy's vectorized operations instead of explicit loops
        for th in range(n_theta):
            result[:, th] = np.interp(ARG_flat[:, th], ka_norm_np, J1_np)
        
        return result
    
    J1_interp = awk.Array([
        _fast_bessel_interp(
            np.array(ARG[i]).ravel(order="F").reshape((n_k[i] * n_segments[i], n_theta), order="F"),
            np.array(ka_f_norm[i]), np.array(J1[i]), n_theta
        )
        for i in range(len(n_k))
    ])

    # Normalize the interpolated Bessel function
    t7 = time.time()
    print(f"[run_copilot_1] Bessel interpolation: {t7-t6:.4f}s")
    J1_norm = awk.Array([
        (
            np.array(J1_interp[i]) / 
            np.array(ARG[i]).ravel(order="F").reshape((n_k[i] * n_segments[i], n_theta), order="F")
        ).reshape((n_k[i], n_segments[i], n_theta), order="F")
        for i in range(len(n_k))
    ])
    t8 = time.time()
    print(f"[run_my_original] Bessel normalization: {t8-t7:.4f}s")
    f_j = ka_f_tapered ** 2 * J1_norm * phase
    f_bs = awk.Array([np.einsum("ijk, j->ik", f_j[i], M2[i]) for i in range(len(n_k))])
    t9 = time.time()
    print(f"[run_my_original] DWBA and einsum: {t9-t8:.4f}s")
    # ---- Get interval
    orientation_interval = np.diff(theta_values).mean()
    # ---- Compute the PDF
    PDF = (
        orientation_interval
        * np.exp(-0.5 * (theta_values - theta_mean) ** 2 / theta_sd**2)
        / (sqrt_2pi * theta_sd)
    )
    f_bs_orientation = awk.Array([
        np.sqrt(np.matmul((np.array(f_bs[i]).real ** 2 + 
                            np.array(f_bs[i]).imag ** 2), PDF))
        for i in range(len(f_bs))
    ])
    t10 = time.time()
    print(f"[run_my_original] Orientation averaging: {t10-t9:.4f}s")
    length_values = np.linspace(
        length_mean - 3 * (length_sd_norm * length_mean),
        length_mean + 3 * (length_sd_norm * length_mean),
        length_distribution["bins"],
    )


    length_deviation = length_mean * length_sd_norm

    # Normalize the length values, if needed
    length_norm = length_values / length_mean
    # ---- Also normalize the standard deviation
    length_sd_norm = length_deviation / length_mean

    # ---- Get the interval
    length_interval = np.diff(length_norm).mean()
    # ---- Compute the PDF
    PDF = (
        length_interval
        * np.exp(-0.5 * (length_norm - 1) ** 2 / length_sd_norm**2)
        / (sqrt_2pi * length_sd_norm)
    )

    # Vectorized computation - compute length-weighted ka for all frequencies
    ka_weighted = length_norm * ka_c.reshape(-1, 1)

    # Vectorized computation using awkward arrays
    sigma_bs_length = awk.Array([
        (
            length_norm**2
            * PDF
            * np.interp(ka_weighted[i], np.array(ka_f[i]), np.array(f_bs_orientation[i]) ** 2)
        ).sum()
        for i in range(len(f_bs_orientation))
    ])

    # ---- Convert to sigma_bs (linear backscattering cross-section)
    sigma_bs = sigma_bs_length * (length_mean) ** 2

    # Switch to logarithmic domain to compute S_V (volumetric backscattering strength)
    t11 = time.time()
    print(f"[run_my_original] Length integration: {t11-t10:.4f}s")
    Sv_prediction = 10 * np.log10(number_density * sigma_bs_length)
    ed = time.time()
    print(f"[run_my_original] Total: {ed-t0:.4f}s")
    return Sv_prediction

params = {
    k: v for k, v in locals().items()
    if k in [
        "center_frequencies", "length_mean", "length_sd_norm",
        "length_radius_ratio", "taper_order", "radius_of_curvature_ratio",
        "theta_mean", "theta_sd", "orientation_distribution",
        "g", "h", "sound_speed_sw", "frequency_interval",
        "n_integration", "n_wavelength", "number_density", "length_distribution"
    ]
}


run_my_original(**params)
run_copilot_2(**params)
