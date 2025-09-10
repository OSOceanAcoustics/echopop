import cProfile
import pstats
import os
import time
from echopop.inversion.inversion_matrix_krill import InversionMatrix
from echopop.typing import InvParameters
import echopop.nwfsc_feat.utils as utils
from pathlib import Path
from echopop.nwfsc_feat import ingest_sv

# ==================================================================================================
# ==================================================================================================
# DEFINE DATA ROOT DIRECTORY
# --------------------------
DATA_ROOT = Path("C:/Users/Brandyn/Documents/GitHub/EchoPro_data/echopop_inversion/")

# ==================================================================================================
# ==================================================================================================
# DATA INGESTION
# ==================================================================================================

sv_path = DATA_ROOT / "raw_files/2023/"
transect_pattern = r"x(\d+)"
impute_coordinates = True
center_frequencies = {18e3: {"min": -90., "max": -50.},
                      38e3: {"min": -90., "max": -50.},
                      70e3: {"min": -90., "max": -50.},
                      120e3: {"min": -90., "max": -50.},
                      200e3: {"min": -90., "max": -50.}}
method="transects"

sv_data = ingest_sv.ingest_echoview_sv(sv_path=sv_path, 
                                       center_frequencies=center_frequencies, 
                                       transect_pattern=transect_pattern, 
                                       aggregate_method=method,
                                       impute_coordinates=True)

# ==================================================================================================
# ==================================================================================================
# DATA SUBSET
# ==================================================================================================

sv_data_sub = utils.apply_filters(sv_data)

# ==================================================================================================
# ==================================================================================================
# DEFINE PARAMETERS
# ==================================================================================================

MODEL_PARAMETERS = {
    "number_density": {"value": 500., "min": 10., "max": 10000., "vary": True},
    "theta_mean": {"value": 10., "min": 0., "max": 90., "vary": True},
    "theta_sd": {"value": 20., "min": 0., "max": 90., "vary": False},
    "length_mean": {"value": 0.030, "min": 0.008, "max": 0.040, "vary": True},
    "length_sd_norm": {"value": 0.15, "min": 0.05, "max": 0.15, "vary": False},
    "g": {"value": 1.015, "min": 1.015, "max": 1.060, "vary": False},
    "h": {"value": 1.020, "min": 1.015, "max": 1.060, "vary":False},
    "radius_of_curvature_ratio": {"value": 3.0, "min": 0.5, "max": 100.0, "vary": False},
    "length_radius_ratio": {"value": 18.2, "min": 14.0, "max": 20.0, "vary": False},
}

# Model-specific settings, including distributions
MODEL_SETTINGS = {
    "type": "pcdwba",
    "taper_order": 10.,
    "frequency_interval": 2000.,
    "n_integration": 50,
    "n_wavelength": 10,
    "orientation_distribution": {"family": "gaussian", "bins": 60},
    "length_distribution": {"family": "gaussian", "bins": 100},
}

ENVIRONMENT = {
    "sound_speed_sw": 1500.0,
    "density_sw": 1026.9,    
}

SIMULATION_SETTINGS = {
    "environment": ENVIRONMENT,
    "monte_carlo": True,
    "mc_realizations": 2,
    "scale_parameters": True, 
    "minimum_frequency_count": 2,
    "reference_frequency": 120e3,
}

OPTIMIZATION_KWARGS = {
    "max_nfev": 1000,
    "method": "least_squares",
    "loss": "huber",
    "xtol": 1e-8,
    "ftol": 1e-8,
    "gtol": 1e-8,
    "diff_step": 1e-7,
    "verbose": 2,
    "burnin": {
        "method": "nelder",
        "max_nfev": 100,
        "tol": 1e-3,
    }
}

def run_workflow():
# ==================================================================================================
# ==================================================================================================
# FORMAT SCATTERING MODEL PARAMETERS
# ==================================================================================================

    scattering_parameters = InvParameters(MODEL_PARAMETERS)

# ==================================================================================================
# ==================================================================================================
# INITIALIZE
# ==================================================================================================

    INVERSION = InversionMatrix(sv_data_sub, SIMULATION_SETTINGS)

# ==================================================================================================
# ==================================================================================================
# BUILD AND FORMAT SCATTERING MODEL OPTIMIZERS
# ==================================================================================================

    INVERSION.build_scattering_model(scattering_parameters, MODEL_SETTINGS)

# ==================================================================================================
# ==================================================================================================
# BUILD AND FORMAT SCATTERING MODEL OPTIMIZERS
# ==================================================================================================
    res = INVERSION.invert(optimization_kwargs=OPTIMIZATION_KWARGS)

if __name__ == "__main__":
    prof_file = os.path.join(os.getcwd(), "inversion_profile.prof")
    pr = cProfile.Profile()
    pr.enable()
    t0 = time.perf_counter()
    run_workflow()
    elapsed = time.perf_counter() - t0
    pr.disable()
    print(f"Total elapsed: {elapsed:.2f}s")
    pr.dump_stats(prof_file)
    ps = pstats.Stats(pr).strip_dirs().sort_stats("tottime")
    ps.print_stats(40)            # top 40 by total time
    # Save readable text
    with open("inversion_profile_top.txt", "w") as fh:
        ps.stream = fh
        ps.print_stats(200)
    print(f"Profile saved to: {prof_file} and inversion_profile_top.txt")

# snakeviz "c:\Users\Brandyn\Documents\GitHub\echopop\inversion_profile.prof"

# import pstats
# p = pstats.Stats(r"c:\Users\Brandyn\Documents\GitHub\echopop\inversion_profile.prof").strip_dirs().sort_stats("cumtime")
# # show top entries with file:lineno(func) — this prints the file and line info

# # # 2) helper: find Stats keys by substring (file path or function name)
# def find_keys_by_substr(p, substr):
#     return [k for k in p.stats.keys() if substr in (k[0] + k[2])]

# # # Example: find all keys that mention 'highlevel' or your repo 'echopop'
# print(find_keys_by_substr(p, 'highlevel'))      # hot lib entry
# print(find_keys_by_substr(p, r'echopop'))       # your code files

# project_tag = 'echopop'
# agg = {}
# for key in p.stats:
#     callers = p.stats[key][4]
#     for ck in callers:
#         if project_tag in ck[0]:
#             agg[ck] = agg.get(ck, 0.0) + p.stats[ck][3]   # accumulate caller cumulative time

# # print top 20 user-callers
# for ck, cum in sorted(agg.items(), key=lambda kv: kv[1], reverse=True)[:20]:
#     print(ck, "cumtime:", cum)
    
# def heavy_call_chain(p, key, depth=12, project_tag='echopop'):
#     cur = key
#     chain = [cur]
#     for _ in range(depth):
#         callers = p.stats[cur][4]
#         if not callers:
#             break
#         # pick caller with largest cumulative time
#         cur = max(callers.keys(), key=lambda ck: p.stats[ck][3])
#         chain.append(cur)
#         if project_tag in cur[0]:
#             break
#     return chain

# # pick a single hot key (copy one tuple from p.print_stats output)
# hot_key = find_keys_by_substr(p, 'highlevel')[0]
# chain = heavy_call_chain(p, hot_key)
# for item in chain:
#     print(item, "cumtime:", p.stats[item][3])

# hot_key = find_keys_by_substr(p, r'echopop')[1]
# chain = heavy_call_chain(p, hot_key)
# for item in chain:
#     print(item, "cumtime:", p.stats[item][3])


# import inspect
# import echopop.inversion.operations as ops
# print(inspect.getsource(ops.length_average))
    
# # Example: build chain for the first hot key you found
# if hot_keys:
#     chain = heavy_call_chain(p, hot_keys[0])
#     print("Call chain (hot -> caller -> ...):")
#     for item in chain:
#         print(item, "  cumtime:", p.stats[item][3])

# python "c:\Users\Brandyn\Documents\GitHub\echopop\echopop\scratch\profile_inversion.py"

# import line_profiler
# from echopop.inversion.scattering_models import pcdwba, pcdwba_fbs

# profile = line_profiler.LineProfiler()
# profile.add_function(pcdwba)
# profile.add_function(pcdwba_fbs)
# profile.enable_by_count()

# # Call your model function with representative arguments
# res = pcdwba(
#     center_frequencies=np.array([18e3, 38e3, 70e3, 120e3, 200e3]),
#     length_mean=0.030,
#     length_sd_norm=0.15,
#     length_radius_ratio=18.2,
#     taper_order=10.,
#     radius_of_curvature_ratio=3.0,
#     theta_mean=10.,
#     theta_sd=20.,
#     orientation_distribution={"family": "gaussian", "bins": 60},
#     g=1.015,
#     h=1.020,
#     sound_speed_sw=1500.0,
#     frequency_interval=2000.,
#     n_integration=50,
#     n_wavelength=10,
#     number_density=500.,
#     length_distribution={"family": "gaussian", "bins": 100}
# )

# profile.print_stats()

# from echopop.inversion.inversion_matrix_krill import optim

# profile = line_profiler.LineProfiler()
# profile.add_function(optim)
# profile.enable_by_count()

# # Call optim with representative arguments (fill in with actual objects)
# best_fit_set = optim(
#     Sv_measured=self.measurements[["sv_mean", "minimizer", "label"]].loc[1],  # replace with actual pd.Series
#     scattering_parameters=self.model_params,  # replace with actual InvParameters
#     simulation_settings=SIMULATION_SETTINGS,
#     optimization_kwargs=OPTIMIZATION_KWARGS,
#     verbose=True
# )

# profile.print_stats()