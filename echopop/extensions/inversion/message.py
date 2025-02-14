import numpy as np
import pandas as pd

def inversion_kriging_results_msg(kriging_results_dict: pd.DataFrame, settings_dict: dict) -> None:

    # Extract dictionary results
    kriging_mesh_results = kriging_results_dict

    # Break down strings
    # ---- Mesh cropping
    if not settings_dict["extrapolate"]:
        crop_str = (
            f"Mesh cropping method: "
            f"{settings_dict['cropping_parameters']['crop_method'].capitalize().replace('_', ' ')}"
        )
    else:
        crop_str = "Extrapolated over uncropped grid"

    # Generate message output
    return print(
        f"--------------------------------\n"
        f"(INVERTED) KRIGING RESULTS (MESH)\n"
        f"--------------------------------\n"
        f"| Kriged variable: {
            ('Inverted ' 
            + settings_dict['variable'].replace('_', ' ')).capitalize()
        } "
        f"(kg/nmi^2)\n"
        f"| Mesh extrapolation: {settings_dict['extrapolate']}\n"
        f"    {crop_str}\n"
        f"| Mesh and transect coordinate standardization: "
        f"{settings_dict['standardize_coordinates']}\n"
        f"--------------------------------\n"
        f"GENERAL RESULTS\n"
        f"--------------------------------\n"
        f"| Mean {settings_dict['variable'].replace('_', ' ')}: "
        f"{np.round(kriging_mesh_results['survey_mean'], 2)} kg/nmi^2\n"
        f"| Total survey biomass estimate: "
        f"{np.round(kriging_mesh_results['survey_estimate'] * 1e-6, 2)} kmt\n"
        f"| Mean mesh sample CV: "
        f"{np.round(kriging_mesh_results['mesh_results_df']['sample_cv'].mean(), 4)}\n"
        f"| Overall survey CV: {np.round(kriging_mesh_results['survey_cv'], 4)}\n"
        f"| Total area coverage: "
        f"{np.round(kriging_mesh_results['mesh_results_df']['area'].sum(), 1)} nmi^2\n"
        f"--------------------------------"
    )
