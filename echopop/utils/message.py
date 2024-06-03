import numpy as np
import pandas as pd


def transect_results_msg(transect_results_dict: dict, settings_dict: dict) -> None:

    # Create copy of biomass summary dataframe
    biomass_message = transect_results_dict["biomass_summary_df"].copy()

    # Define unit string
    units = " kmt"

    # Convert entire summary dataframe to units kmt
    # ---- Convert from kg to kmt
    biomass_message.iloc[:, 1:] = biomass_message.iloc[:, 1:] * 1e-6
    # ---- Format (round to nearest tenth)
    biomass_message.iloc[:, 1:] = biomass_message.iloc[:, 1:].round(1)
    # ---- Index by sex
    biomass_message.set_index(["sex"], inplace=True)

    # Generate message output
    return print(
        f"""
    --------------------------------
    TRANSECT RESULTS
    --------------------------------
    | Variable: Biomass (kmt)
    | Age-1 fish excluded: {settings_dict['exclude_age1']}
    | Stratum definition: {settings_dict['stratum'].upper()}
    --------------------------------
    GENERAL RESULTS
    --------------------------------
    Total biomass: {biomass_message.loc['all']['biomass_all']}{units}
        Age-1: {biomass_message.loc['all']['biomass_age1']}{units}
        Age-2+: {biomass_message.loc['all']['biomass_adult']}{units}
    Total female biomass: {biomass_message.loc['female']['biomass_all']}{units}
        Age-1: {biomass_message.loc['female']['biomass_age1']}{units}
        Age-2+: {biomass_message.loc['female']['biomass_adult']}{units}
    Total male biomass: {biomass_message.loc['male']['biomass_all']}{units}
        Age-1: {biomass_message.loc['male']['biomass_age1']}{units}
        Age-2+: {biomass_message.loc['male']['biomass_adult']}{units}
    Total unsexed biomass: {biomass_message.loc['unsexed']['biomass_all']}{units}
    Total mixed biomass: {biomass_message.loc['mixed']['biomass_all']}{units}
    --------------------------------"""
    )


def stratified_results_msg(stratified_results_dict: pd.DataFrame, settings_dict: dict) -> None:

    # Assign the max number of results per line for multi-line entries
    max_results_per_row = 8

    # Determine variable units
    if settings_dict["variable"] == "biomass":
        units = "kmt"
        bio_units_round = 1
        bio_units_mult = 1e-6
    elif settings_dict["variable"] == "abundance":
        units = "number"
        bio_units_round = 0
        bio_units_mult = 1
    else:
        units = "m^2 nmi^-2"
        bio_units_round = 1
        bio_units_mult = 1

    # Breakdown strings
    # ---- Transect variable name
    transect_name = "virtual transects" if settings_dict["dataset"] == "kriging" else "transects"
    # ---- Number of strata
    strata_num = stratified_results_dict["estimate"]["strata"]["total"].size
    # ---- Strata area
    strata_area = " | ".join(str(x) for x in stratified_results_dict["stratum_area"].round())
    # ---- Strata densities
    # -------- Format the estimates
    strata_densities = (
        stratified_results_dict["estimate"]["strata"]["density"] * bio_units_mult
    ).round(bio_units_round + 2)
    # -------- Format the CIs
    strata_density_cis = stratified_results_dict["ci"]["strata"]["density"]
    strata_density_cis_scaled = [
        np.round(ci * bio_units_mult, bio_units_round + 2) for ci in strata_density_cis
    ]
    # -------- Create the string
    strata_density_cis_str = [
        f"{mean} [{ci[0]}, {ci[1]}]"
        for mean, ci in zip(strata_densities, strata_density_cis_scaled)
    ]
    # -------- Join
    density_rows = []
    for i in range(0, len(strata_density_cis_str), max_results_per_row):
        row = " | ".join(strata_density_cis_str[i : i + max_results_per_row])
        density_rows.append(row)
    # ---- Strata totals
    # -------- Format the estimates
    strata_totals = (stratified_results_dict["estimate"]["strata"]["total"] * bio_units_mult).round(
        bio_units_round
    )
    # -------- Format the CIs
    strata_total_cis = stratified_results_dict["ci"]["strata"]["total"]
    strata_total_cis_scaled = [
        np.round(ci * bio_units_mult, bio_units_round) for ci in strata_total_cis
    ]
    # -------- Create the string
    strata_total_cis_str = [
        f"{mean} [{ci[0]}, {ci[1]}]" for mean, ci in zip(strata_totals, strata_total_cis_scaled)
    ]
    # -------- Join
    total_rows = []
    for i in range(0, len(strata_total_cis_str), max_results_per_row):
        row = " | ".join(strata_total_cis_str[i : i + max_results_per_row])
        total_rows.append(row)
    # ---- Survey densities
    # -------- Format the estimates
    survey_densities = (
        stratified_results_dict["estimate"]["survey"]["density"] * bio_units_mult
    ).round(bio_units_round + 2)
    # -------- Format the CIs
    survey_density_cis = stratified_results_dict["ci"]["survey"]["density"]
    survey_density_cis_scaled = [
        np.round(ci * bio_units_mult, bio_units_round + 2) for ci in survey_density_cis
    ]
    # -------- Create the string
    survey_density_cis_str = (
        f"{survey_densities} [{survey_density_cis_scaled[0]}, {survey_density_cis_scaled[1]}]"
    )
    # ---- Survey totals
    # -------- Format the estimates
    survey_total = (stratified_results_dict["estimate"]["survey"]["total"] * bio_units_mult).round(
        bio_units_round
    )
    # -------- Format the CIs
    survey_total_cis = stratified_results_dict["ci"]["survey"]["total"]
    survey_total_cis_scaled = [
        np.round(ci * bio_units_mult, bio_units_round) for ci in survey_total_cis
    ]
    # -------- Create the string
    survey_total_cis_str = (
        f"{survey_total} [{survey_total_cis_scaled[0]}, {survey_total_cis_scaled[1]}]"
    )
    # ---- Survey CVs
    # -------- Format the estimates
    survey_cv = (stratified_results_dict["estimate"]["survey"]["cv"]).round(4)
    # -------- Format the CIs
    survey_cv_cis = stratified_results_dict["ci"]["survey"]["cv"]
    survey_cv_cis_scaled = [np.round(ci, 4) for ci in survey_cv_cis]
    # -------- Create the string
    survey_cv_cis_str = f"{survey_cv} [{survey_cv_cis_scaled[0]}, {survey_cv_cis_scaled[1]}]"
    # Fix strings
    density_str = "\n".join(density_rows)
    total_str = "\n".join(total_rows)
    # Generate message output
    return print(
        f"--------------------------------\n"
        f" STRATIFIED RESULTS ({settings_dict['dataset'].upper()})\n"
        f"--------------------------------\n"
        f"| Stratified variable: {settings_dict['variable'].title()} ({units})\n"
        f"| Number of {transect_name}: {stratified_results_dict['num_transects']}\n"
        f"| Number of strata ({settings_dict['stratum'].upper()}): {strata_num}\n"
        f"| Total area coverage: {stratified_results_dict['total_area'].round()} nmi^2\n"
        f"| Age-1 fish excluded: {settings_dict['exclude_age1']}\n"
        f"| Bootstrap replicates: {settings_dict['transect_replicates']} samples\n"
        f"| Resampling proportion: {settings_dict['transect_sample']}\n"
        f"| Bootstrap interval method: {settings_dict['bootstrap_ci_method']} (CI: "
        f"{settings_dict['bootstrap_ci'] * 1e2}%)\n"
        f"--------------------------------\n"
        f"STRATUM-SPECIFIC ESTIMATES\n"
        f"--------------------------------\n"
        f"| Stratum area coverage:\n"
        f"{strata_area} nmi^2\n"
        f"| Stratum mean {settings_dict['variable']} density ({units}/nmi^2):\n"
        f"{density_str}\n"
        f"| Stratum mean {settings_dict['variable']} ({units}):\n"
        f"{total_str}\n"
        f"--------------------------------\n"
        f"SURVEY RESULTS\n"
        f"--------------------------------\n"
        f"| Survey mean {settings_dict['variable']} density ({units}/nmi^2):"
        f" {survey_density_cis_str}\n"
        f"| Survey mean {settings_dict['variable']} ({units}):"
        f" {survey_total_cis_str}\n"
        f"| Survey CV: {survey_cv_cis_str}"
    )


def kriging_results_msg(kriging_results_dict: pd.DataFrame, settings_dict: dict) -> None:

    # Extract dictionary results
    kriging_mesh_results = kriging_results_dict

    # Break down strings
    # ---- Mesh cropping
    mesh_crop = (
        settings_dict["crop_method"].capitalize() if not settings_dict["extrapolate"] else None
    )

    # Generate message output
    return print(
        f"""
    --------------------------------
    KRIGING RESULTS (MESH)
    --------------------------------
    | Kriged variable: {settings_dict['variable'].replace('_', ' ').capitalize()} (kg/nmi^2)
    | Age-1 fish excluded: {settings_dict['exclude_age1']}
    | Stratum definition: {settings_dict['stratum'].upper()}
    | Mesh extrapolation: {settings_dict['extrapolate']}
    --- Mesh cropping method: {mesh_crop}
    | Mesh and transect coordinate standardization: {settings_dict['standardize_coordinates']}"""
        """
    --------------------------------
    GENERAL RESULTS
    --------------------------------\n"""
        f"""    Mean {settings_dict['variable'].replace("_", " ")}: """
        f"""{np.round(kriging_mesh_results['survey_mean'], 2)} kg/nmi^2\n"""
        f"""    Total survey {settings_dict['variable'].replace("_density", "")} estimate: """
        f"""{np.round(kriging_mesh_results['survey_estimate'] * 1e-6, 2)} kmt\n"""
        f"""    Mean mesh sample CV: """
        f"""{np.round(kriging_mesh_results['mesh_results_df']['sample_cv'].mean(), 4)}\n"""
        f"""    Overall survey CV: {np.round(kriging_mesh_results['survey_cv'], 4)}\n"""
        f"""    Total area coverage:"""
        f""" {np.round(kriging_mesh_results['mesh_results_df']['area'].sum(), 1)} nmi^2\n"""
        """   --------------------------------"""
    )
