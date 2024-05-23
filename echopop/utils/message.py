import copy

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
    | Stratum definition: {settings_dict[ 'stratum' ].upper() }
    --------------------------------
    GENERAL RESULTS
    --------------------------------
    Total biomass: {biomass_message.loc[ 'all' ][ 'biomass_all' ]}{units}
        Age-1: {biomass_message.loc[ 'all' ][ 'biomass_age1' ]}{units}
        Age-2+: {biomass_message.loc[ 'all' ][ 'biomass_adult' ]}{units}
    Total female biomass: {biomass_message.loc[ 'female' ][ 'biomass_all' ]}{units}
        Age-1: {biomass_message.loc[ 'female' ][ 'biomass_age1' ]}{units}
        Age-2+: {biomass_message.loc[ 'female' ][ 'biomass_adult' ]}{units}
    Total male biomass: {biomass_message.loc[ 'male' ][ 'biomass_all' ]}{units}
        Age-1: {biomass_message.loc[ 'male' ][ 'biomass_age1' ]}{units}
        Age-2+: {biomass_message.loc[ 'male' ][ 'biomass_adult' ]}{units}
    Total unsexed biomass: {biomass_message.loc[ 'unsexed' ][ 'biomass_all' ]}{units}
    Total mixed biomass: {biomass_message.loc[ 'mixed' ][ 'biomass_all' ]}{units}
    --------------------------------"""
    )


def stratified_results_msg(stratified_results_dict: pd.DataFrame, settings_dict: dict) -> None:

    # Create copy
    stratified_results = copy.deepcopy(stratified_results_dict)

    # Apply back-transformation of total
    # ---- Initialize
    stratified_results.update({"back_transform": copy.deepcopy(stratified_results["total"])})
    # ---- Adjust the estimate and confidence interval values
    stratified_results["back_transform"].update(
        {
            "estimate": (
                stratified_results["back_transform"]["estimate"] / settings_dict["transect_sample"]
            ),
            "confidence_interval": (
                stratified_results["back_transform"]["confidence_interval"]
                / settings_dict["transect_sample"]
            ),
        }
    )

    # Create confidence interval string
    # ---- Helper function
    def format_values(value, ci_pct, variable, statistic, metric):
        if "confidence_interval" in metric:
            if statistic == "cv":
                return f"""[{', '.join(map(lambda x: str(x.round(4)), value))}; \
{int(ci_pct * 100)}% CI]"""
            elif (
                statistic in ["back_transform", "total", "mean", "variance"]
                and variable == "biomass"
            ):
                return f"""[{', '.join(map(lambda x: str(np.round(x* 1e-6,1)), value))}; \
{int(ci_pct * 100)}% CI]"""
            elif (
                statistic in ["back_transform", "total", "mean", "variance"]
                and variable == "abundance"
            ):
                return f"""[{', '.join(map(lambda x: str(int(np.round(x,0))), value))}; \
{int(ci_pct * 100)}% CI]"""
            else:
                return f"""[{', '.join(map(lambda x: str(np.round(x,2)), value))}; \
{int(ci_pct * 100)}% CI]"""
        elif "estimate" in metric:
            if statistic == "cv":
                return f"{str(np.round(value, 4))}"
            elif statistic in ["back_transform", "total", "mean"] and variable == "biomass":
                return f"{str(np.round(value*1e-6,1))} kmt"
            elif statistic == "variance" and variable == "biomass":
                return f"{str(np.round(value*1e-6,1))} kmt^2"
            elif (
                statistic in ["back_transform", "total", "mean", "variance"]
                and variable == "abundance"
            ):
                return f"{str(int(np.round(value,0)))} fish"
            elif statistic == "variance" and variable == "abundance":
                return f"{str(int(np.round(value,0)))} fish^2"
            elif (
                statistic in ["back_transform", "total", "mean", "variance"]
                and variable == "biomass_density"
            ):
                return f"{str(np.round(value*1e-6,1))} kmt/nmi^2"
            elif statistic == "variance" and variable == "biomass_density":
                return f"{str(np.round(value*1e-6,1))} (kmt/nmi^2)^2"
            elif (
                statistic in ["back_transform", "total", "mean", "variance"]
                and variable == "number_density"
            ):
                return f"{str(int(np.round(value,0)))} fish/nmi^2"
            elif statistic == "variance" and variable == "number_density":
                return f"{str(int(np.round(value,0)))} (fish/nmi^2)^2"

    # ---- Initialize the message reference dictionary
    stratified_message = {}
    # ---- Format the CI's (and initialize the updated stratified results)
    for statistic, statistic_data in stratified_results.items():
        if isinstance(statistic_data, dict):
            for metric, metric_value in statistic_data.items():
                stratified_message[f"{statistic}:{metric}"] = format_values(
                    metric_value,
                    stratified_results["ci_percentile"],
                    stratified_results["variable"],
                    statistic,
                    metric,
                )

    # TODO: Add to the `stratified_analysis(...)` output: estimates for each of the strata
    # Generate message output
    return print(
        f"""
    --------------------------------
    STRATIFIED RESULTS ({settings_dict[ 'dataset' ].upper()})
    --------------------------------
    | Stratified variable: {settings_dict[ 'variable' ].title() } (kmt)
    | Number of { 'virtual transects'
if settings_dict[ 'dataset' ] == 'kriging' else 'transects' }: \
{ stratified_results[ 'num_transects' ] }
    | Total strata area coverage: { np.round(stratified_results[ 'total_area' ],1)} nmi^2
    | Age-1 fish excluded: {settings_dict['exclude_age1']}
    | Stratum definition: {settings_dict[ 'stratum' ].upper() }
    | Bootstrap replicates: {settings_dict[ 'transect_replicates' ] } samples
    | Resampling proportion: {settings_dict[ 'transect_sample' ] }
    --------------------------------
    STRATUM-SPECIFIC ESTIMATES
    --------------------------------
    _coming soon_
    --------------------------------
    GENERAL RESULTS
    --------------------------------
    CV: {stratified_message['cv:estimate']} {stratified_message['cv:confidence_interval']}\n"""
        f"""    Total (across sub-sampled transects): {stratified_message['total:estimate']}"""
        f""" {stratified_message['total:confidence_interval']}\n"""
        f"""          (back-transformed): {stratified_message['back_transform:estimate']}"""
        f""" {stratified_message['back_transform:confidence_interval']}\n"""
        f"""    Mean (across sub-sampled transects):\
        {stratified_message['mean:unweighted_estimate']}"""
        f""" {stratified_message['mean:unweighted_confidence_interval']}
    --------------------------------"""
    )


def kriging_results_msg(kriging_results_dict: pd.DataFrame, settings_dict: dict) -> None:

    # Extract dictionary results
    kriging_mesh_results = kriging_results_dict

    # Generate message output
    return print(
        f"""
    --------------------------------
    KRIGING RESULTS (MESH)
    --------------------------------
    | Kriged variable: {settings_dict[ 'variable' ].replace("_"," ").capitalize()} (kg/nmi^2)
    | Age-1 fish excluded: {settings_dict['exclude_age1']}
    | Stratum definition: {settings_dict[ 'stratum' ].upper() }
    | Mesh extrapolation: {settings_dict[ 'extrapolate' ] }
    --- Mesh cropping method: { settings_dict[ 'crop_method' ].capitalize() if not
    settings_dict[ 'extrapolate' ] else None }
    | Mesh and transect coordinate standardization: {settings_dict[ 'standardize_coordinates' ] }"""
        """
    --------------------------------
    GENERAL RESULTS
    --------------------------------\n"""
        f"""    Mean {settings_dict[ 'variable' ].replace("_"," ")}: """
        f"""{np.round(kriging_mesh_results['survey_mean'],2)} kg/nmi^2\n"""
        f"""    Total survey {settings_dict[ 'variable' ].replace("_density","")} estimate: """
        f"""{np.round(kriging_mesh_results['survey_estimate']*1e-6,2)} kmt\n"""
        f"""    Mean mesh sample CV: """
        f"""{np.round(kriging_mesh_results[ 'mesh_results_df' ][ 'sample_cv' ].mean(),4)}\n"""
        f"""    Overall survey CV: {np.round(kriging_mesh_results['survey_cv'], 4)}\n"""
        f"""    Total area coverage:"""
        f""" {np.round(kriging_mesh_results[ 'mesh_results_df' ]['area'].sum(),1)} nmi^2\n"""
        """   --------------------------------"""
    )