from typing import Union, Optional
import numpy as np
import pandas as pd

from echopop.acoustics import ts_length_regression, to_linear, to_dB

# TODO: Documentation
def configure_transmit_frequency(frequency_values: pd.Series,
                                 transmit_settings: dict, 
                                 current_units: str):
    
    # Extract transmit frequency units defined in configuration file
    configuration_units = transmit_settings["units"]
    
    # Transform the units, if necessary
    # ---- Hz to kHz
    if current_units == "Hz" and configuration_units == "kHz":
        return frequency_values * 1e-3
    # ---- kHz to Hz
    elif current_units == "kHz" and configuration_units == "Hz":
        return frequency_values * 1e3
    # ---- No change
    else:
        return frequency_values
    
# TODO: Documentation
def preprocess_acoustic_data(prc_nasc_df: pd.DataFrame,
                             file_configuration: dict) -> pd.DataFrame:

    # Get acoustic processing settings
    acoustic_analysis_settings = file_configuration["acoustics"]
    # ---- Extract the fined acoustic frequency
    transmit_settings = acoustic_analysis_settings["transmit"]

    # Filter the dataset
    # ---- Configure `frequency_nominal`, if necessary
    prc_nasc_df["frequency_nominal"] = (
        configure_transmit_frequency(prc_nasc_df["frequency_nominal"],
                                     transmit_settings,
                                     acoustic_analysis_settings["dataset_units"]["frequency"])
    )
    # ---- Filter out any unused frequency coordinates
    prc_nasc_df_filtered = (
        prc_nasc_df[prc_nasc_df["frequency_nominal"] == transmit_settings["frequency"]]
    )

    # Remaining adjustments to the acoustic data prior to being passed to the `LiveSurvey` object
    # ---- Replace NASC `NaN` values with `0.0`
    prc_nasc_df_filtered.loc[:, "NASC"] = prc_nasc_df_filtered.loc[:, "NASC"].fillna(0.0)
    # ---- Drop the `frequency_nominal` column and return the output 
    return prc_nasc_df_filtered.drop(columns = ["frequency_nominal"])

# TODO: Documentation
def average_sigma_bs(length: Union[pd.DataFrame, float, int], 
                     weights: Optional[Union[float, int, str]] = None):

    # Function approach for dataframe input
    if isinstance(length, pd.DataFrame):
        if "length" not in length.columns: 
            raise ValueError(
                "Column [`length`] missing from dataframe input `length`."
            )
        elif "TS_L_slope" not in length.columns:
            raise ValueError(
                "Column [`TS_L_slope`] missing from dataframe input `length`."
            )
        elif "TS_L_slope" not in length.columns:
            raise ValueError(
                "Column [`TS_L_intercept`] missing from dataframe input `length`."
            )
        else:           
            # ---- Compute the TS (as an array)
            target_strength = ts_length_regression(length["length"], length["TS_L_slope"], 
                                                   length["TS_L_intercept"])
            # ---- Convert to `sigma_bs`
            sigma_bs_value = to_linear(target_strength)
            # ---- Weighted or arithmetic avveraging
            if weights is None:
                return sigma_bs_value.mean()
            elif weights not in length.columns:
                raise ValueError(
                    f"Defined `weights` column, {weights}, missing from dataframe input "
                    f"`length`."
                )               
            else:
                return (sigma_bs_value * length[weights]).sum() / length[weights].sum()

# TODO: Documentation
# TODO: Refactor
def estimate_echometrics(acoustic_data_df: pd.DataFrame):

    # Create copy
    acoustic_df = acoustic_data_df.copy().reset_index(drop=True)

    # Compute ABC
    # ---- Convert NASC to ABC
    acoustic_df["ABC"] = acoustic_df["NASC"] / (4 * np.pi * 1852 ** 2)

    # Pre-compute the change in depth
    acoustic_df["dz"] = acoustic_df["depth"].diff()

    # Initialize echometrics dictionary
    echometrics = {}

    # Compute the metrics center-of-mass
    if acoustic_df["NASC"].sum() == 0.0:
        echometrics.update({
            "n_layers": 0,
            "mean_Sv": -999,
            "max_Sv": -999,
            "nasc_db": np.nan,
            "center_of_mass": np.nan,
            "dispersion": np.nan,
            "evenness": np.nan,
            "aggregation_index": np.nan,    
            "occupied_area": 0.0,        
        })
    else:
        
        # Create the `echometrics` dictionary 
        echometrics.update({
            # ---- Number of layers
            "n_layers": int(acoustic_df["depth"][acoustic_df["NASC"] > 0.0].size),
            # ---- Mean Sv (back-calculated)
            "mean_Sv": float(
                10.0 * np.log10(acoustic_df["ABC"].sum() / acoustic_df["depth"].max())
            ),
            # ---- Max Sv (back-calculated)
            "max_Sv": float(
                10 * np.log10(acoustic_df["ABC"].max() 
                              / acoustic_df.loc[np.argmax(acoustic_df["ABC"]), "dz"])
            ),
            # ---- (Logarithmic) acoustic abundance
            "nasc_db": float(10 * np.log10(acoustic_df["ABC"].sum())),
            # ---- Center-of-mass
            "center_of_mass": float(
                (acoustic_df["depth"] * acoustic_df["NASC"]).sum() / (acoustic_df["NASC"]).sum()
            ),
            # ---- Evenness
            "evenness": float(
                (acoustic_df["NASC"] **2).sum() / ((acoustic_df["NASC"]).sum()) ** 2
            ),
            # ---- Occupied area
            "occupied_area": float(
                acoustic_df["dz"][acoustic_df["ABC"] > 0.0].sum() / acoustic_df["depth"].max()
            )
        })

        # Update variable-dependent metrics
        echometrics.update({
            # ---- Dispersion
            "dispersion": float(
                ((acoustic_df["depth"] - echometrics["center_of_mass"]) ** 2 
                * acoustic_df["NASC"]).sum() / (acoustic_df["NASC"]).sum()                
            ),
            # ---- Index of aggregation
            "aggregation_index": float(1 / echometrics["evenness"]), 
        })

    # Return the dictionary
    return echometrics

def integrate_nasc(acoustic_data_df: pd.DataFrame, echometrics: bool = True):

    # Vertically integrate PRC NASC
    nasc_dict = {"nasc": acoustic_data_df["NASC"].sum()}
    
    # Horizontally concatenate `echometrics`, if `True`
    if echometrics:
        # ---- Compute values
        # NOTE: This uses NASC instead of linear `sv`
        echometrics_dict = estimate_echometrics(acoustic_data_df)
        # ---- Merge
        nasc_dict.update(echometrics_dict)

    # Convert `nasc_dict` to a DataFrame and return the output
    return pd.Series(nasc_dict)

def compute_nasc(acoustic_data_df: pd.DataFrame, echometrics: bool = True):

    # Integrate NASC (and compute the echometrics, if necessary)
    nasc_data_df = (
        acoustic_data_df.groupby(["longitude", "latitude", "ping_time"])
        .apply(integrate_nasc, echometrics, include_groups=False)
        .unstack().reset_index()
    )
    # ---- Amend the dtypes if echometrics were computed
    if echometrics:
        # ---- Set dtypes
        nasc_data_df = (
            nasc_data_df
            .astype({"n_layers": int, "mean_Sv": float, "max_Sv": float, "nasc_db": float,
                    "center_of_mass": float, "dispersion": float, "evenness": float,
                    "aggregation_index": float, "occupied_area": float})
        )
        # ---- Reorder columns
        nasc_data_df = nasc_data_df[[
            "longitude", "latitude", "ping_time", "nasc", "n_layers", "nasc_db", "mean_Sv", 
            "max_Sv", "aggregation_index", "center_of_mass", "dispersion", "evenness", 
            "occupied_area"
        ]]
