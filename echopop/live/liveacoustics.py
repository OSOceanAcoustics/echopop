from typing import Union, Optional

import pandas as pd

from echopop.acoustics import ts_length_regression, to_linear, to_dB

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
            "aggregation": np.nan,    
            "occupied_area": 0.0,        
        })
    else:
        
        # Compute the number of layers
        echometrics.update({
            "n_layers": acoustic_df["depth"][acoustic_df["NASC"] > 0.0].size
        })

        # Compute ABC
        # ---- Convert NASC to ABC
        acoustic_df["ABC"] = acoustic_df["NASC"] / (4 * np.pi * 1852 ** 2)
        # ---- Estimate mean Sv
        echometrics.update({
            "mean_Sv": 10.0 * np.log10(acoustic_df["ABC"].sum() / acoustic_df["depth"].max())
        })
        # --- Estimate max Sv (i.e. )
        echometrics.update({
            "max_Sv": 10 * np.log10(acoustic_df["ABC"].max() 
                                    / acoustic_df.loc[np.argmax(acoustic_df["ABC"]), "dz"])
        })

        # Compute (acoustic) abundance
        echometrics.update({
            "nasc_db": 10 * np.log10(acoustic_df["ABC"].sum())
        })

        # Compute center of mass
        echometrics.update({
            "center_of_mass": (
                (acoustic_df["depth"] * acoustic_df["NASC"]).sum()
                / (acoustic_df["NASC"]).sum()
            )
        })

        # Compute the dispersion
        echometrics.update({
            "dispersion": (
                ((acoustic_df["depth"] - echometrics["center_of_mass"]) ** 2 
                * acoustic_df["NASC"]).sum() / (acoustic_df["NASC"]).sum()                
            )
        })

        # Compute the evenness
        echometrics.update({
            "evenness": (acoustic_df["NASC"] **2).sum() / ((acoustic_df["NASC"]).sum()) ** 2
        })

        # Compute the index of aggregation
        echometrics.update({
            "aggregation": 1 / echometrics["evenness"]
        })

        # Get the occupied area
        echometrics.update({
            "occupied_area": (
                acoustic_df["dz"][acoustic_df["ABC"] > 0.0].sum() / acoustic_df["depth"].max()
            )
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
