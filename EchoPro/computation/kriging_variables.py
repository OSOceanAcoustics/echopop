import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr


class ComputeKrigingVariables:
    # TODO: see if we want to change the name of this class
    # TODO: this class may not correctly account for bootstrapping!
    """
    A class that computes key variables corresponding to each
    Kriging mesh point, using output produced by Kriging.

    Parameters
    ----------
    krig : ComputeKrigingVariables
        An initialized ComputeKrigingVariables object. Note that any change to
        self.krig will also change this object.
    """

    def __init__(self, krig=None):

        self.krig = krig

    def _set_gender_biomass(self, ds: xr.Dataset) -> None:
        """
        Calculates the biomass for males and females at
        each mesh point. Additionally, adds the corresponding
        variables to ``kriging_results_gdf``.

        Parameters
        ----------
        ds: xr.Dataset
            A Dataset containing all parameters necessary to compute
            desired variables
        """

        # obtain stratum column from Kriging results
        stratum_vals = self.krig.survey.bio_calc.kriging_results_gdf[
            "stratum_num"
        ].values

        # create variables to improve readability
        biomass_density_adult_mean = self.krig.survey.bio_calc.kriging_results_gdf[
            "biomass_density_adult_mean"
        ]
        cell_area_nmi2 = self.krig.survey.bio_calc.kriging_results_gdf["cell_area_nmi2"]

        # calculate the aged biomass for males and females
        for sex in ["M", "F"]:
            dist_weight_sum = (
                (
                    ds[f"len_age_weight_dist_{sex}_normalized"]
                    * ds[f"len_age_weight_prop_{sex}"]
                )
                .sum(dim=["len_bin", "age_bin"])
                .sel(stratum=stratum_vals)
                .values
            )

            biomass_aged = biomass_density_adult_mean * dist_weight_sum * cell_area_nmi2

            # calculate the unaged biomass
            unaged_wgt_prop_expan = (
                ds[f"unaged_{sex}_wgt_proportion"].sel(stratum=stratum_vals).values
            )

            biomass_unaged = (
                biomass_density_adult_mean * unaged_wgt_prop_expan * cell_area_nmi2
            )

            # calculate and assign the total biomass
            total_biomass = biomass_aged + biomass_unaged
            if sex == "M":
                self.krig.survey.bio_calc.kriging_results_male_gdf[
                    "biomass_adult"
                ] = total_biomass
            else:
                self.krig.survey.bio_calc.kriging_results_female_gdf[
                    "biomass_adult"
                ] = total_biomass

    def _set_abundance(self) -> None:
        """
        Calculates the abundance for males, females, and all sexes at
        each mesh point. Additionally, adds the corresponding
        variables to ``kriging_results_gdf``.
        """

        # expand the bio parameters dataframe so that it corresponds to mesh points
        averaged_weight_expanded = (
            self.krig.survey.bio_calc.bio_param_df.averaged_weight.loc[
                self.krig.survey.bio_calc.kriging_results_gdf["stratum_num"].values
            ]
        )

        # calculate and add abundance to Kriging results
        self.krig.survey.bio_calc.kriging_results_gdf["abundance_adult"] = (
            self.krig.survey.bio_calc.kriging_results_gdf["biomass_adult"]
            / averaged_weight_expanded.values
        )

        self.krig.survey.bio_calc.kriging_results_male_gdf["abundance_adult"] = (
            self.krig.survey.bio_calc.kriging_results_male_gdf["biomass_adult"]
            / averaged_weight_expanded.values
        )

        self.krig.survey.bio_calc.kriging_results_female_gdf["abundance_adult"] = (
            self.krig.survey.bio_calc.kriging_results_female_gdf["biomass_adult"]
            / averaged_weight_expanded.values
        )

    def _set_biomass_cell_CV(self):
        """
        Compute the coefficient of Variation (CV) of biomass
        at each grid cell using Kriging output. Additionally,
        assigns the created variable to ``kriging_results_gdf``.
        """

        # create variables to improve readability
        biomass_density_adult = self.krig.survey.bio_calc.transect_results_gdf[
            "biomass_density_adult"
        ].values
        biomass_density_adult_mean = self.krig.survey.bio_calc.kriging_results_gdf[
            "biomass_density_adult_mean"
        ]
        biomass_density_adult_var = self.krig.survey.bio_calc.kriging_results_gdf[
            "biomass_density_adult_var"
        ]
        cell_area_nmi2 = self.krig.survey.bio_calc.kriging_results_gdf["cell_area_nmi2"]
        kriging_A0 = self.krig.survey.params["kriging_A0"]

        C0 = np.std(biomass_density_adult, ddof=1) ** 2
        Bn = np.nansum(biomass_density_adult_mean * cell_area_nmi2) * 1e-9
        self.krig.survey.bio_calc.kriging_results_gdf["biomass_adult_cell_CV"] = (
            kriging_A0
            * np.sqrt(biomass_density_adult_var * C0)
            * 1e-9
            / Bn
            * np.sqrt(len(biomass_density_adult_var))
        )

    def _compute_biomass_all_ages(
        self,
        weight_fraction_all_ages_df: pd.DataFrame,
        biomass_column: pd.DataFrame,
        results_gdf: gpd.GeoDataFrame,
    ) -> None:
        """
        Compute the biomass for each age bin based off the provided input.
        Additionally, add computed variables to the input ``results_gdf``.

        Parameters
        ----------
        weight_fraction_all_ages_df: pd.DataFrame
            A DataFrame containing the weight fraction at each age
            bin for all strata (corresponds to ``biomass_column``)
        biomass_column: pd.DataFrame
            A DataFrame corresponding to a biomass column from the Kriging results
        results_gdf: gpd.GeoDataFrame
            A GeoDataFrame where the biomass at each age bin should be stored
        """

        # obtain stratum column from Kriging results
        stratum_vals = results_gdf["stratum_num"].values

        for bin_str in weight_fraction_all_ages_df:

            # expand the weight fraction for all ages
            expanded_age_bin = (
                weight_fraction_all_ages_df[bin_str].loc[stratum_vals].values
            )

            results_gdf["biomass_" + bin_str] = expanded_age_bin * biomass_column

    def set_variables(self) -> None:
        """
        Calculates variables over Kriging mesh points that are useful
        for analysis (e.g. abundance, NASC, CV). Additionally, assigns
        these variables to the appropriate Kriging results GeoDataFrame.
        """

        # calculate and add the male and female biomass to Kriging results
        self._set_gender_biomass(self.krig.survey.bio_calc.param_ds)

        # calculate and add abundance variables to Kriging results
        self._set_abundance()

        # add sig_b values to Kriging results
        self.krig.survey.bio_calc.kriging_results_gdf[
            "sig_b"
        ] = self.krig.survey.bio_calc.strata_sig_b.loc[
            self.krig.survey.bio_calc.kriging_results_gdf["stratum_num"].values
        ].values

        # calculate and add NASC to Kriging results
        self.krig.survey.bio_calc.kriging_results_gdf["NASC"] = (
            self.krig.survey.bio_calc.kriging_results_gdf["abundance_adult"]
            * self.krig.survey.bio_calc.kriging_results_gdf["sig_b"]
        )

        # calculate and add biomass CV at each grid cell to Kriging results
        self._set_biomass_cell_CV()

        # calculate and add male biomass for all ages to Kriging results
        self._compute_biomass_all_ages(
            self.krig.survey.bio_calc.weight_fraction_all_ages_male_df,
            self.krig.survey.bio_calc.kriging_results_male_gdf["biomass_adult"],
            self.krig.survey.bio_calc.kriging_results_male_gdf,
        )

        # calculate and add female biomass for all ages to Kriging results
        self._compute_biomass_all_ages(
            self.krig.survey.bio_calc.weight_fraction_all_ages_female_df,
            self.krig.survey.bio_calc.kriging_results_female_gdf["biomass_adult"],
            self.krig.survey.bio_calc.kriging_results_female_gdf,
        )

        # calculate and add biomass for all ages to Kriging results
        self._compute_biomass_all_ages(
            self.krig.survey.bio_calc.weight_fraction_all_ages_df,
            self.krig.survey.bio_calc.kriging_results_gdf["biomass_adult"],
            self.krig.survey.bio_calc.kriging_results_gdf,
        )
