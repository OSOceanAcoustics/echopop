import numpy as np

from ..utils.input_checks_read import check_and_read


class LoadStrataData:  # TODO: Does it make sense for this to be a class?
    """
    A Class that loads and processes all
    stratification data. Additionally, it
    calculates the associated mean backscattering
    cross-section for each stratum.

    Parameters
    ----------
    survey : Survey
        An initialized Survey object. Note that any change to
        self.survey will also change this object.
    """

    def __init__(self, survey=None):
        self.survey = survey

        # expected columns for strata Dataframe
        self.strata_cols_types = {
            "stratum_num": int,
            "haul_num": int,
            "fraction_hake": np.float64
        }

        # expected columns for geo strata Dataframe
        self.geo_strata_cols_types = {"stratum_num": int, "northlimit_latitude": np.float64}

        self._load_stratification_file()
        self._load_geographic_stratification()

    def _load_stratification_file(self) -> None:
        """
        Loads and checks the stratification file associated
        with relating the stratification to the Haul.

        Returns
        -------
        strata_df : pd.Dataframe
            Dataframe representation of stratification file.
        """

        if self.survey.params["stratification"]["strata"]["sheetname"] in ["Base KS", "INPFC"]:
            strata_df = check_and_read(
                "stratification/strata",
                self.strata_cols_types,
                self.survey.params
            )

            # set index of dataframe
            strata_df.set_index(["haul_num", "stratum_num"], inplace=True)
            strata_df.sort_index(inplace=True)

            self.survey.strata_df = strata_df
        else:
            raise NotImplementedError("strata_filename has unknown sheet name!")

    def _load_geographic_stratification(self) -> None:
        """
        Loads and checks the geographic stratification
        file defining the stratum and the associated
        Latitude upper limit.

        Returns
        -------
        geo_strata_df : pd.Dataframe
            Dataframe representation of geographic stratification file.
        """

        if self.survey.params["stratification"]["geo_strata"]["sheetname"] in ["stratification1", "INPFC"]: # noqa
            geo_strata_df = check_and_read(
                "stratification/geo_strata",
                self.geo_strata_cols_types,
                self.survey.params
            )

            self.survey.geo_strata_df = geo_strata_df
        else:
            raise NotImplementedError("geo_strata_filename has unknown sheet name!")
