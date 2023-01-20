import pathlib
from typing import Union

import numpy as np
import pandas as pd


class Reports:
    """
    A Class that writes requested variables to
    consolidated files.

    Parameters
    ----------
    survey : Survey
        Survey object that contains all necessary parameters

    """

    def __init__(self, survey):

        self.survey = survey

    def __process_gear_data(self, df):
        """
        Parameters
        ----------
        df : Pandas Dataframe
            Dataframe holding the gear data
        Returns
        -------
        Processed Dataframe
        """

        df = df[
            [
                "Haul",
                "Average_Footrope_Depth",
                "Surface_Temperature",
                "Gear_Temperature",
                "Average_Wireout",
                "Transect",
                "Net_Height",
            ]
        ].copy()

        # set data types of dataframe
        df = df.astype(
            {
                "Haul": int,
                "Average_Footrope_Depth": np.float64,
                "Surface_Temperature": np.float64,
                "Gear_Temperature": np.float64,
                "Average_Wireout": np.float64,
                "Transect": np.float64,
            }
        )

        if self.EPro.params["exclude_age1"] is False:
            raise NotImplementedError("Including age 1 data has not been implemented!")

        df["Net_Height"] = (
            df["Net_Height"]
            .apply(lambda x: np.nan if type(x) == str else x)
            .astype(float)
        )

        df.set_index("Haul", inplace=True)
        df.sort_index(inplace=True)

        return df

    def __load_gear_data(self):

        # TODO: if this function ends up being used, we need to document and create a check for df

        if self.EPro.params["source"] == 3:

            if self.EPro.params["filename_gear_US"]:
                gear_us_df = pd.read_excel(
                    self.EPro.params["data_root_dir"]
                    + self.EPro.params["filename_gear_US"],
                    sheet_name="biodata_gear",
                )
                gear_us_df = self.__process_gear_data(gear_us_df)
            else:
                gear_us_df = None

            if self.EPro.params["filename_gear_CAN"]:
                gear_can_df = pd.read_excel(
                    self.EPro.params["data_root_dir"]
                    + self.EPro.params["filename_gear_CAN"],
                    sheet_name="biodata_gear_CAN",
                )
                gear_can_df = self.__process_gear_data(gear_can_df)

                # transect offset is set to zero since we will not have overlap transects
                # TODO: Is this a necessary variable? If so, we should make it an input to the function.  # noqa
                CAN_Transect_offset = 0

                # combine US & CAN gear files
                gear_can_df.index = (
                    gear_can_df.index + self.EPro.params["CAN_haul_offset"]
                )  # add haul_offset
                gear_can_df["Transect"] = gear_can_df["Transect"] + CAN_Transect_offset
            else:
                gear_can_df = None

            # combine US & CAN trawl files
            if isinstance(gear_us_df, pd.DataFrame) and isinstance(
                gear_can_df, pd.DataFrame
            ):
                self.EPro.gear_df = pd.concat([gear_us_df, gear_can_df])
            else:
                raise SystemError(
                    "Cannot construct gear_df for source = 3, "
                    "since either US or CAN is not available."
                )

        else:
            raise NotImplementedError(
                f"Source of {self.EPro.params['source']} not implemented yet."
            )

    def __process_trawl_data(self, df):
        """
        Parameters
        ----------
        df : Pandas Dataframe
            Dataframe holding the trawl data
        Returns
        -------
        Processed Dataframe
        """

        df = df[
            [
                "Haul",
                "Haul_Type",
                "Performance_Code",
                "Duration",
                "Distance_Fished",
                "Stratum",
                "EQ_Latitude",
                "EQ_Longitude",
                "Average_Bottom_Depth",
                "Vessel_Log_Start",
                "Vessel_Log_Stop",
            ]
        ].copy()

        if df.dtypes["Average_Bottom_Depth"] == object:
            # remove any instances of >, <, or m from the column Average_Bottom_Depth
            df["Average_Bottom_Depth"] = (
                df["Average_Bottom_Depth"]
                .apply(lambda x: x.strip("><m") if type(x) == str else x)
                .astype(float)
            )

        df["Duration"] = pd.to_timedelta(
            df["Duration"]
        )  # change datatype from string to timedelta

        # set data types of dataframe
        df = df.astype(
            {
                "Haul": int,
                "Haul_Type": int,
                "Performance_Code": np.float64,
                "Distance_Fished": np.float64,
                "Stratum": int,
                "EQ_Latitude": np.float64,
                "EQ_Longitude": np.float64,
                "Average_Bottom_Depth": np.float64,
                "Vessel_Log_Start": np.float64,
                "Vessel_Log_Stop": np.float64,
            }
        )

        df.set_index("Haul", inplace=True)
        df.sort_index(inplace=True)

        if self.EPro.params["exclude_age1"] is False:
            raise NotImplementedError("Including age 1 data has not been implemented!")

        if self.EPro.params["hemisphere"][0] == "N":
            df["EQ_Latitude"] = df["EQ_Latitude"].abs()
        elif self.EPro.params["hemisphere"][0] == "S":
            df["EQ_Latitude"] = -1.0 * (df["EQ_Latitude"].abs())
        else:
            raise ValueError("Wrong N/S Hemisphere provided! \n")

        if self.EPro.params["hemisphere"][1] == "W":
            df["EQ_Longitude"] = -1.0 * (df["EQ_Longitude"].abs())
        elif self.EPro.params["hemisphere"][1] == "E":
            df["EQ_Longitude"] = df["EQ_Longitude"].abs()
        else:
            raise ValueError("Wrong E/W Hemisphere provided! \n")

        return df

    def __load_trawl_data(self):

        # TODO: if this function ends up being used, we need to document and create a check for df

        if self.EPro.params["source"] == 3:

            if self.EPro.params["filename_trawl_US"]:
                trawl_us_df = pd.read_excel(
                    self.EPro.params["data_root_dir"]
                    + self.EPro.params["filename_trawl_US"],
                    converters={"Duration": str},
                    sheet_name="biodata_haul",
                )
                trawl_us_df = self.__process_trawl_data(trawl_us_df)
            else:
                trawl_us_df = None

            if self.EPro.params["filename_trawl_CAN"]:
                trawl_can_df = pd.read_excel(
                    self.EPro.params["data_root_dir"]
                    + self.EPro.params["filename_trawl_CAN"],
                    converters={"Duration": str},
                    sheet_name="biodata_haul_CAN",
                )
                trawl_can_df = self.__process_trawl_data(trawl_can_df)

                trawl_can_df.index = (
                    trawl_can_df.index + self.EPro.params["CAN_haul_offset"]
                )  # add haul_offset

                # change lon to -lon
                trawl_can_df["EQ_Longitude"] = trawl_can_df["EQ_Longitude"].apply(
                    lambda x: -1.0 * x if x > 0 else x
                )
            else:
                trawl_can_df = None

            # combine US & CAN trawl files
            if isinstance(trawl_us_df, pd.DataFrame) and isinstance(
                trawl_can_df, pd.DataFrame
            ):
                self.EPro.trawl_df = pd.concat([trawl_us_df, trawl_can_df])
            else:
                raise SystemError(
                    "Cannot construct trawl_df for source = 3, "
                    "since either US or CAN is not available."
                )

        else:
            raise NotImplementedError(
                f"Source of {self.EPro.params['source']} not implemented yet."
            )

    def _aged_len_haul_counts_report(self):
        """
        Creates aged length-haul-counts table, which specifies the aged hake
        counts at each length bin from all hauls for male, female, and all.
        Additionally, writes this report to an Excel file.

        Returns
        -------

        """

        pass

    def _transect_based_biomass_ages_report(self):
        """
        Generates a transect based biomass at different age bins
        report and writes it to an Excel file.

        Returns
        -------

        """

        pass

    def _transect_based_core_variables_report(self):
        """
        Generates a report containing core transect based variables
        and writes it to an Excel file.

        Returns
        -------

        """

        pass

    def _total_len_haul_counts_report(self):
        """
        Generates a report for total hake counts at length with
        actually measured length (sampled at both biological sampling
        stations for US/CAN) for all hauls. Additionally, writes this
        report to an Excel file.

        Returns
        -------

        """

        pass

    def _transect_based_len_age_abundance_report(self):
        """
        Generates a report that is a 40 x 21 matrix, with 40 length bins (1st column)
        (L =2, 4, 6, … 80 cm, and ith length bin include length i-1 and i cm) and 21
        age bins (1st row), 1, 2, 3, … 20 year old, and the 21st age bin is un-aged
        abundance. The abundance values are based on the transect results. There are
        three tabs with male, female, and all (male+female). Additionally, writes this
        report to an Excel file.

        Returns
        -------

        """

        pass

    def _transect_based_len_age_biomass_report(self):
        """
        Generates a report that is a 40 x 21 matrix, with 40 length bins (1st column)
        (L =2, 4, 6, … 80 cm, and ith length bin include length i-1 and i cm) and 21
        age bins (1st row), 1, 2, 3, … 20 year old. The biomass values are based on the
        transect results. There are three tabs with male, female, and all (male+female).
        Additionally, writes this report to an Excel file.

        Returns
        -------

        """

        pass

    def create_and_write_reports(self, output_path: Union[str, pathlib.Path]) -> None:
        """
        Constructs Kriging mesh and Transect report DataFrames and writes
        them to Excel files.

        Parameters
        ----------
        output_path: str or pathlib.Path
            The output path where all Excel files should be saved

        # TODO: maybe include an option to overwrite files (default to False)
        """

        if not isinstance(output_path, (str, pathlib.Path)):
            raise TypeError("output_path must be a string or pathlib.Path!")

        if isinstance(output_path, str):
            # convert string path to pathlib.Path
            output_path = pathlib.Path(output_path)

        # check if path exists, if it doesn't create it
        output_path.mkdir(parents=True, exist_ok=True)

        self._aged_len_haul_counts_report()

        self._transect_based_biomass_ages_report()

        self._transect_based_core_variables_report()

        self._total_len_haul_counts_report()

        self._transect_based_len_age_abundance_report()

        self._transect_based_len_age_biomass_report()
