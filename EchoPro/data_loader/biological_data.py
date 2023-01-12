from pathlib import Path

import numpy as np
import pandas as pd

from ..utils.input_checks import check_column_names, check_existence_of_file


class LoadBioData:  # TODO: Does it make sense for this to be a class?
    """
    This class loads and checks those files
    associated with the biological data.

    Parameters
    ----------
    survey : Survey
        An initialized Survey object. Note that any change to
        self.survey will also change this object.
    """

    def __init__(self, survey=None):

        self.survey = survey

        # expected columns for length Dataframe
        self.len_cols = {"haul_num", "species_id", "sex", "length", "length_count"}

        # expected columns for specimen Dataframe
        self.spec_cols = {"haul_num", "species_id", "sex", "length", "weight", "age"}

        # expected columns for catch Dataframe
        self.catch_cols = {"haul_num", "species_id", "haul_count", "haul_weight"}

        # expected columns for haul_to_transect_mapping Dataframe
        self.haul_to_transect_mapping_cols = {"haul_num", "transect_num"}

        self._load_length_data()
        self._load_specimen_data()
        self._load_catch_data()
        self._load_haul_to_transect_mapping_data()

    def _check_length_df(self, len_df: pd.DataFrame, df_path: Path) -> None:
        """
        Ensures that the appropriate columns are
        contained in the length Dataframe.

        Parameters
        ----------
        len_df: pd.DataFrame
            The constructed Length DataFrame
        df_path: Path
            The path to the Excel file used to construct the DataFrame
        """

        # TODO: should we add more in-depth checks here?

        check_column_names(df=len_df, expected_names=self.len_cols, path_for_df=df_path)

    def _check_specimen_df(self, spec_df: pd.DataFrame, df_path: Path) -> None:
        """
        Ensures that the appropriate columns are
        contained in the specimen Dataframe.

        Parameters
        ----------
        spec_df: pd.DataFrame
            The constructed Specimen DataFrame
        df_path: Path
            The path to the Excel file used to construct the DataFrame
        """

        # TODO: should we add more in-depth checks here?

        check_column_names(
            df=spec_df, expected_names=self.spec_cols, path_for_df=df_path
        )

    def _check_catch_df(self, catch_df: pd.DataFrame, df_path: Path) -> None:
        """
        Ensures that the appropriate columns are
        contained in the Catch Dataframe.

        Parameters
        ----------
        catch_df: pd.DataFrame
            The constructed Catch DataFrame
        df_path: Path
            The path to the Excel file used to construct the DataFrame
        """

        # TODO: should we add more in-depth checks here?

        check_column_names(
            df=catch_df, expected_names=self.catch_cols, path_for_df=df_path
        )

    def _check_haul_to_transect_mapping_df(
        self, haul_to_transect_mapping_df: pd.DataFrame, df_path: Path
    ) -> None:
        """
        Ensures that the appropriate columns are
        contained in the haul to transect mapping Dataframe.

        Parameters
        ----------
        haul_to_transect_mapping_df: pd.DataFrame
            The constructed haul to transect mapping DataFrame
        df_path: Path
            The path to the Excel file used to construct the DataFrame
        """

        # TODO: should we add more in-depth checks here?

        check_column_names(
            df=haul_to_transect_mapping_df,
            expected_names=self.haul_to_transect_mapping_cols,
            path_for_df=df_path,
        )

    def _process_length_data_df(
        self, df: pd.DataFrame, haul_num_offset: int
    ) -> pd.DataFrame:
        """
        Processes the length dataframe by:
        * Obtaining the required columns from the dataframe
        * Extracting only the target species
        * Setting the data type of each column
        * Applying a haul offset, if necessary
        * Setting the index required for downstream processes

        Parameters
        ----------
        df : Pandas Dataframe
            Dataframe holding the length data
        haul_num_offset : int
            The offset that should be applied to the ``haul_num`` column

        Returns
        -------
        Processed Dataframe
        """

        # obtaining those columns that are required
        df = df[["haul_num", "species_id", "sex", "length", "length_count"]].copy()

        # set data types of dataframe
        df = df.astype(
            {
                "haul_num": int,
                "species_id": int,
                "sex": int,
                "length": np.float64,
                "length_count": np.float64,
            }
        )

        # extract target species
        df = df.loc[df["species_id"] == self.survey.params["species_id"]]

        # Apply haul offset
        df["haul_num"] = df["haul_num"] + haul_num_offset

        if self.survey.params["exclude_age1"] is False:
            raise NotImplementedError("Including age 1 data has not been implemented!")

        # remove species_id column
        df.drop(columns=["species_id"], inplace=True)

        df.set_index("haul_num", inplace=True)

        return df

    def _process_specimen_data(
        self, df: pd.DataFrame, haul_num_offset: int
    ) -> pd.DataFrame:
        """
        Processes the specimen dataframe by:
        * Obtaining the required columns from the dataframe
        * Extracting only the target species
        * Setting the data type of each column
        * Applying a haul offset, if necessary
        * Setting the index required for downstream processes

        Parameters
        ----------
        df : Pandas Dataframe
            Dataframe holding the specimen data
        haul_num_offset : int
            The offset that should be applied to the ``haul_num`` column

        Returns
        -------
        Processed Dataframe
        """

        # obtaining those columns that are required
        df = df[["haul_num", "species_id", "sex", "length", "weight", "age"]].copy()

        # set data types of dataframe
        df = df.astype(
            {
                "haul_num": int,
                "species_id": int,
                "sex": int,
                "length": np.float64,
                "weight": np.float64,
                "age": np.float64,
            }
        )

        # extract target species
        df = df.loc[df["species_id"] == self.survey.params["species_id"]]

        # Apply haul_num_offset
        df["haul_num"] = df["haul_num"] + haul_num_offset

        if self.survey.params["exclude_age1"] is False:
            raise NotImplementedError("Including age 1 data has not been implemented!")

        # remove species_id column
        df.drop(columns=["species_id"], inplace=True)

        # set and organize index
        df.set_index("haul_num", inplace=True)
        df.sort_index(inplace=True)

        # perform check on data
        if len(df["age"]) - df["age"].isna().sum() < 0.1 * len(df["age"]):
            raise RuntimeWarning("Aged data are less than 10%!\n")

        return df

    def _process_catch_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Processes the catch dataframe by:
        * Obtaining the required columns from the dataframe
        * Extracting only the target species
        * Setting the data type of each column
        * Setting the index required for downstream processes

        Parameters
        ----------
        df : Pandas Dataframe
            Dataframe holding the catch data

        Returns
        -------
        Processed Dataframe
        """

        # obtaining those columns that are required
        df = df[["haul_num", "species_id", "haul_count", "haul_weight"]].copy()

        # set data types of dataframe
        df = df.astype(
            {
                "haul_num": int,
                "species_id": int,
                "haul_count": np.float64,
                "haul_weight": np.float64,
            }
        )

        # extract target species
        df = df.loc[df["species_id"] == self.survey.params["species_id"]]

        # remove species_id column
        df.drop(columns=["species_id"], inplace=True)

        df.set_index("haul_num", inplace=True)

        df.sort_index(inplace=True)

        return df

    def _load_length_data(self) -> None:
        """
        Loads and prepares data associated with a station
        that records the length and sex of the animal.
        Additionally, it sets survey.length_df using the
        final processed dataframe.
        """

        if self.survey.params["source"] == 3:

            # check existence of the files
            file_path_us = (
                self.survey.params["data_root_dir"]
                / self.survey.params["length_US_filename"]
            )
            file_path_can = (
                self.survey.params["data_root_dir"]
                / self.survey.params["length_CAN_filename"]
            )
            check_existence_of_file(file_path_us)
            check_existence_of_file(file_path_can)

            # read in and check US and Canada Excel files
            df_us = pd.read_excel(
                file_path_us, sheet_name=self.survey.params["length_US_sheet"]
            )
            self._check_length_df(df_us, file_path_us)

            df_can = pd.read_excel(
                file_path_can, sheet_name=self.survey.params["length_CAN_sheet"]
            )
            self._check_length_df(df_can, file_path_can)

            # process US and Canada dataframes
            length_us_df = self._process_length_data_df(df_us, 0)
            length_can_df = self._process_length_data_df(
                df_can, self.survey.params["CAN_haul_offset"]
            )

            # Construct full length dataframe from US and Canada sections
            self.survey.length_df = pd.concat([length_us_df, length_can_df])

        else:
            raise NotImplementedError(
                f"Source of {self.survey.params['source']} not implemented yet."
            )

    def _load_specimen_data(self) -> None:
        """
        Loads and prepares data associated with a station
        that records the length, weight, age, and sex of
        the animal. Additionally, it sets survey.specimen_df
        using the final processed dataframe.
        """

        if self.survey.params["source"] == 3:

            # check existence of the files
            file_path_us = (
                self.survey.params["data_root_dir"]
                / self.survey.params["specimen_US_filename"]
            )
            file_path_can = (
                self.survey.params["data_root_dir"]
                / self.survey.params["specimen_CAN_filename"]
            )
            check_existence_of_file(file_path_us)
            check_existence_of_file(file_path_can)

            # read in and check US and Canada Excel files
            specimen_us_df = pd.read_excel(
                file_path_us, sheet_name=self.survey.params["specimen_US_sheet"]
            )
            self._check_specimen_df(specimen_us_df, file_path_us)

            specimen_can_df = pd.read_excel(
                file_path_can, sheet_name=self.survey.params["specimen_CAN_sheet"]
            )
            self._check_specimen_df(specimen_can_df, file_path_can)

            # process US and Canada dataframes
            specimen_us_df = self._process_specimen_data(specimen_us_df, 0)
            specimen_can_df = self._process_specimen_data(
                specimen_can_df, self.survey.params["CAN_haul_offset"]
            )

            # Construct full specimen dataframe from US and Canada sections
            self.survey.specimen_df = pd.concat([specimen_us_df, specimen_can_df])

        else:
            raise NotImplementedError(
                f"Source of {self.survey.params['source']} not implemented yet."
            )

    def _load_catch_data(self):
        """
        Loads and prepares data associated with a station
        that records the catch data of a particular haul.
        Additionally, it sets survey.catch_df using the
        final processed dataframe.
        """

        if self.survey.params["source"] == 3:

            # check existence of the files
            file_path_us = (
                self.survey.params["data_root_dir"]
                / self.survey.params["catch_US_filename"]
            )
            file_path_can = (
                self.survey.params["data_root_dir"]
                / self.survey.params["catch_CAN_filename"]
            )
            check_existence_of_file(file_path_us)
            check_existence_of_file(file_path_can)

            # read in and check US and Canada Excel files
            catch_us_df = pd.read_excel(
                file_path_us, sheet_name=self.survey.params["catch_US_sheet"]
            )
            self._check_catch_df(catch_us_df, file_path_us)

            catch_can_df = pd.read_excel(
                file_path_can, sheet_name=self.survey.params["catch_CAN_sheet"]
            )
            self._check_catch_df(catch_can_df, file_path_can)

            # process US and Canada dataframes
            catch_us_df = self._process_catch_data(catch_us_df)
            catch_can_df = self._process_catch_data(catch_can_df)

            # apply haul offset to Canada data
            catch_can_df.index = (
                catch_can_df.index + self.survey.params["CAN_haul_offset"]
            )

            # Construct full catch dataframe from US and Canada sections
            self.survey.catch_df = pd.concat([catch_us_df, catch_can_df])

        else:

            raise NotImplementedError(
                f"Source of {self.survey.params['source']} not implemented yet."
            )

    def _process_haul_to_transect_mapping_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Processes the haul to transect mapping data by
        * selecting the haul and transect columns
        * ensuring the dataframe has the appropriate data types
        * setting the ``haul_num`` column as the Dataframe index
        * sorting the ``haul_num`` index in ascending order

        Parameters
        ----------
        df : pd. Dataframe
            Dataframe holding the haul to transect mapping data

        Returns
        -------
        pd.DataFrame
            Processed haul to transect mapping Dataframe
        """

        # obtain those columns necessary for core downstream processes
        df = df[["haul_num", "transect_num"]].copy()

        # set data types of dataframe
        df = df.astype({"haul_num": int, "transect_num": np.float64})

        if self.survey.params["exclude_age1"] is False:
            raise NotImplementedError("Including age 1 data has not been implemented!")

        # set haul_num as index and sort it
        df.set_index("haul_num", inplace=True)
        df.sort_index(inplace=True)

        return df

    def _load_haul_to_transect_mapping_data(self) -> None:
        """
        Loads and prepares the data that maps each ``haul_num`` to a ``transect_num``.
        Additionally, it sets survey.haul_to_transect_mapping_df using the final
        processed dataframe.

        Notes
        -----
        This data is necessary for obtaining subsets of data during transect selection.
        """

        if self.survey.params["source"] == 3:

            # check existence of the files
            file_path_us = (
                self.survey.params["data_root_dir"]
                / self.survey.params["filename_haul_to_transect_US"]
            )
            file_path_can = (
                self.survey.params["data_root_dir"]
                / self.survey.params["filename_haul_to_transect_CAN"]
            )
            check_existence_of_file(file_path_us)
            check_existence_of_file(file_path_can)

            # read in, check, and process the US haul to transect mapping file
            haul_to_transect_mapping_us_df = pd.read_excel(
                file_path_us,
                sheet_name=self.survey.params["haul_to_transect_US_sheetname"],
            )
            self._check_haul_to_transect_mapping_df(
                haul_to_transect_mapping_us_df, file_path_us
            )
            haul_to_transect_mapping_us_df = (
                self._process_haul_to_transect_mapping_data(
                    haul_to_transect_mapping_us_df
                )
            )

            # read in, check, and process the Canada haul to transect mapping file
            haul_to_transect_mapping_can_df = pd.read_excel(
                file_path_can,
                sheet_name=self.survey.params["haul_to_transect_CAN_sheetname"],
            )
            self._check_haul_to_transect_mapping_df(
                haul_to_transect_mapping_can_df, file_path_can
            )
            haul_to_transect_mapping_can_df = (
                self._process_haul_to_transect_mapping_data(
                    haul_to_transect_mapping_can_df
                )
            )

            # transect offset is set to zero since we will not have overlapping transects
            # TODO: Is this a necessary variable? If so, we should make it an input to the function.
            CAN_Transect_offset = 0

            # add Canada transect and haul offset
            haul_to_transect_mapping_can_df.index = (
                haul_to_transect_mapping_can_df.index
                + self.survey.params["CAN_haul_offset"]
            )
            haul_to_transect_mapping_can_df["transect_num"] = (
                haul_to_transect_mapping_can_df["transect_num"] + CAN_Transect_offset
            )

            # combine US & CAN data
            self.survey.haul_to_transect_mapping_df = pd.concat(
                [haul_to_transect_mapping_us_df, haul_to_transect_mapping_can_df]
            )

        else:
            raise NotImplementedError(
                f"Source of {self.survey.params['source']} not implemented yet."
            )
