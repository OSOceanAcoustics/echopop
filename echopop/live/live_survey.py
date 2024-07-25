from typing import Union
from pathlib import Path
import copy

from .live_core import(
    LIVE_DATA_STRUCTURE,
)

from ..acoustics import (
    ts_length_regression,
    to_dB,
    to_linear
)

from . import live_data_processing as eldp
from . import live_data_loading as eldl
class LiveSurvey:
    """
    A real-time processing version of the `echopop` base `Survey` class that ingests biological, 
    acoustic, and event meta data to provide population estimates when generated.
    """

    def __init__(
        self,
        live_init_config_path: Union[str, Path], 
        live_file_config_path: Union[str, Path],
        verbose: bool = True,
    ):
        # Initialize `meta` attribute
        self.meta = copy.deepcopy(LIVE_DATA_STRUCTURE["meta"])

        # Loading the configuration settings and definitions that are used to
        # initialize the Survey class object
        self.config = eldl.live_configuration(Path(live_init_config_path), 
                                              Path(live_file_config_path))
        # ---- Initialize config key for database files
        self.config.update(
            {"database": {key: None for key in self.config["input_directories"].keys()}}
        )
        
        # Initialize input attribute
        self.input = copy.deepcopy(LIVE_DATA_STRUCTURE["input"])

        # Initialize database attribute
        self.database = copy.deepcopy(LIVE_DATA_STRUCTURE["database"])

        # Initialize the results attribute
        self.results = copy.deepcopy(LIVE_DATA_STRUCTURE["results"])

        # TODO: Replace Tuple output by appending the "database" key to the respective dataset dict
        # Ingest data
        # ---- Acoustics
        self.input["acoustics"]["prc_nasc_df"] = eldl.load_acoustic_data(self.config)
        # ---- Biology
        self.input["biology"] = eldp.load_biology_data(self.config)
        
        # TODO: Add verbosity for printing database filepaths/connections 
        if verbose: 
            pass