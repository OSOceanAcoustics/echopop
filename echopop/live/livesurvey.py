from typing import Union
from pathlib import Path
import copy
import yaml

from .core import(
    DATA_STRUCTURE
)

from ..acoustics import (
    ts_length_regression,
    to_dB,
    to_linear
)

class LiveSurvey:
    """
    A real-time processing version of the `echopop` base 
    `Survey` class that ingests biological, acoustic, and
    event meta data to provide population estimates when 
    generated.
    """

    def __init__(
        self
    ):
        # Initialize `meta` attribute
        self.meta = copy.deepcopy(DATA_STRUCTURE["meta"])

        # Loading the configuration settings and definitions that are used to
        # initialize the Survey class object
        self.config = el.load_configuration(Path(init_config_path), Path(survey_year_config_path))

        # Loading the datasets defined in the configuration files
        self.input = el.load_survey_data(self.config)

        # Initialize the `analysis` data attribute
        self.analysis = copy.deepcopy(DATA_STRUCTURE["analysis"])

        # Initialize the `results` data attribute
        self.results = copy.deepcopy(DATA_STRUCTURE["results"])