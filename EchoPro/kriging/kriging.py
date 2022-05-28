import numpy as np
import pandas as pd
import warnings
import sys
import folium
from matplotlib.colors import to_hex


class Kriging:
    """
    This class constructs all data necessary
    for kriging and performs kriging

    Parameters
    ----------
    EPro : EchoPro object
        An initialized EchoPro object. Note that any change to
        self.EPro will also change this object.
    """

    def __init__(self, EPro = None):

        self.EPro = EPro
