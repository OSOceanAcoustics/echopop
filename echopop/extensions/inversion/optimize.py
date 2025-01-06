import numpy as np
from lmfit import Minimizer, Parameters
from typing import Any, Callable, Dict, Literal, Union
from numpy.typing import ArrayLike

def mad(
    prediction: ArrayLike[float],
    measurement: ArrayLike[float],
):
    """
    Mean absolute deviation (MAD) in logarithmic space (dB)
    """
    # == functions/cost_functionALL.m
    pass

def rmse(
    prediction: ArrayLike[float],
    measurement: ArrayLike[float],
):
    """
    Root mean square deviation (RMSE) in logarithmic space (dB)
    """
    # == functions/cost_functionALL.m
    pass


def scattering_model_optimizer(
    prediction: ArrayLike[float],
    measurement: ArrayLike[float],
    parameters: Parameters,
    cost_function: Callable = mad,
):
    """
    Optimize scattering model parameters
    """
    # == functions/SVpredictionALL.m
    # == KrillSvInversion_simu_data_2020_05_01.m
    pass
