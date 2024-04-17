import numpy as np


def ts_length_regression(x, slope, intercept):
    return slope * np.log10(x) + intercept


def to_linear(x):
    return 10.0 ** (x / 10.0)


def to_dB(x):
    return 10 * np.log10(x)
