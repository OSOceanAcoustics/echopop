"""Inherited validators for model-specific validators."""

from ..inversion.pcdwba import pcdwba
from .scattering import ValidatePCDWBAParams, ValidatePCDWBASettings

SCATTERING_MODEL_PARAMETERS = {
    "pcdwba": {
        "function": pcdwba,
        "parameters": ValidatePCDWBAParams,
        "settings": ValidatePCDWBASettings,
    }
}
