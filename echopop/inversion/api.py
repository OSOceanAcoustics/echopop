from __future__ import annotations
from ..validators.scattering import ValidatePCDWBAParams, ValidatePCDWBASettings
from .scattering_models import pcdwba


SCATTERING_MODEL_PARAMETERS = {
    "pcdwba": {
        "function": pcdwba,
        "parameters": ValidatePCDWBAParams,
        "settings": ValidatePCDWBASettings,
    }
}