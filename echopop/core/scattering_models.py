from __future__ import annotations

from ..inversion.scattering_models import pcdwba
from ..validators.scattering import ValidatePCDWBAParams, ValidatePCDWBASettings

SCATTERING_MODEL_PARAMETERS = {
    "pcdwba": {
        "function": pcdwba,
        "parameters": ValidatePCDWBAParams,
        "settings": ValidatePCDWBASettings,
    }
}
