"""Pydantic validators for workflow-based UID definitions and other utility functions."""

from typing import Any

from pydantic import ConfigDict, Field

from ..core.validators import BaseDictionary


class ShipID(BaseDictionary):
    """Country-keyed ship identifier mapping."""

    US: Any = Field(default=1)
    CAN: Any = Field(default=2)


class SurveyID(BaseDictionary):
    """Country-keyed survey identifier mapping."""

    US: Any = Field(default=1)
    CAN: Any = Field(default=2)


class ValidateHaulUID(BaseDictionary):
    """Validation model for haul unique-identifier construction parameters."""

    ship_id: ShipID = Field(default_factory=ShipID)
    survey_id: SurveyID = Field(default_factory=SurveyID)
    species_id: Any = Field(default=99999)
    haul_offset: int | float = Field(default=0)
    single_country: bool = Field(default=False)
    model_config = ConfigDict(title="haul-based identifier")
