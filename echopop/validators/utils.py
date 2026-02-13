from typing import Any, Union

from pydantic import ConfigDict, Field

from ..core.validators import BaseDictionary


class ShipID(BaseDictionary):
    US: Any = Field(default=1)
    CAN: Any = Field(default=2)


class SurveyID(BaseDictionary):
    US: Any = Field(default=1)
    CAN: Any = Field(default=2)


class ValidateHaulUID(BaseDictionary):
    ship_id: ShipID = Field(default_factory=ShipID)
    survey_id: SurveyID = Field(default_factory=SurveyID)
    species_id: Any = Field(default=99999)
    haul_offset: Union[int, float] = Field(default=0)
    single_country: bool = Field(default=False)
    model_config = ConfigDict(title="haul-based identifier")
