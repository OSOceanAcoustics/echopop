from typing import Optional

from pydantic import Field, model_validator
from typing_extensions import Self

from ..core.validators import BaseDictionary


class ValidateSelectivityParams(BaseDictionary):
    intercept: Optional[float] = Field(default=None, allow_inf_nan=False)
    slope: Optional[float] = Field(default=None, allow_inf_nan=False)
    l50: Optional[float] = Field(default=None, allow_inf_nan=False, gt=0.0)
    sr: Optional[float] = Field(default=None, allow_inf_nan=False, gt=0.0)
    minimum_selectivity: float = Field(default=1e-12, allow_inf_nan=False)

    @model_validator(mode="after")
    def validate_parameter_sets(self) -> Self:
        # Define the two valid sets
        regression_set = {self.intercept, self.slope}
        metrics_set = {self.l50, self.sr}

        # Check if a set is "complete" (no Nones)
        has_regression = all(x is not None for x in regression_set)
        has_metrics = all(x is not None for x in metrics_set)

        # Check if both are provided
        if has_regression and has_metrics:
            raise ValueError(
                "Ambiguous parameters: Provide either 'intercept'/'slope' or 'l50'/'sr'."
            )

        # Check if neither is complete
        if not has_regression and not has_metrics:
            raise ValueError(
                "Missing parameters: Must provide a complete set of 'intercept'/'slope' or "
                "'l50'/'sr'."
            )

        # Prevent "mixing and matching" (e.g., providing slope and l50)
        # ---- This catches if they provided one part of one set and one part of the other
        provided_fields = {
            k for k, v in self.model_dump().items() if v is not None and k != "minimum_selectivity"
        }

        valid_regression = {"intercept", "slope"}
        valid_metrics = {"l50", "sr"}

        # Raise Error
        if provided_fields != valid_regression and provided_fields != valid_metrics:
            # ---- Format
            fields_msg = "', '".join(provided_fields)
            raise ValueError(f"Invalid parameter combination: '{fields_msg}'.")

        return self
