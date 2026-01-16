from .inversion import (
    ValidateBuildModelArgs,
    ValidateInversionMatrix,
    ValidateLengthTS,
)
from .kriging import ValidateKrigingClass
from .spatial import ValidateHullCropArgs
from .variogram import (
    ValidateEmpiricalVariogramArgs,
    ValidateFitVariogramArgs,
    ValidateVariogramClass,
)


# Rebuild ValidateBuildModelArgs after InvParameters is available
def _rebuild_on_import():
    try:
        from ..inversion import InvParameters  # noqa: F401

        ValidateBuildModelArgs.model_rebuild()
    except ImportError:
        pass


_rebuild_on_import()

__all__ = [
    "ValidateBuildModelArgs",
    "ValidateHullCropArgs",
    "ValidateEmpiricalVariogramArgs",
    "ValidateFitVariogramArgs",
    "ValidateInversionMatrix",
    "ValidateKrigingClass",
    "ValidateLengthTS",
    "ValidateVariogramClass",
]
