import functools
from typing import List

from ..survey import Survey
from .feat_report import FEATReports

####################################################################################################
# PATCH METHODS
# --------------------------------------------------------------------------------------------------


def _generate_reports(self: Survey, **kwargs) -> None:

    # Create `FEATReports` instance
    reports_instance = FEATReports(self, **kwargs)

    # Return the result for `FEATReports` class methods
    return reports_instance.generate()


def _report_options(self) -> List[str]:

    # Get report options
    return FEATReports.report_options()


####################################################################################################
# PATCHERS
# --------------------------------------------------------------------------------------------------
def patch_generate_reports():
    """
    Patch `generate_reports` to `Survey`
    """

    # Copy all hidden attributes
    # ---- `generate()`
    functools.update_wrapper(_generate_reports, FEATReports)
    # ---- `report_options()`
    functools.update_wrapper(_report_options, FEATReports.report_options)

    # Assign method
    # ---- `generate()`
    Survey.generate_reports = _generate_reports
    # ---- `report_options()`
    Survey.report_options = _report_options
