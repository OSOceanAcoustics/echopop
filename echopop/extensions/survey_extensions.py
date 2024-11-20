from ..survey import Survey
from .feat_report import FEATReports
import functools

####################################################################################################
# PATCH METHODS
# --------------------------------------------------------------------------------------------------

def _generate_reports(self: Survey, **kwargs) -> None:

    # Create `FEATReports` insance
    reports_instance = FEATReports(self, **kwargs)
    
    # Return the result for `FEATReports` class methods
    return reports_instance.generate()

####################################################################################################
# PATCHERS
# --------------------------------------------------------------------------------------------------
def patch_generate_reports():
    """
    Patch `generate_reports` to `Survey`
    """
    
    # Copy all hidden attributes
    functools.update_wrapper(_generate_reports, FEATReports.generate)
    
    # Assign method
    Survey.generate_reports = _generate_reports

