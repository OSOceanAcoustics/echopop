import functools
from pathlib import Path
from typing import List, Union

from ..survey import Survey
from .diagnostics import DiagnosticPlot
from .feat_report import FEATReports
from .inversion.invert import AcousticInversion

####################################################################################################
# PATCH METHODS
####################################################################################################

# --------------------------------------------------------------------------------------------------
# .feat_report
# --------------------------------------------------------------------------------------------------


def _generate_reports(self: Survey, **kwargs) -> None:

    # Create `FEATReports` instance
    reports_instance = FEATReports(self, **kwargs)

    # Return the result for `FEATReports` class methods
    return reports_instance.generate()


def _report_options(self) -> List[str]:

    # Get report options
    return FEATReports.report_options()


# --------------------------------------------------------------------------------------------------
# .diagnostics
# --------------------------------------------------------------------------------------------------


def _mesh_cropping_results(self: Survey) -> None:

    # Plot cropped mesh results
    return DiagnosticPlot(self).mesh_cropping_results()


def _mesh_regions(self: Survey) -> None:

    # Plot mesh transect regions
    return DiagnosticPlot(self).mesh_regions()


def _nasc_map(self: Survey) -> None:

    # Plot high-leverage NASC data
    return DiagnosticPlot(self).nasc_map()


def _stratified_results(self: Survey) -> None:

    # Plot stratified results
    return DiagnosticPlot(self).stratified_results()


# --------------------------------------------------------------------------------------------------
# .inversion
# --------------------------------------------------------------------------------------------------


def _initialize_inversion(A_cls):

    @classmethod
    def _patch_initialize_inversion(cls, filepath) -> None:

        # Create `Inversion` instance
        instance = cls.__new__(cls)

        # Initialize
        AcousticInversion.__init__(instance, filepath)

        # Return
        return instance

    # Patch the method
    A_cls.initialize_inversion = _patch_initialize_inversion


def _invert_population(self: Survey, **kwargs) -> None:

    # Run inversion
    return AcousticInversion.invert_population(self, **kwargs)


def _fit_inversion_variogram(self: Survey, **kwargs) -> None:

    # Run variogram fitting
    return AcousticInversion.fit_inversion_variogram(self, **kwargs)


def _inversion_kriging_analysis(self: Survey, **kwargs) -> None:

    # Run kriging interpolation
    return AcousticInversion.inversion_kriging_analysis(self, **kwargs)


####################################################################################################
# PATCHERS
####################################################################################################

# --------------------------------------------------------------------------------------------------
# .feat_report
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


# --------------------------------------------------------------------------------------------------
# .diagnostics
# --------------------------------------------------------------------------------------------------


def patch_diagnostic_plots():
    """
    Patch diagnostic plotting functions to `Survey`
    """

    # Copy all hidden attributes
    # ---- `mesh_cropping_results()`
    functools.update_wrapper(_mesh_cropping_results, DiagnosticPlot.mesh_cropping_results)
    # ---- `mesh_regions()`
    functools.update_wrapper(_mesh_regions, DiagnosticPlot.mesh_regions)
    # ---- `nasc_map()`
    functools.update_wrapper(_nasc_map, DiagnosticPlot.nasc_map)
    # ---- `stratified_results()`
    functools.update_wrapper(_stratified_results, DiagnosticPlot.stratified_results)

    # Assign method
    # ---- `generate()`
    Survey.mesh_cropping_results = _mesh_cropping_results
    # ---- `mesh_regions()`
    Survey.mesh_regions = _mesh_regions
    # ---- `nasc_map()`
    Survey.nasc_map = _nasc_map
    # ---- `stratified_results()`
    Survey.stratified_results = _stratified_results


# --------------------------------------------------------------------------------------------------
# .inversion
# --------------------------------------------------------------------------------------------------


def patch_invert_population():
    """
    Patch inversion placeholder method to `Survey`
    """

    # Copy all hidden attributes
    # ---- `initialize_inversion()`
    functools.update_wrapper(_initialize_inversion, AcousticInversion)
    # ---- `invert_population()`
    functools.update_wrapper(_invert_population, AcousticInversion.invert_population)
    # ---- `fit_inversion_variogram()`
    functools.update_wrapper(_fit_inversion_variogram, AcousticInversion.fit_inversion_variogram)
    # ---- `inversion_kriging_analysis()`
    functools.update_wrapper(
        _inversion_kriging_analysis, AcousticInversion.inversion_kriging_analysis
    )

    # Assign method
    # ---- `initialize_inversion()`
    _initialize_inversion(Survey)
    # ---- `invert_population()`
    Survey.invert_population = _invert_population
    # ---- `fit_inversion_variogram()`
    Survey.fit_inversion_variogram = _fit_inversion_variogram
    # ---- `inversion_kriging_analysis()`
    Survey.inversion_kriging_analysis = _inversion_kriging_analysis
