from .survey_extensions import (
    patch_diagnostic_plots as diagnostic_plots,
    patch_generate_reports as generate_reports,
)


# Flex patch import
def import_all_patches():
    diagnostic_plots()
    generate_reports()


# ---- Automatic import for `import echopop.extensions`
import_all_patches()

# Exposed individual patches for selective importing
__all__ = ["diagnostic_plots", "generate_reports"]
