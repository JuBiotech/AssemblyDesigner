"""Public exports for the `assembly_designer.utils` subpackage."""

from __future__ import annotations

# Plotting / visualization
# Mastermix calculators (new snake_case API)
from .utils import (
    FZcmaps,
    FZcolors,
    calculate_dreamtaq_mastermix_components,
    calculate_gg_mastermix_components,
    calculate_pcr_mastermix_components,
    check_and_log_volumes,
    find_planning_file,
    fn_plot,
    fn_plot2,
    mtpshow,
    plot_gif,
    setup_logging,
    # --- NEW: logging helpers (file logging + DF volume checks)
    setup_run_logging_files,
    to_colormap,
    transparentify,
    values_in_2d,
)

__all__ = [
    # core
    "setup_logging",
    "setup_run_logging_files",
    "check_and_log_volumes",
    "find_planning_file",
    # plotting
    "fn_plot",
    "fn_plot2",
    "plot_gif",
    "mtpshow",
    "values_in_2d",
    "to_colormap",
    "transparentify",
    "FZcolors",
    "FZcmaps",
    # mastermix calculators (snake_case API)
    "calculate_pcr_mastermix_components",
    "calculate_gg_mastermix_components",
    "calculate_dreamtaq_mastermix_components",
    # legacy (camelCase) names are exported lazily via __getattr__ below
    "calculate_GGmastermix_components",
    "calculate_DreamTaqMastermix_components",
    "calculate_GibsonMastermix_components",
]


def __getattr__(name: str):
    """
    Lazy compatibility re-exports for legacy names without importing at module import time.
    Avoids Ruff F401 ('imported but unused') in this __init__.
    """
    if name in {
        "calculate_GGmastermix_components",
        "calculate_DreamTaqMastermix_components",
        "calculate_GibsonMastermix_components",
    }:
        from . import utils as _u

        return getattr(_u, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
