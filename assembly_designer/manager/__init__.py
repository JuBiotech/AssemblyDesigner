# assembly_designer/manager/__init__.py
"""Public exports for the manager package."""
from assembly_designer.utils import fn_plot, fn_plot2, plot_gif  # optional

from .manager import AssemblyWorkflowManager

__all__ = ["AssemblyWorkflowManager", "fn_plot", "fn_plot2", "plot_gif"]
