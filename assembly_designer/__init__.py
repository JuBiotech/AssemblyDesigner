# assembly_designer/__init__.py
"""
Assembly Designer — in-silico plasmid construction & liquid-handling simulation.

This package exposes two pillars under one roof:

1) Liquid-handling design & simulation
   - Build and validate automated workflows
   - Generate worklists for robots (Opentrons / Tecan)
   - Plot and inspect steps

2) In-silico plasmid design (“Plasmidio”)
   - Run DnaCauldron assemblies (combinatorial & recipe-based)
   - Extract and clean GenBank outputs
   - Simulate restriction digests & gels
   - Align reads, visualize, and export reports
   - PCR/Gibson helpers and build provenance (DAGs)

What this module does
---------------------
This top-level ``__init__`` provides **lazy re-exports** of the most used API so
you can write::

    import assembly_designer as ad
    from assembly_designer import LiquidHandlingWorkflow
    from assembly_designer import generate_assembly_reports  # (Plasmidio shortcut)

without importing all submodules up front. Imports are resolved on first use.

Submodule namespace
-------------------
- The full in-silico toolkit is available as
  ``assembly_designer.plasmidio``.
- For backward compatibility, ``assembly_designer.insilico`` is an alias of
  the same namespace.

Plasmidio convenience re-exports (shortcuts)
--------------------------------------------
- **Reports:**
  - ``generate_assembly_reports``  (combinatorial/cartesian builder)
  - ``generate_recipe_assemblies`` (explicit recipe builder)
  - ``organize_assembly_reports``  (unpack/arrange and optionally delete ZIPs)
  - ``run_all_functions``          (convenience pipeline wrapper)

- **I/O & features:**
  - ``snapgene_to_seqrecord``, ``read_plasmid_file``, ``load_dna_file``
  - ``filter_features``, ``remove_near_duplicate_features``
  - ``load_parts_from_folders``

- **Digest & gel:**
  - ``perform_restriction_digest``, ``get_ladder``, ``plot_gel``,
    ``digest_and_plot_plasmids``

- **Batch & export:**
  - ``prepare_plasmids``, ``load_reads``, ``align_batch``,
    ``results_to_dataframe``, ``export_results``, ``export_feature_inventory``
  - **Parallel alignment:** ``align_batch_parallel`` (threads/processes; optional)

- **Visualization:**
  - ``plot_alignment``, ``plot_alignment_by_ref``
  - ``view_alignment*`` (text/HTML helpers)

- **PCR/Gibson helpers:**
  - ``PCRResult``, ``simulate_pcr``, ``simulate_pcr_overhangs``,
    ``merge_with_gibson_features``, ``longest_overlap``

- **Build naming & provenance:**
  - ``build_construct_label``
  - ``AssemblyHistory``, ``HistoryNode``, ``HistoryEdge``,
    ``build_histories_for_all_constructs``

- **3G batch (Golden Gate → PCR → Gibson):**
  - ``run_3g_batch_safe``, ``BuildStatus``

Core (liquid-handling) API
--------------------------
- ``LiquidHandlingWorkflow``
- ``AssemblyWorkflowManager``
- ``MTPManager``

Backward compatibility
----------------------
- ``generate_combinatorial_assemblies`` → alias of ``generate_assembly_reports``.
- ``insilico`` → alias of the ``plasmidio`` namespace.

Notes
-----
- The actual implementations live in their respective modules (see
  ``assembly_designer.plasmidio`` for in-silico utilities and
  ``assembly_designer.workflow`` / ``assembly_designer.manager`` / ``assembly_designer.mtp_manager`` for
  liquid-handling).
- This top-level file focuses on keeping imports **lightweight** and the public
  API **stable** for users and notebooks.
"""

from __future__ import annotations

from importlib import import_module
from typing import Any

# --------------------------------------------------------------------------- #
# Version (from installed distribution; safe fallback for editable installs)  #
# --------------------------------------------------------------------------- #
try:
    from importlib.metadata import version as _version  # Python 3.8+

    __version__ = _version("assembly-designer")
except Exception:  # pragma: no cover
    __version__ = "0.0.0"

# --------------------------------------------------------------------------- #
# Public API                                                                  #
# --------------------------------------------------------------------------- #
__all__ = [
    # Core classes
    "LiquidHandlingWorkflow",
    "AssemblyWorkflowManager",
    "MTPManager",
    # Utilities
    "find_planning_file",
    "mtpshow",
    "values_in_2d",
    # Plotting helpers
    "fn_plot",
    "fn_plot2",
    "plot_gif",
    # Mastermix calculators
    "calculate_pcr_mastermix_components",
    "calculate_gg_mastermix_components",
    "calculate_dreamtaq_mastermix_components",
    "calculate_gibson_mastermix_components",
    # Submodule namespaces
    "plasmidio",
    "insilico",
    # Plasmidio re-exports (reports / I/O / digestion)
    "generate_assembly_reports",  # combinatorial
    "generate_recipe_assemblies",  # explicit recipes
    "generate_combinatorial_assemblies",  # alias for backward-compat
    "organize_assembly_Reports",  # keep legacy-cased name
    "organize_assembly_reports",  # <-- NEW: canonical casing added (non-breaking)
    "run_all_functions",
    "snapgene_to_seqrecord",
    "read_plasmid_file",
    "load_dna_file",  # <-- NEW: was in docs, now exported
    "perform_restriction_digest",
    "plot_gel",
    "digest_and_plot_plasmids",
    # Plasmidio re-exports (batch + export)
    "prepare_plasmids",
    "load_reads",
    "align_batch",
    "align_batch_parallel",  # shortcut
    "results_to_dataframe",
    "summarize_matches" "find_unbuilt_plasmids" "export_results",
    "export_feature_inventory",
    # Plasmidio re-exports (viewers)
    "plot_alignment",
    "plot_alignment_by_ref",
    "view_alignment",
    "view_alignment_by_ref",
    "view_alignment_html",
    "view_alignment_html_by_ref",
    # Plasmidio re-exports (PCR/Gibson + naming)
    "PCRResult",
    "simulate_pcr",
    "simulate_pcr_overhangs",
    "merge_with_gibson_features",
    "longest_overlap",
    "build_construct_label",
    # NEW: convenience wrappers (optional in plasmidio.py)
    "dealign_reads_easy",  # <-- NEW optional shortcut
    "read_anchored_align",  # <-- NEW optional low-level aligner
    # NEW: 3G batch runner
    "run_3g_batch_safe",
    "BuildStatus",
    # NEW: History / provenance
    "AssemblyHistory",
    "HistoryNode",
    "HistoryEdge",
    "PCRResultLike",
    "build_histories_for_all_constructs",
    # Meta
    "__version__",
]

# Map attribute -> (module, symbol) for lazy loading
_PUBLIC_MAP: dict[str, tuple[str, str]] = {
    # Core
    "LiquidHandlingWorkflow": ("assembly_designer.workflow", "LiquidHandlingWorkflow"),
    "AssemblyWorkflowManager": ("assembly_designer.manager", "AssemblyWorkflowManager"),
    "MTPManager": ("assembly_designer.mtp_manager", "MTPManager"),
    # Plotting helpers
    "fn_plot": ("assembly_designer.utils", "fn_plot"),
    "fn_plot2": ("assembly_designer.utils", "fn_plot2"),
    "plot_gif": ("assembly_designer.utils", "plot_gif"),
    # Utilities
    "find_planning_file": ("assembly_designer.utils", "find_planning_file"),
    "mtpshow": ("assembly_designer.utils", "mtpshow"),
    "values_in_2d": ("assembly_designer.utils", "values_in_2d"),
    # Mastermix calculators
    "calculate_pcr_mastermix_components": (
        "assembly_designer.utils",
        "calculate_pcr_mastermix_components",
    ),
    "calculate_gg_mastermix_components": (
        "assembly_designer.utils",
        "calculate_gg_mastermix_components",
    ),
    "calculate_dreamtaq_mastermix_components": (
        "assembly_designer.utils",
        "calculate_dreamtaq_mastermix_components",
    ),
    "calculate_gibson_mastermix_components": (
        "assembly_designer.utils",
        "calculate_gibson_mastermix_components",
    ),
    # Plasmidio: reports / I/O / digestion
    "generate_assembly_reports": (
        "assembly_designer.plasmidio",
        "generate_assembly_reports",
    ),
    "generate_recipe_assemblies": (
        "assembly_designer.plasmidio",
        "generate_recipe_assemblies",
    ),
    # Back-compat alias → same target as generate_assembly_reports
    "generate_combinatorial_assemblies": (
        "assembly_designer.plasmidio",
        "generate_assembly_reports",
    ),
    # Keep both names mapped to the same target (non-breaking)
    "organize_assembly_Reports": (
        "assembly_designer.plasmidio",
        "organize_assembly_reports",
    ),
    "organize_assembly_reports": (
        "assembly_designer.plasmidio",
        "organize_assembly_reports",
    ),
    "run_all_functions": ("assembly_designer.plasmidio", "run_all_functions"),
    "snapgene_to_seqrecord": ("assembly_designer.plasmidio", "snapgene_to_seqrecord"),
    "read_plasmid_file": ("assembly_designer.plasmidio", "read_plasmid_file"),
    "load_dna_file": ("assembly_designer.plasmidio", "load_dna_file"),  # <-- NEW
    "perform_restriction_digest": (
        "assembly_designer.plasmidio",
        "perform_restriction_digest",
    ),
    "plot_gel": ("assembly_designer.plasmidio", "plot_gel"),
    "digest_and_plot_plasmids": (
        "assembly_designer.plasmidio",
        "digest_and_plot_plasmids",
    ),
    # Plasmidio: batch + export
    "prepare_plasmids": ("assembly_designer.plasmidio", "prepare_plasmids"),
    "load_reads": ("assembly_designer.plasmidio", "load_reads"),
    "align_batch": ("assembly_designer.plasmidio", "align_batch"),
    "align_batch_parallel": (
        "assembly_designer.plasmidio",
        "align_batch_parallel",
    ),
    "results_to_dataframe": ("assembly_designer.plasmidio", "results_to_dataframe"),
    "export_results": ("assembly_designer.plasmidio", "export_results"),
    "export_feature_inventory": (
        "assembly_designer.plasmidio",
        "export_feature_inventory",
    ),
    # Plasmidio: viewers
    "plot_alignment": ("assembly_designer.plasmidio", "plot_alignment"),
    "plot_alignment_by_ref": ("assembly_designer.plasmidio", "plot_alignment_by_ref"),
    "view_alignment": ("assembly_designer.plasmidio", "view_alignment"),
    "view_alignment_by_ref": ("assembly_designer.plasmidio", "view_alignment_by_ref"),
    "view_alignment_html": ("assembly_designer.plasmidio", "view_alignment_html"),
    "view_alignment_html_by_ref": (
        "assembly_designer.plasmidio",
        "view_alignment_html_by_ref",
    ),
    # Plasmidio: PCR/Gibson + naming
    "PCRResult": ("assembly_designer.plasmidio", "PCRResult"),
    "simulate_pcr": ("assembly_designer.plasmidio", "simulate_pcr"),
    "simulate_pcr_overhangs": ("assembly_designer.plasmidio", "simulate_pcr_overhangs"),
    "merge_with_gibson_features": (
        "assembly_designer.plasmidio",
        "merge_with_gibson_features",
    ),
    "longest_overlap": ("assembly_designer.plasmidio", "longest_overlap"),
    "build_construct_label": ("assembly_designer.plasmidio", "build_construct_label"),
    # NEW: convenience wrappers (optional in plasmidio.py)
    "dealign_reads_easy": (
        "assembly_designer.plasmidio",
        "dealign_reads_easy",
    ),  # <-- NEW
    "read_anchored_align": (
        "assembly_designer.plasmidio",
        "read_anchored_align",
    ),  # <-- NEW
    # NEW: 3G batch runner
    "run_3g_batch_safe": ("assembly_designer.plasmidio", "run_3g_batch_safe"),
    "BuildStatus": ("assembly_designer.plasmidio", "BuildStatus"),
    # NEW: History / provenance
    "AssemblyHistory": ("assembly_designer.plasmidio", "AssemblyHistory"),
    "HistoryNode": ("assembly_designer.plasmidio", "HistoryNode"),
    "HistoryEdge": ("assembly_designer.plasmidio", "HistoryEdge"),
    "PCRResultLike": ("assembly_designer.plasmidio", "PCRResultLike"),
    "build_histories_for_all_constructs": (
        "assembly_designer.plasmidio",
        "build_histories_for_all_constructs",
    ),
}


def __getattr__(name: str) -> Any:
    """Lazy attribute loader so `import assembly_designer` stays light-weight.

    - For names present in `_PUBLIC_MAP`, import the module and return the symbol.
    - For 'plasmidio' load and return the submodule namespace.
    - For 'insilico' (backward-compat alias), return the 'plasmidio' namespace.
    """
    if name in ("plasmidio", "insilico"):
        mod = import_module("assembly_designer.plasmidio")
        # Cache both names to the same module object (alias)
        globals()["plasmidio"] = mod
        globals()["insilico"] = mod
        return mod

    try:
        mod_name, attr = _PUBLIC_MAP[name]
    except KeyError as exc:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from exc

    module = import_module(mod_name)
    value = getattr(module, attr)
    globals()[name] = value  # cache for subsequent lookups
    return value


def __dir__() -> list[str]:
    """Make IDE autocompletion aware of our lazy attributes."""
    base = set(globals().keys())
    base.update(__all__)
    base.update(_PUBLIC_MAP.keys())
    base.update({"plasmidio", "insilico"})
    return sorted(base)
