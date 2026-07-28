# assembly_designer/plasmidio/__init__.py
"""
Plasmidio — in-silico plasmid design & analysis toolkit.

This package-level ``__init__`` re-exports the public API from the
per-concern implementation modules (``io_.py``, ``assembly.py``, ``pcr.py``,
``batch.py``, ``plotting.py``, ``history.py``, ``s1_doc.py``, ...) so you can
write::

    from assembly_designer.plasmidio import <function/class>

without depending on internal layout. The goal is to keep the import surface
stable and explicit.

Highlights
----------
- **I/O & feature hygiene**
  - Read GenBank / SnapGene (``load_dna_file``, ``read_plasmid_file``,
    ``snapgene_to_seqrecord``)
  - Feature filtering / de-duplication (``filter_features``,
    ``remove_near_duplicate_features``)
  - Part loaders by folder (``load_parts_from_folders``)

- **Assembly reports (DnaCauldron)**
  - Combinatorial builder: ``generate_assembly_reports`` (cartesian product)
  - Recipe builder: ``generate_recipe_assemblies`` (explicit stems only)
  - Report organization / extraction: ``organize_assembly_reports``,
    ``run_all_functions``, ``delete_all_zip_files``

- **Restriction digest & gel simulation**
  - ``perform_restriction_digest``, ``get_ladder``, ``plot_gel``,
    ``digest_and_plot_plasmids``

- **Batch analysis & exports**
  - Prepare and align reads: ``prepare_plasmids``, ``load_reads``, ``align_batch``
  - **Parallel alignment** (threads/processes, optional): ``align_batch_parallel``
  - Tabular outputs: ``results_to_dataframe``, ``export_results``,
    ``export_feature_inventory``
  - Convenience pickers: ``pick_rows_by_ref``, ``pick_row_by_ref``
  - **Easy wrapper (optional):** ``dealign_reads_easy``  # ← ADDED (optional)

- **Visualization**
  - Read-anchored alignment plots: ``plot_alignment``, ``plot_alignment_by_ref``
  - HTML/inline helpers: ``view_alignment*``
  - **Note:** ``plot_alignment`` supports
    ``feature_mode={"both","full","window","none"}`` to avoid stacked or
    overlapping feature tracks; use ``feature_mode="window"`` for a clean zoom.

- **PCR / Gibson helpers**
  - ``simulate_pcr``, ``simulate_pcr_overhangs``, ``merge_with_gibson_features``,
    ``longest_overlap``
  - ``PCRResult`` dataclass for rich PCR metadata

- **3G batch pipeline**
  - ``run_3g_batch_safe`` and enum ``BuildStatus`` (Golden Gate → PCR → Gibson)

- **History / provenance**
  - ``AssemblyHistory``, ``HistoryNode``, ``HistoryEdge``,
    ``build_histories_for_all_constructs`` (build & plot assembly DAGs)

- **Template resolver**
  - ``make_template_lookup`` to resolve TU templates from UNS context

Backward compatibility
----------------------
- ``generate_combinatorial_assemblies`` is an alias of
  ``generate_assembly_reports``.

This docstring documents what is re-exported via ``__all__``; for details,
see the underlying implementations in the per-concern modules listed above.
"""


from __future__ import annotations

# --------------------------------------------------------------------------- #
# Public API, re-exported from the per-concern implementation modules         #
# --------------------------------------------------------------------------- #
from .alignment import read_anchored_align as read_anchored_align  # noqa: F401
from .alignment import read_anchored_align_edlib  # low-level aligner
from .assembly import (
    delete_all_zip_files,
    generate_assembly_reports,
    generate_recipe_assemblies,  # explicit recipe builder
    organize_assembly_reports,
    run_all_functions,
)
from .batch import (
    PlasmidRef,
    ReadItem,
    ResultRow,
    align_batch,
    export_feature_inventory,
    export_results,
    find_unbuilt_plasmids,
    load_reads,
    prepare_plasmids,
    results_to_dataframe,
    summarize_matches,
)
from .batch import align_batch_parallel as align_batch_parallel  # noqa: F401
from .batch import dealign_reads_easy as dealign_reads_easy  # noqa: F401
from .features import filter_features
from .gel import (
    digest_and_plot_plasmids,
    get_ladder,
    perform_restriction_digest,
    plot_gel,
    read_plasmid_file,
)
from .history import (
    AssemblyHistory,
    HistoryEdge,
    HistoryNode,
    PCRResultLike,
    build_histories_for_all_constructs,
)
from .io_ import (
    _safe_filename,  # Intentional private export (used in tests/tools)
    load_dna_file,
    load_parts_from_folders,
    remove_near_duplicate_features,
    snapgene_to_seqrecord,
)
from .lexicon import (
    Lexicon,
    compile_lexicon,
    default_lexicon,
    match_token,
    normalize_ref_id,
    revcomp,
)
from .mapping_summary import plot_alignment_results_summary
from .naming import build_construct_label
from .pcr import (
    BuildStatus,
    PCRResult,
    longest_overlap,
    merge_with_gibson_features,
    run_3g_batch_safe,  # 3G batch runner (Golden Gate → PCR → Gibson)
    simulate_pcr,
    simulate_pcr_overhangs,
)
from .plotting import (
    pick_row_by_ref,
    pick_rows_by_ref,
    plot_alignment,
    plot_alignment_by_ref,
    view_alignment,
    view_alignment_by_ref,
    view_alignment_html,
    view_alignment_html_by_ref,
)
from .template import make_template_lookup

# --- Optional: S1 document helpers (re-export only if their extra deps are installed) ---
try:
    from .s1_doc import build_plasmid_catalog_docx as build_plasmid_catalog_docx
    from .s1_doc import build_s1_documentation as build_s1_documentation

    _HAS_S1 = True
except Exception:  # pragma: no cover
    _HAS_S1 = False

_HAS_EASY = True
_HAS_PARALLEL = True


# --------------------------------------------------------------------------- #
# Backward-compat alias (old name → combinatorial builder)                    #
# --------------------------------------------------------------------------- #
generate_combinatorial_assemblies = generate_assembly_reports

# --------------------------------------------------------------------------- #
# Export surface
# --------------------------------------------------------------------------- #
__all__ = [
    # I/O & features
    "load_dna_file",
    "snapgene_to_seqrecord",
    "read_plasmid_file",
    "filter_features",
    "remove_near_duplicate_features",
    "load_parts_from_folders",
    # Reports / DnaCauldron pipeline
    "generate_assembly_reports",
    "generate_recipe_assemblies",
    "generate_combinatorial_assemblies",
    "organize_assembly_reports",
    "run_all_functions",
    "delete_all_zip_files",
    # Restriction digest & gel simulation
    "get_ladder",
    "perform_restriction_digest",
    "plot_gel",
    "digest_and_plot_plasmids",
    # Batch orchestration
    "PlasmidRef",
    "ReadItem",
    "ResultRow",
    "prepare_plasmids",
    "load_reads",
    "align_batch",
    "results_to_dataframe",
    "export_results",
    "export_feature_inventory",
    # Pickers
    "pick_rows_by_ref",
    "pick_row_by_ref",
    # Viewers/Plotters
    "plot_alignment",
    "plot_alignment_by_ref",
    "view_alignment",
    "view_alignment_by_ref",
    "view_alignment_html",
    "view_alignment_html_by_ref",
    "read_anchored_align_edlib",  # low-level aligner (optional)
    "summarize_matches",
    "find_unbuilt_plasmids",
    "plot_alignment_results_summary",
    # --- PCR & Gibson assembly API ---
    "PCRResult",
    "simulate_pcr",
    "simulate_pcr_overhangs",
    "merge_with_gibson_features",
    "longest_overlap",
    # --- Naming helper ---
    "build_construct_label",
    # --- NEW: 3G batch runner ---
    "run_3g_batch_safe",
    "BuildStatus",
    # --- NEW: history / provenance ---
    "PCRResultLike",
    "HistoryNode",
    "HistoryEdge",
    "AssemblyHistory",
    "build_histories_for_all_constructs",
    # --- NEW: template resolver factory ---
    "make_template_lookup",
    # Private (intentional) exports
    "_safe_filename",
    "normalize_ref_id",
    "revcomp",
    "Lexicon",
    "compile_lexicon",
    "default_lexicon",
    "match_token",
]

# Append S1-exports only if available
if _HAS_S1:
    __all__.extend(
        [
            "build_s1_documentation",
            "build_plasmid_catalog_docx",
        ]
    )

# Only append parallel-export if available
if _HAS_PARALLEL:
    __all__.append("align_batch_parallel")

# Only attach the easy wrapper and low-level aligner if available
if _HAS_EASY:
    __all__.extend(["dealign_reads_easy", "read_anchored_align"])

__version__ = "0.9.1"
