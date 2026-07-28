"""DNAcauldron-based assembly report generation & organization."""

from __future__ import annotations

import logging
import os
import re
import zipfile
from collections.abc import Callable, Iterable, Mapping, Sequence
from itertools import product
from pathlib import Path
from typing import Optional, Union

from Bio.SeqRecord import SeqRecord

from .io_ import (
    _coerce_to_repo_keys,
    _default_sequence_loader,
    _is_placeholder,
    load_dna_file,  # noqa: F401 - looked up via globals() below
)

try:
    import dnacauldron as dc
except Exception:  # pragma: no cover
    dc = None

LOGGER = logging.getLogger(__name__)


def generate_assembly_reports(
    *,
    category_dirs: Mapping[str, Union[str, Path]],
    category_order: Sequence[str],
    output_dir: Union[str, Path],
    sequence_loader: Optional[Callable[[Path], SeqRecord]] = None,
    valid_suffixes: Sequence[str] = (".gb", ".gbk", ".genbank", ".dna"),
) -> list[tuple[Path, list[str]]]:
    """Create Golden-Gate assemblies as the Cartesian product across categories.

    This function builds combinatorial assemblies (e.g., **Promoter × RBS × Gene ×
    Terminator × Backbone**). For each **category** (folder) it loads all sequence
    files, normalizes their identifiers, and keys them by the **file stem**. It then
    enumerates all combinations in ``category_order``, creates a per-combination
    repository, simulates the assembly with **DNAcauldron**, and writes a
    ``*_report.zip`` for each result.

    **Hard normalization to avoid ``<unknown id>`` errors**
    After loading each file, the corresponding :class:`Bio.SeqRecord.SeqRecord` is
    normalized so that **``record.id`` and ``record.name`` equal the repository key
    (the file stem)**. Placeholder or missing ``description`` is also set to the stem,
    and minimal GenBank annotations are ensured (``accessions``, ``sequence_version``,
    ``molecule_type``). This guarantees that the DNAcauldron report writer cannot
    look up a placeholder like ``<unknown id>`` in the repository.

    Parameters
    ----------
    category_dirs
        Mapping from category name to folder path. Example:
        ``{"Promoter": ".../Promoter_parts", "RBS": ".../RBS_parts", ...}``.
    category_order
        Order of categories in each assembly (also used for the assembly name),
        e.g., ``["Promoter", "RBS", "Gene", "Terminator", "Backbone"]``.
    output_dir
        Directory where ``*_report.zip`` files will be written. Created if needed.
    sequence_loader
        Optional callable to load a file into a :class:`SeqRecord`. If ``None``,
        a project-level ``load_dna_file(Path)`` is used if available; otherwise a
        fallback loader is used (GenBank + SnapGene if installed).
    valid_suffixes
        File extensions to consider as sequence inputs.

    Returns
    -------
    list of (pathlib.Path, list of str)
        A list of tuples ``(zip_path, part_names_in_order)`` for each assembly.

    Raises
    ------
    FileNotFoundError
        If a provided category folder does not exist or is not a directory.
    RuntimeError
        If a category has no valid sequence files, or if parts cannot be resolved
        to repository keys.
    """
    loader = sequence_loader or _default_sequence_loader
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1) Load all files per category (keyed by file stem)
    catalogs: dict[str, dict[str, SeqRecord]] = {}
    for cat, d in category_dirs.items():
        dpath = Path(d)
        if not dpath.exists() or not dpath.is_dir():
            raise FileNotFoundError(f"Category folder not found: {cat} -> {dpath}")

        files = sorted(p for p in dpath.iterdir() if p.suffix.lower() in valid_suffixes)
        if not files:
            raise RuntimeError(f"No sequence files in category '{cat}' at {dpath}")

        cat_records: dict[str, SeqRecord] = {}
        for fp in files:
            rec = loader(fp)
            key = fp.stem  # repository key = filename stem

            # --- Hard normalization: make the writer <unknown id>-proof ---
            rec.id = key
            rec.name = key
            if _is_placeholder(getattr(rec, "description", None)) or not getattr(
                rec, "description", ""
            ):
                rec.description = key
            accs = rec.annotations.get("accessions")
            if not accs:
                rec.annotations["accessions"] = [key]
            ver = rec.annotations.get("sequence_version")
            if not isinstance(ver, int) or ver <= 0:
                rec.annotations["sequence_version"] = 1
            rec.annotations.setdefault("molecule_type", "DNA")
            # ----------------------------------------------------------------

            cat_records[key] = rec

        catalogs[cat] = cat_records

    # 2) Cartesian product across categories in the specified order
    choices: list[list[str]] = [sorted(catalogs[c].keys()) for c in category_order]
    if any(len(lst) == 0 for lst in choices):
        raise RuntimeError("At least one category is empty after loading.")

    results: list[tuple[Path, list[str]]] = []
    for combo in product(*choices):
        parts_in_order: list[str] = list(combo)

        # 3) Build per-combination repo with just the selected records
        combined: dict[str, SeqRecord] = {
            stem: catalogs[cat][stem]
            for cat, stem in zip(category_order, parts_in_order, strict=False)
        }

        # 4) Optional remap (defensive); should be a no-op with hard normalization
        resolved = _coerce_to_repo_keys(parts_in_order, combined.keys())
        missing = [n for n in resolved if n not in combined]
        if missing:
            sample = ", ".join(list(combined.keys())[:10])
            raise RuntimeError(
                "Parts not found in repository for this combination:\n"
                f"  Missing: {missing}\n"
                f"  Available keys (sample): {sample} ..."
            )

        # 5) Simulate with DNAcauldron
        repo = dc.SequenceRepository(collections={"parts": combined})
        assembly_name = "_".join(resolved)
        assembly = dc.Type2sRestrictionAssembly(parts=resolved, name=assembly_name)
        sim = assembly.simulate(sequence_repository=repo)

        # 6) Write report ZIP
        zip_path = out_dir / f"{assembly_name}_report.zip"
        sim.write_report(str(zip_path))
        LOGGER.info("Wrote report: %s", zip_path.name)

        results.append((zip_path, resolved))

    return results


def generate_recipe_assemblies(
    *,
    category_dirs: Mapping[str, Union[str, Path]],
    designs: Sequence[Mapping[str, str]],
    category_order: Sequence[str],
    output_dir: Union[str, Path],
    sequence_loader: Optional[Callable[[Path], SeqRecord]] = None,
    valid_suffixes: Sequence[str] = (
        ".gb",
        ".gbk",
        ".genbank",
        ".fa",
        ".fasta",
        ".fna",
        ".dna",
    ),
    allow_extra_keys: bool = True,
    category_aliases: Optional[Mapping[str, str]] = None,
) -> list[tuple[Path, list[str]]]:
    """
    Build DNA Cauldron assemblies for the **exact designs provided** (no combinatorics).

    This function loads sequence parts from the folders in ``category_dirs`` using
    a caller-provided loader (``sequence_loader``) or a global ``load_dna_file(Path)``,
    normalizes their identifiers to the **file stem**, resolves each requested part
    name in ``designs`` (exact → case-insensitive → canonical match), and runs one
    DNA Cauldron simulation per design. It writes a single ``*_report.zip`` per design.

    The function is tolerant to *documentation-only* keys in ``designs``:
    when ``allow_extra_keys=True`` (default), unknown keys are ignored. You may also
    supply ``category_aliases`` (e.g., ``{"UNS_Context": "Backbone"}``) to map
    alternate column names to canonical categories.

    Parameters
    ----------
    category_dirs : Mapping[str, str | Path]
        Mapping from **category name** to a directory with sequence files,
        e.g. ``{"Promoter": ".../Promoters", "Backbone": ".../Backbones"}``.
        The directory contents are indexed by **file stem** (filename without extension).
    designs : Sequence[Mapping[str, str]]
        A list of dictionaries, **one per assembly**, that specify which file stem to
        use for each category. Extra keys (IDs, notes, etc.) are ignored when
        ``allow_extra_keys`` is True. Example item::

            {
                "Promoter": "P45_AB_(A)",
                "RBS": "RBS_01_BC_(A)",
                "Gene": "E0030_CD",
                "Terminator": "TrrnB_DE_(A)",
                "Backbone": "UNS1A_UNS3_E"
            }

        You can also use aliases via ``category_aliases``, e.g. the key
        ``"UNS_Context"`` will be treated as ``"Backbone"`` if
        ``{"UNS_Context": "Backbone"}`` is provided.
    category_order : Sequence[str]
        The exact order of categories to use when building each assembly and to form
        the assembly name (``"_".join(resolved_stems)``).
    output_dir : str | Path
        Folder where per-design ZIP reports will be written.
    sequence_loader : Callable[[Path], SeqRecord] | None
        Sequence file loader. If ``None``, a global ``load_dna_file(Path)`` is used
        if present. If neither is available, a ``RuntimeError`` is raised.
    valid_suffixes : Sequence[str]
        File extensions to accept when scanning each category folder.
    allow_extra_keys : bool
        If True (default), keys in each design dict that are **not** in
        ``category_order`` (after aliasing) are ignored. If False, such keys raise.
    category_aliases : Mapping[str, str] | None
        Optional mapping of **incoming design keys** to **canonical category names**,
        e.g. ``{"UNS_Context": "Backbone"}``.

    Returns
    -------
    list[tuple[Path, list[str]]]
        A list with one entry per design: ``(zip_path, resolved_part_names_in_order)``.
        ``zip_path`` is the path to the written report ZIP. The second element is the
        list of **resolved file stems** ordered by ``category_order``.

    Raises
    ------
    FileNotFoundError
        If any category folder is missing.
    RuntimeError
        If no sequence loader is available; if a category folder is empty; if a design
        is missing a required category; if an unknown key is present while
        ``allow_extra_keys`` is False; or if a requested stem cannot be resolved.
    ImportError
        If the **DNA Cauldron** package (``dnacauldron``) is not installed.

    Notes
    -----
    - **Identifier normalization**: for every loaded record, ``id`` and ``name`` are
      set to the file stem, and minimal annotations are ensured to avoid
      ``<unknown id>`` in reports.
    - **Name resolution**: a requested stem is matched in this order:
      *exact match* → *case-insensitive* → *canonical* (lowercased alphanumerics).
    - **DNA Cauldron**: this function uses a simple
      :class:`dnacauldron.Type2sRestrictionAssembly` over the provided parts.
    """
    # ------------------------ helpers ------------------------

    def _is_placeholder(val: Optional[str]) -> bool:
        return val is None or str(val).strip().lower() in {
            "",
            ".",
            "unknown",
            "<unknown id>",
            "<unknown>",
            "nan",
            "null",
            "none",
        }

    def _canonical(text: str) -> str:
        return re.sub(r"[^a-z0-9]", "", text.lower())

    def _resolve_stem(requested: str, keys: Iterable[str]) -> Optional[str]:
        """Resolve ``requested`` in ``keys`` via exact → lower → canonical match."""
        key_list = list(keys)
        if requested in key_list:
            return requested
        lower_map = {k.lower(): k for k in key_list}
        canon_map = {_canonical(k): k for k in key_list}
        req_l = requested.lower()
        if req_l in lower_map:
            return lower_map[req_l]
        req_c = _canonical(requested)
        return canon_map.get(req_c)

    def _load_dnacauldron():
        """Import dnacauldron; raise a clear error if missing."""
        try:
            # pylint: disable=import-error
            import dnacauldron as dc
        except Exception as exc:  # pylint: disable=broad-except
            raise ImportError(
                "dnacauldron is required for this function. "
                "Install with: pip install dnacauldron"
            ) from exc
        return dc

    # ------------------------ setup -------------------------

    # Choose loader: explicit argument → global load_dna_file → error.
    if sequence_loader is not None:
        loader: Callable[[Path], SeqRecord] = sequence_loader
    else:
        global_loader = globals().get("load_dna_file")
        if callable(global_loader):
            loader = global_loader
        else:
            raise RuntimeError(
                "No sequence loader available. Provide `sequence_loader=` or define a "
                "global `load_dna_file(Path) -> SeqRecord`."
            )

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1) Build per-category catalogs (keyed by file stem) with hard normalization.
    catalogs: dict[str, dict[str, SeqRecord]] = {}
    valid_set = {s.lower() for s in valid_suffixes}

    for cat, folder in category_dirs.items():
        dpath = Path(folder)
        if not dpath.is_dir():
            raise FileNotFoundError(f"Category folder not found: {cat} -> {dpath}")

        files = sorted(p for p in dpath.iterdir() if p.suffix.lower() in valid_set)
        if not files:
            raise RuntimeError(f"No sequence files in category '{cat}' at {dpath}")

        cat_records: dict[str, SeqRecord] = {}
        for fp in files:
            rec = loader(fp)
            stem = fp.stem

            # Hard normalization to stable IDs.
            rec.id = stem
            rec.name = stem
            if _is_placeholder(getattr(rec, "description", None)):
                rec.description = stem
            if not rec.annotations.get("accessions"):
                rec.annotations["accessions"] = [stem]
            if not isinstance(rec.annotations.get("sequence_version"), int):
                rec.annotations["sequence_version"] = 1
            rec.annotations.setdefault("molecule_type", "DNA")

            cat_records[stem] = rec

        catalogs[cat] = cat_records

    # 2) Resolve each design and run a DNA Cauldron simulation.
    results: list[tuple[Path, list[str]]] = []
    aliases = dict(category_aliases or {})

    dc = _load_dnacauldron()

    for idx, raw_design in enumerate(designs, start=1):
        # Apply aliases and (optionally) ignore documentation keys.
        normalized: dict[str, str] = {}
        unknown_keys: list[str] = []

        for key, value in raw_design.items():
            target = aliases.get(key, key)
            if target in category_order:
                normalized[target] = value
            else:
                unknown_keys.append(key)

        if unknown_keys and not allow_extra_keys:
            raise RuntimeError(
                f"Design #{idx} contains unknown keys not in category_order: {unknown_keys}"
            )
        if unknown_keys and allow_extra_keys:
            LOGGER.info(
                "Design #%d: ignoring extra keys: %s",
                idx,
                ", ".join(sorted(unknown_keys)),
            )

        # Ensure all required categories are present.
        missing = [c for c in category_order if c not in normalized]
        if missing:
            raise RuntimeError(f"Design #{idx} is missing categories: {missing}")

        # Resolve requested stems to actual catalog keys.
        resolved: list[str] = []
        for cat in category_order:
            requested_stem = normalized[cat]
            available = catalogs[cat].keys()
            match = _resolve_stem(requested_stem, available)
            if match is None:
                examples = ", ".join(list(available)[:10])
                raise RuntimeError(
                    f"Design #{idx}: requested stem '{requested_stem}' not found in "
                    f"category '{cat}'. Examples: {examples} ..."
                )
            resolved.append(match)

        # Build a small repository of the selected parts.
        repo_parts: dict[str, SeqRecord] = {
            stem: catalogs[cat][stem]
            for cat, stem in zip(category_order, resolved, strict=False)
        }
        repository = dc.SequenceRepository(collections={"parts": repo_parts})

        assembly_name = "_".join(resolved)
        assembly = dc.Type2sRestrictionAssembly(parts=resolved, name=assembly_name)
        simulation = assembly.simulate(sequence_repository=repository)

        zip_path = out_dir / f"{assembly_name}_report.zip"
        simulation.write_report(str(zip_path))
        LOGGER.info("Wrote report: %s", zip_path)

        results.append((zip_path, resolved))

    return results


def organize_assembly_reports(
    report_dir: str | Path,
    reports: Optional[Sequence[tuple[Path, list[str]]]] = None,
    *,
    extract_subdir: str = "Assembly",
    delete_zip: bool = False,
    final_only: bool = True,
) -> list[Path]:
    """Extract GenBank files from DNAcauldron ``*_report.zip`` archives.

    By default, the function extracts exactly **one final construct per ZIP**
    (``final_only=True``). The final construct is identified by matching the
    ZIP's assembly name (ZIP stem without the ``_report`` suffix) to a file
    named ``<assembly_name>.gb`` or ``<assembly_name>.gbk`` inside the archive.
    If no exact match is present, it falls back to any GenBank file located in a
    path containing the substring ``"construct"``. As a last resort (to avoid
    over-extraction), it picks the first ``.gb``/``.gbk`` it finds.

    Parameters
    ----------
    report_dir : str | pathlib.Path
        Directory containing the report ZIPs (e.g., ``"reports"``).
    reports : Sequence[tuple[pathlib.Path, list[str]]] | None, optional
        Optional result of an assembly-generation function; if provided, its ZIP
        paths are used instead of scanning ``report_dir``.
    extract_subdir : str, optional
        Subdirectory (under ``report_dir``) where extracted files will be saved.
        Default is ``"Assembly"``.
    delete_zip : bool, optional
        If ``True``, delete each ZIP after successful extraction. Default ``False``.
    final_only : bool, optional
        If ``True`` (default), extract only the final construct per ZIP using the
        identification strategy described above. If ``False``, extract **all**
        ``.gb``/``.gbk`` files found in each ZIP.

    Returns
    -------
    list[pathlib.Path]
        Paths to the extracted GenBank files.

    Raises
    ------
    FileNotFoundError
        If ``report_dir`` does not exist.

    Notes
    -----
    - When both ``.gb`` and ``.gbk`` candidates exist for the chosen file, the
      ``.gb`` variant is preferred.
    - Filenames on disk are normalized to ``<assembly_name>.<ext>`` when
      ``final_only=True`` to avoid collisions and keep outputs tidy.
    """
    base = Path(report_dir)
    if not base.exists():
        raise FileNotFoundError(f"Report folder not found: {base}")

    # Determine ZIPs to process: passed-in reports or scan the directory.
    zip_paths: list[Path]
    if reports is not None:
        zip_paths = [p for (p, _parts) in reports]
    else:
        zip_paths = sorted(base.glob("*_report.zip"))

    out_dir = base / extract_subdir
    out_dir.mkdir(parents=True, exist_ok=True)

    extracted: list[Path] = []

    for zpath in zip_paths:
        if not zpath.exists():
            LOGGER.warning("Skipping missing ZIP: %s", zpath)
            continue

        with zipfile.ZipFile(zpath) as zf:
            members = [m for m in zf.infolist() if not m.is_dir()]

            candidates: list[zipfile.ZipInfo] = []
            if final_only:
                # Assembly name = ZIP stem without trailing "_report"
                asm_stem = (
                    zpath.stem[:-7] if zpath.stem.endswith("_report") else zpath.stem
                )
                wanted = {f"{asm_stem}.gb".lower(), f"{asm_stem}.gbk".lower()}

                # 1) Exact filename match
                candidates = [
                    m for m in members if Path(m.filename).name.lower() in wanted
                ]

                # 2) Fallback: any GenBank file under a "construct" path segment
                if not candidates:
                    candidates = [
                        m
                        for m in members
                        if "construct" in m.filename.lower()
                        and m.filename.lower().endswith((".gb", ".gbk"))
                    ]

                # 3) Last resort: first GenBank file (keeps output to one per ZIP)
                if not candidates:
                    for m in members:
                        if m.filename.lower().endswith((".gb", ".gbk")):
                            candidates = [m]
                            break
            else:
                candidates = [
                    m for m in members if m.filename.lower().endswith((".gb", ".gbk"))
                ]

            if not candidates:
                LOGGER.warning("No GenBank candidates in %s", zpath.name)
                if delete_zip:
                    zpath.unlink(missing_ok=True)
                continue

            # Prefer .gb over .gbk (and shorter paths if otherwise equal)
            candidates.sort(
                key=lambda m: (
                    Path(m.filename).suffix.lower() != ".gb",
                    len(m.filename),
                )
            )
            chosen = candidates[0]

            # Destination filename
            asm_stem = zpath.stem[:-7] if zpath.stem.endswith("_report") else zpath.stem
            ext = Path(chosen.filename).suffix.lower()
            out_name = (
                f"{asm_stem}{ext}"
                if final_only
                else f"{asm_stem}__{Path(chosen.filename).name}"
            )
            dest = out_dir / out_name

            with zf.open(chosen) as src, dest.open("wb") as dst:
                dst.write(src.read())

            extracted.append(dest)
            LOGGER.info("Extracted %s -> %s", chosen.filename, dest)

        if delete_zip:
            zpath.unlink(missing_ok=True)

    return extracted


def delete_all_zip_files(folder_path: str | os.PathLike[str]) -> None:
    """Deleting all `.zip` files in a folder.

    Parameters
    ----------
    folder_path : str | os.PathLike
        Path to the folder where ZIP files are to be deleted.

    Returns
    -------
    None
    """
    folder = Path(folder_path)
    for p in folder.glob("*.zip"):
        try:
            p.unlink()
            print(f"Deleted: {p.name}")
        except OSError as exc:
            print(f"Failed to delete '{p.name}': {exc}")


def run_all_functions(
    folders: Mapping[str, Union[str, Path]] | Iterable[Union[str, Path]],
    *,
    report_dir: Union[str, Path] = ".",
    delete_zip: bool = False,
    category_order: Sequence[str] | None = None,
):
    """
    Build DNA Cauldron reports and organize them into `report_dir`.

    Parameters
    ----------
    folders
        Either a mapping of {category_name -> folder} or a plain iterable of
        folders. When an iterable is given, category names are inferred from the
        folder basenames and `category_order` defaults to that order.
    report_dir
        Output directory for the reports.
    delete_zip
        If True, remove the original *_report.zip files after organizing.
    category_order
        Order of categories when generating the combinatorial assemblies.
    """
    # Normalize to a mapping
    if isinstance(folders, Mapping):
        category_dirs = {str(k): Path(v) for k, v in folders.items()}
        if category_order is None:
            category_order = list(category_dirs.keys())
    else:
        seq = [Path(p) for p in folders]
        category_dirs = {p.name: p for p in seq}
        if category_order is None:
            category_order = [p.name for p in seq]

    report_dir = Path(report_dir)
    report_dir.mkdir(parents=True, exist_ok=True)

    reports = generate_assembly_reports(
        category_dirs=category_dirs,
        category_order=category_order,
        output_dir=report_dir,
    )
    return organize_assembly_reports(report_dir, reports, delete_zip=delete_zip)
