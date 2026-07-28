"""I/O & feature hygiene: GenBank/SnapGene readers, header normalization, part loaders."""

from __future__ import annotations

import hashlib
import re
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from .features import dedupe_record_features

try:
    from snapgene_reader import snapgene_file_to_dict
except Exception:  # pragma: no cover
    snapgene_file_to_dict = None


def snapgene_to_seqrecord(path: str | Path) -> SeqRecord:
    """
    Convert a SnapGene ``.dna`` file to a :class:`Bio.SeqRecord.SeqRecord`.

    The function normalizes ``id``/``name`` (GenBank LOCUS-safe),
    preserves topology (circular/linear) when available, copies qualifiers,
    and removes near-duplicate features via :func:`dedupe_record_features`.

    Parameters
    ----------
    path
        File path to a SnapGene ``.dna`` file.

    Returns
    -------
    Bio.SeqRecord.SeqRecord
        Parsed, normalized record with deduplicated features.

    Raises
    ------
    ImportError
        If ``snapgene-reader`` is not installed.
    ValueError
        If the SnapGene file does not contain a sequence string.

    Notes
    -----
    - LOCUS/``id``/``name`` are sanitized to avoid spaces and overlong values.
      This improves downstream GenBank writing and repo handling.
    - Unknown/empty names fall back to the filename stem.
    """
    if snapgene_file_to_dict is None:
        raise ImportError(
            "snapgene-reader is required for .dna files. "
            "Install via `pip install snapgene-reader`."
        )

    p = Path(path)
    data = snapgene_file_to_dict(str(p))

    seq_str = data.get("seq") or ""
    if not isinstance(seq_str, str) or not seq_str:
        raise ValueError(f"SnapGene file has no 'seq' content: {p}")

    seq = Seq(seq_str)

    # Build features
    strand_map: dict[str | None, int] = {"+": 1, "-": -1, None: 0}
    feats: list[SeqFeature] = []
    for f in data.get("features", []) or []:
        try:
            start_i = int(f["start"])
            end_i = int(f["end"])
        except (KeyError, TypeError, ValueError) as exc:
            # Skip malformed feature, keep robust.
            _ = exc
            continue
        feats.append(
            SeqFeature(
                location=FeatureLocation(
                    start=start_i,
                    end=end_i,
                    strand=strand_map.get(f.get("strand"), 0),
                ),
                type=str(f.get("type", "feature")),
                qualifiers=f.get("qualifiers", {}) or {},
            )
        )

    # Safe, LOCUS-friendly identifiers
    stem = p.stem
    raw_id = str(data.get("id") or data.get("name") or stem) or stem
    raw_name = str(data.get("name") or data.get("id") or stem) or stem

    def _safe_locus(text: str, maxlen: int = 20) -> str:
        txt = re.sub(r"\s+", "_", text.strip())
        return txt[:maxlen] or "record"

    rec_id = _safe_locus(raw_id)
    rec_name = _safe_locus(raw_name)

    # Topology and annotations
    dna_meta = data.get("dna", {}) or {}
    topo = dna_meta.get("topology")
    if not topo:
        topo = "circular" if bool(data.get("is_circular")) else "linear"

    annotations: dict[str, object] = {
        "molecule_type": "DNA",
        "topology": topo,
    }
    notes = data.get("notes") or {}
    if isinstance(notes, dict):
        annotations.update(notes)

    rec = SeqRecord(
        seq=seq,
        id=rec_id,
        name=rec_name,
        features=feats,
        annotations=annotations,
    )

    return dedupe_record_features(rec)


_PLACEHOLDERS = {"", ".", "unknown", "<unknown id>", "<unknown>", "none", "null", "nan"}


def _is_placeholder(name: Optional[str]) -> bool:
    """Return True if `name` is None/empty or a known placeholder."""
    return name is None or name.strip().lower() in _PLACEHOLDERS


def _canonical(s: str) -> str:
    """Canonicalize a part name for fuzzy lookup (lowercased, non-alnum removed)."""
    return re.sub(r"[^a-z0-9]", "", s.lower())


def _coerce_to_repo_keys(parts: Sequence[str], repo_keys: Iterable[str]) -> list[str]:
    """Map requested part names to repository keys using several fallbacks.

    Strategy:
    1) Exact match
    2) Case-insensitive match
    3) Canonical match (lowercase + strip non-alphanumeric)

    Parameters
    ----------
    parts
        Requested part names (may include placeholders like ``"<unknown id>"``).
    repo_keys
        Keys available in the SequenceRepository.

    Returns
    -------
    list of str
        Resolved repository keys in the same order as ``parts``. Unresolved
        names are returned unchanged so a subsequent missing-check can fail
        with a precise error.
    """
    keys = list(repo_keys)
    exact = {k: k for k in keys}
    lower = {k.lower(): k for k in keys}
    canon = {_canonical(k): k for k in keys}

    resolved: list[str] = []
    for p in parts:
        if p in exact:
            resolved.append(p)
            continue
        pl = p.lower()
        if pl in lower:
            resolved.append(lower[pl])
            continue
        pc = _canonical(p)
        if pc in canon:
            resolved.append(canon[pc])
            continue
        resolved.append(p)  # unchanged; will be flagged later
    return resolved


def _needs_fix(value: Optional[str]) -> bool:
    """Return True if a string is empty/None or equals 'Unknown' or '.'."""
    return value is None or value.strip() in {"", "Unknown", "."}


def normalize_genbank_header_from_filename(
    record: SeqRecord, *, stem: str
) -> SeqRecord:
    """Normalize critical GenBank header fields using the file stem.

    Fills these when missing/placeholder:
    - LOCUS name  -> record.id (and record.name)
    - DEFINITION  -> record.description
    - ACCESSION   -> record.annotations['accessions']  (list)
    - VERSION     -> record.annotations['sequence_version'] (int; set to 1)
    - KEYWORDS    -> record.annotations['keywords'] (list)
    Also ensures molecule_type='DNA' (helps some toolchains).
    """
    # LOCUS (name)
    if _needs_fix(getattr(record, "id", None)):
        record.id = stem
    if _needs_fix(getattr(record, "name", None)):
        record.name = stem

    # DEFINITION
    if _needs_fix(getattr(record, "description", None)):
        record.description = stem

    # ACCESSION (list)
    accs = record.annotations.get("accessions")
    if not accs or all(_needs_fix(a) for a in accs):
        record.annotations["accessions"] = [stem]

    # VERSION (integer in Biopython; use 1 as a sane default)
    ver = record.annotations.get("sequence_version")
    if ver in (None, 0) or (isinstance(ver, str) and _needs_fix(ver)):
        record.annotations["sequence_version"] = 1

    # KEYWORDS (list)
    kws = record.annotations.get("keywords")
    if not kws or (isinstance(kws, list) and all(_needs_fix(k) for k in kws)):
        record.annotations["keywords"] = [stem]

    # Molecule type
    if "molecule_type" not in record.annotations:
        record.annotations["molecule_type"] = "DNA"

    return record


def load_dna_file(path: Path) -> SeqRecord:
    """Load a DNA file (GenBank or SnapGene) and return a normalized SeqRecord.

    Supported formats:
    - GenBank: ``.gb``, ``.gbk``, ``.genbank``
    - SnapGene: ``.dna`` (requires ``snapgene_reader``)

    Normalization (prevents ``<unknown id>``/``.`` issues):
    - For **GenBank**, the following placeholders (empty/``"."``/``"Unknown"``)
      are replaced with the **file stem**:
        * LOCUS name → ``record.id`` and ``record.name``
        * DEFINITION → ``record.description``
        * ACCESSION  → ``record.annotations["accessions"]`` (list)
        * VERSION    → ``record.annotations["sequence_version"]`` (set to ``1``)
        * KEYWORDS   → ``record.annotations["keywords"]`` (list)
    - For **SnapGene**, placeholder values in ``id``/``name``/``description`` are
      also replaced with the **file stem**.
    - Ensures ``record.annotations["molecule_type"] == "DNA"`` when missing.

    Parameters
    ----------
    path : pathlib.Path
        Path to the input file.

    Returns
    -------
    Bio.SeqRecord.SeqRecord
        Parsed and normalized sequence record.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    ValueError
        If the file extension is unsupported.
    RuntimeError
        If reading a ``.dna`` file but ``snapgene_reader`` is unavailable.
    """
    if not path.exists():
        raise FileNotFoundError(f"Sequence file not found: {path}")

    suffix = path.suffix.lower()
    stem = path.stem

    # ------------------------- GenBank -------------------------
    if suffix in {".gb", ".gbk", ".genbank"}:
        record: SeqRecord = SeqIO.read(str(path), "genbank")

        # LOCUS → id/name
        if _needs_fix(getattr(record, "id", None)):
            record.id = stem
        if _needs_fix(getattr(record, "name", None)):
            record.name = stem

        # DEFINITION → description
        if _needs_fix(getattr(record, "description", None)):
            record.description = stem

        # ACCESSION → annotations['accessions'] (list)
        accs = record.annotations.get("accessions")
        if not accs or all(_needs_fix(a) for a in accs):
            record.annotations["accessions"] = [stem]

        # VERSION → annotations['sequence_version'] (int)
        ver = record.annotations.get("sequence_version")
        if not isinstance(ver, int) or ver <= 0:
            record.annotations["sequence_version"] = 1

        # KEYWORDS → annotations['keywords'] (list)
        kws = record.annotations.get("keywords")
        if not kws or (isinstance(kws, list) and all(_needs_fix(k) for k in kws)):
            record.annotations["keywords"] = [stem]

        # Molecule type
        if "molecule_type" not in record.annotations:
            record.annotations["molecule_type"] = "DNA"

        return record

    # ------------------------- SnapGene ------------------------
    if suffix == ".dna":
        try:
            try:
                # Preferred: directly returns a SeqRecord
                from snapgene_reader import (
                    snapgene_file_to_seqrecord,
                )

                record = snapgene_file_to_seqrecord(str(path))
            except Exception:
                # Fallback: dict API → minimal SeqRecord
                from snapgene_reader import (
                    snapgene_file_to_dict,
                )

                data = snapgene_file_to_dict(str(path))
                seq_str = str(data.get("seq", ""))
                record = SeqRecord(Seq(seq_str), id=stem, name=stem, description=stem)

            # Normalize placeholders for id/name/description
            if _needs_fix(getattr(record, "id", None)):
                record.id = stem
            if _needs_fix(getattr(record, "name", None)):
                record.name = stem
            if _needs_fix(getattr(record, "description", None)):
                record.description = stem
            if "molecule_type" not in record.annotations:
                record.annotations["molecule_type"] = "DNA"
            return record
        except Exception as exc:  # pragma: no cover
            raise RuntimeError(
                "Reading SnapGene '.dna' requires 'snapgene_reader'. "
                "Install with: pip install snapgene-reader"
            ) from exc

    # ----------------------- Unsupported -----------------------
    raise ValueError(
        f"Unsupported file extension '{suffix}' for {path.name}. "
        "Supported: .gb, .gbk, .genbank, .dna"
    )


def load_parts_from_folders(
    folders: Iterable[str | Path],
) -> dict[str, dict[str, SeqRecord]]:
    """Load all DNA parts from one or more folders into a nested dictionary.

    Supported file types per folder:
    - SnapGene: ``.dna``
    - GenBank: ``.gb``, ``.gbk``, ``.genbank`` (features are de-duplicated)
    - FASTA: ``.fa``, ``.fasta``

    The return structure is:
    ``{ <folder_name>: { <part_name>: SeqRecord, ... }, ... }``,
    where ``part_name`` is the file stem (filename without extension).

    Parameters
    ----------
    folder_paths
        Iterable of directory paths to scan (non-recursive).

    Returns
    -------
    dict[str, dict[str, Bio.SeqRecord.SeqRecord]]
        A mapping from folder *basename* to a mapping of part names to
        :class:`~Bio.SeqRecord.SeqRecord` objects.

    Raises
    ------
    FileNotFoundError
        If any of the given folders does not exist.
    NotADirectoryError
        If any given path exists but is not a directory.

    Notes
    -----
    - Unknown file extensions are ignored.
    - If two files in the same folder share the same stem (e.g., ``x.gb`` and
      ``x.fasta``), the **last** one seen will overwrite the earlier entry.
    - Parsing errors for individual files are skipped; you can add logging
      around :func:`load_dna_file` if you want visibility.

    Examples
    --------
    >>> parts = load_parts_from_folders(["vectors", "promoters"])
    >>> "vectors" in parts and isinstance(parts["vectors"], dict)
    True
    >>> any(hasattr(rec, "seq") for rec in parts["vectors"].values())
    True
    """
    result: dict[str, dict[str, SeqRecord]] = {}
    for folder in folders:
        fpath = Path(folder)
        if not fpath.is_dir():
            raise FileNotFoundError(f"Not a directory: {fpath}")
        container: dict[str, SeqRecord] = {}
        for file in sorted(fpath.iterdir()):
            if not file.is_file():
                continue
            if file.suffix.lower() not in {
                ".dna",
                ".gb",
                ".gbk",
                ".genbank",
                ".fa",
                ".fasta",
            }:
                continue
            rec = load_dna_file(file)
            container[file.stem] = rec
        result[fpath.name] = container
    return result


def _safe_filename(stem: str, maxlen: int = 120) -> str:
    """Return a cross-platform file-system safe filename stem.

    Sanitizes characters that are invalid on Windows (``<>:\"/\\|?*`` and
    whitespace), trims leading/trailing dots/underscores, and shortens overly
    long names by appending an 8-character SHA1 suffix to preserve uniqueness.

    Parameters
    ----------
    stem : str
        The unsanitized file name (without extension).
    maxlen : int, optional
        Maximum length for the sanitized stem (default: 120).

    Returns
    -------
    str
        A safe, normalized file name stem suitable for all platforms.
    """
    safe = re.sub(r'[<>:"/\\|?*\s]+', "_", stem).strip("._")
    if len(safe) > maxlen:
        digest = hashlib.sha1(stem.encode("utf-8")).hexdigest()[:8]
        safe = f"{safe[: maxlen - 9]}-{digest}"
    return safe


def _default_sequence_loader(path: Path) -> SeqRecord:
    """Load a sequence file into a normalized :class:`Bio.SeqRecord.SeqRecord`.

    Prefers a project-level ``load_dna_file(Path)`` if present in the global
    namespace. Otherwise falls back to:
    - GenBank: ``.gb``, ``.gbk``, ``.genbank`` (via Biopython)
    - SnapGene: ``.dna`` (requires ``snapgene_reader``; dict fallback builds a minimal SeqRecord)

    Normalization (prevents ``<unknown id>``/``.`` issues):
    - Placeholder values in ``id``/``name``/``description`` (empty/``"."``/``"Unknown"``/``"<unknown id>"``)
      are replaced with the **file stem**.
    - Ensures ``record.annotations["molecule_type"] == "DNA"``.
    """
    # Prefer your robust project loader if available
    loader = globals().get("load_dna_file")
    if callable(loader):
        return loader(path)

    suffix = path.suffix.lower()
    stem = path.stem

    if suffix in {".gb", ".gbk", ".genbank"}:
        rec: SeqRecord = SeqIO.read(str(path), "genbank")
        if _is_placeholder(getattr(rec, "id", None)):
            rec.id = stem
        if _is_placeholder(getattr(rec, "name", None)):
            rec.name = stem
        if _is_placeholder(getattr(rec, "description", None)):
            rec.description = stem
        rec.annotations.setdefault("molecule_type", "DNA")
        return rec

    if suffix == ".dna":
        try:
            try:
                from snapgene_reader import (
                    snapgene_file_to_seqrecord,
                )

                rec = snapgene_file_to_seqrecord(str(path))
            except Exception:
                from snapgene_reader import (
                    snapgene_file_to_dict,
                )

                data = snapgene_file_to_dict(str(path))
                seq_str = str(data.get("seq", ""))
                rec = SeqRecord(Seq(seq_str), id=stem, name=stem, description=stem)
            if _is_placeholder(getattr(rec, "id", None)):
                rec.id = stem
            if _is_placeholder(getattr(rec, "name", None)):
                rec.name = stem
            if _is_placeholder(getattr(rec, "description", None)):
                rec.description = stem
            rec.annotations.setdefault("molecule_type", "DNA")
            return rec
        except Exception as exc:
            raise RuntimeError(
                "Reading SnapGene '.dna' requires 'snapgene_reader'. "
                "Install with: pip install snapgene-reader"
            ) from exc

    raise RuntimeError(
        f"No loader for '{path.suffix}'. Provide 'load_dna_file(Path)' to support this format."
    )


def remove_near_duplicate_features(file_path: str | Path, tolerance: int = 3) -> Path:
    """
    Remove near-duplicate features from a GenBank file based on (start, end)
    coordinates within a positional tolerance.

    Two features are considered duplicates if both their start and end positions
    fall within ``±tolerance`` bases of an already-seen feature. The cleaned
    record is written **to the same path** (in-place overwrite).

    Parameters
    ----------
    file_path : str or pathlib.Path
        Path to the GenBank file (``.gb``/``.gbk``/``.genbank``).
    tolerance : int, optional
        Allowed deviation (in bases) for start and end positions, by default 3.

    Returns
    -------
    pathlib.Path
        Path of the cleaned (overwritten) GenBank file.

    Raises
    ------
    FileNotFoundError
        If ``file_path`` does not exist.
    ValueError
        If ``tolerance`` is negative.
    OSError
        If writing the output file fails.

    Notes
    -----
    - Only the (start, end) coordinates are compared; strand, type, and
      qualifiers are **not** considered.
    - Complex locations (e.g., ``join(...)``) are treated via their extremal
      coordinates through ``int()`` conversion.
    """
    path_obj = Path(file_path)
    if not path_obj.exists():
        raise FileNotFoundError(f"File not found: {path_obj}")
    if tolerance < 0:
        raise ValueError("tolerance must be >= 0")

    record: SeqRecord = SeqIO.read(str(path_obj), "genbank")

    unique_positions: list[tuple[int, int]] = []
    unique_features: list[SeqFeature] = []

    for feature in record.features:
        feature_start = int(feature.location.start)
        feature_end = int(feature.location.end)

        is_duplicate = any(
            abs(feature_start - start) <= tolerance
            and abs(feature_end - end) <= tolerance
            for start, end in unique_positions
        )
        if not is_duplicate:
            unique_positions.append((feature_start, feature_end))
            unique_features.append(feature)

    record.features = unique_features

    with path_obj.open("w") as output_handle:
        SeqIO.write(record, output_handle, "genbank")

    print(f"✅ Cleaned file saved: {path_obj}")
    return path_obj
