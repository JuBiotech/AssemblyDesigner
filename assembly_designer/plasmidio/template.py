"""Template resolver factory for PCR/3G templates."""

from __future__ import annotations

import re
from collections.abc import Callable
from functools import cache
from pathlib import Path
from typing import Any, Optional

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

_UNS_NUM = re.compile(r"uns\s*([0-9]+)", flags=re.I)


def _as_bool(value: Any, default: bool = True) -> bool:
    """Convert a fuzzy value to ``bool``.

    Parameters
    ----------
    value
        Arbitrary input (bool/int/str/None).
    default
        Fallback when conversion is ambiguous.

    Returns
    -------
    bool
        Parsed truth value.

    Examples
    --------
    >>> _as_bool("yes")
    True
    >>> _as_bool("0")
    False
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return default
    s = str(value).strip().lower()
    if s in {"1", "true", "yes", "y"}:
        return True
    if s in {"0", "false", "no", "n"}:
        return False
    return default


def _get_int(value: Any, default: int) -> int:
    """Convert a value to ``int`` with a fallback.

    Parameters
    ----------
    value
        Arbitrary input.
    default
        Fallback if conversion fails.

    Returns
    -------
    int
        Parsed integer or ``default``.
    """
    try:
        return int(value)
    except Exception:  # pylint: disable=broad-except
        return default


def make_template_lookup(
    *,
    reports_dir: Path,
    tus_df: Optional[pd.DataFrame],
    vector_rec: SeqRecord,
) -> Callable[[str, str], SeqRecord]:
    """
    Creates a resolver function to retrieve template `SeqRecord` objects for a given construct ID and role.
    Parameters
    ----------
    reports_dir : Path
        Path to the directory containing assembly reports. The function expects a subdirectory named "Assembly" with GenBank files.
    tus_df : Optional[pd.DataFrame]
        DataFrame containing TU (Transcription Unit) information, including UNS tags in the `UNS_Context` column.
    vector_rec : SeqRecord
        The backbone sequence record to be returned when the role is "Backbone".
    Returns
    -------
    Callable[[str, str], SeqRecord]
        A function that takes a construct ID (`cid`) and a role (`role`) and returns the corresponding `SeqRecord`:
        - If `role` is "Backbone", returns `vector_rec`.
        - For TU roles (e.g., "TU1", "TU2"), selects the best matching TU-final sequence from the assembly directory based on UNS tags in the DataFrame.
    - GenBank files in `<reports_dir>/Assembly` are loaded and cached for efficient lookup.
    - The resolver matches TU roles using UNS tags extracted from the DataFrame to select the most appropriate sequence.
    - Raises `RuntimeError` if the assembly directory or required records are missing.
    """
    assembly_dir = Path(reports_dir) / "Assembly"

    @cache
    def _load_all() -> list[tuple[Path, SeqRecord]]:
        out: list[tuple[Path, SeqRecord]] = []
        if assembly_dir.is_dir():
            for p in sorted(assembly_dir.glob("*.gb*")):
                try:
                    rec = SeqIO.read(str(p), "genbank")
                    rec.id = rec.name = p.stem
                    rec.annotations.setdefault("molecule_type", "DNA")
                    out.append((p, rec))
                except Exception:
                    pass
        return out

    def _uns_tags_from_context(ctx: str) -> list[str]:
        return [f"UNS{n}" for n in _UNS_NUM.findall(str(ctx))]

    def _text_pool(rec: SeqRecord) -> list[str]:
        pool = [rec.id, rec.name, getattr(rec, "description", "")]
        for f in getattr(rec, "features", []):
            for k in (
                "label",
                "note",
                "gene",
                "product",
                "locus_tag",
                "name",
                "ApEinfo_label",
            ):
                v = f.qualifiers.get(k, [])
                if isinstance(v, list):
                    pool.extend(map(str, v))
                elif v:
                    pool.append(str(v))
            pool.append(getattr(f, "type", ""))
        return [s.lower() for s in pool if s]

    def _has_uns(rec: SeqRecord, tag: str) -> bool:
        t = tag.lower()
        return any(t in s for s in _text_pool(rec))

    def _lookup(cid: str, role: str) -> SeqRecord:
        if role.strip().lower() == "backbone":
            return vector_rec

        if not assembly_dir.is_dir():
            raise RuntimeError(
                f"Assembly folder not found: {assembly_dir}. "
                "Run Golden Gate first so TU finals exist."
            )

        loaded = _load_all()
        if not loaded:
            raise RuntimeError(f"No TU finals found in {assembly_dir}.")

        tags: list[str] = []
        if tus_df is not None and not tus_df.empty:
            row = tus_df[
                (tus_df["ConstructID"].astype(str) == str(cid))
                & (tus_df["TU"].astype(str).str.upper() == str(role).upper())
            ]
            if not row.empty:
                tags = _uns_tags_from_context(str(row.iloc[0].get("UNS_Context", "")))

        best_rec: Optional[SeqRecord] = None
        best_score = -1
        for _, rec in loaded:
            score = sum(_has_uns(rec, t) for t in tags) if tags else 0
            if score > best_score:
                best_rec, best_score = rec, score

        if best_rec is None:
            raise RuntimeError(
                f"No TU record found for ({cid}, {role}) in {assembly_dir}."
            )

        return best_rec

    return _lookup
