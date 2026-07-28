"""Small helpers shared by two or more plasmidio submodules.

Kept separate (rather than duplicated or owned by one peer module) to avoid
circular imports between ``batch.py``, ``plotting.py``, ``template.py`` and
``history.py``.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from typing import Any, Optional

from .alignment import summarize_feature_coverage

_UNS_NUMBER = re.compile(r"uns\s*([0-9]+)", flags=re.IGNORECASE)
_ADAPTER_2L = re.compile(r"(?:_|[\(])([A-Z]{2})(?:_|[\)])")


def _format_feature_hits_from_alignment(
    pref: Any,
    start_ref: int,
    end_ref: int,
    read_seq: str,
    strand: str,
    *,
    match: float,
    mismatch: float,
    gap_open: float,
    gap_extend: float,
) -> str:
    """
    Build the `features_hit` string with **PID per feature** using Biopython's
    pairwise2 only (no edlib).

    The function semiglobally aligns the reference window `[start_ref:end_ref]`
    against the read (reverse-complementing the read if needed), obtains aligned
    strings `(aligned_ref, aligned_read)`, and passes them to your internal
    `_feature_pid_stats(...)` helper to compute per-feature PID values. If the
    alignment cannot be produced, it falls back to coverage-based formatting via
    `summarize_feature_coverage(...)`.

    Parameters
    ----------
    pref : Any
        Plasmid reference object exposing at least `concat_ref` (str-like) and
        `feature_map`.
    start_ref : int
        Start (inclusive) of the reference slice in `pref.concat_ref`.
    end_ref : int
        End (exclusive) of the reference slice in `pref.concat_ref`.
    read_seq : str
        Read sequence (5'→3').
    strand : str
        `"F"`/`"forward"` or `"R"`/`"reverse"`. `"R"` triggers reverse complement
        of `read_seq` before alignment.
    match : float
        Match reward for `pairwise2.align.globalms`.
    mismatch : float
        Mismatch penalty for `pairwise2.align.globalms`.
    gap_open : float
        Gap-open penalty for `pairwise2.align.globalms`.
    gap_extend : float
        Gap-extend penalty for `pairwise2.align.globalms`.

    Returns
    -------
    str
        Comma-separated per-feature PID like
        ``"promoter:J23106 (100%), 5'UTR:B0033m (98%), CDS:PanD (64%)"``.
        Empty string if no overlapping features exist.

    Notes
    -----
    - Requires your helpers: `_feature_pid_stats(...)` and
      `summarize_feature_coverage(...)`.
    """
    # Reverse-complement read if necessary
    if strand.upper().startswith("R"):
        tbl = str.maketrans("ACGTNacgtn", "TGCANtgcan")
        read_seq = read_seq.translate(tbl)[::-1]

    ref_seg = str(pref.concat_ref)[int(start_ref) : int(end_ref)]
    read = read_seq.upper()

    aligned_ref: Optional[str] = None
    aligned_read: Optional[str] = None

    # Semiglobal-ish alignment using pairwise2 (free end gaps on the READ)
    try:
        from Bio import pairwise2

        alns = pairwise2.align.globalms(
            ref_seg,
            read,
            match,
            mismatch,
            gap_open,
            gap_extend,
            penalize_end_gaps=(False, True),  # free gaps at read ends
            one_alignment_only=True,
        )
        if alns:
            a = alns[0]
            aligned_ref = str(a.seqA)
            aligned_read = str(a.seqB)
    except Exception:
        aligned_ref = None
        aligned_read = None

    # Fallback: coverage-based summary if alignment failed
    if aligned_ref is None or aligned_read is None:
        try:
            hits = list(
                summarize_feature_coverage(
                    pref.feature_map, (int(start_ref), int(end_ref))
                )
            )
        except Exception:
            hits = []
        if not hits:
            return ""
        return ", ".join(
            f"{getattr(h, 'ftype', 'feat')}:{getattr(h, 'fname', '')} "
            f"({int(round(float(getattr(h, 'coverage_pct', 0))))}%)"
            for h in hits
        )

    # Compute PID per feature via your helper
    parts: list[str] = []
    try:
        stats = _feature_pid_stats(
            pref, int(start_ref), int(end_ref), aligned_ref, aligned_read
        )
        for _fs, _fe, ftype, fname, pid, denom, covered in stats:
            if denom > 0:
                parts.append(f"{ftype}:{fname} ({float(pid):.0f}%)")
            elif covered > 0:
                parts.append(f"{ftype}:{fname} (PID n/a)")
        # If no overlapping features, return empty string
        return ", ".join(parts)
    except Exception:
        return ""


def _feature_pid_stats(
    pref: Any,
    start_ref: int,
    end_ref: int,
    aln_ref: str,
    aln_read: str,
) -> list[tuple[int, int, str, str, float, int, int]]:
    """
    Compute per-feature identity inside a reference window given gapped strings.

    Parameters
    ----------
    pref : Any
        Prepared plasmid reference that exposes `feature_map` as a list of
        (start, end, ftype, fname) in concatenated coordinates.
    start_ref, end_ref : int
        Window bounds on the concatenated reference (half-open).
    aln_ref, aln_read : str
        Gapped alignment strings of equal length for the chosen window,
        using '-' as gap. `aln_ref` corresponds to ref[start_ref:end_ref].

    Returns
    -------
    list of tuple
        Each tuple is (fs, fe, ftype, fname, pid, denom_core, covered_cols) where:

        - `pid` is % identity across **core columns** (both sides consume a base)
          whose absolute reference position falls into that feature.
        - `denom_core` counts only the number of core columns within the feature
          (0 means no identity could be computed for that feature).
        - `covered_cols` counts **all** alignment columns (including gaps on one
          side) whose reference position lies in the feature.

    Notes
    -----
    - This measures *identity*, not coverage. Use `covered_cols` if you need
      coverage-related logic.
    """
    if len(aln_ref) != len(aln_read):
        raise ValueError("aln_ref and aln_read must have the same length.")

    # Map each alignment column to an absolute reference position or None
    ref_pos: list[Optional[int]] = []
    rp = int(start_ref)
    for ch in aln_ref:
        if ch == "-":
            ref_pos.append(None)
        else:
            ref_pos.append(rp)
            rp += 1

    out: list[tuple[int, int, str, str, float, int, int]] = []
    fmap: Sequence[tuple[int, int, str, str]] = getattr(pref, "feature_map", ())

    for fs, fe, ftype, fname in fmap:
        matches = 0
        denom = 0
        covered = 0
        for i, (ra, rb) in enumerate(zip(aln_ref, aln_read, strict=False)):
            p = ref_pos[i]
            if p is None or not (fs <= p < fe):
                continue
            covered += 1
            if ra != "-" and rb != "-":
                denom += 1
                if ra == rb:
                    matches += 1
        pid = 0.0 if denom == 0 else 100.0 * matches / float(denom)
        out.append(
            (
                int(fs),
                int(fe),
                str(ftype),
                str(fname),
                float(pid),
                int(denom),
                int(covered),
            )
        )

    return out
