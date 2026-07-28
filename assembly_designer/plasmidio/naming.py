"""Naming helper for final constructs."""

from __future__ import annotations

import re
from collections.abc import Sequence

import pandas as pd


def _sanitize_token(s: object) -> str:
    """
    Convert any object to a filesystem- and GenBank-id-friendly token.

    Parameters
    ----------
    s : object
        Value to sanitize (often a string; may be NaN/None).

    Returns
    -------
    str
        Sanitized token (ASCII-ish, compact underscores).
    """
    text = "" if s is None else str(s)
    text = re.sub(r"\s+", "_", text)
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text or "NA"


def _tu_label(row: pd.Series, *, include_uns: bool = True) -> str:
    """
    Build a TU label from common columns if present.

    Parameters
    ----------
    row : pd.Series
        Row from the TUs sheet.
    include_uns : bool, optional
        Include UNS/backbone context if available, by default True.

    Returns
    -------
    str
        Label like 'P45_AB_(A)_RBS_01_BC_(A)_E0030_CD_TrrnB_DE_(A)[_UNS1A_UNS3_E]'.
    """
    cols: list[str] = ["Promoter", "RBS", "Gene", "Terminator"]
    parts: list[str] = []
    for col in cols:
        if col in row and pd.notna(row[col]):
            parts.append(_sanitize_token(row[col]))
    if include_uns and "UNS_Context" in row and pd.notna(row["UNS_Context"]):
        parts.append(_sanitize_token(row["UNS_Context"]))
    return "_".join(parts) if parts else "TU"


def _role_label(
    cid: str,
    role: str,
    constructs_df: pd.DataFrame,
    tus_df: pd.DataFrame,
    *,
    include_uns: bool = True,
) -> str:
    """
    Resolve a human-readable label for one assembly role within a construct.

    Parameters
    ----------
    cid : str
        ConstructID.
    role : str
        Role name ('Backbone', 'TU1', 'TU2', ...).
    constructs_df : pd.DataFrame
        Constructs sheet.
    tus_df : pd.DataFrame
        TUs sheet.
    include_uns : bool, optional
        Include UNS context in TU labels, by default True.

    Returns
    -------
    str
        Resolved label for the role.
    """
    r = role.upper()
    if r == "BACKBONE":
        if not constructs_df.empty and "Vector" in constructs_df.columns:
            rowc = constructs_df[constructs_df["ConstructID"] == cid]
            if not rowc.empty and pd.notna(rowc["Vector"].iloc[0]):
                return _sanitize_token(rowc["Vector"].iloc[0])
        return "Backbone"

    if r.startswith("TU") and not tus_df.empty:
        mask = (tus_df["ConstructID"] == cid) & (
            tus_df["TU"].astype(str).str.upper() == r
        )
        if mask.any():
            return _tu_label(tus_df.loc[mask].iloc[0], include_uns=include_uns)

    # Fallback to the role name itself
    return _sanitize_token(role)


def build_construct_label(
    cid: str,
    order: Sequence[str],
    constructs_df: pd.DataFrame,
    tus_df: pd.DataFrame,
    *,
    prefix_with_cid: bool = True,
    include_uns_in_tu: bool = True,
) -> tuple[str, list[str]]:
    """
    Build a descriptive label for an assembled construct.

    The label format is:
        '<CID>__<BackboneLabel>__<TU1Label>__<TU2Label>...'
    (CID prefix is optional). Each per-role label is derived from sheet fields.

    Parameters
    ----------
    cid : str
        ConstructID.
    order : sequence of str
        Fragment roles in assembly order (e.g., ['Backbone', 'TU1', 'TU2']).
    constructs_df : pd.DataFrame
        Constructs sheet (should contain 'ConstructID' and optionally 'Vector').
    tus_df : pd.DataFrame
        TUs sheet (should contain 'ConstructID', 'TU', and part columns).
    prefix_with_cid : bool, optional
        Prepend the ConstructID, by default True.
    include_uns_in_tu : bool, optional
        Include 'UNS_Context' in TU labels when present, by default True.

    Returns
    -------
    tuple[str, list[str]]
        (full_label, per_role_labels). `full_label` is safe for filenames.
    """
    per_role = [
        _role_label(cid, role, constructs_df, tus_df, include_uns=include_uns_in_tu)
        for role in order
    ]
    head = [_sanitize_token(cid)] if prefix_with_cid else []
    full = "__".join(head + per_role)
    return full, per_role
