"""Read-anchored alignment viewers/plotters (ASCII, matplotlib, HTML)."""

from __future__ import annotations

import html
import math
import re
import textwrap
from collections.abc import Iterable, Sequence
from typing import Any, Literal, Optional

import matplotlib.pyplot as plt
import pandas as pd
from Bio.Align import PairwiseAligner
from IPython.display import HTML, display
from matplotlib.figure import Figure

from ._shared import _feature_pid_stats


def norm_ref(s: str) -> str:
    """Return first contiguous digit block found in *s* ('' if none)."""
    m = re.search(r"[0-9]+", str(s))
    return m.group(0) if m else ""


def pick_rows_by_ref(
    df: pd.DataFrame, ref: str | int, *, min_pid: float = 0.0
) -> pd.DataFrame:
    """Return all rows whose `sequence_name` contains the numeric *ref*."""
    key = df["sequence_name"].astype(str).map(norm_ref)
    hits = df[key.eq(str(ref))]
    if min_pid > 0:
        hits = hits[hits["pid"] >= float(min_pid)]
    return hits.sort_values(["pid", "score"], ascending=[False, False])


def pick_row_by_ref(
    df: pd.DataFrame, ref: str | int, *, nth: int = 0, min_pid: float = 0.0
) -> Optional[pd.Series]:
    """Return the nth best row (by PID then score) for *ref* or None."""
    hits = pick_rows_by_ref(df, ref, min_pid=min_pid)
    return None if hits.empty else hits.iloc[nth]


def _pick_pref(
    plasmids: Sequence[Any], *, file_name: Optional[str], plasmid_id: Optional[str]
) -> Any:
    """Return a PlasmidRef by exact **file basename** first, else by plasmid_id.

    Notes
    -----
    - `prepare_plasmids` stores `PlasmidRef.file` as a basename.
      We therefore match `os.path.basename(file_name)`.
    - If `file_name` is provided but not found, we raise a clear error
      instead of silently falling back to a different plasmid.
    """
    import os

    # Normalize inputs
    fn = (file_name or "").strip()
    pid = (plasmid_id or "").strip()

    # 1) Prefer exact file basename match when provided
    if fn:
        base = os.path.basename(fn)
        for p in plasmids:
            if str(getattr(p, "file", "")) == base:
                return p
        avail = ", ".join(str(getattr(p, "file", "")) for p in plasmids)
        raise KeyError(f"Plasmid file not found: {base!r}. " f"Available: {avail}")

    # 2) Fallback: plasmid_id (only if no file was provided)
    if pid:
        for p in plasmids:
            if str(getattr(p, "plasmid_id", "")) == pid:
                return p
        raise KeyError(f"Plasmid with plasmid_id={pid!r} not found.")

    raise KeyError(f"Plasmid not found (file={file_name!r}, id={plasmid_id!r}).")


def _pick_read(reads: Sequence[Any], name: str) -> Any:
    """Return a read object by its `name` attribute (exact match)."""
    target = str(name)
    for r in reads:
        if str(getattr(r, "name", "")) == target:
            return r
    avail = ", ".join(str(getattr(r, "name", "")) for r in reads[:10])
    raise KeyError(f"Read not found: {target!r}. Examples: {avail}")


def plot_alignment_only(
    row: Any,  # pandas.Series-like (mapping with string keys)
    plasmids: Sequence[Any],
    reads: Sequence[Any],
    *,
    span: Literal["window", "full"] = "window",
    strand: Literal["F", "R", "auto"] = "auto",
    wrap: int = 90,
    fontsize: int = 10,
    title_width: int = 100,
    fig_width: float = 14.0,
    mid_per_block: float = 0.18,  # height per 3-line wrapped block (inches)
    line_gap: float = 0.035,  # gap between ref|match|read lines
    block_gap: float = 0.06,  # gap between wrapped blocks
) -> Figure:
    """
    Render a standalone, read-anchored ASCII alignment figure.

    The function recomputes a local alignment (Smith–Waterman) against either a
    window or the full concatenated reference, then draws wrapped monospace
    lines: reference, match bars, and read (the read is padded with gaps to
    span the chosen reference segment).

    Parameters
    ----------
    row : pandas.Series-like
        Result row with at least: ``plasmid_file`` or ``plasmid_id``,
        ``sequence_name``, ``start_ref``, ``end_ref``, ``pid``, ``core_len``,
        ``score``, ``snps``.
    plasmids : Sequence[Any]
        Prepared plasmid objects exposing: ``file``, ``plasmid_id``,
        ``concat_ref``.
    reads : Sequence[Any]
        Read objects exposing: ``name``, ``seq``.
    span : {"window", "full"}, optional
        Reference span for the alignment. Default ``"window"``.
    strand : {"F", "R", "auto"}, optional
        Orientation to use; ``"auto"`` tries both and picks the best. Default
        ``"auto"``.
    wrap : int, optional
        Characters per wrapped alignment line. Default 90.
    fontsize : int, optional
        Font size for alignment text. Default 10.
    title_width : int, optional
        Truncation width for labels in the title. Default 100.
    fig_width : float, optional
        Figure width in inches. Default 14.0.
    mid_per_block : float, optional
        Vertical space in inches per wrapped 3-line block. Default 0.18.
    line_gap : float, optional
        Vertical gap between ref/match/read lines. Default 0.035.
    block_gap : float, optional
        Vertical gap between wrapped blocks. Default 0.06.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The alignment panel as a Matplotlib figure.

    Notes
    -----
    - ``align_read_semiglobal`` should uppercase sequences internally to avoid
      case-sensitive mismatches in Biopython's aligner.
    """
    # Resolve objects without requiring pandas imports
    pref = _pick_pref(
        plasmids,
        file_name=row.get("plasmid_file"),
        plasmid_id=row.get("plasmid_id"),
    )
    read = _pick_read(reads, str(row["sequence_name"]))

    start_ref = int(row["start_ref"])
    end_ref = int(row["end_ref"])
    ref_full = str(pref.concat_ref)
    read_seq = str(read.seq)

    # Compute gapped strings via the alignment core
    a_ref, a_read, s_eff, e_eff = align_read_semiglobal(
        ref_full,
        start_ref=start_ref,
        end_ref=end_ref,
        read_seq=read_seq,
        strand=strand,
        span=span,
    )

    # Build match line and prepare wrapping
    match_line = "".join(
        "|" if (x == y and x != "-" and y != "-") else " "
        for x, y in zip(a_ref, a_read, strict=False)
    )

    def _chunks(txt: str, n: int) -> Iterable[str]:
        """Yield non-overlapping chunks of size n."""
        step = max(1, n)
        for i in range(0, len(txt), step):
            yield txt[i : i + step]

    n_blocks = max(1, math.ceil(len(a_ref) / max(1, wrap)))
    mid_height = max(1.4, n_blocks * mid_per_block)  # enforce a usable minimum

    fig = plt.figure(figsize=(fig_width, mid_height))
    ax = fig.add_subplot(111)

    y = 1.0
    for rline, mline, qline in zip(
        _chunks(a_ref, wrap),
        _chunks(match_line, wrap),
        _chunks(a_read, wrap),
        strict=False,
    ):
        ax.text(0.01, y, rline, family="monospace", va="top", fontsize=fontsize)
        y -= line_gap
        ax.text(0.01, y, mline, family="monospace", va="top", fontsize=fontsize)
        y -= line_gap
        ax.text(0.01, y, qline, family="monospace", va="top", fontsize=fontsize)
        y -= block_gap

    ax.set_ylim(0, 1.02)
    ax.set_axis_off()

    # Concise header
    name_pl = textwrap.shorten(
        str(getattr(pref, "file", "")), width=title_width, placeholder="…"
    )
    name_rd = textwrap.shorten(
        str(getattr(read, "name", "")), width=title_width, placeholder="…"
    )
    title = (
        f"{name_pl} [{s_eff}–{e_eff}]  vs  {name_rd}  (strand {strand}; span={span})\n"
        f"PID={float(row['pid']):.2f}%  core={int(row['core_len'])}  "
        f"score={float(row['score']):.1f}  SNPs={int(row['snps'])}"
    )
    ax.set_title(title, pad=6)

    return fig


def plot_feature_track_only(
    pref: Any,
    *,
    start_ref: int,
    end_ref: int,
    fstats: list[tuple[int, int, str, str, float, int, int]],
    fig_width: float = 12.0,
    height: float = 0.95,
    track_height: float = 0.34,
    rect_height: float = 0.20,
    lane_gap: float = 0.10,
    label_thresh: float = 0.06,
    fontsize: int = 10,
    label_bbox: bool = True,
    label_bbox_alpha: float = 0.9,
    zoom_to_features: bool = False,
    pad_bp: int = 120,
    verbose: bool = False,
) -> Figure:
    """
    Render a standalone feature-track panel for a reference window.

    Colors encode identity:
      • PID ≥ 99.5% → green (treated as 100%)
      • otherwise   → yellow
      • only zero-width artifacts are gray

    When no coverage stats are available, the function falls back to
    ``pref.feature_map`` (drawn as yellow boxes since PID is unknown).

    Parameters
    ----------
    pref : Any
        Plasmid reference exposing ``feature_map`` as
        ``[(start, end, type, name), ...]`` in concatenated coordinates.
    start_ref, end_ref : int
        Window bounds (half-open) within the concatenated reference.
    fstats : list of tuple
        Feature stats from ``_feature_pid_stats``:
        ``(fs, fe, ftype, fname, pid, denom, covered_bp)``.
    fig_width : float, optional
        Figure width in inches. Default 12.0.
    height : float, optional
        Figure height in inches. Default 0.95.
    track_height : float, optional
        Vertical space for the feature lanes (axis fraction). Default 0.34.
    rect_height : float, optional
        Rectangle height for a feature (axis fraction). Default 0.20.
    lane_gap : float, optional
        Vertical gap between lanes. Default 0.10.
    label_thresh : float, optional
        Minimum relative width to draw a label. Default 0.06.
    fontsize : int, optional
        Font size for labels. Default 10.
    label_bbox : bool, optional
        Draw a white rounded label background. Default True.
    label_bbox_alpha : float, optional
        Opacity for the label background. Default 0.9.
    zoom_to_features : bool, optional
        If True, x-limits are zoomed to the union of visible features with
        padding ``pad_bp``. Default False.
    pad_bp : int, optional
        Padding (bp) used only when ``zoom_to_features`` is True. Default 120.
    verbose : bool, optional
        Print placement diagnostics. Default False.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Feature-track figure.
    """

    def _color_from_pid(pid: float, denom: int, width_bp: int) -> str:
        """
        PID-based color rule:
          - gray only if nothing drawable (width<=0)
          - yellow if no PID data (denom==0) or PID < 99.5
          - green if PID >= 99.5
        """
        if width_bp <= 0:
            return "#bdbdbd"  # gray: degenerate box
        if denom == 0:
            return "#fdd835"  # yellow: unknown PID
        return "#2e7d32" if pid >= 99.5 else "#fdd835"

    win_l, win_r = int(start_ref), int(end_ref)
    span_bp = max(1, win_r - win_l)

    # Keep only stats overlapping the window
    def _visible(stats: Iterable[tuple[int, int, str, str, float, int, int]]):
        out: list[tuple[int, int, str, str, float, int, int]] = []
        for item in stats:
            fs, fe = item[0], item[1]
            if min(win_r, fe) > max(win_l, fs):
                out.append(item)
        return out

    items = _visible(fstats)

    # Fallback: raw feature_map as "unknown PID" (yellow)
    if not items:
        raw: list[tuple[int, int, str, str, float, int, int]] = []
        for fs, fe, ftype, fname in getattr(pref, "feature_map", []):
            if min(win_r, fe) > max(win_l, fs):
                raw.append((fs, fe, ftype, fname, 0.0, 0, 0))  # denom=0 => unknown PID
        items = raw

    # Optional zoom to visible features (+ padding)
    x0, x1 = win_l, win_r
    if zoom_to_features and items:
        x0 = max(win_l, min(fs for fs, *_ in items) - int(pad_bp))
        x1 = min(win_r, max(fe for _, fe, *_ in items) + int(pad_bp))
        span_bp = max(1, x1 - x0)

    def _rel(u: int) -> float:
        return (u - x0) / span_bp

    # Greedy lane packing
    lanes: list[list[tuple[float, float]]] = []
    placed: list[tuple[float, float, int, int, str, str, float, int, int, int]] = []
    for fs, fe, ftype, fname, pid, denom, covered in sorted(items, key=lambda t: t[0]):
        left, right = max(x0, fs), min(x1, fe)
        xr0, xr1 = _rel(left), _rel(right)
        if xr1 <= xr0:
            continue
        lane_idx = 0
        for li, lane in enumerate(lanes):
            if not lane or xr0 >= lane[-1][1]:
                lane_idx = li
                lane.append((xr0, xr1))
                break
        else:
            lane_idx = len(lanes)
            lanes.append([(xr0, xr1)])
        placed.append(
            (
                xr0,
                xr1,
                fs,
                fe,
                ftype,
                fname,
                float(pid),
                int(denom),
                int(covered),
                lane_idx,
            )
        )

    if verbose:
        print(
            f"[feature-panel] window {win_l}-{win_r} | "
            f"coverage_items={len(fstats)} visible={len(_visible(fstats))} "
            f"raw_map={len(getattr(pref, 'feature_map', []))} placed={len(placed)}"
        )

    # --- plotting -------------------------------------------------------------
    fig = plt.figure(figsize=(fig_width, height))
    ax = fig.add_subplot(111)

    ax.hlines(0.5, 0, 1, color="black", linewidth=0.8)  # baseline

    if placed:
        nlanes = max(1, len(lanes))
        lane_step = track_height / nlanes
        h = min(rect_height, max(0.02, lane_step - lane_gap))
        base_y = 0.5 - 0.5 * track_height + (lane_step - h) / 2.0

        for (
            xr0,
            xr1,
            fs,
            fe,
            ftype,
            fname,
            pid,
            denom,
            covered,  # noqa: B007
            li,
        ) in placed:  # noqa: B007
            y = base_y + (nlanes - 1 - li) * lane_step
            width_rel = max(1e-6, xr1 - xr0)
            width_bp = max(0, min(x1, fe) - max(x0, fs))

            ax.add_patch(
                plt.Rectangle(
                    (xr0, y),
                    width_rel,
                    h,
                    facecolor=_color_from_pid(pid, denom, width_bp),
                    edgecolor="black",
                    linewidth=0.8,
                )
            )

            if width_rel >= float(label_thresh):
                label = (
                    f"{ftype}:{fname}"
                    if denom == 0
                    else f"{ftype}:{fname} (PID {pid:.0f}%)"
                )
                bbox_cfg: Optional[dict[str, Any]] = (
                    dict(
                        boxstyle="round,pad=0.2",
                        fc="white",
                        ec="none",
                        alpha=float(label_bbox_alpha),
                    )
                    if label_bbox
                    else None
                )
                ax.text(
                    xr0 + width_rel / 2.0,
                    y + h * 0.9,
                    label,
                    ha="center",
                    va="bottom",
                    fontsize=fontsize,
                    bbox=bbox_cfg,
                )

        ax.set_title("Feature context (green=PID 100%; yellow<100%)", pad=6)
    else:
        ax.text(
            0.5,
            0.65,
            "No features overlap this window",
            ha="center",
            va="center",
            fontsize=fontsize + 1,
        )
        ax.set_title("Feature context (no overlapping features)", pad=6)

    ax.set(xlim=(0, 1), ylim=(0, 1), yticks=[], xticks=[0, 1])
    ax.set_xticklabels([str(x0), str(x1)])
    return fig


def plot_alignment(
    row: Any,  # pandas.Series-like (mapping with string keys)
    plasmids: Sequence[Any],
    reads: Sequence[Any],
    *,
    # --- alignment (same knobs as before) ---
    span: Literal["window", "full"] = "window",
    strand: Literal["F", "R", "auto"] = "auto",
    wrap: int = 90,
    fontsize: int = 10,
    title_width: int = 100,
    fig_width: float = 14.0,
    mid_per_block: float = 0.18,
    line_gap: float = 0.035,
    block_gap: float = 0.06,
    # --- feature track control (separate second figure) ---
    feature_span: Optional[Literal["window", "full"]] = None,
    feature_track_kwargs: Optional[dict[str, Any]] = None,
) -> tuple[Figure, Figure]:
    """
    Render two figures for a result row: an ASCII alignment and a feature track.

    The function first computes a local alignment (Smith–Waterman) for the chosen
    `span` (window or full concat) to draw the ASCII alignment. If the feature
    panel span differs (`feature_span != span`), the alignment is recomputed for
    the feature span so per-feature PID/coverage are derived from the correct
    sequences (avoids gray boxes).

    Parameters
    ----------
    row : pandas.Series-like
        One row from your results table. Must provide:
        ``plasmid_file/plasmid_id``, ``sequence_name``, ``start_ref``, ``end_ref``,
        ``pid``, ``core_len``, ``score``, ``snps``.
    plasmids : Sequence[Any]
        Prepared plasmid references exposing at least ``file``, ``plasmid_id``,
        ``concat_ref``, ``feature_map``.
    reads : Sequence[Any]
        Read objects exposing at least ``name``, ``seq``.
    span : {"window", "full"}, optional
        Alignment span for the ASCII panel. Default ``"window"``.
    strand : {"F", "R", "auto"}, optional
        Orientation to use; default ``"auto"`` chooses the best.
    wrap, fontsize, title_width, fig_width, mid_per_block, line_gap, block_gap
        Layout controls passed to the ASCII plotter.
    feature_span : {"window", "full"}, optional
        Span used for the feature panel. If ``None``, falls back to ``span``.
        When ``"window"``, the reported window ``[start_ref:end_ref]`` is used.
    feature_track_kwargs : dict, optional
        Additional kwargs for ``plot_feature_track_only``.

    Returns
    -------
    (fig_alignment, fig_features) : tuple[Figure, Figure]
        The ASCII alignment figure and the feature-track figure.

    Notes
    -----
    - Recomputing the alignment when spans differ ensures PID coloring uses
      sequences aligned over the same coordinates shown in the feature panel.
    """
    # --- resolve inputs (no pandas import required) --------------------------
    pref = _pick_pref(
        plasmids,
        file_name=row.get("plasmid_file"),
        plasmid_id=row.get("plasmid_id"),
    )
    rd = _pick_read(reads, str(row["sequence_name"]))

    start_ref = int(row["start_ref"])
    end_ref = int(row["end_ref"])
    ref_concat = str(pref.concat_ref)
    read_seq = str(rd.seq)

    # --- 1) alignment for the ASCII panel ------------------------------------
    a_ref, a_read, s_eff, e_eff = align_read_semiglobal(
        ref_concat,
        start_ref=start_ref,
        end_ref=end_ref,
        read_seq=read_seq,
        strand=strand,
        span=span,
    )

    fig_align: Figure = plot_alignment_only(
        row,
        plasmids,
        reads,
        span=span,
        strand=strand,
        wrap=wrap,
        fontsize=fontsize,
        title_width=title_width,
        fig_width=fig_width,
        mid_per_block=mid_per_block,
        line_gap=line_gap,
        block_gap=block_gap,
    )

    # --- 2) feature-track alignment (recompute if spans differ) --------------
    span_for_feat: Literal["window", "full"] = feature_span or span
    if span_for_feat == "full":
        s_feat, e_feat = 0, len(ref_concat)
    else:
        s_feat, e_feat = start_ref, end_ref

    # If the panel spans differ, recompute alignment for the feature span
    a_ref_feat, a_read_feat = a_ref, a_read
    if span_for_feat != span:
        a_ref_feat, a_read_feat, _, _ = align_read_semiglobal(
            ref_concat,
            start_ref=start_ref,
            end_ref=end_ref,
            read_seq=read_seq,
            strand=strand,
            span=span_for_feat,
        )

    # Compute per-feature PID/coverage for that exact span
    fstats = _feature_pid_stats(pref, int(s_feat), int(e_feat), a_ref_feat, a_read_feat)

    # Defaults for the feature panel; allow overrides
    ft_kwargs: dict[str, Any] = {
        "fig_width": fig_width,
        "height": 0.95,
        "track_height": 0.34,
        "rect_height": 0.20,
        "lane_gap": 0.10,
        "label_thresh": 0.06,
        "label_bbox": True,
        "label_bbox_alpha": 0.90,
        "zoom_to_features": False,
        "pad_bp": 120,
    }
    if feature_track_kwargs:
        ft_kwargs.update(feature_track_kwargs)

    fig_feat: Figure = plot_feature_track_only(
        pref,
        start_ref=int(s_feat),
        end_ref=int(e_feat),
        fstats=fstats,
        **ft_kwargs,
    )

    _ = (s_eff, e_eff)  # quieten linters if not otherwise used
    return fig_align, fig_feat


def plot_alignment_by_ref(
    ref: str | int,
    df: pd.DataFrame,
    plasmids: Sequence[Any],
    reads: Sequence[Any],
    *,
    nth: int = 0,
    min_pid: float = 0.0,
    **plot_kwargs: Any,
) -> Optional[pd.Series]:
    """
    Plot the alignment for the *nth* best row matching a numeric reference ID.

    This is a thin convenience wrapper around :func:`plot_alignment`. It looks up
    a row in the results table by a numeric identifier contained in
    ``sequence_name`` (e.g., ``2034`` for ``EF73802034_...``), then forwards all
    remaining keyword arguments to :func:`plot_alignment`.

    Parameters
    ----------
    ref : str or int
        Numeric identifier to search for in ``df['sequence_name']``.
    df : pandas.DataFrame
        Results table with the columns required by ``plot_alignment``.
    plasmids : sequence
        Prepared plasmid refs (from ``prepare_plasmids``).
    reads : sequence
        Loaded reads (from ``load_reads``).
    nth : int, optional
        0-based rank among matching rows (``0`` = best), by default 0.
    min_pid : float, optional
        Minimum PID (%) to consider when selecting candidates, by default 0.0.
    **plot_kwargs
        Forwarded to :func:`plot_alignment`. For backward compatibility,
        the legacy aliases ``feature_mode`` and ``draw_features`` are accepted
        and mapped to the new ``span`` parameter (``"window"`` or ``"full"``).

    Returns
    -------
    pandas.Series or None
        The plotted row, or ``None`` if nothing matched.
    """
    # Find the row (helper expected to exist in your module)
    row = pick_row_by_ref(df, ref, nth=nth, min_pid=min_pid)
    if row is None:
        print(f"No rows for ref #{ref} (min_pid={min_pid}, nth={nth}).")
        return None

    # ── Backward-compat: map legacy kwargs to the new `span` ─────────────────
    # Accept `feature_mode` (old) or `draw_features` (intermediate) as aliases.
    if "span" not in plot_kwargs:
        if "draw_features" in plot_kwargs:
            plot_kwargs["span"] = plot_kwargs.pop("draw_features")
        elif "feature_mode" in plot_kwargs:
            plot_kwargs["span"] = plot_kwargs.pop("feature_mode")

    # Sensible default if the caller didn’t specify a span
    plot_kwargs.setdefault("span", "window")

    print(f"Showing ref #{ref} — rank {nth + 1}, PID={row['pid']:.2f}%")
    plot_alignment(row, plasmids, reads, **plot_kwargs)
    return row


def _find_plasmid_by_file(plasmids: Sequence[Any], file_name: str) -> Any:
    """Return the prepared plasmid whose `.file` matches *file_name*."""
    for p in plasmids:
        if getattr(p, "file", None) == file_name:
            return p
    raise KeyError(f"Plasmid with file={file_name!r} not found.")


def _find_read_by_name(reads: Sequence[Any], name: str) -> Any:
    """Return the prepared read whose `.name` matches *name*."""
    for r in reads:
        if getattr(r, "name", None) == name:
            return r
    raise KeyError(f"Read with name={name!r} not found.")


def _revcomp_local(seq: str) -> str:
    """Reverse-complement A/C/G/T/N (case-insensitive)."""
    tbl = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(tbl)[::-1]


def _resolve_span(
    ref_full: str,
    start_ref: Optional[int],
    end_ref: Optional[int],
    span: Literal["full", "window"] | tuple[str, dict[str, Any]],
) -> tuple[str, int, int]:
    """Return (reference_segment, s_eff, e_eff) for window/full (+ optional padding)."""
    pad_bp = 0
    if isinstance(span, tuple):
        mode, opts = span
        pad_bp = int(opts.get("pad_bp", 0))
    else:
        mode = span

    if mode == "full":
        return ref_full, 0, len(ref_full)

    s0 = int(start_ref or 0)
    e0 = int(end_ref or len(ref_full))
    s_eff = max(0, s0 - pad_bp)
    e_eff = min(len(ref_full), e0 + pad_bp)
    return ref_full[s_eff:e_eff], s_eff, e_eff


def _gapped_from_blocks(
    ref_seg: str,
    query: str,
    ref_blocks: list[tuple[int, int]],
    qry_blocks: list[tuple[int, int]],
) -> tuple[str, str]:
    """Rebuild gapped strings from PairwiseAligner block coordinates."""
    out_r: list[str] = []
    out_q: list[str] = []
    r_cur = 0
    q_cur = 0
    for (rs, re), (qs, qe) in zip(ref_blocks, qry_blocks, strict=False):  # noqa: F402
        # gaps before next block
        if rs > r_cur:
            out_r.append(ref_seg[r_cur:rs])
            out_q.append("-" * (rs - r_cur))
        if qs > q_cur:
            out_r.append("-" * (qs - q_cur))
            out_q.append(query[q_cur:qs])
        # block
        out_r.append(ref_seg[rs:re])
        out_q.append(query[qs:qe])
        r_cur, q_cur = re, qe
    # tails
    if r_cur < len(ref_seg):
        out_r.append(ref_seg[r_cur:])
        out_q.append("-" * (len(ref_seg) - r_cur))
    if q_cur < len(query):
        out_r.append("-" * (len(query) - q_cur))
        out_q.append(query[q_cur:])
    a_ref = "".join(out_r)
    a_read = "".join(out_q)
    L = min(len(a_ref), len(a_read))  # noqa: N806
    return a_ref[:L], a_read[:L]


AlignSpan = Literal["window", "full"]


def align_read_semiglobal(
    ref_full: str,
    *,
    start_ref: Optional[int],
    end_ref: Optional[int],
    read_seq: str,
    strand: Literal["F", "R", "auto"] = "auto",
    span: AlignSpan | tuple[str, dict[str, Any]] = "window",
    match: float = 2.0,
    mismatch: float = -1.0,
    gap_open: float = -6.0,
    gap_extend: float = -1.0,
    anchor_pad_bp: int = 300,  # kept for backwards-compatibility; not used here
) -> tuple[str, str, int, int]:
    """
    Perform a robust, read-anchored alignment using **local** (Smith–Waterman) scoring,
    then pad to the requested reference span for a continuous, readable ASCII layout.

    This function normalizes both reference and read to uppercase (case-sensitive
    matching in Biopython), selects the best read orientation automatically by
    default, reconstructs gapped strings from `Alignment.aligned` blocks, and pads
    the read with gaps so that the full chosen reference segment is covered.

    Parameters
    ----------
    ref_full : str
        The full reference sequence (e.g., feature-concatenated or full plasmid).
        It will be uppercased internally to avoid case-mismatch artefacts.
    start_ref : int, optional
        Start coordinate of the hit on the reference (0-based, half-open) used
        when `span="window"`. Ignored for `span="full"`.
    end_ref : int, optional
        End coordinate of the hit on the reference (0-based, half-open) used
        when `span="window"`. Ignored for `span="full"`.
    read_seq : str
        Raw read sequence to align. Reverse complement will be considered if
        `strand="auto"` (default).
    strand : {"F", "R", "auto"}, default "auto"
        Orientation to use:
        - "F": forward as provided,
        - "R": reverse complement,
        - "auto": try both and pick the higher-scoring alignment.
    span : {"window", "full"} or tuple, default "window"
        Reference segment to align against:
        - "window": use `[start_ref:end_ref]` (with optional padding),
        - "full": use the entire `ref_full`.
        You may also pass a tuple `("window", {"pad_bp": int})` to add symmetric
        padding (in bp) around the window to provide local context.
    match : float, default 2.0
        Match score for the PairwiseAligner.
    mismatch : float, default -1.0
        Mismatch score for the PairwiseAligner.
    gap_open : float, default -6.0
        Gap open penalty for the PairwiseAligner.
    gap_extend : float, default -1.0
        Gap extension penalty for the PairwiseAligner.
    anchor_pad_bp : int, default 300
        Deprecated / reserved for compatibility. Not used by this function.
        Use `span=("window", {"pad_bp": ...})` instead.

    Returns
    -------
    a_ref : str
        Gapped reference string covering the chosen segment.
    a_read : str
        Gapped read string, padded with '-' so it spans the same length as `a_ref`.
    s_eff : int
        Effective start coordinate of the reference segment within `ref_full`.
    e_eff : int
        Effective end coordinate of the reference segment within `ref_full`.

    Notes
    -----
    - Biopython's `PairwiseAligner` is case-sensitive. Both sequences are uppercased
      here to ensure that identical bases match as expected.
    - If the local alignment yields no blocks (extremely rare), a zero-width
      dummy block is injected to produce well-formed padded strings.

    Examples
    --------
    >>> a_ref, a_read, s, e = align_read_semiglobal(
    ...     ref_full, start_ref=1000, end_ref=1200, read_seq=read, span=("window", {"pad_bp": 25})
    ... )
    >>> print(a_ref)   # doctest: +ELLIPSIS
    >>> print(a_read)  # doctest: +ELLIPSIS
    """
    # Normalize reference and read sequences to uppercase to avoid case-mismatch artefacts.
    ref_full = str(ref_full).upper()
    qF = str(read_seq).upper()  # noqa: N806
    qR = _revcomp_local(qF)  # noqa: N806

    # Resolve the reference segment to align against (window/full + optional padding).
    ref_seg, s_eff, e_eff = _resolve_span(ref_full, start_ref, end_ref, span)

    # Prepare candidate strands based on requested orientation.
    candidates: list[tuple[str, str]] = (
        [("F", qF)]
        if strand == "F"
        else [("R", qR)] if strand == "R" else [("F", qF), ("R", qR)]  # auto: try both
    )

    # Configure a local (Smith–Waterman) aligner.
    al = PairwiseAligner()
    al.mode = "local"
    al.match_score = float(match)
    al.mismatch_score = float(mismatch)
    al.open_gap_score = float(gap_open)
    al.extend_gap_score = float(gap_extend)

    # Evaluate candidates and keep the highest-scoring alignment.
    best: Optional[tuple[float, str, Any]] = None  # (score, strand, alignment)
    for st, q in candidates:
        aln = next(iter(al.align(ref_seg, q)))
        sc = float(aln.score)
        if best is None or sc > best[0]:
            best = (sc, st, aln)
    assert best is not None
    _, chosen, aln = best
    q_seq = qF if chosen == "F" else qR

    # Reconstruct gapped strings from alignment blocks; pad to cover the segment.
    ref_blocks = [(int(s), int(e)) for (s, e) in aln.aligned[0]]
    qry_blocks = [(int(s), int(e)) for (s, e) in aln.aligned[1]]

    # Edge case: if no blocks (very rare), inject a zero-width block to get clean padding.
    if not ref_blocks or not qry_blocks:
        ref_blocks = [(0, 0)]
        qry_blocks = [(0, 0)]

    a_ref, a_read = _gapped_from_blocks(ref_seg, q_seq, ref_blocks, qry_blocks)

    # Ensure both strings have identical length (truncate to the shorter just in case).
    L = min(len(a_ref), len(a_read))  # noqa: N806
    return a_ref[:L], a_read[:L], s_eff, e_eff


def view_alignment(
    row: pd.Series,
    plasmids: Sequence[Any],
    reads: Sequence[Any],
    *,
    span: Literal["window", "full"] = "window",
    wrap: int = 100,
    match: float = 1.0,
    mismatch: float = -1.0,
    gap_open: float = -1.5,
    gap_extend: float = -0.5,
) -> None:
    """
    Pretty-print a read-anchored alignment for a single result row.

    This uses `align_read_semiglobal` (PairwiseAligner only) to reconstruct
    gapped strings, then prints a wrapped ASCII view.

    Parameters
    ----------
    row : pandas.Series
        One row from your results table (needs keys:
        plasmid_file, sequence_name, strand, start_ref, end_ref, pid, core_len, score, snps).
    plasmids : list[Any]
        Prepared plasmid objects (must provide `.file`, `.concat_ref`).
    reads : list[Any]
        Loaded reads (must provide `.name`, `.seq`).
    span : {"window","full"}, default "window"
        Align only the reported hit ("window") or the full concat-band ("full").
    wrap : int, default 100
        Width of each wrapped alignment line.
    match, mismatch, gap_open, gap_extend : float
        Scoring parameters forwarded to `align_read_semiglobal`.

    Returns
    -------
    None
        Prints to stdout.
    """

    def _find_plasmid_by_file(objs: Sequence[Any], fname: str) -> Any:
        for p in objs:
            if str(getattr(p, "file", "")) == fname:
                return p
        raise ValueError(f"Plasmid with file '{fname}' not found.")

    def _find_read_by_name(objs: Sequence[Any], name: str) -> Any:
        for r in objs:
            if str(getattr(r, "name", "")) == name:
                return r
        raise ValueError(f"Read with name '{name}' not found.")

    # Resolve objects and fields
    pref = _find_plasmid_by_file(plasmids, str(row["plasmid_file"]))
    rd = _find_read_by_name(reads, str(row["sequence_name"]))
    s0, e0 = int(row["start_ref"]), int(row["end_ref"])
    strand: Literal["F", "R"] = "R" if str(row.get("strand", "F")) == "R" else "F"

    # Compute gapped strings via the core
    a_ref, a_read, s_eff, e_eff = align_read_semiglobal(
        str(pref.concat_ref),
        start_ref=s0,
        end_ref=e0,
        read_seq=str(rd.seq),
        strand=strand,
        span=span,
        match=match,
        mismatch=mismatch,
        gap_open=gap_open,
        gap_extend=gap_extend,
    )

    # Build match line
    match_line = "".join(
        "|" if (x == y and x != "-" and y != "-") else " "
        for x, y in zip(a_ref, a_read, strict=False)
    )

    # Header
    left = textwrap.shorten(str(getattr(pref, "file", "")), width=96, placeholder="…")
    right = textwrap.shorten(str(getattr(rd, "name", "")), width=96, placeholder="…")
    print(
        f"{left} [{s_eff}–{e_eff}]  vs  {right} (strand {strand}; span={span})\n"
        f"PID={float(row['pid']):.2f}%  core={int(row['core_len'])}  "
        f"score={float(row['score']):.1f}  SNPs={int(row['snps'])}"
    )

    # Wrapped output
    def _chunks(t: str, n: int) -> Iterable[str]:
        for i in range(0, len(t), n):
            yield t[i : i + n]

    for rline, mline, qline in zip(
        _chunks(a_ref, wrap),
        _chunks(match_line, wrap),
        _chunks(a_read, wrap),
        strict=False,
    ):
        print(rline)
        print(mline)
        print(qline)
        print()


def view_alignment_by_ref(
    ref: str | int,
    df: pd.DataFrame,
    plasmids: Sequence[Any],
    reads: Sequence[Any],
    *,
    nth: int = 0,
    min_pid: float = 0.0,
    **kwargs: Any,
) -> Optional[pd.Series]:
    """
    Convenience wrapper: pick the *nth* best row for a numeric ref and print it.

    Parameters
    ----------
    ref
        Numeric ID contained in ``sequence_name`` (e.g. ``2034`` for
        ``EF73802034_*``).
    df
        Results dataframe from ``results_to_dataframe``.
    nth
        Rank among matches (0 = best). Default is 0.
    min_pid
        Filter out rows with ``pid < min_pid`` before ranking. Default 0.0.
    **kwargs
        Forwarded to ``view_alignment`` (e.g. custom scoring).

    Returns
    -------
    pd.Series or None
        The row that was displayed, or ``None`` if no match.
    """
    key = (
        df["sequence_name"]
        .astype(str)
        .str.extract(r"([0-9]+)", expand=False)
        .fillna("")
    )
    hits = df[key.eq(str(ref))]
    if min_pid > 0:
        hits = hits[hits["pid"] >= float(min_pid)]
    hits = hits.sort_values(["pid", "score"], ascending=[False, False])

    if hits.empty:
        display(HTML(f"<b>No rows for ref #{html.escape(str(ref))}.</b>"))
        return None

    row = hits.iloc[min(nth, len(hits) - 1)]
    view_alignment(row, plasmids, reads, **kwargs)
    return row


def view_alignment_html(
    row,  # pandas.Series (deliberately typed in a casual style)
    plasmids: Sequence[Any],
    reads: Sequence[Any],
    *,
    span: Literal["window", "full"] | tuple[str, dict[str, Any]] = "window",
    wrap: int = 100,
    theme: str = "dark",  # "dark" or "light"
    font_px: int = 13,
    # Scores are passed on to the core (Biopython PairwiseAligner)
    match: float = 1.0,
    mismatch: float = -1.0,
    gap_open: float = -1.5,
    gap_extend: float = -0.5,
    # --- static figure export (optional) ---
    save_prefix: Optional[str] = None,
    save_font_size: float = 7.0,
    save_chars_per_row: int = 60,
    save_fig_width_mm: float = 180.0,
    save_dpi: int = 300,
) -> None:
    """
    Render a colored, read-anchored HTML alignment (with a simple position ruler)
    using the **Biopython PairwiseAligner** via our helper
    :func:`align_read_semiglobal`.

    Parameters
    ----------
    row, plasmids, reads
        Row from the results dataframe plus prepared plasmids / loaded reads.
        The row must provide: ``plasmid_file`` (or ``plasmid_id``),
        ``sequence_name``, ``strand``, ``start_ref``, ``end_ref``,
        ``pid``, ``core_len``, ``score``, ``snps``.
    span : {"window","full"} or ("window", {"pad_bp":int, ...})
        Alignment span. ``"window"`` anchors the read to the reported
        ``[start_ref:end_ref]``; ``"full"`` aligns against the full concat band.
        As a tuple you can pass padding options when using ``"window"``.
    wrap : int
        Number of alignment columns per printed row.
    theme : {"dark","light"}
        Color scheme for the HTML block.
    font_px : int
        Base font size (CSS px).
    match, mismatch, gap_open, gap_extend : float
        Scoring parameters forwarded to the PairwiseAligner inside the core.
    save_prefix : str, optional
        If given, export a static Matplotlib figure as
        ``{save_prefix}.png`` and ``{save_prefix}.pdf``.
        The HTML output in Jupyter is unaffected.
    save_font_size : float
        Font size in pt for the exported figure (default 7.0).
        Increase for larger, more readable text; decrease to fit more
        alignment columns on a single page.
    save_chars_per_row : int
        Alignment columns per row in the exported figure (default 60).
        Fewer columns yield larger characters.
    save_fig_width_mm : float
        Figure width in millimetres (default 180 = full A4 text width).
    save_dpi : int
        Raster resolution for PNG export (default 300).

    Notes
    -----
    - This function **no longer** uses edlib/pairwise2 directly. All gapped
      strings are produced by :func:`align_read_semiglobal` to ensure the same
      anchoring semantics as the rest of the pipeline.
    - When ``save_prefix`` is set, a Matplotlib figure is rendered and saved
      alongside the HTML display. The export uses only standard library and
      Matplotlib; no external binaries are required.
    """

    # ---- resolve prepared objects (support both helper names) ----------------
    try:
        pref = _pick_pref(
            plasmids,
            file_name=row.get("plasmid_file"),
            plasmid_id=row.get("plasmid_id"),
        )
    except NameError:
        pref = _find_plasmid_by_file(plasmids, str(row["plasmid_file"]))

    try:
        rd = _pick_read(reads, str(row["sequence_name"]))
    except NameError:
        rd = _find_read_by_name(reads, str(row["sequence_name"]))

    # ---- window + orientation ------------------------------------------------
    s0, e0 = int(row["start_ref"]), int(row["end_ref"])
    strand: Literal["F", "R"] = "R" if str(row.get("strand", "F")) == "R" else "F"
    ref_full = str(pref.concat_ref)
    read_seq = str(rd.seq)

    # ---- align via the shared core ------------------------------------------
    a_ref, a_read, s_eff, e_eff = align_read_semiglobal(
        ref_full,
        start_ref=s0,
        end_ref=e0,
        read_seq=read_seq,
        strand=strand,
        span=span,
        match=match,
        mismatch=mismatch,
        gap_open=gap_open,
        gap_extend=gap_extend,
    )

    # ---- theme colors --------------------------------------------------------
    if theme.lower() == "light":
        css_bg, css_fg = "#ffffff", "#222222"
        c_match, c_mis, c_gap, c_tick = "#81c784", "#ff8a80", "#bdbdbd", "#555555"
    else:
        css_bg, css_fg = "#111111", "#eeeeee"
        c_match, c_mis, c_gap, c_tick = "#1b5e20", "#b71c1c", "#616161", "#bbbbbb"

    # Small helpers to paint per character
    def _span_ref(x: str, y: str) -> str:
        bg = c_gap if (x == "-" or y == "-") else (c_match if x == y else c_mis)
        return f"<span style='background:{bg};padding:0 1px'>{html.escape(x)}</span>"

    def _span_read(x: str, y: str) -> str:
        bg = c_gap if (x == "-" or y == "-") else (c_match if x == y else c_mis)
        return f"<span style='background:{bg};padding:0 1px'>{html.escape(y)}</span>"

    # Cumulative counts of non-gaps → absolute coordinates
    ref_cum: list[int] = [0]
    read_cum: list[int] = [0]
    rc_ref, rc_read = 0, 0
    for i in range(len(a_ref)):
        if a_ref[i] != "-":
            rc_ref += 1
        if a_read[i] != "-":
            rc_read += 1
        ref_cum.append(rc_ref)
        read_cum.append(rc_read)

    # Build HTML lines with simple ruler (· every base, | at 5, number at 10)
    rows_html: list[str] = []
    for i0 in range(0, len(a_ref), wrap):
        r = a_ref[i0 : i0 + wrap]
        q = a_read[i0 : i0 + wrap]

        # absolute ref start/end visible in this chunk
        ref_start_abs: Optional[int] = None
        ref_end_abs: Optional[int] = None
        for j, ch in enumerate(r):
            if ch != "-":
                ref_start_abs = s_eff + ref_cum[i0 + j]
                break
        for j in range(len(r) - 1, -1, -1):
            if r[j] != "-":
                ref_end_abs = s_eff + ref_cum[i0 + j]
                break

        # read coordinates within aligned region (0-based)
        read_start: Optional[int] = None
        read_end: Optional[int] = None
        for j, ch in enumerate(q):
            if ch != "-":
                read_start = read_cum[i0 + j]
                break
        for j in range(len(q) - 1, -1, -1):
            if q[j] != "-":
                read_end = read_cum[i0 + j]
                break

        # position ruler over reference
        ruler_cells: list[str] = []
        for j, ch in enumerate(r):
            if ch == "-":
                ruler_cells.append("&nbsp;")
                continue
            pos = s_eff + ref_cum[i0 + j]
            if pos % 10 == 0:
                lab = str(pos)
                ruler_cells.append(f"<span style='color:{c_tick}'>{lab}</span>")
                # occupy characters so layout stays monospace
                for _ in range(len(lab) - 1):
                    ruler_cells.append("")
            elif pos % 5 == 0:
                ruler_cells.append(f"<span style='color:{c_tick}'>|</span>")
            else:
                ruler_cells.append("&middot;")
        ruler_line = "".join(ruler_cells)

        line_ref = "".join(_span_ref(x, y) for x, y in zip(r, q, strict=False))
        line_read = "".join(_span_read(x, y) for x, y in zip(r, q, strict=False))

        pos_ref = f"ref [{ref_start_abs if ref_start_abs is not None else '-'}–{ref_end_abs if ref_end_abs is not None else '-'}]"
        pos_read = f"read[{read_start if read_start is not None else '-'}–{read_end if read_end is not None else '-'}]"

        rows_html.append(
            f"<div style='font-family:monospace'>{pos_ref} &nbsp; {pos_read}</div>"
            f"<pre style='margin:2px 0;font-family:monospace'>{ruler_line}</pre>"
            f"<pre style='margin:2px 0;font-family:monospace'>{line_ref}</pre>"
            f"<pre style='margin:2px 0;font-family:monospace'>{line_read}</pre>"
            "<hr style='border:none;border-top:1px dashed #888;margin:6px 0'>"
        )

    # Header + final wrapper
    header = (
        "<div style='font-family:system-ui'>"
        f"<b>{html.escape(str(pref.file))}</b> [{s_eff}–{e_eff}] &nbsp;vs&nbsp; "
        f"<b>{html.escape(str(rd.name))}</b> "
        f"(strand {html.escape(str(strand))}; span={html.escape(str(span))})<br>"
        f"PID={row['pid']:.2f}% &nbsp; core={row['core_len']} &nbsp; "
        f"score={row['score']:.1f} &nbsp; SNPs={row['snps']}</div>"
    )

    html_out = (
        f"<div style='background:{css_bg};color:{css_fg};"
        f"padding:8px;border-radius:8px;font-size:{font_px}px'>"
        f"{header}<hr>" + "".join(rows_html) + "</div>"
    )
    display(HTML(html_out))

    # ── optional static PNG/PDF export (Matplotlib, no external binaries) ─────
    if save_prefix is None:
        return

    import math as _math
    import os as _os

    import matplotlib as _mpl

    _mpl.use("Agg")
    import matplotlib.patches as _mpatch
    import matplotlib.pyplot as _plt

    _MONO = "DejaVu Sans Mono"  # noqa: N806
    _SANS = "Liberation Sans"  # noqa: N806
    _MM = 1 / 25.4  # noqa: N806
    _CPR = save_chars_per_row  # noqa: N806
    _FS = save_font_size  # noqa: N806

    # export color palette (matches theme)
    if theme.lower() == "light":
        _col_bg, _col_fg = "#ffffff", "#222222"
        _col_match = "#a5d6a7"
        _col_mis = "#ef9a9a"
        _col_gap = "#e0e0e0"
        _col_tick = "#888888"
    else:
        _col_bg, _col_fg = "#111111", "#eeeeee"
        _col_match = "#1b5e20"
        _col_mis = "#b71c1c"
        _col_gap = "#424242"
        _col_tick = "#bbbbbb"

    # cumulative non-gap counters for export (independent of HTML counters above)
    _rc: list[int] = [0]
    _qc: list[int] = [0]
    for _cr, _cq in zip(a_ref, a_read, strict=False):
        _rc.append(_rc[-1] + (_cr != "-"))
        _qc.append(_qc[-1] + (_cq != "-"))

    # figure height: 4 text lines per block + 2 header lines
    _n_blocks = _math.ceil(len(a_ref) / _CPR)
    _n_lines = _n_blocks * 4 + 2
    _line_h_mm = max(4.5, _FS * 0.55)
    _fig_h_mm = _n_lines * _line_h_mm + 10

    _fig = _plt.figure(
        figsize=(save_fig_width_mm * _MM, _fig_h_mm * _MM),
        dpi=save_dpi,
    )
    _fig.patch.set_facecolor(_col_bg)
    _ax = _fig.add_axes([0.01, 0.0, 0.99, 1.0])
    _ax.set_facecolor(_col_bg)
    _ax.set_xlim(0, 1)
    _ax.set_ylim(0, 1)
    _ax.axis("off")

    _n_wide = _CPR + 16
    _cw = 1.0 / _n_wide  # width of one character in axes units
    _lh = 1.0 / (_n_lines + 1)
    _cy = 1.0 - _lh * 0.6  # current y position

    def _txt(
        x, y, s, color=None, fs=None, mono=True, bold=False
    ):  # pylint: disable=too-many-arguments
        _ax.text(
            x,
            y,
            s,
            color=color or _col_fg,
            fontsize=fs or _FS,
            fontfamily=_MONO if mono else _SANS,
            va="center",
            ha="left",
            fontweight="bold" if bold else "normal",
            transform=_ax.transAxes,
        )

    def _box(x, y, w, h, color):
        _ax.add_patch(
            _mpatch.FancyBboxPatch(
                (x, y),
                w,
                h,
                boxstyle="square,pad=0",
                linewidth=0,
                facecolor=color,
                zorder=2,
                transform=_ax.transAxes,
            )
        )

    # header line
    _txt(
        0.0,
        _cy,
        f"{str(pref.file)}  vs  {str(rd.name)}"
        f"  |  strand={strand}"
        f"  |  PID={row['pid']:.2f}%"
        f"  |  core={row['core_len']}"
        f"  |  SNPs={row['snps']}",
        fs=_FS * 0.9,
        mono=False,
        bold=True,
    )
    _cy -= _lh * 1.8

    for _i0 in range(0, len(a_ref), _CPR):
        _r = a_ref[_i0 : _i0 + _CPR]
        _q = a_read[_i0 : _i0 + _CPR]

        _rs = next((s_eff + _rc[_i0 + j] for j, c in enumerate(_r) if c != "-"), None)
        _re = next(
            (s_eff + _rc[_i0 + j] for j, c in enumerate(reversed(_r), 1) if c != "-"),
            None,
        )
        _qs = next((_qc[_i0 + j] for j, c in enumerate(_q) if c != "-"), None)
        _qe = next(
            (_qc[_i0 + j] for j, c in enumerate(reversed(_q), 1) if c != "-"), None
        )

        _lr = f"ref [{_rs}-{_re}]"
        _lq = f"read[{_qs}-{_qe}]"
        _mg = max(len(_lr), len(_lq)) + 1

        # ruler
        _cx = _mg * _cw
        _skip = 0
        for _j, _ch in enumerate(_r):
            if _ch == "-":
                _cx += _cw
                continue
            if _skip:
                _skip -= 1
                _cx += _cw
                continue
            _p = s_eff + _rc[_i0 + _j]
            if _p % 10 == 0:
                _lab = str(_p)
                _txt(_cx, _cy, _lab, color=_col_tick, fs=_FS * 0.8)
                _cx += len(_lab) * _cw
                _skip = len(_lab) - 1
            elif _p % 5 == 0:
                _txt(_cx, _cy, "|", color=_col_tick, fs=_FS * 0.8)
                _cx += _cw
            else:
                _txt(_cx, _cy, "·", color=_col_tick, fs=_FS * 0.65)
                _cx += _cw
        _cy -= _lh

        # reference row
        _txt(0.0, _cy, _lr, color=_col_tick, fs=_FS * 0.85)
        _cx = _mg * _cw
        for _xc, _yc in zip(_r, _q, strict=False):
            _bc = (
                _col_gap
                if (_xc == "-" or _yc == "-")
                else (_col_match if _xc == _yc else _col_mis)
            )
            _box(_cx - _cw * 0.05, _cy - _lh * 0.44, _cw * 1.05, _lh * 0.88, _bc)
            _txt(_cx, _cy, _xc)
            _cx += _cw
        _cy -= _lh

        # read row
        _txt(0.0, _cy, _lq, color=_col_tick, fs=_FS * 0.85)
        _cx = _mg * _cw
        for _xc, _yc in zip(_r, _q, strict=False):
            _bc = (
                _col_gap
                if (_xc == "-" or _yc == "-")
                else (_col_match if _xc == _yc else _col_mis)
            )
            _box(_cx - _cw * 0.05, _cy - _lh * 0.44, _cw * 1.05, _lh * 0.88, _bc)
            _txt(_cx, _cy, _yc)
            _cx += _cw
        _cy -= _lh

        # block separator
        _ax.axhline(
            _cy + _lh * 0.35,
            color="#555555",
            lw=0.5,
            linestyle="--",
            xmin=0.0,
            xmax=1.0,
        )
        _cy -= _lh * 0.6

    # save PNG and PDF
    _os.makedirs(_os.path.dirname(_os.path.abspath(save_prefix)), exist_ok=True)
    for _ext in ("png", "pdf"):
        _fig.savefig(
            f"{save_prefix}.{_ext}",
            dpi=save_dpi,
            bbox_inches="tight",
            facecolor=_col_bg,
        )
    _plt.close(_fig)
    print(f"Saved: {save_prefix}.png  |  {save_prefix}.pdf")


def view_alignment_html_by_ref(
    ref: str | int,
    df,  # pandas.DataFrame (locker getypt)
    plasmids: Sequence[Any],
    reads: Sequence[Any],
    *,
    nth: int = 0,
    min_pid: float = 0.0,
    **kwargs: Any,
) -> Optional[Any]:
    """
    Pick the *nth* best row for a numeric reference (by PID then score)
    and render the HTML alignment via :func:`view_alignment_html`.

    Parameters
    ----------
    ref : str or int
        Numeric token in ``sequence_name`` (e.g. ``2034`` for ``EF...2034_*``).
    df : pandas.DataFrame
        Results table containing the usual columns.
    nth : int, default 0
        Rank among matches (0 = best).
    min_pid : float, default 0.0
        Drop rows with ``pid < min_pid`` before ranking.
    **kwargs
        Forwarded to :func:`view_alignment_html` (e.g. ``wrap=90``,
        ``theme="light"``, ``span="full"``).

    Returns
    -------
    pandas.Series or None
        The displayed row or ``None`` if nothing matched.
    """
    key = df["sequence_name"].astype(str).map(norm_ref)
    hits = df[key.eq(str(ref))]
    if min_pid > 0:
        hits = hits[hits["pid"] >= float(min_pid)]
    hits = hits.sort_values(["pid", "score"], ascending=[False, False])

    if hits.empty:
        display(HTML(f"<b>No rows for ref #{html.escape(str(ref))}.</b>"))
        return None

    row = hits.iloc[min(nth, len(hits) - 1)]
    view_alignment_html(row, plasmids, reads, **kwargs)
    return row
