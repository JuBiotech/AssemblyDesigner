"""QC summary visualization: binary pass/fail sequence-to-plasmid mapping dotplot.

Builds a single overview figure from a results table (as produced by
:func:`results_to_dataframe` / :func:`export_results` / :func:`dealign_reads_easy`)
showing, for every read × construct pair, whether the alignment passed a PID
threshold. Accepts either an in-memory DataFrame or a path to a saved
Excel file.
"""

from __future__ import annotations

import re
from collections.abc import Callable
from pathlib import Path
from typing import Optional, Union

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.figure import Figure

_FONT = "Liberation Sans"
_MM = 1 / 25.4

# Known RBS/5'UTR tokens present in the construct naming scheme.
# Adjust this list if your construct library uses other RBS names.
RBS_TOKENS = ["BCD12", "BCD2", "BCD8", "B0032m", "B0033m", "B0034m"]


def _prep_pid(df: pd.DataFrame, pid_col: str = "pid") -> pd.DataFrame:
    d = df.copy()
    if pid_col in d.columns:
        d[pid_col] = pd.to_numeric(
            d[pid_col].astype(str).str.replace(",", ".", regex=False),
            errors="coerce",
        )
        # auto-detect 0..1 scale -> convert to percent
        if d[pid_col].notna().any() and d[pid_col].max() <= 1.5:
            d[pid_col] = d[pid_col] * 100.0
    return d


def extract_rbs(construct: str) -> str:
    """
    Return the RBS/5'UTR token found in a construct name, by matching
    against the known ``RBS_TOKENS`` list (order-independent, robust to
    missing promoter prefix, e.g. 'B0032m_PanD_B0015').
    """
    parts = str(construct).split("_")
    for tok in parts:
        if tok in RBS_TOKENS:
            return tok
    return "unknown"


def default_x_label(construct: str) -> str:
    """Remove the constant suffix '_PanD_B0015' (or '_ecPanD_B0015')."""
    return re.sub(r"_(cg|ec)?PanD_B0015$", "", str(construct))


def default_y_label(seq_name: str) -> str:
    """Extract the longest digit block from a FASTA ID."""
    m = re.search(r"(\d{5,})", str(seq_name))
    return m.group(1) if m else str(seq_name)[-10:]


def _darken(hex_color: str, factor: float = 0.7) -> str:
    h = hex_color.lstrip("#")
    r, g, b = (int(h[i : i + 2], 16) for i in (0, 2, 4))
    return f"#{int(r * factor):02x}{int(g * factor):02x}{int(b * factor):02x}"


def plot_mapping_dotplot_binary(
    df: pd.DataFrame,
    *,
    x_col: str = "construct",
    y_col: str = "sequence_name",
    pid_col: str = "pid",
    pid_pass: float = 99.999,
    x_label_fn: Callable[[str], str] = default_x_label,
    y_label_fn: Callable[[str], str] = default_y_label,
    group_col: Optional[str] = "rbs",
    color_pass: str = "#2ca02c",
    color_fail: str = "#d62728",
    dot_size: float = 32,
    fig_width_mm: float = 180,
    fig_height_mm: Optional[float] = None,
    font_size_tick: float = 6.5,
    font_size_label: float = 9.0,
    font_size_title: float = 10.0,
    x_rotation: int = 50,
    dpi: int = 300,
    theme: str = "white",
    save_prefix: Optional[str] = None,
    title: str = "Sequence-to-plasmid mapping",
) -> Figure:
    """Binary pass/fail sequence-to-plasmid mapping dotplot.

    Parameters
    ----------
    df : pandas.DataFrame
        Results table (one row per read × best-hit construct), typically from
        :func:`results_to_dataframe`.
    x_col, y_col, pid_col : str
        Column names for the construct, read identifier, and percent identity.
    pid_pass : float
        Minimum PID (in the same scale as ``pid_col``, percent by default) to
        count as a pass.
    x_label_fn, y_label_fn : callable
        Functions to derive shortened axis tick labels from the raw values.
    group_col : str or None
        Column used to group/label constructs along the x-axis (e.g. by RBS).
        If missing from ``df``, it is derived via :func:`extract_rbs`.
    save_prefix : str or None
        If given, save ``{save_prefix}.png/.svg/.pdf`` alongside returning the figure.

    Returns
    -------
    matplotlib.figure.Figure
    """
    bg, fg = ("#111111", "#eeeeee") if theme == "dark" else ("#ffffff", "#222222")

    d = _prep_pid(df, pid_col)

    if group_col is not None and group_col not in d.columns:
        d[group_col] = d[x_col].apply(extract_rbs)

    if group_col is not None:
        d = d.sort_values([group_col, x_col])

    x_cats = d[x_col].astype(str).unique().tolist()
    y_cats = sorted(d[y_col].astype(str).unique(), reverse=True)

    x_map = {v: i for i, v in enumerate(x_cats)}
    y_map = {v: i for i, v in enumerate(y_cats)}
    n_x, n_y = len(x_cats), len(y_cats)

    x_labels = [x_label_fn(v) for v in x_cats]
    y_labels = [y_label_fn(v) for v in y_cats]

    xs = d[x_col].astype(str).map(x_map).values
    ys = d[y_col].astype(str).map(y_map).values

    d["pass"] = d[pid_col] >= pid_pass
    m_pass = d["pass"].values
    m_fail = ~m_pass

    if fig_height_mm is None:
        fig_height_mm = min(240.0, max(80.0, n_y * 3.8 + 40))

    fig, ax = plt.subplots(figsize=(fig_width_mm * _MM, fig_height_mm * _MM), dpi=dpi)
    fig.patch.set_facecolor(bg)
    ax.set_facecolor(bg)

    # ── scatter ────────────────────────────────────────────────────────────
    if m_pass.any():
        ax.scatter(
            xs[m_pass],
            ys[m_pass],
            s=dot_size,
            color=color_pass,
            linewidths=0.4,
            edgecolors=_darken(color_pass),
            zorder=4,
            alpha=0.92,
        )

    if m_fail.any():
        ax.scatter(
            xs[m_fail],
            ys[m_fail],
            s=dot_size,
            color=color_fail,
            marker="X",
            linewidths=0.4,
            edgecolors=_darken(color_fail),
            zorder=4,
            alpha=0.92,
        )

    # ── legend ─────────────────────────────────────────────────────────────
    handles = []
    if m_pass.any():
        handles.append(
            mpatches.Patch(
                color=color_pass,
                label=f"PID = 100%  (verified, n={int(m_pass.sum())})",
            )
        )
    if m_fail.any():
        handles.append(
            mpatches.Patch(
                color=color_fail,
                label=f"PID < 100%  (suspicious, n={int(m_fail.sum())})",
            )
        )

    leg = ax.legend(
        handles=handles,
        frameon=False,
        fontsize=font_size_tick + 0.5,
        loc="upper right",
        labelcolor=fg,
    )
    for t in leg.get_texts():
        t.set_fontfamily(_FONT)

    # ── group separators + top-axis labels ──────────────────────────────────
    if group_col is not None:
        group_seq = d.drop_duplicates(subset=[x_col])[group_col].tolist()
        boundaries = []
        for i in range(1, len(group_seq)):
            if group_seq[i] != group_seq[i - 1]:
                boundaries.append(i - 0.5)
        for b in boundaries:
            ax.axvline(b, color="#bbbbbb", lw=0.7, linestyle="--", zorder=2)

        starts = [0] + [int(b + 0.5) for b in boundaries]
        ends = [int(b + 0.5) for b in boundaries] + [n_x]
        for s, e in zip(starts, ends, strict=False):
            label = group_seq[s]
            ax.text(
                (s + e - 1) / 2,
                n_y - 0.3 + n_y * 0.02,
                label,
                ha="center",
                va="bottom",
                fontsize=font_size_label,
                fontweight="bold",
                color=fg,
                fontfamily=_FONT,
            )

    # ── grid ───────────────────────────────────────────────────────────────
    ax.set_axisbelow(True)
    ax.grid(
        True,
        color="#e0e0e0" if theme == "white" else "#333333",
        linewidth=0.35,
        linestyle=":",
        zorder=1,
    )

    # ── axes ───────────────────────────────────────────────────────────────
    ax.set_xlim(-0.5, n_x - 0.5)
    ax.set_ylim(-0.5, n_y - 0.5)
    ax.set_xticks(range(n_x))
    ax.set_xticklabels(
        x_labels,
        rotation=x_rotation,
        ha="right",
        fontsize=font_size_tick,
        fontfamily=_FONT,
        color=fg,
    )
    ax.set_yticks(range(n_y))
    ax.set_yticklabels(y_labels, fontsize=font_size_tick, fontfamily=_FONT, color=fg)
    ax.set_xlabel("Construct", fontsize=font_size_label, fontfamily=_FONT, color=fg)
    ax.set_ylabel(
        "Sequencing read (FASTA ID)",
        fontsize=font_size_label,
        fontfamily=_FONT,
        color=fg,
    )
    ax.set_title(title, fontsize=font_size_title, fontfamily=_FONT, color=fg)

    for sp in ax.spines.values():
        sp.set_color(fg)
        sp.set_linewidth(0.5)
    ax.tick_params(colors=fg)

    info = (
        f"n reads = {len(d)}   verified = {int(m_pass.sum())}   "
        f"suspicious = {int(m_fail.sum())}"
    )
    ax.text(
        0.01,
        0.01,
        info,
        transform=ax.transAxes,
        fontsize=font_size_tick,
        fontfamily=_FONT,
        color="#888888",
        va="bottom",
    )

    fig.tight_layout()

    if save_prefix:
        for ext in ("png", "svg", "pdf"):
            fig.savefig(
                f"{save_prefix}.{ext}", dpi=dpi, bbox_inches="tight", facecolor=bg
            )
        print(f"Saved: {save_prefix}.png / .svg / .pdf")

    return fig


def _load_results(
    data: Union[str, Path, pd.DataFrame],
    *,
    sheet_name: int | str,
) -> pd.DataFrame:
    """Resolve `data` into a DataFrame, whether it's a path or already a df."""
    if isinstance(data, pd.DataFrame):
        return data
    return pd.read_excel(data, sheet_name=sheet_name)


def plot_alignment_results_summary(
    data: Union[str, Path, pd.DataFrame],
    *,
    sheet_name: int | str = 0,
    max_reads: Optional[int] = 42,
    save_prefix: Optional[str] = None,
    **plot_kwargs,
) -> Figure:
    """
    Build the mapping dotplot from alignment results.

    Parameters
    ----------
    data : str, Path, or pandas.DataFrame
        Either a path to an exported results Excel file, or a DataFrame
        already in memory (e.g. the ``df`` returned by
        ``results_to_dataframe(rows)`` / ``dealign_reads_easy(...)``) —
        no need to export to Excel first.
    sheet_name : int or str
        Sheet to read when `data` is a file path (default: first sheet).
        Ignored when `data` is already a DataFrame.
    max_reads : int or None
        Only use the first ``max_reads`` rows (in file/frame order) so the
        figure does not grow too tall. Default 42. Set to None to use
        all reads.
    save_prefix : str or None
        Base path for output files (.png/.svg/.pdf appended). ``None`` by
        default (figure is returned but not written to disk).
    **plot_kwargs
        Passed through to :func:`plot_mapping_dotplot_binary`.

    Returns
    -------
    matplotlib.figure.Figure
    """
    df_best = _load_results(data, sheet_name=sheet_name)
    if max_reads is not None:
        df_best = df_best.head(max_reads)
    return plot_mapping_dotplot_binary(df_best, save_prefix=save_prefix, **plot_kwargs)
