from __future__ import annotations

import io
import logging
import os
import warnings
from collections.abc import Callable, Iterable, Sequence
from pathlib import Path
from typing import (
    Any,
    Optional,
    Union,
)

import fastprogress
import matplotlib
import mpl_toolkits.axes_grid1
import numpy as np
import pandas
from matplotlib import cm, colors, pyplot
from PIL import Image as PILImage

Data = Union[pandas.Series, dict[str, int | float]]
LOGGER = logging.getLogger(__name__)


def setup_logging(
    level: int | str = "INFO",
    name: str = "assembly_designer",
    propagate: bool = False,
) -> logging.Logger:
    """Configure and return the package logger.

    This function sets the log level, installs a stream handler if none is
    present, and applies a concise formatter. It is safe to call multiple
    times (no duplicate handlers will be added).

    Parameters
    ----------
    level : int or str, optional
        Log level (e.g., ``"DEBUG"``, ``"INFO"``, ``"WARNING"`` or
        ``logging.DEBUG``). Default is ``"INFO"``.
    name : str, optional
        Logger name to configure. Default is ``"assembly_designer"``.
    propagate : bool, optional
        Whether the logger should propagate records to the root logger.
        Default is ``False``.

    Returns
    -------
    logging.Logger
        The configured logger instance.

    Examples
    --------
    >>> logger = setup_logging("DEBUG")
    >>> logger.info("CLI started")
    """
    # Resolve string levels like "INFO" to numeric
    if isinstance(level, str):
        numeric = logging.getLevelName(level.upper())
        if not isinstance(numeric, int):
            raise ValueError(f"Unknown log level: {level!r}")
        level = numeric

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.propagate = propagate

    # Add a single StreamHandler if none present
    if not any(isinstance(h, logging.StreamHandler) for h in logger.handlers):
        handler = logging.StreamHandler()
        handler.setLevel(level)
        handler.setFormatter(
            logging.Formatter(
                fmt="%(asctime)s %(levelname)s %(name)s: %(message)s",
                datefmt="%H:%M:%S",
            )
        )
        logger.addHandler(handler)

    return logger


def setup_run_logging_files(
    log_dir: str,
    *,
    logfile_name: str = "assembly.log",
    warnfile_name: str = "assembly_warnerr.log",
    logger_name: str = "assembly_designer",  # passt zu deinem setup_logging()
    console_level: int | str = "INFO",
    file_level: int | str = "DEBUG",
    warn_file_level: int | str = "WARNING",
) -> dict[str, logging.Handler]:
    """Attach file handlers for the current run (full + warn/error) to the package logger.

    Creates/attaches two FileHandlers (deduplicated):
      - full log (DEBUG+):  assembly.log
      - warn/error only:    assembly_warnerr.log
    Keeps your existing console handler from setup_logging() untouched.
    """

    # -- helpers ---------------------------------------------------------------
    def _norm(level: int | str) -> int:
        """Normalize int/str log-level to an int (mypy-friendly)."""
        if isinstance(level, str):
            val = logging.getLevelName(level.upper())
        else:
            val = level
        if isinstance(val, str):
            raise ValueError(f"Unknown log level: {level!r}")
        return int(val)

    # resolve levels if strings
    console_level = _norm(console_level)
    file_level = _norm(file_level)
    warn_file_level = _norm(warn_file_level)

    os.makedirs(log_dir or ".", exist_ok=True)
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)  # low floor; handlers filter

    # remember existing file paths to avoid duplicates (common in notebooks)
    existing_paths: set[str | None] = {
        getattr(h, "baseFilename", None) for h in logger.handlers
    }

    fmt = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(name)s:%(lineno)d | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # console (only ensure level; setup_logging already added a StreamHandler)
    for h in logger.handlers:
        if isinstance(h, logging.StreamHandler) and not hasattr(h, "baseFilename"):
            h.setLevel(console_level)

    # full log file
    full_path = os.path.join(log_dir, logfile_name)
    if full_path not in existing_paths:
        fh: logging.Handler = logging.FileHandler(full_path, mode="w", encoding="utf-8")
        fh.setLevel(file_level)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    else:
        matches = [
            h for h in logger.handlers if getattr(h, "baseFilename", "") == full_path
        ]
        if matches:
            fh = matches[0]
        else:
            # Fallback (should practically never happen)
            fh = logging.FileHandler(full_path, mode="w", encoding="utf-8")
            fh.setLevel(file_level)
            fh.setFormatter(fmt)
            logger.addHandler(fh)

    # warn/error file
    warn_path = os.path.join(log_dir, warnfile_name)
    if warn_path not in existing_paths:
        wh: logging.Handler = logging.FileHandler(warn_path, mode="w", encoding="utf-8")
        wh.setLevel(warn_file_level)
        wh.setFormatter(fmt)
        logger.addHandler(wh)
    else:
        matches = [
            h for h in logger.handlers if getattr(h, "baseFilename", "") == warn_path
        ]
        if matches:
            wh = matches[0]
        else:
            # Fallback (should practically never happen)
            wh = logging.FileHandler(warn_path, mode="w", encoding="utf-8")
            wh.setLevel(warn_file_level)
            wh.setFormatter(fmt)
            logger.addHandler(wh)

    logger.debug("File logging ready | full=%s | warn=%s", full_path, warn_path)
    return {"file_full": fh, "file_warn": wh}


def check_and_log_volumes(
    df: pandas.DataFrame,
    volume_column: str,
    *,
    min_volume_ul: float = 1.0,
    group_columns: Optional[Sequence[str]] = None,
    context: str = "",
    logger: Optional[logging.Logger] = None,
) -> None:
    """Validate volumes and log warnings/errors (NumPy-style docstring).

    Logs:
      - ERROR: missing/NaN volumes
      - WARNING: volume < `min_volume_ul`
      - INFO: summary (totals)

    Parameters
    ----------
    df : pandas.DataFrame
        Input table containing the `volume_column`.
    volume_column : str
        Name of the volume column (µL).
    min_volume_ul : float, optional
        Threshold (µL) below which a WARNING is logged. Default is 1.0.
    group_columns : sequence of str or None, optional
        Columns printed to identify rows (e.g., ["Mastermix", "Well Nr. destination plate"]).
        If None, a sensible subset is auto-selected when available.
    context : str, optional
        Free-form label included in messages (e.g., a step name).
    logger : logging.Logger or None, optional
        Logger to use; falls back to the package LOGGER.
    """
    log = logger or LOGGER

    if volume_column not in df.columns:
        log.error(
            "Volume check failed: column %r not found. %s", volume_column, context
        )
        return

    # Auto-pick identifiers if not provided
    if group_columns is None:
        candidates = [
            "Mastermix",
            "Plasmidmix",
            "Part",
            "Well",
            "Well Nr. destination plate",
            "Plate",
        ]
        group_columns = [c for c in candidates if c in df.columns]

    # Compute masks
    is_na = df[volume_column].isna()
    vol_numeric = pandas.to_numeric(df[volume_column], errors="coerce")
    below = vol_numeric < float(min_volume_ul)

    # Log missing volumes (ERROR)
    if bool(is_na.any()):
        cols = (list(group_columns) if group_columns else []) + [volume_column]
        for _, row in df.loc[is_na, cols].iterrows():
            ident = (
                ", ".join(f"{c}={row[c]}" for c in group_columns)
                if group_columns
                else "(row)"
            )
            log.error(
                "Missing volume: %s | %s | %s=%r",
                context,
                ident,
                volume_column,
                row[volume_column],
            )

    # Log low volumes (WARNING)
    if bool(below.fillna(False).any()):
        cols = (list(group_columns) if group_columns else []) + [volume_column]
        for _, row in df.loc[below.fillna(False), cols].iterrows():
            ident = (
                ", ".join(f"{c}={row[c]}" for c in group_columns)
                if group_columns
                else "(row)"
            )
            # row[volume_column] is possible here str/obj; therefore vol_numeric for output
            v = float(pandas.to_numeric(row[volume_column], errors="coerce"))
            log.warning(
                "Low volume: %s | %s | %s=%.3f µL < min=%.3f µL",
                context,
                ident,
                volume_column,
                v,
                float(min_volume_ul),
            )

    # Summary (INFO)
    n_na = int(is_na.sum())
    n_low = int(below.fillna(False).sum())
    total = int(len(df))
    log.info(
        "Volume check summary: total=%d, missing=%d, low(<%.3f µL)=%d%s",
        total,
        n_na,
        float(min_volume_ul),
        n_low,
        f" | {context}" if context else "",
    )


# ---------------------------------------------------------------------------


def get_crop_slices(well_ids: Iterable[str]) -> tuple[slice, slice]:
    """Compute minimal row/column slices that contain all given well IDs.

    Parameters
    ----------
    well_ids : iterable of str
        Alphanumeric well IDs (e.g., ``"A01"``, ``"H12"``). Case-insensitive.
        Must not be empty.

    Returns
    -------
    rslice : slice
        Slice covering 0-based rows from the topmost to the bottommost well.
    cslice : slice
        Slice covering 0-based columns from the leftmost to the rightmost well.

    Raises
    ------
    ValueError
        If ``well_ids`` is empty or contains malformed IDs.
    """
    well_ids = list(well_ids)
    if not well_ids:
        raise ValueError("well_ids must not be empty.")

    rows: list[int] = []
    cols: list[int] = []

    for w in well_ids:
        if not isinstance(w, str) or len(w) < 2:
            raise ValueError(f"Malformed well ID: {w!r}")
        w = w.strip().upper()

        row_char = w[0]
        col_part = w[1:]

        if not ("A" <= row_char <= "Z") or not col_part.isdigit():
            raise ValueError(f"Malformed well ID: {w!r}")

        row_idx = ord(row_char) - ord("A")
        col_idx = int(col_part) - 1  # to 0-based

        if row_idx < 0 or col_idx < 0:
            raise ValueError(f"Malformed well ID: {w!r}")

        rows.append(row_idx)
        cols.append(col_idx)

    rslice = slice(min(rows), max(rows) + 1)
    cslice = slice(min(cols), max(cols) + 1)
    return rslice, cslice


def make_well_array(rows: int, cols: int) -> np.ndarray:
    """Create a 2D array of well IDs like ``A01``, ``B12`` for a plate.

    Parameters
    ----------
    rows : int
        Number of plate rows (max. 26, i.e., ``A``–``Z``).
    cols : int
        Number of plate columns (positive integer).

    Returns
    -------
    np.ndarray
        Array of shape ``(rows, cols)`` with alphanumeric well IDs.

    Raises
    ------
    ValueError
        If ``rows`` or ``cols`` are not positive, or if ``rows`` > 26.
    """
    if rows <= 0 or cols <= 0:
        raise ValueError("rows and cols must be positive integers.")
    if rows > 26:
        raise ValueError("rows must be ≤ 26 (A–Z).")

    return np.array(
        [
            [f"{row_letter}{col:02d}" for col in range(1, cols + 1)]
            for row_letter in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[:rows]
        ]
    )


def values_in_2d(
    data: pandas.Series | dict[str, int | float],
    shape: tuple[int, int] | None = None,
    *,
    missing_values: float = np.nan,
) -> np.ndarray:
    """Convert well-indexed data into a 2D array aligned to a plate.

    Parameters
    ----------
    data : pandas.Series or dict
        Mapping from alphanumeric well IDs (e.g., ``"A01"``) to values.
        If a dict is provided it is converted to a ``Series``.
    shape : tuple of (int, int), optional
        Full plate shape ``(rows, cols)``. If ``None``, the minimal shape
        covering all wells in ``data`` is inferred automatically.
    missing_values : float, optional
        Fill value for wells missing from ``data`` (default: ``np.nan``).

    Returns
    -------
    np.ndarray
        A 2D array with values placed at positions specified by well IDs.

    Raises
    ------
    ValueError
        If ``data`` contains malformed well IDs or ``shape`` is invalid.
    """
    series = pandas.Series(data=data) if isinstance(data, dict) else data
    if not isinstance(series, pandas.Series):
        raise ValueError("data must be a pandas.Series or a dict[str, number].")

    # Determine cropping slices and full shape
    if shape is None:
        rslice, cslice = get_crop_slices(series.index.astype(str))
        full_shape = (rslice.stop, cslice.stop)
    else:
        if len(shape) != 2 or shape[0] <= 0 or shape[1] <= 0:
            raise ValueError("shape must be a tuple of two positive integers.")
        full_shape = (shape[0], shape[1])
        rslice = slice(0, shape[0])
        cslice = slice(0, shape[1])

    wells = make_well_array(*full_shape)
    flat_vals = [
        series.get(well_id, missing_values)  # Series.get handles missing keys
        for well_id in wells.ravel()
    ]
    values = np.array(flat_vals).reshape(full_shape)

    # Return just the relevant crop
    return values[rslice, cslice]


def annotate_heatmap(
    im: matplotlib.image.AxesImage,
    *,
    valfmt: str | matplotlib.ticker.Formatter = "{x:.2f}",
    textcolors: tuple[str, str] = ("black", "white"),
    threshold: float | None = None,
    **text_kwargs: Any,
) -> list[matplotlib.text.Text]:
    """Add value labels on top of a heatmap.

    Parameters
    ----------
    im : matplotlib.image.AxesImage
        The AxesImage to annotate.
    valfmt : str or matplotlib.ticker.Formatter, optional
        Format string or formatter for values (default ``"{x:.2f}"``).
    textcolors : tuple of str, optional
        Colors for values below/above the threshold (default: (``"black"``, ``"white"``)).
    threshold : float, optional
        Data value used to switch between the two `textcolors`.
        If ``None``, use the midpoint of the colormap.
    **text_kwargs
        Extra keyword arguments forwarded to ``Axes.text``.

    Returns
    -------
    list of matplotlib.text.Text
        The created text artists.
    """
    data = im.get_array()
    thr = im.norm(threshold) if threshold is not None else im.norm(np.nanmax(data) / 2)
    formatter = (
        matplotlib.ticker.StrMethodFormatter(valfmt)
        if isinstance(valfmt, str)
        else valfmt
    )

    base_kwargs: dict[str, Any] = {
        "horizontalalignment": "center",
        "verticalalignment": "center",
        "fontsize": 8,
    }
    base_kwargs.update(text_kwargs)

    texts: list[matplotlib.text.Text] = []
    nrows, ncols = data.shape
    for i in range(nrows):
        for j in range(ncols):
            val = data[i, j]
            if isinstance(val, np.ma.core.MaskedConstant):
                s = "--"
                color = "black"
            else:
                s = formatter(val)
                color = textcolors[int(im.norm(val) > thr)]
            txt = im.axes.text(j, i, s=s, color=color, **base_kwargs)
            texts.append(txt)
    return texts


def to_colormap(col: tuple[int, int, int]) -> colors.ListedColormap:
    """Create a transparent→color colormap from an RGB triple.

    Parameters
    ----------
    col : tuple of int
        RGB values in the range 0–255 (e.g., ``(191, 21, 33)``).

    Returns
    -------
    matplotlib.colors.ListedColormap
        Colormap that fades from fully transparent white to the given color.
    """
    n = 256
    rgba = np.array((*((np.array(col) / 255.0)[:3]), 1.0))
    white = np.ones(4)
    ramp = np.array([(1 - t) * white + t * rgba for t in np.linspace(0, 1, n)])
    ramp[:, 3] = np.linspace(0, 1, n)  # alpha
    return colors.ListedColormap(ramp)


def mtpshow(
    ax: matplotlib.axes.Axes,
    data: Data | np.ndarray,
    *,
    shape: tuple[int, int] | None = None,
    cbar_kwargs: dict[str, Any] | None = None,
    cbar_label: str | None = None,
    grid: bool = True,
    annotate: bool = True,
    **imshow_kwargs: Any,
) -> matplotlib.image.AxesImage:
    """Display microtiter-plate data similar to ``ax.imshow``.

    Works best with a default figure size of roughly (6, 4) per Axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    data : np.ndarray or mapping-like
        Either a 2D array of values or a mapping/Series from well IDs
        (e.g., ``"A01"``) to values.
    shape : tuple of (int, int), optional
        Full plate shape. If ``None``, infer minimal shape from well IDs.
    cbar_kwargs : dict, optional
        Extra kwargs forwarded to ``Figure.colorbar``.
    cbar_label : str, optional
        Label for the colorbar. If ``None``, no label is set.
    grid : bool, optional
        Draw thin grid lines between wells (default ``True``).
    annotate : bool, optional
        If ``True``, overlay numeric values using :func:`annotate_heatmap`.
    **imshow_kwargs
        Forwarded to ``ax.imshow`` (e.g., ``vmin``, ``vmax``, ``cmap``).

    Returns
    -------
    matplotlib.image.AxesImage
        The image artist returned by ``ax.imshow``.

    Raises
    ------
    ValueError
        If ``data`` is neither a 2D array nor a well-indexed mapping.
    """
    if isinstance(data, np.ndarray):
        values = data
    elif isinstance(data, (dict, pandas.Series)):
        values = values_in_2d(data, shape)
    else:
        raise ValueError(
            "data must be a 2D np.ndarray or a pandas.Series/dict keyed by well IDs."
        )

    # Determine tick crops for labels
    if shape is None and isinstance(data, pandas.Series):
        rslice, cslice = get_crop_slices(data.index.astype(str))
    elif shape is None and isinstance(data, dict):
        rslice, cslice = get_crop_slices([str(k) for k in data.keys()])
    else:
        rslice = slice(0, values.shape[0])
        cslice = slice(0, values.shape[1])

    imshow_kw: dict[str, Any] = {"cmap": "binary"}
    imshow_kw.update(imshow_kwargs)
    im = ax.imshow(values, **imshow_kw)

    # Colorbar alongside the image
    divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    bar = ax.figure.colorbar(im, cax=cax, **(cbar_kwargs or {}))
    if cbar_label is not None:
        bar.ax.set_ylabel(cbar_label)

    # Axis ticks & labels
    ax.xaxis.tick_top()
    ax.set(
        yticks=np.arange(values.shape[0]),
        yticklabels=[chr(ord("A") + r) for r in range(rslice.start, rslice.stop)],
        xticks=np.arange(values.shape[1]),
        xticklabels=[str(c + 1) for c in range(cslice.start, cslice.stop)],
    )

    # Optional grid
    if grid:
        ax.set_xticks(np.arange(0, values.shape[1]) - 0.5, minor=True)
        ax.set_yticks(np.arange(0, values.shape[0]) - 0.5, minor=True)
        ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.2)
        ax.tick_params(which="minor", top=False, left=False)

    if annotate:
        annotate_heatmap(im)

    return im


def fn_plot(
    concentrations: np.ndarray | pandas.Series | dict[str, int | float], title: str
) -> None:
    """Plot a microtiter-plate heatmap for concentration data in opaque blue.

    Parameters
    ----------
    concentrations : numpy.ndarray or pandas.Series or dict[str, int | float]
        2D array or well-indexed data (alphanumeric well IDs).
    title : str
        Title for the plot.

    Notes
    -----
    - Uses an **opaque** blue colormap (``cm.Blues``) to avoid grey/black
      appearance when values are near zero or when transparency is undesired.
    - Sets ``vmin=0`` so the colorbar starts at zero.
    """
    fig, ax = pyplot.subplots()
    mtpshow(ax, concentrations, cmap=cm.Blues, vmin=0)
    ax.set(ylabel="row", xlabel="column", title=title)


def fn_plot2(label_vols):
    """
    Wrapper function to plot concentration data with a label as the title.

    Parameters
    ----------
    label_vols : tuple
        A tuple containing:
        - label (str): The title of the plot.
        - vols (array-like): A 2D array or matrix representing the concentration data for the plot.

    Returns
    -------
    None
        Displays the heatmap plot of the concentration data using `fn_plot`.

    Notes
    -----
    - This function extracts the label and volumes from the input tuple and calls `fn_plot`.
    - The `label` is used as the title of the plot.
    """
    label, vols = label_vols
    fn_plot(vols, label)


def plot_gif(
    fn_plot: Callable[[Any], None],
    fp_out: str | os.PathLike[str],
    *,
    data: Iterable[Any],
    fps: float | int = 3,
    delay_frames: int = 3,
    close: bool = True,
) -> Path:
    """Create an animated GIF from matplotlib figures.

    Parameters
    ----------
    fn_plot : callable
        Function that takes one element from ``data`` and draws/updates the current figure.
        If the same figure is reused between frames, set ``close=False``.
    fp_out : str or os.PathLike
        Output filename or path (e.g. ``"workflow.gif"``).
    data : iterable
        Elements to iterate over and pass to ``fn_plot``.
    fps : float or int, optional
        Frames per second of the GIF (default: 3).
    delay_frames : int, optional
        Number of extra copies of the last frame to append (default: 3).
    close : bool, optional
        Close the current figure after each frame (default: True).

    Returns
    -------
    pathlib.Path
        Path to the written GIF.

    Raises
    ------
    ValueError
        If ``fps`` <= 0 or ``data`` yields no frames.
    """
    out_path = Path(fp_out)
    if fps <= 0:
        raise ValueError("fps must be > 0")

    frames: list[PILImage.Image] = []
    for dat in fastprogress.progress_bar(data):
        fn_plot(dat)
        fig = pyplot.gcf()
        with io.BytesIO() as buf:
            # Save raw RGBA buffer (avoids transparency artifacts in GIF quantization)
            pyplot.savefig(buf, format="raw", facecolor="w")
            buf.seek(0)
            arr = np.frombuffer(buf.getvalue(), dtype="uint8").reshape(
                *fig.canvas.get_width_height()[::-1], 4
            )
            # Remove alpha (white background), then quantize to paletted image for GIF
            frames.append(PILImage.fromarray(arr[..., :3], "RGB").quantize())
        if close:
            pyplot.close()

    if not frames:
        raise ValueError("No frames were produced from the provided data.")

    # Add delay on the last frame
    if delay_frames > 0:
        frames.extend([frames[-1]] * delay_frames)

    out_path.parent.mkdir(parents=True, exist_ok=True)

    duration_ms = int(round(1000.0 / float(fps)))
    frames[0].save(
        out_path,
        save_all=True,
        append_images=frames[1:],
        duration=duration_ms,
        loop=0,
    )
    return out_path


def find_planning_file(
    directory: str | os.PathLike[str],
    type_of_assembly: str,
    *,
    case_sensitive: bool = False,
) -> Path:
    """Find exactly one planning Excel file matching a substring.

    Searches ``directory`` for a **single** ``.xlsx`` file whose
    filename contains the text from ``type_of_assembly`` (e.g., ``"PCR"`` or
    ``"Golden Gate"``). Case sensitivity can be controlled with ``case_sensitive``.
    No interactive input is used.

    Parameters
    ----------
    directory : str or os.PathLike
        Folder in which to search.
    type_of_assembly : str
        Substring that must be present in the filename (e.g., ``"PCR"``).
    case_sensitive : bool, optional
        Whether the comparison should be case-sensitive. Default: ``False``.

    Returns
    -------
    pathlib.Path
        Path to the found file.

    Raises
    ------
    FileNotFoundError
        If the folder does not exist or no matching file was found.
    ValueError
        If **more than one** matching file was found.
    """
    base = Path(directory).resolve()
    if not base.exists():
        raise FileNotFoundError(f"Directory does not exist: {base}")

    if not isinstance(type_of_assembly, str) or not type_of_assembly.strip():
        raise ValueError("type_of_assembly must be a non-empty string.")

    # collect candidates
    candidates: list[Path] = []
    for p in base.iterdir():
        if not p.is_file() or p.suffix.lower() != ".xlsx":
            continue
        name = p.name if case_sensitive else p.name.lower()
        needle = type_of_assembly if case_sensitive else type_of_assembly.lower()
        if needle in name:
            candidates.append(p)

    if not candidates:
        raise FileNotFoundError(
            f"No *.xlsx file containing {type_of_assembly!r} found in {base}."
        )
    if len(candidates) > 1:
        # Liste kurzhalten
        listed = ", ".join(c.name for c in candidates[:5])
        more = " …" if len(candidates) > 5 else ""
        raise ValueError(f"More than one matching file found in {base}: {listed}{more}")

    LOGGER.info("Found planning file: %s", candidates[0])
    return candidates[0]


def transparentify(
    cmap: colors.Colormap | str,
    *,
    start_alpha: float = 0.0,
    end_alpha: float = 1.0,
) -> colors.ListedColormap:
    """Return a copy of `cmap` with a linear alpha ramp.

    Parameters
    ----------
    cmap : matplotlib.colors.Colormap or str
        Colormap instance or name, e.g. ``cm.Greys`` or ``"Greys"``.
    start_alpha : float, optional
        Alpha at the first color (default 0.0).
    end_alpha : float, optional
        Alpha at the last color (default 1.0).

    Returns
    -------
    matplotlib.colors.ListedColormap
        A ListedColormap with the same RGB values and updated alpha.

    Raises
    ------
    ValueError
        If alpha values are not in ``[0.0, 1.0]``.
    """
    if isinstance(cmap, str):
        base = cm.get_cmap(cmap)
    else:
        base = cmap

    if not (0.0 <= start_alpha <= 1.0 and 0.0 <= end_alpha <= 1.0):
        raise ValueError("start_alpha and end_alpha must be within [0.0, 1.0].")

    rgba = np.array(base(np.arange(base.N)))
    rgba[:, 3] = np.linspace(start_alpha, end_alpha, base.N)
    name = getattr(base, "name", "transparentified")
    return colors.ListedColormap(rgba, name=f"{name}_alpha")


class FZcolors:
    """Colors from the FZJ corporate design."""

    red = np.array((191, 21, 33)) / 255
    green = np.array((0, 153, 102)) / 255
    blue = np.array((2, 61, 107)) / 255
    orange = np.array((220, 110, 0)) / 255


class FZcmaps:
    """Color maps correspoding to FZJ corporate design colors."""

    Reds = to_colormap(FZcolors.red)
    Greens = to_colormap(FZcolors.green)
    Blues = to_colormap(FZcolors.blue)
    Oranges = to_colormap(FZcolors.orange)
    Greys = transparentify(cm.Greys)


def calculate_gg_mastermix_components(
    t4_concentration: int = 500,
    rest_enzyme_concentration: int = 20,
    t4_units: int = 10,
    rest_enzyme_units: int = 10,
    mastermix_volume: int = 700,
    volume_per_reaction: float = 7.5,
    bsa_100x: bool = False,
) -> pandas.DataFrame:
    """Compute Golden Gate mastermix components per reaction and for the total mix.

    The per-reaction recipe is:
    - 10x ligase buffer: 0.5 × reaction volume
    - T4 DNA ligase: ``t4_units / t4_concentration`` (µL)
    - 100x BSA: 0.02 × reaction volume (if ``bsa_100x``)
    - 10 U BsaI: ``rest_enzyme_units / rest_enzyme_concentration`` (µL)
    - H2O: to reach ``volume_per_reaction``

    Parameters
    ----------
    t4_concentration : int, optional
        Concentration of T4 DNA ligase (U/µL), default 500.
    rest_enzyme_concentration : int, optional
        Concentration of restriction enzyme (U/µL), default 20.
    t4_units : int, optional
        Desired units of T4 ligase per reaction, default 10.
    rest_enzyme_units : int, optional
        Desired units of restriction enzyme per reaction, default 10.
    mastermix_volume : int, optional
        Total mastermix volume to prepare (µL), default 700.
    volume_per_reaction : float, optional
        Reaction volume used for scaling (µL), default 7.5.
    bsa_100x : bool, optional
        Include 100× BSA at 1:100 (i.e., 0.02 × reaction volume), default False.
        It is not the enzyme bsaI that is meant, but bsa, the protein that is important for freezing.

    Returns
    -------
    pandas.DataFrame
        Table with columns: ``Component``, ``Per Reaction (µL)``, ``Total Volume (µL)``.

    Raises
    ------
    ValueError
        If any concentration or ``volume_per_reaction`` is non-positive.
    """
    if t4_concentration <= 0 or rest_enzyme_concentration <= 0:
        raise ValueError("Enzyme concentrations must be positive.")
    if volume_per_reaction <= 0:
        raise ValueError("volume_per_reaction must be positive.")
    if mastermix_volume < 0:
        raise ValueError("mastermix_volume must be non-negative.")

    per_reaction: dict[str, float] = {
        "10x ligase buffer": round(volume_per_reaction * 0.5, 2),
        "T4 DNA ligase": round(t4_units / float(t4_concentration), 2),
        "100x BSA": round(volume_per_reaction * 0.02, 2) if bsa_100x else 0.0,
        "10 U BsaI": round(rest_enzyme_units / float(rest_enzyme_concentration), 2),
        "H2O": 0.0,  # computed below
    }

    total_reactions = mastermix_volume / float(volume_per_reaction)

    # Water to top up: ensure non-negative after rounding errors
    water = volume_per_reaction - sum(per_reaction.values())
    per_reaction["H2O"] = round(max(0.0, water), 2)

    total_volume = {
        comp: round(val * total_reactions, 2) for comp, val in per_reaction.items()
    }

    return pandas.DataFrame(
        {
            "Component": list(per_reaction.keys()),
            "Per Reaction (µL)": list(per_reaction.values()),
            "Total Volume (µL)": list(total_volume.values()),
        }
    )


def calculate_GGmastermix_components(*args, **kwargs):  # noqa: N802
    """Compatibility wrapper for the old camelCase name.

    Prefer :func:`calculate_gg_mastermix_components`.
    """
    warnings.warn(
        "calculate_GGmastermix_components() is deprecated; "
        "use calculate_gg_mastermix_components() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return calculate_gg_mastermix_components(*args, **kwargs)


def calculate_pcr_mastermix_components(
    volume_per_reaction: float = 50.0,
    *,
    primer_final_um: float = 0.5,
    primer_stock_um: float = 10.0,
    template_volume: float = 1.0,
    dmso_fraction: float = 0.0,  # e.g. 0.05 for 5% v/v DMSO (Coryne cPCR)
    mastermix_volume: float = 700.0,
    mastermix_2x_fraction: float = 0.5,  # 2X mix ⇒ 50% of reaction volume
    polymerase_name: str = "Q5 High-Fidelity 2X Master Mix",
    extra_components_per_rxn: dict[str, float] | None = None,
) -> pandas.DataFrame:
    """Compute PCR mastermix components per reaction and for the total mix.

    Defaults are set for Q5 2× Master Mix (50 µL Reaction Volume, 0.5 µM Primer final, 10 µM Primer-Stock).
    DMSO can optionally be specified as a v/v fraction (e.g. 0.05 for 5%).

    Parameters
    ----------
    volume_per_reaction : float, optional
        Total volume per reaction in µL. Default 50.0.
    primer_final_uM : float, optional
        Desired final concentration per primer (µM), Default 0.5.
    primer_stock_uM : float, optional
        Stock solution concentration primer (µM), Default 10.0.
    template_volume : float, optional
        Volume of template DNA per reaction (µL), Default 1.0.
    dmso_fraction : float, optional
        DMSO v/v fraction (0.0–0.2 sensible). 0.05 equals 5 %. Default 0.0.
    mastermix_volume : float, optional
        Total volume of the master mix to be prepared (µL), Default 700.0.
    mastermix_2x_fraction : float, optional
        Proportion of the 2× master mix in the reaction volume (typical 0.5), Default 0.5.
    polymerase_name : str, optional
        Display name of the 2× mixture in the output. Default
        "Q5 High-Fidelity 2X Master Mix".
    extra_components_per_rxn : dict[str, float], optional
        Optional additional components (µL per reaction), e.g. {"MgCl2 (25 mM)": 0.5}.

    Returns
    -------
    pandas.DataFrame
        Table with columns: ``Component``, ``Per Reaction (µL)``, ``Total Volume (µL)``.

    Raises
    ------
    ValueError
        In the event of invalid (non-positive) input parameters,
        or if the sum of the partial volumes exceeds the reaction volume.
    """
    # --- Validation ---
    if volume_per_reaction <= 0:
        raise ValueError("volume_per_reaction must be positive.")
    if primer_stock_um <= 0:
        raise ValueError("primer_stock_uM must be positive.")
    if not (0.0 <= dmso_fraction < 1.0):
        raise ValueError("dmso_fraction must be in [0.0, 1.0).")
    if not (0.0 < mastermix_2x_fraction <= 1.0):
        raise ValueError("mastermix_2x_fraction must be in (0.0, 1.0].")
    if mastermix_volume < 0:
        raise ValueError("mastermix_volume must be non-negative.")
    if template_volume < 0:
        raise ValueError("template_volume must be non-negative.")

    # --- Per-Reaction Calculation ---
    v_mastermix_2x = round(volume_per_reaction * mastermix_2x_fraction, 2)
    v_primer_fwd = round((primer_final_um * volume_per_reaction) / primer_stock_um, 2)
    v_primer_rev = v_primer_fwd
    v_template = round(template_volume, 2)
    v_dmso = round(volume_per_reaction * dmso_fraction, 2)

    extra = extra_components_per_rxn or {}
    for k, v in extra.items():
        if v < 0:
            raise ValueError(f"Extra component {k!r} volume must be non-negative.")

    known_sum = (
        v_mastermix_2x
        + v_primer_fwd
        + v_primer_rev
        + v_template
        + v_dmso
        + sum(extra.values())
    )
    v_water = round(volume_per_reaction - known_sum, 2)

    if v_water < -1e-6:
        raise ValueError(
            f"Component volumes ({known_sum} µL) exceed reaction volume "
            f"({volume_per_reaction} µL). Decrease template/DMSO/extra components."
        )
    # truncate small negative values ​​to 0
    v_water = max(0.0, v_water)

    per_rxn: dict[str, float] = {
        polymerase_name: v_mastermix_2x,
        "Primer Fwd": v_primer_fwd,
        "Primer Rev": v_primer_rev,
        "Template DNA": v_template,
        "DMSO": v_dmso,
        **extra,
        "Nuclease-free water": v_water,
    }

    # --- Scaling to total volume ---
    n_rxn = (
        mastermix_volume / float(volume_per_reaction) if volume_per_reaction else 0.0
    )
    total = {k: round(v * n_rxn, 2) for k, v in per_rxn.items()}

    return pandas.DataFrame(
        {
            "Component": list(per_rxn.keys()),
            "Per Reaction (µL)": list(per_rxn.values()),
            "Total Volume (µL)": list(total.values()),
        }
    )


def calculate_PCRmastermix_components(*args, **kwargs):  # noqa: N802
    """Compatibility wrapper for the old camelCase name (prefer snake_case)."""
    warnings.warn(
        "calculate_PCRmastermix_components() is deprecated; "
        "use calculate_pcr_mastermix_components() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return calculate_pcr_mastermix_components(*args, **kwargs)


def calculate_dreamtaq_mastermix_components(
    volume_per_reaction: float = 50.0,
    *,
    primer_final_um: float = 0.5,
    primer_stock_um: float = 10.0,
    template_volume: float = 1.0,
    dmso_fraction: float = 0.0,  # 0.05 ⇒ 5% v/v DMSO (optional)
    mastermix_volume: float = 700.0,
    mastermix_2x_fraction: float = 0.5,  # 2X mix ⇒ 50% of reaction volume
    polymerase_name: str = "DreamTaq PCR Master Mix (2X)",
    extra_components_per_rxn: dict[str, float] | None = None,
) -> pandas.DataFrame:
    """Compute DreamTaq 2× PCR mastermix per reaction and total.

    Defaults match the vendor recipe for a 50 µL reaction:
    25 µL 2× mix, 0.5 µM Primer (from 10 µM Stock), X µL Template, water to volume.
    DMSO can optionally be specified as a v/v fraction (e.g., 0.05 for 5%).

    Parameters
    ----------
    volume_per_reaction : float, optional
        Total reaction volume (µL). Default 50.0.
    primer_final_uM : float, optional
        Target final concentration per primer (µM). Default 0.5.
    primer_stock_uM : float, optional
        Primer stock concentration (µM). Default 10.0.
    template_volume : float, optional
        Template DNA volume per reaction (µL). Default 1.0.
    dmso_fraction : float, optional
        DMSO v/v fraction (0.0–1.0). 0.05 = 5 %. Default 0.0.
    mastermix_volume : float, optional
        Total mastermix to prepare (µL). Default 700.0.
    mastermix_2x_fraction : float, optional
        Fraction of 2× mix in the reaction volume. Default 0.5.
    polymerase_name : str, optional
        Display name for the 2× mix. Default "DreamTaq PCR Master Mix (2X)".
    extra_components_per_rxn : dict[str, float], optional
        Additional components in µL per reaction (e.g. {"MgCl2 (25 mM)": 0.5}).

    Returns
    -------
    pandas.DataFrame
        Columns: ``Component``, ``Per Reaction (µL)``, ``Total Volume (µL)``.
    """
    return calculate_pcr_mastermix_components(
        volume_per_reaction=volume_per_reaction,
        primer_final_um=primer_final_um,
        primer_stock_um=primer_stock_um,
        template_volume=template_volume,
        dmso_fraction=dmso_fraction,
        mastermix_volume=mastermix_volume,
        mastermix_2x_fraction=mastermix_2x_fraction,
        polymerase_name=polymerase_name,
        extra_components_per_rxn=extra_components_per_rxn,
    )


def calculate_DreamTaqMastermix_components(*args, **kwargs):  # noqa: N802
    """Compatibility wrapper (prefer snake_case)."""
    warnings.warn(
        "calculate_DreamTaqMastermix_components() is deprecated; "
        "use calculate_dreamtaq_mastermix_components() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return calculate_dreamtaq_mastermix_components(*args, **kwargs)


def calculate_gibson_mastermix_components(
    total_mix_volume: float = 1200.0,
) -> pandas.DataFrame:
    """
    Return the Gibson Assembly master mix scaled to ``total_mix_volume`` (µL).

    The scaling is derived from the following baseline recipe which totals **1200 µL**:
      - 5× isothermal reaction buffer: **320.00 µL**
      - T5 Exonuclease (10 U/µL): **0.64 µL**
      - Phusion DNA Polymerase (2 U/µL): **20.00 µL**
      - Taq DNA ligase (40 U/µL): **160.00 µL**
      - H₂O: **699.36 µL** (computed as the remainder)

    Parameters
    ----------
    total_mix_volume : float, optional
        Desired total volume of the master mix in µL. Must be positive.
        Default is ``1200.0``.

    Returns
    -------
    pandas.DataFrame
        A two-column table with:
        - ``Component``: reagent name
        - ``Total Volume (µL)``: scaled volume in µL (rounded to 2 decimals)

    Raises
    ------
    ValueError
        If ``total_mix_volume`` is not positive.

    Notes
    -----
    - Water is calculated as the remainder so that all components sum to
      exactly ``total_mix_volume``.
    - Values are rounded to two decimal places for pipetting practicality.

    Examples
    --------
    >>> df = calculate_gibson_mastermix_components(1200)
    >>> df
                      Component  Total Volume (µL)
    0  5× isothermal reaction buffer            320.00
    1     T5 Exonuclease (10 U/µL)               0.64
    2  Phusion DNA Polymerase (2 U/µL)          20.00
    3        Taq DNA ligase (40 U/µL)          160.00
    4                          H2O             699.36
    """
    if total_mix_volume <= 0:
        raise ValueError("total_mix_volume must be positive.")

    baseline_total = 1200.0
    baseline = {
        "5× isothermal reaction buffer": 320.00,
        "T5 Exonuclease (10 U/µL)": 0.64,
        "Phusion DNA Polymerase (2 U/µL)": 20.00,
        "Taq DNA ligase (40 U/µL)": 160.00,
    }

    scale = total_mix_volume / baseline_total
    per_total = {k: round(v * scale, 2) for k, v in baseline.items()}
    per_total["H2O"] = round(total_mix_volume - sum(per_total.values()), 2)

    return pandas.DataFrame(
        {
            "Component": list(per_total.keys()),
            "Total Volume (µL)": list(per_total.values()),
        }
    )


# Optional: Legacy alias (CamelCase), in case you want consistency, like with GG/Taq/Q5
def calculate_GibsonMastermix_components(*args, **kwargs):  # noqa: N802
    """Deprecated alias; use :func:`calculate_gibson_mastermix_components`."""
    import warnings

    warnings.warn(
        "calculate_GibsonMastermix_components() is deprecated; "
        "use calculate_gibson_mastermix_components() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return calculate_gibson_mastermix_components(*args, **kwargs)
