from __future__ import annotations

import datetime
import logging
import os
import re
import textwrap
from collections import defaultdict
from collections.abc import Callable, Iterable, Sequence
from pathlib import Path
from typing import Any, Optional, TypeAlias

import matplotlib.pyplot as plt
import numpy as np
import openpyxl
import pandas as pd
import robotools
from IPython.display import Image
from matplotlib.patches import Rectangle

from assembly_designer.utils import (
    check_and_log_volumes,
    fn_plot2,
    plot_gif,
    setup_logging,
    setup_run_logging_files,
)

LOGGER = logging.getLogger(__name__)
setup_logging(level="INFO", name="assembly_designer", propagate=False)

# File-Logging for this run (full + Warn/Error)
_RUN_DIR = os.path.join("logs", datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S"))
setup_run_logging_files(log_dir=_RUN_DIR)  # writes assembly.log + assembly_warnerr.log

WellValue: TypeAlias = str | float | int | None


def normalize_well_id(value: WellValue) -> WellValue:
    """Normalize Opentrons well IDs like 'A01' -> 'A1'.

    Values that are empty, NaN, or not valid 96-well IDs are returned unchanged.
    """
    if value is None:
        return value

    if isinstance(value, float) and pd.isna(value):
        return value

    text = str(value).strip()
    match: re.Match[str] | None = re.fullmatch(
        r"([A-Ha-h])0*([1-9]|1[0-2])",
        text,
    )
    if match is None:
        return value

    row, col = match.groups()
    return f"{row.upper()}{int(col)}"


def normalize_worklist_well_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with known worklist well columns normalized.

    The function only modifies columns if they are present.
    """
    out = df.copy()

    candidate_columns: list[str] = [
        "Well Nr. source plate",
        "Well Nr. destination plate",
        "Well",
        "Source Well",
        "Destination Well",
        "source_well",
        "destination_well",
    ]

    for col in candidate_columns:
        if col in out.columns:
            out[col] = out[col].apply(normalize_well_id)

    return out


def build_master_worklist(worklist_data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Combine all individual worklists into one master worklist.

    Each row gets a leading ``Source_Worklist`` column so the Opentrons script
    can reconstruct the original worklists later.
    """
    frames: list[pd.DataFrame] = []

    for worklist_name, df in worklist_data.items():
        if not isinstance(df, pd.DataFrame) or df.empty:
            continue

        part = normalize_worklist_well_columns(df)
        part.insert(0, "Source_Worklist", str(worklist_name))
        frames.append(part)

    if not frames:
        return pd.DataFrame()

    return pd.concat(frames, ignore_index=True)


def find_worklist_sheets(
    file_path: str | os.PathLike[str],
    pattern: str = "Worklist",
    case_sensitive: bool = False,
    regex: bool = False,
) -> list[str]:
    """Find sheet names in an Excel workbook matching a pattern.

    Parameters
    ----------
    file_path : str or os.PathLike
        Path to the Excel file.
    pattern : str, optional
        Substring or regular expression used to match sheet names.
        Default is ``"Worklist"``.
    case_sensitive : bool, optional
        If ``True``, matching is case-sensitive. Default ``False``.
    regex : bool, optional
        If ``True``, interpret `pattern` as a regular expression.
        Default ``False`` (simple substring match).

    Returns
    -------
    list of str
        Sheet names that match `pattern`.

    Raises
    ------
    FileNotFoundError
        If `file_path` does not exist.
    ValueError
        If the file cannot be opened as a valid Excel workbook.
    """
    path = os.fspath(file_path)
    if not os.path.exists(path):
        raise FileNotFoundError(f"The file '{path}' does not exist.")

    try:
        wb = openpyxl.load_workbook(path, read_only=True)
    except Exception as err:  # pylint: disable=broad-except
        raise ValueError(f"Failed to open Excel file: {err}") from err

    names = wb.sheetnames
    if regex:
        flags = 0 if case_sensitive else re.IGNORECASE
        matched = [s for s in names if re.search(pattern, s, flags)]
    else:
        if case_sensitive:
            matched = [s for s in names if pattern in s]
        else:
            p = pattern.lower()
            matched = [s for s in names if p in s.lower()]

    LOGGER.info("Found %d matching sheet(s): %s", len(matched), matched)
    return matched


def load_worklist_sheets(
    file_name: str,
    sheet_names: list[str],
    usecols: Iterable[int] | None = (0, 1, 2, 3),
    dropna: bool = True,
) -> dict[str, pd.DataFrame]:
    """Load specific sheets from an Excel file into DataFrames.

    The function validates that the file exists and that all requested
    sheet names are present. Each sheet is read (optionally restricted
    to `usecols`) and rows with missing values can be dropped.

    Parameters
    ----------
    file_name : str
        Path to the Excel file.
    sheet_names : list of str
        Sheet names to load.
    usecols : iterable of int or None, optional
        Column indices to read (e.g., ``(0, 1, 2, 3)``). If ``None``,
        all columns are read. Default is ``(0, 1, 2, 3)``.
    dropna : bool, optional
        If ``True``, drop rows containing any NaN values. Default ``True``.

    Returns
    -------
    dict[str, pd.DataFrame]
        Mapping of sheet name to loaded DataFrame.

    Raises
    ------
    FileNotFoundError
        If the Excel file does not exist.
    ValueError
        If the file cannot be opened, required sheets are missing,
        or a sheet fails to load.
    """
    if not os.path.exists(file_name):
        raise FileNotFoundError(f"The file '{file_name}' does not exist.")

    try:
        xls = pd.ExcelFile(file_name)
    except Exception as err:  # pylint: disable=broad-except
        raise ValueError(f"Failed to open Excel file: {err}") from err

    missing = [s for s in sheet_names if s not in xls.sheet_names]
    if missing:
        raise ValueError(f"Unknown sheet(s): {missing}. Available: {xls.sheet_names}")

    result: dict[str, pd.DataFrame] = {}
    for sheet_name in sheet_names:
        try:
            if usecols is None:
                df = pd.read_excel(
                    file_name, sheet_name=sheet_name, header=0, index_col=None
                )
            else:
                df = pd.read_excel(
                    file_name,
                    sheet_name=sheet_name,
                    usecols=list(usecols),
                    header=0,
                    index_col=None,
                )
            if dropna:
                df = df.dropna()
            result[sheet_name] = df
            LOGGER.info(
                "Loaded sheet '%s' (%d rows, %d cols).",
                sheet_name,
                df.shape[0],
                df.shape[1],
            )
        except Exception as err:  # pylint: disable=broad-except
            raise ValueError(f"Failed to load sheet '{sheet_name}': {err}") from err

    return result


def create_empty_gwl_files(file_names: list[str]) -> None:
    """Create empty ``.gwl`` files for each provided name.

    The function appends the ``.gwl`` suffix if missing, ensures parent
    directories exist, and creates (or truncates) each file. It logs progress
    and raises on OS-level errors.

    Parameters
    ----------
    file_names : list of str
        File names with or without the ``.gwl`` extension.

    Returns
    -------
    None
        The function creates files but does not return a value.

    Raises
    ------
    OSError
        If any file cannot be created due to an OS-level error.
    """
    for name in file_names:
        path = name if name.endswith(".gwl") else f"{name}.gwl"
        directory = os.path.dirname(path) or "."
        os.makedirs(directory, exist_ok=True)

        try:
            # Create or truncate the file; avoid binding an unused variable.
            with open(path, "w", encoding="utf-8"):
                pass
            LOGGER.info("Created GWL file: %s", path)
        except OSError as err:
            LOGGER.error("Failed to create GWL file %s: %s", path, err)
            # Preserve original traceback (Ruff B904).
            raise OSError(f"Failed to create GWL file {path}") from err

    LOGGER.info("All .gwl files were successfully created.")


# def load_metadata_dataframe(file_name, sheet_name):
#     """
#     Automatically searches for the "Meta data" column in an Excel sheet and loads
#     all rows beneath it, along with its associated value column, into a DataFrame.

#     Parameters
#     ----------
#     file_name : str
#         The path to the Excel file.
#     sheet_name : str
#         The name of the sheet to search.

#     Returns
#     -------
#     pd.DataFrame
#         A DataFrame containing the "Meta data" and associated value columns.

#     Raises
#     ------
#     ValueError
#         If the "Meta data" column cannot be found.
#     FileNotFoundError
#         If the file or sheet does not exist.
#     """
#     try:
#         # Load the entire sheet
#         sheet_df = pd.read_excel(file_name, sheet_name=sheet_name, header=None)

#         # Find the "Meta data" header
#         meta_data_location = sheet_df.apply(lambda row: row.astype(str).str.contains("Meta data", na=False)).any(axis=1)

#         if not meta_data_location.any():
#             raise ValueError("'Meta data' column not found in the sheet.")

#         # Get the row and column index of "Meta data"
#         start_row = meta_data_location.idxmax()
#         meta_data_col = sheet_df.iloc[start_row].tolist().index("Meta data")

#         # Assume the next column is the value column
#         value_col = meta_data_col + 1

#         # Extract the relevant rows and columns
#         metadata_df = sheet_df.iloc[start_row + 1:, [meta_data_col, value_col]].dropna()

#         # Rename columns for clarity
#         metadata_df.columns = ["Meta data", "Value"]

#         return metadata_df

#     except FileNotFoundError:
#         raise FileNotFoundError(f"File not found: {file_name}")
#     except Exception as e:
#         raise ValueError(f"An error occurred while processing the file: {e}")


def col_letter_to_index(col_letters: str) -> int:
    """Convert Excel-style column letters to a zero-based integer index.

    Examples
    --------
    >>> col_letter_to_index("A")
    0
    >>> col_letter_to_index("Z")
    25
    >>> col_letter_to_index("AA")
    26
    >>> col_letter_to_index("AB")
    27

    Parameters
    ----------
    col_letters : str
        One or multiple uppercase/lowercase ASCII letters representing
        an Excel column (e.g., ``"A"``, ``"Z"``, ``"AA"``).

    Returns
    -------
    int
        Zero-based column index.

    Raises
    ------
    ValueError
        If `col_letters` is empty or contains non-letter characters.
    """
    if not col_letters or not col_letters.strip():
        raise ValueError("Column letters must be a non-empty string.")

    letters = col_letters.strip().upper()
    idx = 0
    for ch in letters:
        if not ("A" <= ch <= "Z"):
            raise ValueError(
                f"Invalid column letter {col_letters!r}: only A–Z are allowed."
            )
        idx = idx * 26 + (ord(ch) - ord("A") + 1)

    # Convert 1-based to 0-based
    return idx - 1


def load_metadata_dataframe(file_name: str, sheet_name: str) -> pd.DataFrame:
    """Extract a two-column \"Meta data\" block from an Excel sheet.

    The function searches the given sheet **case-insensitively** for a cell
    whose normalized content equals ``\"meta data\"``. It then returns the
    two columns starting at that cell: the label column (``\"Meta data\"``)
    and the value column immediately to its right (``\"Value\"``),
    reading all subsequent rows until the first fully-empty row is reached.

    Examples
    --------
    >>> # df = load_metadata_dataframe("design.xlsx", "Sheet1")
    >>> # df.head()
    >>> #        Meta data           Value
    >>> # 0      Project ID           AD-42
    >>> # 1     Organism ID  C. glutamicum
    >>> # 2       Designer      T. Stoltmann

    Parameters
    ----------
    file_name : str
        Path to the Excel file.
    sheet_name : str
        Name of the sheet to read.

    Returns
    -------
    pd.DataFrame
        DataFrame with exactly two columns: ``\"Meta data\"`` and ``\"Value\"``.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the sheet cannot be read, the \"Meta data\" header cannot be found,
        or the value column is missing.
    """
    try:
        sheet_df = pd.read_excel(file_name, sheet_name=sheet_name, header=None)
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"File not found: {file_name}") from exc
    except Exception as exc:  # pylint: disable=broad-except
        raise ValueError(f"Failed to read sheet {sheet_name!r}: {exc}") from exc

    # Normalize for robust matching
    norm = sheet_df.applymap(lambda x: str(x).strip().lower() if pd.notna(x) else "")

    # Find the first row containing an exact "meta data" cell
    start_row = None
    start_col = None
    for r_idx in range(norm.shape[0]):
        for c_idx in range(norm.shape[1]):
            if norm.iat[r_idx, c_idx] == "meta data":
                start_row = r_idx
                start_col = c_idx
                break
        if start_row is not None:
            break

    if start_row is None or start_col is None:
        raise ValueError("Could not locate a 'Meta data' header in the sheet.")

    value_col = start_col + 1
    if value_col >= sheet_df.shape[1]:
        raise ValueError(
            "'Meta data' header is the last column; a value column to the right is required."
        )

    # Slice rows below the header and the two columns (label + value)
    out = sheet_df.iloc[start_row + 1 :, [start_col, value_col]].copy()

    # Stop at the first fully empty row (optional; keeps tidy tables)
    first_empty = None
    for r_idx in range(out.shape[0]):
        if out.iloc[r_idx].isna().all():
            first_empty = r_idx
            break
    if first_empty is not None:
        out = out.iloc[:first_empty]

    # Final cleanup
    out.columns = ["Meta data", "Value"]
    out = out.dropna(how="all")
    LOGGER.debug("Loaded %d metadata rows from %s[%s]", len(out), file_name, sheet_name)
    return out.reset_index(drop=True)


def plot_mikrotiterplatte(
    df: pd.DataFrame,
    title: str = "Mikrotiterplatte",
    highlight_column: Optional[str] = None,
    cmap: str = "tab10",
    color_palette: Optional[Sequence[Any]] = None,
    max_lines: int = 3,
    font_size: int = 10,
    wrap_width: int = 20,
    save_path: Optional[str] = None,
    exclude_keywords: Optional[Sequence[str]] = None,
) -> None:
    """Plot a 96-well microtiter plate (8×12) with well annotations.

    Optionally highlights wells by unique values in `highlight_column`, accepts a
    custom `color_palette`, and can exclude columns before rendering.

    Parameters
    ----------
    df : pandas.DataFrame
        Input table with at least a ``"Well"`` column containing alphanumeric
        well IDs (e.g., ``"A1"`` or ``"A01"``). Other columns are rendered as
        per-well text (wrapped/limited).
    title : str, optional
        Title for the figure. Default is ``"Mikrotiterplatte"``.
    highlight_column : str or None, optional
        Column used to color wells by unique values. Default ``None``.
    cmap : str, optional
        Matplotlib colormap name used if no `color_palette` is provided.
        Default ``"tab10"``.
    color_palette : Sequence or None, optional
        Custom colors (cycled) to map unique values in `highlight_column`.
        If ``None``, colors are drawn from `cmap`. Default ``None``.
    max_lines : int, optional
        Maximum number of text lines per well. Default ``3``.
    font_size : int, optional
        Font size for well annotations. Default ``10``.
    wrap_width : int, optional
        Soft wrap width (characters) for per-well text. Default ``20``.
    save_path : str or None, optional
        If provided, path to save the plot (parent directory is created).
        If ``None``, the figure is shown. Default ``None``.
    exclude_keywords : Sequence[str] or None, optional
        Case-insensitive list of column *names* to exclude (exact match).
        Default ``None``.

    Returns
    -------
    None

    Raises
    ------
    KeyError
        If the ``"Well"`` column is missing.
    """
    # --- Validate input
    if "Well" not in df.columns:
        raise KeyError("The 'Well' column is missing from the DataFrame.")

    # Filter valid wells only
    df_valid = df[df["Well"].notna()].copy()

    # Exclude columns by exact name (case-insensitive)
    if exclude_keywords:
        exclude_lower = {kw.lower() for kw in exclude_keywords}
        keep_mask = ~df_valid.columns.str.lower().isin(exclude_lower)
        removed = df_valid.columns[~keep_mask].tolist()
        df_valid = df_valid.loc[:, keep_mask]
        if removed:
            LOGGER.info("Excluded columns: %s", removed)

    if df_valid.empty:
        LOGGER.warning("DataFrame is empty after filtering; nothing to plot.")
        return

    # Plate geometry
    rows = list("ABCDEFGH")
    cols = [f"{i:02d}" for i in range(1, 13)]
    plate = pd.DataFrame(index=rows, columns=cols, data="")

    def _parse_well(value: Any) -> Optional[tuple[str, str]]:
        """Return (row, col) as ('A'..'H', '01'..'12') or None if invalid."""
        if not isinstance(value, str) or len(value) < 2:
            return None
        row_id = value[0].upper()
        try:
            col_num = int(value[1:])
        except ValueError:
            return None
        if row_id not in rows or not (1 <= col_num <= 12):
            return None
        return row_id, f"{col_num:02d}"

    # Compose well text (wrapped, limited)
    for _, r in df_valid.iterrows():
        parsed = _parse_well(r["Well"])
        if parsed is None:
            continue
        r_id, c_id = parsed
        contents = [
            str(r[c])
            for c in df_valid.columns
            if c != "Well" and pd.notna(r[c]) and str(r[c]).strip() not in {"", " "}
        ]
        if not contents:
            plate.at[r_id, c_id] = ""
        else:
            wrapped = [textwrap.fill(s, width=wrap_width) for s in contents[:max_lines]]
            plate.at[r_id, c_id] = "\n".join(wrapped)

    # Build color mapping
    well_colors: dict[str, Any] = {}
    if highlight_column and highlight_column in df_valid.columns:
        unique_vals = pd.Series(df_valid[highlight_column]).dropna().unique().tolist()
        if color_palette and len(color_palette) > 0:
            palette = list(color_palette)
        else:
            cmap_obj = plt.get_cmap(cmap)
            # Sample evenly from the colormap
            palette = [
                cmap_obj(i / max(1, len(unique_vals) - 1))
                for i in range(len(unique_vals))
            ]

        value_to_color = {
            val: palette[i % len(palette)] for i, val in enumerate(unique_vals)
        }
        for _, r in df_valid.iterrows():
            parsed = _parse_well(r["Well"])
            if parsed is None:
                continue
            r_id, c_id = parsed
            well_colors[f"{r_id}{c_id}"] = value_to_color.get(
                r.get(highlight_column), "white"
            )

    # --- Plot
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_title(title, fontsize=16)
    ax.axis("off")

    for r_idx, r_id in enumerate(rows):
        for c_idx, c_id in enumerate(cols):
            well_id = f"{r_id}{c_id}"
            rect_color = "white"
            cell_text = plate.at[r_id, c_id]
            if isinstance(cell_text, str) and cell_text.strip():
                rect_color = well_colors.get(well_id, "white")
            ax.add_patch(
                Rectangle((c_idx, r_idx), 1, 1, facecolor=rect_color, edgecolor="black")
            )
            if cell_text:
                ax.text(
                    c_idx + 0.5,
                    r_idx + 0.5,
                    cell_text,
                    ha="center",
                    va="center",
                    fontsize=font_size,
                    wrap=True,
                )

    # Axis cosmetics
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 8)
    ax.set_xticks(np.arange(12) + 0.5)
    ax.set_yticks(np.arange(8) + 0.5)
    ax.set_xticklabels(cols, fontsize=10)
    ax.set_yticklabels(rows, fontsize=10)
    ax.invert_yaxis()

    # Save or show
    if save_path:
        os.makedirs(os.path.dirname(save_path) or ".", exist_ok=True)
        plt.savefig(save_path, bbox_inches="tight")
        LOGGER.info("Plate plot saved to: %s", save_path)
        plt.close(fig)
    else:
        plt.show()


def process_worklist(
    mtp_manager: Any,
    worklist_data: dict[str, pd.DataFrame],
    worklist_name: str,
    source_plate_name: str,
    destination_plate_name: str,
    list_liquidhandling_data: list[tuple[str, Any]],
    worklist_file: str = "Worklist.gwl",
    label: str = "transfer",
    plot: bool = True,
    plot_function: Optional[Callable[[Any, str], Any]] = None,
) -> None:
    """Process a transfer worklist: source → destination, then plot destination volumes.

    The function validates inputs, writes a Fluent `.gwl` worklist via `robotools`,
    updates the `destination_plate` history summary, and (optionally) plots the
    resulting fractional volumes.

    Parameters
    ----------
    mtp_manager : object
        Manager that provides `get_mtp(name)` and returns plate objects compatible
        with `robotools.FluentWorklist`.
    worklist_data : dict[str, pandas.DataFrame]
        Mapping from worklist name to table with columns:
        ``"Well Nr. source plate"``, ``"Well Nr. destination plate"``, ``"Volumen"``.
    worklist_name : str
        Key to select the DataFrame from `worklist_data`.
    source_plate_name : str
        Name of the source plate.
    destination_plate_name : str
        Name of the destination plate.
    list_liquidhandling_data : list[tuple[str, Any]]
        Log/list to append a summary `(title, data)` after execution.
    worklist_file : str, optional
        Path to write the `.gwl` file. Default ``"Worklist.gwl"``.
    label : str, optional
        Label written into the worklist. Default ``"transfer"``.
    plot : bool, optional
        If ``True``, create a plot of destination volumes. Default ``True``.
    plot_function : callable or None, optional
        Function like ``fn_plot(volumes, title)``. If ``None``, no plot is created
        (unless caller injects a function).

    Returns
    -------
    None

    Raises
    ------
    KeyError
        If `worklist_name` is missing in `worklist_data`.
    ValueError
        If plates are missing, required columns are missing, or columns have
        inconsistent lengths after dropping NaNs.
    OSError
        If the worklist file cannot be written.
    """
    if worklist_name not in worklist_data:
        raise KeyError(f"Worklist '{worklist_name}' not found in worklist_data.")
    df = worklist_data[worklist_name]

    source_plate = mtp_manager.get_mtp(source_plate_name)
    destination_plate = mtp_manager.get_mtp(destination_plate_name)
    if source_plate is None:
        raise ValueError(f"Source plate '{source_plate_name}' not found in MTPManager.")
    if destination_plate is None:
        raise ValueError(
            f"Destination plate '{destination_plate_name}' not found in MTPManager."
        )

    required = {
        "Well Nr. source plate",
        "Well Nr. destination plate",
        "Volumen",
    }
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required column(s): {missing}")

    src_wells = df["Well Nr. source plate"].dropna().astype(str).tolist()
    dst_wells = df["Well Nr. destination plate"].dropna().astype(str).tolist()
    volumes = df["Volumen"].dropna().astype(float).tolist()

    if not (len(src_wells) == len(dst_wells) == len(volumes)):
        raise ValueError(
            "Source wells, destination wells, and volumes must have equal lengths "
            f"(got {len(src_wells)}, {len(dst_wells)}, {len(volumes)})."
        )

    # Ensure directory for worklist exists
    os.makedirs(os.path.dirname(worklist_file) or ".", exist_ok=True)

    try:
        with robotools.FluentWorklist(worklist_file) as wl:
            wl.transfer(
                source=source_plate,
                source_wells=src_wells,
                destination=destination_plate,
                destination_wells=dst_wells,
                volumes=volumes,
                label=label,
            )
        LOGGER.info("Worklist '%s' written to %s.", worklist_name, worklist_file)
    except OSError as err:
        LOGGER.error(
            "Failed to write worklist '%s' → %s: %s", worklist_name, worklist_file, err
        )
        raise OSError(f"Failed to write worklist to {worklist_file}") from err

    # Append a brief summary for downstream display/log
    summary_title = f"To plate {destination_plate_name} add {label}"
    try:
        # assumes destination_plate.history is list of (title, data)
        summary_data = destination_plate.history[-1][1]
    except Exception:
        # fall back gracefully if history is not populated yet
        summary_data = getattr(destination_plate, "volumes", None)
    list_liquidhandling_data.append((summary_title, summary_data))

    # Optional plotting
    if plot and plot_function is not None and hasattr(destination_plate, "volumes"):
        plot_function(destination_plate.volumes, summary_title)


def process_worklist_sample_distribution(
    worklist_file: str,
    list_items: Sequence[str],
    df_worklist: pd.DataFrame,
    source_plate: Any,
    destination_plate: Any,
    group_by_column: str,
    well_column: str,
    volume_column: str,
    label: str,
    list_liquidhandling_data: list[tuple[str, Any]],
    source_column_prefix: str = "column",
    diti_reuse: int = 12,
    multi_disp: int = 12,
    plot: bool = True,
    plot_function: Optional[Callable[[Any, str], Any]] = None,
) -> None:
    """Distribute grouped items (e.g., mastermixes) to a destination plate.

    For each item in `list_items`, this function selects the corresponding rows
    from `df_worklist`, verifies a unique volume per item, and writes distribute
    commands to a Fluent `.gwl` via `robotools`.

    Parameters
    ----------
    worklist_file : str
        Path/name of the `.gwl` file to create (e.g., ``"Mastermix.gwl"``).
    list_items : sequence of str
        Items to distribute (e.g., ``["MM1", "MM2"]``).
    df_worklist : pandas.DataFrame
        Data with at least `group_by_column`, `well_column`, and `volume_column`.
    source_plate : object
        Reservoir/trough or plate providing the items (compatible with `robotools`).
    destination_plate : object
        Plate to receive the items (compatible with `robotools`).
    group_by_column : str
        Column to group on (e.g., ``"Mastermix"``).
    well_column : str
        Column name holding destination wells (e.g., ``"Well Nr. destination plate"``).
    volume_column : str
        Column name holding volumes to distribute (µL).
    label : str
        Label written into the worklist file for these operations.
    list_liquidhandling_data : list[tuple[str, Any]]
        Log/list to append a summary `(title, data)` after execution.
    source_column_prefix : str, optional
        Text prefix for the (logical) source columns. Only used for messages.
    diti_reuse : int, optional
        How often a DITI can be reused. Default ``12``.
    multi_disp : int, optional
        Number of multi-dispenses per aspiration. Default ``12``.
    plot : bool, optional
        If ``True``, create a plot of destination volumes. Default ``True``.
    plot_function : callable or None, optional
        Function like ``fn_plot(volumes, title)``. If ``None``, no plot is created.

    Returns
    -------
    None

    Raises
    ------
    KeyError
        If required columns are missing from `df_worklist`.
    ValueError
        If an item is missing in `df_worklist` or its volume is not uniquely defined.
    OSError
        If the worklist file cannot be written.
    """
    required = {group_by_column, well_column, volume_column}
    missing = [c for c in required if c not in df_worklist.columns]
    if missing:
        raise KeyError(f"Missing required column(s): {missing}")

    groups = df_worklist.groupby(group_by_column, dropna=False)

    os.makedirs(os.path.dirname(worklist_file) or ".", exist_ok=True)
    try:
        with robotools.FluentWorklist(worklist_file) as wl:
            for idx, item in enumerate(list_items):
                if item not in groups.groups:
                    raise ValueError(
                        f"Item {item!r} not found in '{group_by_column}'. "
                        f"Available: {sorted(map(str, groups.groups.keys()))}"
                    )

                sub = groups.get_group(item)

                dest_wells = sub[well_column].dropna().astype(str).tolist()
                if not dest_wells:
                    raise ValueError(f"No destination wells for item {item!r}.")

                uniq_vols = pd.Series(sub[volume_column]).dropna().unique().tolist()
                if len(uniq_vols) != 1:
                    raise ValueError(
                        f"Expected exactly one unique volume for item {item!r} in "
                        f"column '{volume_column}', got {uniq_vols}."
                    )
                volume = float(uniq_vols[0])

                # Source column: use integer index for robotools; prefix used in logs only
                LOGGER.info(
                    "Distribute %s µL of %s from %s%d to %d destination wells.",
                    volume,
                    item,
                    source_column_prefix,
                    idx + 1,
                    len(dest_wells),
                )

                wl.distribute(
                    label=label,
                    source=source_plate,
                    source_column=idx,  # robotools expects zero-based integer index
                    destination=destination_plate,
                    destination_wells=dest_wells,
                    volume=volume,
                    diti_reuse=diti_reuse,
                    multi_disp=multi_disp,
                )

        LOGGER.info("Worklist '%s' created successfully.", worklist_file)
    except OSError as err:
        LOGGER.error("Failed to write worklist '%s': %s", worklist_file, err)
        raise OSError(f"Failed to write worklist to {worklist_file}") from err

    summary_title = f"To plate {destination_plate.name} add {label}"
    try:
        summary_data = destination_plate.history[-1][1]
    except Exception:
        summary_data = getattr(destination_plate, "volumes", None)
    list_liquidhandling_data.append((summary_title, summary_data))

    if plot and plot_function is not None and hasattr(destination_plate, "volumes"):
        plot_function(destination_plate.volumes, summary_title)


# def generate_coordinates(
#     group_data: pd.DataFrame,
#     group_column: str,
#     offset: int = 1,
#     items: Optional[Sequence[Any]] = None,
# ) -> list[int]:
#     """Generate 1-D coordinates (integer indices + offset) from grouped data.

#     For each distinct value in `group_column` (or in the explicit `items` order,
#     if provided), the function returns the **first positional index** where that
#     value occurs in `group_data[group_column]`, plus an integer `offset`.
#     This is useful when mapping grouped items (e.g., mastermixes) to
#     1-based column indices required by external tools.

#     Examples
#     --------
#     >>> df = pd.DataFrame({"Group": ["MM1", "MM1", "MM2", "MM3"]})
#     >>> generate_coordinates(df, "Group", offset=1)
#     [1, 3, 4]
#     >>> # with explicit order (even if some groups appear later in the table)
#     >>> generate_coordinates(df, "Group", offset=0, items=["MM3", "MM1"])
#     [3, 0]

#     Parameters
#     ----------
#     group_data : pandas.DataFrame
#         DataFrame containing at least the `group_column`.
#     group_column : str
#         Column used for grouping (e.g., ``"Mastermix"``).
#     offset : int, optional
#         Integer added to each positional index (e.g., ``1`` for 1-based indexing).
#         Default is ``1``.
#     items : sequence or None, optional
#         If given, only these values are considered (in the given order). Values
#         not present in the data raise a ``ValueError``. If ``None``, the unique
#         values from the data (order of first appearance) are used.

#     Returns
#     -------
#     list of int
#         Positional indices (0-based) plus `offset` for each value.

#     Raises
#     ------
#     KeyError
#         If `group_column` is missing.
#     ValueError
#         If `items` contains a value not present in the column, or if the column
#         has no valid (non-NA) values.
#     """
#     if group_column not in group_data.columns:
#         raise KeyError(f"Column {group_column!r} not found in DataFrame.")

#     # Determine the sequence of values to map
#     if items is None:
#         # pandas.unique preserves order of first appearance
#         values = pd.unique(group_data[group_column].dropna())
#         values = values.tolist()
#     else:
#         # keep given order, drop NAs explicitly
#         values = [v for v in items if not (pd.isna(v))]

#     if not values:
#         raise ValueError(f"No valid values found in column {group_column!r}.")

#     col = group_data[group_column]
#     coords: list[int] = []

#     # Compute first positional index for each value
#     for val in values:
#         mask = col == val
#         # Use numpy for fast positional lookup
#         positions = np.flatnonzero(mask.to_numpy())
#         if positions.size == 0:
#             raise ValueError(f"Value {val!r} not found in column {group_column!r}.")
#         first_pos = int(positions[0])  # 0-based position
#         coords.append(first_pos + int(offset))

#     return coords


def generate_coordinates(
    group_data: pd.DataFrame,
    group_column: str,
    offset: int = 1,
    items: Optional[Sequence[Any]] = None,
) -> list[int]:
    """Generate 1-D coordinates (integer indices + offset) from grouped data.

    For each distinct value in `group_column` (or in the explicit `items` order,
    if provided), the function returns the **first positional index** where that
    value occurs in `group_data[group_column]`, plus an integer `offset`.
    This is useful when mapping grouped items (e.g., mastermixes) to
    1-based column indices required by external tools.

    Examples
    --------
    >>> df = pd.DataFrame({"Group": ["MM1", "MM1", "MM2", "MM3"]})
    >>> generate_coordinates(df, "Group", offset=1)
    [1, 3, 4]
    >>> # with explicit order (even if some groups appear later in the table)
    >>> generate_coordinates(df, "Group", offset=0, items=["MM3", "MM1"])
    [3, 0]

    Parameters
    ----------
    group_data : pandas.DataFrame
        DataFrame containing at least the `group_column`.
    group_column : str
        Column used for grouping (e.g., ``"Mastermix"``).
    offset : int, optional
        Integer added to each positional index (e.g., ``1`` for 1-based indexing).
        Default is ``1``.
    items : sequence or None, optional
        If given, only these values are considered (in the given order). Values
        not present in the data raise a ``ValueError``. If ``None``, the unique
        values from the data (order of first appearance) are used.

    Returns
    -------
    list of int
        Positional indices (0-based) plus `offset` for each value.

    Raises
    ------
    KeyError
        If `group_column` is missing.
    ValueError
        If `items` contains a value not present in the column, or if the column
        has no valid (non-NA) values.
    """
    # --- logging setup --------------------------------------------------------
    # Use module-level logger; does not change behavior if logging is disabled.
    _logger = logging.getLogger(__name__)
    _logger.debug(
        "generate_coordinates: start | column=%r, offset=%r, items_provided=%s, df_shape=%s",
        group_column,
        offset,
        items is not None,
        getattr(group_data, "shape", None),
    )
    # --------------------------------------------------------------------------

    if group_column not in group_data.columns:
        _logger.error(
            "Column %r not found in DataFrame columns: %s",
            group_column,
            list(group_data.columns),
        )
        raise KeyError(f"Column {group_column!r} not found in DataFrame.")

    # Determine the sequence of values to map
    if items is None:
        # pandas.unique preserves order of first appearance
        values = pd.unique(group_data[group_column].dropna())
        values = values.tolist()
        _logger.debug(
            "Values inferred from data (order of first appearance): %s", values
        )
    else:
        # keep given order, drop NAs explicitly
        values = [v for v in items if not (pd.isna(v))]
        _logger.debug("Values provided explicitly (after dropping NA): %s", values)

    if not values:
        _logger.error("No valid values found in column %r.", group_column)
        raise ValueError(f"No valid values found in column {group_column!r}.")

    col = group_data[group_column]
    coords: list[int] = []

    # Compute first positional index for each value
    for val in values:
        mask = col == val
        # Use numpy for fast positional lookup
        positions = np.flatnonzero(mask.to_numpy())
        if positions.size == 0:
            _logger.error("Value %r not found in column %r.", val, group_column)
            raise ValueError(f"Value {val!r} not found in column {group_column!r}.")
        first_pos = int(positions[0])  # 0-based position
        coord = first_pos + int(offset)
        coords.append(coord)
        _logger.debug(
            "Mapped value=%r -> first_pos=%d (0-based) -> coord=%d (offset=%d)",
            val,
            first_pos,
            coord,
            offset,
        )

    _logger.info(
        "generate_coordinates: done | n_values=%d, offset=%d, coords=%s",
        len(values),
        offset,
        coords,
    )
    return coords


def modify_gwl_file(
    input_file: str | os.PathLike[str],
    output_file: str | os.PathLike[str],
    list_coordinates: Sequence[int | str],
    source_name: str = "mtp_source[001]",
    line_prefix: str = "R;",
) -> None:
    """Rewrite selected lines in a Tecan GWL file.

    Lines beginning with `line_prefix` are updated: the source name (field 1)
    is set to `source_name`, and both coordinate fields (indices 4 and 5) are
    replaced by values from `list_coordinates` in order.

    Parameters
    ----------
    input_file : str or os.PathLike
        Path to the input ``.gwl`` file.
    output_file : str or os.PathLike
        Path to the output ``.gwl`` file to write.
    list_coordinates : sequence of int or str
        Coordinates to apply to matching lines (must be at least as many as
        there are matching lines).
    source_name : str, optional
        Replacement for the source field (default ``"mtp_source[001]"``).
    line_prefix : str, optional
        Line prefix that marks rows to modify (default ``"R;"``).

    Returns
    -------
    None

    Raises
    ------
    FileNotFoundError
        If `input_file` does not exist.
    ValueError
        If there are more matching lines than provided coordinates, or if
        a matching line is malformed (fewer than 6 ``;``-separated fields).
    OSError
        If `output_file` cannot be written.
    """
    in_path = os.fspath(input_file)
    out_path = os.fspath(output_file)

    # Read all lines (preserve newlines)
    try:
        with open(in_path, encoding="utf-8", newline="") as handle:
            lines = handle.readlines()
    except FileNotFoundError as err:
        raise FileNotFoundError(f"Input GWL not found: {in_path}") from err
    except OSError as err:
        raise OSError(f"Failed to read GWL file: {in_path}") from err

    # Count target lines first to validate coordinate availability
    match_indices = [i for i, ln in enumerate(lines) if ln.startswith(line_prefix)]
    if len(list_coordinates) < len(match_indices):
        raise ValueError(
            "Not enough coordinates provided: "
            f"{len(list_coordinates)} < {len(match_indices)} required."
        )

    new_lines: list[str] = []
    coord_idx = 0

    for line in lines:
        if line.startswith(line_prefix):
            parts = line.rstrip("\n").split(";")
            if len(parts) < 6:
                # Malformed R; line: fail fast with context
                raise ValueError(
                    f"Malformed GWL line (expected ≥6 fields): {line.strip()!r}"
                )
            parts[1] = source_name
            coord = str(list_coordinates[coord_idx])
            parts[4] = coord
            parts[5] = coord
            new_lines.append(";".join(parts) + "\n")
            coord_idx += 1
        else:
            new_lines.append(line)

    # Ensure output directory exists
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    try:
        with open(out_path, "w", encoding="utf-8", newline="") as handle:
            handle.writelines(new_lines)
    except OSError as err:
        raise OSError(f"Failed to write modified GWL to: {out_path}") from err

    extra = len(list_coordinates) - coord_idx
    if extra > 0:
        LOGGER.info(
            "GWL modified (%d lines). %d coordinate(s) unused.", coord_idx, extra
        )
    else:
        LOGGER.info("GWL modified (%d lines).", coord_idx)


def process_plan(
    mtp_manager: Any,
    worklist_data: dict[str, pd.DataFrame],
    workflow: Any,  # or: "LiquidHandlingWorkflow"
    stocks_information: pd.DataFrame,
    mastermix_information: pd.DataFrame,
    water_information: pd.DataFrame,
    sourceplate_information: pd.DataFrame,
    list_liquidhandling_data: list[tuple[str, Any]],
    plot: bool = True,
    plot_function: Optional[Callable[[Any, str], Any]] = None,
) -> None:
    """Execute all workflow steps: sample distribution or transfer per step.

    For each step returned by ``workflow.list_steps()``, this function dispatches
    to either :func:`process_worklist_sample_distribution` or
    :func:`process_worklist`, rewrites the generated ``.gwl`` (via
    :func:`modify_gwl_file`) with the correct source coordinates, and appends a
    short summary to ``list_liquidhandling_data``. Optional plotting is supported
    through ``plot_function(volumes, title)``.

    Parameters
    ----------
    mtp_manager : object
        Manager providing ``get_mtp(name)`` for plate objects.
    worklist_data : dict[str, pandas.DataFrame]
        Mapping of worklist-name → table used for the step.
    workflow : LiquidHandlingWorkflow-like
        Object exposing ``list_steps() -> Iterable[dict]``; each dict must at
        least have: ``process_label``, ``sample_handling``, ``liquid_handling``,
        ``source_plate`` (or ``fake_source_plate``), ``destination_plate``.
    stocks_information : pandas.DataFrame
        Table with stock/source information (must include a grouping column used
        by :func:`generate_coordinates`, typically ``"Content"``).
    mastermix_information : pandas.DataFrame
        Source table used when step label indicates a mastermix.
    water_information : pandas.DataFrame
        Source table used when step label indicates water.
    sourceplate_information : pandas.DataFrame
        Fallback source table if the chosen information frame is empty.
    list_liquidhandling_data : list[tuple[str, Any]]
        List that will be appended with (title, data) tuples for downstream
        display/logging.
    plot : bool, optional
        If ``True``, plots are produced via ``plot_function`` when provided.
    plot_function : callable or None, optional
        Function like ``fn_plot(volumes, title)`` injected by the caller.
        If ``None``, no plot is generated.

    Returns
    -------
    None

    Raises
    ------
    KeyError
        If required keys are missing in a step or a needed worklist is absent.
    ValueError
        If mandatory columns/plates are missing or coordinate derivation fails.
    """
    # defensive check: workflow must offer list_steps
    if not hasattr(workflow, "list_steps"):
        raise ValueError("workflow must provide a 'list_steps()' method.")

    # Helper: infer grouping column for sample distributions
    candidate_columns = [
        "Mastermix",
        "Primer paar",
        "Water",
        "Backbone",
        "Promotor (A-B)",
        "RBS",
        "Signal Peptide",
        "GoI",
        "Tag/Term",
        "Plasmidmix",
        "Part",
        "Part_1",
        "Part_2",
        "Part_3",
        "Part_4",
        "Part_5",
        "Part_6",
        "Mix_1",
        "Mix_2",
    ]

    for step in workflow.list_steps():
        label = step.get("process_label", "<unnamed>")
        LOGGER.info("Processing step: %s", label)

        handling = step.get("sample_handling", "").lower()
        if handling == "sample_distribution":
            LOGGER.info("Handling sample distribution.")

            try:
                worklist_df = worklist_data.get(step["process_label"])
                if worklist_df is None:
                    raise KeyError(
                        f"Worklist {step['process_label']!r} not found in worklist_data."
                    )

                # Volume validation (WARN/ERROR go to warn/error log if handlers are attached)
                id_cols: list[str] = [
                    c
                    for c in [
                        "Mastermix",
                        "Plasmidmix",
                        "Part",
                        "Well Nr. destination plate",
                        "Plate",
                    ]
                    if c in worklist_df.columns
                ]
                check_and_log_volumes(
                    df=worklist_df,
                    volume_column="Volumen",
                    min_volume_ul=2.0,  # adjust threshold as needed
                    group_columns=id_cols,
                    context=f"Step={label}",
                )

                # Determine grouping column
                group_column: Optional[str] = next(
                    (c for c in candidate_columns if c in worklist_df.columns), None
                )
                if not group_column:
                    raise ValueError(
                        f"Could not determine group column for worklist {step['process_label']!r}."
                    )

                # Choose source/destination plates (fake overrides real if provided)
                fake_src = step.get("fake_source_plate")
                source_plate = (
                    mtp_manager.get_mtp(fake_src)
                    if fake_src
                    else mtp_manager.get_mtp(step["source_plate"])
                )
                destination_plate = mtp_manager.get_mtp(step["destination_plate"])
                if source_plate is None or destination_plate is None:
                    raise ValueError(
                        "Source or destination plate not found in MTPManager."
                    )

                # Write distribute worklist (keeps Excel order as-is)
                process_worklist_sample_distribution(
                    worklist_file=f"{step['process_label']}.gwl",
                    list_items=worklist_df[group_column].dropna().unique().tolist(),
                    df_worklist=worklist_df,
                    source_plate=source_plate,
                    destination_plate=destination_plate,
                    group_by_column=group_column,
                    well_column="Well Nr. destination plate",
                    volume_column="Volumen",
                    label=step["liquid_handling"],
                    list_liquidhandling_data=list_liquidhandling_data,
                    plot=plot,
                    plot_function=plot_function,
                )

                # Select correct source-information table
                plabel_lower = str(step.get("process_label", "")).lower()
                if "mastermix" in plabel_lower:
                    source_information = mastermix_information
                elif "water" in plabel_lower:
                    source_information = water_information
                else:
                    source_information = stocks_information

                if source_information is None or source_information.empty:
                    source_information = sourceplate_information

                # Determine coordinates from source_information (kept for logging/compat)
                group_col_for_coords = (
                    "Content"
                    if "Content" in source_information.columns
                    else group_column
                )
                if group_col_for_coords not in source_information.columns:
                    raise ValueError(
                        "Cannot determine coordinate grouping column: "
                        f"expected 'Content' or '{group_column}'."
                    )

                coords_global = generate_coordinates(
                    group_data=source_information,
                    group_column=group_col_for_coords,
                    offset=1,
                )
                LOGGER.info(
                    "Global coordinate mapping (may contain gaps): %s", coords_global
                )

                # Compact, step-local coordinates (no gaps), Excel order
                list_items_order = (
                    worklist_df[group_column]
                    .dropna()
                    .astype(str)
                    .drop_duplicates(keep="first")
                    .tolist()
                )
                coordinates = list(range(1, len(list_items_order) + 1))
                LOGGER.info(
                    "Assigned compact source coordinates (per-step): items=%s -> coords=%s",
                    list_items_order,
                    coordinates,
                )

                # Rewrite the GWL with proper source and coordinates (in place)
                modify_gwl_file(
                    input_file=f"{step['process_label']}.gwl",
                    output_file=f"{step['process_label']}.gwl",
                    list_coordinates=coordinates,  # compact per-step mapping
                    source_name=step["source_plate"],
                )

            except (KeyError, ValueError) as err:
                LOGGER.error(
                    "Error in sample distribution for step '%s': %s", label, err
                )
                # continue with next step
            except Exception as err:  # pylint: disable=broad-except
                LOGGER.exception(
                    "Unexpected error in sample distribution for step '%s': %s",
                    label,
                    err,
                )
                # continue with next step

        elif handling == "transfer":
            LOGGER.info("Handling transfer.")

            try:
                process_worklist(
                    mtp_manager=mtp_manager,
                    worklist_data=worklist_data,
                    worklist_name=step["process_label"],
                    source_plate_name=step["source_plate"],
                    destination_plate_name=step["destination_plate"],
                    worklist_file=f"{step['process_label']}.gwl",
                    label=step["liquid_handling"],
                    list_liquidhandling_data=list_liquidhandling_data,
                    plot=plot,
                    plot_function=plot_function,
                )
            except (KeyError, ValueError) as err:
                LOGGER.error("Error in transfer for step '%s': %s", label, err)
            except Exception as err:  # pylint: disable=broad-except
                LOGGER.exception(
                    "Unexpected error in transfer for step '%s': %s", label, err
                )

        else:
            LOGGER.error(
                "Unknown sample handling type: %r", step.get("sample_handling")
            )


def read_and_process_dataframe(
    file_path: str,
    sheet_name: str,
    start_column: str,
    well_column: str,
    empty_fill: str = " ",
    required_columns: Iterable[str] | None = None,
) -> pd.DataFrame:
    """Read a range of columns from an Excel sheet and return a cleaned table.

    The function reads the sheet with header detection, slices column-wise
    from the column *letter* `start_column` up to and including `well_column`,
    removes fully-empty columns, fills internal NaNs with a placeholder,
    and drops rows that are empty in all **non-well** columns.

    This is useful when your design/worklist tables have varying leading
    columns but always contain a trailing ``Well`` column.

    Parameters
    ----------
    file_path : str
        Path to the Excel file.
    sheet_name : str
        Sheet name to read.
    start_column : str
        Excel-style column letters where the relevant table starts, e.g. ``\"C\"`` or ``\"AA\"``.
    well_column : str
        Name of the column that contains well IDs (e.g., ``\"Well\"``).
        This column must exist in the sheet.
    empty_fill : str, optional
        Placeholder used to fill NaNs in non-well columns (default: single space).
    required_columns : Iterable[str] or None, optional
        Additional column names that must be present in the sliced DataFrame.
        If any are missing, a ``ValueError`` is raised.

    Returns
    -------
    pd.DataFrame
        Cleaned DataFrame from `start_column` to `well_column` (inclusive).

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the sheet cannot be read, `well_column` is missing,
        `start_column` resolves to an index beyond the sheet width,
        or required columns are missing.
    """
    try:
        df_full = pd.read_excel(file_path, sheet_name=sheet_name, header=0)
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"File not found: {file_path}") from exc
    except Exception as exc:  # pylint: disable=broad-except
        raise ValueError(f"Failed to read sheet {sheet_name!r}: {exc}") from exc

    if well_column not in df_full.columns:
        raise ValueError(f"The well column {well_column!r} is missing in the sheet.")

    start_idx = col_letter_to_index(start_column)
    if start_idx >= df_full.shape[1]:
        raise ValueError(
            f"start_column {start_column!r} resolves to index {start_idx}, "
            f"but the sheet has only {df_full.shape[1]} columns."
        )

    # Slice from letter-based start to the column named `well_column` (inclusive)
    well_idx = df_full.columns.get_loc(well_column)
    if start_idx > well_idx:
        raise ValueError(
            "start_column is to the right of the 'well_column'. Check your inputs."
        )

    df = df_full.iloc[:, start_idx : well_idx + 1].copy()

    # Drop columns that are entirely empty
    df = df.dropna(how="all", axis=1)

    # Fill NaNs in non-well columns so that row-drop below keeps intended rows
    non_well_cols = [c for c in df.columns if c != well_column]
    if non_well_cols:
        df[non_well_cols] = df[non_well_cols].fillna(empty_fill)

    # Remove rows that are empty in all non-well columns
    if non_well_cols:
        mask_all_empty = (df[non_well_cols] == empty_fill).all(axis=1)
        df = df.loc[~mask_all_empty].copy()

    # Validate required columns if any
    if required_columns:
        missing = [c for c in required_columns if c not in df.columns]
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

    LOGGER.debug(
        "Processed table from %s[%s]: shape=%s (start=%s, well=%s)",
        file_path,
        sheet_name,
        tuple(df.shape),
        start_column,
        well_column,
    )
    return df.reset_index(drop=True)


class AssemblyWorkflowManager:
    """Coordinate loading, planning, worklist generation and simulation for assemblies.

    Notes
    -----
    This manager ties together:
    - metadata & plan loading from an Excel workbook
    - MTP (labware) access via an injected `mtp_manager`
    - workflow steps via an injected `workflow`
    - worklist sheet discovery, `.gwl` creation, and volume simulation
    """

    # Default volume factors (kept as attributes for runtime tweaking)
    SAMPLE_DISTRIBUTION_DEAD_VOLUME_FACTOR: float = 1.3
    TRANSFER_DEAD_VOLUME_MINIMUM: float = 10.0  # µL
    TRANSFER_DEAD_VOLUME_FACTOR: float = 1.1
    DEAD_VOLUME_PER_WELL: float = 8.0  # µL

    def __init__(self) -> None:
        """Initialize empty manager state.

        Attributes set to `None` will be populated by :meth:`initialize`.
        """
        self.file_path: Optional[str] = None
        self.mtp_manager: Optional[Any] = None
        self.workflow: Optional[Any] = None

        self.metadata: Optional[pd.DataFrame] = None
        self.stocks_information: Optional[pd.DataFrame] = None
        self.mastermix_information: Optional[pd.DataFrame] = None
        self.water_information: Optional[pd.DataFrame] = None
        self.sourceplate_information: Optional[pd.DataFrame] = None
        self.df_pm: Optional[pd.DataFrame] = None
        self.df_dest_plate: Optional[pd.DataFrame] = None
        self.df_volumes: Optional[pd.DataFrame] = None
        self.df_mastermix_components: Optional[pd.DataFrame] = None

        self.worklist_data: dict[str, pd.DataFrame] = {}
        self.sheet_name: Optional[str] = None
        self.list_liquidhandling_data: list[tuple[str, Any]] = []
        self.worklist_sheets: Optional[list[str]] = None
        self.parts: Optional[pd.DataFrame] = None

        self.sum_mastermix_volume: float = 7.5
        self.volume_per_reaction: Optional[float] = None

    def initialize(
        self,
        mtp_manager: Any,
        workflow: Any,
        sheet_name: str,
        file_path: str,
        type_of_molbiowork: str,
        plot: bool = True,
        plot_information: bool = True,
    ) -> None:
        """Load metadata/plan, prepare workflow components, and set up MTP volumes.

        This:
        - loads metadata and plan sheets from `file_path`
        - wires `workflow` and `mtp_manager`
        - discovers worklist sheets and creates empty `.gwl` files
        - loads worklist DataFrames
        - simulates initial volumes and (optionally) plots MTP volumes

        Parameters
        ----------
        mtp_manager : object
            Manager for labware (must provide `get_mtp(name)` and plotting support).
        workflow : object
            Workflow object (must provide `list_steps()`).
        sheet_name : str
            Name of the Excel sheet containing the destination plan.
        file_path : str
            Path to the Excel workbook with metadata and planning information.
        type_of_molbiowork : str
            Type of workflow (e.g., ``"PCR"`` or ``"Golden Gate"``).
        plot : bool, optional
            If ``True``, visualize MTP volumes after initialization. Default ``True``.
        plot_information : bool, optional
            If ``True``, include tabular volume info in logs. Default ``True``.

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If required data cannot be loaded or parsed.
        """
        LOGGER.info("Initializing AssemblyWorkflowManager...")
        try:
            # Store wiring & basic params
            self.file_path = file_path
            self.mtp_manager = mtp_manager
            self.workflow = workflow
            self.sheet_name = sheet_name

            # --- Metadata & plan
            self.metadata = self._load_metadata()  # uses self.file_path
            self.df_dest_plate = self._load_plan()

            # --- PM info depends on workflow type
            if type_of_molbiowork == "PCR":
                self.df_pm = self._load_df(
                    usecols=[1, 3, 4], sheet_name="PM_info+"
                ).dropna()
            else:
                self.df_pm = self._load_df(
                    usecols=[1, 3, 4, 5, 6], sheet_name="PM_info+"
                ).dropna()

            # Combine destination & PM to a parts view (non-PCR)
            if type_of_molbiowork != "PCR":
                self.parts = self._combine_and_filter_parts(
                    self.df_dest_plate, self.df_pm
                )

            # --- Source information sheets
            self.stocks_information = self._load_df(usecols=[0, 1])
            self.mastermix_information = self._load_df(usecols=[4, 5]).rename(
                {"Content.1": "Content", "Well.1": "Well"}, axis=1
            )
            self.water_information = self._load_df(usecols=[12, 13]).rename(
                {"Content.3": "Content", "Well.3": "Well"}, axis=1
            )
            self.sourceplate_information = self._load_df(usecols=[8, 9]).rename(
                {"Content.2": "Content", "Well.2": "Well"}, axis=1
            )

            # --- Consistency checks & worklist loading
            self._log_mismatched_parts()
            self.worklist_sheets = self._find_worklist_sheets()
            self._create_empty_worklist_files(self.worklist_sheets)
            self.worklist_data = self._load_worklist_data(self.worklist_sheets)

            # --- Volume simulation & derived values
            self.df_volumes = self.simulate_volumes(plot_information=plot_information)
            self.sum_mastermix_volume = self._sum_mastermix_volumes()

            # Derive volume per reaction if available
            try:
                vols = (
                    pd.Series(self.worklist_data["Worklist_MM"]["Volumen"])
                    .dropna()
                    .unique()
                    .tolist()
                )
                self.volume_per_reaction = float(vols[0]) if vols else None
            except Exception as err:  # pylint: disable=broad-except
                LOGGER.warning("Could not derive 'volume_per_reaction': %s", err)
                self.volume_per_reaction = None

            # --- Optional plots
            self._plot_mtp_volumes(plot)
            LOGGER.info("Initialization complete.")

        except (KeyError, ValueError) as err:
            LOGGER.error("Initialization failed: %s", err)
            raise
        except Exception as err:  # pylint: disable=broad-except
            LOGGER.exception("Unexpected error during initialization: %s", err)
            raise ValueError(f"Initialization failed: {err}") from err

    def simulate_volumes(self, plot_information: bool = True) -> pd.DataFrame:
        """Simulate required volumes and initialize starting plates.

        Parameters
        ----------
        plot_information : bool, optional
            If ``True``, log the computed volume table. Default ``True``.

        Returns
        -------
        pandas.DataFrame
            The computed volume table. If an error occurs, an empty DataFrame
            is returned and the error is logged.
        """
        LOGGER.info("Simulating volumes...")
        try:
            volume_df = self.calculate_volumes()
            self.initialize_volumes(volume_df)
            if plot_information:
                # Compact but readable logging of the table
                LOGGER.info("Volume summary:\n%s", volume_df.to_string(index=False))
            LOGGER.info("Simulation complete.")
            return volume_df
        except Exception as err:  # pylint: disable=broad-except
            LOGGER.exception("Error during simulation: %s", err)
            return pd.DataFrame()

    def _get_step_factor(self, step: dict[str, Any]) -> float:
        """Return per-step additional liquid factor (defaults to 1.0)."""
        try:
            return float(step.get("additional_liquid_factor", 1.0))
        except Exception:  # pylint: disable=broad-except
            return 1.0

    def calculate_volumes(self) -> pd.DataFrame:
        """Calculate required volumes for starting and source plates based on the workflow.

        The procedure:
        1) Process the ``Worklist_DNA`` (if present and `start == "Yes"`) to mark parts
        as already handled.
        2) Process remaining worklists (with `start == "Yes"`) and accumulate volumes
        per component for the corresponding source plates. Each row volume is scaled
        by the per-step ``additional_liquid_factor``.
        3) Apply dead-volume adjustments and produce a summary table. If a given
        (plate, component) has contributions from both handling types, apply the
        corresponding dead-volume factors **separately** and sum the results.
        Finally, enforce a minimum total volume.

        Returns
        -------
        pandas.DataFrame
            A table with columns: ``Plate``, ``Component``, ``Total Volume (µL)``, ``Wells``.

        Notes
        -----
        - Expected columns in relevant worklists: ``"Part"``, ``"Volumen"``,
        ``"Well Nr. source plate"``.
        - If these columns are missing for a step, that step is skipped.
        """
        if self.workflow is None:
            raise ValueError("Workflow is not initialized.")
        if not isinstance(self.worklist_data, dict):
            raise ValueError("worklist_data must be a dictionary of DataFrames.")

        processed_parts: set[str] = set()
        # Separate buckets per handling type:
        volume_summary_sd: dict[str, dict[str, float]] = defaultdict(
            lambda: defaultdict(float)
        )
        volume_summary_tr: dict[str, dict[str, float]] = defaultdict(
            lambda: defaultdict(float)
        )
        wells_summary: dict[str, dict[str, set[str]]] = defaultdict(
            lambda: defaultdict(set)
        )

        # --- 1) Handle Worklist_DNA first ---------------------------------------
        for step in self.workflow.list_steps():
            process_label = step.get("process_label")
            if not process_label:
                continue
            if (
                process_label == "Worklist_DNA"
                and step.get("start", "No").lower() == "yes"
            ):
                df = self.worklist_data.get(process_label)
                if df is None or df.empty:
                    continue
                required = {"Part", "Volumen", "Well Nr. source plate"}
                if not required.issubset(df.columns):
                    LOGGER.warning(
                        "Skipping Worklist_DNA; missing columns: %s",
                        sorted(required - set(df.columns)),
                    )
                    continue

                handling_type = (
                    str(step.get("sample_handling", "transfer")).lower() or "transfer"
                )
                try:
                    step_factor = float(step.get("additional_liquid_factor", 1.0))
                except Exception:  # pylint: disable=broad-except
                    step_factor = 1.0

                src_plate = step.get("source_plate")
                if not src_plate:
                    continue

                for _, row in df.iterrows():
                    component = row.get("Part")
                    volume = row.get("Volumen", 0)
                    well = row.get("Well Nr. source plate")
                    if not component or pd.isna(component) or pd.isna(well):
                        continue

                    comp_str = str(component)
                    vol_scaled = float(volume) * step_factor
                    wells_summary[src_plate][comp_str].add(str(well))

                    if handling_type == "sample_distribution":
                        volume_summary_sd[src_plate][comp_str] += vol_scaled
                    else:
                        volume_summary_tr[src_plate][comp_str] += vol_scaled

                    processed_parts.add(comp_str)

                    LOGGER.debug(
                        "Accum(DNA): plate=%s part=%s step=%s handling=%s vol=%.3f * factor=%.3f -> %.3f",
                        src_plate,
                        comp_str,
                        process_label,
                        handling_type,
                        float(volume),
                        step_factor,
                        vol_scaled,
                    )

        # --- 2) Handle all other start steps (excluding Worklist_DNA) -----------
        for step in self.workflow.list_steps():
            process_label = step.get("process_label")
            if not process_label or process_label == "Worklist_DNA":
                continue
            if step.get("start", "No").lower() != "yes":
                continue

            df = self.worklist_data.get(process_label)
            if df is None or df.empty:
                continue

            required = {"Part", "Volumen", "Well Nr. source plate"}
            if not required.issubset(df.columns):
                LOGGER.warning(
                    "Skipping '%s'; missing columns: %s",
                    process_label,
                    sorted(required - set(df.columns)),
                )
                continue

            handling_type = (
                str(step.get("sample_handling", "transfer")).lower() or "transfer"
            )
            try:
                step_factor = float(step.get("additional_liquid_factor", 1.0))
            except Exception:  # pylint: disable=broad-except
                step_factor = 1.0

            src_plate = step.get("source_plate")
            if not src_plate:
                continue

            for _, row in df.iterrows():
                component = row.get("Part")
                volume = row.get("Volumen", 0)
                well = row.get("Well Nr. source plate")
                if not component or pd.isna(component):
                    continue
                comp_str = str(component)
                if comp_str in processed_parts:
                    continue

                if not pd.isna(well):
                    wells_summary[src_plate][comp_str].add(str(well))

                vol_scaled = float(volume) * step_factor
                if handling_type == "sample_distribution":
                    volume_summary_sd[src_plate][comp_str] += vol_scaled
                else:
                    volume_summary_tr[src_plate][comp_str] += vol_scaled

                LOGGER.debug(
                    "Accum: plate=%s part=%s step=%s handling=%s vol=%.3f * factor=%.3f -> %.3f",
                    src_plate,
                    comp_str,
                    process_label,
                    handling_type,
                    float(volume),
                    step_factor,
                    vol_scaled,
                )

        # --- 3) Apply dead-volume logic and build final table --------------------
        final_rows: list[dict[str, object]] = []
        all_plates = set(volume_summary_sd.keys()) | set(volume_summary_tr.keys())

        for plate in sorted(all_plates):
            components = set(volume_summary_sd[plate].keys()) | set(
                volume_summary_tr[plate].keys()
            )
            for component in sorted(components):
                wells = sorted(wells_summary[plate][component])
                num_wells = len(wells) if wells else 1

                total_sd = float(volume_summary_sd[plate].get(component, 0.0))
                total_tr = float(volume_summary_tr[plate].get(component, 0.0))

                # apply handling-specific dead-volume scaling separately
                adj_sd = (
                    total_sd * num_wells * self.SAMPLE_DISTRIBUTION_DEAD_VOLUME_FACTOR
                )
                adj_tr = total_tr * num_wells * self.TRANSFER_DEAD_VOLUME_FACTOR

                # keep existing model: one per-plate additive dead volume
                adjusted_volume = adj_sd + adj_tr + self.DEAD_VOLUME_PER_WELL

                # enforce minimum total
                if adjusted_volume < self.TRANSFER_DEAD_VOLUME_MINIMUM:
                    LOGGER.debug(
                        "Enforcing minimum dead volume for %s/%s: %.2f -> %.2f µL",
                        plate,
                        component,
                        adjusted_volume,
                        self.TRANSFER_DEAD_VOLUME_MINIMUM,
                    )
                    adjusted_volume = self.TRANSFER_DEAD_VOLUME_MINIMUM

                final_rows.append(
                    {
                        "Plate": plate,
                        "Component": component,
                        "Total Volume (µL)": round(adjusted_volume, 2),
                        "Wells": ", ".join(wells),
                    }
                )

        # Optional: detailed per-component logging (debug level)
        for row in final_rows:
            LOGGER.debug(
                "Plate=%s Component=%s Total=%.2fµL Wells=%s",
                row["Plate"],
                row["Component"],
                row["Total Volume (µL)"],
                row["Wells"],
            )

        return pd.DataFrame(final_rows)

    def initialize_volumes(self, volume_df: pd.DataFrame) -> None:
        """Initialize source plates with calculated volumes.

        For each plate in `volume_df`, distribute the listed total volume evenly
        across its wells and add these volumes to the corresponding MTP object.

        Parameters
        ----------
        volume_df : pandas.DataFrame
            Must contain columns: ``"Plate"``, ``"Wells"``, ``"Total Volume (µL)"``.

        Returns
        -------
        None
        """
        required_cols = {"Plate", "Wells", "Total Volume (µL)"}
        if not required_cols.issubset(volume_df.columns):
            missing = required_cols - set(volume_df.columns)
            raise ValueError(
                f"volume_df is missing required columns: {sorted(missing)}"
            )

        if self.mtp_manager is None:
            raise ValueError("MTP manager is not initialized.")

        for plate_name, plate_df in volume_df.groupby("Plate"):
            plate = self.mtp_manager.get_mtp(plate_name)
            if plate is None:
                LOGGER.warning(
                    "Plate '%s' not found in MTPManager; skipping.", plate_name
                )
                continue

            wells: list[str] = []
            volumes: list[float] = []

            for _, row in plate_df.iterrows():
                wells_str = str(row["Wells"]).strip()
                if not wells_str or wells_str.lower() in {"nan", "none"}:
                    continue
                component_wells = [w.strip() for w in wells_str.split(",") if w.strip()]
                if not component_wells:
                    continue
                total_vol = float(row["Total Volume (µL)"])
                per_well = total_vol / len(component_wells)
                wells.extend(component_wells)
                volumes.extend([per_well] * len(component_wells))

            if wells:
                LOGGER.info(
                    "Initializing %.1f µL in %s for %d wells.",
                    sum(volumes),
                    plate_name,
                    len(wells),
                )
                # Assumes your plate object supports .add(wells, volumes)
                plate.add(wells, volumes)
            else:
                LOGGER.info("No wells to initialize for plate '%s'.", plate_name)

    def identify_start_plates(self) -> set[str]:
        """Return names of source plates whose steps are marked as start=='Yes'."""
        if self.workflow is None:
            raise ValueError("Workflow is not initialized.")
        start_plates: set[str] = set()
        for step in self.workflow.list_steps():
            if step.get("start", "No").lower() == "yes":
                sp = step.get("source_plate")
                if isinstance(sp, str):
                    start_plates.add(sp)
        start_plates.discard("water")
        return start_plates

    def process_workflow(
        self,
        plot: bool = True,
        plot_function: Optional[Callable[[Any, str], Any]] = None,
        *,
        export_worklists: bool = True,
        export_dir: str = "worklist",
        overwrite_worklists: bool = True,
    ) -> None:
        """Process all workflow steps and optionally export worklists as CSV files.

        This method executes all liquid-handling steps defined in the workflow.
        Each step triggers the corresponding processing function (e.g. transfer,
        distribution, mixing) and optionally generates plots.

        If enabled, all generated worklists (stored in ``self.worklist_data``)
        are exported as individual CSV files. Additionally, a combined
        ``Master_Worklist.csv`` is created for use with Opentrons. This master
        file contains all worklists merged into a single table and includes a
        ``Source_Worklist`` column to trace the origin of each row.

        Well identifiers are normalized to Opentrons-compatible format
        (e.g. ``A01`` → ``A1``), and CSV files are written using comma
        separation.

        Parameters
        ----------
        plot : bool, optional
            If ``True``, plotting is enabled for workflow steps that support it.
            Default is ``True``.
        plot_function : Callable or None, optional
            Optional plotting function (e.g. ``fn_plot(volumes, title)``) that is
            passed to workflow steps. If ``None``, no plots are generated.
            Default is ``None``.
        export_worklists : bool, optional
            If ``True``, export each worklist DataFrame in ``self.worklist_data``
            as an individual CSV file and generate a combined
            ``Master_Worklist.csv`` for Opentrons.
            Default is ``True``.
        export_dir : str, optional
            Output directory where CSV files are written. The directory is
            created if it does not exist.
            Default is ``"worklist"``.
        overwrite_worklists : bool, optional
            If ``True``, existing CSV files are overwritten. If ``False``,
            numeric suffixes are appended to avoid overwriting existing files.
            Default is ``True``.

        Returns
        -------
        None
            This method processes the workflow in-place and writes output files
            if export is enabled.

        Notes
        -----
        - Individual worklists are preserved as separate CSV files.
        - The combined ``Master_Worklist.csv`` is required for Opentrons and
        contains all worklists in a single table.
        - The ``Source_Worklist`` column allows reconstruction of the original
        worklists from the master file.
        - CSV files use comma (`,`) as separator for compatibility with
        Opentrons.
        - Invalid or non-standard well identifiers are left unchanged.

        Raises
        ------
        ValueError
            If workflow steps are malformed or required data is missing.
        RuntimeError
            If execution of a workflow step fails.
        """
        if any(
            x is None for x in (self.mtp_manager, self.workflow, self.worklist_data)
        ):
            raise ValueError("Manager not initialized. Call 'initialize' first.")

        LOGGER.info("Processing workflow...")

        # Optional export of the loaded worklist DataFrames as Excel files.
        if export_worklists:
            try:
                self.export_worklists_to_csvs(
                    directory=export_dir,
                    include_empty=False,
                    overwrite=overwrite_worklists,
                    export_master=True,
                    master_filename="Master_Worklist.csv",
                )
            except Exception as err:  # pylint: disable=broad-except
                LOGGER.error("Exporting worklists failed: %s", err)

        process_plan(
            mtp_manager=self.mtp_manager,
            worklist_data=self.worklist_data,
            workflow=self.workflow,
            stocks_information=self.stocks_information,
            mastermix_information=self.mastermix_information,
            water_information=self.water_information,
            sourceplate_information=self.sourceplate_information,
            plot=plot,
            list_liquidhandling_data=self.list_liquidhandling_data,
            plot_function=plot_function,
        )

    def plot_destination_plate(
        self,
        save_path: str = "PCR_Plan.png",
        highlight_column: str = "Mastermix",
        font_size: int = 8,
        wrap_width: int = 15,
        cmap: str = "tab10",
        max_lines: int = 5,
        exclude_keywords: Optional[Sequence[str]] = None,
        color_palette: Optional[Sequence[Any]] = None,
    ) -> None:
        """Render the destination plate layout as an image."""
        if self.df_dest_plate is None:
            raise ValueError("No destination plate data loaded.")
        title = "Destination"
        try:
            mtp = (
                self.mtp_manager.get_mtp("mtp_destination")
                if self.mtp_manager
                else None
            )
            if mtp and getattr(mtp, "name", None):
                title = mtp.name
        except Exception as err:  # pylint: disable=broad-except
            LOGGER.warning("Failed to read destination plate name: %s", err)

        plot_mikrotiterplatte(
            df=self.df_dest_plate,
            title=title,
            highlight_column=highlight_column,
            font_size=font_size,
            wrap_width=wrap_width,
            save_path=save_path,
            max_lines=max_lines,
            exclude_keywords=exclude_keywords,
            color_palette=color_palette,
            cmap=cmap,
        )
        LOGGER.info("Destination plate plot saved to %s", save_path)

    def simulate_liquid_handling(
        self, fps: float = 0.5, delay_frames: int = 1
    ) -> Image:
        """Simulate liquid handling workflow and return an animated GIF as Image."""
        out_path = Path("results") / "Workflow_simulation.gif"

        fp_gif = plot_gif(
            fn_plot=fn_plot2,  # injected util in your package
            fp_out=out_path,
            data=self.list_liquidhandling_data,
            fps=fps,
            delay_frames=delay_frames,
        )

        LOGGER.info("Workflow simulation GIF written to %s", fp_gif)
        return Image(filename=str(fp_gif))

    def export_worklists_to_csvs(
        self,
        directory: str = "worklist",
        *,
        include_empty: bool = False,
        overwrite: bool = True,
        export_master: bool = True,
        master_filename: str = "Master_Worklist.csv",
    ) -> Path:
        """Export worklists to individual CSV files and one combined master CSV.

        Parameters
        ----------
        directory : str, optional
            Target directory for exported CSV files. Default is ``"worklist"``.
        include_empty : bool, optional
            If ``True``, empty DataFrames are also exported. Default ``False``.
        overwrite : bool, optional
            If ``True``, existing files are overwritten. Otherwise numeric suffixes
            are appended. Default ``True``.
        export_master : bool, optional
            If ``True``, also export a combined master CSV. Default ``True``.
        master_filename : str, optional
            File name for the combined master CSV. Default is
            ``"Master_Worklist.csv"``.

        Returns
        -------
        pathlib.Path
            Directory where the files were written.

        Raises
        ------
        ValueError
            If no worklist data are available.
        """
        if not isinstance(self.worklist_data, dict) or not self.worklist_data:
            raise ValueError("No worklist data available to export.")

        out_dir = Path(directory)
        out_dir.mkdir(parents=True, exist_ok=True)

        def _sanitize(name: str) -> str:
            safe = "".join(
                ch if ch.isalnum() or ch in {"-", "_"} else "_" for ch in name
            )
            return safe or "Worklist"

        written = 0
        exported_frames: dict[str, pd.DataFrame] = {}

        for key, df in self.worklist_data.items():
            if df is None:
                continue

            if not isinstance(df, pd.DataFrame):
                LOGGER.warning("Skipping non-DataFrame worklist '%s'.", key)
                continue

            if not include_empty and df.empty:
                continue

            clean_df = normalize_worklist_well_columns(df)
            exported_frames[str(key)] = clean_df

            base = _sanitize(str(key))
            file_path = out_dir / f"{base}.csv"

            if not overwrite and file_path.exists():
                index = 1
                while True:
                    candidate = out_dir / f"{base}_{index}.csv"
                    if not candidate.exists():
                        file_path = candidate
                        break
                    index += 1

            try:
                clean_df.to_csv(file_path, index=False, sep=",")
                written += 1
                LOGGER.info("Exported worklist '%s' -> %s", key, file_path)
            except Exception as err:  # pylint: disable=broad-except
                LOGGER.error(
                    "Failed to export worklist '%s' to %s: %s",
                    key,
                    file_path,
                    err,
                )

        if export_master:
            try:
                master_df = build_master_worklist(exported_frames)

                if include_empty or not master_df.empty:
                    master_path = out_dir / master_filename

                    if not overwrite and master_path.exists():
                        stem = Path(master_filename).stem
                        suffix = Path(master_filename).suffix or ".csv"
                        index = 1
                        while True:
                            candidate = out_dir / f"{stem}_{index}{suffix}"
                            if not candidate.exists():
                                master_path = candidate
                                break
                            index += 1

                    master_df.to_csv(master_path, index=False, sep=",")
                    LOGGER.info(
                        "Exported master worklist with %d rows -> %s",
                        len(master_df),
                        master_path,
                    )
                else:
                    LOGGER.info(
                        "Master worklist not written because no rows were available."
                    )
            except Exception as err:  # pylint: disable=broad-except
                LOGGER.error("Failed to export master worklist CSV: %s", err)

        LOGGER.info(
            "Worklist CSV export done: %d individual file(s) -> %s",
            written,
            out_dir,
        )
        return out_dir

    # private methods

    def _normalize_part_name(self, value: object) -> str:
        """Normalize part names for robust comparison."""
        if value is None or pd.isna(value):
            return ""
        return " ".join(str(value).strip().split())

    def _check_for_mismatched_parts(self) -> pd.DataFrame:
        """Check if all planned parts exist EXACTLY in source tables."""

        # --- Collect source parts ---
        source_frames = []

        def _extract(df: pd.DataFrame, label: str) -> pd.DataFrame:
            col = df.columns[0]
            out = df[[col]].copy()
            out = out.rename(columns={col: "Part"})
            out["Source"] = label
            out["Part_norm"] = out["Part"].map(self._normalize_part_name)
            return out

        if isinstance(self.stocks_information, pd.DataFrame):
            source_frames.append(_extract(self.stocks_information, "Stocks"))

        if isinstance(self.mastermix_information, pd.DataFrame):
            source_frames.append(_extract(self.mastermix_information, "Mastermix"))

        if isinstance(self.water_information, pd.DataFrame):
            source_frames.append(_extract(self.water_information, "Water"))

        if isinstance(self.sourceplate_information, pd.DataFrame):
            source_frames.append(_extract(self.sourceplate_information, "Source Plate"))

        if not source_frames:
            LOGGER.warning("No source data available.")
            return pd.DataFrame()

        source_df = pd.concat(source_frames, ignore_index=True)
        source_set = set(source_df["Part_norm"])

        # --- Collect planned parts ---
        dest_frames = []

        if isinstance(self.df_dest_plate, pd.DataFrame):
            for col in self.df_dest_plate.columns:
                if col == "Well":
                    continue

                df = (
                    self.df_dest_plate[["Well", col]]
                    .dropna()
                    .rename(columns={col: "Part"})
                )
                df["Source"] = f"Dest::{col}"
                dest_frames.append(df)

        if isinstance(self.df_pm, pd.DataFrame):
            melted = self.df_pm.melt(id_vars=["Well"], value_name="Part").dropna(
                subset=["Part"]
            )
            melted["Source"] = "PM"
            dest_frames.append(melted[["Well", "Part", "Source"]])

        if not dest_frames:
            return pd.DataFrame()

        dest_df = pd.concat(dest_frames, ignore_index=True)
        dest_df["Part_norm"] = dest_df["Part"].map(self._normalize_part_name)

        # --- Check mismatches ---
        mismatches = dest_df[~dest_df["Part_norm"].isin(source_set)].copy()

        if mismatches.empty:
            return pd.DataFrame()

        return mismatches[["Part", "Well", "Source"]]

    def _log_mismatched_parts(self, *, raise_on_mismatch: bool = True) -> None:
        """Log and optionally stop on mismatched parts."""

        mismatches = self._check_for_mismatched_parts()

        if mismatches.empty:
            LOGGER.info("All parts matched correctly.")
            return

        LOGGER.error(
            "Mismatched parts detected:\n%s",
            mismatches.to_string(index=False),
        )

        if raise_on_mismatch:
            raise ValueError(
                "Unmatched parts found. Fix spelling in Plan/PM before continuing."
            )

    def _sum_mastermix_volumes(self) -> float:
        """Sum all volumes for components containing 'mastermix' (case-insensitive).

        Returns
        -------
        float
            Total volume over rows whose ``Component`` contains ``mastermix``.
            Returns ``0.0`` if no volume table is available.
        """
        if not isinstance(self.df_volumes, pd.DataFrame) or self.df_volumes.empty:
            return 0.0
        mask = self.df_volumes["Component"].str.contains(
            "mastermix", case=False, na=False
        )
        return float(self.df_volumes.loc[mask, "Total Volume (µL)"].sum())

    def _combine_and_filter_parts(
        self,
        df_dest_plate: pd.DataFrame,
        df_pm: pd.DataFrame,
    ) -> pd.DataFrame:
        """Combine destination-plate parts with plasmid mix information.

        Removes mix-specific columns (``Mastermix``/``Plasmidmix``), joins with
        ``df_pm`` on ``Well``, sets ``Well`` as index, and drops all-NaN columns.

        Parameters
        ----------
        df_dest_plate : pandas.DataFrame
            Destination plate data; must include ``'Mastermix'`` and ``'Well'``.
        df_pm : pandas.DataFrame
            PM information; must include ``'Well'``.

        Returns
        -------
        pandas.DataFrame
            Combined table with ``Well`` as index and no all-NaN columns.

        Raises
        ------
        ValueError
            If required columns are missing.
        """
        if not {"Mastermix", "Well"}.issubset(df_dest_plate.columns):
            raise ValueError(
                "df_dest_plate is missing required columns 'Mastermix' or 'Well'."
            )
        if "Well" not in df_pm.columns:
            raise ValueError("df_pm is missing required column 'Well'.")

        filtered = df_dest_plate.drop(
            columns=["Mastermix", "Plasmidmix"], errors="ignore"
        )
        combined = pd.merge(filtered, df_pm, on="Well", how="left")
        combined.set_index("Well", inplace=True)
        combined.dropna(axis=1, how="all", inplace=True)
        return combined

    def _load_metadata(self) -> pd.DataFrame:
        """Load the metadata sheet from the workbook.

        Uses the helper :func:`load_metadata_dataframe` with sheet
        ``"Plate organization+"`` and drops empty rows.

        Returns
        -------
        pandas.DataFrame

        Raises
        ------
        ValueError
            If `file_path` is not set or the sheet cannot be loaded.
        """
        if not self.file_path:
            raise ValueError("file_path is not set.")
        try:
            df = load_metadata_dataframe(
                self.file_path, sheet_name="Plate organization+"
            )
            return df.dropna(how="all")
        except Exception as err:  # pylint: disable=broad-except
            raise ValueError(f"Failed to load metadata: {err}") from err

    def _load_plan(self) -> pd.DataFrame:
        """Load and sanitize the destination plan sheet.

        Reads the plan via :func:`read_and_process_dataframe`, limits rows to a 96-well
        plate, converts empty strings to ``None``, and drops rows where all non-``Well``
        columns are empty.

        Returns
        -------
        pandas.DataFrame

        Raises
        ------
        ValueError
            If reading or processing fails.
        """
        if not self.file_path or not self.sheet_name:
            raise ValueError("file_path and sheet_name must be set.")
        try:
            df = read_and_process_dataframe(
                self.file_path,
                self.sheet_name,
                start_column="B",
                well_column="Well",
            )
            df = df.iloc[:96].copy()
            df = df.applymap(
                lambda x: None if isinstance(x, str) and x.strip() == "" else x
            )
            non_well = [c for c in df.columns if c != "Well"]
            df = df.dropna(how="all", subset=non_well)
            return df
        except Exception as err:  # pylint: disable=broad-except
            raise ValueError(f"Failed to load plan: {err}") from err

    def _plot_mtp_volumes(self, plot: bool) -> None:
        """Plot MTP volumes for steps with ``start == 'Yes'`` if `plot` is True.

        Parameters
        ----------
        plot : bool
            If ``True``, calls ``mtp_manager.display_mtp_volumes`` for each
            starting source plate.
        """
        if not plot or self.workflow is None or self.mtp_manager is None:
            return
        for step in self.workflow.list_steps():
            if str(step.get("start", "")).lower() == "yes":
                try:
                    mtp = self.mtp_manager.get_mtp(step["source_plate"])
                    if mtp is not None:
                        self.mtp_manager.display_mtp_volumes(list_of_mtps=[mtp.name])
                except Exception as err:  # pylint: disable=broad-except
                    LOGGER.warning("Failed to plot volumes for step %s: %s", step, err)

    def _load_df(
        self,
        usecols: Sequence[int] | Sequence[str],
        sheet_name: str = "Source_Plate",
    ) -> pd.DataFrame:
        """Load a sheet from the Excel workbook with selected columns.

        Parameters
        ----------
        usecols : sequence of int or str
            Columns to read (as indices or names).
        sheet_name : str, optional
            Excel sheet name. Default ``"Source_Plate"``.

        Returns
        -------
        pandas.DataFrame
            DataFrame with empty rows dropped.

        Raises
        ------
        ValueError
            If `file_path` is not set or reading fails.
        """
        if not self.file_path:
            raise ValueError("file_path is not set.")
        try:
            df = pd.read_excel(
                self.file_path,
                sheet_name=sheet_name,
                usecols=list(usecols),
                header=0,
                index_col=None,
            )
            return df.dropna(how="all")
        except Exception as err:  # pylint: disable=broad-except
            raise ValueError(
                f"Failed to read sheet '{sheet_name}' with usecols {list(usecols)}: {err}"
            ) from err

    def _find_worklist_sheets(self) -> list[str]:
        """Return all Excel sheet names that look like worklists.

        Uses :func:`find_worklist_sheets` with default pattern (``"Worklist"``).

        Returns
        -------
        list of str

        Raises
        ------
        ValueError
            If `file_path` is not set or the workbook cannot be inspected.
        """
        if not self.file_path:
            raise ValueError("file_path is not set.")
        try:
            sheets = find_worklist_sheets(self.file_path)
            LOGGER.info("Found %d worklist sheet(s): %s", len(sheets), sheets)
            return sheets
        except Exception as err:  # pylint: disable=broad-except
            raise ValueError(f"Failed to find worklist sheets: {err}") from err

    def _create_empty_worklist_files(self, worklist_sheets: Sequence[str]) -> None:
        """Create empty ``.gwl`` files for each sheet name in `worklist_sheets`.

        Parameters
        ----------
        worklist_sheets : sequence of str
            Each entry is used as a base name for an empty ``.gwl`` file.
        """
        if not worklist_sheets:
            return
        try:
            create_empty_gwl_files(list(worklist_sheets))
            LOGGER.info(
                "Created empty .gwl files for %d sheet(s).", len(worklist_sheets)
            )
        except Exception as err:  # pylint: disable=broad-except
            raise ValueError(f"Failed to create empty GWL files: {err}") from err

    def _load_worklist_data(
        self,
        worklist_sheets: Sequence[str],
    ) -> dict[str, pd.DataFrame]:
        """Load all worklist sheets into a dict of DataFrames.

        Parameters
        ----------
        worklist_sheets : sequence of str
            Sheet names to load.

        Returns
        -------
        dict[str, pandas.DataFrame]
            Mapping of sheet name → DataFrame.

        Raises
        ------
        ValueError
            If `file_path` is not set or loading fails.
        """
        if not self.file_path:
            raise ValueError("file_path is not set.")
        try:
            data = load_worklist_sheets(self.file_path, list(worklist_sheets))
            LOGGER.info("Loaded %d worklist table(s).", len(data))
            return data
        except Exception as err:  # pylint: disable=broad-except
            raise ValueError(f"Failed to load worklist sheets: {err}") from err
