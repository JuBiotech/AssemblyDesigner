"""Batch orchestration (many plasmids × many reads).

This module wires together:
- Plasmid loading (.gb/.gbk), feature-concatenation, K-mer indexing
- Read loading (FASTA/FASTQ file OR a folder with multiple files)
- Seeding to select top-K candidate plasmids per read
- Read-anchored alignment on candidate windows (both strands)
- Result table assembly (one best-hit row per read), plus optional SNP dump

Notes
-----
- PID uses the **read** as denominator → stable across different read lengths.
- For performance: reuse per-plasmid concat/k-index across all reads.
- The caller can parallelize externally if needed; this module keeps things simple/reliable.
"""

from __future__ import annotations

import math
import os
import time
import traceback
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import (
    Any,
    Literal,
    Optional,
    Protocol,
    Union,
    cast,
)

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ._shared import _feature_pid_stats, _format_feature_hits_from_alignment
from .alignment import (
    _format_feature_hits,
    align_read_to_concat,
    build_kmer_index,
    seed_window,
)
from .features import (
    DEFAULT_ALLOWED_TYPES,
    build_feature_concat_sequence,
    derive_construct_from_record,
    feature_display_name,
)

try:
    from snapgene_reader import snapgene_file_to_dict
except Exception:  # pragma: no cover
    snapgene_file_to_dict = None

try:

    _HAS_TQDM = True
except Exception:
    _HAS_TQDM = False


@dataclass(slots=True)
class PlasmidRef:
    """Prepared plasmid reference.

    Attributes
    ----------
    file :
        Basename of the GenBank file (e.g., ``"pFoo.gb"``).
    plasmid_id :
        Record ID from the GenBank file.
    construct :
        Human-readable construct name derived from the record (project-specific).
    concat_ref :
        Feature-concatenated reference sequence (uppercase).
    feature_map :
        List of tuples ``(start, end, ftype, fname)`` in **concat_ref** coordinates.
    k_index :
        K-mer index over **concat_ref** for seeding (implementation specific).
    full_ref :
        **Full continuous plasmid sequence** (uppercase) for whole-band alignment.
    """

    file: str
    plasmid_id: str
    construct: str
    concat_ref: str
    feature_map: list[tuple[int, int, str, str]]
    k_index: dict[str, list[int]]
    full_ref: str


@dataclass(frozen=True)
class ReadItem:
    """
    Single input read.

    Attributes
    ----------
    name : str
        Identifier (e.g., FASTA/FASTQ record id).
    seq : str
        Sequence (uppercase recommended).
    """

    name: str
    seq: str


class PrefLike(Protocol):
    """Minimal interface of a plasmid reference expected by helpers."""

    concat_ref: str  # full linearized reference (concatenated)
    feature_map: Any  # your feature mapping structure


@dataclass(slots=True)
class ResultRow:
    """
    One summarized alignment result for a read (best hit kept).

    Parameters
    ----------
    sequence_name : str
        Read identifier.
    read_len : int
        Length of the read sequence in bases.
    plasmid_file : str
        Source file/path of the plasmid reference.
    plasmid_id : str
        Plasmid identifier (e.g., record name or stable id).
    construct : str
        Construct/assembly name this plasmid belongs to.
    strand : str
        Orientation used for the best alignment, e.g., '+' or '-'.
    pid : float
        Percent identity for the chosen PID mode.
    core_len : int
        Length of the aligned core (bp).
    score : float
        Raw alignment score from the backend.
    start_ref : int
        Start coordinate on the concatenated reference (inclusive).
    end_ref : int
        End coordinate on the concatenated reference (inclusive).
    snps : int
        Number of single-nucleotide mismatches detected.
    features_hit : str
        Compact per-feature coverage string (human-readable).
    concat_pos : str
        Preformatted reference window string, e.g. "123–456".

    min_feature_pid : float, optional
        (Compatibility) Minimum per-feature PID observed for this hit.
        Set to None if not computed.
    rbs_pid : float, optional
        (Compatibility) PID of the RBS feature, if present and computed.
        Set to None if not computed.
    n_features_below_95 : int, optional
        (Compatibility) Count of features with PID < 95%.
        Set to None if not computed.
    suspicious : bool, optional
        (Compatibility) Heuristic flag (e.g., high read PID but low feature PID).
        Set to None if not computed.
    """

    sequence_name: str
    read_len: int
    plasmid_file: str
    plasmid_id: str
    construct: str
    strand: str
    pid: float
    core_len: int
    score: float
    start_ref: int
    end_ref: int
    snps: int
    features_hit: str
    concat_pos: str

    # Optional/compat fields (safe defaults so older code still works)
    min_feature_pid: Optional[float] = None
    rbs_pid: Optional[float] = None
    n_features_below_95: Optional[int] = None
    suspicious: Optional[bool] = None


def _is_genbank(path: str, file_exts: tuple[str, ...]) -> bool:
    """Return True if `path` looks like a GenBank file by extension."""
    low = path.lower()
    return any(low.endswith(ext) for ext in file_exts)


def _guess_seq_format(path: str) -> Optional[str]:
    """
    Detect sequence file format by extension.

    Returns
    -------
    str or None
        "fasta", "fastq", "dna" (SnapGene) or None if unknown.
    """
    low = path.lower()
    if low.endswith((".fa", ".fasta", ".fna")):
        return "fasta"
    if low.endswith((".fq", ".fastq")):
        return "fastq"
    if low.endswith(".dna"):
        return "dna"
    return None


def prepare_plasmids(
    plasmids_folder: str,
    *,
    file_exts: tuple[str, ...] = (".gb", ".gbk"),
    include_tokens: Optional[Sequence[str]] = None,
    exclude_tokens: Optional[Sequence[str]] = ("backbone", "vector", "ori", "dvl1"),
    allowed_types: Sequence[str] = DEFAULT_ALLOWED_TYPES,
    k: int = 11,
    include_terminators: bool = True,
    terminator_tokens: Sequence[str] = ("terminator", "B0015"),
    use_full_sequence: bool = False,
) -> list[PlasmidRef]:
    """
    Load GenBank plasmids and prepare references.

    When `use_full_sequence=False` (default), build a feature-concatenated
    reference (`concat_ref`) plus a K-mer index over it. When
    `use_full_sequence=True`, set `concat_ref` to the *full plasmid sequence*
    (i.e., whole-band alignment target) and build the K-mer index over that.

    In both cases, `full_ref` is stored alongside for downstream use.

    Returns
    -------
    list[PlasmidRef]
        Each item has fields: file, plasmid_id, construct, concat_ref,
        feature_map (coordinates match `concat_ref`), k_index, full_ref.
    """

    # --- effective exclude list (optionally expand for terminators) ----------
    if not include_terminators:
        base = list(exclude_tokens or ())
        seen = {t.lower() for t in base}
        for tok in terminator_tokens:
            if tok.lower() not in seen:
                base.append(tok)
        eff_exclude: tuple[str, ...] = tuple(base)
    else:
        eff_exclude = tuple(exclude_tokens or ())

    def _name_ok(name: str) -> bool:
        n = name.lower()
        if include_tokens and not any(t.lower() in n for t in include_tokens):
            return False
        if eff_exclude and any(t.lower() in n for t in eff_exclude):
            return False
        return True

    def _feature_map_full(
        rec: SeqRecord,
        allowed: Sequence[str],
    ) -> list[tuple[int, int, str, str]]:
        fmap: list[tuple[int, int, str, str]] = []
        for feat in rec.features:
            if allowed and getattr(feat, "type", "") not in allowed:
                continue
            name = feature_display_name(feat)
            if not _name_ok(name):
                continue
            try:
                s = int(feat.location.start)  # 0-based, half-open
                e = int(feat.location.end)
            except Exception:
                continue
            fmap.append(
                (max(0, s), max(0, e), getattr(feat, "type", "misc_feature"), name)
            )
        return fmap

    out: list[PlasmidRef] = []

    for fname in sorted(os.listdir(plasmids_folder)):
        path = os.path.join(plasmids_folder, fname)
        if not os.path.isfile(path) or not _is_genbank(path, file_exts):
            continue
        if os.path.getsize(path) <= 0:
            continue

        try:
            rec: SeqRecord = SeqIO.read(path, "genbank")
        except Exception:
            continue

        # Always keep the full continuous plasmid sequence
        full_ref = str(rec.seq).upper()

        if use_full_sequence:
            # Use FULL plasmid as concat_ref; feature map in genomic coords
            fmap = _feature_map_full(rec, allowed_types)
            if not fmap:
                # no selectable features -> synthetic catch-all
                fmap = [
                    (0, len(full_ref), "misc_feature", feature_display_name(f))
                    for f in rec.features
                ] or [(0, len(full_ref), "misc_feature", "record")]
            concat_ref = full_ref
        else:
            # Build FEATURE-CONCAT reference + feature map in concat coords
            concat_ref, fmap = build_feature_concat_sequence(
                rec,
                include_tokens=include_tokens,
                exclude_tokens=eff_exclude,
                allowed_types=allowed_types,
            )
            if not concat_ref:
                # fallback: full sequence also as concat_ref
                concat_ref = full_ref
                fmap = [
                    (0, len(concat_ref), "misc_feature", feature_display_name(f))
                    for f in rec.features
                ] or [(0, len(concat_ref), "misc_feature", "record")]

        # K-mer index over whatever we decided `concat_ref` is
        k_index = build_kmer_index(concat_ref, k=k)
        construct = derive_construct_from_record(rec)

        out.append(
            PlasmidRef(
                file=os.path.basename(path),
                plasmid_id=str(rec.id),
                construct=construct,
                concat_ref=concat_ref,
                feature_map=fmap,
                k_index=k_index,
                full_ref=full_ref,  # always available
            )
        )

    return out


def dealign_reads_easy(
    reads_input: Union[str, Path],
    plasmids_dir: Union[str, Path],
    *,
    # --- Feature selection ----------------------------------------------------
    include_tokens: Optional[Sequence[str]] = None,
    exclude_tokens: Optional[Sequence[str]] = ("dvl1", "backbone", "vector", "ori"),
    allowed_types: Optional[Sequence[str]] = None,
    # --- Seeding / alignment knobs -------------------------------------------
    k: int = 11,
    step: int = 4,
    margin: int = 120,
    within: Optional[float] = 0.90,
    top_k: Optional[int] = 5,
    try_both_strands: bool = True,
    pid_mode: Literal["read", "feature", "feature_len"] = "read",
    match: float = 1.0,
    mismatch: float = -1.0,
    gap_open: float = -1.5,
    gap_extend: float = -0.5,
    # --- Parallel execution ---------------------------------------------------
    backend: Literal["threads", "processes", "none"] = "threads",
    workers: Optional[int] = None,
    per_read_timeout_s: Optional[float] = 10.0,
    max_retries: int = 2,
    fallback_kwargs: Optional[dict[str, Any]] = None,
    # NEW: light-weight candidate prefilter (parallel path only)
    min_jaccard: float = 0.10,
    # --- Read trimming --------------------------------------------------------
    trim_left: int = 0,
    trim_right: int = 0,
    # --- Outputs --------------------------------------------------------------
    out_table: Optional[Union[str, Path]] = None,
    out_format: Literal["xlsx", "csv"] = "xlsx",
    snp_report_path: Optional[Union[str, Path]] = None,
    # --- UX -------------------------------------------------------------------
    verbose: bool = True,
    progress: bool = True,
    # --- Aligner control & full-pass fallback --------------------------------
    read_aligner: Union[str, Callable[..., Any], None] = "auto",
    fallback_read_aligner: Union[str, Callable[..., Any], None] = "pairwise2",
    auto_fallback_on_exception: bool = True,
    auto_fallback_on_zero_hits: bool = True,
    min_matched_reads_to_accept: int = 1,
) -> tuple[list[ResultRow], pd.DataFrame]:
    """
    One-call wrapper to align reads against plasmid references, with an optional
    full-pass fallback between aligners.

    The function prepares references, loads reads, runs one full alignment pass
    (threads or processes), and optionally re-runs with a fallback aligner if
    the first pass raises or yields too few matches.

    Parameters
    ----------
    reads_input : str or pathlib.Path
        FASTA/FASTQ/SnapGene(.dna) file **or** directory containing such files.
    plasmids_dir : str or pathlib.Path
        Directory with GenBank plasmids used to build feature-concatenated references.
    include_tokens, exclude_tokens : Sequence[str] or None, optional
        Case-insensitive substrings that filter features by *display name*.
        Excludes typical backbone tokens by default.
    allowed_types : Sequence[str] or None, optional
        GenBank feature types to retain. If ``None``, a safe default is used.
    k : int, default 11
        K-mer length for seeding and index construction.
    step : int, default 4
        Seeding stride in bp (slide the read k-mer every ``step`` bases).
    margin : int, default 120
        Extra bases added left/right to the seeded window (robustness vs. speed).
    within : float or None, default 0.90
        Dynamic gating: keep candidates with ``hits ≥ ceil(within * best_hits)``.
        ``None`` disables the gate.
    top_k : int or None, default 5
        Keep at most this many candidates per read (after gating).
    try_both_strands : bool, default True
        Evaluate forward and reverse-complement; keep the better.
    pid_mode : {"read","feature","feature_len"}, default "read"
        Returned PID definition (read-based by default). Feature modes require
        a low-level aligner that accepts ``pid_mode``, ``feature_map``,
        and ``ref_abs_start``; this is provided by :func:`align_read_to_concat`.
    match, mismatch, gap_open, gap_extend : float
        Scoring knobs for pairwise fallbacks (edlib ignores these).
    backend : {"threads","processes","none"}, default "threads"
        Execution backend. ``"none"`` runs in the current thread (no pool).
    workers : int or None, optional
        Number of worker threads/processes. Defaults to ``min(8, cpu_count())``.
    per_read_timeout_s : float or None, default 10.0
        Per-read timeout; only *hard* for ``backend="processes"``.
    max_retries : int, default 2
        Maximum retries for timed-out reads (parallel backends).
    fallback_kwargs : dict or None, optional
        Parameter overrides applied on retries (e.g., smaller windows).
    min_jaccard : float, default 0.10
        K-mer Jaccard threshold (read vs. seeded ref segment). Only used in the
        *parallel* path via ``_align_one_read``. Set lower (e.g., 0.02–0.05) for
        noisy Sanger reads.
    trim_left, trim_right : int
        Fixed-base trimming applied to reads before alignment.
    out_table : str or pathlib.Path or None
        Optional output path for the result table.
    out_format : {"xlsx","csv"}, default "xlsx"
        Format for ``out_table``.
    snp_report_path : str or pathlib.Path or None
        Optional path to write a human-readable SNP report (single-thread use).
    verbose : bool, default True
        Print a header, progress, and a final summary.
    progress : bool, default True
        Show a tqdm bar if available; else periodic prints.
    read_aligner : {"auto","edlib","pairwise2"} or callable or None, default "auto"
        Primary per-read aligner. Tags are resolved downstream; use a small
        shim to map "pairwise2" to your PairwiseAligner-based implementation.
    fallback_read_aligner : same as ``read_aligner``, default "pairwise2"
        Fallback aligner used for a **full second pass** if enabled.
    auto_fallback_on_exception : bool, default True
        If the primary pass raises, re-run once with the fallback aligner.
    auto_fallback_on_zero_hits : bool, default True
        If primary matches are fewer than ``min_matched_reads_to_accept``,
        re-run once with the fallback aligner.
    min_matched_reads_to_accept : int, default 1
        Acceptance threshold for the primary pass.

    Returns
    -------
    rows : list[ResultRow]
        One best-hit row per read. Unmatched reads are omitted.
    df : pandas.DataFrame
        Result table (sorted by ``sequence_name``, ``pid``, ``score``) with a
        column ``pid_mode`` indicating the denominator used.

    Notes
    -----
    - For ``backend="processes"``, aligner specs are transported as *strings*
      (picklable) to workers. For ``threads``/``none``, callables are fine.
    - ``align_batch`` (single-thread) does **not** apply the Jaccard prefilter.
      The prefilter is only active in the parallel path through
      ``align_batch_parallel`` → ``_align_one_read``.
    """
    # Resolve allowed_types default robustly (do not crash if global is missing)
    if allowed_types is None:
        allowed_types = globals().get(
            "DEFAULT_ALLOWED_TYPES",
            ("promoter", "RBS", "CDS", "terminator", "misc_feature", "5'UTR", "3'UTR"),
        )

    # -- Helpers to normalize/transport aligner specs --------------------------
    def _name_for_callable(fn: Callable[..., Any]) -> str:
        nm = getattr(fn, "__name__", None) or str(fn)
        low = nm.lower()
        if "edlib" in low:
            return "edlib"
        if "pairwise" in low or "read_anchored_align" in low:
            return "pairwise2"
        return nm

    def _transport_spec(
        spec: Union[str, Callable[..., Any], None],
        is_processes: bool,
    ) -> Union[str, Callable[..., Any], None]:
        if spec is None:
            return None
        if isinstance(spec, str):
            return spec
        return _name_for_callable(spec) if is_processes else spec

    is_proc = backend == "processes"
    primary_spec = _transport_spec(read_aligner, is_proc)
    fallback_spec = _transport_spec(fallback_read_aligner, is_proc)

    if verbose:
        print("New aligner:", primary_spec is not None)
        print("Using:", primary_spec)

    # 1) Prepare plasmids
    plasmids = prepare_plasmids(
        str(plasmids_dir),
        include_tokens=include_tokens,
        exclude_tokens=exclude_tokens,
        allowed_types=allowed_types,
        k=k,
    )
    if verbose:
        print(f"Prepared {len(plasmids)} plasmid references.")

    # 2) Load reads
    reads = load_reads(
        str(reads_input),
        trim_left=trim_left,
        trim_right=trim_right,
    )
    if verbose:
        print(f"Loaded {len(reads)} reads.")

    # 3) Inner runner (so we can re-run with a different aligner spec)
    def _run_with(spec: Union[str, Callable[..., Any], None]) -> list[ResultRow]:
        """Run a single pass with the given aligner spec."""
        # Common knobs used by both single-thread and parallel paths
        common_kwargs: dict[str, Any] = {
            "k": k,
            "step": step,
            "margin": margin,
            "within": within,
            "top_k": top_k,
            "try_both_strands": try_both_strands,
            "pid_mode": pid_mode,
            "match": match,
            "mismatch": mismatch,
            "gap_open": gap_open,
            "gap_extend": gap_extend,
            # Pass aligner choice down to the workers/single-thread path
            "read_aligner": spec,
            "fallback_read_aligner": None,
        }

        if backend == "none":
            # Single-threaded path (no Jaccard prefilter)
            return align_batch(
                plasmids,
                reads,
                snp_report_path=str(snp_report_path) if snp_report_path else None,
                show_progress=progress,
                progress_desc="Aligning reads",
                **common_kwargs,
            )

        # Parallel path: add the Jaccard prefilter (used by _align_one_read)
        parallel_kwargs = dict(common_kwargs)
        parallel_kwargs["min_jaccard"] = float(min_jaccard)

        return align_batch_parallel(
            plasmids,
            reads,
            backend=("threads" if backend == "threads" else "processes"),
            workers=workers,
            verbose=verbose,
            progress=progress,
            per_read_timeout_s=per_read_timeout_s,
            max_retries=max_retries,
            fallback_kwargs=fallback_kwargs,
            **parallel_kwargs,
        )

    rows: list[ResultRow] = []

    # Primary run
    try:
        rows = _run_with(primary_spec)
    except Exception as exc:  # pylint: disable=broad-except
        if verbose:
            print(f"[plasmidio] Primary aligner failed: {type(exc).__name__}: {exc}")
        if auto_fallback_on_exception and fallback_spec is not None:
            if verbose:
                print(f"[plasmidio] Retrying with fallback aligner… {fallback_spec}")
            rows = _run_with(fallback_spec)
        else:
            raise

    # Low-yield fallback: full second pass with fallback aligner
    if (
        auto_fallback_on_zero_hits
        and (len(rows) < int(min_matched_reads_to_accept))
        and fallback_spec is not None
    ):
        if verbose:
            print(
                f"[plasmidio] Low yield ({len(rows)} matches). "
                f"Retrying with fallback aligner… {fallback_spec}"
            )
        rows = _run_with(fallback_spec)

    # 4) DataFrame (+ optional export)
    df = results_to_dataframe(rows).sort_values(
        ["sequence_name", "pid", "score"], ascending=[True, False, False]
    )
    df["pid_mode"] = pid_mode  # record which denominator was used

    if out_table:
        out_path = Path(out_table)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        print(rows[0].__dict__.keys())
        export_results(rows, str(out_path), fmt=out_format)
        if verbose:
            print(f"Wrote results → {out_path.resolve()}")

    return rows, df


def _iter_reads_from_file(
    path: str,
    *,
    trim_left: int = 0,
    trim_right: int = 0,
) -> Iterable[ReadItem]:
    """Yield ReadItem(s) from a FASTA/FASTQ/SnapGene (.dna) file."""
    fmt = _guess_seq_format(path)
    if fmt is None:
        return  # silently skip unknown formats

    def _apply_trim(seq: str) -> str:
        if not (trim_left or trim_right):
            return seq
        L = len(seq)  # noqa: N806
        left = max(0, min(trim_left, L))
        right = max(0, min(trim_right, max(0, L - left)))
        return seq[left : L - right] if right else seq[left:]

    if fmt in ("fasta", "fastq"):
        for rec in SeqIO.parse(path, fmt):
            seq = _apply_trim(str(rec.seq).upper())
            yield ReadItem(name=str(rec.id), seq=seq)
        return

    if fmt == "dna":
        if snapgene_file_to_dict is None:
            raise ImportError(
                "Reading .dna requires 'snapgene_reader'. Install via: pip install snapgene-reader"
            )
        data = snapgene_file_to_dict(path)
        raw = data.get("seq", b"")
        if isinstance(raw, (bytes, bytearray)):
            seq = raw.decode("ascii", errors="ignore").upper()
        else:
            seq = str(raw).upper()
        name = str(data.get("name") or os.path.splitext(os.path.basename(path))[0])
        seq = _apply_trim(seq)
        # Treat a .dna file as a single 'read' item
        yield ReadItem(name=name, seq=seq)
        return


def load_reads(
    input_path: str,
    *,
    trim_left: int = 0,
    trim_right: int = 0,
) -> list[ReadItem]:
    """
    Load reads from a FASTA/FASTQ/SnapGene (.dna) file **or** from a folder of such files.

    Parameters
    ----------
    input_path : str
        Path to a .fasta/.fastq/.dna file or a folder containing multiple files.
    trim_left, trim_right : int, optional
        Fixed-base trimming at the ends (before alignment).

    Returns
    -------
    list[ReadItem]
        All reads in the given input (order preserved within each file).
    """
    if os.path.isdir(input_path):
        items: list[ReadItem] = []
        for fname in sorted(os.listdir(input_path)):
            path = os.path.join(input_path, fname)
            if not os.path.isfile(path) or _guess_seq_format(path) is None:
                continue
            items.extend(
                list(
                    _iter_reads_from_file(
                        path, trim_left=trim_left, trim_right=trim_right
                    )
                )
            )
        return items

    # Single file
    if not os.path.isfile(input_path):
        raise FileNotFoundError(f"Not a file or folder: {input_path}")
    return list(
        _iter_reads_from_file(input_path, trim_left=trim_left, trim_right=trim_right)
    )


def _topk_candidates_for_read(
    read: ReadItem,
    plasmids: Sequence[PlasmidRef],
    *,
    k: int,
    margin: int,
    step: int,
    top_k: Optional[int] = 4,
) -> list[tuple[int, PlasmidRef, tuple[int, int]]]:
    """
    Rank plasmids for a read using k-mer seeding and return the top-K windows.

    Parameters
    ----------
    read : ReadItem
        Input read (expects attribute ``seq``).
    plasmids : Sequence[PlasmidRef]
        Prepared references (each exposes ``concat_ref`` and ``k_index``).
    k : int
        K-mer length used for seeding. **Must match** how ``k_index`` was built.
    margin : int
        Extra bases to include left/right of the seeded median offset.
    step : int
        Stride for sliding the read k-mer window (every ``step`` bases).
    top_k : int or None, optional
        Keep at most this many candidates (after sorting by seed hits).
        If ``None``, all candidates are returned.

    Returns
    -------
    list of tuple
        Each entry is ``(hits, plasmid_ref, (start, end))`` where
        ``(start, end)`` is the seed window on the concatenated reference
        (0-based, end-exclusive). The list is sorted by ``hits`` descending.

    Notes
    -----
    - Windows come directly from :func:`seed_window`; sizing is controlled
      by ``margin`` and the read length.
    - If a plasmid yields no seed hits, it is skipped.
    - The caller is responsible for any **additional gating** (e.g., Jaccard
      prefilter, max window length), because those checks might differ between
      single-thread und parallelem Pfad.
    """
    cands: list[tuple[int, PlasmidRef, tuple[int, int]]] = []
    for pref in plasmids:
        win = seed_window(
            pref.concat_ref, read.seq, pref.k_index, k=k, margin=margin, step=step
        )
        if win is None:
            continue
        w_start, w_end, hits = win
        cands.append((hits, pref, (int(w_start), int(w_end))))
    cands.sort(key=lambda x: x[0], reverse=True)
    return cands if top_k is None else cands[: int(top_k)]


def _kmerset(s: str, k: int) -> set[str]:
    """
    Return the set of all k-mers in `s`.

    Parameters
    ----------
    s
        Input sequence.
    k
        k-mer length.

    Returns
    -------
    set of str
        Unique k-mers (empty if `len(s) < k`).
    """
    n = len(s) - k + 1
    if n <= 0:
        return set()
    return {s[i : i + k] for i in range(n)}


def _jaccard(a: str, b: str, k: int) -> float:
    """
    Jaccard index of k-mer sets between two strings.

    Parameters
    ----------
    a, b
        Sequences to compare.
    k
        k-mer length.

    Returns
    -------
    float
        |A∩B| / |A∪B| in [0, 1]. Returns 0.0 if either side has no k-mers.
    """
    A = _kmerset(a, k)  # noqa: N806
    B = _kmerset(b, k)  # noqa: N806
    if not A or not B:
        return 0.0
    inter = len(A & B)
    uni = len(A | B)
    return inter / uni


def align_batch(
    plasmids: Sequence[PlasmidRef],
    reads: Sequence[ReadItem],
    *,
    k: int = 11,
    margin: int = 200,
    step: int = 2,
    top_k: int = 4,
    try_both_strands: bool = True,
    match: float = 1.0,
    mismatch: float = -1.0,
    gap_open: float = -1.5,
    gap_extend: float = -0.5,
    snp_report_path: Optional[str] = None,
    pid_mode: Literal["read", "feature", "feature_len"] = "read",
    show_progress: bool = True,
    progress_desc: str = "Aligning reads",
    progress_print_every: Optional[int] = None,
    # Aligner selection
    read_aligner: Union[str, Callable[..., Any], None] = None,
    fallback_read_aligner: Union[str, Callable[..., Any], None] = None,
    # Direct low-level aligner (overrides read_aligner when set)
    read_align_impl: Optional[Callable[..., Any]] = None,
) -> list[ResultRow]:
    """
    Align a batch of reads against a set of plasmids (sequential, single-thread).

    The function uses k-mer seeding to pick top-K candidate plasmids per read,
    then performs a read-anchored semiglobal alignment within a seeded window.
    For each read, the best result is selected by highest percent identity (PID),
    with ties broken by the raw score.

    Parameters
    ----------
    plasmids : Sequence[PlasmidRef]
        Prepared references produced by :func:`prepare_plasmids`. Each item must
        expose ``file``, ``plasmid_id``, ``construct``, ``concat_ref``,
        ``feature_map``, and ``k_index``.
    reads : Sequence[ReadItem]
        Input reads. Each item must provide ``name`` (str) and ``seq`` (str).
    k : int, default=11
        K-mer length for seeding (must match reference preparation).
    margin : int, default=200
        Extra bases to add around the seeded window before alignment.
    step : int, default=2
        Stride for read k-mers during seeding.
    top_k : int, default=4
        Maximum number of candidate plasmids kept per read after seeding.
    try_both_strands : bool, default=True
        If True, evaluate both strands and keep the better alignment.
    match : float, default=1.0
        Match score (used by pairwise aligners).
    mismatch : float, default=-1.0
        Mismatch penalty (used by pairwise aligners).
    gap_open : float, default=-1.5
        Gap-open penalty (used by pairwise aligners).
    gap_extend : float, default=-0.5
        Gap-extend penalty (used by pairwise aligners).
    snp_report_path : str or None, optional
        If provided, append a simple SNP report for all reads to this path.
    pid_mode : {"read", "feature", "feature_len"}, default="read"
        PID calculation mode forwarded to :func:`align_read_to_concat`.
    show_progress : bool, default=True
        Show a tqdm progress bar if available; otherwise periodic prints.
    progress_desc : str, default="Aligning reads"
        Label for the tqdm bar or log prints.
    progress_print_every : int or None, optional
        When tqdm is unavailable, print once every N reads. If None, an
        interval is chosen to produce about 20 lines.
    read_aligner : {"auto","edlib","pairwise2"} or callable or None, optional
        Primary per-candidate aligner specification. If string, it is resolved
        from globals; if callable, passed directly. When None, a reasonable
        default (edlib→pairwise2) is used.
    fallback_read_aligner : same as ``read_aligner``, optional
        Fallback per-candidate aligner when the primary returns ``None``.
    read_align_impl : callable or None, optional
        Direct low-level aligner injected into :func:`align_read_to_concat`.
        When provided, it overrides ``read_aligner``.

    Returns
    -------
    list of ResultRow
        One ``ResultRow`` per read that yielded a successful alignment.
        Unmatched reads are omitted.

    Notes
    -----
    - This function is single-threaded and does not spawn workers.
    - It attempts to format the feature coverage string via
      :func:`_format_feature_hits_from_alignment`. If that helper is not
      available, it falls back to :func:`_format_feature_hits`, and as a last
      resort creates a minimal string from ``AlignmentResult.feature_hits``.
    """

    def _resolve_aligner(
        spec: Union[str, Callable[..., Any], None]
    ) -> Optional[Callable[..., Any]]:
        """Map a user-specified aligner spec to a callable from globals()."""
        if spec is None:
            return (
                globals().get("read_anchored_align_edlib")
                or globals().get("read_anchored_align_pairwise2")
                or globals().get("read_anchored_align")
            )
        if callable(spec):
            return spec
        name = str(spec).lower()
        if name == "auto":
            return (
                globals().get("read_anchored_align_edlib")
                or globals().get("read_anchored_align_pairwise2")
                or globals().get("read_anchored_align")
            )
        if name in ("ed", "edlib"):
            return globals().get("read_anchored_align_edlib")
        if name in ("pw2", "pairwise2", "biopython"):
            # Back-compat alias; may map to a PairwiseAligner-based path.
            return globals().get("read_anchored_align_pairwise2") or globals().get(
                "read_anchored_align"
            )
        # Last resort: exact global symbol name
        return globals().get(spec)

    primary_impl: Optional[Callable[..., Any]] = read_align_impl or _resolve_aligner(
        read_aligner
    )
    fallback_impl: Optional[Callable[..., Any]] = _resolve_aligner(
        fallback_read_aligner
    )

    rows: list[ResultRow] = []
    snp_buf: list[str] = []

    total = len(reads)

    # --- Progress handling (tqdm optional) -------------------------------------
    tqdm = globals().get("tqdm")
    has_tqdm_flag = bool(globals().get("_HAS_TQDM", False))
    if show_progress and total > 0 and (tqdm is None or not has_tqdm_flag):
        try:
            from tqdm import tqdm as _tqdm

            tqdm = _tqdm
            has_tqdm_flag = True
        except Exception:
            has_tqdm_flag = False

    use_tqdm = bool(show_progress and has_tqdm_flag and total > 0)

    if not use_tqdm and show_progress:
        if progress_print_every is None:
            progress_print_every = max(1, total // 20) if total > 0 else 1
        print(
            f"[plasmidio] {progress_desc}… (printing every {progress_print_every} reads)"
        )

    if use_tqdm:
        assert tqdm is not None  # implied by use_tqdm/has_tqdm_flag above
        iterable: Iterable[tuple[int, ReadItem]] = enumerate(
            tqdm(
                reads, total=total, unit="read", dynamic_ncols=True, desc=progress_desc
            ),
            start=1,
        )
    else:
        iterable = enumerate(reads, start=1)

    # --- Main loop --------------------------------------------------------------
    for i, rd in iterable:
        # 1) Seed & select top-K candidate plasmids for this read
        cands = _topk_candidates_for_read(
            rd, plasmids, k=k, margin=margin, step=step, top_k=top_k
        )
        if not cands:
            if (
                not use_tqdm
                and show_progress
                and progress_print_every is not None
                and (i % progress_print_every == 0 or i == total)
            ):
                pct = 100.0 * i / float(max(1, total))
                print(f"[plasmidio] {i}/{total} reads ({pct:.1f}%)")
            continue

        best_row: Optional[ResultRow] = None
        best_pid = -1.0
        best_score = float("-inf")

        # 2) Align read to each candidate and keep best by PID, then score
        for _hits, pref, _win in cands:
            # Primary attempt
            res = align_read_to_concat(
                concat_ref=pref.concat_ref,
                feature_map=pref.feature_map,
                read_seq=rd.seq,
                k_index=pref.k_index,
                k=k,
                margin=margin,
                step=step,
                try_both_strands=try_both_strands,
                match=match,
                mismatch=mismatch,
                gap_open=gap_open,
                gap_extend=gap_extend,
                pid_mode=pid_mode,
                read_align_impl=primary_impl,
            )

            # Fallback per-candidate
            if (
                res is None
                and fallback_impl is not None
                and fallback_impl is not primary_impl
            ):
                res = align_read_to_concat(
                    concat_ref=pref.concat_ref,
                    feature_map=pref.feature_map,
                    read_seq=rd.seq,
                    k_index=pref.k_index,
                    k=k,
                    margin=margin,
                    step=step,
                    try_both_strands=try_both_strands,
                    match=match,
                    mismatch=mismatch,
                    gap_open=gap_open,
                    gap_extend=gap_extend,
                    pid_mode=pid_mode,
                    read_align_impl=fallback_impl,
                )

            if res is None:
                continue

            # Optional SNP aggregation
            if snp_report_path is not None and res.snps:
                snp_buf.append(
                    f"[{rd.name}] {pref.file} ({pref.plasmid_id}) {res.strand} "
                    f"PID={res.pid:.2f}% core={res.core_len} ref={res.start_ref}-{res.end_ref}"
                )
                for s in res.snps:
                    snp_buf.append(
                        f"  SNP ref@{s.ref_pos} {s.ref_base}>{s.read_base}  read@{s.read_pos}"
                    )
                snp_buf.append("")

            # Feature coverage string
            try:
                feats = _format_feature_hits_from_alignment(
                    pref,
                    res.start_ref,
                    res.end_ref,
                    rd.seq,
                    res.strand,
                    match=match,
                    mismatch=mismatch,
                    gap_open=gap_open,
                    gap_extend=gap_extend,
                )
            except NameError:
                try:
                    feats = _format_feature_hits(res)
                except NameError:
                    if getattr(res, "feature_hits", None):
                        feats = ", ".join(
                            f"{getattr(h, 'ftype', 'feat')}:{getattr(h, 'fname', '')} "
                            f"({int(round(float(getattr(h, 'coverage_pct', 0))))}%)"
                            for h in res.feature_hits
                        )
                    else:
                        feats = ""

            row = ResultRow(
                sequence_name=rd.name,
                read_len=len(rd.seq),
                plasmid_file=pref.file,
                plasmid_id=pref.plasmid_id,
                construct=pref.construct,
                strand=res.strand,
                pid=res.pid,
                core_len=res.core_len,
                score=res.score,
                start_ref=res.start_ref,
                end_ref=res.end_ref,
                snps=len(res.snps),
                features_hit=feats,
                concat_pos=f"{res.start_ref}–{res.end_ref}",
            )

            if (row.pid > best_pid) or (row.pid == best_pid and row.score > best_score):
                best_pid = row.pid
                best_score = row.score
                best_row = row

        if best_row is not None:
            rows.append(best_row)

        if (
            not use_tqdm
            and show_progress
            and progress_print_every is not None
            and (i % progress_print_every == 0 or i == total)
        ):
            pct = 100.0 * i / float(max(1, total))
            print(f"[plasmidio] {i}/{total} reads ({pct:.1f}%)")

    # 3) Write SNP report if requested
    if snp_report_path is not None:
        try:
            with open(snp_report_path, "w", encoding="utf-8") as fh:
                fh.write("\n".join(snp_buf))
        except Exception:
            # best-effort reporting only
            pass

    return rows


def feature_quality_metrics_from_alignment(
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
    prefer: Union[
        Literal["auto", "edlib", "pairwise2"], str, Callable[..., Any]
    ] = "auto",
) -> tuple[float, float, int]:
    """
    Convenience wrapper: compute `min_feature_pid`, `rbs_pid`, `n_features_below_95`
    using the same span/engine as the pretty-printer.

    Returns
    -------
    (min_feature_pid, rbs_pid, n_features_below_95) : tuple
        `rbs_pid` is NaN if no RBS-like feature had denom>0.
    """
    # Reuse the alignment done inside the formatter by recomputing here as well.
    # For code clarity we re-run the same internal engine logic.
    # (If you want to avoid double-work, factor the engine-selection into a shared helper.)
    # Reverse complement if needed
    if strand == "R":
        tbl = str.maketrans("ACGTNacgtn", "TGCANtgcan")
        read_seq = read_seq.translate(tbl)[::-1]

    ref_seg = str(pref.concat_ref)[int(start_ref) : int(end_ref)]
    read = read_seq.upper()

    # Try engines
    def _normalize_prefer(p: Union[str, Callable[..., Any]]) -> str:
        nm = (getattr(p, "__name__", str(p)) if callable(p) else str(p)).lower()
        if "edlib" in nm:
            return "edlib"
        if "pairwise" in nm or "pairwise2" in nm or "biopython" in nm:
            return "pairwise2"
        return "auto" if nm not in {"edlib", "pairwise2"} else nm

    tag = _normalize_prefer(prefer)

    a_ref: Optional[str] = None
    a_read: Optional[str] = None

    try:
        import edlib
    except Exception:
        edlib = None

    if tag in ("edlib", "auto") and edlib is not None:
        try:
            res = edlib.align(read, ref_seg, mode="HW", task="path")
            cigar = res.get("cigar", "")
            locs = res.get("locations", [])
            if cigar and locs:
                start = int(locs[0][0])

                # build gapped strings (same as above)
                def _ops(cig: str) -> Iterable[tuple[int, str]]:
                    n = 0
                    for ch in cig:
                        if ch.isdigit():
                            n = n * 10 + (ord(ch) - 48)
                        else:
                            yield (n if n else 1, ch)
                            n = 0

                i_ref = start
                i_read = 0
                r_buf: list[str] = []
                q_buf: list[str] = []
                if start > 0:
                    r_buf.append(ref_seg[:start])
                    q_buf.append("-" * start)
                for length, op in _ops(cigar):
                    if op == "M":
                        r_buf.append(ref_seg[i_ref : i_ref + length])
                        q_buf.append(read[i_read : i_read + length])
                        i_ref += length
                        i_read += length
                    elif op == "I":
                        r_buf.append("-" * length)
                        q_buf.append(read[i_read : i_read + length])
                        i_read += length
                    elif op == "D":
                        r_buf.append(ref_seg[i_ref : i_ref + length])
                        q_buf.append("-" * length)
                        i_ref += length
                if i_ref < len(ref_seg):
                    r_buf.append(ref_seg[i_ref:])
                    q_buf.append("-" * (len(ref_seg) - i_ref))
                a_ref, a_read = "".join(r_buf), "".join(q_buf)
        except Exception:
            a_ref = a_read = None

    if (a_ref is None or a_read is None) and tag in ("pairwise2", "auto"):
        try:
            from Bio import pairwise2

            alns = pairwise2.align.globalms(
                ref_seg,
                read,
                match,
                mismatch,
                gap_open,
                gap_extend,
                penalize_end_gaps=(False, True),
                one_alignment_only=True,
            )
            if alns:
                a = alns[0]
                a_ref, a_read = str(a.seqA), str(a.seqB)
        except Exception:
            a_ref = a_read = None

    # If no engine available, return safe defaults
    if a_ref is None or a_read is None:
        return 100.0, float("nan"), 0

    stats = _feature_pid_stats(pref, int(start_ref), int(end_ref), a_ref, a_read)

    min_pid = 100.0
    rbs_pid = float("nan")
    below95 = 0

    for _fs, _fe, ftype, fname, pid, denom, _cov in stats:
        if denom == 0:
            continue
        if pid < min_pid:
            min_pid = pid
        if pid < 95.0:
            below95 += 1
        if (ftype.lower() == "rbs") or ("bcd" in fname.lower()):
            rbs_pid = max(pid, rbs_pid) if not math.isnan(rbs_pid) else pid

    if min_pid == 100.0 and all(denom == 0 for *_, denom, _ in stats):
        # no measurable features; keep consistent defaults
        min_pid = 100.0

    return float(min_pid), float(rbs_pid), int(below95)


def _resolve_read_aligner_tag(
    tag: Optional[Union[str, Callable[..., Any]]]
) -> Optional[Callable[..., Any]]:
    """
    Resolve an aligner *tag* (string alias or callable) to a callable.

    The returned callable is expected to implement the per-segment API::

        (ref_segment: str, read_seq: str, **kwargs) ->
            (pid, core_len, start_ref, end_ref, start_read, end_read, score, snps)

    Parameters
    ----------
    tag : {None, str, callable}
        - ``None``: no preference; return ``None`` (caller may decide a default).
        - ``str``: alias for a known implementation:
            * ``"auto"``         → prefer Pairwise-first wrapper, else pairwise2 shim, else edlib
            * ``"pairwise"``     → pairwise2 shim (Biopython PairwiseAligner under the hood)
            * ``"pairwise2"``    → same as ``"pairwise"``
            * ``"biopython"``    → same as ``"pairwise"``
            * ``"ed"``/``"edlib"`` → edlib implementation
          Any other string is resolved via ``globals()[name]`` if present.
        - ``callable``: returned as-is.

    Returns
    -------
    callable or None
        Resolved implementation, or ``None`` if no suitable symbol exists.

    Notes
    -----
    - This function does **not** import modules; it only inspects ``globals()``.
      Make sure the target functions are imported into the worker/module namespace:
      ``read_anchored_align``, ``read_anchored_align_pairwise2``, ``read_anchored_align_edlib``.
    - For process pools, prefer passing **string tags** so specs are picklable.
    """
    if tag is None:
        return None
    if callable(tag):
        return tag

    t = str(tag).lower().strip()

    # Known symbols in module globals
    pairwise_first = cast(
        Optional[Callable[..., Any]],
        globals().get("read_anchored_align"),  # our PairwiseAligner-first wrapper
    )
    pairwise2_impl = cast(
        Optional[Callable[..., Any]],
        globals().get(
            "read_anchored_align_pairwise2"
        ),  # back-compat shim → same wrapper
    )
    edlib_impl = cast(
        Optional[Callable[..., Any]], globals().get("read_anchored_align_edlib")
    )

    if t == "auto":
        # prefer the modern Pairwise-first wrapper, then the shim, then edlib
        return pairwise_first or pairwise2_impl or edlib_impl
    if t in ("pairwise", "pairwise2", "biopython", "pw2"):
        return pairwise2_impl or pairwise_first
    if t in ("ed", "edlib"):
        return edlib_impl

    # last resort: try to resolve an arbitrary symbol name from globals
    return cast(Optional[Callable[..., Any]], globals().get(tag))


def _normalize_aligner_kwargs(kwargs: dict[str, Any]) -> dict[str, Any]:
    """
    Ensure kwargs carry a *callable* under 'read_align_impl' for the child process.

    - If 'read_align_impl' is missing, but 'read_aligner' is present (string or callable),
      resolve it to a callable and store under 'read_align_impl'.
    - Keep 'read_aligner' in kwargs for logging/debug if you want; it’s ignored by
      the low-level aligners.
    """
    out = dict(kwargs)
    if "read_align_impl" not in out:
        impl = _resolve_read_aligner_tag(out.get("read_aligner"))
        if impl is not None:
            out["read_align_impl"] = impl
    # optional fallback tag → callable (only if your pipeline uses it)
    if "fallback_read_aligner" in out and "fallback_read_align_impl" not in out:
        fb_impl = _resolve_read_aligner_tag(out.get("fallback_read_aligner"))
        if fb_impl is not None:
            out["fallback_read_align_impl"] = fb_impl
    return out


def _parallel_align_one_read_worker(
    read: Any,
    plasmids: Sequence[Any],
    *,
    k: int,
    margin: int,
    step: int,
    top_k: int,
    try_both_strands: bool,
    pid_mode: Literal["read", "feature", "feature_len"],
    match: float,
    mismatch: float,
    gap_open: float,
    gap_extend: float,
    read_aligner: Union[str, Callable[..., Any], None],
    fallback_read_aligner: Union[str, Callable[..., Any], None],
    # --- NEW: accept legacy/extra kwargs safely (ignored unless wired in)
    within: Optional[float] = None,
    min_jaccard: Optional[float] = None,
    fallback_kwargs: Optional[dict[str, Any]] = None,
    **_ignored: Any,
) -> Optional[ResultRow]:
    """
    Worker: seed → select candidates → align → build a ResultRow with QA fields.

    Parameters
    ----------
    read : Any
        Read object exposing `.name` and `.seq`.
    plasmids : sequence
        Prepared references exposing `.concat_ref`, `.feature_map`, `.k_index`,
        `.file`, `.plasmid_id`, `.construct`.
    k, margin, step, top_k : int
        Seeding parameters.
    try_both_strands : bool
        Evaluate both orientations in alignment.
    pid_mode : {'read','feature','feature_len'}
        PID definition returned by `align_read_to_concat`.
    match, mismatch, gap_open, gap_extend : float
        Pairwise scoring knobs (ignored by edlib paths).
    read_aligner, fallback_read_aligner : str | callable | None
        Low-level aligner(s) selection; strings resolved in globals().

    Notes
    -----
    - Accepts legacy kwargs (e.g., `within`, `fallback_kwargs`) for
      backward compatibility. They are ignored here unless you wire them in.
    """

    # --- resolve aligner tags to callables (kept exactly like before) ----------
    def _resolve(
        spec: Union[str, Callable[..., Any], None]
    ) -> Optional[Callable[..., Any]]:
        if spec is None:
            return (
                globals().get("read_anchored_align_edlib")
                or globals().get("read_anchored_align_pairwise2")
                or globals().get("read_anchored_align")
            )
        if callable(spec):
            return spec
        name = str(spec).lower()
        if name == "auto":
            return (
                globals().get("read_anchored_align")
                or globals().get("read_anchored_align_pairwise2")
                or globals().get("read_anchored_align_edlib")
            )
        if name in ("ed", "edlib"):
            return globals().get("read_anchored_align_edlib")
        if name in ("pw2", "pairwise2", "biopython", "pairwise"):
            return globals().get("read_anchored_align") or globals().get(
                "read_anchored_align_pairwise2"
            )
        return globals().get(spec)

    primary_impl = _resolve(read_aligner)
    fallback_impl = _resolve(fallback_read_aligner)

    # --- candidate selection via seeding (kept local like before) --------------
    def _topk_candidates_for_read(
        rd: Any,
        plasmids_: Sequence[Any],
        *,
        k: int,
        margin: int,
        step: int,
        top_k: int,
    ) -> list[tuple[int, Any, tuple[int, int]]]:
        cands: list[tuple[int, Any, tuple[int, int]]] = []
        for pref in plasmids_:
            win = seed_window(
                pref.concat_ref, rd.seq, pref.k_index, k=k, margin=margin, step=step
            )
            if win is None:
                continue
            s, e, hits = win
            cands.append((hits, pref, (s, e)))
        cands.sort(key=lambda x: x[0], reverse=True)
        return cands[:top_k]

    cands = _topk_candidates_for_read(
        read, plasmids, k=k, margin=margin, step=step, top_k=top_k
    )
    if not cands:
        return None

    best_row: Optional[ResultRow] = None
    best_pid = -1.0
    best_score = float("-inf")

    for _hits, pref, _win in cands:
        res = align_read_to_concat(
            concat_ref=pref.concat_ref,
            feature_map=pref.feature_map,
            read_seq=read.seq,
            k_index=pref.k_index,
            k=k,
            margin=margin,
            step=step,
            try_both_strands=try_both_strands,
            match=match,
            mismatch=mismatch,
            gap_open=gap_open,
            gap_extend=gap_extend,
            pid_mode=pid_mode,
            read_align_impl=primary_impl,
        )
        if (
            res is None
            and (fallback_impl is not None)
            and (fallback_impl is not primary_impl)
        ):
            res = align_read_to_concat(
                concat_ref=pref.concat_ref,
                feature_map=pref.feature_map,
                read_seq=read.seq,
                k_index=pref.k_index,
                k=k,
                margin=margin,
                step=step,
                try_both_strands=try_both_strands,
                match=match,
                mismatch=mismatch,
                gap_open=gap_open,
                gap_extend=gap_extend,
                pid_mode=pid_mode,
                read_align_impl=fallback_impl,
            )
        if res is None:
            continue

        # Pretty per-feature PID text (your existing function)
        feats_text = _format_feature_hits_from_alignment(
            pref,
            res.start_ref,
            res.end_ref,
            read.seq,
            res.strand,
            match=match,
            mismatch=mismatch,
            gap_open=gap_open,
            gap_extend=gap_extend,
        )

        # QA metrics (unchanged)
        min_feat_pid, rbs_pid, n_below_95 = feature_quality_metrics_from_alignment(
            pref,
            res.start_ref,
            res.end_ref,
            read.seq,
            res.strand,
            match=match,
            mismatch=mismatch,
            gap_open=gap_open,
            gap_extend=gap_extend,
        )
        suspicious = (float(res.pid) >= 98.0) and (float(min_feat_pid) < 95.0)

        row = ResultRow(
            sequence_name=str(read.name),
            read_len=len(read.seq),
            plasmid_file=str(pref.file),
            plasmid_id=str(pref.plasmid_id),
            construct=str(pref.construct),
            strand=str(res.strand),
            pid=float(res.pid),
            core_len=int(res.core_len),
            score=float(res.score),
            start_ref=int(res.start_ref),
            end_ref=int(res.end_ref),
            snps=int(len(res.snps)),
            features_hit=feats_text,
            concat_pos=f"{int(res.start_ref)}–{int(res.end_ref)}",
            min_feature_pid=float(min_feat_pid),
            rbs_pid=float(rbs_pid),
            n_features_below_95=int(n_below_95),
            suspicious=bool(suspicious),
        )

        if (row.pid > best_pid) or (row.pid == best_pid and row.score > best_score):
            best_row, best_pid, best_score = row, row.pid, row.score

    return best_row


def align_batch_parallel(
    plasmids: Sequence[PrefLike],
    reads: Sequence[Any],
    *,
    backend: Literal["threads", "processes"] = "threads",
    workers: Optional[int] = None,
    # Progress / logging
    verbose: bool = True,
    progress: bool = True,
    progress_print_every: Optional[int] = None,  # kept for API compatibility (unused)
    # Timeout & retries
    per_read_timeout_s: Optional[float] = 20.0,
    max_retries: int = 1,
    # Worker/aligner knobs (forwarded 1:1 to the worker)
    k: int = 11,
    step: int = 2,
    margin: int = 200,
    top_k: Optional[int] = 4,
    within: Optional[float] = 0.92,
    try_both_strands: bool = True,
    pid_mode: Literal["read", "feature", "feature_len"] = "read",
    match: float = 1.0,
    mismatch: float = -1.0,
    gap_open: float = -1.5,
    gap_extend: float = -0.5,
    read_aligner: Union[str, Callable[..., Any], None] = "auto",
    fallback_read_aligner: Union[str, Callable[..., Any], None] = None,
    min_jaccard: float = 0.10,
    debug: bool = False,
    # API compatibility: override kwargs for retries/timeouts
    fallback_kwargs: Optional[dict[str, Any]] = None,
) -> list[ResultRow]:
    """
    Run per-read alignments against plasmids in parallel and return best hit per read.

    Each read is processed by `_parallel_align_one_read_worker(read, plasmids, **kwargs)`,
    which performs seeding (top-K) and read-anchored alignment. Arguments such as
    `within`, `min_jaccard`, scoring, and aligner selection are forwarded unchanged
    to the worker. Failed/long-running tasks can be retried with `fallback_kwargs`.

    Parameters
    ----------
    plasmids : sequence of PrefLike
        Prepared plasmid references (must provide `.concat_ref`, `.feature_map`,
        `.k_index`, `.file`, `.plasmid_id`, `.construct`).
    reads : sequence
        Read objects; the worker uses at least `.name` and `.seq`.
    backend : {'threads', 'processes'}, default 'threads'
        Parallel backend to use.
    workers : int or None, optional
        Number of workers. If None/<=0: 4 for threads, or `max(1, cpu_count()-1)` for processes.
    verbose : bool, default True
        Print high-level status and retry messages.
    progress : bool, default True
        Print lightweight progress counters.
    progress_print_every : int or None, optional
        Kept for API compatibility (not used).
    per_read_timeout_s : float or None, default 20.0
        If not None, reads that take longer than this (and return None) are retried.
    max_retries : int, default 1
        Number of retry attempts per read on exception/timeout/None-result.
    k, step, margin, top_k : int
        Seeding and windowing parameters forwarded to the worker.
    within : float or None, default 0.92
        Candidate gating ratio in the worker; None disables gating.
    try_both_strands : bool, default True
        Allow worker to consider both orientations.
    pid_mode : {'read', 'feature', 'feature_len'}, default 'read'
        PID definition the worker should compute.
    match, mismatch, gap_open, gap_extend : float
        Scoring parameters forwarded to the worker (ignored by edlib).
    read_aligner : {'auto','edlib','pairwise2'} or callable or None, default 'auto'
        Primary low-level aligner selector for the worker.
    fallback_read_aligner : same as `read_aligner`, optional
        Secondary aligner the worker may try if primary returns None.
    min_jaccard : float, default 0.10
        Minimum Jaccard threshold used by the worker prefilter.
    debug : bool, default False
        Extra diagnostics printed on failures/retries.
    fallback_kwargs : dict or None, optional
        Extra kwargs merged into worker kwargs **only on retries** (e.g., relax thresholds).

    Returns
    -------
    list of ResultRow
        One best row per successfully aligned read. Reads failing all retries are omitted.

    Notes
    -----
    Requires `_parallel_align_one_read_worker(read, plasmids, **kwargs)` to accept:
    `k, step, margin, top_k, within, try_both_strands, pid_mode, match, mismatch,
    gap_open, gap_extend, read_aligner, fallback_read_aligner, min_jaccard, debug`.
    """
    import concurrent.futures as _fut

    # Resolve worker count
    if workers is None or workers <= 0:
        if backend == "threads":
            workers = 4
        else:
            cpu = os.cpu_count() or 4
            workers = max(1, cpu - 1)

    # Base kwargs forwarded to the worker
    worker_kwargs: dict[str, Any] = {
        "k": int(k),
        "step": int(step),
        "margin": int(margin),
        "top_k": (None if top_k is None else int(top_k)),
        "within": (None if within is None else float(within)),
        "try_both_strands": bool(try_both_strands),
        "pid_mode": pid_mode,  # Literal is fine to forward as-is
        "match": float(match),
        "mismatch": float(mismatch),
        "gap_open": float(gap_open),
        "gap_extend": float(gap_extend),
        "read_aligner": read_aligner,
        "fallback_read_aligner": fallback_read_aligner,
        "min_jaccard": float(min_jaccard),
        "debug": bool(debug),
    }

    rows: list[ResultRow] = []

    Executor = (  # noqa: N806
        _fut.ThreadPoolExecutor if backend == "threads" else _fut.ProcessPoolExecutor
    )
    with Executor(max_workers=workers) as ex:
        pending: set[_fut.Future[Optional[ResultRow]]] = set()
        meta: dict[_fut.Future[Optional[ResultRow]], dict[str, Any]] = {}

        # Submit initial tasks
        for rd in reads:
            fut = ex.submit(
                _parallel_align_one_read_worker, rd, plasmids, **worker_kwargs
            )
            pending.add(fut)
            meta[fut] = {"t0": time.perf_counter(), "retries": 0, "read": rd}

        done_n = 0
        last_print = time.perf_counter()

        while pending:
            done, pending = _fut.wait(
                pending, timeout=0.25, return_when=_fut.FIRST_COMPLETED
            )

            # Lightweight progress
            if progress:
                now = time.perf_counter()
                if now - last_print > 1.0:
                    if verbose:
                        print(f"[align] finished={done_n} pending={len(pending)}")
                    last_print = now

            for fut in done:
                info = meta.pop(
                    fut, {"t0": time.perf_counter(), "retries": 0, "read": None}
                )
                rd = info.get("read")
                retries = int(info.get("retries", 0))
                t0 = float(info.get("t0", time.perf_counter()))

                try:
                    res = fut.result(timeout=0)  # already finished
                except Exception as e:
                    # Exception → retry
                    if retries < int(max_retries) and rd is not None:
                        if verbose:
                            print(
                                f"[align] read '{getattr(rd, 'name', '?')}' failed → retry {retries+1}/{max_retries}"
                            )
                            if debug:
                                traceback.print_exc()
                        kw = dict(worker_kwargs)
                        if fallback_kwargs:
                            kw.update(fallback_kwargs)
                        fut2 = ex.submit(
                            _parallel_align_one_read_worker, rd, plasmids, **kw
                        )
                        pending.add(fut2)
                        meta[fut2] = {
                            "t0": time.perf_counter(),
                            "retries": retries + 1,
                            "read": rd,
                        }
                        continue
                    else:
                        if verbose:
                            print(
                                f"[align] read '{getattr(rd, 'name', '?')}' permanently failed: {e}"
                            )
                        res = None

                # Timeout-based retry: result is None and took too long
                took = time.perf_counter() - t0
                if (
                    res is None
                    and per_read_timeout_s is not None
                    and took > float(per_read_timeout_s)
                    and retries < int(max_retries)
                    and rd is not None
                ):
                    if verbose:
                        print(
                            f"[align] read '{getattr(rd, 'name', '?')}' timeout → retry {retries+1}/{max_retries}"
                        )
                    kw = dict(worker_kwargs)
                    if fallback_kwargs:
                        kw.update(fallback_kwargs)
                    fut2 = ex.submit(
                        _parallel_align_one_read_worker, rd, plasmids, **kw
                    )
                    pending.add(fut2)
                    meta[fut2] = {
                        "t0": time.perf_counter(),
                        "retries": retries + 1,
                        "read": rd,
                    }
                    continue

                if res is not None:
                    rows.append(res)
                done_n += 1

    return rows


def summarize_matches(rows: Sequence[ResultRow], min_pid: float = 85.0) -> pd.DataFrame:
    """
    Summarize per-plasmid matches: number of reads, best PID, and mean PID.

    Parameters
    ----------
    rows : Sequence[ResultRow]
        Best-hit rows (e.g., from `align_batch`, `align_batch_parallel`,
        or `dealign_reads_easy`).
    min_pid : float, default 85.0
        Only count matches with PID ≥ `min_pid`.

    Returns
    -------
    pandas.DataFrame
        One row per plasmid with columns:
        - ``plasmid_id``
        - ``construct``
        - ``plasmid_file``
        - ``n_reads``    (unique read names counted)
        - ``best_pid``   (max PID)
        - ``mean_pid``   (mean PID)

    Notes
    -----
    - If there are no rows or no rows meeting the threshold, an empty
      DataFrame with the expected columns is returned.
    - ``n_reads`` counts unique read names (``sequence_name``) per plasmid.

    Examples
    --------
    >>> summary = summarize_matches(rows, min_pid=90.0)
    >>> summary.head()
    """
    columns = [
        "plasmid_id",
        "construct",
        "plasmid_file",
        "n_reads",
        "best_pid",
        "mean_pid",
    ]
    if not rows:
        return pd.DataFrame(columns=columns)

    df = pd.DataFrame(
        {
            "plasmid_id": [r.plasmid_id for r in rows],
            "construct": [r.construct for r in rows],
            "plasmid_file": [r.plasmid_file for r in rows],
            "sequence_name": [r.sequence_name for r in rows],
            "pid": [float(r.pid) for r in rows],
        }
    )

    df = df[df["pid"] >= float(min_pid)]
    if df.empty:
        return pd.DataFrame(columns=columns)

    out = (
        df.groupby(["plasmid_id", "construct", "plasmid_file"], as_index=False)
        .agg(
            n_reads=("sequence_name", "nunique"),
            best_pid=("pid", "max"),
            mean_pid=("pid", "mean"),
        )
        .sort_values(["construct", "plasmid_id"], kind="stable")
    )
    # Spaltenreihenfolge sicherstellen
    return out[columns]


def find_unbuilt_plasmids(
    plasmids: Sequence[PlasmidRef],
    rows: Sequence[ResultRow],
    min_pid: float = 85.0,
) -> pd.DataFrame:
    """
    List plasmids with **no** matching read at or above the given PID threshold.

    Parameters
    ----------
    plasmids : Sequence[PlasmidRef]
        Prepared references returned by `prepare_plasmids(...)`.
    rows : Sequence[ResultRow]
        Best-hit rows (e.g., from `align_batch*` / `dealign_reads_easy`).
    min_pid : float, default 85.0
        Threshold: a plasmid is considered "built" if any row for it has
        PID ≥ `min_pid`.

    Returns
    -------
    pandas.DataFrame
        One row per plasmid **without** any qualifying match, with columns:
        - ``plasmid_id``
        - ``construct``
        - ``plasmid_file``

    Examples
    --------
    >>> missing = find_unbuilt_plasmids(plasmids, rows, min_pid=90.0)
    >>> missing
    """
    threshold = float(min_pid)
    matched_ids = {r.plasmid_id for r in rows if float(r.pid) >= threshold}
    missing = [p for p in plasmids if p.plasmid_id not in matched_ids]

    if not missing:
        return pd.DataFrame(columns=["plasmid_id", "construct", "plasmid_file"])

    out = pd.DataFrame(
        {
            "plasmid_id": [p.plasmid_id for p in missing],
            "construct": [p.construct for p in missing],
            "plasmid_file": [p.file for p in missing],
        }
    ).sort_values(["construct", "plasmid_id"], kind="stable")

    return out[["plasmid_id", "construct", "plasmid_file"]]


def results_to_dataframe(rows: Sequence[ResultRow]):
    """
    Convert a sequence of ResultRow objects to a pandas DataFrame.

    Works with both regular and slotted dataclasses (i.e., no __dict__).
    Falls back to `vars()` for non-dataclass records.

    Parameters
    ----------
    rows : sequence of ResultRow
        Result records to convert.

    Returns
    -------
    pandas.DataFrame
        A DataFrame containing one row per result. Columns are ordered using
        the annotations of `ResultRow` when available; unexpected/extra keys
        are appended to the end in stable order.

    Notes
    -----
    - If `rows` is empty, returns an empty DataFrame with the annotated
      columns of `ResultRow`.
    - Uses `dataclasses.asdict` for dataclass instances to support `slots=True`.
    """
    from dataclasses import asdict, is_dataclass
    from dataclasses import fields as dc_fields

    # Preferred column order taken from the dataclass annotations (if present)
    try:
        base_cols: list[str] = list(ResultRow.__annotations__.keys())
    except Exception:
        # Fallback: derive from dataclass fields, or leave empty
        try:
            base_cols = [f.name for f in dc_fields(ResultRow)]
        except Exception:
            base_cols = []

    if not rows:
        return pd.DataFrame(columns=base_cols)

    recs: list[dict[str, Any]] = []
    extra_cols_order: list[str] = []

    for r in rows:
        if is_dataclass(r):
            d: dict[str, Any] = asdict(r)
        else:
            # Non-dataclass fallback (may still fail if truly slot-only non-dataclass)
            try:
                d = dict(vars(r))
            except TypeError:
                # Last resort: attempt attribute extraction using annotated fields
                d = {k: getattr(r, k) for k in base_cols if hasattr(r, k)}

        # Track extra keys not in the canonical set
        for k in d.keys():
            if k not in base_cols and k not in extra_cols_order:
                extra_cols_order.append(k)

        recs.append(d)

    cols = base_cols + extra_cols_order
    df = pd.DataFrame.from_records(recs)

    # Reorder columns; keep any missing keys safely handled
    df = df.reindex(columns=cols)

    return df


def _pick_xlsx_engine() -> str:
    """Return a working Excel writer engine name or raise if none is available."""
    try:
        import xlsxwriter  # noqa: F401

        return "xlsxwriter"
    except Exception:
        pass
    try:
        import openpyxl  # noqa: F401

        return "openpyxl"
    except Exception as exc:  # pylint: disable=broad-except
        raise RuntimeError(
            "No Excel engine found. Install 'xlsxwriter' or 'openpyxl' for XLSX export."
        ) from exc


def _autosize_columns(df: pd.DataFrame, worksheet, workbook_engine: str) -> None:
    """
    Best-effort column autosize for Excel.

    Parameters
    ----------
    df : pd.DataFrame
        Data to size against.
    worksheet : Any
        Worksheet object from the engine.
    workbook_engine : {"xlsxwriter", "openpyxl"}
        Selected engine.
    """
    # Compute naive width: max(len(header), len(value))
    widths: list[int] = []
    for col in df.columns:
        max_len = len(str(col))
        for val in df[col].astype(str).values:
            if len(val) > max_len:
                max_len = len(val)
        # Add a little padding
        widths.append(min(max_len + 2, 80))

    if workbook_engine == "xlsxwriter":
        for i, w in enumerate(widths):
            worksheet.set_column(i, i, float(w))
    else:
        # openpyxl: has a 'width' attribute in column_dimensions, but it's approximate.
        from openpyxl.utils import get_column_letter

        for i, w in enumerate(widths, start=1):
            worksheet.column_dimensions[get_column_letter(i)].width = float(w)


def export_results(
    rows: Sequence[ResultRow],
    out_path: str,
    *,
    fmt: Literal["csv", "xlsx"] = "xlsx",
    include_sheet_name: str = "Summary",
    text_wrap_columns: tuple[str, ...] = ("features_hit",),
) -> str:
    """
    Export batch alignment results to CSV or XLSX.

    Parameters
    ----------
    rows : sequence of ResultRow
        Rows produced by `align_batch`.
    out_path : str
        Target file path. Extension is not required; `fmt` decides.
    fmt : {"csv", "xlsx"}, optional
        Output format, by default "xlsx".
    include_sheet_name : str, optional
        Sheet name for the results in XLSX, by default "Summary".
    text_wrap_columns : tuple of str, optional
        Columns that should be wrapped in Excel, by default ("features_hit",).

    Returns
    -------
    str
        Absolute path of the written file.

    Notes
    -----
    - PID is written as a percentage in Excel (e.g., 97.5%).
    - If `rows` is empty, a minimal header-only sheet/file is written.
    """
    if not out_path:
        raise ValueError("out_path must not be empty")
    if fmt not in ("csv", "xlsx"):
        raise ValueError("fmt must be 'csv' or 'xlsx'")

    df = results_to_dataframe(rows)

    # Ensure stable column order (matching ResultRow)
    cols = list(getattr(ResultRow, "__annotations__", {}).keys())
    if df.empty:
        df = pd.DataFrame(columns=cols)
    else:
        df = df.loc[:, cols]

    # Pick final path
    base, ext = os.path.splitext(out_path)
    if fmt == "csv":
        final_path = base + ".csv" if ext.lower() != ".csv" else out_path
        df.to_csv(final_path, index=False)
        return os.path.abspath(final_path)

    # XLSX
    engine = _pick_xlsx_engine()
    final_path = base + ".xlsx" if ext.lower() != ".xlsx" else out_path

    with pd.ExcelWriter(
        final_path, engine=engine, datetime_format="YYYY-MM-DD HH:MM:SS"
    ) as xw:
        df.to_excel(xw, sheet_name=include_sheet_name, index=False)
        workbook = xw.book
        worksheet = xw.sheets[include_sheet_name]

        # Formats
        if engine == "xlsxwriter":
            pid_fmt = workbook.add_format({"num_format": "0.00%"})

            int_fmt = workbook.add_format({"num_format": "0"})

            wrap_fmt = workbook.add_format({"text_wrap": True})

        else:
            # openpyxl formats are applied differently via styles; use pandas Styler as fallback:
            pid_fmt = None
            int_fmt = None
            wrap_fmt = None

        # Locate columns
        col_index: dict[str, int] = {name: i for i, name in enumerate(df.columns)}

        # Apply number formats (xlsxwriter path; for openpyxl we at least autosize)
        if engine == "xlsxwriter":
            if "pid" in col_index:
                j = col_index["pid"]
                worksheet.set_column(j, j, None, pid_fmt)
                # Convert 0..100 to 0..1 for Excel percentage format
                if not df.empty:
                    pct = df["pid"].astype(float) / 100.0
                    df2 = df.copy()
                    df2["pid"] = pct
                    worksheet.write_column(1, j, df2["pid"].tolist())
            for int_col in ("core_len", "snps", "start_ref", "end_ref", "read_len"):
                if int_col in col_index:
                    j = col_index[int_col]
                    worksheet.set_column(j, j, None, int_fmt)
            for name in text_wrap_columns:
                if name in col_index:
                    j = col_index[name]
                    worksheet.set_column(j, j, None, wrap_fmt)

        # Freeze header
        worksheet.freeze_panes(1, 0)
        # Auto filter
        worksheet.autofilter(0, 0, max(1, len(df)), max(0, len(df.columns) - 1))

        # Autosize columns (best-effort)
        _autosize_columns(df, worksheet, engine)

    return os.path.abspath(final_path)


def export_feature_inventory(
    plasmids: Sequence[PlasmidRef],
    out_path: str,
    *,
    fmt: Literal["csv", "xlsx"] = "xlsx",
    sheet_name: str = "Inventory",
) -> str:
    """
    Export a slot-wise unique-parts inventory over prepared plasmids.

    Parameters
    ----------
    plasmids : sequence of PlasmidRef
        Prepared references from `prepare_plasmids`.
    out_path : str
        Target file path. Extension is decided by `fmt`.
    fmt : {"csv", "xlsx"}, optional
        Output format, by default "xlsx".
    sheet_name : str, optional
        Worksheet name for XLSX, by default "Inventory".

    Returns
    -------
    str
        Absolute path of the written file.
    """
    # Collect unique tokens per slot
    slots = ("promoter", "rbs", "gene", "term")
    seen: dict[str, dict[str, None]] = {s: {} for s in slots}

    # Reconstruct SeqRecords is not needed; PlasmidRef already has feature_map and construct,
    # but we want the *tokens per slot* across records; we can derive it from constructs.
    # Fast path: parse constructs; fallback to empty lists.

    for pref in plasmids:
        tokens = [t.strip() for t in pref.construct.split("_") if t.strip()]
        for s, tok in zip(slots, tokens, strict=False):
            if tok:
                seen[s][tok] = None

    # <-- patched: use Series so columns may have different lengths; fill blanks
    inv_df = pd.DataFrame(
        {slot: pd.Series(sorted(seen[slot].keys(), key=str.lower)) for slot in slots}
    ).fillna("")

    base, ext = os.path.splitext(out_path)
    if fmt == "csv":
        final = base + ".csv" if ext.lower() != ".csv" else out_path
        inv_df.to_csv(final, index=False)
        return os.path.abspath(final)

    engine = _pick_xlsx_engine()
    final = base + ".xlsx" if ext.lower() != ".xlsx" else out_path
    with pd.ExcelWriter(final, engine=engine) as xw:
        inv_df.to_excel(xw, sheet_name=sheet_name, index=False)
        ws = xw.sheets[sheet_name]
        _autosize_columns(inv_df, ws, engine)
        ws.freeze_panes(1, 0)
    return os.path.abspath(final)
