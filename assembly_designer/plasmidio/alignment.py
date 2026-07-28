"""Core read-anchored alignment engine (k-mer seeding + per-segment aligners).

Leaf module: does not depend on any other plasmidio submodule except
:mod:`.lexicon` (for ``revcomp``).
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any, Literal, Optional

from .lexicon import revcomp

KmerIndex = dict[str, list[int]]


@dataclass(frozen=True)
class SNP:
    """
    Single-nucleotide difference within the core alignment.

    Attributes
    ----------
    ref_pos : int
        Position on the reference concat-band (0-based).
    read_pos : int
        Position on the read (0-based, after trimming performed by caller).
    ref_base : str
        Base on the reference.
    read_base : str
        Base on the read.
    """

    ref_pos: int
    read_pos: int
    ref_base: str
    read_base: str


def _format_feature_hits(res: AlignmentResult) -> str:
    """
    Convert `res.feature_hits` into a compact, human-readable string.

    Examples
    --------
    >>> # FeatureHit(ftype='promoter', fname='J23106', coverage_pct=100)
    >>> _format_feature_hits(res)
    "promoter:J23106 (100%), 5'UTR:B0033m (100%)"

    Notes
    -----
    - Returns an empty string if there are no feature hits.
    - Robust against both dataclass objects and dict-like items.
    """
    hits = getattr(res, "feature_hits", None)
    if not hits:
        return ""

    parts: list[str] = []
    for h in hits:
        # Accept dataclass FeatureHit or dict-like records
        try:
            ftype = getattr(h, "ftype", None) or (
                h.get("ftype") if isinstance(h, dict) else "feat"
            )
            fname = getattr(h, "fname", None) or (
                h.get("fname") if isinstance(h, dict) else ""
            )
            pct = getattr(h, "coverage_pct", None) or (
                h.get("coverage_pct") if isinstance(h, dict) else 0
            )
            pct_i = int(round(float(pct))) if pct is not None else 0
            parts.append(f"{ftype}:{fname} ({pct_i}%)")
        except Exception:
            # Never break the whole render on a single odd item
            parts.append(str(h))

    return ", ".join(parts)


@dataclass(frozen=True)
class FeatureHit:
    """
    Coverage summary for a single feature block on the concat-band.

    Attributes
    ----------
    start : int
        Block start on concat-band (0-based, inclusive).
    end : int
        Block end on concat-band (0-based, exclusive).
    ftype : str
        GenBank feature type (e.g., "promoter", "RBS", "CDS", "terminator", ...).
    fname : str
        Display name (e.g., "J23100", "BCD12", "EcPanD", "B0015").
    covered_bp : int
        Overlap length between aligned reference span and this block.
    coverage_pct : float
        100 * covered_bp / (end - start).
    """

    start: int
    end: int
    ftype: str
    fname: str
    covered_bp: int
    coverage_pct: float


@dataclass(frozen=True)
class AlignmentResult:
    """
    Result of read-anchored alignment against a concat-band segment.

    Attributes
    ----------
    pid : float
        Percent identity in the core, **denominator = read bases in core**.
    core_len : int
        Number of read bases contributing to PID (≈ aligned read length).
    score : float
        Raw alignment score from Biopython (for tie-breaking).
    strand : str
        'F' for forward, 'R' for reverse-complement chosen.
    start_ref : int
        Start coordinate on the **full** concat-band (0-based, inclusive).
    end_ref : int
        End coordinate on the **full** concat-band (0-based, exclusive).
    start_read : int
        First consumed base on the read within the alignment (0-based).
    end_read : int
        Past-the-end index on the read (0-based).
    snps : list[SNP]
        List of SNP differences (core mismatches).
    feature_hits : list[FeatureHit]
        Feature coverage summary intersecting [start_ref, end_ref).
    """

    pid: float
    core_len: int
    score: float
    strand: str
    start_ref: int
    end_ref: int
    start_read: int
    end_read: int
    snps: list[SNP]
    feature_hits: list[FeatureHit]


def build_kmer_index(seq: str, k: int = 11) -> KmerIndex:
    """
    Build a simple K-mer → positions index for a reference string.

    Parameters
    ----------
    seq : str
        Upper/lowercase accepted; it will be uppercased internally.
    k : int, optional
        K-mer length, by default 11.

    Returns
    -------
    dict[str, list[int]]
        Map from K-mer (uppercase) to a sorted list of start positions (0-based).

    Notes
    -----
    - For Sanger-length reads, k=9..13 work well; k=11 is a good default.
    """
    s = seq.upper()
    idx: KmerIndex = {}
    n = len(s) - k + 1
    if n < 1:
        return idx
    # Use a dict of lists; append in-order yields sorted positions per key.
    for i in range(n):
        key = s[i : i + k]
        idx.setdefault(key, []).append(i)
    return idx


def seed_window(
    ref_seq: str,
    read_seq: str,
    k_index: KmerIndex,
    *,
    k: int = 11,
    margin: int = 200,
    step: int = 2,
) -> Optional[tuple[int, int, int]]:
    """
    Estimate a search window on the reference using K-mer offsets (median).

    Parameters
    ----------
    ref_seq : str
        Full reference (concat-band).
    read_seq : str
        Read sequence to seed (will be uppercased internally).
    k_index : dict[str, list[int]]
        Prebuilt K-mer index for `ref_seq`.
    k : int, optional
        K-mer length, by default 11. Must match the index.
    margin : int, optional
        Extra bases left/right of the median offset, by default 200.
    step : int, optional
        Subsample the read K-mers (every `step`-th k-mer), by default 2.

    Returns
    -------
    tuple or None
        (start, end, hits) window in the reference, or None if no seeds found.

    Notes
    -----
    - We only use **one** position per read K-mer (the first in the index) to keep it fast.
    - For very short reads (<k), returns None.
    """
    r = read_seq.upper()
    if len(r) < k:
        return None

    offsets: list[int] = []
    hits = 0
    for i in range(0, len(r) - k + 1, max(1, step)):
        mer = r[i : i + k]
        pos_list = k_index.get(mer)
        if not pos_list:
            continue
        hits += 1
        offsets.append(pos_list[0] - i)

    if not offsets:
        return None

    offsets.sort()
    median_off = offsets[len(offsets) // 2]
    start = max(0, median_off - margin)
    end = min(len(ref_seq), median_off + len(r) + margin)
    if end <= start:
        return None
    return start, end, hits


def _core_metrics_and_snps(
    aligned_ref: str,
    aligned_read: str,
    *,
    read_pid_core_only: bool = False,
) -> tuple[float, int, int, int, int, int, list]:
    """
    Compute read-anchored PID, core coordinates, and SNPs from two gapped strings.

    Parameters
    ----------
    aligned_ref, aligned_read : str
        Gapped alignment strings of equal length ('-' as gap).
    read_pid_core_only : bool, optional
        If False (default), PID denominator = #columns where the **read** consumes a base
        anywhere in the alignment (penalizes read-only tails).
        If True, denominator = #core columns only (core = both sides consume).

    Returns
    -------
    pid : float
        Percent identity (0..100) using the chosen denominator.
    core_len : int
        Number of core columns (both sides consume a base).
    start_ref, end_ref : int
        Ungapped reference consumption bounds of the core (0-based, end-exclusive).
    start_read, end_read : int
        Ungapped read   consumption bounds of the core (0-based, end-exclusive).
    snps : list
        One entry per mismatch **inside the core**.
        Each entry is a dict with keys:
        ``ref_pos``, ``read_pos``, ``ref_base``, ``read_base``.
        (The calling code can convert these dicts into your `SNP` dataclass using `snp_factory(**rec)`)

    Notes
    -----
    - Coordinates are measured in consumed (ungapped) units on each side.
    - If there is no overlapping core, returns zeros and an empty list.
    """
    if len(aligned_ref) != len(aligned_read):
        raise ValueError("aligned_ref and aligned_read must have the same length.")

    n = len(aligned_ref)
    if n == 0:
        return 0.0, 0, 0, 0, 0, 0, []

    # locate first/last column where both sides consume a base (the core window)
    first_core: int | None = None
    last_core: int | None = None
    for i in range(n):
        if aligned_ref[i] != "-" and aligned_read[i] != "-":
            first_core = i
            break
    for i in range(n - 1, -1, -1):
        if aligned_ref[i] != "-" and aligned_read[i] != "-":
            last_core = i
            break

    # no overlapping core
    if first_core is None or last_core is None or last_core < first_core:
        read_consumed = sum(1 for i in range(n) if aligned_read[i] != "-")
        denom = 0 if read_pid_core_only else read_consumed
        pid = 0.0 if denom == 0 else 0.0
        return pid, 0, 0, 0, 0, 0, []

    # walk once to compute coords + SNPs
    ref_pos = 0
    read_pos = 0
    start_ref: int | None = None
    start_read: int | None = None
    end_ref = 0
    end_read = 0

    matches = 0
    core_len = 0
    snps: list[dict[str, int | str]] = []

    # pre-count read-consumed columns (for denominator when read_pid_core_only=False)
    read_consumed_all = sum(1 for i in range(n) if aligned_read[i] != "-")

    for i in range(n):
        a = aligned_ref[i]
        b = aligned_read[i]
        ref_consumes = a != "-"
        read_consumes = b != "-"

        if i == first_core and ref_consumes and read_consumes:
            start_ref = ref_pos
            start_read = read_pos

        if ref_consumes and read_consumes:
            # inside core
            core_len += 1
            if a == b:
                matches += 1
            else:
                snps.append(
                    {
                        "ref_pos": ref_pos,
                        "read_pos": read_pos,
                        "ref_base": a,
                        "read_base": b,
                    }
                )
            ref_pos += 1
            read_pos += 1
            end_ref = ref_pos
            end_read = read_pos
        elif ref_consumes:
            ref_pos += 1
            if i <= last_core:
                end_ref = ref_pos
        elif read_consumes:
            read_pos += 1
            if i <= last_core:
                end_read = read_pos

        if i == last_core:
            break

    if start_ref is None:
        start_ref = 0
    if start_read is None:
        start_read = 0

    denom = core_len if read_pid_core_only else read_consumed_all
    pid = 0.0 if denom == 0 else (100.0 * matches / float(denom))
    return (
        float(pid),
        int(core_len),
        int(start_ref),
        int(end_ref),
        int(start_read),
        int(end_read),
        snps,
    )


def read_anchored_align(
    ref_segment: str,
    read_seq: str,
    *,
    # scoring (only for pairwise2)
    match: float = 1.0,
    mismatch: float = -1.0,
    gap_open: float = -1.5,
    gap_extend: float = -0.5,
    # Feature-PID
    pid_mode: Literal["read", "feature", "feature_len"] = "read",
    feature_map: Optional[
        Sequence[tuple[int, int, str, str]]
    ] = None,  # (s,e,type,name)
    ref_abs_start: Optional[int] = None,
    # SNP-Objects
    snp_factory: Optional[Callable[..., Any]] = None,
) -> tuple[float, int, int, int, int, int, float, list[Any]]:
    """
    Align read to a reference segment (pairwise2 → edlib fallback),
    compute read-anchored PID + core coords via `_core_metrics_and_snps`,
    and optionally compute a feature-basiertes PID (pid_mode).

    Returns
    -------
    (pid, core_len, start_ref, end_ref, start_read, end_read, score, snps)
    """
    r = read_seq.upper()
    ref = ref_segment.upper()

    # ---- 1) obtain gapped Alignment ----
    alnA = alnB = None  # noqa: N806
    score = float("-inf")
    try:
        from Bio import pairwise2

        alns = pairwise2.align.globalms(
            ref,
            r,
            match,
            mismatch,
            gap_open,
            gap_extend,
            one_alignment_only=True,
            penalize_end_gaps=False,
        )
        if alns:
            a_ref, a_read, sc, *_ = alns[0]
            alnA, alnB, score = a_ref, a_read, float(sc)  # noqa: N806
    except Exception:
        alnA = alnB = None  # noqa: N806

    if alnA is None or alnB is None:
        # ---- edlib Fallback → reconstruct gapped ----
        try:
            import edlib

            ed = edlib.align(r, ref, mode="NW", task="path")
            if ed and "cigar" in ed:
                cigar = ed["cigar"]
                a_ref, a_read = [], []
                i = j = 0
                num = ""
                for ch in cigar:
                    if ch.isdigit():
                        num += ch
                        continue
                    n = int(num) if num else 1
                    num = ""
                    if ch in ("=", "M", "X"):
                        a_ref.append(ref[i : i + n])
                        a_read.append(r[j : j + n])
                        i += n
                        j += n
                    elif ch == "I":
                        a_ref.append("-" * n)
                        a_read.append(r[j : j + n])
                        j += n
                    elif ch == "D":
                        a_ref.append(ref[i : i + n])
                        a_read.append("-" * n)
                        i += n
                alnA, alnB = "".join(a_ref), "".join(a_read)  # noqa: N806
                score = float(-ed.get("editDistance", 0))
            else:
                return (0.0, 0, 0, 0, 0, 0, float("-inf"), [])
        except Exception:
            return (0.0, 0, 0, 0, 0, 0, float("-inf"), [])

    # ---- 2) Core metrics & SNPs centralized from Helper ----
    pid_read, core_len, sA, eA, sB, eB, snp_recs = _core_metrics_and_snps(alnA, alnB)  # noqa: N806  # fmt: skip

    # ---- 3) Optional: Calculate Feature-PID ----
    pid_final = pid_read
    if pid_mode != "read" and feature_map and (ref_abs_start is not None):
        # mappe gapped REF-Columns → absolute REF-coordinates
        abs_pos: list[Optional[int]] = []
        rp = ref_abs_start
        for ch in alnA:
            if ch == "-":
                abs_pos.append(None)
            else:
                abs_pos.append(rp)
                rp += 1

        # Determine core mask on gapped indices
        # (We only know sA..eA is ungapped; run it again and keep count)
        ref_pos = 0
        core_idx: list[int] = []
        for i, a in enumerate(alnA):
            if a != "-":
                if sA <= ref_pos < eA:
                    core_idx.append(i)
                ref_pos += 1

        def _in_feature(p: int) -> bool:
            for fs, fe, *_ in feature_map:
                if fs <= p < fe:
                    return True
            return False

        matches_core = 0
        denom_core = 0  # genuine core feature columns
        for i in core_idx:
            ap = abs_pos[i]
            if ap is None:
                continue
            a = alnA[i]
            b = alnB[i]
            if a != "-" and b != "-" and _in_feature(ap):
                denom_core += 1
                if a == b:
                    matches_core += 1

        # Length-Denominator (Overlap of features ∩ Core in REF-Coordinates)
        denom_len = 0
        abs_s = ref_abs_start + sA
        abs_e = ref_abs_start + eA
        for fs, fe, *_ in feature_map:
            left = max(abs_s, fs)
            right = min(abs_e, fe)
            if right > left:
                denom_len += right - left

        pid_feat_core = (
            (100.0 * matches_core / float(denom_core)) if denom_core else 0.0
        )
        pid_feat_len = (100.0 * matches_core / float(denom_len)) if denom_len else 0.0
        pid_final = pid_feat_core if pid_mode == "feature" else pid_feat_len

    # ---- 4) Create SNP-Objects (if factory) ----
    if snp_factory is not None:
        if snp_recs:
            first = snp_recs[0]
            if hasattr(first, "__dict__"):
                snps = [
                    snp_factory(
                        ref_pos=s.ref_pos,
                        read_pos=s.read_pos,
                        ref_base=s.ref_base,
                        read_base=s.read_base,
                    )
                    for s in snp_recs
                ]
            else:
                snps = [snp_factory(**rec) for rec in snp_recs]
        else:
            snps = []
    else:
        snps = list(snp_recs)

    return (
        float(pid_final),
        int(core_len),
        int(sA),
        int(eA),
        int(sB),
        int(eB),
        float(score),
        snps,
    )


def summarize_feature_coverage(
    feature_map: Sequence[tuple[int, int, str, str]],
    ref_span: tuple[int, int],
) -> list[FeatureHit]:
    """
    Summarize coverage per feature block given an aligned reference span.

    Parameters
    ----------
    feature_map : sequence of (start, end, ftype, fname)
        Mapping returned by `build_feature_concat_sequence`.
    ref_span : tuple of int
        (start_ref, end_ref) on concat-band (0-based, exclusive end).

    Returns
    -------
    list[FeatureHit]
        Coverage report for blocks overlapping the span (sorted by start).
    """
    start_ref, end_ref = ref_span
    hits: list[FeatureHit] = []
    for blk_start, blk_end, ftype, fname in feature_map:
        # Overlap with [start_ref, end_ref)
        left = max(start_ref, blk_start)
        right = min(end_ref, blk_end)
        if right <= left:
            continue
        covered = right - left
        pct = 100.0 * covered / float(blk_end - blk_start)
        hits.append(
            FeatureHit(
                start=blk_start,
                end=blk_end,
                ftype=str(ftype),
                fname=str(fname),
                covered_bp=covered,
                coverage_pct=pct,
            )
        )
    hits.sort(key=lambda h: (h.start, h.end))
    return hits


def feature_pid_from_alignment(
    aln_ref: str,
    aln_read: str,
    *,
    ref_abs_start: int,
    feature_map: list[tuple[int, int, str, str]],
    ref_abs_end: int | None = None,
) -> tuple[float, float]:
    """
    Compute PID restricted to feature regions (two flavors).

    Parameters
    ----------
    aln_ref, aln_read
        Gapped alignment strings (same length), '-' as gap.
    ref_abs_start
        Absolute concat-band start (0-based) of the aligned reference segment.
    feature_map
        List of (start, end, ftype, fname) blocks on the concat-band.
    ref_abs_end
        Absolute end (exclusive) of the aligned reference segment. If None,
        inferred from `aln_ref` length and gaps.

    Returns
    -------
    pid_feature_core : float
        Matches / aligned columns that are inside any feature (both sides consume).
    pid_feature_len : float
        Matches / reference feature length intersecting [ref_abs_start, ref_abs_end).

    Notes
    -----
    - This does not inspect read-only tails; it purely evaluates feature overlap.
    """
    # Map alignment columns to absolute ref positions (or None for ref gaps).
    ref_pos: list[int | None] = []
    rp = ref_abs_start
    for ch in aln_ref:
        if ch == "-":
            ref_pos.append(None)
        else:
            ref_pos.append(rp)
            rp += 1

    if ref_abs_end is None:
        ref_abs_end = ref_abs_start + sum(1 for ch in aln_ref if ch != "-")

    def in_feat(p: int) -> bool:
        for s, e, *_ in feature_map:
            if s <= p < e:
                return True
        return False

    matches_core = 0
    denom_core = 0
    for i, (a, b) in enumerate(zip(aln_ref, aln_read, strict=False)):
        p = ref_pos[i]
        if p is None:
            continue
        if a != "-" and b != "-" and in_feat(p):
            denom_core += 1
            if a == b:
                matches_core += 1

    denom_len = 0
    for s, e, *_ in feature_map:
        left = max(ref_abs_start, s)
        right = min(ref_abs_end, e)
        if right > left:
            denom_len += right - left

    pid_core = (100.0 * matches_core / float(denom_core)) if denom_core else 0.0
    pid_len = (100.0 * matches_core / float(denom_len)) if denom_len else 0.0
    return float(pid_core), float(pid_len)


def align_read_to_concat(
    concat_ref: str,
    feature_map: Sequence[tuple[int, int, str, str]],
    read_seq: str,
    *,
    k_index: Optional[KmerIndex] = None,
    k: int = 11,
    margin: int = 200,
    step: int = 2,
    try_both_strands: bool = True,
    match: float = 1.0,
    mismatch: float = -1.0,
    gap_open: float = -1.5,
    gap_extend: float = -0.5,
    pid_mode: Literal["read", "feature", "feature_len"] = "read",
    read_align_impl: Optional[Callable[..., Any]] = None,
    debug: bool = False,
) -> Optional[AlignmentResult]:
    """
    Align a single read to a *feature-concatenated* reference and return a structured result.

    Workflow
    --------
    1) **Seeding**: build/reuse a k-mer index over the concat-band and select a narrow window
       via :func:`seed_window`.
    2) **Read-anchored alignment** inside that window (reference ends free, read ends penalized).
       A low-level per-segment aligner is selected (in this order if not injected):
       ``read_anchored_align`` → ``read_anchored_align_pairwise2`` → ``read_anchored_align_edlib``.
    3) **Strand selection**: optionally evaluate reverse complement and keep the better strand
       by PID (ties by score).
    4) **Feature coverage**: summarize overlap of the chosen reference span with annotated
       feature blocks.
    5) **PID mode**: return PID as requested:
       - ``"read"``        → matches / read bases in the aligned core
       - ``"feature"``     → matches / aligned columns **inside features** only
       - ``"feature_len"`` → matches / **reference feature length** inside the final span

    Parameters
    ----------
    concat_ref : str
        Feature-concatenated reference sequence (uppercase recommended).
    feature_map : sequence of (start, end, ftype, fname)
        Feature blocks on the concat-band, 0-based half-open.
    read_seq : str
        Read sequence (uppercased internally).
    k_index : dict or None, optional
        Optional prebuilt k-mer index for ``concat_ref``. If ``None``, it is built on demand.
    k : int, default 11
        K-mer length for seeding (must match the index).
    margin : int, default 200
        Extra bases added left/right to the seeded window (robustness vs speed).
    step : int, default 2
        Seeding stride on the read (every ``step``-th k-mer).
    try_both_strands : bool, default True
        If True, evaluate forward and reverse complement and keep the better.
    match, mismatch, gap_open, gap_extend : float
        Scoring knobs passed to PairwiseAligner-based paths (ignored by edlib).
    pid_mode : {"read","feature","feature_len"}, default "read"
        Mode of the **returned** percent identity.
    read_align_impl : callable or None, optional
        Per-segment aligner with signature::

            (ref_segment: str, read_seq: str, **kwargs) ->
                (pid, core_len, start_ref, end_ref, start_read, end_read, score, snps)

        If ``None``, a best-effort default is chosen (see *Workflow*).
    debug : bool, default False
        If True, prints which low-level aligner was chosen for easier diagnosis.

    Returns
    -------
    AlignmentResult or None
        Structured result (PID, score, chosen strand, absolute reference span, read span,
        SNP list, feature coverage). Returns ``None`` if no alignment was produced.

    Notes
    -----
    - Feature-based PID requires absolute coordinates; this function passes
      ``pid_mode``, ``feature_map`` and ``ref_abs_start`` down to the aligner.
    - Kwarg filtering ensures compatibility with older helper implementations
      (e.g., mapping ``snp_factory`` → ``snp_cls`` if needed).
    """
    from typing import cast

    # -- local helpers ---------------------------------------------------------
    # (a) robust revcomp if global not present
    try:
        _revcomp = revcomp
    except NameError:
        _RC = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}  # noqa: N806

        def _revcomp(seq: str) -> str:
            return "".join(_RC.get(b, "N") for b in reversed(seq.upper()))

    # (b) safe call into arbitrary low-level aligner
    def _call_aligner_safely(read_align_impl, ref_seg, r_used, align_kwargs):
        """
        Calls the concrete read aligner and normalizes the result to an 8-tuple:
        (pid, core_len, sA, eA, sB, eB, score, snps_raw).

        Returns:
            tuple[float, int, int, int, int, int, float, list] or None
        """
        try:
            res = read_align_impl(ref_seg, r_used, **align_kwargs)
        except Exception:
            # Swallow aligner-specific exceptions and signal a soft failure
            return None

        if res is None:
            return None

        # --- Normalization ---
        # Support dict-like outputs as well as tuple/list/namedtuple.
        if isinstance(res, dict):
            pid_ret = float(res.get("pid"))
            core_len = int(res.get("core_len"))
            sA = int(res.get("sA"))  # noqa: N806
            eA = int(res.get("eA"))  # noqa: N806
            sB = int(res.get("sB"))  # noqa: N806
            eB = int(res.get("eB"))  # noqa: N806
            score = float(res.get("score"))
            snps_raw = res.get("snps") or res.get("snps_raw") or []
            if not isinstance(snps_raw, (list, tuple)):
                snps_raw = [snps_raw] if snps_raw is not None else []
            return (pid_ret, core_len, sA, eA, sB, eB, score, list(snps_raw))

        # Sequence-like outputs (tuple/list/namedtuple)
        try:
            n = len(res)
        except TypeError:
            raise ValueError(  # noqa: B904
                f"Unexpected aligner return type: {type(res)!r}"
            )  # noqa: B904

        if n == 8:
            pid_ret, core_len, sA, eA, sB, eB, score, snps_raw = res  # noqa: N806
        elif n == 7:
            # Aligner did not provide SNPs → default to empty list
            pid_ret, core_len, sA, eA, sB, eB, score = res  # noqa: N806
            snps_raw = []
        else:
            raise ValueError(f"Unexpected aligner return length: {n}")

        # Type safety
        pid_ret = float(pid_ret)
        core_len = int(core_len)
        sA, eA = int(sA), int(eA)  # noqa: N806
        sB, eB = int(sB), int(eB)  # noqa: N806
        score = float(score)
        if snps_raw is None:
            snps_raw = []
        elif not isinstance(snps_raw, (list, tuple)):
            snps_raw = [snps_raw]

        return (pid_ret, core_len, sA, eA, sB, eB, score, list(snps_raw))

    # -- choose aligner --------------------------------------------------------
    if read_align_impl is None:
        read_align_impl = (
            cast(Optional[Callable[..., Any]], globals().get("read_anchored_align"))
            or cast(
                Optional[Callable[..., Any]],
                globals().get("read_anchored_align_pairwise2"),
            )
            or cast(
                Optional[Callable[..., Any]], globals().get("read_anchored_align_edlib")
            )
        )
    if read_align_impl is None:
        if debug:
            print("[align] No low-level aligner available (None).")
        return None
    if debug:
        nm = getattr(read_align_impl, "__name__", str(read_align_impl))
        print(f"[align] Using low-level aligner: {nm}")

    # -- build/reuse k-index ---------------------------------------------------
    kidx: KmerIndex = k_index or build_kmer_index(concat_ref, k=k)

    # NOTE: extended tuple to include feature-only PID for ranking (last element)
    best: Optional[tuple[str, float, int, int, int, int, float, list[Any], float]] = (
        None
    )
    read_f = read_seq.upper()
    strands: tuple[str, ...] = ("F", "R") if try_both_strands else ("F",)

    for strand in strands:
        r_used = read_f if strand == "F" else _revcomp(read_f)

        win = seed_window(concat_ref, r_used, kidx, k=k, margin=margin, step=step)
        if win is None:
            continue
        w_s, w_e, _hits = win
        ref_seg = concat_ref[w_s:w_e]

        # build kwargs for the aligner
        align_kwargs: dict[str, Any] = {
            "pid_mode": pid_mode,
            "feature_map": feature_map,
            "ref_abs_start": w_s,  # absolute offset of this segment
            "match": match,
            "mismatch": mismatch,
            "gap_open": gap_open,
            "gap_extend": gap_extend,
        }
        # SNP factory → dataclass (if available); else dicts
        try:
            align_kwargs["snp_factory"] = lambda **rec: SNP(**rec)
        except NameError:
            align_kwargs["snp_factory"] = lambda **rec: dict(**rec)

        res = _call_aligner_safely(read_align_impl, ref_seg, r_used, align_kwargs)
        if res is None:
            continue

        pid_ret, core_len, sA, eA, sB, eB, score, snps_raw = res  # noqa: N806
        sA_abs, eA_abs = w_s + sA, w_s + eA  # noqa: N806

        # --- Feature-only PID for ranking robustness -------------------------
        # We approximate feature-only PID by scaling the reported PID with the
        # fraction of core reference columns that overlap annotated features.
        # If the low-level aligner can provide per-column masks/CIGAR-expansion,
        # consider computing exact in-feature matches instead of this heuristic.
        core_bp = max(eA_abs - sA_abs, 1)
        feat_bp = 0
        if feature_map:
            span_start, span_end = sA_abs, eA_abs
            for f_s, f_e, ftype, fname in feature_map:  # noqa: B007
                ovl_s = max(span_start, f_s)
                ovl_e = min(span_end, f_e)
                if ovl_e > ovl_s:
                    feat_bp += ovl_e - ovl_s
        coverage_frac = (feat_bp / core_bp) if core_bp else 0.0
        matched_bp_est = (pid_ret / 100.0) * core_bp
        feat_matches_est = matched_bp_est * coverage_frac
        feat_pid_for_rank = (
            (100.0 * feat_matches_est / max(feat_bp, 1)) if feat_bp > 0 else 0.0
        )

        # Build candidate; keep reported PID as-is (based on pid_mode)
        candidate = (
            strand,  # 0
            float(pid_ret),  # 1 reported/displayed PID
            int(sA_abs),
            int(eA_abs),  # 2-3 absolute ref span
            int(sB),
            int(eB),  # 4-5 read span
            float(score),  # 6 score
            list(snps_raw),  # 7 SNPs
            float(feat_pid_for_rank),  # 8 ranking PID (feature-only)
        )

        # Prefer higher feature-only PID; tie-break by reported PID, then by score
        if (
            best is None
            or candidate[8] > best[8]
            or (candidate[8] == best[8] and candidate[1] > best[1])
            or (
                candidate[8] == best[8]
                and candidate[1] == best[1]
                and candidate[6] > best[6]
            )
        ):
            best = candidate

    if best is None:
        return None

    strand, pid_final, sA_abs, eA_abs, sB, eB, score, snps_raw, _feat_rank = best  # noqa: N806  # fmt: skip
    cov = summarize_feature_coverage(feature_map, (sA_abs, eA_abs))

    if debug:
        print(
            "[align] feature_hits type:",
            type(cov),
            "sample:",
            (cov[:2] if isinstance(cov, (list, tuple)) else cov),
        )

    try:
        snps_cast = cast(list[SNP], snps_raw)
    except Exception:
        snps_cast = []

    return AlignmentResult(
        pid=pid_final,
        core_len=(eB - sB),
        score=score,
        strand=strand,
        start_ref=sA_abs,
        end_ref=eA_abs,
        start_read=sB,
        end_read=eB,
        snps=snps_cast,
        feature_hits=cov,
    )


def read_anchored_align_pairwise2(ref_segment: str, read_seq: str, **kwargs):
    """
    Backward-compatible shim: maps 'pairwise2' to the new implementation and
    coerces the return value to the legacy 7-tuple:
    (ref_aln, read_aln, start_ref, end_ref, strand, score, anchors)
    """
    res = read_anchored_align(ref_segment, read_seq, **kwargs)

    # If the new class is used:
    try:
        ReadAlignResult = kwargs.get("ReadAlignResult")  # noqa: N806
    except Exception:
        ReadAlignResult = None  # noqa: N806

    if ReadAlignResult and isinstance(res, ReadAlignResult):
        # Class → 7-tuple
        return (
            res.ref_aln,
            res.read_aln,
            res.start_ref,
            res.end_ref,
            res.strand,
            res.score,
            getattr(res, "anchors", None),
        )

    # Tuple → truncate to 7 fields (e.g., if 8 including debug is returned)
    if isinstance(res, tuple):
        if len(res) >= 7:
            return res[:7]
        raise ValueError(f"aligner returned {len(res)} fields; expected ≥7")

    raise TypeError(f"unexpected return type from read_anchored_align: {type(res)!r}")


def read_anchored_align_edlib(
    ref_segment: str,
    read_seq: str,
    *,
    snp_factory: Optional[Callable[..., Any]] = None,
    max_error_rate: float = 0.40,
    **kwargs: Any,
) -> Optional[tuple[float, int, int, int, int, int, float, list[Any]]]:
    """
    Read-anchored Semiglobal-Alignment with edlib (HW-Mode).

    Semantics
    --------
    - Anchoring: READ must be fully covered (query=read); REF may be a substring.
      → edlib mode="HW". Equivalent to pairwise2 with penalize_end_gaps=(False, True).
    - PID: (Matches / len(read)) * 100 (read-based denominator, stable across read lengths).
    - core_len: # Cleavage sites where both strands consume a base (including mismatches).
    - score: -EditDistance (bigger = better).

    Parameter
    ---------
    ref_segment : str
        Reference-Slice (e.g. Seed-Window + Edge).
    read_seq : str
        Read-Sequence. Internally uppercased.
    snp_factory : callable | None
        Factory/functor for SNP objects, signature:
            snp_factory(ref_pos=..., ref_base=..., read_pos=..., read_base=...)
        If None, dictionaries are created.
    max_error_rate : float
        Early rejection, if editDistance > rate * len(read).
    **kwargs : Any
        Unexpected or additional arguments are ignored (robust API).

    Return value
    ------------
    (pid, core_len, start_ref, end_ref, start_read, end_read, score, snps) or None
    """

    try:
        import edlib
    except Exception as exc:  # pragma: no cover
        raise ImportError(
            "edlib is required for read_anchored_align_edlib requires (pip install edlib)."
        ) from exc

    read = read_seq.upper()
    ref = ref_segment.upper()
    read_len = len(read)
    if read_len == 0 or len(ref) == 0:
        return None

    result = edlib.align(read, ref, mode="HW", task="path")
    edit_distance = int(result.get("editDistance", -1))
    if edit_distance < 0:
        return None

    if (max_error_rate is not None) and (
        edit_distance > int(max_error_rate * read_len)
    ):
        return None

    locations = result.get("locations") or []
    if not locations:
        return None

    start_ref = int(locations[0][0])
    end_ref = int(locations[0][1]) + 1

    cigar = result.get("cigar", "")
    if not cigar:
        return None

    def _make_snp(ref_pos: int, ref_base: str, read_pos: int, read_base: str) -> Any:
        if snp_factory is None:
            return {
                "ref_pos": ref_pos,
                "ref_base": ref_base,
                "read_pos": read_pos,
                "read_base": read_base,
            }
        return snp_factory(
            ref_pos=ref_pos, ref_base=ref_base, read_pos=read_pos, read_base=read_base
        )

    def _ops(cig: str):
        n = 0
        for ch in cig:
            if ch.isdigit():
                n = n * 10 + (ord(ch) - 48)
            else:
                yield (n if n else 1, ch)
                n = 0

    i_ref = start_ref  # Index in ref
    i_read = 0  # Index in read
    matches = 0
    core_len = 0
    snps: list[Any] = []

    for length, op in _ops(cigar):
        if op == "M":
            for j in range(length):
                rb = ref[i_ref + j]
                qb = read[i_read + j]
                core_len += 1
                if rb == qb:
                    matches += 1
                else:
                    snps.append(_make_snp(i_ref + j, rb, i_read + j, qb))
            i_ref += length
            i_read += length

        elif op == "I":
            for j in range(length):
                qb = read[i_read + j]
                snps.append(_make_snp(i_ref, "-", i_read + j, qb))
            i_read += length

        elif op == "D":
            for j in range(length):
                rb = ref[i_ref + j]
                snps.append(_make_snp(i_ref + j, rb, i_read, "-"))
            i_ref += length

        else:
            continue

    pid = 100.0 * (matches / read_len)
    score = float(-edit_distance)
    start_read = 0
    end_read = read_len

    return (
        float(pid),
        int(core_len),
        int(start_ref),
        int(end_ref),
        int(start_read),
        int(end_read),
        float(score),
        snps,
    )
