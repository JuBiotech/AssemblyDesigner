"""PCR / Gibson assembly simulation, plus the 3G batch pipeline (Golden Gate → PCR → Gibson)."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping, MutableMapping
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, NamedTuple, Optional, Union

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from .assembly import generate_recipe_assemblies, organize_assembly_reports
from .template import _as_bool, _get_int


class PCRResult(NamedTuple):
    """Container for PCR simulation output."""

    product: SeqRecord
    fwd_start: int
    rev_end: int
    # Defaults keep backwards compatibility with 3-field returns.
    fwd_anneal_len: int = 0
    rev_anneal_len: int = 0


def strip_terminal_ambiguous_bases(
    record: SeqRecord, *, ambiguous: str = "N", min_run: int = 1
) -> SeqRecord:
    """
    Trim leading/trailing runs of an ambiguous base (e.g., 'N').

    Parameters
    ----------
    record : Bio.SeqRecord.SeqRecord
        Input sequence record (linear or circular).
    ambiguous : str, optional
        Base to strip from the ends, by default "N".
    min_run : int, optional
        Minimum run length to consider for stripping, by default 1.

    Returns
    -------
    Bio.SeqRecord.SeqRecord
        New record with terminal runs removed (internal ambiguous bases are kept).
    """
    s = str(record.seq)
    a = ambiguous.upper()

    # Left scan
    i = 0
    while i < len(s) and s[i].upper() == a:
        i += 1

    # Right scan
    j = len(s)
    while j > i and s[j - 1].upper() == a:
        j -= 1

    # Thresholds
    if i < min_run:
        i = 0
    if (len(s) - j) < min_run:
        j = len(s)

    out = record[i:j]
    out.id = record.id
    out.name = record.name
    if not getattr(out, "description", ""):
        out.description = record.description
    out.annotations.setdefault("molecule_type", "DNA")
    return out


def _shift_location(
    loc: Union[FeatureLocation, CompoundLocation], offset: int
) -> Union[FeatureLocation, CompoundLocation]:
    """
    Shift a feature location by a fixed offset.

    Parameters
    ----------
    loc : FeatureLocation or CompoundLocation
        Original location to shift.
    offset : int
        Offset in bases to add to start/end.

    Returns
    -------
    FeatureLocation or CompoundLocation
        Shifted location with the same structure/strand.
    """
    if isinstance(loc, CompoundLocation):
        return CompoundLocation(
            [_shift_location(p, offset) for p in loc.parts],
            operator=loc.operator,
        )
    return FeatureLocation(
        int(loc.start) + offset, int(loc.end) + offset, strand=loc.strand
    )


def _copy_feature_with_shift(feat: SeqFeature, offset: int) -> SeqFeature:
    """Copy a feature and shift its location by `offset` bases."""
    return SeqFeature(
        location=_shift_location(feat.location, offset),
        type=feat.type,
        qualifiers={
            k: (list(v) if isinstance(v, list) else v)
            for k, v in feat.qualifiers.items()
        },
    )


def _concat_records_with_features(
    parts: list[SeqRecord], new_id: str = "CONCAT"
) -> SeqRecord:
    """
    Concatenate SeqRecords while preserving and offsetting features.

    Parameters
    ----------
    parts : list of Bio.SeqRecord.SeqRecord
        Parts to concatenate in order.
    new_id : str, optional
        ID/name for the resulting record, by default "CONCAT".

    Returns
    -------
    Bio.SeqRecord.SeqRecord
        Concatenated record with merged, shifted features.
    """
    seq_chunks: list[str] = []
    feats: list[SeqFeature] = []
    offset = 0

    for rec in parts:
        s = str(rec.seq)
        seq_chunks.append(s)
        for f in rec.features:
            feats.append(_copy_feature_with_shift(f, offset))
        offset += len(s)

    out = SeqRecord(
        Seq("".join(seq_chunks)), id=new_id, name=new_id, description="concatenated"
    )
    out.features = feats
    out.annotations.setdefault("molecule_type", "DNA")
    return out


def _tag_overhang(seq: str, label: str) -> SeqRecord:
    """
    Create a mini-SeqRecord for an overhang and tag it as a misc_feature.

    Parameters
    ----------
    seq : str
        Overhang sequence.
    label : str
        Label shown in the feature's qualifiers.

    Returns
    -------
    Bio.SeqRecord.SeqRecord
        Record containing `seq` and a single misc_feature covering it.
    """
    r = SeqRecord(Seq(seq), id=label, name=label, description=label)
    if seq:
        r.features = [
            SeqFeature(
                FeatureLocation(0, len(seq)),
                type="misc_feature",
                qualifiers={"label": [label], "note": ["PCR overhang"]},
            )
        ]
    return r


def simulate_pcr(
    template: SeqRecord,
    *,
    fwd_primer: str,  # 5'→3'
    rev_primer: str,  # 5'→3' (search RC)
    include_primers: bool = True,
    trim_terminal_N: bool = True,  # noqa: N803
    enforce_unique: bool = True,
    circular: bool = False,
) -> PCRResult:
    """
    Simulate a PCR on a linear view of `template` using exact primer matches.

    Parameters
    ----------
    template : Bio.SeqRecord.SeqRecord
        Template sequence; treated as linear for indexing.
    fwd_primer : str
        Forward primer (5'→3'), matched exactly.
    rev_primer : str
        Reverse primer (5'→3'); its reverse complement is matched exactly.
    include_primers : bool, optional
        If True, include primer sequences in the amplicon, by default True.
    trim_terminal_N : bool, optional
        If True, strip terminal 'N' runs from the product, by default True.
    enforce_unique : bool, optional
        If True, require exactly one FWD and one REV site, by default True.
    circular : bool, optional
        If True, allow wrap-around amplicons on circular templates, by default False.

    Returns
    -------
    PCRResult
        Tuple-like result with product record and site indices.

    Raises
    ------
    ValueError
        If sites are not found (or not unique when `enforce_unique=True`), or if
        the reverse site is upstream of the forward site on linear templates.

    Notes
    -----
    - Slicing (`record[start:end]`) preserves features; wrap-around uses a
      feature-preserving concat.
    """
    seq_str = str(template.seq).upper()
    fwd = fwd_primer.upper()
    rev_rc = str(Seq(rev_primer).reverse_complement()).upper()

    if enforce_unique:
        f_hits = seq_str.count(fwd)
        r_hits = seq_str.count(rev_rc)
        if f_hits != 1 or r_hits != 1:
            raise ValueError(
                f"Primer sites must be unique (found fwd={f_hits}, rev={r_hits})."
            )

    fpos = seq_str.find(fwd)
    rpos = seq_str.find(rev_rc)
    if fpos < 0 or rpos < 0:
        raise ValueError("Primer binding site not found (exact match required).")

    if not circular:
        if rpos <= fpos:
            raise ValueError(
                "Reverse primer site must be downstream of forward primer."
            )
        start = fpos if include_primers else (fpos + len(fwd))
        end = (rpos + len(rev_rc)) if include_primers else rpos
        product = template[start:end]  # slicing preserves features
    else:
        # Circular: either standard span (rpos > fpos) or wrap-around (rpos < fpos)
        if rpos > fpos:
            start = fpos if include_primers else (fpos + len(fwd))
            end = (rpos + len(rev_rc)) if include_primers else rpos
            product = template[start:end]
        else:
            part1_start = fpos if include_primers else (fpos + len(fwd))
            part1 = template[part1_start:]
            part2_end = (rpos + len(rev_rc)) if include_primers else rpos
            part2 = template[:part2_end]
            product = _concat_records_with_features(
                [part1, part2], new_id=f"{template.id}_PCR"
            )

    product.id = product.name = f"{template.id}_PCR"
    if not getattr(product, "description", ""):
        product.description = "PCR product"

    if trim_terminal_N:
        product = strip_terminal_ambiguous_bases(product, ambiguous="N", min_run=1)

    rev_end = (rpos + len(rev_rc)) % len(template) if circular else (rpos + len(rev_rc))
    return PCRResult(product=product, fwd_start=fpos, rev_end=rev_end)


def simulate_pcr_overhangs(
    template: SeqRecord,
    *,
    fwd_primer: str,  # 5'→3' including overhang
    rev_primer: str,  # 5'→3' including overhang (RC searched)
    min_anneal: int = 14,  # minimal annealing length to accept
    circular: bool = False,
    include_overhangs: bool = True,
    product_id: Optional[str] = None,
    trim_terminal_N: bool = False,  # noqa: N803
) -> PCRResult:
    """
    Simulate a PCR where primers carry 5' overhangs; only the annealing parts bind.

    The function automatically finds the **longest unique annealing** suffix (forward
    primer) and prefix (reverse primer in RC) of length ≥ `min_anneal`. It returns the
    product including the annealing regions and, optionally, the 5' overhangs as
    sequence flanks annotated with `misc_feature`. Features from the template are
    preserved via slicing/feature-aware concatenation.

    Parameters
    ----------
    template : Bio.SeqRecord.SeqRecord
        Template sequence; treated as linear for indexing (wrap-around handled when
        `circular=True`).
    fwd_primer : str
        Forward primer (5'→3') possibly with a 5' overhang.
    rev_primer : str
        Reverse primer (5'→3') possibly with a 5' overhang; RC is searched.
    min_anneal : int, optional
        Minimal length for the annealing part to be considered, by default 14.
    circular : bool, optional
        If True, allow wrap-around amplicons on circular templates, by default False.
    include_overhangs : bool, optional
        If True, include 5' overhang sequences flanking the product and tag them
        as `misc_feature`, by default True.
    product_id : str | None, optional
        ID for the resulting product; if None, defaults to `<template.id>_PCR`.
    trim_terminal_N : bool, optional
        If True, strip terminal 'N' runs from the final product, by default False.

    Returns
    -------
    PCRResult
        Tuple-like result with product record and site indices; also encodes the
        annealing lengths in `fwd_anneal_len` and `rev_anneal_len`.

    Raises
    ------
    ValueError
        If annealing sites cannot be found uniquely for the given `min_anneal`, or if
        the reverse site is upstream of the forward site on linear templates.

    Notes
    -----
    - Overhangs are defined as the non-annealing 5' segments of each primer.
      The right overhang is appended in template orientation (i.e., suffix of RC).
    - If `include_overhangs=False`, only the annealed core product is returned.
    """

    def _longest_unique_suffix_hit(
        seq: str, primer: str, min_len: int
    ) -> tuple[int, int]:
        p = primer.upper()
        for L in range(len(p), min_len - 1, -1):  # noqa: N806
            motif = p[-L:]
            if seq.count(motif) == 1:
                return seq.find(motif), L
        return -1, 0

    def _longest_unique_prefix_hit(
        seq: str, rc_primer: str, min_len: int
    ) -> tuple[int, int]:
        p = rc_primer.upper()
        for L in range(len(p), min_len - 1, -1):  # noqa: N806
            motif = p[:L]
            if seq.count(motif) == 1:
                return seq.find(motif), L
        return -1, 0

    seq_str = str(template.seq).upper()
    fwd = fwd_primer.upper()
    rev_rc = str(Seq(rev_primer).reverse_complement()).upper()

    fpos, fL = _longest_unique_suffix_hit(seq_str, fwd, min_anneal)  # noqa: N806
    rpos, rL = _longest_unique_prefix_hit(seq_str, rev_rc, min_anneal)  # noqa: N806
    if fpos < 0 or rpos < 0:
        raise ValueError(
            f"Anneal not found uniquely (fwd_len={fL}, rev_len={rL}). "
            "Check primers or increase `min_anneal`."
        )

    # Core product (includes the annealing segments)
    if not circular:
        if rpos <= fpos:
            raise ValueError(
                "Reverse site must be downstream of forward site (linear)."
            )
        core = template[fpos : rpos + rL]
    else:
        if rpos >= fpos:
            core = template[fpos : rpos + rL]
        else:
            core = _concat_records_with_features(
                [template[fpos:], template[: rpos + rL]],
                new_id="PCR_core",
            )

    # Optional 5' overhangs
    left_ov = fwd[:-fL]  # 5' segment before the annealing suffix
    right_ov = rev_rc[
        rL:
    ]  # suffix of RC (5' overhang of reverse in template orientation)

    parts: list[SeqRecord] = []
    if include_overhangs:
        parts.append(_tag_overhang(left_ov, "left_overhang"))
    parts.append(core)
    if include_overhangs:
        parts.append(_tag_overhang(right_ov, "right_overhang"))

    product = _concat_records_with_features(
        parts, new_id=product_id or f"{template.id}_PCR"
    )
    product.description = "PCR product with 5' overhangs (features preserved)"

    if trim_terminal_N:
        product = strip_terminal_ambiguous_bases(product, ambiguous="N", min_run=1)

    rev_end = (rpos + rL) % len(template) if circular else (rpos + rL)
    return PCRResult(
        product=product,
        fwd_start=fpos,
        rev_end=rev_end,
        fwd_anneal_len=fL,
        rev_anneal_len=rL,
    )


def longest_overlap(a, b, min_overlap: int = 1):
    """
    Longest 3' overlap length between two sequences (strings or SeqRecords).

    Returns the largest L >= min_overlap such that:
      - str(a).upper().endswith(str(b)[:L].upper())  OR
      - str(b).upper().endswith(str(a)[:L].upper())
    """
    sa = str(getattr(a, "seq", a)).upper()
    sb = str(getattr(b, "seq", b)).upper()
    if not sa or not sb:
        return 0
    maxlen = min(len(sa), len(sb))
    for L in range(maxlen, max(min_overlap - 1, 0), -1):  # noqa: N806
        if sa.endswith(sb[:L]):
            return L
    return 0


# Back-compat alias used elsewhere in the file
_longest_overlap = longest_overlap


def merge_with_gibson_features(
    frags: list[SeqRecord],
    *,
    min_overlap: int = 20,
    circularize: bool = True,
) -> SeqRecord:
    """
    Assemble fragments by Gibson-style overlaps while preserving features.

    Fragments are concatenated in order using the longest suffix/prefix overlap
    (length >= `min_overlap`) between the current assembled sequence and the next
    fragment. Only the non-overlapping tail of the next fragment is appended,
    and the concatenation is performed via `_concat_records_with_features(...)`
    so that features (and their coordinates) are preserved.

    If `circularize=True`, the terminal overlap between the final linear assembly
    and the first fragment is removed from the end (to avoid duplication),
    yielding a linear representation of the closed circle.

    Parameters
    ----------
    frags
        Input fragments in assembly order (already oriented correctly).
    min_overlap
        Minimum overlap length at each junction.
    circularize
        If True, enforce and trim the terminal (last→first) overlap.

    Returns
    -------
    SeqRecord
        Assembled (linear) record with features preserved.

    Raises
    ------
    ValueError
        If a junction has insufficient overlap.
    ImportError
        If `_concat_records_with_features` is not available.
    """
    if not frags:
        raise ValueError("No fragments provided.")
    if min_overlap <= 0:
        raise ValueError("min_overlap must be a positive integer.")

    # Ensure the feature-preserving concat helper exists
    if not callable(globals().get("_concat_records_with_features")):
        raise ImportError(
            "Feature-preserving concatenation requires `_concat_records_with_features(parts, new_id=...)`."
        )

    assembled = frags[0]

    # Append each next fragment's non-overlapping tail
    for i in range(len(frags) - 1):
        a_rec = assembled
        b_rec = frags[i + 1]
        ov = _longest_overlap(str(a_rec.seq), str(b_rec.seq), min_overlap)
        if ov < min_overlap:
            raise ValueError(
                f"No sufficient overlap between fragment {i} and {i + 1} "
                f"(found {ov}, need ≥ {min_overlap})."
            )
        b_tail = b_rec[ov:]  # slicing preserves B's features
        assembled = _concat_records_with_features(
            [a_rec, b_tail], new_id="GIBSON_linear"
        )

    # Enforce and trim the terminal overlap when circularizing
    if circularize:
        terminal_ov = _longest_overlap(
            str(assembled.seq), str(frags[0].seq), min_overlap
        )
        if terminal_ov < min_overlap:
            raise ValueError(
                f"No sufficient terminal overlap to circularize "
                f"(found {terminal_ov}, need ≥ {min_overlap})."
            )
        if terminal_ov:
            assembled = assembled[:-terminal_ov]
        assembled.id = "GIBSON_3G"
        assembled.name = "GIBSON_3G"
        assembled.description = "Assembled (Gibson), features preserved"

    return assembled


@dataclass(slots=True)
class BuildStatus:
    """Per-construct status after Golden Gate → PCR → Gibson."""

    cid: str
    gg_ok: bool = False
    gg_msg: str = ""
    pcr_ok: bool = False
    pcr_msg: str = ""
    gibson_ok: bool = False
    gibson_msg: str = ""
    final_path: Path | None = None
    final_bp: int | None = None


class BatchResult(NamedTuple):
    """Return container for high-throughput builds."""

    status_df: pd.DataFrame
    products: dict[tuple[str, str], SeqRecord]
    pcr_results: dict[tuple[str, str], object]
    finals: dict[str, SeqRecord]


def run_3g_batch_safe(
    *,
    category_dirs: Mapping[str, Path],
    designs: Iterable[Mapping[str, Any]],
    category_order: Iterable[str],
    reports_dir: Path,
    pcr_plan_df: pd.DataFrame,
    assembly_df: pd.DataFrame,
    tus_df: pd.DataFrame,
    constructs_df: pd.DataFrame,
    template_lookup: Callable[[str, str], SeqRecord],
) -> BatchResult:
    """Run Golden Gate → PCR → Gibson across many constructs.

    Parameters
    ----------
    category_dirs
        Mapping ``{category: folder}`` containing the part libraries.
    designs
        Iterable of explicit TU recipes (usually rows from the ``TUs`` sheet).
    category_order
        Order of categories (e.g., ``["Promoter","RBS","Gene","Terminator","UNS_Context"]``).
    reports_dir
        Output directory for report ZIPs and extracted GenBank files.
    pcr_plan_df
        ``PCR`` sheet.
    assembly_df
        ``Assembly`` sheet (defines fragment order per construct).
    tus_df
        ``TUs`` sheet (may be empty; used by downstream history tools).
    constructs_df
        ``Constructs`` sheet (``MinOverlap`` / ``Circularize`` optional).
    template_lookup
        Callable that resolves the PCR template for ``(cid, role)``.

    Returns
    -------
    BatchResult
        Tuple containing:
        * ``status_df`` : per-CID status table
        * ``products``  : PCR products keyed by ``(cid, role)``
        * ``pcr_results`` : raw PCR result objects (opaque)
        * ``finals``    : Gibson finals keyed by ``cid``

    Notes
    -----
    * Golden Gate report ZIPs are written to ``reports_dir``.
    * TU finals (GenBank) are extracted to ``reports_dir/Assembly``.
    * Gibson finals (GenBank) are written to ``reports_dir/Finals/<CID>_final.gb``.
    """
    reports_dir = Path(reports_dir)
    reports_dir.mkdir(parents=True, exist_ok=True)

    # Initialize per-CID status
    statuses: dict[str, BuildStatus] = {
        str(cid): BuildStatus(cid=str(cid))
        for cid in sorted(assembly_df["ConstructID"].astype(str).unique())
    }

    # ---- 1) Golden Gate (DNAcauldron) --------------------------------
    try:
        gg_reports = generate_recipe_assemblies(
            category_dirs=category_dirs,
            designs=list(designs),
            category_order=list(category_order),
            output_dir=reports_dir,
        )
    except Exception as exc:  # pylint: disable=broad-except
        raise RuntimeError(
            f"Golden Gate stage failed before per-construct split: {exc}"
        ) from exc

    try:
        _ = organize_assembly_reports(
            report_dir=reports_dir,
            reports=gg_reports,
            delete_zip=False,
            final_only=True,
        )
    except Exception as exc:  # pylint: disable=broad-except
        # Continue; some reports may still have finals
        print(f"⚠️ organize_assembly_reports: {exc}")

    assembly_out_dir = reports_dir / "Assembly"
    tu_final_paths = (
        list(assembly_out_dir.glob("*.gb*")) if assembly_out_dir.is_dir() else []
    )
    if tu_final_paths:
        for s in statuses.values():
            s.gg_ok = True
    else:
        for s in statuses.values():
            s.gg_ok = False
            s.gg_msg = f"No TU finals found in {assembly_out_dir}"

    # ---- 2) PCR -------------------------------------------------------
    products: dict[tuple[str, str], SeqRecord] = {}
    pcr_results: dict[tuple[str, str], object] = {}

    pcr_ok_map: dict[str, bool] = {cid: True for cid in statuses}
    pcr_msgs_map: dict[str, list[str]] = {cid: [] for cid in statuses}

    for _, row in pcr_plan_df.iterrows():
        cid = str(row["ConstructID"])
        role = str(row["Role"])
        if cid not in statuses or not statuses[cid].gg_ok:
            continue

        try:
            fwd = str(row["FwdPrimer"]).replace(" ", "")
            rev = str(row["RevPrimer"]).replace(" ", "")
            circ = _as_bool(row.get("Circular", True), True)
            use_ovh = _as_bool(
                row.get("IncludeOverhangs", role.upper() == "BACKBONE"),
                role.upper() == "BACKBONE",
            )
            min_anneal = _get_int(row.get("MinAnneal", 18), 18)

            template = template_lookup(cid, role)
            pid = f"{cid}_{role}_PCR"

            if use_ovh:
                res = simulate_pcr_overhangs(
                    template=template,
                    fwd_primer=fwd,
                    rev_primer=rev,
                    circular=circ,
                    min_anneal=min_anneal,
                    include_overhangs=True,
                    product_id=pid,
                )
            else:
                res = simulate_pcr(
                    template=template,
                    fwd_primer=fwd,
                    rev_primer=rev,
                    circular=circ,
                )

            prod = res.product
            meta: MutableMapping[str, Any] = prod.annotations.setdefault("pcr_meta", {})
            meta.setdefault("fwd_primer", fwd)
            meta.setdefault("rev_primer", rev)

            products[(cid, role)] = prod
            pcr_results[(cid, role)] = res

        except Exception as exc:  # pylint: disable=broad-except
            pcr_ok_map[cid] = False
            pcr_msgs_map[cid].append(f"{role}: {exc!s}")

    for cid, status in statuses.items():
        status.pcr_ok = pcr_ok_map[cid] and status.gg_ok
        status.pcr_msg = "; ".join(pcr_msgs_map[cid])

    # ---- 3) Gibson ----------------------------------------------------
    finals: dict[str, SeqRecord] = {}
    finals_dir = reports_dir / "Finals"
    finals_dir.mkdir(parents=True, exist_ok=True)

    for cid, grp in assembly_df.groupby("ConstructID"):
        scid = str(cid)
        status = statuses[scid]
        if not (status.gg_ok and status.pcr_ok):
            continue

        try:
            roles = grp.sort_values("Order")["FragmentRole"].astype(str).tolist()
            try:
                frags = [products[(scid, r)] for r in roles]
            except KeyError as missing:
                raise KeyError(
                    f"Missing PCR product for role {missing!s} (CID {scid})."
                ) from None

            rowc = constructs_df[constructs_df["ConstructID"].astype(str) == scid]
            min_ov = (
                _get_int(rowc["MinOverlap"].iloc[0], 40)
                if (not rowc.empty and "MinOverlap" in rowc)
                else 40
            )
            circ = (
                _as_bool(rowc["Circularize"].iloc[0], True)
                if (not rowc.empty and "Circularize" in rowc)
                else True
            )

            final = merge_with_gibson_features(
                frags, min_overlap=min_ov, circularize=circ
            )
            finals[scid] = final

            status.gibson_ok = True
            status.final_bp = len(final)
            out_path = finals_dir / f"{scid}_final.gb"
            SeqIO.write(final, out_path, "genbank")
            status.final_path = out_path

        except Exception as exc:  # pylint: disable=broad-except
            status.gibson_ok = False
            status.gibson_msg = str(exc)

    status_df = pd.DataFrame([asdict(s) for s in statuses.values()]).set_index("cid")
    return BatchResult(
        status_df=status_df, products=products, pcr_results=pcr_results, finals=finals
    )
