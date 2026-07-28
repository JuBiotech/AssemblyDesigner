"""Feature utilities (display name, filtering, concat band, constructs).

This module provides robust helpers to work with GenBank features:
- `feature_display_name`: pick a readable name from GenBank qualifiers.
- `classify_feature_slot`: map a feature to a parts slot (promoter/rbs/gene/term) + token.
- `filter_features`: filter by feature type and include/exclude tokens.
- `build_feature_concat_sequence`: concatenate selected feature sequences and return a back-map.
- Construct derivation:
  * `derive_construct_from_record`: make "Promoter_RBS_Gene_Term" from a SeqRecord.
  * `derive_construct_from_labels`: same from a list of labels.
  * `derive_construct_from_filename`: fallback: infer parts from a filename.
  * `parts_key_from_construct`: reduce a construct to selected parts (e.g., promoter+RBS).
- `inventory_features`: list unique tokens per slot across records.
"""

from __future__ import annotations

import os
from collections.abc import Iterable, Sequence
from typing import Optional

from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from .lexicon import Lexicon, default_lexicon, match_token

# Include common GenBank feature types you want to consider for MoClo parts.
DEFAULT_ALLOWED_TYPES: tuple[str, ...] = (
    "promoter",
    "RBS",
    "CDS",
    "terminator",
    "misc_feature",
    "5'UTR",
    "five_prime_UTR",
)

# Loose mapping of obvious type→slot (used as a strong hint).
# Anything unknown falls back to name-based token matching.
_TYPE_TO_SLOT: dict[str, str] = {
    "promoter": "promoter",
    "RBS": "rbs",
    "CDS": "gene",
    "terminator": "term",
    "5'UTR": "rbs",
    "five_prime_UTR": "rbs",
    # "misc_feature": no direct mapping → resolve by tokens/heuristics
}


def feature_display_name(feat: SeqFeature) -> str:
    """
    Choose a human-readable display name for a GenBank feature.

    The function checks these qualifiers in order: ``label``, ``gene``,
    ``product``, ``note``, ``locus_tag``. If none exist, the feature
    ``type`` is returned.

    Parameters
    ----------
    feat : SeqFeature
        A Biopython feature from a GenBank record.

    Returns
    -------
    str
        Best-effort display name for the feature.
    """
    qualifiers = getattr(feat, "qualifiers", {}) or {}
    for key in ("label", "gene", "product", "note", "locus_tag"):
        values = qualifiers.get(key)
        if values:
            return str(values[0])
    ftype = str(getattr(feat, "type", "feature"))
    return ftype


def _is_5utr_label(text: str) -> bool:
    """Return True if a label string looks like a 5'UTR."""
    low = text.lower()
    return "utr" in low or "5'utr" in low or "five_prime" in low


def classify_feature_slot(
    feat: SeqFeature,
    lex: Optional[Lexicon] = None,
    *,
    treat_5utr_as_rbs: bool = True,
) -> tuple[Optional[str], Optional[str]]:
    """
    Classify a feature into a parts slot and extract a token.

    Strategy:
    1) If the feature type is a known MoClo slot (`_TYPE_TO_SLOT`), use that slot hint.
    2) Use the feature display name and the lexicon's regex patterns to find a token.
    3) For ``misc_feature`` or unknown types, purely rely on token matching.
    4) If nothing matches, return ``(None, None)``.

    Parameters
    ----------
    feat : SeqFeature
        Feature to classify.
    lex : Lexicon or None
        Parts lexicon; defaults to `default_lexicon()`.
    treat_5utr_as_rbs : bool, optional
        If True, an explicit 5'UTR type or label is mapped to the RBS slot.

    Returns
    -------
    (slot, token) : tuple of (str or None, str or None)
        Slot is one of {"promoter", "rbs", "gene", "term"} or None.
        Token is the canonical token for that slot (e.g., "J23100", "BCD12", "ecPanD").
    """
    use_lex = lex or default_lexicon()

    ftype = str(getattr(feat, "type", "") or "")
    slot_hint = _TYPE_TO_SLOT.get(ftype)
    name = feature_display_name(feat).strip()

    # 5'UTR labeled features should count as RBS (common in GB annotations).
    if treat_5utr_as_rbs and (
        ftype in ("5'UTR", "five_prime_UTR") or _is_5utr_label(name)
    ):
        slot_hint = "rbs"

    # Try to extract tokens for all slots from the name.
    tok_p = match_token(name, "promoter", use_lex)
    tok_r = match_token(name, "rbs", use_lex)
    tok_g = match_token(name, "gene", use_lex)
    tok_t = match_token(name, "term", use_lex)

    # Prefer the token consistent with the slot hint, if present.
    if slot_hint == "promoter" and tok_p:
        return "promoter", tok_p
    if slot_hint == "rbs" and tok_r:
        return "rbs", tok_r
    if slot_hint == "gene" and tok_g:
        return "gene", tok_g
    if slot_hint == "term" and tok_t:
        return "term", tok_t

    # Otherwise, decide by strongest name-based evidence.
    if tok_p:
        return "promoter", tok_p
    if tok_r:
        return "rbs", tok_r
    if tok_g:
        return "gene", tok_g
    if tok_t:
        return "term", tok_t

    return None, None


def _token_hit(name_lower: str, tokens: Sequence[str]) -> bool:
    """Return True if any token (case-insensitive) is contained in `name_lower`."""
    return any(tok.lower() in name_lower for tok in tokens)


def filter_features(
    features: Iterable[SeqFeature],
    *,
    include_tokens: Optional[Sequence[str]] = None,
    exclude_tokens: Optional[Sequence[str]] = None,
    allowed_types: Sequence[str] = DEFAULT_ALLOWED_TYPES,
) -> list[SeqFeature]:
    """
    Filter features by allowed types and optional include/exclude token lists.

    Parameters
    ----------
    features : Iterable[SeqFeature]
        Input features (e.g., ``record.features``).
    include_tokens : sequence of str or None, optional
        If provided, only keep features whose *display name* contains at least
        one of these tokens (case-insensitive).
    exclude_tokens : sequence of str or None, optional
        If provided, drop features whose *display name* contains any of these tokens.
    allowed_types : sequence of str, optional
        Whitelisted feature types (as they appear in the GenBank file).

    Returns
    -------
    list of SeqFeature
        Filtered features, sorted by genomic start position.
    """
    out: list[SeqFeature] = []
    allowed = set(allowed_types)
    for feat in features:
        ftype = str(getattr(feat, "type", ""))
        if ftype not in allowed:
            continue
        name = feature_display_name(feat)
        name_l = name.lower()
        if include_tokens is not None and not _token_hit(name_l, include_tokens):
            continue
        if exclude_tokens is not None and _token_hit(name_l, exclude_tokens):
            continue
        out.append(feat)

    out.sort(key=lambda f: (int(f.location.start), int(f.location.end)))
    return out


def dedupe_record_features(
    record: SeqRecord,
    *,
    tolerance: int = 0,
    include_type: set[str] | None = None,
) -> SeqRecord:
    """
    Remove near-duplicate features **in-place** on a SeqRecord.

    Two features are considered duplicates if their **type**, **strand**, and **label**
    match, and both coordinates (start, end) differ by at most ``tolerance`` bases.

    Parameters
    ----------
    record : Bio.SeqRecord.SeqRecord
        Input record whose ``.features`` will be filtered.
    tolerance : int, optional
        Maximum allowed deviation (bp) for start/end, default is 0.
    include_type : set[str] | None, optional
        If set, only features of these types are considered/retained.

    Returns
    -------
    Bio.SeqRecord.SeqRecord
        The same record with deduplicated features (``record.features`` replaced).
    """
    seen: list[tuple[int, int, int, str, str]] = []
    uniq: list[SeqFeature] = []

    def _key(f: SeqFeature) -> tuple[int, int, int, str, str]:
        start = int(f.location.start)
        end = int(f.location.end)
        strand = int(getattr(f.location, "strand", 0) or 0)
        ftype = str(getattr(f, "type", "feature"))
        label = (f.qualifiers.get("label") or f.qualifiers.get("gene") or [""])[0]
        return start, end, strand, ftype, label

    for feat in record.features:
        if include_type is not None and feat.type not in include_type:
            continue
        s, e, st, ft, lb = _key(feat)
        is_dup = False
        for s0, e0, st0, ft0, lb0 in seen:
            if (
                ft0 == ft
                and st0 == st
                and lb0 == lb
                and abs(s - s0) <= tolerance
                and abs(e - e0) <= tolerance
            ):
                is_dup = True
                break
        if not is_dup:
            seen.append((s, e, st, ft, lb))
            uniq.append(feat)

    record.features = uniq
    return record


def build_feature_concat_sequence(
    record: SeqRecord,
    *,
    include_tokens: Optional[Sequence[str]] = None,
    exclude_tokens: Optional[Sequence[str]] = None,
    allowed_types: Sequence[str] = DEFAULT_ALLOWED_TYPES,
) -> tuple[str, list[tuple[int, int, str, str]]]:
    """
    Build a feature-concatenated reference string and a back-map to features.

    For each filtered feature, the sequence is extracted via
    ``feature.extract(record.seq)`` and appended in genomic order. The function
    returns both the concatenated sequence and a map that records the
    concatenation coordinates and metadata for each block.

    Parameters
    ----------
    record : SeqRecord
        GenBank record with features.
    include_tokens : sequence of str or None, optional
        Keep only features whose name contains at least one token.
    exclude_tokens : sequence of str or None, optional
        Drop features whose name contains any of these tokens (e.g., "backbone", "vector", "ori").
    allowed_types : sequence of str, optional
        Feature types to consider.

    Returns
    -------
    concat_seq : str
        Concatenated feature sequence (uppercase).
    feature_map : list of tuple
        List of tuples ``(start, end, feature_type, feature_name)`` where
        ``start`` and ``end`` are half-open coordinates on the concatenation band.

    Notes
    -----
    - Strand direction is respected by Biopython's extractor.
    - If the filter excludes everything, an empty string and empty map are returned.
    """
    concat_chunks: list[str] = []
    fmap: list[tuple[int, int, str, str]] = []

    feats = filter_features(
        record.features,
        include_tokens=include_tokens,
        exclude_tokens=exclude_tokens,
        allowed_types=allowed_types,
    )

    cursor = 0
    for feat in feats:
        try:
            frag = str(feat.extract(record.seq)).upper()
        except Exception as exc:  # pylint: disable=broad-except
            # Be robust against odd/partial features; skip on failure.
            _ = exc  # keep variable for potential logging hooks
            continue
        start = cursor
        end = cursor + len(frag)
        concat_chunks.append(frag)
        fmap.append(
            (start, end, str(getattr(feat, "type", "")), feature_display_name(feat))
        )
        cursor = end

    return "".join(concat_chunks), fmap


def derive_construct_from_record(
    record: SeqRecord,
    lex: Optional[Lexicon] = None,
    *,
    include_slots: Sequence[str] = ("promoter", "rbs", "gene", "term"),
    treat_5utr_as_rbs: bool = True,
) -> str:
    """
    Derive a "Promoter_RBS_Gene_Term" string from a GenBank record.

    The function iterates features in genomic order and grabs at most one token
    per slot. If a slot is missing (e.g., no terminator was annotated), it is
    simply omitted from the construct string.

    Parameters
    ----------
    record : SeqRecord
        GenBank record.
    lex : Lexicon or None
        Parts lexicon; defaults to `default_lexicon()`.
    include_slots : sequence of str, optional
        Which slots to include and in which order.
    treat_5utr_as_rbs : bool, optional
        If True, treat 5'UTR annotations as RBS.

    Returns
    -------
    str
        Construct string containing the tokens in the given slot order.
    """
    use_lex = lex or default_lexicon()
    feats: list[SeqFeature] = sorted(
        record.features, key=lambda f: (int(f.location.start), int(f.location.end))
    )
    found: dict[str, str] = {}

    for feat in feats:
        slot, token = classify_feature_slot(
            feat, use_lex, treat_5utr_as_rbs=treat_5utr_as_rbs
        )
        if slot is None or token is None:
            continue
        if slot not in include_slots:
            continue
        # Keep the first occurrence per slot (typical MoClo layout is unique per slot).
        if slot not in found:
            found[slot] = token

    parts = [found[s] for s in include_slots if s in found]
    return "_".join(parts)


def derive_construct_from_labels(
    labels: Sequence[str],
    lex: Optional[Lexicon] = None,
    *,
    treat_5utr_as_rbs: bool = True,
) -> str:
    """
    Derive a "Promoter_RBS_Gene_Term" string from arbitrary labels.

    Labels can be "type:token" (e.g., "misc_feature:J23100") or just a token.
    5'UTR-like labels are mapped to the RBS slot if requested.

    Parameters
    ----------
    labels : sequence of str
        Free-form labels such as from feature qualifiers.
    lex : Lexicon or None
        Parts lexicon; defaults to `default_lexicon()`.
    treat_5utr_as_rbs : bool, optional
        If True, treat UTR-like labels as RBS.

    Returns
    -------
    str
        Construct string with the discovered tokens, joined by underscores.
    """
    use_lex = lex or default_lexicon()
    prom: Optional[str] = None
    rbs: Optional[str] = None
    gene: Optional[str] = None
    term: Optional[str] = None

    for lab in labels:
        raw = lab.split(":", 1)[-1].strip() if ":" in lab else lab.strip()
        # Name-based matching
        prom = prom or match_token(raw, "promoter", use_lex)
        gene = gene or match_token(raw, "gene", use_lex)
        term = term or match_token(raw, "term", use_lex)

        if treat_5utr_as_rbs and _is_5utr_label(lab):
            rbs = rbs or match_token(raw, "rbs", use_lex)
        else:
            rbs = rbs or match_token(raw, "rbs", use_lex)

    parts = [p for p in (prom, rbs, gene, term) if p]
    return "_".join(parts)


def derive_construct_from_filename(
    fname: str,
    lex: Optional[Lexicon] = None,
) -> Optional[str]:
    """
    Fallback: derive a construct string from a filename.

    Parameters
    ----------
    fname : str
        File name such as "J23100_AB_BCD8_ecPanD_B0015.gb".
    lex : Lexicon or None
        Parts lexicon; defaults to `default_lexicon()`.

    Returns
    -------
    str or None
        Construct if any parts were matched; otherwise None.
    """
    use_lex = lex or default_lexicon()
    base = os.path.basename(fname)
    prom = match_token(base, "promoter", use_lex)
    rbs = match_token(base, "rbs", use_lex)
    gene = match_token(base, "gene", use_lex)
    term = match_token(base, "term", use_lex)
    parts = [p for p in (prom, rbs, gene, term) if p]
    return "_".join(parts) if parts else None


def parts_key_from_construct(
    construct: str,
    compare_parts: Sequence[str] = ("promoter", "rbs", "gene", "term"),
    lex: Optional[Lexicon] = None,
) -> str:
    """
    Build a comparison key from a construct using selected slots only.

    This is useful, e.g., to compare only promoter+RBS while ignoring gene/term.

    Parameters
    ----------
    construct : str
        Full construct (e.g., "J23100_BCD12_ecPanD_B0015").
    compare_parts : sequence of str, optional
        Which slots to include in the key.
    lex : Lexicon or None
        Parts lexicon; defaults to `default_lexicon()`.

    Returns
    -------
    str
        Reduced key such as "J23100_BCD12".
    """
    if not construct:
        return ""
    use_lex = lex or default_lexicon()

    tokens = [t.strip() for t in construct.split("_") if t.strip()]
    slots: dict[str, Optional[str]] = {
        "promoter": None,
        "rbs": None,
        "gene": None,
        "term": None,
    }

    # Try to assign each segment to a slot (first matching slot wins).
    for seg in tokens:
        for slot in ("promoter", "rbs", "gene", "term"):
            if slots[slot] is None:
                tok = match_token(seg, slot, use_lex)
                if tok is not None:
                    slots[slot] = tok

    return "_".join(v for p in compare_parts if (v := slots[p]))


def inventory_features(
    records: Iterable[SeqRecord],
    lex: Optional[Lexicon] = None,
) -> dict[str, list[str]]:
    """
    Build an inventory of unique tokens per slot across multiple records.

    Parameters
    ----------
    records : iterable of SeqRecord
        One or more GenBank records.
    lex : Lexicon or None
        Parts lexicon; defaults to `default_lexicon()`.

    Returns
    -------
    dict
        Mapping: slot → sorted list of unique tokens (e.g., {"promoter": ["J23100", ...], ...}).
    """
    use_lex = lex or default_lexicon()
    seen: dict[str, dict[str, None]] = {
        "promoter": {},
        "rbs": {},
        "gene": {},
        "term": {},
    }

    for rec in records:
        for feat in rec.features:
            slot, token = classify_feature_slot(feat, use_lex)
            if slot is None or token is None:
                continue
            seen[slot][token] = None

    return {slot: sorted(tokens.keys(), key=str.lower) for slot, tokens in seen.items()}
